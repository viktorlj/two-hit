"""MAF and SEG file parsers using polars.

Provides both eager (small files / web uploads) and lazy (GENIE-scale) modes.
"""

from __future__ import annotations

import logging
from io import BytesIO
from pathlib import Path

import polars as pl

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Required columns
# ---------------------------------------------------------------------------

REQUIRED_MAF_COLUMNS = {
    "Hugo_Symbol",
    "Chromosome",
    "Start_Position",
    "End_Position",
    "Variant_Classification",
    "HGVSp_Short",
    "Tumor_Sample_Barcode",
}

MAF_SELECT_COLUMNS = [
    "Hugo_Symbol",
    "Chromosome",
    "Start_Position",
    "End_Position",
    "Variant_Classification",
    "HGVSp_Short",
    "Tumor_Sample_Barcode",
    "t_alt_count",
    "t_ref_count",
]

REQUIRED_SEG_COLUMNS = {"ID", "chrom", "loc.start", "loc.end", "seg.mean"}

# Column renaming for SEG normalization
SEG_RENAME: dict[str, str] = {
    "ID": "sample_id",
    "chrom": "chromosome",
    "loc.start": "start",
    "loc.end": "end",
    "num.mark": "num_mark",
    "seg.mean": "seg_mean",
}

# Case-insensitive alias map for SEG columns
_SEG_ALIASES: dict[str, str] = {
    "id": "ID",
    "sample": "ID",
    "sample_id": "ID",
    "chromosome": "chrom",
    "chr": "chrom",
    "loc.start": "loc.start",
    "start": "loc.start",
    "loc_start": "loc.start",
    "loc.end": "loc.end",
    "end": "loc.end",
    "loc_end": "loc.end",
    "num.mark": "num.mark",
    "num_mark": "num.mark",
    "markers": "num.mark",
    "seg.mean": "seg.mean",
    "seg_mean": "seg.mean",
    "segmean": "seg.mean",
}


def _strip_comments(content: str) -> str:
    """Remove lines starting with # from content."""
    lines = [line for line in content.splitlines(keepends=True) if not line.startswith("#")]
    return "".join(lines)


def _validate_columns(df: pl.DataFrame, required: set[str], file_type: str) -> None:
    """Raise ValueError if required columns are missing."""
    present = set(df.columns)
    missing = required - present
    if missing:
        raise ValueError(f"Missing required {file_type} columns: {sorted(missing)}")


def _normalize_seg_columns(df: pl.DataFrame) -> pl.DataFrame:
    """Normalize SEG column names to standard format, then rename."""
    rename_map: dict[str, str] = {}
    lower_to_orig = {c.lower(): c for c in df.columns}

    for alias, canonical in _SEG_ALIASES.items():
        if alias in lower_to_orig and canonical not in rename_map.values():
            rename_map[lower_to_orig[alias]] = canonical

    if rename_map:
        df = df.rename(rename_map)

    # Now apply the standard rename to friendly names
    final_rename = {k: v for k, v in SEG_RENAME.items() if k in df.columns}
    return df.rename(final_rename)


def _coerce_maf_types(df: pl.DataFrame) -> pl.DataFrame:
    """Ensure MAF numeric columns have correct types."""
    casts = {}
    for col in ("Start_Position", "End_Position"):
        if col in df.columns and df[col].dtype != pl.Int64:
            casts[col] = pl.Int64
    for col in ("t_alt_count", "t_ref_count"):
        if col in df.columns and df[col].dtype != pl.Int64:
            casts[col] = pl.Int64
    if casts:
        df = df.cast(casts, strict=False)
    return df


def _coerce_seg_types(df: pl.DataFrame) -> pl.DataFrame:
    """Ensure SEG numeric columns have correct types."""
    casts = {}
    for col in ("start", "end"):
        if col in df.columns and df[col].dtype != pl.Int64:
            casts[col] = pl.Int64
    if "seg_mean" in df.columns and df["seg_mean"].dtype != pl.Float64:
        casts["seg_mean"] = pl.Float64
    if casts:
        df = df.cast(casts, strict=False)
    return df


# ---------------------------------------------------------------------------
# Eager mode — small files / web uploads
# ---------------------------------------------------------------------------


def parse_maf(source: Path | str) -> pl.DataFrame:
    """Parse a MAF file from a file path or string content.

    Handles comment lines starting with #, validates required columns,
    and coerces numeric types.
    """
    if isinstance(source, Path):
        content = source.read_text()
    else:
        content = source

    cleaned = _strip_comments(content)
    if not cleaned.strip():
        raise ValueError("Empty MAF content after removing comments")

    # Detect separator
    first_line = cleaned.split("\n", 1)[0]
    separator = "\t" if "\t" in first_line else ","

    df = pl.read_csv(BytesIO(cleaned.encode()), separator=separator, infer_schema_length=0)
    _validate_columns(df, REQUIRED_MAF_COLUMNS, "MAF")

    # Select available columns from our desired set
    available = [c for c in MAF_SELECT_COLUMNS if c in df.columns]
    df = df.select(available)
    return _coerce_maf_types(df)


def parse_seg(source: Path | str) -> pl.DataFrame:
    """Parse a SEG file from a file path or string content.

    Normalizes column names (case-insensitive), validates required columns,
    and coerces numeric types.
    """
    if isinstance(source, Path):
        content = source.read_text()
    else:
        content = source

    cleaned = _strip_comments(content)
    if not cleaned.strip():
        raise ValueError("Empty SEG content after removing comments")

    first_line = cleaned.split("\n", 1)[0]
    separator = "\t" if "\t" in first_line else ","

    df = pl.read_csv(BytesIO(cleaned.encode()), separator=separator, infer_schema_length=0)
    df = _normalize_seg_columns(df)

    required_normalized = {"sample_id", "chromosome", "start", "end", "seg_mean"}
    _validate_columns(df, required_normalized, "SEG")

    return _coerce_seg_types(df)


# ---------------------------------------------------------------------------
# Lazy mode — large files (GENIE-scale)
# ---------------------------------------------------------------------------


def scan_maf(path: Path) -> pl.LazyFrame:
    """Lazily scan a MAF file, selecting only needed columns."""
    lf = pl.scan_csv(
        path,
        separator="\t",
        infer_schema_length=10000,
        schema_overrides={
            "Chromosome": pl.Utf8,
            "Start_Position": pl.Int64,
            "End_Position": pl.Int64,
            "t_alt_count": pl.Int64,
            "t_ref_count": pl.Int64,
        },
    )
    schema_cols = lf.collect_schema().names()
    available = [c for c in MAF_SELECT_COLUMNS if c in schema_cols]
    return lf.select(available)


def scan_seg(path: Path) -> pl.LazyFrame:
    """Lazily scan a SEG file."""
    return pl.scan_csv(
        path,
        separator="\t",
        infer_schema_length=10000,
        schema_overrides={"chrom": pl.Utf8},
    )


def filter_maf_for_sample(lf: pl.LazyFrame, sample_id: str) -> pl.DataFrame:
    """Filter a lazy MAF frame to a single sample and collect."""
    df = lf.filter(pl.col("Tumor_Sample_Barcode") == sample_id).collect()
    return _coerce_maf_types(df)


def filter_seg_for_sample(lf: pl.LazyFrame, sample_id: str) -> pl.DataFrame:
    """Filter a lazy SEG frame to a single sample, normalize, and collect."""
    df = lf.filter(pl.col("ID") == sample_id).collect()
    df = _normalize_seg_columns(df)
    return _coerce_seg_types(df)
