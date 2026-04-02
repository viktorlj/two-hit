"""Core two-hit analysis logic.

Integrates somatic mutations (MAF) with copy number segments (SEG) to identify
biallelic inactivation in tumor suppressors and mutation+amplification in oncogenes.
"""

from __future__ import annotations

import logging

import polars as pl

from .gene_list import GeneRole, load_gene_roles
from .models import (
    BIALLELIC_STATUSES,
    AnalysisParams,
    CnaStatus,
    Confidence,
    GeneResult,
    MutationRecord,
    SampleResult,
    TwoHitStatus,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# CNA classification
# ---------------------------------------------------------------------------


def classify_cna(seg_mean: float | None, params: AnalysisParams) -> CnaStatus:
    """Classify a segment's copy number status from its seg.mean value."""
    if seg_mean is None:
        return CnaStatus.NO_DATA
    if seg_mean < params.deep_del_threshold:
        return CnaStatus.HOMDEL
    if seg_mean < params.del_threshold:
        return CnaStatus.HEMIDEL
    if seg_mean > params.amp_threshold:
        return CnaStatus.AMP
    if seg_mean > params.gain_threshold:
        return CnaStatus.GAIN
    return CnaStatus.NEUTRAL


# ---------------------------------------------------------------------------
# Mutation-to-segment mapping
# ---------------------------------------------------------------------------


def _normalize_chrom(chrom: str) -> str:
    """Normalize chromosome name: strip 'chr' prefix, uppercase."""
    c = str(chrom).strip().upper()
    if c.startswith("CHR"):
        c = c[3:]
    return c


def map_mutations_to_segments(
    maf_df: pl.DataFrame,
    seg_df: pl.DataFrame,
    params: AnalysisParams,
) -> list[MutationRecord]:
    """Map each mutation to its overlapping CNA segment.

    For each mutation, finds the segment on the same chromosome where
    Start_Position falls within [start, end]. If multiple segments overlap,
    keeps the one with the most extreme |seg_mean|.

    Returns a list of MutationRecord with seg_mean and cna_status populated.
    """
    if maf_df.is_empty():
        return []

    # Normalize chromosomes in both dataframes
    maf_with_chrom = maf_df.with_columns(
        pl.col("Chromosome")
        .cast(pl.Utf8)
        .map_elements(_normalize_chrom, return_dtype=pl.Utf8)
        .alias("chrom_norm")
    )

    records: list[MutationRecord] = []

    if seg_df.is_empty():
        # No segments — all mutations get NO_DATA
        for row in maf_with_chrom.iter_rows(named=True):
            records.append(_row_to_mutation(row, seg_mean=None, cna_status=CnaStatus.NO_DATA))
        return records

    seg_with_chrom = seg_df.with_columns(
        pl.col("chromosome")
        .cast(pl.Utf8)
        .map_elements(_normalize_chrom, return_dtype=pl.Utf8)
        .alias("chrom_norm")
    )

    # Group segments by chromosome for efficient lookup
    seg_by_chrom: dict[str, pl.DataFrame] = {}
    for chrom_val in seg_with_chrom["chrom_norm"].unique().to_list():
        seg_by_chrom[chrom_val] = seg_with_chrom.filter(pl.col("chrom_norm") == chrom_val)

    for row in maf_with_chrom.iter_rows(named=True):
        chrom = row["chrom_norm"]
        pos = row["Start_Position"]

        seg_mean = None
        cna_status = CnaStatus.NO_DATA

        if chrom in seg_by_chrom and pos is not None:
            chrom_segs = seg_by_chrom[chrom]
            overlapping = chrom_segs.filter((pl.col("start") <= pos) & (pl.col("end") >= pos))
            if not overlapping.is_empty():
                # Pick segment with most extreme |seg_mean|
                seg_means = overlapping["seg_mean"].to_list()
                seg_mean = max(seg_means, key=lambda x: abs(x) if x is not None else 0)
                cna_status = classify_cna(seg_mean, params)

        records.append(_row_to_mutation(row, seg_mean=seg_mean, cna_status=cna_status))

    return records


def _row_to_mutation(row: dict, seg_mean: float | None, cna_status: CnaStatus) -> MutationRecord:
    """Convert a MAF row dict to a MutationRecord."""
    t_alt = row.get("t_alt_count")
    t_ref = row.get("t_ref_count")
    return MutationRecord(
        hugo_symbol=row.get("Hugo_Symbol", ""),
        chromosome=str(row.get("Chromosome", "")),
        start_position=row.get("Start_Position") or 0,
        end_position=row.get("End_Position") or 0,
        variant_classification=row.get("Variant_Classification", ""),
        hgvsp_short=row.get("HGVSp_Short") or "",
        t_alt_count=int(t_alt) if t_alt is not None else None,
        t_ref_count=int(t_ref) if t_ref is not None else None,
        seg_mean=seg_mean,
        cna_status=cna_status,
    )


# ---------------------------------------------------------------------------
# Gene-level classification
# ---------------------------------------------------------------------------


def _classify_tsg(
    mutations: list[MutationRecord],
    params: AnalysisParams,
) -> tuple[TwoHitStatus, Confidence, list[str]]:
    """Classify a tumor suppressor gene."""
    notes: list[str] = []
    truncating = [m for m in mutations if m.is_truncating(params)]
    has_deletion = any(m.cna_status in (CnaStatus.HEMIDEL, CnaStatus.HOMDEL) for m in mutations)
    has_homdel = any(m.cna_status == CnaStatus.HOMDEL for m in mutations)

    # 1. Homozygous deletion
    if has_homdel:
        notes.append("Deep deletion detected")
        return TwoHitStatus.HOMDEL, Confidence.HIGH, notes

    # 2. Compound heterozygous (>=2 truncating mutations)
    if len(truncating) >= 2:
        notes.append(f"{len(truncating)} truncating mutations")
        return TwoHitStatus.BIALLELIC_COMPOUND_HET, Confidence.HIGH, notes

    # 3. Mutation + LOH with high VAF
    if mutations and has_deletion:
        vafs = [m.vaf for m in mutations if m.vaf is not None]
        max_vaf = max(vafs) if vafs else None
        if max_vaf is not None and max_vaf > params.vaf_threshold:
            notes.append(f"Mutation + deletion with high VAF ({max_vaf:.2f})")
            return TwoHitStatus.BIALLELIC_MUT_LOH_HIGH_VAF, Confidence.HIGH, notes
        notes.append("Mutation + deletion (LOH)")
        return TwoHitStatus.BIALLELIC_MUT_LOH, Confidence.MEDIUM, notes

    # 4. Single hit
    if mutations:
        if has_deletion:
            # Shouldn't reach here, but defensive
            return TwoHitStatus.SINGLE_HIT_DEL, Confidence.LOW, notes
        notes.append("Single mutation, no CNA support")
        return TwoHitStatus.SINGLE_HIT_MUT, Confidence.LOW, notes

    # 5. Deletion only (no mutation)
    if has_deletion:
        return TwoHitStatus.SINGLE_HIT_DEL, Confidence.LOW, notes

    return TwoHitStatus.UNCLASSIFIED, Confidence.LOW, notes


def _classify_oncogene(
    mutations: list[MutationRecord],
    params: AnalysisParams,
) -> tuple[TwoHitStatus, Confidence, list[str]]:
    """Classify an oncogene."""
    notes: list[str] = []
    has_amp = any(m.cna_status == CnaStatus.AMP for m in mutations)
    has_gain = any(m.cna_status in (CnaStatus.AMP, CnaStatus.GAIN) for m in mutations)

    if mutations and has_amp:
        notes.append("Mutation + amplification")
        return TwoHitStatus.MUT_AMP, Confidence.HIGH, notes

    if mutations and has_gain:
        notes.append("Mutation + copy number gain")
        return TwoHitStatus.MUT_GAIN, Confidence.MEDIUM, notes

    if mutations:
        notes.append("Mutation without CNA")
        return TwoHitStatus.MUT_ONLY, Confidence.LOW, notes

    # Amplification only (no mutation) — shouldn't normally reach here
    # since we group by genes with mutations, but handle defensively
    if has_amp:
        return TwoHitStatus.AMP_ONLY, Confidence.MEDIUM, notes

    return TwoHitStatus.UNCLASSIFIED, Confidence.LOW, notes


_CONFIDENCE_ORDER = {Confidence.HIGH: 3, Confidence.MEDIUM: 2, Confidence.LOW: 1}


def classify_gene(
    hugo_symbol: str,
    role: GeneRole,
    mutations: list[MutationRecord],
    params: AnalysisParams,
) -> GeneResult:
    """Classify a single gene based on its mutations and CNA context."""
    # Compute gene-level summary values
    seg_means = [m.seg_mean for m in mutations if m.seg_mean is not None]
    gene_seg_mean = max(seg_means, key=abs) if seg_means else None
    gene_cna = classify_cna(gene_seg_mean, params)
    vafs = [m.vaf for m in mutations if m.vaf is not None]
    vaf_max = max(vafs) if vafs else None

    if role == GeneRole.UNKNOWN:
        return GeneResult(
            hugo_symbol=hugo_symbol,
            role=role,
            mutations=mutations,
            seg_mean=gene_seg_mean,
            cna_status=gene_cna,
            two_hit_status=TwoHitStatus.UNCLASSIFIED,
            confidence=Confidence.LOW,
            vaf_max=vaf_max,
            notes=["Gene not in OncoKB — role unknown"],
        )

    if role == GeneRole.TSG:
        status, conf, notes = _classify_tsg(mutations, params)
    elif role == GeneRole.ONCOGENE:
        status, conf, notes = _classify_oncogene(mutations, params)
    elif role == GeneRole.BOTH:
        # Run both classifiers, pick the one with higher confidence
        tsg_status, tsg_conf, tsg_notes = _classify_tsg(mutations, params)
        onc_status, onc_conf, onc_notes = _classify_oncogene(mutations, params)

        tsg_score = _CONFIDENCE_ORDER[tsg_conf]
        onc_score = _CONFIDENCE_ORDER[onc_conf]

        if tsg_score > onc_score:
            status, conf, notes = tsg_status, tsg_conf, tsg_notes
        elif onc_score > tsg_score:
            status, conf, notes = onc_status, onc_conf, onc_notes
        else:
            # Tie-break: prefer biallelic
            if tsg_status in BIALLELIC_STATUSES:
                status, conf, notes = tsg_status, tsg_conf, tsg_notes
            else:
                status, conf, notes = onc_status, onc_conf, onc_notes
        notes.append("Dual-role gene (oncogene + TSG)")
    else:
        status, conf, notes = TwoHitStatus.UNCLASSIFIED, Confidence.LOW, []

    return GeneResult(
        hugo_symbol=hugo_symbol,
        role=role,
        mutations=mutations,
        seg_mean=gene_seg_mean,
        cna_status=gene_cna,
        two_hit_status=status,
        confidence=conf,
        vaf_max=vaf_max,
        notes=notes,
    )


# ---------------------------------------------------------------------------
# Sample-level orchestrator
# ---------------------------------------------------------------------------

# Sort priority: biallelic/mut_amp first, then by confidence, then alphabetical
_STATUS_PRIORITY: dict[TwoHitStatus, int] = {
    TwoHitStatus.BIALLELIC_MUT_LOH_HIGH_VAF: 10,
    TwoHitStatus.BIALLELIC_COMPOUND_HET: 10,
    TwoHitStatus.BIALLELIC_MUT_LOH: 9,
    TwoHitStatus.HOMDEL: 9,
    TwoHitStatus.MUT_AMP: 8,
    TwoHitStatus.MUT_GAIN: 6,
    TwoHitStatus.SINGLE_HIT_MUT: 4,
    TwoHitStatus.SINGLE_HIT_DEL: 4,
    TwoHitStatus.AMP_ONLY: 5,
    TwoHitStatus.MUT_ONLY: 3,
    TwoHitStatus.UNCLASSIFIED: 1,
}


def analyze_sample(
    maf_df: pl.DataFrame,
    seg_df: pl.DataFrame,
    sample_id: str | None = None,
    params: AnalysisParams | None = None,
    gene_roles: dict[str, GeneRole] | None = None,
) -> SampleResult:
    """Run the full two-hit analysis for a single sample.

    Args:
        maf_df: Parsed MAF dataframe (single sample).
        seg_df: Parsed SEG dataframe (single sample).
        sample_id: Sample identifier. If None, inferred from MAF.
        params: Analysis parameters. Defaults to AnalysisParams().
        gene_roles: Gene role mapping. If None, loaded from OncoKB.

    Returns:
        SampleResult with per-gene classifications.
    """
    params = params or AnalysisParams()
    gene_roles = gene_roles if gene_roles is not None else load_gene_roles()

    # Infer sample ID
    if sample_id is None:
        if "Tumor_Sample_Barcode" in maf_df.columns and not maf_df.is_empty():
            sample_id = maf_df["Tumor_Sample_Barcode"][0]
        else:
            sample_id = "unknown"

    logger.info(
        "Analyzing sample %s: %d mutations, %d segments", sample_id, len(maf_df), len(seg_df)
    )

    # Map mutations to segments
    mutation_records = map_mutations_to_segments(maf_df, seg_df, params)

    # Group by gene
    genes: dict[str, list[MutationRecord]] = {}
    for mut in mutation_records:
        genes.setdefault(mut.hugo_symbol, []).append(mut)

    # Classify each gene
    gene_results: list[GeneResult] = []
    for hugo_symbol, muts in genes.items():
        role = gene_roles.get(hugo_symbol, GeneRole.UNKNOWN)
        result = classify_gene(hugo_symbol, role, muts, params)
        gene_results.append(result)

    # Sort: highest priority status first, then confidence, then alphabetical
    gene_results.sort(
        key=lambda g: (
            -_STATUS_PRIORITY.get(g.two_hit_status, 0),
            -_CONFIDENCE_ORDER.get(g.confidence, 0),
            g.hugo_symbol,
        )
    )

    biallelic_count = sum(1 for g in gene_results if g.is_biallelic)
    mut_amp_count = sum(1 for g in gene_results if g.two_hit_status == TwoHitStatus.MUT_AMP)

    return SampleResult(
        sample_id=sample_id,
        params=params,
        gene_results=gene_results,
        total_mutations=len(mutation_records),
        total_genes_with_events=len(gene_results),
        biallelic_count=biallelic_count,
        mut_amp_count=mut_amp_count,
    )
