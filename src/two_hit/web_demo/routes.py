"""Routes for the two-hit demo: landing, preset cases, bounded upload."""

from __future__ import annotations

import logging
from pathlib import Path

from fastapi import APIRouter, File, Request, UploadFile
from fastapi.responses import HTMLResponse

logger = logging.getLogger(__name__)

router = APIRouter()

CASES_DIR = Path(__file__).parent / "cases"
MAX_BYTES = 1_000_000  # 1 MB per file

# name -> (human label, short description)
CASES: dict[str, tuple[str, str]] = {
    "tp53_loh": (
        "TP53 biallelic (mutation + LOH)",
        "Truncating TP53 mutation with high VAF over a chr17p deletion.",
    ),
    "kras_mut_only": (
        "KRAS mutation only",
        "Activating KRAS G12D with no supporting copy-number change.",
    ),
    "egfr_amp_mut": (
        "EGFR mutation + amplification",
        "EGFR L858R on an amplified chr7 segment.",
    ),
    "brca_compound_het": (
        "BRCA1 compound heterozygous",
        "Two truncating BRCA1 mutations (compound het).",
    ),
}


def _render_report(request: Request, maf_src, seg_src, note: str | None = None):
    """Run the panel-restricted analysis and render the report."""
    from ..analysis import analyze_sample
    from ..io import parse_maf, parse_seg
    from ..models import AnalysisParams
    from .panel import PANEL, filter_maf_to_panel
    from .svgplot import genome_svg

    maf_df = filter_maf_to_panel(parse_maf(maf_src))
    seg_df = parse_seg(seg_src)

    sample_id = None
    if not maf_df.is_empty() and "Tumor_Sample_Barcode" in maf_df.columns:
        ids = maf_df["Tumor_Sample_Barcode"].unique().to_list()
        sample_id = ids[0]
        if len(ids) > 1:
            maf_df = maf_df.filter(pl_eq("Tumor_Sample_Barcode", sample_id))
            note = (note + " " if note else "") + f"Multiple samples found; showing {sample_id}."

    result = analyze_sample(
        maf_df, seg_df, sample_id=sample_id, params=AnalysisParams(), gene_roles=PANEL
    )
    svg = genome_svg(seg_df, result)

    templates = request.app.state.templates
    return templates.TemplateResponse(
        request, "report.html", {"r": result, "plot_svg": svg, "note": note}
    )


def pl_eq(col: str, val):
    """Tiny helper to avoid importing polars at module load."""
    import polars as pl

    return pl.col(col) == val


def _index(request: Request, error: str | None = None, status_code: int = 200):
    templates = request.app.state.templates
    return templates.TemplateResponse(
        request, "index.html", {"cases": CASES, "error": error}, status_code=status_code
    )


@router.get("/", response_class=HTMLResponse)
async def index(request: Request):
    return _index(request)


@router.get("/case/{name}", response_class=HTMLResponse)
async def case(request: Request, name: str):
    if name not in CASES:
        return _index(request, error=f"Unknown case: {name}", status_code=404)
    maf_path = CASES_DIR / f"{name}.maf"
    seg_path = CASES_DIR / f"{name}.seg"
    if not maf_path.exists() or not seg_path.exists():
        return _index(request, error="Case data missing.", status_code=404)
    try:
        return _render_report(request, maf_path, seg_path)
    except Exception as e:  # noqa: BLE001
        logger.error("Case %s failed: %s", name, e, exc_info=True)
        return _index(request, error=str(e), status_code=400)


@router.post("/analyze", response_class=HTMLResponse)
async def analyze(
    request: Request,
    maf_file: UploadFile = File(...),
    seg_file: UploadFile = File(...),
):
    maf_bytes = await maf_file.read()
    seg_bytes = await seg_file.read()
    if len(maf_bytes) > MAX_BYTES or len(seg_bytes) > MAX_BYTES:
        return _index(
            request,
            error=(
                "File too large — the demo accepts files up to 1 MB each."
                " Run the full app locally for larger inputs."
            ),
            status_code=400,
        )
    try:
        maf_src = maf_bytes.decode("utf-8")
        seg_src = seg_bytes.decode("utf-8")
        return _render_report(request, maf_src, seg_src)
    except Exception as e:  # noqa: BLE001
        logger.error("Upload analysis failed: %s", e, exc_info=True)
        return _index(request, error=str(e), status_code=400)
