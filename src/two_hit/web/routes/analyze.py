"""Analysis routes: upload, analyze, report, PDF download."""

from __future__ import annotations

import logging
from pathlib import Path

from fastapi import APIRouter, File, Form, Request, UploadFile
from fastapi.responses import HTMLResponse, JSONResponse, StreamingResponse

from ...analysis import analyze_sample
from ...io import parse_maf, parse_seg
from ...models import AnalysisParams
from ...plot import create_cna_plot, fig_to_html_div
from ...report import result_to_dict
from ..pdf import generate_pdf

logger = logging.getLogger(__name__)

router = APIRouter()


def _find_demo_dir() -> Path:
    """Locate demo/ directory in development or container."""
    candidates = [
        Path(__file__).resolve().parent.parent.parent.parent.parent / "demo",
        Path("/app/demo"),
        Path.cwd() / "demo",
    ]
    for p in candidates:
        if p.is_dir():
            return p
    return candidates[0]


@router.get("/", response_class=HTMLResponse)
async def index(request: Request):
    """Upload page."""
    templates = request.app.state.templates
    return templates.TemplateResponse(request, "index.html")


@router.get("/methods", response_class=HTMLResponse)
async def methods(request: Request):
    """Methods page."""
    templates = request.app.state.templates
    return templates.TemplateResponse(request, "methods.html")


@router.post("/demo", response_class=HTMLResponse)
async def demo(request: Request):
    """Run analysis on the built-in demo case."""
    templates = request.app.state.templates
    results_store = request.app.state.results

    demo_dir = _find_demo_dir()
    maf_path = demo_dir / "demo_sample.maf"
    seg_path = demo_dir / "demo_sample.seg"

    if not maf_path.exists() or not seg_path.exists():
        return templates.TemplateResponse(
            request,
            "index.html",
            {"error": "Demo files not found."},
            status_code=404,
        )

    try:
        maf_df = parse_maf(maf_path)
        seg_df = parse_seg(seg_path)
        params = AnalysisParams()
        result = analyze_sample(maf_df, seg_df, params=params)

        try:
            fig = create_cna_plot(seg_df, result)
            plot_html = fig_to_html_div(fig)
        except Exception:
            logger.warning("Could not generate CNA plot", exc_info=True)
            plot_html = "<p>CNA plot unavailable.</p>"

        result_data = result_to_dict(result)
        result_data["_plot_html"] = plot_html
        rid = results_store.put(result_data)

        return templates.TemplateResponse(
            request,
            "report.html",
            {"r": result_data, "result_id": rid, "plot_html": plot_html},
        )
    except Exception as e:
        logger.error("Demo analysis failed: %s", e, exc_info=True)
        return templates.TemplateResponse(
            request,
            "index.html",
            {"error": str(e)},
            status_code=400,
        )


@router.post("/analyze", response_class=HTMLResponse)
async def analyze(
    request: Request,
    maf_file: UploadFile = File(...),
    seg_file: UploadFile = File(...),
    sample_id: str = Form(""),
    del_threshold: float = Form(-0.3),
    deep_del_threshold: float = Form(-0.5),
    gain_threshold: float = Form(0.3),
    amp_threshold: float = Form(0.7),
    vaf_threshold: float = Form(0.7),
    include_splice_region: bool = Form(False),
):
    """Run analysis on uploaded files and render report."""
    templates = request.app.state.templates
    results_store = request.app.state.results

    try:
        maf_content = (await maf_file.read()).decode("utf-8")
        seg_content = (await seg_file.read()).decode("utf-8")

        maf_df = parse_maf(maf_content)
        seg_df = parse_seg(seg_content)

        params = AnalysisParams(
            del_threshold=del_threshold,
            deep_del_threshold=deep_del_threshold,
            gain_threshold=gain_threshold,
            amp_threshold=amp_threshold,
            vaf_threshold=vaf_threshold,
            include_splice_region=include_splice_region,
        )

        sid = sample_id.strip() or None
        result = analyze_sample(maf_df, seg_df, sample_id=sid, params=params)

        # Generate CNA plot HTML
        try:
            fig = create_cna_plot(seg_df, result)
            plot_html = fig_to_html_div(fig)
        except Exception:
            logger.warning("Could not generate CNA plot", exc_info=True)
            plot_html = "<p>CNA plot unavailable.</p>"

        result_data = result_to_dict(result)
        result_data["_plot_html"] = plot_html
        rid = results_store.put(result_data)

        return templates.TemplateResponse(
            request,
            "report.html",
            {"r": result_data, "result_id": rid, "plot_html": plot_html},
        )
    except Exception as e:
        logger.error("Analysis failed: %s", e, exc_info=True)
        return templates.TemplateResponse(
            request,
            "index.html",
            {"error": str(e)},
            status_code=400,
        )


@router.get("/results/{result_id}/pdf")
async def download_pdf(request: Request, result_id: str):
    """Download PDF report for a stored result."""
    result_data = request.app.state.results.get(result_id)
    if result_data is None:
        return HTMLResponse("Result not found or expired.", status_code=404)

    # Re-generate plot as PNG for PDF
    plot_png = None
    try:
        # We can't easily re-create the plot without the original seg_df,
        # so we embed the plot PNG only if available
        pass
    except Exception:
        pass

    pdf_bytes = generate_pdf(result_data, plot_png=plot_png)
    sid = result_data.get("sample_id", "sample")
    return StreamingResponse(
        iter([pdf_bytes]),
        media_type="application/pdf",
        headers={"Content-Disposition": f'attachment; filename="{sid}_two_hit.pdf"'},
    )


@router.get("/results/{result_id}/json")
async def result_json(request: Request, result_id: str):
    """Return stored result as JSON."""
    result_data = request.app.state.results.get(result_id)
    if result_data is None:
        return JSONResponse({"error": "Result not found or expired."}, status_code=404)
    # Strip internal keys
    clean = {k: v for k, v in result_data.items() if not k.startswith("_")}
    return JSONResponse(clean)


@router.post("/api/analyze")
async def api_analyze(
    maf_file: UploadFile = File(...),
    seg_file: UploadFile = File(...),
    sample_id: str = Form(""),
    del_threshold: float = Form(-0.3),
    deep_del_threshold: float = Form(-0.5),
    gain_threshold: float = Form(0.3),
    amp_threshold: float = Form(0.7),
    vaf_threshold: float = Form(0.7),
    include_splice_region: bool = Form(False),
):
    """JSON API endpoint for programmatic access."""
    try:
        maf_content = (await maf_file.read()).decode("utf-8")
        seg_content = (await seg_file.read()).decode("utf-8")

        maf_df = parse_maf(maf_content)
        seg_df = parse_seg(seg_content)

        params = AnalysisParams(
            del_threshold=del_threshold,
            deep_del_threshold=deep_del_threshold,
            gain_threshold=gain_threshold,
            amp_threshold=amp_threshold,
            vaf_threshold=vaf_threshold,
            include_splice_region=include_splice_region,
        )

        sid = sample_id.strip() or None
        result = analyze_sample(maf_df, seg_df, sample_id=sid, params=params)
        return JSONResponse(result_to_dict(result))
    except Exception as e:
        return JSONResponse({"error": str(e)}, status_code=400)
