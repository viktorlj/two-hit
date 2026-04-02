"""CLI entry point using typer."""

from __future__ import annotations

import logging
from pathlib import Path

import typer

from . import __version__

app = typer.Typer(
    name="two-hit",
    help="Integrated mutation + copy number biallelic inactivation analysis.",
    no_args_is_help=True,
)

logger = logging.getLogger("two_hit")


def _setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(levelname)s  %(message)s",
    )


def _build_params(
    del_threshold: float,
    deep_del_threshold: float,
    gain_threshold: float,
    amp_threshold: float,
    vaf_threshold: float,
    include_splice_region: bool,
):
    from .models import AnalysisParams

    return AnalysisParams(
        del_threshold=del_threshold,
        deep_del_threshold=deep_del_threshold,
        gain_threshold=gain_threshold,
        amp_threshold=amp_threshold,
        vaf_threshold=vaf_threshold,
        include_splice_region=include_splice_region,
    )


@app.command()
def report(
    maf: Path = typer.Argument(..., help="Path to MAF file"),
    seg: Path = typer.Argument(..., help="Path to SEG file"),
    sample_id: str | None = typer.Option(None, "--sample-id", "-s", help="Sample identifier"),
    output: Path | None = typer.Option(None, "--output", "-o", help="Output file path"),
    fmt: str = typer.Option("text", "--format", "-f", help="Output format: text, json, tsv, pdf"),
    del_threshold: float = typer.Option(-0.3, "--del-threshold", help="LOH deletion threshold"),
    deep_del_threshold: float = typer.Option(
        -0.5, "--deep-del-threshold", help="Homozygous deletion threshold"
    ),
    gain_threshold: float = typer.Option(0.3, "--gain-threshold", help="CN gain threshold"),
    amp_threshold: float = typer.Option(0.7, "--amp-threshold", help="Amplification threshold"),
    vaf_threshold: float = typer.Option(0.7, "--vaf-threshold", help="High VAF threshold"),
    include_splice_region: bool = typer.Option(
        False, "--include-splice-region", help="Count Splice_Region as truncating"
    ),
    verbose: bool = typer.Option(False, "--verbose", "-v"),
) -> None:
    """Analyze a MAF + SEG file pair and generate a two-hit report."""
    _setup_logging(verbose)

    from .analysis import analyze_sample
    from .io import parse_maf, parse_seg
    from .plot import create_cna_plot, fig_to_png_bytes
    from .report import format_human_readable, gene_results_to_tsv, result_to_dict, result_to_json
    from .web.pdf import generate_pdf

    params = _build_params(
        del_threshold,
        deep_del_threshold,
        gain_threshold,
        amp_threshold,
        vaf_threshold,
        include_splice_region,
    )

    typer.echo(f"Reading MAF: {maf}")
    maf_df = parse_maf(maf)
    typer.echo(f"Reading SEG: {seg}")
    seg_df = parse_seg(seg)
    typer.echo(f"Loaded {len(maf_df)} mutations, {len(seg_df)} segments")

    result = analyze_sample(maf_df, seg_df, sample_id=sample_id, params=params)

    if fmt == "json":
        text = result_to_json(result)
    elif fmt == "tsv":
        text = gene_results_to_tsv(result.gene_results)
    elif fmt == "pdf":
        result_dict = result_to_dict(result)
        try:
            fig = create_cna_plot(seg_df, result)
            plot_png = fig_to_png_bytes(fig)
        except Exception:
            logger.warning("Could not generate CNA plot for PDF")
            plot_png = None
        pdf_bytes = generate_pdf(result_dict, plot_png=plot_png)
        out_path = output or Path(f"{result.sample_id}_two_hit.pdf")
        out_path.write_bytes(pdf_bytes)
        typer.echo(f"PDF report written to {out_path}")
        return
    else:
        text = format_human_readable(result)

    if output:
        output.write_text(text)
        typer.echo(f"Report written to {output}")
    else:
        typer.echo(text)


@app.command()
def genie(
    genie_dir: Path = typer.Argument(..., help="Path to GENIE data directory"),
    sample: str = typer.Argument(..., help="Tumor_Sample_Barcode / ID"),
    output: Path | None = typer.Option(None, "--output", "-o", help="Output file path"),
    fmt: str = typer.Option("text", "--format", "-f", help="Output format: text, json, tsv, pdf"),
    del_threshold: float = typer.Option(-0.3, "--del-threshold"),
    deep_del_threshold: float = typer.Option(-0.5, "--deep-del-threshold"),
    gain_threshold: float = typer.Option(0.3, "--gain-threshold"),
    amp_threshold: float = typer.Option(0.7, "--amp-threshold"),
    vaf_threshold: float = typer.Option(0.7, "--vaf-threshold"),
    include_splice_region: bool = typer.Option(False, "--include-splice-region"),
    verbose: bool = typer.Option(False, "--verbose", "-v"),
) -> None:
    """Analyze a sample from GENIE data using lazy scanning."""
    _setup_logging(verbose)

    from .analysis import analyze_sample
    from .io import filter_maf_for_sample, filter_seg_for_sample, scan_maf, scan_seg
    from .plot import create_cna_plot, fig_to_png_bytes
    from .report import format_human_readable, gene_results_to_tsv, result_to_dict, result_to_json
    from .web.pdf import generate_pdf

    params = _build_params(
        del_threshold,
        deep_del_threshold,
        gain_threshold,
        amp_threshold,
        vaf_threshold,
        include_splice_region,
    )

    # Find MAF and SEG files
    maf_path = genie_dir / "data_mutations_extended.txt"
    seg_candidates = ["data_cna_hg19.seg", "genie_data_cna_hg19.seg"]
    seg_path = None
    for name in seg_candidates:
        candidate = genie_dir / name
        if candidate.exists():
            seg_path = candidate
            break
    if seg_path is None:
        typer.echo(f"Error: No SEG file found in {genie_dir}", err=True)
        raise typer.Exit(1)

    typer.echo(f"Scanning MAF: {maf_path.name} (lazy)")
    maf_lf = scan_maf(maf_path)
    typer.echo(f"Scanning SEG: {seg_path.name} (lazy)")
    seg_lf = scan_seg(seg_path)

    typer.echo(f"Filtering for sample: {sample}")
    maf_df = filter_maf_for_sample(maf_lf, sample)
    seg_df = filter_seg_for_sample(seg_lf, sample)
    typer.echo(f"Found {len(maf_df)} mutations, {len(seg_df)} segments")

    if maf_df.is_empty() and seg_df.is_empty():
        typer.echo(f"Error: No data found for sample {sample}", err=True)
        raise typer.Exit(1)

    result = analyze_sample(maf_df, seg_df, sample_id=sample, params=params)

    if fmt == "json":
        text = result_to_json(result)
    elif fmt == "tsv":
        text = gene_results_to_tsv(result.gene_results)
    elif fmt == "pdf":
        result_dict = result_to_dict(result)
        try:
            fig = create_cna_plot(seg_df, result)
            plot_png = fig_to_png_bytes(fig)
        except Exception:
            logger.warning("Could not generate CNA plot for PDF")
            plot_png = None
        pdf_bytes = generate_pdf(result_dict, plot_png=plot_png)
        out_path = output or Path(f"{result.sample_id}_two_hit.pdf")
        out_path.write_bytes(pdf_bytes)
        typer.echo(f"PDF report written to {out_path}")
        return
    else:
        text = format_human_readable(result)

    if output:
        output.write_text(text)
        typer.echo(f"Report written to {output}")
    else:
        typer.echo(text)


@app.command()
def serve(
    host: str = typer.Option("0.0.0.0", "--host"),
    port: int | None = typer.Option(None, "--port"),
    reload: bool = typer.Option(False, "--reload"),
) -> None:
    """Start the web interface."""
    import os

    import uvicorn

    from .web.app import create_app

    if port is None:
        port = int(os.environ.get("PORT", "8000"))

    web_app = create_app()
    typer.echo(f"Starting two-hit web interface at http://{host}:{port}")
    uvicorn.run(web_app, host=host, port=port)


@app.callback(invoke_without_command=True)
def main(
    version: bool = typer.Option(False, "--version", "-V", help="Show version"),
) -> None:
    """Two-hit biallelic inactivation analysis tool."""
    if version:
        typer.echo(f"two-hit {__version__}")
        raise typer.Exit()
