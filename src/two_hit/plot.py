"""Genome-wide CNA plot with mutation overlay using plotly."""

from __future__ import annotations

from collections import OrderedDict

import plotly.graph_objects as go
import polars as pl

from .models import (
    SampleResult,
    TwoHitStatus,
)

# ---------------------------------------------------------------------------
# Chromosome layout (hg19 / GRCh37)
# ---------------------------------------------------------------------------

# Chromosome sizes in base pairs (GRCh37)
CHROM_SIZES: OrderedDict[str, int] = OrderedDict(
    [
        ("1", 249250621),
        ("2", 243199373),
        ("3", 198022430),
        ("4", 191154276),
        ("5", 180915260),
        ("6", 171115067),
        ("7", 159138663),
        ("8", 146364022),
        ("9", 141213431),
        ("10", 135534747),
        ("11", 135006516),
        ("12", 133851895),
        ("13", 115169878),
        ("14", 107349540),
        ("15", 102531392),
        ("16", 90354753),
        ("17", 81195210),
        ("18", 78077248),
        ("19", 59128983),
        ("20", 63025520),
        ("21", 48129895),
        ("22", 51304566),
        ("X", 155270560),
        ("Y", 59373566),
    ]
)


def get_chrom_offsets() -> dict[str, int]:
    """Compute cumulative genomic offsets for sequential chromosome layout."""
    offsets: dict[str, int] = {}
    cumulative = 0
    for chrom, size in CHROM_SIZES.items():
        offsets[chrom] = cumulative
        cumulative += size
    return offsets


def _normalize_chrom(chrom: str) -> str:
    """Normalize chromosome name for offset lookup."""
    c = str(chrom).strip().upper()
    if c.startswith("CHR"):
        c = c[3:]
    return c


# ---------------------------------------------------------------------------
# Color palette
# ---------------------------------------------------------------------------

COLORS = {
    "del": "#2962ff",
    "amp": "#d32f2f",
    "neutral": "#bdbdbd",
    "biallelic": "#c62828",
    "oncogene_hit": "#1565c0",
    "single_hit": "#757575",
    "chrom_alt": "rgba(0,0,0,0.03)",
    "threshold": "rgba(0,0,0,0.2)",
}

STATUS_COLORS: dict[TwoHitStatus, str] = {
    TwoHitStatus.BIALLELIC_COMPOUND_HET: COLORS["biallelic"],
    TwoHitStatus.BIALLELIC_MUT_LOH: COLORS["biallelic"],
    TwoHitStatus.BIALLELIC_MUT_LOH_HIGH_VAF: COLORS["biallelic"],
    TwoHitStatus.HOMDEL: COLORS["biallelic"],
    TwoHitStatus.MUT_AMP: COLORS["oncogene_hit"],
    TwoHitStatus.MUT_GAIN: COLORS["oncogene_hit"],
    TwoHitStatus.SINGLE_HIT_MUT: COLORS["single_hit"],
    TwoHitStatus.SINGLE_HIT_DEL: COLORS["single_hit"],
    TwoHitStatus.MUT_ONLY: COLORS["single_hit"],
    TwoHitStatus.AMP_ONLY: COLORS["oncogene_hit"],
    TwoHitStatus.UNCLASSIFIED: COLORS["single_hit"],
}


# ---------------------------------------------------------------------------
# Plot construction
# ---------------------------------------------------------------------------


def create_cna_plot(
    seg_df: pl.DataFrame,
    sample_result: SampleResult,
    height: int = 450,
) -> go.Figure:
    """Create genome-wide CNA plot with mutation overlay.

    Args:
        seg_df: Parsed SEG dataframe for the sample.
        sample_result: Analysis result containing gene classifications.
        height: Plot height in pixels.

    Returns:
        A plotly Figure.
    """
    params = sample_result.params
    offsets = get_chrom_offsets()
    genome_size = sum(CHROM_SIZES.values())

    fig = go.Figure()

    # --- Chromosome background bands ---
    for i, (chrom, size) in enumerate(CHROM_SIZES.items()):
        if i % 2 == 1:
            fig.add_vrect(
                x0=offsets[chrom],
                x1=offsets[chrom] + size,
                fillcolor=COLORS["chrom_alt"],
                line_width=0,
                layer="below",
            )

    # --- Segments ---
    if not seg_df.is_empty():
        for row in seg_df.iter_rows(named=True):
            chrom = _normalize_chrom(row.get("chromosome", ""))
            if chrom not in offsets:
                continue
            start = offsets[chrom] + (row.get("start") or 0)
            end = offsets[chrom] + (row.get("end") or 0)
            seg_mean = row.get("seg_mean")
            if seg_mean is None:
                continue

            if seg_mean < params.del_threshold:
                color = COLORS["del"]
            elif seg_mean > params.gain_threshold:
                color = COLORS["amp"]
            else:
                color = COLORS["neutral"]

            fig.add_trace(
                go.Scatter(
                    x=[start, end],
                    y=[seg_mean, seg_mean],
                    mode="lines",
                    line={"color": color, "width": 3},
                    showlegend=False,
                    hoverinfo="text",
                    text=f"chr{chrom}:{row.get('start')}-{row.get('end')}<br>"
                    f"seg.mean: {seg_mean:.3f}",
                )
            )

    # --- Threshold lines ---
    for val, label, dash in [
        (params.deep_del_threshold, "Deep del", "dash"),
        (params.del_threshold, "LOH", "dot"),
        (params.gain_threshold, "Gain", "dot"),
        (params.amp_threshold, "Amp", "dash"),
    ]:
        fig.add_hline(
            y=val,
            line_dash=dash,
            line_color=COLORS["threshold"],
            line_width=1,
            annotation_text=label,
            annotation_position="right",
            annotation_font_size=9,
            annotation_font_color=COLORS["threshold"],
        )

    # --- Mutation overlay ---
    # Collect mutations for scatter plot
    mut_x, mut_y, mut_text, mut_color, mut_labels = [], [], [], [], []

    for gene_result in sample_result.gene_results:
        for mut in gene_result.mutations:
            chrom = _normalize_chrom(mut.chromosome)
            if chrom not in offsets:
                continue
            x = offsets[chrom] + mut.start_position
            y = mut.seg_mean if mut.seg_mean is not None else 0.0
            color = STATUS_COLORS.get(gene_result.two_hit_status, COLORS["single_hit"])

            vaf_str = f"{mut.vaf:.2f}" if mut.vaf is not None else "N/A"
            seg_str = f"{mut.seg_mean:.3f}" if mut.seg_mean is not None else "N/A"
            hover = (
                f"<b>{mut.hugo_symbol}</b> {mut.hgvsp_short}<br>"
                f"{mut.variant_classification}<br>"
                f"VAF: {vaf_str}<br>"
                f"seg.mean: {seg_str}<br>"
                f"Status: {gene_result.status_label}"
            )

            mut_x.append(x)
            mut_y.append(y)
            mut_text.append(hover)
            mut_color.append(color)

            # Label notable events
            if gene_result.is_biallelic or gene_result.two_hit_status in (
                TwoHitStatus.MUT_AMP,
                TwoHitStatus.MUT_GAIN,
            ):
                mut_labels.append((x, y, mut.hugo_symbol))

    if mut_x:
        fig.add_trace(
            go.Scatter(
                x=mut_x,
                y=mut_y,
                mode="markers",
                marker={
                    "symbol": "diamond",
                    "size": 9,
                    "color": mut_color,
                    "line": {"width": 1, "color": "white"},
                },
                hovertext=mut_text,
                hoverinfo="text",
                showlegend=False,
            )
        )

    # --- Gene labels for key events ---
    for x, y, label in mut_labels:
        fig.add_annotation(
            x=x,
            y=y,
            text=f"<b>{label}</b>",
            showarrow=True,
            arrowhead=0,
            arrowwidth=1,
            arrowcolor="#555",
            ax=0,
            ay=-25,
            font={"size": 10, "color": "#213448"},
            bgcolor="rgba(255,255,255,0.8)",
            borderpad=2,
        )

    # --- Layout ---
    # Chromosome tick labels at midpoints
    tickvals = []
    ticktext = []
    for chrom, size in CHROM_SIZES.items():
        mid = offsets[chrom] + size // 2
        tickvals.append(mid)
        ticktext.append(chrom)

    fig.update_layout(
        template="plotly_white",
        height=height,
        margin={"l": 50, "r": 30, "t": 40, "b": 50},
        title={
            "text": f"Genome-wide CNA — {sample_result.sample_id}",
            "font": {"size": 13},
        },
        xaxis={
            "range": [0, genome_size],
            "tickvals": tickvals,
            "ticktext": ticktext,
            "tickangle": 0,
            "tickfont": {"size": 9},
            "title": None,
            "showgrid": False,
        },
        yaxis={
            "range": [-2.5, 2.5],
            "title": {"text": "seg.mean (log₂ ratio)", "font": {"size": 11}},
            "zeroline": True,
            "zerolinecolor": "rgba(0,0,0,0.1)",
        },
        hovermode="closest",
    )

    return fig


# ---------------------------------------------------------------------------
# Export helpers
# ---------------------------------------------------------------------------


def fig_to_html_div(fig: go.Figure) -> str:
    """Export plotly figure as an HTML div (no full page wrapper)."""
    return fig.to_html(full_html=False, include_plotlyjs="cdn")


def fig_to_png_bytes(fig: go.Figure, width: int = 1200, height: int = 450) -> bytes:
    """Export plotly figure as PNG bytes (for PDF embedding)."""
    return fig.to_image(format="png", width=width, height=height, scale=2)
