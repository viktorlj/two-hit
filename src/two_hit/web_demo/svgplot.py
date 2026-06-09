"""Server-side inline SVG genome-wide CNA plot.

Replaces the plotly figure used by the full app. Produces a self-contained
<svg> string (no JS, no CDN) with copy-number segments and mutation markers
laid out across hg19 chromosomes.
"""

from __future__ import annotations

from html import escape

import polars as pl

from ..models import BIALLELIC_STATUSES, AnalysisParams, SampleResult, TwoHitStatus

# hg19 chromosome lengths (bp)
HG19_CHROM_LENGTHS: dict[str, int] = {
    "1": 249250621,
    "2": 243199373,
    "3": 198022430,
    "4": 191154276,
    "5": 180915260,
    "6": 171115067,
    "7": 159138663,
    "8": 146364022,
    "9": 141213431,
    "10": 135534747,
    "11": 135006516,
    "12": 133851895,
    "13": 115169878,
    "14": 107349540,
    "15": 102531392,
    "16": 90354753,
    "17": 81195210,
    "18": 78077248,
    "19": 59128983,
    "20": 63025520,
    "21": 48129895,
    "22": 51304566,
    "X": 155270560,
    "Y": 59373566,
}

_CHROM_ORDER = list(HG19_CHROM_LENGTHS)
_YMAX = 1.5  # clamp seg.mean to [-1.5, 1.5]

def _norm_chrom(c: str) -> str:
    c = str(c).strip().upper()
    return c[3:] if c.startswith("CHR") else c


def _seg_color(seg_mean: float | None, params: AnalysisParams) -> str:
    # Branches assume del_threshold < gain_threshold < amp_threshold (AnalysisParams defaults)
    if seg_mean is None:
        return "#9aa0a6"
    if seg_mean < params.del_threshold:
        return "#3a6ea5"  # deletion (blue)
    if seg_mean > params.amp_threshold:
        return "#c0392b"  # amplification (red)
    if seg_mean > params.gain_threshold:
        return "#e08e0b"  # gain (orange)
    return "#9aa0a6"  # neutral (grey)


def _mut_color(status: TwoHitStatus) -> str:
    if status in BIALLELIC_STATUSES:
        return "#6c3483"  # biallelic (purple)
    if status == TwoHitStatus.MUT_AMP:
        return "#c0392b"  # mut + amp (red)
    return "#444"


def genome_svg(
    seg_df: pl.DataFrame,
    result: SampleResult,
    width: int = 920,
    height: int = 320,
) -> str:
    """Render a genome-wide CNA plot with mutation overlay as an inline SVG string.

    Args:
        seg_df: Parsed SEG dataframe with columns chromosome, start, end, seg_mean.
        result: SampleResult from analyze_sample, containing gene_results and params.
        width: SVG canvas width in pixels.
        height: SVG canvas height in pixels.

    Returns:
        A self-contained SVG string starting with ``<svg`` and ending with ``</svg>``.
    """
    params = result.params
    margin_l, margin_r, margin_t, margin_b = 44, 12, 14, 26
    plot_w = width - margin_l - margin_r
    plot_h = height - margin_t - margin_b

    total = sum(HG19_CHROM_LENGTHS.values())
    offsets: dict[str, int] = {}
    acc = 0
    for c in _CHROM_ORDER:
        offsets[c] = acc
        acc += HG19_CHROM_LENGTHS[c]

    def gx(chrom: str, pos: int) -> float:
        return margin_l + (offsets[chrom] + pos) / total * plot_w

    def gy(val: float) -> float:
        v = max(-_YMAX, min(_YMAX, val))
        return margin_t + (1 - (v + _YMAX) / (2 * _YMAX)) * plot_h

    parts: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" '
        f'width="100%" role="img" aria-label="Genome-wide copy number plot" '
        f'style="font-family:var(--mp-mono,monospace);font-size:9px">'
    ]

    # zero baseline
    y0 = gy(0.0)
    parts.append(
        f'<line x1="{margin_l}" y1="{y0:.1f}" x2="{margin_l + plot_w}" y2="{y0:.1f}" '
        f'stroke="#d0d0d0" stroke-width="1"/>'
    )
    # y-axis ticks at -1, 0, +1
    for tick in (-1.0, 0.0, 1.0):
        yt = gy(tick)
        parts.append(
            f'<text x="{margin_l - 6}" y="{yt + 3:.1f}" text-anchor="end" fill="#888">'
            f"{tick:+.0f}</text>"
        )

    # chromosome boundaries + labels
    for c in _CHROM_ORDER:
        x = gx(c, 0)
        parts.append(
            f'<line x1="{x:.1f}" y1="{margin_t}" x2="{x:.1f}" y2="{margin_t + plot_h}" '
            f'stroke="#efefef" stroke-width="1"/>'
        )
        xc = gx(c, HG19_CHROM_LENGTHS[c] // 2)
        parts.append(
            f'<text x="{xc:.1f}" y="{height - 8}" text-anchor="middle" fill="#aaa">{c}</text>'
        )

    # segments as colored horizontal lines
    if not seg_df.is_empty():
        for row in seg_df.iter_rows(named=True):
            chrom = _norm_chrom(row.get("chromosome", ""))
            if chrom not in HG19_CHROM_LENGTHS:
                continue
            sm = row.get("seg_mean")
            if sm is None:
                continue
            x1 = gx(chrom, int(row["start"]))
            x2 = gx(chrom, int(row["end"]))
            yv = gy(float(sm))
            color = _seg_color(float(sm), params)
            parts.append(
                f'<line x1="{x1:.1f}" y1="{yv:.1f}" x2="{x2:.1f}" y2="{yv:.1f}" '
                f'stroke="{color}" stroke-width="2.5">'
                f"<title>chr{chrom} seg.mean={float(sm):.2f}</title></line>"
            )

    # mutation markers
    for g in result.gene_results:
        color = _mut_color(g.two_hit_status)
        for m in g.mutations:
            chrom = _norm_chrom(m.chromosome)
            if chrom not in HG19_CHROM_LENGTHS:
                continue
            x = gx(chrom, m.start_position)
            yv = gy(m.seg_mean if m.seg_mean is not None else 0.0)
            label = escape(
                f"{g.hugo_symbol} {m.protein_change or m.variant_classification}"
                f" [{g.status_label}]"
            )
            parts.append(
                f'<circle cx="{x:.1f}" cy="{yv:.1f}" r="3.5" fill="{color}" '
                f'stroke="#fff" stroke-width="0.8"><title>{label}</title></circle>'
            )

    parts.append("</svg>")
    return "".join(parts)
