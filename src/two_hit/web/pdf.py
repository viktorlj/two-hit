"""PDF report generation using fpdf2."""

from __future__ import annotations

import io
import tempfile
from datetime import datetime

from fpdf import FPDF


def _sanitize(text: str) -> str:
    """Replace Unicode characters that built-in PDF fonts can't render."""
    return (
        text.replace("\u2014", "--")
        .replace("\u2013", "-")
        .replace("\u2018", "'")
        .replace("\u2019", "'")
        .replace("\u201c", '"')
        .replace("\u201d", '"')
        .replace("\u2265", ">=")
        .replace("\u2264", "<=")
        .replace("\u2192", "->")
        .replace("\u2082", "2")
        .replace("\u2080", "0")
    )


STATUS_COLORS: dict[str, tuple[int, int, int]] = {
    "biallelic_compound_het": (198, 40, 40),
    "biallelic_mut_loh": (198, 40, 40),
    "biallelic_mut_loh_high_vaf": (198, 40, 40),
    "homdel": (198, 40, 40),
    "mut_amp": (21, 101, 192),
    "mut_gain": (21, 101, 192),
    "single_hit_mut": (97, 97, 97),
    "single_hit_del": (97, 97, 97),
    "mut_only": (97, 97, 97),
    "amp_only": (21, 101, 192),
    "unclassified": (150, 150, 150),
}

STATUS_LABELS: dict[str, str] = {
    "biallelic_compound_het": "Biallelic (compound het)",
    "biallelic_mut_loh": "Biallelic (mut + LOH)",
    "biallelic_mut_loh_high_vaf": "Biallelic (mut + LOH, high VAF)",
    "homdel": "Homozygous deletion",
    "mut_amp": "Mutation + amplification",
    "mut_gain": "Mutation + gain",
    "single_hit_mut": "Single hit (mutation)",
    "single_hit_del": "Single hit (deletion)",
    "mut_only": "Mutation only",
    "amp_only": "Amplification only",
    "unclassified": "Unclassified",
}


def generate_pdf(result_data: dict, plot_png: bytes | None = None) -> bytes:
    """Generate a PDF two-hit analysis report."""
    pdf = FPDF()
    pdf.set_auto_page_break(auto=True, margin=20)
    pdf.add_page()
    pw = pdf.w - pdf.l_margin - pdf.r_margin

    sample_id = result_data.get("sample_id", "unknown")
    version = result_data.get("version", "0.1.0")

    # --- Header ---
    pdf.set_font("Helvetica", "B", 16)
    pdf.cell(pw, 10, "Two-Hit Biallelic Inactivation Report", new_x="LMARGIN", new_y="NEXT")
    pdf.set_font("Helvetica", "", 9)
    pdf.set_text_color(100, 100, 100)
    pdf.cell(
        pw,
        5,
        f"Generated {datetime.now().strftime('%Y-%m-%d %H:%M')}  |  two-hit v{version}",
        new_x="LMARGIN",
        new_y="NEXT",
    )
    pdf.ln(4)

    # --- Sample info ---
    pdf.set_font("Helvetica", "", 11)
    pdf.set_text_color(0, 0, 0)
    pdf.cell(pw, 7, f"Sample: {sample_id}", new_x="LMARGIN", new_y="NEXT")

    biallelic = result_data.get("biallelic_count", 0)
    mut_amp = result_data.get("mut_amp_count", 0)
    total_muts = result_data.get("total_mutations", 0)
    total_genes = result_data.get("total_genes_with_events", 0)

    pdf.set_font("Helvetica", "", 10)
    pdf.cell(
        pw,
        6,
        f"Mutations: {total_muts}  |  Genes: {total_genes}  |  "
        f"Biallelic: {biallelic}  |  Mut+Amp: {mut_amp}",
        new_x="LMARGIN",
        new_y="NEXT",
    )
    pdf.ln(3)

    # --- Thresholds ---
    params = result_data.get("params", {})
    _section_header(pdf, pw, "Analysis Parameters")
    pdf.set_font("Courier", "", 8)
    pdf.set_text_color(60, 60, 60)
    for key in (
        "del_threshold",
        "deep_del_threshold",
        "gain_threshold",
        "amp_threshold",
        "vaf_threshold",
    ):
        val = params.get(key, "N/A")
        label = key.replace("_", " ").title()
        pdf.cell(pw, 4, f"  {label}: {val}", new_x="LMARGIN", new_y="NEXT")
    splice = "Yes" if params.get("include_splice_region") else "No"
    pdf.cell(pw, 4, f"  Include Splice Region: {splice}", new_x="LMARGIN", new_y="NEXT")
    pdf.set_text_color(0, 0, 0)
    pdf.ln(3)

    # --- Gene results table ---
    _section_header(pdf, pw, "Gene Results")

    # Column widths in mm, positioned using absolute x coordinates
    # to guarantee alignment regardless of fpdf2 cursor behavior
    col_defs = [
        # (width, header, align)
        (18, "Gene", "L"),
        (20, "Role", "L"),
        (28, "Mutation", "L"),
        (16, "seg.mean", "R"),
        (50, "Status", "L"),
        (18, "Conf.", "L"),
    ]
    row_h = 6
    font_sz = 8
    x_start = pdf.l_margin

    # Compute absolute x positions
    col_x = []
    x = x_start
    for w, _, _ in col_defs:
        col_x.append(x)
        x += w

    # Header row
    pdf.set_font("Helvetica", "B", font_sz)
    pdf.set_fill_color(235, 235, 235)
    y = pdf.get_y()
    for i, (w, header, _) in enumerate(col_defs):
        pdf.set_xy(col_x[i], y)
        pdf.cell(w, row_h, header, border="B", fill=True)
    pdf.set_y(y + row_h)

    # Table rows
    gene_results = result_data.get("gene_results", [])
    for g in gene_results:
        status = g.get("two_hit_status", "unclassified")
        r, gb, b = STATUS_COLORS.get(status, (0, 0, 0))

        seg_str = f"{g['seg_mean']:.2f}" if g.get("seg_mean") is not None else "N/A"
        best_mut = ""
        muts = g.get("mutations", [])
        for m in muts:
            if m.get("hgvsp_short"):
                best_mut = _sanitize(m["hgvsp_short"])
                break
        if not best_mut and muts:
            best_mut = muts[0].get("variant_classification", "")

        status_label = STATUS_LABELS.get(status, status)
        conf = g.get("confidence", "")
        role = g.get("role", "")

        row_data = [
            g.get("hugo_symbol", ""),
            role,
            best_mut[:20],
            seg_str,
            status_label,
            conf,
        ]

        y = pdf.get_y()
        # Draw bottom border line across full table width
        table_w = sum(w for w, _, _ in col_defs)
        pdf.set_draw_color(220, 220, 220)
        pdf.line(x_start, y + row_h, x_start + table_w, y + row_h)

        for i, (w, _, align) in enumerate(col_defs):
            pdf.set_xy(col_x[i], y)

            # Gene name bold, status colored
            if i == 0:
                pdf.set_font("Helvetica", "B", font_sz)
                pdf.set_text_color(0, 0, 0)
            elif i == 4:  # Status
                pdf.set_font("Helvetica", "B", font_sz)
                pdf.set_text_color(r, gb, b)
            else:
                pdf.set_font("Helvetica", "", font_sz)
                pdf.set_text_color(0, 0, 0)

            pdf.cell(w, row_h, row_data[i], align=align)

        pdf.set_text_color(0, 0, 0)
        pdf.set_y(y + row_h)

    pdf.ln(4)

    # --- CNA plot ---
    if plot_png:
        _section_header(pdf, pw, "Genome-wide CNA Plot")
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
            tmp.write(plot_png)
            tmp.flush()
            pdf.image(tmp.name, x=pdf.l_margin, w=pw)
        pdf.ln(4)

    # --- Methods ---
    _section_header(pdf, pw, "Methods")
    pdf.set_font("Helvetica", "", 8)
    pdf.set_text_color(60, 60, 60)
    methods_text = (
        "This report integrates somatic mutations (MAF) with copy number segments (SEG) to "
        "identify biallelic inactivation in tumor suppressor genes and mutation+amplification "
        "in oncogenes. Gene roles are derived from the OncoKB cancer gene list. "
        "For tumor suppressors, biallelic events are identified as: compound heterozygous "
        "truncating mutations, single mutation with loss of heterozygosity (LOH from copy "
        "number deletion), or homozygous deletion. For oncogenes, mutation co-occurring with "
        "amplification or gain is flagged. Confidence is HIGH for biallelic events with "
        "strong evidence (compound het, high VAF + LOH), MEDIUM for mutation + deletion, "
        "and LOW for single-hit events."
    )
    pdf.multi_cell(pw, 4, _sanitize(methods_text))
    pdf.ln(3)

    # --- Disclaimer ---
    pdf.set_font("Helvetica", "I", 8)
    pdf.set_text_color(120, 120, 120)
    pdf.multi_cell(
        pw,
        4,
        "DISCLAIMER: This report is for research and clinical decision support only. "
        "Biallelic inactivation calls should be confirmed by orthogonal methods. "
        "Copy number thresholds are configurable and should be validated for the "
        "specific sequencing assay used.",
    )

    buf = io.BytesIO()
    pdf.output(buf)
    return buf.getvalue()


def _section_header(pdf: FPDF, pw: float, title: str) -> None:
    pdf.set_font("Helvetica", "B", 11)
    pdf.set_text_color(0, 0, 0)
    pdf.cell(pw, 7, title, new_x="LMARGIN", new_y="NEXT")
    pdf.set_draw_color(200, 200, 200)
    pdf.line(pdf.l_margin, pdf.get_y(), pdf.l_margin + pw, pdf.get_y())
    pdf.ln(2)
