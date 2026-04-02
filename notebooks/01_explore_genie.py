# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "marimo",
#     "polars>=1.0",
#     "plotly>=5.18",
#     "two-hit",
# ]
# ///

import marimo

__generated_with = "0.22.0"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell
def _(mo):
    mo.md(
        """
        # GENIE v18 — Two-Hit Exploration

        Explore the AACR GENIE v18 dataset using the **two-hit** analysis pipeline
        to identify biallelic inactivation events in tumor suppressors and
        mutation+amplification events in oncogenes.
        """
    )
    return


@app.cell
def _():
    from pathlib import Path

    import polars as pl

    from two_hit.analysis import analyze_sample
    from two_hit.gene_list import GeneRole, load_gene_roles
    from two_hit.io import filter_maf_for_sample, filter_seg_for_sample, scan_maf, scan_seg
    from two_hit.models import AnalysisParams
    from two_hit.plot import create_cna_plot
    from two_hit.report import format_human_readable

    PROJECT_ROOT = Path(__file__).resolve().parent.parent
    GENIE_DIR = PROJECT_ROOT / "data" / "genie_v18"
    gene_roles = load_gene_roles()
    params = AnalysisParams()
    return (
        GENIE_DIR,
        GeneRole,
        AnalysisParams,
        PROJECT_ROOT,
        Path,
        analyze_sample,
        filter_maf_for_sample,
        filter_seg_for_sample,
        format_human_readable,
        gene_roles,
        create_cna_plot,
        load_gene_roles,
        params,
        pl,
        scan_maf,
        scan_seg,
    )


@app.cell
def _(mo):
    mo.md("## 1. Lazy scan of MAF + SEG")
    return


@app.cell
def _(GENIE_DIR, mo, pl, scan_maf, scan_seg):
    maf_lf = scan_maf(GENIE_DIR / "data_mutations_extended.txt")
    seg_lf = scan_seg(GENIE_DIR / "data_cna_hg19.seg")

    maf_schema = maf_lf.collect_schema()
    seg_schema = seg_lf.collect_schema()

    mo.md(
        f"""
        **MAF schema** ({len(maf_schema)} columns): `{', '.join(maf_schema.names())}`

        **SEG schema** ({len(seg_schema)} columns): `{', '.join(seg_schema.names())}`
        """
    )
    return maf_lf, maf_schema, seg_lf, seg_schema


@app.cell
def _(mo):
    mo.md("## 2. Find interesting samples")
    return


@app.cell
def _(GENIE_DIR, mo, pl):
    _maf_full = pl.scan_csv(
        GENIE_DIR / "data_mutations_extended.txt",
        separator="\t",
        schema_overrides={"Chromosome": pl.Utf8},
    )

    truncating = ["Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site"]

    # TP53 with mutation + potential LOH (look for samples with TP53 mutations)
    tp53_samples = (
        _maf_full.filter(
            (pl.col("Hugo_Symbol") == "TP53")
            & pl.col("Variant_Classification").is_in(truncating)
        )
        .group_by("Tumor_Sample_Barcode")
        .agg(pl.len().alias("n_trunc"))
        .filter(pl.col("n_trunc") >= 2)
        .sort("n_trunc", descending=True)
        .head(5)
        .collect()
    )

    # BRCA1/2 with >=2 mutations
    brca_samples = (
        _maf_full.filter(
            pl.col("Hugo_Symbol").is_in(["BRCA1", "BRCA2"])
            & pl.col("Variant_Classification").is_in(truncating)
        )
        .group_by("Tumor_Sample_Barcode")
        .agg(
            pl.len().alias("n_trunc"),
            pl.col("Hugo_Symbol").unique().alias("genes"),
        )
        .filter(pl.col("n_trunc") >= 2)
        .sort("n_trunc", descending=True)
        .head(5)
        .collect()
    )

    # EGFR/ERBB2 mutations (potential mut+amp)
    egfr_samples = (
        _maf_full.filter(pl.col("Hugo_Symbol").is_in(["EGFR", "ERBB2"]))
        .group_by("Tumor_Sample_Barcode")
        .agg(pl.len().alias("n_muts"))
        .sort("n_muts", descending=True)
        .head(5)
        .collect()
    )

    mo.md(
        f"""
        ### TP53 compound het candidates (>=2 truncating)
        {mo.as_html(tp53_samples)}

        ### BRCA1/2 compound het candidates
        {mo.as_html(brca_samples)}

        ### EGFR/ERBB2 mutation candidates
        {mo.as_html(egfr_samples)}
        """
    )
    return brca_samples, egfr_samples, tp53_samples, truncating


@app.cell
def _(mo):
    mo.md("## 3. Analyze selected samples")
    return


@app.cell
def _(mo):
    sample_input = mo.ui.text(
        value="GENIE-MSK-P-0033859-T05-IH4",
        label="Sample barcode",
        full_width=True,
    )
    sample_input
    return (sample_input,)


@app.cell
def _(
    analyze_sample,
    filter_maf_for_sample,
    filter_seg_for_sample,
    gene_roles,
    maf_lf,
    mo,
    params,
    sample_input,
    seg_lf,
):
    _sid = sample_input.value.strip()
    _maf_df = filter_maf_for_sample(maf_lf, _sid)
    _seg_df = filter_seg_for_sample(seg_lf, _sid)

    mo.stop(
        _maf_df.is_empty() and _seg_df.is_empty(),
        mo.md(f"**No data found for sample `{_sid}`**"),
    )

    result = analyze_sample(_maf_df, _seg_df, sample_id=_sid, params=params, gene_roles=gene_roles)
    seg_df_for_plot = _seg_df

    mo.md(
        f"""
        ### {result.sample_id}
        - **Mutations:** {result.total_mutations}
        - **Genes:** {result.total_genes_with_events}
        - **Biallelic events:** {result.biallelic_count}
        - **Mut+Amp events:** {result.mut_amp_count}
        """
    )
    return result, seg_df_for_plot


@app.cell
def _(mo):
    mo.md("## 4. Gene results table")
    return


@app.cell
def _(mo, result):
    _rows = []
    for _g in result.gene_results:
        _muts_str = ", ".join(
            m.hgvsp_short for m in _g.mutations if m.hgvsp_short
        )
        _rows.append(
            {
                "Gene": _g.hugo_symbol,
                "Role": _g.role.value,
                "#Mut": _g.mutation_count,
                "Mutations": _muts_str,
                "CNA": _g.cna_status.value,
                "seg.mean": round(_g.seg_mean, 3) if _g.seg_mean is not None else None,
                "VAF": round(_g.vaf_max, 2) if _g.vaf_max is not None else None,
                "Status": _g.two_hit_status.value,
                "Confidence": _g.confidence.value,
            }
        )
    mo.ui.table(_rows, label="Gene results")
    return


@app.cell
def _(mo):
    mo.md("## 5. Genome-wide CNA plot")
    return


@app.cell
def _(create_cna_plot, mo, result, seg_df_for_plot):
    fig = create_cna_plot(seg_df_for_plot, result)
    mo.ui.plotly(fig)
    return (fig,)


@app.cell
def _(mo):
    mo.md("## 6. Human-readable report")
    return


@app.cell
def _(format_human_readable, mo, result):
    mo.md(f"```\n{format_human_readable(result)}\n```")
    return


@app.cell
def _(mo):
    mo.md(
        """
        ## 7. Batch validation

        Run the analysis on a sample of barcodes to check the distribution
        of two-hit statuses across the cohort.
        """
    )
    return


@app.cell
def _(
    analyze_sample,
    filter_maf_for_sample,
    filter_seg_for_sample,
    gene_roles,
    maf_lf,
    mo,
    params,
    pl,
    seg_lf,
):
    # Pick 20 random samples that have both MAF and SEG data
    _sample_ids = (
        maf_lf.select("Tumor_Sample_Barcode")
        .unique()
        .head(100)
        .collect()["Tumor_Sample_Barcode"]
        .sample(20, seed=42)
        .to_list()
    )

    _batch_rows = []
    for _sid in _sample_ids:
        _m = filter_maf_for_sample(maf_lf, _sid)
        _s = filter_seg_for_sample(seg_lf, _sid)
        if _m.is_empty():
            continue
        _r = analyze_sample(_m, _s, sample_id=_sid, params=params, gene_roles=gene_roles)
        _batch_rows.append(
            {
                "sample": _sid,
                "mutations": _r.total_mutations,
                "genes": _r.total_genes_with_events,
                "biallelic": _r.biallelic_count,
                "mut_amp": _r.mut_amp_count,
            }
        )

    batch_df = pl.DataFrame(_batch_rows)

    mo.md(
        f"""
        ### Batch results ({len(_batch_rows)} samples)

        {mo.as_html(batch_df)}

        **Summary:**
        - Samples with biallelic events: {batch_df.filter(pl.col("biallelic") > 0).height}
        - Samples with mut+amp events: {batch_df.filter(pl.col("mut_amp") > 0).height}
        - Mean mutations per sample: {batch_df["mutations"].mean():.1f}
        - Mean biallelic per sample: {batch_df["biallelic"].mean():.2f}
        """
    )
    return (batch_df,)


if __name__ == "__main__":
    app.run()
