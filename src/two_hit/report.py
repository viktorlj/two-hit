"""Output formatters: JSON, TSV, and human-readable text."""

from __future__ import annotations

from .models import GeneResult, SampleResult


def result_to_dict(result: SampleResult) -> dict:
    """Convert SampleResult to a JSON-serializable dict."""
    return result.model_dump(mode="json")


def result_to_json(result: SampleResult) -> str:
    """Render SampleResult as formatted JSON."""
    return result.model_dump_json(indent=2)


def gene_results_to_tsv(gene_results: list[GeneResult]) -> str:
    """Render gene results as a TSV table."""
    header = "Gene\tRole\t#Mutations\tBest Mutation\tCNA\tseg.mean\tVAF\tStatus\tConfidence"
    lines = [header]
    for g in gene_results:
        seg_str = f"{g.seg_mean:.3f}" if g.seg_mean is not None else "N/A"
        vaf_str = f"{g.vaf_max:.2f}" if g.vaf_max is not None else "N/A"
        lines.append(
            f"{g.hugo_symbol}\t{g.role.value}\t{g.mutation_count}\t{g.best_mutation}\t"
            f"{g.cna_status.value}\t{seg_str}\t{vaf_str}\t{g.two_hit_status.value}\t"
            f"{g.confidence.value}"
        )
    return "\n".join(lines) + "\n"


def format_human_readable(result: SampleResult) -> str:
    """Render SampleResult as human-readable terminal output."""
    lines = [
        "=" * 70,
        f"  Two-Hit Analysis Report — {result.sample_id}",
        "=" * 70,
        "",
        f"  Total mutations:     {result.total_mutations}",
        f"  Genes with events:   {result.total_genes_with_events}",
        f"  Biallelic events:    {result.biallelic_count}",
        f"  Mutation + amp:      {result.mut_amp_count}",
        "",
        "  Thresholds:",
        f"    Deletion (LOH):    {result.params.del_threshold}",
        f"    Deep deletion:     {result.params.deep_del_threshold}",
        f"    Gain:              {result.params.gain_threshold}",
        f"    Amplification:     {result.params.amp_threshold}",
        f"    High VAF:          {result.params.vaf_threshold}",
        "",
        "-" * 70,
    ]

    for g in result.gene_results:
        lines.append(
            f"  {g.hugo_symbol:<12s}  {g.role.value:<10s}  {g.status_label:<35s}  "
            f"[{g.confidence.value}]"
        )
        for m in g.mutations:
            m_vaf = f"{m.vaf:.2f}" if m.vaf is not None else "N/A"
            lines.append(
                f"    {m.variant_classification:<25s} {m.hgvsp_short:<15s} "
                f"CNA={m.cna_status.value:<8s} VAF={m_vaf}"
            )
        if g.notes:
            lines.append(f"    Notes: {'; '.join(g.notes)}")
        lines.append("")

    lines.append("=" * 70)
    return "\n".join(lines)
