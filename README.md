# two-hit

Integrated mutation + copy number analysis tool for identifying **biallelic inactivation** in tumor suppressors and **mutation-with-amplification** in oncogenes. Part of the [molpath.tools](https://molpath.tools) ecosystem.

## What it does

Given a **MAF** file (somatic mutations) and a **SEG** file (copy number segments), two-hit produces a per-gene report that classifies each event:

| Pattern | Example | Confidence |
|---------|---------|------------|
| Biallelic — compound het | 2 truncating TP53 mutations | HIGH |
| Biallelic — mutation + LOH | PTEN missense + hemizygous deletion (VAF > 0.7) | HIGH |
| Homozygous deletion | CDKN2A deep deletion (seg.mean < −0.5) | HIGH |
| Mutation + amplification | EGFR L858R + high-level amp | HIGH |
| Mutation + gain | KRAS G12D + copy gain | MEDIUM |
| Single hit | One mutation, no CNA support | LOW |

Gene roles (oncogene / TSG / both) are derived from the [OncoKB cancer gene list](https://www.oncokb.org/cancerGenes) (1,236 genes).

## Output

- **Per-gene summary table** with classification, confidence, VAF, and CNA status
- **Genome-wide CNA plot** (plotly) with mutations overlaid and key events labeled
- **PDF report** with embedded plot, methods, and parameter summary
- **JSON / TSV** for programmatic use

## Install

```bash
git clone https://github.com/viktorlj/two-hit.git
cd two-hit
uv sync
source .venv/bin/activate
```

## Usage

### CLI

```bash
# Analyze a MAF + SEG pair
two-hit report sample.maf sample.seg

# PDF output
two-hit report sample.maf sample.seg -f pdf -o report.pdf

# JSON output
two-hit report sample.maf sample.seg -f json

# Custom thresholds
two-hit report sample.maf sample.seg --del-threshold -0.4 --amp-threshold 0.8

# Analyze a sample from AACR GENIE (lazy-scans large files)
two-hit genie /path/to/genie/v18 GENIE-MSK-P-0001296-T01-IM3
```

### Web interface

```bash
two-hit serve
# Open http://localhost:8000
```

Upload MAF + SEG files or click **Run demo case** to try with a synthetic example.

### API

```bash
curl -X POST http://localhost:8000/api/analyze \
  -F "maf_file=@sample.maf" \
  -F "seg_file=@sample.seg"
```

## Input formats

**MAF** — Standard AACR GENIE MAF (tab-separated). Required columns: `Hugo_Symbol`, `Chromosome`, `Start_Position`, `End_Position`, `Variant_Classification`, `HGVSp_Short`, `Tumor_Sample_Barcode`. Optional: `t_alt_count`, `t_ref_count` (for VAF).

**SEG** — CBS-style segment file (tab-separated). Required columns: `ID`, `chrom`, `loc.start`, `loc.end`, `seg.mean`. Column matching is case-insensitive.

## Default thresholds

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--del-threshold` | −0.3 | Hemizygous deletion / LOH |
| `--deep-del-threshold` | −0.5 | Homozygous deletion |
| `--gain-threshold` | 0.3 | Copy number gain |
| `--amp-threshold` | 0.7 | High-level amplification |
| `--vaf-threshold` | 0.7 | High VAF supporting LOH |

## Tech stack

Python 3.12+ · [polars](https://pola.rs) · [pydantic](https://docs.pydantic.dev) v2 · [typer](https://typer.tiangolo.com) · [FastAPI](https://fastapi.tiangolo.com) + Jinja2 · [plotly](https://plotly.com/python/) · [fpdf2](https://py-pdf.github.io/fpdf2/) · [ruff](https://docs.astral.sh/ruff/)

## License

MIT
