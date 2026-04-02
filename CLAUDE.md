# two-hit

Integrated mutation + copy number analysis tool for identifying biallelic
inactivation events (tumor suppressors) and mutation-with-amplification events
(oncogenes). Part of the [molpath.tools](https://molpath.tools) ecosystem.

## What this tool does

Given a MAF file (somatic mutations) and a SEG file (copy number segments),
two-hit produces a per-gene report that integrates both data types to identify:

1. **TSG biallelic inactivation** — two hits on a tumor suppressor gene:
   - Two truncating/pathogenic mutations (compound heterozygous)
   - One mutation + LOH (deletion on the other allele, from SEG data)
   - Supported by VAF ≈ 1.0 in the context of a deletion
2. **Oncogene mutation + amplification** — activating mutation on an amplified
   allele, suggesting clonal advantage
3. **Single-hit events** — mutations or CNAs alone, for completeness

Gene classification (oncogene / TSG / both) comes from the OncoKB cancer gene
list, which is freely downloadable without a license.

## Tech stack

- **Python 3.12+**, managed with `uv`
- **polars** for all dataframe operations (NEVER pandas)
- **pydantic** v2 for data models and validation
- **typer** for the CLI
- **FastAPI** + Jinja2 for the web interface
- **weasyprint** or **fpdf2** for PDF generation (evaluate which works better)
- **plotly** for the genome-wide CNA plot with mutation overlay
- **ruff** for linting
- Marimo notebooks for development/exploration (NEVER Jupyter)

## Project structure

```
two-hit/
├── CLAUDE.md                 ← you are here
├── BOOTSTRAP_PROMPT.md       ← initial prompt for Claude Code
├── pyproject.toml            ← project metadata + dependencies (uv)
├── src/
│   └── two_hit/
│       ├── __init__.py
│       ├── models.py         ← pydantic models for MAF rows, SEG rows, gene reports
│       ├── io.py             ← MAF and SEG file parsers (polars-based)
│       ├── gene_list.py      ← OncoKB cancer gene list loader + caching
│       ├── analysis.py       ← core two-hit integration logic
│       ├── plot.py           ← genome-wide CNA plot with mutations overlaid
│       ├── report.py         ← PDF report generation
│       ├── cli.py            ← typer CLI entry point
│       └── web/
│           ├── app.py        ← FastAPI application
│           ├── static/
│           │   └── molpath.css  ← shared design system (sync from molpath-design)
│           └── templates/
│               ├── base.html
│               ├── upload.html
│               └── report.html  ← interactive HTML report (also used for PDF)
├── data/
│   └── oncokb_cancer_genes.tsv  ← cached OncoKB cancer gene list
├── tests/
│   └── ...
├── notebooks/
│   └── 01_explore_genie.py   ← Marimo notebook for dev/exploration with GENIE data
├── Dockerfile
└── .gitignore
```

## Input formats

### MAF file (mutations)
Standard AACR GENIE MAF format. Required columns:
- `Hugo_Symbol` — gene name
- `Chromosome`, `Start_Position`, `End_Position` — genomic coordinates
- `Variant_Classification` — e.g. Missense_Mutation, Nonsense_Mutation, Frame_Shift_Del, etc.
- `HGVSp_Short` — protein change (e.g. p.V600E)
- `t_alt_count`, `t_ref_count` — for VAF calculation
- `Tumor_Sample_Barcode` — sample identifier

Optional but useful: `FILTER`, `Variant_Type`, `Reference_Allele`, `Tumor_Seq_Allele2`

### SEG file (copy number)
Standard CBS-style `.seg` format. Required columns:
- `ID` — sample identifier (must match Tumor_Sample_Barcode in MAF)
- `chrom` — chromosome
- `loc.start`, `loc.end` — segment boundaries
- `num.mark` — number of probes/markers
- `seg.mean` — log2 ratio segment mean

Column name matching should be case-insensitive and tolerant of common
variations (e.g., `Chromosome` vs `chrom`, `Start` vs `loc.start`).

## Analysis logic

### CNA thresholds (configurable via CLI flags and web UI)
- **Deep deletion**: seg.mean < -0.5 (homozygous deletion)
- **Hemizygous deletion / LOH**: seg.mean < -0.3
- **Gain**: seg.mean > 0.3
- **Amplification**: seg.mean > 0.7

### Mapping mutations to segments
For each mutation, find the overlapping segment by chromosome and position.
A mutation "overlaps" a segment if its Start_Position falls within
[loc.start, loc.end] on the same chromosome.

### VAF calculation
VAF = t_alt_count / (t_alt_count + t_ref_count)
High VAF (> 0.7) in the context of a hemizygous deletion strongly supports
LOH / biallelic inactivation.

### Two-hit classification logic

For each gene with at least one mutation:

1. **Lookup gene role** in OncoKB (oncogene / TSG / both / unknown)
2. **Count mutations** in the gene for this sample
3. **Find overlapping CNA segment** for each mutation
4. **Classify**:

   **TSG two-hit patterns:**
   - `BIALLELIC_COMPOUND_HET`: ≥2 truncating/pathogenic mutations
   - `BIALLELIC_MUT_LOH`: 1 mutation + overlapping deletion (seg.mean < -0.3)
   - `BIALLELIC_MUT_LOH_HIGH_VAF`: same as above but VAF > 0.7 (high confidence)
   - `HOMDEL`: deep deletion (seg.mean < -0.5), no mutation needed
   - `SINGLE_HIT_MUT`: 1 mutation only, no supporting CNA
   - `SINGLE_HIT_DEL`: deletion only, no mutation

   **Oncogene patterns:**
   - `MUT_AMP`: mutation + overlapping amplification (seg.mean > 0.7)
   - `MUT_GAIN`: mutation + overlapping gain (seg.mean > 0.3)
   - `MUT_ONLY`: mutation without CNA
   - `AMP_ONLY`: amplification without mutation

   **Both (e.g. TP53 classified as both):**
   - Apply both TSG and oncogene logic, report the most informative pattern

### Truncating mutation types
These Variant_Classification values count as truncating for compound het logic:
`Nonsense_Mutation`, `Frame_Shift_Del`, `Frame_Shift_Ins`, `Splice_Site`,
`Splice_Region` (debatable — make configurable), `Translation_Start_Site`

## Output

### Per-gene summary table
| Gene | Role | Mutations | CNA | VAF | Two-Hit Status | Confidence |
Each row is one gene with at least one event. Confidence is
HIGH / MEDIUM / LOW based on the evidence pattern.

### Genome-wide CNA plot
Plotly figure showing:
- X-axis: genomic position (chromosomes laid out sequentially)
- Y-axis: seg.mean (log2 ratio)
- Segments colored by CNA status (blue=deletion, red=amplification, grey=neutral)
- Mutations overlaid as points/markers at their genomic position
- Color-coded by two-hit status
- Gene labels for key events

### PDF report
Contains: sample ID, date, parameter summary, per-gene table, CNA plot,
and a methods section describing the algorithm.

## CLI interface

```bash
# Basic usage
two-hit report --maf sample.maf --seg sample.seg --output report.pdf

# With custom thresholds
two-hit report --maf sample.maf --seg sample.seg \
    --del-threshold -0.3 \
    --amp-threshold 0.7 \
    --vaf-threshold 0.7 \
    --output report.pdf

# HTML output instead of PDF
two-hit report --maf sample.maf --seg sample.seg --format html --output report.html

# Process GENIE data (for testing/development)
two-hit genie --genie-dir /path/to/genie/v18 --sample GENIE-MSK-P-0000001
```

## Web interface

FastAPI app with:
1. **Upload page** — drag-and-drop MAF + SEG files, configure thresholds
2. **Interactive report page** — same content as PDF but with:
   - Sortable/filterable gene table
   - Interactive Plotly CNA plot (hover for details)
   - "Download PDF" button that generates and serves the PDF
3. **API endpoint** — POST /api/analyze accepting MAF + SEG files, returns JSON

## Design system

Follow the molpath.tools design conventions:
- Use `molpath.css` custom properties (`--mp-*` tokens)
- Typography: Space Grotesk for body, IBM Plex Mono for code/labels
- Accent color: `--mp-accent` (#1f6f61 teal)
- Panel background: `--mp-panel` (#fffef9 warm white)
- No Tailwind, no CSS-in-JS, no frameworks
- Consistent with other tools: alascca-classify, ec-molsubtype, tumor-likelihood

## Development workflow

1. Use Marimo notebooks in `notebooks/` for exploration and algorithm dev
2. Once logic is solid, refactor into `src/two_hit/` modules
3. CLI wraps the core logic for batch use
4. Web interface wraps the same logic for interactive use
5. Test with GENIE v18 data (symlinked in `data/raw/`)

## Commands

```bash
uv sync                        # install dependencies
uv run two-hit --help          # CLI help
uv run two-hit report ...      # generate report
uv run fastapi dev src/two_hit/web/app.py  # dev server
uv run ruff check src/         # lint
uv run ruff format src/        # format
uv run marimo edit notebooks/01_explore_genie.py  # open notebook
```

## Deployment

Same pattern as other molpath.tools apps:
- Dockerfile → Google Cloud Run
- Custom domain: `two-hit.molpath.tools`
- Cloudflare DNS: CNAME → ghs.googlehosted.com

## GENIE v18 test data

For development and testing, use AACR GENIE v18 data:
- MAF: `data_mutations_extended.txt`
- SEG: `data_cna_hg19.seg` or `genie_data_cna_hg19.seg`
- Clinical: `data_clinical_sample.txt` (for sample metadata)

Pick a few well-characterized samples to test with:
- A sample with known TP53 biallelic inactivation (mutation + LOH)
- A sample with KRAS mutation + no CNA (expected for oncogene)
- A sample with EGFR amplification + mutation
- A sample with BRCA1/2 compound heterozygous mutations

## Important conventions

- **NEVER use pandas** — always polars
- **NEVER use Jupyter** — always Marimo
- **Use `uv`** for all package management
- **Use pydantic v2** for all data models
- **Use ruff** for linting and formatting
- All code in `src/two_hit/` with proper `__init__.py` exports
- Type hints on all functions
- Docstrings on all public functions
