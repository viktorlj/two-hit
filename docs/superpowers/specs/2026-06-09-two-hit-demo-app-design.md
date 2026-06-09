# two-hit lightweight demo app — design

**Date:** 2026-06-09
**Status:** Approved (pre-implementation)

## Problem

The full `two-hit` web app (`src/two_hit/web/`) is too resource-intensive to run
reliably on Railway. The measured request-time memory spike is dominated by
**plotly figure generation (+84 MB per request**, 54 → 138 MB even for a tiny
sample), with secondary cost from `fpdf2` (+20 MB). Prior attempts to survive on
Railway (lazy imports, a 30-min idle watchdog that exits the process for a clean
restart, `railway.toml` restart-retry policy) made the deployment fragile rather
than light.

## Goal

Keep the full app intact for local use, and add a **separate, lightweight demo
web app** that runs comfortably on Railway by shedding the heavy dependencies and
bounding its inputs.

## Decisions (locked with user)

1. **Separate lightweight app** — new `src/two_hit/web_demo/`, reusing the core
   analysis modules, importing **no plotly and no fpdf2**.
2. **Curated cases + bounded uploads** — 3–4 preset example cases plus
   panel-restricted, size-capped uploads.
3. **Server-side inline SVG** for the genome-wide CNA plot (no plotly.js, no
   client JS).
4. **~25-gene curated panel**, hardcoded (no OncoKB download for the demo).
5. **Replace the existing Railway service** with the demo; drop the
   idle-watchdog / restart hacks.

## Footprint budget

Measured import cost (this machine):

| Stack | RSS delta |
|-------|-----------|
| FastAPI + Jinja + Uvicorn | ~34 MB |
| + polars | ~17 MB |
| + plotly (figure gen) | **+84 MB at request time** (dropped) |
| + fpdf2 | +20 MB (dropped) |

Demo target: **~50–55 MB idle, ~65 MB per request**. No watchdog needed.

## Architecture

New package `src/two_hit/web_demo/`:

```
src/two_hit/web_demo/
├── __init__.py
├── app.py          # minimal FastAPI app factory; NO idle watchdog, NO lifespan hacks
├── routes.py       # GET / (landing), POST /analyze (upload), GET /case/{name}
├── panel.py        # ~25-gene curated panel -> {Hugo_Symbol: GeneRole}
├── svgplot.py      # server-side SVG genome plot generator
├── cases/          # preset MAF + SEG files for the demo archetypes
│   ├── tp53_loh.maf / .seg
│   ├── kras_mut_only.maf / .seg
│   ├── egfr_amp_mut.maf / .seg
│   └── brca_compound_het.maf / .seg
├── static/
│   ├── molpath.css   # synced from molpath-design
│   ├── demo.css      # small demo-specific styles
│   └── favicon.svg
└── templates/
    ├── base.html
    ├── index.html    # landing: case buttons + bounded upload form
    └── report.html   # results table + inline SVG plot
```

**Reused unchanged from core:** `analysis.analyze_sample`, `io.parse_maf`,
`io.parse_seg`, `models.*`. These depend on polars (acceptable, +17 MB). The
full app under `src/two_hit/web/` is **not modified**.

## Components

### `panel.py`
- Exposes `PANEL: dict[str, GeneRole]` for ~25 core genes:
  `TP53, RB1, PTEN, APC, BRCA1, BRCA2, CDKN2A, NF1, VHL, STK11, SMAD4, KEAP1,
  ATM, MLH1` (TSG); `KRAS, NRAS, BRAF, EGFR, PIK3CA, ERBB2, MET, MYC, CCND1,
  MDM2, CTNNB1` (oncogene/both as appropriate per OncoKB).
- Roles hardcoded so the demo needs no OncoKB file or network at runtime.
- Helper `filter_maf_to_panel(maf_df) -> maf_df` drops rows whose `Hugo_Symbol`
  is not in `PANEL`, before analysis.

### `routes.py`
- `GET /` → `index.html` with the list of preset cases.
- `GET /case/{name}` → load preset MAF+SEG from `cases/`, run analysis, render
  `report.html`. Unknown `name` → 404 back to index with an error message.
- `POST /analyze` → accept `maf_file` + `seg_file` uploads:
  - Enforce **≤1 MB per file** (reject with a clear message otherwise).
  - Decode UTF-8, parse, **filter MAF to panel**, run `analyze_sample` with
    `gene_roles=PANEL` and default `AnalysisParams()`.
  - Single sample only (the parser/analysis already operate per-sample; if
    multiple `Tumor_Sample_Barcode` values are present, use the first and note
    it).
  - Render `report.html`.
- No result store, no `/api`, no PDF route.

### `svgplot.py`
- `genome_svg(seg_df, result) -> str` returns a self-contained `<svg>` element.
- Genome-wide X layout using **fixed hg19 chromosome lengths/offsets** (embedded
  constant dict, chr1–22 + X, Y).
- Y axis = `seg.mean` (log2 ratio), clamped to a sensible range.
- Segments drawn as colored horizontal lines: blue = deletion
  (`seg.mean < del_threshold`), red = amplification (`> amp_threshold`),
  grey = neutral.
- Mutations overlaid as markers at genomic position, colored by two-hit status.
- Hover detail via SVG `<title>` children (no JS).

### Templates
- `report.html` receives the **`SampleResult` object directly** (not the JSON
  dict), so Jinja can use computed properties (`status_label`, `mutation_count`,
  `best_mutation`, `is_biallelic`). Renders: sample id, summary counts, per-gene
  table, and the embedded SVG plot.
- Shares `molpath.css` design tokens + a small `demo.css`.

### `app.py`
- `create_app()` mounts `/static`, configures Jinja2, includes the router.
- No `lifespan`, no idle watchdog, no `IDLE_TIMEOUT`, no `os._exit`.

## Data flow

```
preset case OR upload
   │
   ├─ parse_maf / parse_seg            (polars)
   ├─ filter_maf_to_panel              (panel.py)
   ├─ analyze_sample(gene_roles=PANEL) (core analysis, unchanged)
   │       └─ SampleResult
   ├─ genome_svg(seg_df, result)       (svgplot.py -> <svg> string)
   └─ render report.html(result, svg)
```

## Error handling

- File too large (>1 MB) → re-render `index.html` with an error banner, HTTP 400.
- Parse failure / malformed input → error banner, HTTP 400, exception logged.
- Empty result after panel filtering → render report with an explanatory note
  ("no panel genes mutated in this sample").
- Unknown preset case name → 404 back to index.

## Testing

- `tests/test_panel.py` — `PANEL` contains expected genes with correct roles;
  `filter_maf_to_panel` drops off-panel rows and keeps panel rows.
- `tests/test_web_demo.py` — FastAPI `TestClient`:
  - `GET /` returns 200 and lists cases.
  - each `GET /case/{name}` returns 200 and contains expected gene(s)
    (e.g. TP53 case report mentions biallelic TP53).
  - `POST /analyze` with an oversized file returns 400.
  - `POST /analyze` with a valid small MAF+SEG returns 200 and a report.
- `tests/test_svgplot.py` — `genome_svg` returns a string starting with `<svg`
  and containing segment + mutation elements for a known input.

## Deployment

- **Dockerfile:** change `CMD` to launch `two_hit.web_demo.app:app`. Continue to
  COPY `src/`; `data/` (OncoKB file) no longer required by the demo but harmless
  to keep. `cases/` ship inside the package.
- **Remove** `railway.toml` restart-retry hack and the `IDLE_TIMEOUT` env
  reliance (demo has no watchdog).
- Verify locally (`uvicorn two_hit.web_demo.app:app`), push to GitHub, then
  confirm a successful deploy and low memory via the Railway MCP.
- Domain `two-hit.molpath.tools` stays pointed at the same (now demo) service.

## Out of scope (YAGNI)

PDF export, JSON API, in-memory result store, full OncoKB (1237-gene) list,
plotly, multi-sample batches, UI threshold configuration. Users needing these
run the full app locally (`uvicorn two_hit.web.app:app`).

## Notes / follow-ups

- README should gain a short "Demo vs. full app" section explaining the split
  and how to run each locally.
- The full app's `web/` Railway-survival hacks (idle watchdog) become dead code
  for deployment but stay for anyone who still wants to host the full app; not
  removed in this change.
