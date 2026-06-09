# two-hit Lightweight Demo App Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a separate, dependency-light `two-hit` demo web app (`src/two_hit/web_demo/`) that runs comfortably on Railway by dropping plotly + fpdf2, restricting analysis to a ~25-gene panel, and bounding uploads — while the full app under `src/two_hit/web/` stays unchanged for local use.

**Architecture:** New FastAPI package reuses the core science (`io.parse_maf/parse_seg`, `analysis.analyze_sample`, `models`) but renders the genome plot as a server-side inline SVG (`svgplot.py`) instead of plotly, and serves a small set of preset cases plus size-capped uploads filtered to a hardcoded gene panel (`panel.py`). No result store, no PDF, no JSON API, no idle watchdog.

**Tech Stack:** Python 3.12, FastAPI + Uvicorn, Jinja2, polars (core analysis only), pydantic v2. Tests via pytest + `fastapi.testclient`. No plotly, no fpdf2 imported anywhere in `web_demo`.

**Spec:** `docs/superpowers/specs/2026-06-09-two-hit-demo-app-design.md`

---

## File Structure

```
src/two_hit/web_demo/
├── __init__.py        # exports create_app
├── panel.py           # PANEL dict + filter_maf_to_panel()
├── svgplot.py         # genome_svg() server-side SVG generator
├── app.py             # FastAPI app factory; no watchdog
├── routes.py          # GET / , GET /case/{name}, POST /analyze
├── cases/             # preset MAF+SEG files
│   ├── tp53_loh.maf / tp53_loh.seg
│   ├── kras_mut_only.maf / kras_mut_only.seg
│   ├── egfr_amp_mut.maf / egfr_amp_mut.seg
│   └── brca_compound_het.maf / brca_compound_het.seg
├── static/
│   ├── molpath.css    # copied from web/static/molpath.css
│   ├── demo.css        # small demo styles
│   └── favicon.svg     # copied from web/static/favicon.svg
└── templates/
    ├── base.html
    ├── index.html
    └── report.html

tests/
├── test_panel.py
├── test_svgplot.py
└── test_web_demo.py
```

**Template note:** `report.html` receives the live `SampleResult` object (not a dict), so use `.value` on enums (`g.role.value`, `g.cna_status.value`, `g.two_hit_status.value`, `g.confidence.value`) and the `g.status_label` property for the human label.

---

## Task 1: Gene panel module

**Files:**
- Create: `src/two_hit/web_demo/__init__.py`
- Create: `src/two_hit/web_demo/panel.py`
- Test: `tests/test_panel.py`

- [ ] **Step 1: Write the failing test**

Create `tests/test_panel.py`:

```python
import polars as pl

from two_hit.gene_list import GeneRole
from two_hit.web_demo.panel import PANEL, filter_maf_to_panel


def test_panel_has_25_genes_with_expected_roles():
    assert len(PANEL) == 25
    assert PANEL["TP53"] == GeneRole.TSG
    assert PANEL["BRCA1"] == GeneRole.TSG
    assert PANEL["KRAS"] == GeneRole.ONCOGENE
    assert PANEL["EGFR"] == GeneRole.ONCOGENE
    # 14 TSG + 11 oncogene
    assert sum(1 for r in PANEL.values() if r == GeneRole.TSG) == 14
    assert sum(1 for r in PANEL.values() if r == GeneRole.ONCOGENE) == 11


def test_filter_maf_to_panel_drops_off_panel_genes():
    df = pl.DataFrame(
        {
            "Hugo_Symbol": ["TP53", "NOTCH1", "KRAS"],
            "Chromosome": ["17", "9", "12"],
            "Start_Position": [7577120, 139399365, 25398284],
        }
    )
    out = filter_maf_to_panel(df)
    assert out["Hugo_Symbol"].to_list() == ["TP53", "KRAS"]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_panel.py -v`
Expected: FAIL — `ModuleNotFoundError: two_hit.web_demo.panel`.

- [ ] **Step 3: Write minimal implementation**

Create `src/two_hit/web_demo/__init__.py`:

```python
"""Lightweight two-hit demo web app (Railway-friendly).

Drops plotly and fpdf2; restricts analysis to a curated gene panel and
bounded uploads. The full app lives in two_hit.web.
"""

from __future__ import annotations


def create_app():
    """Lazily build the demo FastAPI app (avoids import cost at package import)."""
    from .app import create_app as _create

    return _create()
```

Create `src/two_hit/web_demo/panel.py`:

```python
"""Curated ~25-gene panel for the two-hit demo.

Roles match the OncoKB cancer gene list. Hardcoded so the demo needs no
OncoKB file or network access at runtime.
"""

from __future__ import annotations

import polars as pl

from ..gene_list import GeneRole

PANEL: dict[str, GeneRole] = {
    # Tumor suppressors (14)
    "TP53": GeneRole.TSG,
    "RB1": GeneRole.TSG,
    "PTEN": GeneRole.TSG,
    "APC": GeneRole.TSG,
    "BRCA1": GeneRole.TSG,
    "BRCA2": GeneRole.TSG,
    "CDKN2A": GeneRole.TSG,
    "NF1": GeneRole.TSG,
    "VHL": GeneRole.TSG,
    "STK11": GeneRole.TSG,
    "SMAD4": GeneRole.TSG,
    "KEAP1": GeneRole.TSG,
    "ATM": GeneRole.TSG,
    "MLH1": GeneRole.TSG,
    # Oncogenes (11)
    "KRAS": GeneRole.ONCOGENE,
    "NRAS": GeneRole.ONCOGENE,
    "BRAF": GeneRole.ONCOGENE,
    "EGFR": GeneRole.ONCOGENE,
    "PIK3CA": GeneRole.ONCOGENE,
    "ERBB2": GeneRole.ONCOGENE,
    "MET": GeneRole.ONCOGENE,
    "MYC": GeneRole.ONCOGENE,
    "CCND1": GeneRole.ONCOGENE,
    "MDM2": GeneRole.ONCOGENE,
    "CTNNB1": GeneRole.ONCOGENE,
}


def filter_maf_to_panel(maf_df: pl.DataFrame) -> pl.DataFrame:
    """Keep only rows whose Hugo_Symbol is in the demo panel."""
    if maf_df.is_empty():
        return maf_df
    return maf_df.filter(pl.col("Hugo_Symbol").is_in(list(PANEL)))
```

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run pytest tests/test_panel.py -v`
Expected: PASS (2 tests).

- [ ] **Step 5: Commit**

```bash
git add src/two_hit/web_demo/__init__.py src/two_hit/web_demo/panel.py tests/test_panel.py
git commit -m "feat: add demo gene panel module"
```

---

## Task 2: Preset case files

**Files:**
- Create: `src/two_hit/web_demo/cases/tp53_loh.maf`
- Create: `src/two_hit/web_demo/cases/tp53_loh.seg`
- Create: `src/two_hit/web_demo/cases/kras_mut_only.maf`
- Create: `src/two_hit/web_demo/cases/kras_mut_only.seg`
- Create: `src/two_hit/web_demo/cases/egfr_amp_mut.maf`
- Create: `src/two_hit/web_demo/cases/egfr_amp_mut.seg`
- Create: `src/two_hit/web_demo/cases/brca_compound_het.maf`
- Create: `src/two_hit/web_demo/cases/brca_compound_het.seg`

> All files are TAB-separated. Preserve the tabs exactly.

- [ ] **Step 1: Create the four MAF files**

`tp53_loh.maf`:
```
Hugo_Symbol	Chromosome	Start_Position	End_Position	Variant_Classification	HGVSp_Short	Tumor_Sample_Barcode	t_alt_count	t_ref_count
TP53	17	7577120	7577120	Nonsense_Mutation	p.R213*	CASE-TP53-LOH	90	10
```

`kras_mut_only.maf`:
```
Hugo_Symbol	Chromosome	Start_Position	End_Position	Variant_Classification	HGVSp_Short	Tumor_Sample_Barcode	t_alt_count	t_ref_count
KRAS	12	25398284	25398284	Missense_Mutation	p.G12D	CASE-KRAS	90	110
```

`egfr_amp_mut.maf`:
```
Hugo_Symbol	Chromosome	Start_Position	End_Position	Variant_Classification	HGVSp_Short	Tumor_Sample_Barcode	t_alt_count	t_ref_count
EGFR	7	55259515	55259515	Missense_Mutation	p.L858R	CASE-EGFR	120	80
```

`brca_compound_het.maf`:
```
Hugo_Symbol	Chromosome	Start_Position	End_Position	Variant_Classification	HGVSp_Short	Tumor_Sample_Barcode	t_alt_count	t_ref_count
BRCA1	17	41244000	41244000	Frame_Shift_Ins	p.E1257fs	CASE-BRCA	60	140
BRCA1	17	41251931	41251931	Splice_Site	p.X1493_splice	CASE-BRCA	35	65
```

- [ ] **Step 2: Create the four SEG files**

`tp53_loh.seg` (chr17p hemizygous deletion overlaps TP53 at 7,577,120):
```
ID	chrom	loc.start	loc.end	num.mark	seg.mean
CASE-TP53-LOH	17	1	25000000	200	-0.40
CASE-TP53-LOH	17	25000001	81195210	300	0.02
CASE-TP53-LOH	7	1	159138663	500	0.01
CASE-TP53-LOH	12	1	133851895	450	-0.03
```

`kras_mut_only.seg` (neutral over KRAS):
```
ID	chrom	loc.start	loc.end	num.mark	seg.mean
CASE-KRAS	12	1	133851895	450	0.05
CASE-KRAS	17	1	81195210	300	0.01
CASE-KRAS	7	1	159138663	500	0.00
```

`egfr_amp_mut.seg` (chr7 amplification overlaps EGFR at 55,259,515):
```
ID	chrom	loc.start	loc.end	num.mark	seg.mean
CASE-EGFR	7	1	159138663	500	1.20
CASE-EGFR	12	1	133851895	450	0.02
CASE-EGFR	17	1	81195210	300	0.01
```

`brca_compound_het.seg` (neutral; two truncating mutations carry the call):
```
ID	chrom	loc.start	loc.end	num.mark	seg.mean
CASE-BRCA	17	1	81195210	300	0.01
CASE-BRCA	7	1	159138663	500	0.02
```

- [ ] **Step 3: Verify the cases produce the expected calls**

Run:
```bash
uv run python -c "
from pathlib import Path
from two_hit.io import parse_maf, parse_seg
from two_hit.analysis import analyze_sample
from two_hit.models import AnalysisParams
from two_hit.web_demo.panel import PANEL, filter_maf_to_panel
d = Path('src/two_hit/web_demo/cases')
expect = {
  'tp53_loh': ('TP53','biallelic_mut_loh_high_vaf'),
  'kras_mut_only': ('KRAS','mut_only'),
  'egfr_amp_mut': ('EGFR','mut_amp'),
  'brca_compound_het': ('BRCA1','biallelic_compound_het'),
}
for name,(gene,status) in expect.items():
    maf = filter_maf_to_panel(parse_maf(d/f'{name}.maf'))
    seg = parse_seg(d/f'{name}.seg')
    r = analyze_sample(maf, seg, params=AnalysisParams(), gene_roles=PANEL)
    g = next(x for x in r.gene_results if x.hugo_symbol==gene)
    print(name, g.hugo_symbol, g.two_hit_status.value)
    assert g.two_hit_status.value == status, (name, g.two_hit_status.value)
print('ALL CASES OK')
"
```
Expected: prints each case and `ALL CASES OK`.

- [ ] **Step 4: Commit**

```bash
git add src/two_hit/web_demo/cases/
git commit -m "feat: add preset demo cases (TP53 LOH, KRAS, EGFR amp, BRCA compound het)"
```

---

## Task 3: Server-side SVG genome plot

**Files:**
- Create: `src/two_hit/web_demo/svgplot.py`
- Test: `tests/test_svgplot.py`

- [ ] **Step 1: Write the failing test**

Create `tests/test_svgplot.py`:

```python
from pathlib import Path

from two_hit.analysis import analyze_sample
from two_hit.io import parse_maf, parse_seg
from two_hit.models import AnalysisParams
from two_hit.web_demo.panel import PANEL, filter_maf_to_panel
from two_hit.web_demo.svgplot import genome_svg

CASES = Path("src/two_hit/web_demo/cases")


def _run(name):
    maf = filter_maf_to_panel(parse_maf(CASES / f"{name}.maf"))
    seg = parse_seg(CASES / f"{name}.seg")
    result = analyze_sample(maf, seg, params=AnalysisParams(), gene_roles=PANEL)
    return seg, result


def test_genome_svg_is_well_formed_svg():
    seg, result = _run("egfr_amp_mut")
    svg = genome_svg(seg, result)
    assert svg.lstrip().startswith("<svg")
    assert svg.rstrip().endswith("</svg>")
    # amplified segment should be drawn (a <line ...> element exists)
    assert "<line" in svg
    # mutation marker(s) drawn as circles with hover titles
    assert "<circle" in svg
    assert "EGFR" in svg  # mutation hover title mentions the gene


def test_genome_svg_handles_empty_segments():
    seg, result = _run("egfr_amp_mut")
    empty = seg.clear()  # same schema, zero rows
    svg = genome_svg(empty, result)
    assert svg.lstrip().startswith("<svg")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_svgplot.py -v`
Expected: FAIL — `ModuleNotFoundError: two_hit.web_demo.svgplot`.

- [ ] **Step 3: Write minimal implementation**

Create `src/two_hit/web_demo/svgplot.py`:

```python
"""Server-side inline SVG genome-wide CNA plot.

Replaces the plotly figure used by the full app. Produces a self-contained
<svg> string (no JS, no CDN) with copy-number segments and mutation markers
laid out across hg19 chromosomes.
"""

from __future__ import annotations

from html import escape

import polars as pl

from ..models import SampleResult

# hg19 chromosome lengths (bp)
HG19_CHROM_LENGTHS: dict[str, int] = {
    "1": 249250621, "2": 243199373, "3": 198022430, "4": 191154276,
    "5": 180915260, "6": 171115067, "7": 159138663, "8": 146364022,
    "9": 141213431, "10": 135534747, "11": 135006516, "12": 133851895,
    "13": 115169878, "14": 107349540, "15": 102531392, "16": 90354753,
    "17": 81195210, "18": 78077248, "19": 59128983, "20": 63025520,
    "21": 48129895, "22": 51304566, "X": 155270560, "Y": 59373566,
}

_CHROM_ORDER = list(HG19_CHROM_LENGTHS)
_YMAX = 1.5  # clamp seg.mean to [-1.5, 1.5]

_BIALLELIC = {
    "biallelic_compound_het", "biallelic_mut_loh",
    "biallelic_mut_loh_high_vaf", "homdel",
}


def _norm_chrom(c: str) -> str:
    c = str(c).strip().upper()
    return c[3:] if c.startswith("CHR") else c


def _seg_color(seg_mean: float | None, params) -> str:
    if seg_mean is None:
        return "#9aa0a6"
    if seg_mean < params.del_threshold:
        return "#3a6ea5"  # deletion (blue)
    if seg_mean > params.amp_threshold:
        return "#c0392b"  # amplification (red)
    if seg_mean > params.gain_threshold:
        return "#e08e0b"  # gain (orange)
    return "#9aa0a6"  # neutral (grey)


def _mut_color(status: str) -> str:
    if status in _BIALLELIC:
        return "#6c3483"  # biallelic (purple)
    if status == "mut_amp":
        return "#c0392b"  # mut + amp (red)
    return "#444"


def genome_svg(seg_df: pl.DataFrame, result: SampleResult, width: int = 920, height: int = 320) -> str:
    """Render a genome-wide CNA plot with mutation overlay as an inline SVG string."""
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
            f'<text x="{margin_l - 6}" y="{yt + 3:.1f}" text-anchor="end" fill="#888">{tick:+.0f}</text>'
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
            parts.append(
                f'<line x1="{x1:.1f}" y1="{yv:.1f}" x2="{x2:.1f}" y2="{yv:.1f}" '
                f'stroke="{_seg_color(float(sm), params)}" stroke-width="2.5">'
                f'<title>chr{chrom} seg.mean={float(sm):.2f}</title></line>'
            )

    # mutation markers
    for g in result.gene_results:
        status = g.two_hit_status.value
        color = _mut_color(status)
        for m in g.mutations:
            chrom = _norm_chrom(m.chromosome)
            if chrom not in HG19_CHROM_LENGTHS:
                continue
            x = gx(chrom, m.start_position)
            yv = gy(m.seg_mean if m.seg_mean is not None else 0.0)
            label = escape(f"{g.hugo_symbol} {m.protein_change or m.variant_classification} [{g.status_label}]")
            parts.append(
                f'<circle cx="{x:.1f}" cy="{yv:.1f}" r="3.5" fill="{color}" '
                f'stroke="#fff" stroke-width="0.8"><title>{label}</title></circle>'
            )

    parts.append("</svg>")
    return "".join(parts)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run pytest tests/test_svgplot.py -v`
Expected: PASS (2 tests).

- [ ] **Step 5: Commit**

```bash
git add src/two_hit/web_demo/svgplot.py tests/test_svgplot.py
git commit -m "feat: add server-side SVG genome plot for demo"
```

---

## Task 4: FastAPI app + routes

**Files:**
- Create: `src/two_hit/web_demo/app.py`
- Create: `src/two_hit/web_demo/routes.py`
- Test: `tests/test_web_demo.py`

> Templates do not exist yet — Task 5 creates them. This task's tests that render
> HTML are written here but will only pass after Task 5. Step 4 runs just the
> non-rendering test; the full suite is re-run at the end of Task 5.

- [ ] **Step 1: Write the failing test**

Create `tests/test_web_demo.py`:

```python
import io

import pytest
from fastapi.testclient import TestClient

from two_hit.web_demo import create_app


@pytest.fixture
def client():
    return TestClient(create_app())


def test_index_lists_cases(client):
    r = client.get("/")
    assert r.status_code == 200
    assert "TP53" in r.text
    assert "/case/tp53_loh" in r.text


@pytest.mark.parametrize(
    "name,gene",
    [
        ("tp53_loh", "TP53"),
        ("kras_mut_only", "KRAS"),
        ("egfr_amp_mut", "EGFR"),
        ("brca_compound_het", "BRCA1"),
    ],
)
def test_case_renders_report(client, name, gene):
    r = client.get(f"/case/{name}")
    assert r.status_code == 200
    assert gene in r.text
    assert "<svg" in r.text


def test_unknown_case_404(client):
    r = client.get("/case/does_not_exist")
    assert r.status_code == 404


def test_analyze_rejects_oversized_file(client):
    big = io.BytesIO(b"x" * (1_000_001))
    small = io.BytesIO(b"ID\tchrom\tloc.start\tloc.end\tseg.mean\n")
    r = client.post(
        "/analyze",
        files={
            "maf_file": ("big.maf", big, "text/plain"),
            "seg_file": ("s.seg", small, "text/plain"),
        },
    )
    assert r.status_code == 400
    assert "1 MB" in r.text or "too large" in r.text.lower()


def test_analyze_valid_upload_returns_report(client):
    maf = (
        "Hugo_Symbol\tChromosome\tStart_Position\tEnd_Position\t"
        "Variant_Classification\tHGVSp_Short\tTumor_Sample_Barcode\t"
        "t_alt_count\tt_ref_count\n"
        "KRAS\t12\t25398284\t25398284\tMissense_Mutation\tp.G12D\tUP\t90\t110\n"
    )
    seg = (
        "ID\tchrom\tloc.start\tloc.end\tnum.mark\tseg.mean\n"
        "UP\t12\t1\t133851895\t450\t0.05\n"
    )
    r = client.post(
        "/analyze",
        files={
            "maf_file": ("u.maf", io.BytesIO(maf.encode()), "text/plain"),
            "seg_file": ("u.seg", io.BytesIO(seg.encode()), "text/plain"),
        },
    )
    assert r.status_code == 200
    assert "KRAS" in r.text


def test_analyze_off_panel_only_shows_note(client):
    # NOTCH1 is not in the panel -> nothing to report
    maf = (
        "Hugo_Symbol\tChromosome\tStart_Position\tEnd_Position\t"
        "Variant_Classification\tHGVSp_Short\tTumor_Sample_Barcode\t"
        "t_alt_count\tt_ref_count\n"
        "NOTCH1\t9\t139399365\t139399365\tNonsense_Mutation\tp.R1598*\tUP\t40\t60\n"
    )
    seg = "ID\tchrom\tloc.start\tloc.end\tnum.mark\tseg.mean\nUP\t9\t1\t141213431\t450\t0.0\n"
    r = client.post(
        "/analyze",
        files={
            "maf_file": ("u.maf", io.BytesIO(maf.encode()), "text/plain"),
            "seg_file": ("u.seg", io.BytesIO(seg.encode()), "text/plain"),
        },
    )
    assert r.status_code == 200
    assert "no panel" in r.text.lower()
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_web_demo.py -v`
Expected: FAIL — `ModuleNotFoundError` / `create_app` import error.

- [ ] **Step 3: Write minimal implementation**

Create `src/two_hit/web_demo/app.py`:

```python
"""FastAPI application factory for the lightweight two-hit demo.

No idle watchdog, no PDF, no JSON API, no result store — kept small enough to
run on a modest Railway container.
"""

from __future__ import annotations

from pathlib import Path

from fastapi import FastAPI
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

from .. import __version__

WEB_DIR = Path(__file__).parent
TEMPLATE_DIR = WEB_DIR / "templates"
STATIC_DIR = WEB_DIR / "static"


def create_app() -> FastAPI:
    app = FastAPI(
        title="two-hit (demo)",
        description="Lightweight demo of integrated mutation + copy number analysis",
        version=__version__,
    )
    app.mount("/static", StaticFiles(directory=str(STATIC_DIR)), name="static")

    templates = Jinja2Templates(directory=str(TEMPLATE_DIR))
    templates.env.globals["version"] = __version__
    app.state.templates = templates

    from .routes import router

    app.include_router(router)
    return app


app = create_app()
```

Create `src/two_hit/web_demo/routes.py`:

```python
"""Routes for the two-hit demo: landing, preset cases, bounded upload."""

from __future__ import annotations

import logging
from pathlib import Path

from fastapi import APIRouter, File, Request, UploadFile
from fastapi.responses import HTMLResponse

logger = logging.getLogger(__name__)

router = APIRouter()

CASES_DIR = Path(__file__).parent / "cases"
MAX_BYTES = 1_000_000  # 1 MB per file

# name -> (human label, short description)
CASES: dict[str, tuple[str, str]] = {
    "tp53_loh": ("TP53 biallelic (mutation + LOH)", "Truncating TP53 mutation with high VAF over a chr17p deletion."),
    "kras_mut_only": ("KRAS mutation only", "Activating KRAS G12D with no supporting copy-number change."),
    "egfr_amp_mut": ("EGFR mutation + amplification", "EGFR L858R on an amplified chr7 segment."),
    "brca_compound_het": ("BRCA1 compound heterozygous", "Two truncating BRCA1 mutations (compound het)."),
}


def _render_report(request: Request, maf_src, seg_src, note: str | None = None):
    """Run the panel-restricted analysis and render the report."""
    from ..analysis import analyze_sample
    from ..io import parse_maf, parse_seg
    from ..models import AnalysisParams
    from .panel import PANEL, filter_maf_to_panel
    from .svgplot import genome_svg

    maf_df = filter_maf_to_panel(parse_maf(maf_src))
    seg_df = parse_seg(seg_src)

    sample_id = None
    if not maf_df.is_empty() and "Tumor_Sample_Barcode" in maf_df.columns:
        ids = maf_df["Tumor_Sample_Barcode"].unique().to_list()
        sample_id = ids[0]
        if len(ids) > 1:
            maf_df = maf_df.filter(pl_eq("Tumor_Sample_Barcode", sample_id))
            note = (note + " " if note else "") + f"Multiple samples found; showing {sample_id}."

    result = analyze_sample(maf_df, seg_df, sample_id=sample_id, params=AnalysisParams(), gene_roles=PANEL)
    svg = genome_svg(seg_df, result)

    templates = request.app.state.templates
    return templates.TemplateResponse(
        request, "report.html", {"r": result, "plot_svg": svg, "note": note}
    )


def pl_eq(col: str, val):
    """Tiny helper to avoid importing polars at module load."""
    import polars as pl

    return pl.col(col) == val


def _index(request: Request, error: str | None = None, status_code: int = 200):
    templates = request.app.state.templates
    return templates.TemplateResponse(
        request, "index.html", {"cases": CASES, "error": error}, status_code=status_code
    )


@router.get("/", response_class=HTMLResponse)
async def index(request: Request):
    return _index(request)


@router.get("/case/{name}", response_class=HTMLResponse)
async def case(request: Request, name: str):
    if name not in CASES:
        return _index(request, error=f"Unknown case: {name}", status_code=404)
    maf_path = CASES_DIR / f"{name}.maf"
    seg_path = CASES_DIR / f"{name}.seg"
    if not maf_path.exists() or not seg_path.exists():
        return _index(request, error="Case data missing.", status_code=404)
    try:
        return _render_report(request, maf_path, seg_path)
    except Exception as e:  # noqa: BLE001
        logger.error("Case %s failed: %s", name, e, exc_info=True)
        return _index(request, error=str(e), status_code=400)


@router.post("/analyze", response_class=HTMLResponse)
async def analyze(
    request: Request,
    maf_file: UploadFile = File(...),
    seg_file: UploadFile = File(...),
):
    maf_bytes = await maf_file.read()
    seg_bytes = await seg_file.read()
    if len(maf_bytes) > MAX_BYTES or len(seg_bytes) > MAX_BYTES:
        return _index(
            request,
            error="File too large — the demo accepts files up to 1 MB each. Run the full app locally for larger inputs.",
            status_code=400,
        )
    try:
        maf_src = maf_bytes.decode("utf-8")
        seg_src = seg_bytes.decode("utf-8")
        return _render_report(request, maf_src, seg_src)
    except Exception as e:  # noqa: BLE001
        logger.error("Upload analysis failed: %s", e, exc_info=True)
        return _index(request, error=str(e), status_code=400)
```

- [ ] **Step 4: Run the non-rendering test to verify the app imports**

Run: `uv run pytest tests/test_web_demo.py::test_analyze_rejects_oversized_file -v`
Expected: PASS (this path returns the error before touching templates — but it renders `index.html`, which does not exist yet, so it will FAIL with a template error). **If it fails on a missing template**, that is expected; proceed to Task 5 and re-run the full suite there.

- [ ] **Step 5: Commit**

```bash
git add src/two_hit/web_demo/app.py src/two_hit/web_demo/routes.py tests/test_web_demo.py
git commit -m "feat: add demo FastAPI app and routes"
```

---

## Task 5: Templates and static assets

**Files:**
- Create: `src/two_hit/web_demo/templates/base.html`
- Create: `src/two_hit/web_demo/templates/index.html`
- Create: `src/two_hit/web_demo/templates/report.html`
- Create: `src/two_hit/web_demo/static/demo.css`
- Create: `src/two_hit/web_demo/static/molpath.css` (copied)
- Create: `src/two_hit/web_demo/static/favicon.svg` (copied)

- [ ] **Step 1: Copy shared static assets**

Run:
```bash
mkdir -p src/two_hit/web_demo/static
cp src/two_hit/web/static/molpath.css src/two_hit/web_demo/static/molpath.css
cp src/two_hit/web/static/favicon.svg src/two_hit/web_demo/static/favicon.svg
```

- [ ] **Step 2: Create `base.html`**

```html
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{% block title %}two-hit demo{% endblock %}</title>
  <link rel="icon" type="image/svg+xml" href="/static/favicon.svg">
  <link href="https://fonts.googleapis.com/css2?family=JetBrains+Mono:ital,wght@0,300;0,400;0,500;0,600;0,700;1,400&display=swap" rel="stylesheet">
  <link rel="stylesheet" href="/static/molpath.css">
  <link rel="stylesheet" href="/static/demo.css">
  {% block head %}{% endblock %}
</head>
<body>
  <div class="mp-background-glow"></div>
  <header class="site-header">
    <div class="mp-layout" style="margin-top: var(--mp-space-sm); margin-bottom: 0;">
      <nav class="mp-flex mp-between">
        <a href="/" class="brand" style="text-decoration: none;">
          <span class="brand-mark">2H</span>
          <span class="brand-name">two-hit <span class="demo-tag">demo</span></span>
        </a>
      </nav>
    </div>
  </header>
  <main class="mp-layout mp-animate-in">
    {% block content %}{% endblock %}
  </main>
  <footer class="site-footer">
    <div class="mp-layout" style="margin-top: 0;">
      <p class="mp-soft" style="font-size: var(--mp-fs-xs); line-height: 1.6;">
        two-hit demo v{{ version }} — restricted to a 25-gene panel and ≤1 MB uploads. For research use and clinical decision support only. Not a standalone diagnostic device. Outputs may be incomplete or incorrect and must be reviewed by qualified experts in the full clinical, laboratory, and histopathological context. The authors assume no responsibility for clinical decisions or other consequences of use.<br>
        Part of <a href="https://molpath.tools" style="color: var(--mp-accent);">molpath.tools</a>
      </p>
    </div>
  </footer>
</body>
</html>
```

- [ ] **Step 3: Create `index.html`**

```html
{% extends "base.html" %}
{% block title %}two-hit demo{% endblock %}
{% block content %}
<div class="mp-panel">
  <p class="mp-eyebrow">Biallelic inactivation analysis — demo</p>
  <h1 style="margin-top: var(--mp-space-xs);">two-hit <span class="demo-tag">demo</span></h1>
  <p class="mp-subtitle">A lightweight demo restricted to a 25-gene panel. Try a preset case, or upload your own MAF + SEG (≤1 MB each). For the full tool — all OncoKB genes, PDF export, larger inputs — run it locally.</p>

  {% if error %}
  <div class="mp-alert mp-alert-error" style="margin-top: var(--mp-space-md);">
    <strong>Error:</strong> {{ error }}
  </div>
  {% endif %}

  <h2 style="margin-top: var(--mp-space-lg);">Preset cases</h2>
  <div class="case-grid">
    {% for name, info in cases.items() %}
    <a class="case-card" href="/case/{{ name }}">
      <span class="case-title">{{ info[0] }}</span>
      <span class="case-desc">{{ info[1] }}</span>
    </a>
    {% endfor %}
  </div>

  <h2 style="margin-top: var(--mp-space-lg);">Upload your own</h2>
  <form action="/analyze" method="post" enctype="multipart/form-data" class="upload-form">
    <div class="upload-grid">
      <label class="file-field">
        <span class="mp-section-title">MAF file</span>
        <input type="file" name="maf_file" accept=".maf,.tsv,.txt" required>
        <span class="file-hint">Somatic mutations, ≤1 MB</span>
      </label>
      <label class="file-field">
        <span class="mp-section-title">SEG file</span>
        <input type="file" name="seg_file" accept=".seg,.tsv,.txt" required>
        <span class="file-hint">Copy number segments, ≤1 MB</span>
      </label>
    </div>
    <button type="submit" class="mp-btn-primary" style="margin-top: var(--mp-space-md);">Analyze</button>
  </form>

  <p class="mp-soft" style="font-size: var(--mp-fs-xs); margin-top: var(--mp-space-md);">
    Panel: TP53, RB1, PTEN, APC, BRCA1, BRCA2, CDKN2A, NF1, VHL, STK11, SMAD4, KEAP1, ATM, MLH1, KRAS, NRAS, BRAF, EGFR, PIK3CA, ERBB2, MET, MYC, CCND1, MDM2, CTNNB1. Off-panel genes are ignored.
  </p>
</div>
{% endblock %}
```

- [ ] **Step 4: Create `report.html`**

```html
{% extends "base.html" %}
{% block title %}Report — {{ r.sample_id }} — two-hit demo{% endblock %}
{% block content %}
<div class="mp-panel">
  <p class="mp-eyebrow">Analysis report (demo)</p>
  <h1 style="margin-top: var(--mp-space-xs);">{{ r.sample_id }}</h1>

  {% if note %}
  <div class="mp-alert" style="margin-top: var(--mp-space-sm);">{{ note }}</div>
  {% endif %}

  <div class="summary-stats" style="margin-top: var(--mp-space-md);">
    <div class="stat"><span class="stat-value">{{ r.total_mutations }}</span><span class="stat-label">Mutations</span></div>
    <div class="stat"><span class="stat-value">{{ r.total_genes_with_events }}</span><span class="stat-label">Genes</span></div>
    <div class="stat"><span class="stat-value biallelic">{{ r.biallelic_count }}</span><span class="stat-label">Biallelic</span></div>
    <div class="stat"><span class="stat-value oncogene">{{ r.mut_amp_count }}</span><span class="stat-label">Mut+Amp</span></div>
  </div>
  <div class="mp-flex" style="margin-top: var(--mp-space-md); gap: var(--mp-space-sm);">
    <a href="/" class="mp-btn-secondary" style="text-decoration: none; font-size: var(--mp-fs-sm);">New analysis</a>
  </div>
</div>

<div class="mp-panel">
  <h2>Gene results</h2>
  {% if r.gene_results %}
  <div style="overflow-x: auto; margin-top: var(--mp-space-sm);">
    <table class="gene-table">
      <thead>
        <tr><th>Gene</th><th>Role</th><th>#Mut</th><th>Best mutation</th><th>CNA</th><th>seg.mean</th><th>VAF</th><th>Status</th><th>Confidence</th></tr>
      </thead>
      <tbody>
        {% for g in r.gene_results %}
        <tr>
          <td><strong>{{ g.hugo_symbol }}</strong></td>
          <td>{{ g.role.value }}</td>
          <td>{{ g.mutations | length }}</td>
          <td>{{ g.best_mutation }}</td>
          <td>{{ g.cna_status.value }}</td>
          <td>{% if g.seg_mean is not none %}{{ "%.3f" | format(g.seg_mean) }}{% else %}N/A{% endif %}</td>
          <td>{% if g.vaf_max is not none %}{{ "%.2f" | format(g.vaf_max) }}{% else %}N/A{% endif %}</td>
          <td><span class="status-badge status-{{ g.two_hit_status.value }}">{{ g.status_label }}</span></td>
          <td><span class="conf-badge conf-{{ g.confidence.value }}">{{ g.confidence.value }}</span></td>
        </tr>
        {% endfor %}
      </tbody>
    </table>
  </div>
  {% else %}
  <p class="mp-soft" style="margin-top: var(--mp-space-sm);">No panel genes mutated in this sample — nothing to report. (The demo only analyzes its 25-gene panel.)</p>
  {% endif %}
</div>

<div class="mp-panel">
  <h2>Genome-wide CNA plot</h2>
  <div class="plot-container" style="margin-top: var(--mp-space-sm); overflow-x:auto;">
    {{ plot_svg | safe }}
  </div>
  <p class="mp-soft" style="font-size: var(--mp-fs-xs); margin-top: var(--mp-space-xs);">
    Blue = deletion · Red = amplification · Orange = gain · Grey = neutral. Circles are mutations (purple = biallelic, red = mut+amp). Hover for details.
  </p>
</div>
{% endblock %}
```

- [ ] **Step 5: Create `demo.css`**

```css
/* two-hit demo — small layer on top of molpath.css */

.demo-tag {
  font-family: var(--mp-mono, monospace);
  font-size: 0.55em;
  text-transform: uppercase;
  letter-spacing: 0.08em;
  color: var(--mp-accent);
  border: 1px solid var(--mp-accent);
  border-radius: 4px;
  padding: 1px 5px;
  vertical-align: middle;
}

.case-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(240px, 1fr));
  gap: var(--mp-space-sm);
  margin-top: var(--mp-space-sm);
}
.case-card {
  display: flex;
  flex-direction: column;
  gap: 4px;
  padding: var(--mp-space-md);
  border: 1px solid var(--mp-border);
  border-radius: 8px;
  background: var(--mp-bg-1);
  text-decoration: none;
  transition: border-color 0.15s ease, transform 0.15s ease;
}
.case-card:hover { border-color: var(--mp-accent); transform: translateY(-1px); }
.case-title { font-weight: 600; color: var(--mp-ink); }
.case-desc { font-size: var(--mp-fs-sm); color: var(--mp-soft, #666); }

.upload-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
  gap: var(--mp-space-md);
}
.file-field { display: flex; flex-direction: column; gap: 4px; }
.file-hint { font-size: var(--mp-fs-xs); color: var(--mp-soft, #888); }

.summary-stats { display: flex; gap: var(--mp-space-lg); flex-wrap: wrap; }
.stat { display: flex; flex-direction: column; }
.stat-value { font-size: var(--mp-fs-xl); font-weight: 700; }
.stat-value.biallelic { color: #6c3483; }
.stat-value.oncogene { color: #c0392b; }
.stat-label { font-size: var(--mp-fs-xs); text-transform: uppercase; letter-spacing: 0.06em; color: var(--mp-soft, #888); }

.gene-table { width: 100%; border-collapse: collapse; font-size: var(--mp-fs-sm); }
.gene-table th, .gene-table td { text-align: left; padding: 6px 10px; border-bottom: 1px solid var(--mp-border); }
.gene-table th { font-family: var(--mp-mono, monospace); font-size: var(--mp-fs-xs); text-transform: uppercase; letter-spacing: 0.05em; color: var(--mp-soft, #888); }

.status-badge, .conf-badge {
  display: inline-block; padding: 2px 8px; border-radius: 999px;
  font-size: var(--mp-fs-xs); font-family: var(--mp-mono, monospace);
  background: var(--mp-bg-1); border: 1px solid var(--mp-border);
}
.conf-high { color: #1e7e34; border-color: #1e7e34; }
.conf-medium { color: #b8860b; border-color: #b8860b; }
.conf-low { color: #888; }
```

- [ ] **Step 6: Run the full demo test suite**

Run: `uv run pytest tests/test_web_demo.py tests/test_svgplot.py tests/test_panel.py -v`
Expected: PASS (all tests).

- [ ] **Step 7: Commit**

```bash
git add src/two_hit/web_demo/templates/ src/two_hit/web_demo/static/
git commit -m "feat: add demo templates and static assets"
```

---

## Task 6: Packaging, Dockerfile, Railway config, README

**Files:**
- Modify: `pyproject.toml` (ensure `cases/`, `static/`, `templates/` ship in the wheel)
- Modify: `Dockerfile`
- Delete: `railway.toml`
- Modify: `README.md`

- [ ] **Step 1: Ensure non-Python package data is included**

Read `pyproject.toml`. Under `[tool.hatch.build.targets.wheel]` add a `force-include` (or `artifacts`) so cases/templates/static ship. Replace the existing `[tool.hatch.build.targets.wheel]` block with:

```toml
[tool.hatch.build.targets.wheel]
packages = ["src/two_hit"]

[tool.hatch.build.targets.wheel.force-include]
"src/two_hit/web_demo/cases" = "two_hit/web_demo/cases"
"src/two_hit/web_demo/templates" = "two_hit/web_demo/templates"
"src/two_hit/web_demo/static" = "two_hit/web_demo/static"
```

> The full app's `web/templates` and `web/static` are already resolved at runtime via `Path(__file__).parent`; the container COPYs the whole `src/` tree, so this force-include is belt-and-suspenders for wheel installs. Keep it.

- [ ] **Step 2: Point the Dockerfile at the demo app and drop the idle env**

Replace the tail of `Dockerfile` (everything from `ENV PORT=8080` onward) with:

```dockerfile
# Railway sets PORT env var; default to 8080
ENV PORT=8080
EXPOSE 8080

# Serve the lightweight demo app (full app lives in two_hit.web for local use)
CMD exec uvicorn two_hit.web_demo.app:app --host 0.0.0.0 --port $PORT
```

- [ ] **Step 3: Remove the Railway restart-retry hack**

Run:
```bash
git rm railway.toml
```

- [ ] **Step 4: Document the demo vs. full split in README**

Read `README.md`, then add a section (after the existing intro/usage) titled `## Demo vs. full app`:

```markdown
## Demo vs. full app

This repo ships two web frontends over the same analysis core:

- **Full app** (`two_hit.web`) — all OncoKB genes, arbitrary uploads, interactive
  Plotly plot, PDF export, JSON API. Heavier; intended for local use:

  ```bash
  uv run uvicorn two_hit.web.app:app --reload
  ```

- **Demo app** (`two_hit.web_demo`) — restricted to a 25-gene panel and ≤1 MB
  uploads, server-side SVG plot, no PDF/API. Lightweight; this is what is
  deployed at `two-hit.molpath.tools`:

  ```bash
  uv run uvicorn two_hit.web_demo.app:app --reload
  ```

The demo offers preset cases (TP53 mutation+LOH, KRAS mutation-only, EGFR
mutation+amplification, BRCA1 compound heterozygous) so visitors can explore
without uploading data.
```

- [ ] **Step 5: Verify the build/test still passes and the app boots**

Run:
```bash
uv run pytest -q
uv run python -c "from two_hit.web_demo.app import app; print('routes:', [r.path for r in app.routes])"
```
Expected: tests pass; printed routes include `/`, `/case/{name}`, `/analyze`, `/static`.

- [ ] **Step 6: Commit**

```bash
git add pyproject.toml Dockerfile README.md
git commit -m "build: serve lightweight demo on Railway; drop idle-restart hack"
```

---

## Task 7: Local smoke test + Railway deploy

**Files:** none (operational)

- [ ] **Step 1: Boot the demo locally and confirm memory footprint**

Run (background, then probe):
```bash
uv run uvicorn two_hit.web_demo.app:app --host 127.0.0.1 --port 8099 &
sleep 3
curl -s -o /dev/null -w "index:%{http_code}\n" http://127.0.0.1:8099/
curl -s -o /dev/null -w "case:%{http_code}\n" http://127.0.0.1:8099/case/tp53_loh
# memory of the uvicorn worker
ps -o rss= -p $(pgrep -f "two_hit.web_demo.app") | awk '{printf "RSS: %.0f MB\n", $1/1024}'
kill %1
```
Expected: `index:200`, `case:200`, RSS well under ~120 MB (target ~60–90 MB).

- [ ] **Step 2: Merge to main and push**

```bash
git checkout main
git merge --no-ff feature/demo-app -m "feat: lightweight Railway demo app"
git push origin main
```

- [ ] **Step 3: Confirm the Railway deploy via MCP**

Using the Railway MCP tools (do not shell out):
1. `list_projects` / `list_services` → locate the two-hit service.
2. `list_deployments` → confirm a new deployment triggered by the push and that it reaches `SUCCESS`.
3. `get_logs` → confirm uvicorn started on `web_demo.app` with no tracebacks.
4. `service_metrics` → confirm memory sits low and stable (no restart loop).
5. Fetch `https://two-hit.molpath.tools/` and `/case/tp53_loh` → expect HTTP 200.

Expected: deployment `SUCCESS`, low steady memory, both URLs return 200.

- [ ] **Step 4: Report deploy status to the user**

Summarize: deployment id/status, observed memory, and the live URL. If the
deploy failed, capture the relevant log lines and stop for review rather than
retrying blindly.

---

## Self-Review

**Spec coverage:**
- Separate lightweight app, no plotly/fpdf2 → Tasks 1–5 (`web_demo`, SVG plot). ✓
- Curated cases + bounded uploads (≤1 MB, single sample, panel-filtered) → Task 2 (cases), Task 4 (`MAX_BYTES`, multi-sample note, `filter_maf_to_panel`). ✓
- Server-side inline SVG → Task 3. ✓
- ~25-gene panel, hardcoded → Task 1. ✓
- Replace Railway service, drop idle watchdog/restart hacks → Task 6 (Dockerfile, delete `railway.toml`; app.py has no watchdog), Task 7 (deploy). ✓
- Error handling (oversized, parse failure, empty-after-filter, unknown case) → Task 4 routes + Task 5 report template note. ✓
- Testing (panel, web_demo, svgplot) → Tasks 1, 3, 4/5. ✓
- README demo-vs-full section → Task 6. ✓

**Placeholder scan:** No TBD/TODO; every code step shows full code; commands have expected output. ✓

**Type consistency:** `PANEL` (dict[str, GeneRole]) and `filter_maf_to_panel` used identically in Tasks 1/3/4; `genome_svg(seg_df, result)` signature consistent across Task 3 and routes; templates use enum `.value` + `status_label`/`best_mutation` properties that exist on the models (verified against `models.py`); `analyze_sample(..., gene_roles=PANEL)` matches the real signature. ✓
