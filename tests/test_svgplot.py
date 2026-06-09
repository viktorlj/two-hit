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


def test_genome_svg_marks_biallelic_in_purple():
    seg, result = _run("tp53_loh")
    svg = genome_svg(seg, result)
    assert 'fill="#6c3483"' in svg  # biallelic TP53 marker is purple
    assert "TP53" in svg
