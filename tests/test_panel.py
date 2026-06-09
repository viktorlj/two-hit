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
