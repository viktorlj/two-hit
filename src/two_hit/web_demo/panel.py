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
