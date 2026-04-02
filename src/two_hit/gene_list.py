"""OncoKB cancer gene list loader.

Downloads and caches the OncoKB cancer gene list, providing a mapping from
Hugo Symbol to gene role (oncogene / TSG / both / unknown).
"""

from __future__ import annotations

import logging
from enum import Enum
from pathlib import Path

import polars as pl

logger = logging.getLogger(__name__)

ONCOKB_URL = "https://www.oncokb.org/api/v1/utils/cancerGeneList.txt"
DEFAULT_DATA_DIR = Path(__file__).resolve().parent.parent.parent / "data"
GENE_LIST_FILENAME = "oncokb_cancer_genes.tsv"


class GeneRole(str, Enum):
    """Gene role classification from OncoKB."""

    ONCOGENE = "oncogene"
    TSG = "tsg"
    BOTH = "both"
    UNKNOWN = "unknown"


_GENE_TYPE_MAP: dict[str, GeneRole] = {
    "ONCOGENE": GeneRole.ONCOGENE,
    "TSG": GeneRole.TSG,
    "ONCOGENE_AND_TSG": GeneRole.BOTH,
}


def download_gene_list(dest: Path) -> None:
    """Download the OncoKB cancer gene list to dest."""
    import httpx

    logger.info("Downloading OncoKB cancer gene list from %s", ONCOKB_URL)
    resp = httpx.get(ONCOKB_URL, follow_redirects=True, timeout=30)
    resp.raise_for_status()
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(resp.content)
    logger.info("Saved %d bytes to %s", len(resp.content), dest)


def load_gene_roles(data_dir: Path | None = None) -> dict[str, GeneRole]:
    """Load gene roles from the OncoKB cancer gene list.

    Returns a dict mapping Hugo_Symbol -> GeneRole.
    Downloads the gene list if not already cached.
    """
    data_dir = data_dir or DEFAULT_DATA_DIR
    path = data_dir / GENE_LIST_FILENAME

    if not path.exists():
        download_gene_list(path)

    df = pl.read_csv(path, separator="\t", infer_schema_length=0)

    roles: dict[str, GeneRole] = {}
    for row in df.iter_rows(named=True):
        symbol = (row.get("Hugo Symbol") or "").strip()
        gene_type = (row.get("Gene Type") or "").strip()
        if not symbol:
            continue
        roles[symbol] = _GENE_TYPE_MAP.get(gene_type, GeneRole.UNKNOWN)

    logger.info(
        "Loaded %d gene roles (%d oncogenes, %d TSGs, %d both)",
        len(roles),
        sum(1 for v in roles.values() if v == GeneRole.ONCOGENE),
        sum(1 for v in roles.values() if v == GeneRole.TSG),
        sum(1 for v in roles.values() if v == GeneRole.BOTH),
    )
    return roles
