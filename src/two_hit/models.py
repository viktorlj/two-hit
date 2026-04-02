"""Pydantic v2 models for two-hit analysis pipeline."""

from __future__ import annotations

from enum import Enum

from pydantic import BaseModel, Field

from .gene_list import GeneRole

# ---------------------------------------------------------------------------
# Enums
# ---------------------------------------------------------------------------


class CnaStatus(str, Enum):
    """Copy number alteration status derived from seg.mean thresholds."""

    HOMDEL = "homdel"
    HEMIDEL = "hemidel"
    NEUTRAL = "neutral"
    GAIN = "gain"
    AMP = "amp"
    NO_DATA = "no_data"


class TwoHitStatus(str, Enum):
    """Two-hit classification outcome."""

    # TSG patterns
    BIALLELIC_COMPOUND_HET = "biallelic_compound_het"
    BIALLELIC_MUT_LOH = "biallelic_mut_loh"
    BIALLELIC_MUT_LOH_HIGH_VAF = "biallelic_mut_loh_high_vaf"
    HOMDEL = "homdel"
    SINGLE_HIT_MUT = "single_hit_mut"
    SINGLE_HIT_DEL = "single_hit_del"

    # Oncogene patterns
    MUT_AMP = "mut_amp"
    MUT_GAIN = "mut_gain"
    MUT_ONLY = "mut_only"
    AMP_ONLY = "amp_only"

    # Fallback
    UNCLASSIFIED = "unclassified"


class Confidence(str, Enum):
    """Evidence confidence level."""

    HIGH = "high"
    MEDIUM = "medium"
    LOW = "low"


# Status display metadata
STATUS_LABELS: dict[TwoHitStatus, str] = {
    TwoHitStatus.BIALLELIC_COMPOUND_HET: "Biallelic (compound het)",
    TwoHitStatus.BIALLELIC_MUT_LOH: "Biallelic (mut + LOH)",
    TwoHitStatus.BIALLELIC_MUT_LOH_HIGH_VAF: "Biallelic (mut + LOH, high VAF)",
    TwoHitStatus.HOMDEL: "Homozygous deletion",
    TwoHitStatus.SINGLE_HIT_MUT: "Single hit (mutation)",
    TwoHitStatus.SINGLE_HIT_DEL: "Single hit (deletion)",
    TwoHitStatus.MUT_AMP: "Mutation + amplification",
    TwoHitStatus.MUT_GAIN: "Mutation + gain",
    TwoHitStatus.MUT_ONLY: "Mutation only",
    TwoHitStatus.AMP_ONLY: "Amplification only",
    TwoHitStatus.UNCLASSIFIED: "Unclassified",
}

BIALLELIC_STATUSES = {
    TwoHitStatus.BIALLELIC_COMPOUND_HET,
    TwoHitStatus.BIALLELIC_MUT_LOH,
    TwoHitStatus.BIALLELIC_MUT_LOH_HIGH_VAF,
    TwoHitStatus.HOMDEL,
}

# ---------------------------------------------------------------------------
# Truncating mutation types
# ---------------------------------------------------------------------------

TRUNCATING_TYPES_BASE: set[str] = {
    "Nonsense_Mutation",
    "Frame_Shift_Del",
    "Frame_Shift_Ins",
    "Splice_Site",
    "Translation_Start_Site",
}

SPLICE_REGION = "Splice_Region"


# ---------------------------------------------------------------------------
# Analysis parameters
# ---------------------------------------------------------------------------


class AnalysisParams(BaseModel):
    """Configurable thresholds for two-hit analysis."""

    del_threshold: float = Field(-0.3, description="Hemizygous deletion / LOH threshold")
    deep_del_threshold: float = Field(-0.5, description="Homozygous deletion threshold")
    gain_threshold: float = Field(0.3, description="Copy number gain threshold")
    amp_threshold: float = Field(0.7, description="Amplification threshold")
    vaf_threshold: float = Field(0.7, description="High VAF threshold for LOH support")
    include_splice_region: bool = Field(False, description="Count Splice_Region as truncating")

    @property
    def truncating_types(self) -> set[str]:
        """Set of variant classifications considered truncating."""
        types = TRUNCATING_TYPES_BASE.copy()
        if self.include_splice_region:
            types.add(SPLICE_REGION)
        return types


# ---------------------------------------------------------------------------
# Data models
# ---------------------------------------------------------------------------


class MutationRecord(BaseModel):
    """A single somatic mutation enriched with CNA overlap data."""

    hugo_symbol: str
    chromosome: str = ""
    start_position: int = 0
    end_position: int = 0
    variant_classification: str = ""
    hgvsp_short: str = ""
    t_alt_count: int | None = None
    t_ref_count: int | None = None
    seg_mean: float | None = None
    cna_status: CnaStatus = CnaStatus.NO_DATA

    @property
    def vaf(self) -> float | None:
        """Variant allele frequency."""
        if self.t_alt_count is None or self.t_ref_count is None:
            return None
        total = self.t_alt_count + self.t_ref_count
        if total == 0:
            return None
        return self.t_alt_count / total

    def is_truncating(self, params: AnalysisParams) -> bool:
        """Whether this mutation is truncating given the current parameters."""
        return self.variant_classification in params.truncating_types

    @property
    def protein_change(self) -> str:
        """Normalized protein change (strip leading 'p.')."""
        hgvs = self.hgvsp_short.strip()
        if hgvs.startswith("p."):
            return hgvs[2:]
        return hgvs


class GeneResult(BaseModel):
    """Per-gene two-hit analysis result."""

    hugo_symbol: str
    role: GeneRole
    mutations: list[MutationRecord] = Field(default_factory=list)
    seg_mean: float | None = None
    cna_status: CnaStatus = CnaStatus.NO_DATA
    two_hit_status: TwoHitStatus = TwoHitStatus.UNCLASSIFIED
    confidence: Confidence = Confidence.LOW
    vaf_max: float | None = None
    notes: list[str] = Field(default_factory=list)

    @property
    def mutation_count(self) -> int:
        return len(self.mutations)

    @property
    def best_mutation(self) -> str:
        """Most informative mutation label (first truncating, or first overall)."""
        if not self.mutations:
            return ""
        for m in self.mutations:
            if m.hgvsp_short:
                return m.hgvsp_short
        return self.mutations[0].variant_classification

    @property
    def is_biallelic(self) -> bool:
        return self.two_hit_status in BIALLELIC_STATUSES

    @property
    def status_label(self) -> str:
        return STATUS_LABELS.get(self.two_hit_status, self.two_hit_status.value)


class SampleResult(BaseModel):
    """Complete two-hit analysis result for one sample."""

    sample_id: str
    params: AnalysisParams
    gene_results: list[GeneResult] = Field(default_factory=list)
    total_mutations: int = 0
    total_genes_with_events: int = 0
    biallelic_count: int = 0
    mut_amp_count: int = 0
    notes: list[str] = Field(default_factory=list)
    version: str = "0.1.0"
