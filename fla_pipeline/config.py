# fla_pipeline/config.py

from typing import Optional, Tuple, Dict, Any, List
from pydantic import BaseModel, Field, model_validator


class GlobalConfig(BaseModel):
    ploidy: int = 2
    min_peak_height: int = 1000
    min_peak_position: int = 15
    relative_peak_threshold: float = 0.3
    bin_tolerance: int = 1
    bin_extend: int = 2
    intensity_weight: float = 1.0
    ratio_weight: float = 1.0
    stutter_weight: float = 1.0
    stutter_ratio_threshold: float = 0.85
    stutter_tolerance: float = 0.5
    instrument_saturation_value: float = 30_000
    qc_confidence_scale: float = 0.75
    diploid_strategy: str = "probabilistic"  # "strict_stutter" or "lenient_ratio"

    liz_peak_penalty_ratio: float = 1.25
    min_genotype_confidence: float = 0.0
    stutter_step_decay: float = 0.5


class MarkerConfig(BaseModel):
    marker: str
    channel: str
    repeat_unit: Optional[int]
    bins: Tuple[int, int]
    stutter_direction: Optional[str] = "both"
    overrides: Optional[Dict[str, Any]] = None
    binning_enabled: bool = True

    def merged_with_global(self, global_cfg: GlobalConfig) -> Dict[str, Any]:
        """
        Merge marker-specific overrides with the global config to produce
        a dict of effective configuration values for this marker.
        """
        base = global_cfg.model_dump()  # Use Pydantic v2-compatible method
        if self.overrides:
            for key, value in self.overrides.items():
                if key in base:
                    base[key] = value
        return base

    @model_validator(mode="after")
    def validate_required_fields(self) -> "MarkerConfig":
        if self.binning_enabled:
            if not self.bins:
                raise ValueError(f"Marker '{self.marker}': bins are required when binning is enabled.")
            if self.repeat_unit is None:
                raise ValueError(f"Marker '{self.marker}': repeat_unit is required when binning is enabled.")

        if not self.channel:
            raise ValueError(f"Marker '{self.marker}': channel is required.")
        return self