# fla_pipeline/config/global_config.py

from pydantic import BaseModel

class GlobalConfig(BaseModel):
    ploidy: int = 2
    min_peak_height: int = 1000
    min_peak_position: int = 5
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
    diploid_strategy: str = "probabilistic"
    liz_peak_penalty_ratio: float = 1.25
    min_genotype_confidence: float = 0.0
    stutter_step_decay: float = 0.5
