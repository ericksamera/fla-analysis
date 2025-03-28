# fla_pipeline/config/distance_config.py

from pydantic import BaseModel
from typing import Literal, Optional

class DistanceConfig(BaseModel):
    metric: Literal["allele_match", "shared_alleles", "nei"] = "allele_match"
    min_confidence: float = 0.8
    min_sample_count: int = 2
    imputation_strategy: Optional[Literal["mode", "mean"]] = None
    confidence_filter_enabled: bool = True
    marker_exclusion_enabled: bool = True