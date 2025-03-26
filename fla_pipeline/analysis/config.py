# fla_pipeline/analysis/config.py

from pydantic import BaseModel, Field
from typing import Literal

class DistanceConfig(BaseModel):
    metric: Literal["bruvo", "nei", "da", "dst", "fst", "dsw", "dmu2", "gst", "gst_star", "jost_d"] = "bruvo"
    min_confidence: float = Field(0.8, ge=0.0, le=1.0)
    min_sample_count: int = Field(2, ge=1)
    impute_missing: bool = False
    bias_correction: bool = False
    debug_mode: bool = False
