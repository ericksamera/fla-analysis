# fla_pipeline/config/marker_config.py

from pydantic import BaseModel, model_validator
from typing import Optional, Tuple, Dict, Any
from fla_pipeline.config.global_config import GlobalConfig

class MarkerConfig(BaseModel):
    marker: str
    channel: str
    repeat_unit: Optional[int]
    bins: Tuple[int, int]
    stutter_direction: Optional[str] = "both"
    overrides: Optional[Dict[str, Any]] = None
    binning_enabled: bool = True

    def merged_with_global(self, global_cfg: GlobalConfig) -> Dict[str, Any]:
        base = global_cfg.model_dump()
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