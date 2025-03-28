# fla_pipeline/config/session_config.py

from pydantic import BaseModel
from fla_pipeline.config.global_config import GlobalConfig
from fla_pipeline.config.marker_config import MarkerConfig
from fla_pipeline.config.distance_config import DistanceConfig
from typing import List, Dict, Any, Optional

class SessionConfig(BaseModel):
    version: str
    samples: List[Dict[str, Any]]  # Serialized Sample.to_dict()
    marker_list: List[MarkerConfig]
    global_config: GlobalConfig
    distance_config: Optional[DistanceConfig] = None