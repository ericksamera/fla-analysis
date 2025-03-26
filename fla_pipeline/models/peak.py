# fla_pipeline/models/peak.py

from dataclasses import dataclass

@dataclass
class Peak:
    position: float
    intensity: float
    saturated: bool = False
    note: str = ""

