# fla_pipeline/utils/qc_score.py

from typing import List
import numpy as np
from fla_pipeline.config import GlobalConfig

class QCScore:
    def __init__(self, config: GlobalConfig):
        self.config = config
        self.penalties: List[tuple[float, str]] = []
        self.total_penalty = 0.0

    def add(self, value: float, reason: str, weight: float = 1.0):
        weighted = value * weight
        self.total_penalty += weighted
        self.penalties.append((weighted, reason))

    def compute_confidence(self) -> float:
        conf = float(np.exp(-self.config.qc_confidence_scale * self.total_penalty))
        return round(min(1.0, max(0.0, conf)), 3)

    def summary_flags(self) -> List[str]:
        return [f"{reason} (-{penalty:.2f})" for penalty, reason in self.penalties]
