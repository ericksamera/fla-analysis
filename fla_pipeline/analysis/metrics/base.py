# fla_pipeline/analysis/metrics/base.py

from abc import ABC, abstractmethod
from fla_pipeline.models.genotype import GenotypeResult
from fla_pipeline.config.marker_config import MarkerConfig
from typing import Optional

class BaseDistanceMetric(ABC):
    name: str

    def __init__(self, debug: bool = False):
        self.debug = debug
        self.logs = []

    @abstractmethod
    def compute(
        self,
        geno1: GenotypeResult,
        geno2: GenotypeResult,
        marker_cfg: MarkerConfig
    ) -> Optional[float]:
        ...

    def log(self, msg: str):
        if self.debug:
            self.logs.append(msg)

    def get_logs(self):
        return self.logs
