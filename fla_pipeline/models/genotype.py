# fla_pipeline/models/genotype.py

from dataclasses import dataclass, field
from typing import List

@dataclass
class GenotypeResult:
    marker: str
    alleles: List[int]
    confidence: float = 0.0
    qc_flags: List[str] = field(default_factory=list)
    strategy: str = ""

    def dict(self):
        return {
            "marker": self.marker,
            "alleles": self.alleles,
            "confidence": self.confidence,
            "qc_flags": self.qc_flags,
            "strategy": self.strategy,
        }

