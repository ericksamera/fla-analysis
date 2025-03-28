from dataclasses import dataclass, field
from typing import List, Optional, Dict, Any


@dataclass
class GenotypeResult:
    marker: str
    alleles: List[str]
    confidence: float = 1.0
    is_valid: bool = True
    strategy: Optional[str] = None
    qc_flags: List[str] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict:
        return {
            "marker": self.marker,
            "alleles": self.alleles,
            "confidence": float(self.confidence),
            "is_valid": bool(self.is_valid),
            "strategy": self.strategy,
            "qc_flags": self.qc_flags,
            "metadata": self.metadata
        }

    @staticmethod
    def from_dict(data: dict) -> "GenotypeResult":
        return GenotypeResult(
            marker=data["marker"],
            alleles=data["alleles"],
            confidence=float(data.get("confidence", 1.0)),
            is_valid=bool(data.get("is_valid", True)),
            strategy=data.get("strategy"),
            qc_flags=data.get("qc_flags", []),
            metadata=data.get("metadata", {})
        )

    @property
    def dict(self):
        return self.to_dict()
