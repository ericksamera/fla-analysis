from dataclasses import dataclass, field
from fla_pipeline.models.peak import Peak
from fla_pipeline.models.genotype import GenotypeResult
from typing import Optional, Dict, Any
import numpy as np


@dataclass
class Sample:
    sample_id: str
    file_path: str
    fsa_data: Optional[Dict[str, Any]] = None
    peaks: Optional[Dict[str, list]] = None
    suppressed_peaks: Optional[Dict[str, Any]] = field(default_factory=dict)
    marker_results: Optional[Dict[str, GenotypeResult]] = None
    metadata: Dict[str, Any] = field(default_factory=dict)
    run_metrics: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict:
        return {
            "sample_id": self.sample_id,
            "file_path": self.file_path,
            "metadata": self.metadata,
            "run_metrics": self.run_metrics,
            "fsa_data": {
                "name": self.fsa_data.get("name"),
                "smap": self.fsa_data["smap"].tolist(),
                "channels": {
                    k: v.tolist() for k, v in self.fsa_data["channels"].items()
                }
            } if self.fsa_data else None,
            "peaks": {
                k: [p.to_dict() for p in v] for k, v in (self.peaks or {}).items()
            },
            "suppressed_peaks": self.suppressed_peaks or {},
            "marker_results": {
                k: v.to_dict() if hasattr(v, "to_dict") else v
                for k, v in (self.marker_results or {}).items()
            }
        }

    @staticmethod
    def from_dict(data: dict) -> "Sample":
        sample = Sample(
            sample_id=data["sample_id"],
            file_path=data["file_path"],
            metadata=data.get("metadata", {}),
            run_metrics=data.get("run_metrics", {}),
            suppressed_peaks=data.get("suppressed_peaks", {})
        )

        if "fsa_data" in data and data["fsa_data"]:
            sample.fsa_data = {
                "name": data["fsa_data"]["name"],
                "smap": np.array(data["fsa_data"]["smap"]),
                "channels": {
                    k: np.array(v) for k, v in data["fsa_data"]["channels"].items()
                }
            }

        if "peaks" in data:
            sample.peaks = {
                k: [Peak.from_dict(p) for p in v] for k, v in data["peaks"].items()
            }

        if "marker_results" in data:
            sample.marker_results = {
                k: GenotypeResult.from_dict(v) if isinstance(v, dict) else v
                for k, v in data["marker_results"].items()
            }

        return sample

    @property
    def dict(self):
        return self.to_dict()
