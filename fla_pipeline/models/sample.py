# fla_pipeline/models/sample.py

from dataclasses import dataclass, field
from fla_pipeline.models.genotype import GenotypeResult
from typing import Optional, Dict, Any

@dataclass
class Sample:
    sample_id: str
    file_path: str
    fsa_data: Optional[Dict[str, Any]] = None
    peaks: Optional[Dict[str, list]] = None
    marker_results: Optional[Dict[str, GenotypeResult]] = None
    metadata: Dict[str, Any] = field(default_factory=dict)
    run_metrics: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self):
        return {
            "sample_id": self.sample_id,
            "file_path": self.file_path,
            "metadata": self.metadata,
            "run_metrics": self.run_metrics,
            "marker_results": {
                k: v.dict() for k, v in (self.marker_results or {}).items()
            }
        }