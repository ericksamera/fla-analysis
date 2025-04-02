from fla_pipeline.analysis.metrics.base import BaseDistanceMetric
from fla_pipeline.models.genotype import GenotypeResult
from fla_pipeline.config.marker_config import MarkerConfig
from collections import Counter
from typing import Optional
import numpy as np

class NeiDistanceMetric(BaseDistanceMetric):
    name = "nei"

    def compute(self, geno1: GenotypeResult, geno2: GenotypeResult, marker_cfg: MarkerConfig) -> Optional[float]:
        if not geno1.alleles or not geno2.alleles:
            self.log(f"Missing alleles for marker {marker_cfg.marker}")
            return None

        # Count alleles
        c1 = Counter(map(str, geno1.alleles))
        c2 = Counter(map(str, geno2.alleles))

        total1 = sum(c1.values())
        total2 = sum(c2.values())
        if total1 == 0 or total2 == 0:
            return None

        freqs1 = {a: c / total1 for a, c in c1.items()}
        freqs2 = {a: c / total2 for a, c in c2.items()}

        shared_alleles = set(freqs1) | set(freqs2)
        I = sum(np.sqrt(freqs1.get(a, 0) * freqs2.get(a, 0)) for a in shared_alleles)
        if I <= 0:
            return 1.0

        D = -np.log(I)
        self.log(f"{marker_cfg.marker}: Nei distance = {D:.4f}")
        return D
