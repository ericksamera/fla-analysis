# fla_pipeline/analysis/metrics/bruvo.py

from fla_pipeline.analysis.metrics.base import BaseDistanceMetric
from fla_pipeline.models.genotype import GenotypeResult
from fla_pipeline.config.marker_config import MarkerConfig
from typing import Optional
import itertools
import numpy as np

class BruvoDistanceMetric(BaseDistanceMetric):
    name = "bruvo"

    def compute(
        self,
        geno1: GenotypeResult,
        geno2: GenotypeResult,
        marker_cfg: MarkerConfig
    ) -> Optional[float]:
        alleles1 = sorted(geno1.alleles)
        alleles2 = sorted(geno2.alleles)

        if not alleles1 or not alleles2:
            self.log(f"Missing genotype(s) at marker {marker_cfg.marker}")
            return None

        repeat_unit = marker_cfg.repeat_unit or 1
        g1_counts = [round(a / repeat_unit) for a in alleles1]
        g2_counts = [round(a / repeat_unit) for a in alleles2]

        score = self._bruvo_distance(g1_counts, g2_counts)
        self.log(f"{marker_cfg.marker}: Bruvo({g1_counts} vs {g2_counts}) = {score:.4f}")
        return score

    def _bruvo_distance(self, g1, g2) -> float:
        n1, n2 = len(g1), len(g2)
        best = np.inf

        if n1 < n2:
            expansions = itertools.product(g1, repeat=(n2 - n1))
            for exp in expansions:
                full = list(g1) + list(exp)
                best = min(best, self._min_perm_distance(full, g2))
        elif n2 < n1:
            expansions = itertools.product(g2, repeat=(n1 - n2))
            for exp in expansions:
                full = list(g2) + list(exp)
                best = min(best, self._min_perm_distance(g1, full))
        else:
            best = self._min_perm_distance(g1, g2)

        return best if best != np.inf else 1.0

    def _min_perm_distance(self, g1, g2):
        best = np.inf
        for perm in itertools.permutations(g2):
            distances = [1 - 2 ** (-abs(a - b)) for a, b in zip(g1, perm)]
            avg = sum(distances) / len(distances) if distances else 1.0
            best = min(best, avg)
        return best
