# fla_pipeline/analysis/distance.py

from fla_pipeline.analysis.config import DistanceConfig
from fla_pipeline.analysis.metrics.bruvo import BruvoDistanceMetric
from fla_pipeline.analysis.metrics.base import BaseDistanceMetric
from fla_pipeline.models.genotype import GenotypeResult
from fla_pipeline.config.marker_config import MarkerConfig
from typing import Dict, List, Optional
import pandas as pd
import numpy as np
import itertools

# Metric registry
DISTANCE_REGISTRY = {
    "bruvo": BruvoDistanceMetric,
    # Future: "nei": NeiDistanceMetric,
    # etc.
}

class DistanceCalculator:
    def __init__(
        self,
        genotype_matrix: Dict[str, Dict[str, GenotypeResult]],
        marker_configs: Dict[str, MarkerConfig],
        sample_metadata: Optional[Dict[str, Dict]] = None
    ):
        self.genotypes = genotype_matrix  # sample -> marker -> GenotypeResult
        self.marker_configs = marker_configs
        self.metadata = sample_metadata or {}
        self.samples = sorted(genotype_matrix.keys())
        self.debug_logs = []
        self.included_markers = []

    def compute_matrix(self, config: DistanceConfig) -> pd.DataFrame:
        metric_cls = DISTANCE_REGISTRY.get(config.metric)
        if not metric_cls:
            raise ValueError(f"Unknown distance metric: {config.metric}")
        metric: BaseDistanceMetric = metric_cls(debug=config.debug_mode)

        dist_matrix = pd.DataFrame(np.nan, index=self.samples, columns=self.samples)
        self.included_markers.clear()

        for i, s1 in enumerate(self.samples):
            for j, s2 in enumerate(self.samples):
                if j < i:
                    continue

                markers = self._shared_markers(s1, s2, config)

                if not markers:
                    dist_matrix.at[s1, s2] = 1.0
                    dist_matrix.at[s2, s1] = 1.0
                    continue

                distances = []
                for marker in markers:
                    geno1 = self.genotypes[s1].get(marker)
                    geno2 = self.genotypes[s2].get(marker)

                    result = metric.compute(geno1, geno2, self.marker_configs[marker])
                    if result is not None:
                        distances.append(result)

                avg = sum(distances) / len(distances) if distances else 1.0
                dist_matrix.at[s1, s2] = avg
                dist_matrix.at[s2, s1] = avg

        self.debug_logs = metric.get_logs()
        return dist_matrix

    def _shared_markers(self, s1, s2, config: DistanceConfig) -> List[str]:
        g1 = self.genotypes[s1]
        g2 = self.genotypes[s2]

        markers = set(g1.keys()) & set(g2.keys())
        filtered = []

        for marker in markers:
            geno1 = g1.get(marker)
            geno2 = g2.get(marker)

            if not geno1 or not geno2:
                continue
            if geno1.confidence < config.min_confidence:
                continue
            if geno2.confidence < config.min_confidence:
                continue

            sample_count = sum(
                1 for g in self.genotypes.values() if marker in g and g[marker].confidence >= config.min_confidence
            )
            if sample_count < config.min_sample_count:
                continue

            filtered.append(marker)

        # Save for UI preview
        if s1 == self.samples[0] and s2 == self.samples[1]:  # only once
            self.included_markers = filtered

        return filtered

    def perform_pcoa(self, distance_matrix: pd.DataFrame, n_components: int = 2):
        D = distance_matrix.values.astype(float)
        n = D.shape[0]

        if n < n_components:
            raise ValueError(f"Not enough samples for {n_components} components")

        D2 = D ** 2
        J = np.eye(n) - np.ones((n, n)) / n
        B = -0.5 * J @ D2 @ J
        eigvals, eigvecs = np.linalg.eigh(B)
        idx = np.argsort(eigvals)[::-1]
        eigvals = eigvals[idx]
        eigvecs = eigvecs[:, idx]
        pos_idx = eigvals > 0
        coords = eigvecs[:, pos_idx] * np.sqrt(eigvals[pos_idx])
        explained = eigvals / eigvals.sum()

        return coords[:, :n_components], explained[:n_components]

    def get_debug_log(self) -> str:
        return "\n".join(self.debug_logs)

    def get_included_markers(self) -> List[str]:
        return self.included_markers
