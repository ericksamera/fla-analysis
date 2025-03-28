# fla_pipline/utils/evaluator.py

from typing import List, Tuple
from fla_pipeline.models.peak import Peak
from fla_pipeline.models.genotype import GenotypeResult
from fla_pipeline.config.global_config import GlobalConfig
from fla_pipeline.config.marker_config import MarkerConfig
from fla_pipeline.utils.qc_score import QCScore


class GenotypeEvaluator:
    def __init__(self, config: GlobalConfig, marker: MarkerConfig, max_liz_intensity: float = 0.0):
        self.config = config
        self.marker = marker
        self.qc = QCScore(config)
        self.qc_flags: List[str] = []
        self.max_liz_intensity = max_liz_intensity

    def evaluate_pair(self, pk1: Peak, pk2: Peak, all_peaks: List[Peak]) -> GenotypeResult:
        pos1, pos2 = round(pk1.position), round(pk2.position)
        alleles = sorted([pos1, pos2])
        ratio = min(pk1.intensity, pk2.intensity) / max(pk1.intensity, pk2.intensity)
        distance = abs(pk1.position - pk2.position)

        stutter_score, stutter_note = self.stutter_likelihood(distance, ratio)
        self.qc.add(1.0 - stutter_score, "Stutter ambiguity", weight=self.config.stutter_weight)
        self.qc_flags.append(stutter_note)

        if self.is_probable_stutter(pk2, all_peaks):
            self.qc.add(0.50, "Minor peak resembles stutter ladder", weight=self.config.stutter_weight)
            self.qc_flags.append("Minor peak may be part of stutter ladder; suppressing heterozygote call.")

        penalty = self.detect_distracting_peaks(all_peaks, pk1, pk2)
        if penalty > 0:
            self.qc.add(penalty, "Distracting secondary peaks", weight=self.config.ratio_weight)

        confidence = min(self.qc.compute_confidence(), 0.95)

        if self.max_liz_intensity:
            allele_peak_min = min(pk1.intensity, pk2.intensity)
            liz_ratio = allele_peak_min / self.max_liz_intensity
            threshold_ratio = self.config.liz_peak_penalty_ratio

            if liz_ratio < threshold_ratio:
                penalty = (threshold_ratio - liz_ratio) * 0.25  # scale the penalty
                confidence = max(0.0, confidence - penalty)
                self.qc_flags.append(f"Low relative to LIZ: {liz_ratio:.2f} < {threshold_ratio} (-{penalty:.2f})")


        return GenotypeResult(
            marker=self.marker.marker,
            alleles=alleles,
            confidence=confidence,
            strategy="probabilistic",
            qc_flags=self.qc_flags + self.qc.summary_flags()
        )

    def stutter_likelihood(self, distance: float, ratio: float) -> Tuple[float, str]:
        if not self.marker.repeat_unit:
            return 0.5, "Repeat unit unknown; default 0.5 likelihood."

        expected_dist = self.marker.repeat_unit
        tolerance = self.config.stutter_tolerance or 1.0
        ratio_thresh = self.config.stutter_ratio_threshold or 0.15
        decay = self.config.stutter_step_decay or 0.5

        repeat_steps = [expected_dist * i for i in range(1, 10)]
        penalties = []

        for i, step in enumerate(repeat_steps):
            closeness = max(0.0, 1.0 - abs(distance - step) / tolerance)
            decay_weight = decay ** i
            penalties.append(closeness * decay_weight)

        total_closeness = max(penalties, default=0.0)

        ratio_penalty = max(0.0, 1.0 - min(ratio / ratio_thresh, 1.0))
        adjusted_penalty = total_closeness * ratio_penalty

        stutter_score = 1.0 - adjusted_penalty

        note = f"Minor peak distance: {distance:.2f}, ratio: {ratio:.2f}, stutter_score: {stutter_score:.2f}"
        return stutter_score, note

    def is_probable_stutter(self, candidate: Peak, peaks: List[Peak]) -> bool:
        repeat_unit = self.marker.repeat_unit or 1.0
        tolerance = self.config.stutter_tolerance
        intensity_threshold = self.config.stutter_ratio_threshold
        direction = self.marker.stutter_direction or "both"

        for parent in peaks:
            if parent.intensity <= candidate.intensity:
                continue

            offset = parent.position - candidate.position
            if direction == "left" and offset <= 0:
                continue
            elif direction == "right" and offset >= 0:
                continue
            elif direction not in {"left", "right", "both"}:
                continue

            if abs(offset - repeat_unit) <= tolerance:
                ratio = candidate.intensity / parent.intensity
                if ratio > intensity_threshold:
                    return False
                return True
        return False

    def detect_distracting_peaks(self, peaks: List[Peak], major: Peak, minor: Peak) -> float:
        penalty = 0.0
        for p in peaks[2:]:
            if abs(p.position - major.position) < 1 or abs(p.position - minor.position) < 1:
                continue

            if self.is_probable_stutter(p, peaks):
                continue

            base = major.intensity if abs(p.intensity - major.intensity) < abs(p.intensity - minor.intensity) else minor.intensity
            ratio = p.intensity / base if base > 0 else 0

            if ratio < 0.5:
                continue

            scaled_penalty = min((ratio - 0.5) / 0.6, 1.0) * 0.2
            penalty += scaled_penalty

        return min(penalty, 1.0)
