# fla_pipeline/strategies/diploid.py

from fla_pipeline.models.peak import Peak
from fla_pipeline.models.genotype import GenotypeResult
from fla_pipeline.config import GlobalConfig, MarkerConfig
from fla_pipeline.utils.evaluator import GenotypeEvaluator
from fla_pipeline.utils.qc_score import QCScore
from typing import List


DIPLOID_STRATEGY_REGISTRY = {}

def register_diploid_strategy(name: str):
    def decorator(cls):
        DIPLOID_STRATEGY_REGISTRY[name] = cls
        return cls
    return decorator

def get_diploid_caller(strategy_name: str):
    caller_cls = DIPLOID_STRATEGY_REGISTRY.get(strategy_name)
    if not caller_cls:
        raise ValueError(f"Unknown diploid strategy: {strategy_name}")
    return caller_cls()

class BaseDiploidCaller:
    def fallback_homozygote(self, marker: MarkerConfig, peaks: List[Peak], strategy: str, note: str = "") -> GenotypeResult:
        pos = round(peaks[0].position)
        return GenotypeResult(
            marker=marker.marker,
            alleles=[pos, pos],
            confidence=0.95,
            strategy=strategy,
            qc_flags=["Only one peak; homozygous call." + (" " + note if note else "")]
        )


@register_diploid_strategy("probabilistic")
class ProbabilisticCaller(BaseDiploidCaller):

    def call_genotype(self, peaks: List[Peak], marker: MarkerConfig, config: GlobalConfig) -> GenotypeResult:
        if not peaks:
            return GenotypeResult(marker=marker.marker, alleles=[], strategy="probabilistic", qc_flags=["No peaks to call."])

        peaks = sorted(peaks, key=lambda p: p.intensity, reverse=True)

        if len(peaks) == 1:
            return self.fallback_homozygote(marker, peaks, "probabilistic")

        evaluator = GenotypeEvaluator(config, marker)
        top_n = min(10, len(peaks))

        results = []
        for j in range(1, top_n):
            pk1, pk2 = peaks[0], peaks[j]
            if abs(pk1.position - pk2.position) < 0.25:
                continue

            result = evaluator.evaluate_pair(pk1, pk2, peaks)
            results.append(result)

        results.sort(key=lambda r: r.confidence, reverse=True)

        best = results[0] if results else None
        second = results[1] if len(results) > 1 else None

        if best and second and abs(best.confidence - second.confidence) < 0.075:
            best.confidence = round(best.confidence - 0.25, 3)
            best.qc_flags.append("Alternate peak pairing nearly as likely - ambiguous call. (-0.25)")

        return best or self.fallback_homozygote(marker, peaks, "probabilistic", note="No confident pairing found.")




@register_diploid_strategy("strict_stutter")
class StrictStutterCaller(BaseDiploidCaller):
    def call_genotype(self, peaks: List[Peak], marker: MarkerConfig, config: GlobalConfig) -> GenotypeResult:
        if not peaks:
            return GenotypeResult(marker=marker.marker, alleles=[], strategy="strict_stutter", qc_flags=["No peaks to call."])

        peaks = sorted(peaks, key=lambda p: p.intensity, reverse=True)
        major = peaks[0]
        if len(peaks) == 1:
            return self.fallback_homozygote(marker, peaks, "strict_stutter")

        minor = peaks[1]
        ratio = minor.intensity / major.intensity
        distance = abs(major.position - minor.position)
        is_stutter_like = (
            marker.repeat_unit
            and distance == marker.repeat_unit
            and ratio <= config.stutter_ratio_threshold
        )

        if is_stutter_like:
            return GenotypeResult(
                marker=marker.marker,
                alleles=[round(major.position)] * 2,
                confidence=0.90,
                strategy="strict_stutter",
                qc_flags=[f"Minor peak ({minor.position:.1f}) may be stutter; homozygous call."]
            )

        return GenotypeResult(
            marker=marker.marker,
            alleles=sorted([round(major.position), round(minor.position)]),
            confidence=0.80,
            strategy="strict_stutter",
            qc_flags=[f"Minor peak included (ratio = {ratio:.2f})."]
        )

@register_diploid_strategy("lenient_ratio")
class LenientRatioCaller(BaseDiploidCaller):
    def call_genotype(self, peaks: List[Peak], marker: MarkerConfig, config: GlobalConfig) -> GenotypeResult:
        if not peaks:
            return GenotypeResult(marker=marker.marker, alleles=[], strategy="lenient_ratio", qc_flags=["No peaks to call."])

        peaks = sorted(peaks, key=lambda p: p.intensity, reverse=True)
        major = peaks[0]
        if len(peaks) == 1:
            return self.fallback_homozygote(marker, peaks, "lenient_ratio")

        minor = peaks[1]
        ratio = minor.intensity / major.intensity
        threshold = marker.het_ratio_lower_bound or config.relative_peak_threshold

        if ratio >= threshold:
            return GenotypeResult(
                marker=marker.marker,
                alleles=sorted([round(major.position), round(minor.position)]),
                confidence=0.85,
                strategy="lenient_ratio",
                qc_flags=[f"Calling heterozygote (ratio = {ratio:.2f})."]
            )
        else:
            return GenotypeResult(
                marker=marker.marker,
                alleles=[round(major.position)] * 2,
                confidence=0.80,
                strategy="lenient_ratio",
                qc_flags=[f"Minor peak below threshold (ratio = {ratio:.2f}); homozygous call."]
            )
