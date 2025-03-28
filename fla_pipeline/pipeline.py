# fla_pipeline/pipeline.py

from fla_pipeline.config.global_config import GlobalConfig
from fla_pipeline.config.marker_config import MarkerConfig
from fla_pipeline.core.fsa_loader import load_fsa
from fla_pipeline.core.peak_detector import detect_peaks
from fla_pipeline.core.peak_binner import bin_peaks
from fla_pipeline.strategies.diploid import get_diploid_caller
from fla_pipeline.models.genotype import GenotypeResult
from typing import Dict, Any

from fla_pipeline.models.sample import Sample 

def run_pipeline(
    file_path: str,
    config: GlobalConfig,
    markers: list[MarkerConfig],
    sample: Sample | None = None
) -> Dict[str, Any]:
    # Use existing data if available
    if sample and sample.fsa_data and sample.peaks:
        fsa_data = sample.fsa_data
        peak_dict = sample.peaks
        max_liz_intensity = sample.metadata.get("max_liz_intensity", 0.0)
        suppressed_peaks = sample.suppressed_peaks or {}
    else:
        fsa_data = load_fsa(file_path)
        smap = fsa_data["smap"]
        channels = fsa_data["channels"]
        peak_dict, max_liz_intensity, suppressed_peaks = detect_peaks(smap, channels, config)

        # If a Sample object was passed, populate it
        if sample is not None:
            sample.fsa_data = fsa_data
            sample.peaks = peak_dict
            sample.suppressed_peaks = suppressed_peaks
            sample.metadata["max_liz_intensity"] = max_liz_intensity

    smap = fsa_data["smap"]
    marker_results: Dict[str, dict] = {}

    if config.ploidy != 2:
        raise NotImplementedError("Only diploid mode is implemented.")

    for marker in markers:
        peaks = peak_dict.get(marker.channel, [])
        per_marker_cfg = GlobalConfig(**marker.merged_with_global(config))

        if marker.binning_enabled:
            binned_peaks, bin_flags = bin_peaks(peaks, marker, per_marker_cfg)
        else:
            binned_peaks = peaks
            bin_flags = ["Binning disabled; using raw peak positions."]

        if per_marker_cfg.ploidy == 2:
            caller = get_diploid_caller(per_marker_cfg.diploid_strategy)
        elif per_marker_cfg.ploidy > 2:
            raise NotImplementedError("Polyploid genotype calling is not yet supported.")
        else:
            raise ValueError("Invalid ploidy setting.")

        genotype: GenotypeResult = caller.call_genotype(binned_peaks, marker, per_marker_cfg, max_liz_intensity=max_liz_intensity)
        genotype.qc_flags = bin_flags + genotype.qc_flags

        if genotype.confidence < per_marker_cfg.min_genotype_confidence:
            genotype.alleles = []
            genotype.qc_flags.append(
                f"Genotype confidence {genotype.confidence:.2f} below threshold ({per_marker_cfg.min_genotype_confidence}). Call suppressed."
            )

        marker_results[marker.marker] = genotype.to_dict()

    return {
        "fsa_data": {
            "name": fsa_data["name"],
            "smap": smap.tolist(),
        },
        "detected_peaks": {
            ch: [vars(p) for p in pk_list] for ch, pk_list in peak_dict.items()
        },
        "marker_results": marker_results,
        "max_liz_intensity": max_liz_intensity,
        "suppressed_peaks": suppressed_peaks
    }
