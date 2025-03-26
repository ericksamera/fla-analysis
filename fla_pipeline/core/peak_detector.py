# fla_pipeline/core/peak_detector.py

from collections import defaultdict
from fla_pipeline.config import GlobalConfig
from fla_pipeline.models.peak import Peak
from fla_pipeline.utils.saturation import correct_if_saturated
from typing import Dict, List, Tuple
import numpy as np
from scipy.signal import find_peaks

def detect_peaks(
    smap: np.ndarray,
    channels: Dict[str, np.ndarray],
    config: GlobalConfig
) -> Tuple[Dict[str, List[Peak]], float]:
    peak_dict: Dict[str, List[Peak]] = {}

    # Step 1: Detect peaks per channel
    for ch_name, intensity_array in channels.items():
        if ch_name == "LIZ":
            peak_indices, _ = find_peaks(intensity_array, height=100, distance=3)
        else:
            peak_indices, _ = find_peaks(
                intensity_array,
                height=config.min_peak_height,
                distance=3
            )

        peaks: List[Peak] = []
        for idx in peak_indices:
            size = smap[idx]
            if size < config.min_peak_position:
                continue

            raw_intensity = intensity_array[idx]
            saturated = raw_intensity >= config.instrument_saturation_value

            if saturated:
                corrected_intensity, corrected_position, width = correct_if_saturated(idx, smap, intensity_array)
                note = f"Saturated; corrected via Gaussian fit (μ={corrected_position:.2f}, FWHM={width:.2f})"
            else:
                corrected_intensity = raw_intensity
                corrected_position = size
                width = 0.0
                note = ""

            peak = Peak(
                position=round(corrected_position, 2),
                intensity=round(corrected_intensity, 2),
                saturated=saturated,
                note=note
            )
            peak.corrected_intensity = round(corrected_intensity, 2)
            peak.width = width

            peaks.append(peak)


            peak_dict[ch_name] = peaks


    # Step 2: Cross-channel artifact detection
    all_peaks = []
    for channel, peaks in peak_dict.items():
        for pk in peaks:
            all_peaks.append((pk.position, pk.intensity, channel, pk))

    # Group by rounded position
    position_map = defaultdict(list)
    for pos, intensity, ch, pk in all_peaks:
        key = round(pos, 1)  # or 0.2 for tighter tolerance
        position_map[key].append((ch, intensity, pk))


    suppressed_peaks = defaultdict(list)
    for shared_pos, peak_group in position_map.items():
        if len(peak_group) < 2:
            continue

        non_liz_peaks = [item for item in peak_group if item[0] != "LIZ"]

        if non_liz_peaks:
            candidates = sorted(
                non_liz_peaks,
                key=lambda x: getattr(x[2], "corrected_intensity", x[1]),
                reverse=True
            )
            dominant_ch, dom_intensity, dom_pk = candidates[0]

            # Only proceed if dominant peak is saturated or very wide (flat)
            if not dom_pk.saturated and dom_pk.width < 0.2:
                continue
        else:
            continue  # skip suppression entirely if only L




        shared_rounded = round(dom_pk.position, 1)

        for ch, intensity, pk in peak_group[1:]:
            if round(pk.position, 1) != shared_rounded:
                continue

            ratio = intensity / dom_intensity if dom_intensity else 0
            if ratio >= 0.1:
                pk.note += f" Suppressed: pull-up from {dominant_ch}, ratio={ratio:.2f}."
                suppressed_peaks[ch].append((pk.position, intensity))

                peak_dict[ch] = [
                    p for p in peak_dict[ch] if round(p.position, 1) != shared_rounded
                ]

    liz_peaks = peak_dict.get("LIZ", [])
    max_liz_intensity = max((p.intensity for p in liz_peaks), default=0.0)
    return peak_dict, max_liz_intensity, suppressed_peaks
