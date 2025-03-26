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

            intensity = intensity_array[idx]
            saturated = intensity >= config.instrument_saturation_value

            if saturated:
                corrected = correct_if_saturated(idx, smap, intensity_array)
                note = "Saturated; intensity corrected via Gaussian fit."
                intensity = corrected
            else:
                note = ""

            peak = Peak(
                position=round(float(size), 2),
                intensity=round(float(intensity), 2),
                saturated=saturated,
                note=note
            )
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



    for shared_pos, peak_group in position_map.items():
        if len(peak_group) < 2:
            continue

        peak_group.sort(key=lambda x: x[1], reverse=True)
        dominant_ch, dom_intensity, dom_pk = peak_group[0]

        shared_rounded = round(dom_pk.position, 1)

        for ch, intensity, pk in peak_group[1:]:
            if round(pk.position, 1) != shared_rounded:
                continue

            ratio = intensity / dom_intensity if dom_intensity else 0
            if ratio >= 0.1:
                pk.note += f" Suppressed: pull-up from {dominant_ch}, ratio={ratio:.2f}."

                print(f"SUPPRESSING {ch} peak @ {pk.position:.3f} (rounded: {round(pk.position,1)}) due to pull-up from {dominant_ch}")

                # Remove peak by rounded position
                peak_dict[ch] = [
                    p for p in peak_dict[ch] if round(p.position, 1) != shared_rounded
                ]

    print("\n[DEBUG] Final LIZ Peaks:")
    for p in peak_dict.get("LIZ", []):
        print(f"  LIZ @ {p.position:.3f} | Intensity: {p.intensity:.1f} | Note: {p.note}")

    liz_peaks = peak_dict.get("LIZ", [])
    max_liz_intensity = max((p.intensity for p in liz_peaks), default=0.0)
    return peak_dict, max_liz_intensity