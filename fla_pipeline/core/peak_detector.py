# fla_pipeline/core/peak_detector.py

from collections import defaultdict
from fla_pipeline.config import GlobalConfig
from fla_pipeline.models.peak import Peak
from typing import Dict, List
import numpy as np
from scipy.signal import find_peaks

def detect_peaks(
    smap: np.ndarray,
    channels: Dict[str, np.ndarray],
    config: GlobalConfig
) -> Dict[str, List[Peak]]:
    peak_dict: Dict[str, List[Peak]] = {}

    # Step 1: Detect peaks per channel
    for ch_name, intensity_array in channels.items():
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
            peak = Peak(
                position=round(float(size), 2),
                intensity=round(float(intensity_array[idx]), 2),
                saturated=intensity_array[idx] >= config.instrument_saturation_value,
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
            continue  # not across multiple channels

        # Find brightest channel
        peak_group.sort(key=lambda x: x[1], reverse=True)
        dominant_channel = peak_group[0][0]

        for ch, intensity, pk in peak_group[1:]:
            ratio = intensity / peak_group[0][1] if peak_group[0][1] else 0
            if ratio >= 0.1:  # Adjustable threshold
                pk.note += f" Possible pull-up artifact (from {dominant_channel}, ratio={ratio:.2f})."

    return peak_dict

