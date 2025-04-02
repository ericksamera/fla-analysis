from fla_pipeline.models.peak import Peak
from fla_pipeline.config.global_config import GlobalConfig
from fla_pipeline.config.marker_config import MarkerConfig
from typing import List, Tuple, Dict
from collections import defaultdict
import numpy as np


def bin_peaks(peaks: List[Peak], marker: MarkerConfig, config: GlobalConfig) -> Tuple[List[Peak], List[str]]:
    bmin, bmax = marker.bins
    repeat = marker.repeat_unit or 1
    tolerance = repeat * config.bin_tolerance
    qc_flags = []

    # Step 1: Define bins with extended range
    extended_min = bmin - repeat * config.bin_extend
    extended_max = bmax + repeat * config.bin_extend
    bins = np.arange(extended_min, extended_max + 1, repeat)

    # Step 2: Assign each peak to its closest bin (within tolerance)
    bin_to_peaks: Dict[float, List[Peak]] = defaultdict(list)
    for peak in peaks:
        closest_bin = min(bins, key=lambda b: abs(peak.position - b))
        if abs(peak.position - closest_bin) <= tolerance:
            bin_to_peaks[closest_bin].append(peak)

    if not bin_to_peaks:
        return [], [f"No peaks assigned to bins for {marker.marker}."]

    # Step 3: Choose primary bin = one with max total intensity
    anchor_bin = max(bin_to_peaks.items(), key=lambda kv: sum(p.intensity for p in kv[1]))[0]
    anchor_idx = np.where(bins == anchor_bin)[0][0]

    # Step 4: Re-assign peaks based on distance from anchor (relative repeat units)
    reassigned_bins: Dict[float, List[Peak]] = defaultdict(list)
    for peak in peaks:
        offset = round((peak.position - anchor_bin) / repeat)
        rel_bin = anchor_bin + offset * repeat
        if extended_min <= rel_bin <= extended_max and abs(peak.position - rel_bin) <= tolerance:
            reassigned_bins[rel_bin].append(peak)

    # Step 5: From each bin, keep strongest peak, avoiding 1bp overlap
    assigned_peaks = []
    for bin_pos, bin_peaks in reassigned_bins.items():
        bin_peaks = sorted(bin_peaks, key=lambda p: p.intensity, reverse=True)
        filtered = []
        for i, pk in enumerate(bin_peaks):
            too_close = any(
                abs(pk.position - other.position) < 1.0 and pk.intensity < other.intensity
                for other in bin_peaks[:i]
            )
            if not too_close:
                filtered.append(pk)
        if filtered:
            top = max(filtered, key=lambda p: p.intensity)
            assigned_peaks.append(Peak(
                position=bin_pos,
                intensity=top.intensity,
                saturated=top.saturated,
                note=top.note + f" (snapped from {top.position:.2f})"
            ))

    if not assigned_peaks:
        return [], [f"All peaks filtered due to proximity or weakness for {marker.marker}."]

    # Step 6: Intensity threshold relative to strongest retained
    max_intensity = max(p.intensity for p in assigned_peaks)
    threshold = config.relative_peak_threshold * max_intensity
    final_peaks = [p for p in assigned_peaks if p.intensity >= threshold]
    n_removed = len(assigned_peaks) - len(final_peaks)

    if n_removed > 0:
        qc_flags.append(
            f"{n_removed} peak(s) below {int(config.relative_peak_threshold * 100)}% of retained max intensity were removed."
        )

    return final_peaks, qc_flags
