# fla_pipeline/core/peak_binner.py

from fla_pipeline.models.peak import Peak
from fla_pipeline.config import GlobalConfig, MarkerConfig
from typing import List, Tuple, Dict
from collections import defaultdict
import numpy as np


def bin_peaks(peaks: List[Peak], marker: MarkerConfig, config: GlobalConfig) -> Tuple[List[Peak], List[str]]:
    bmin, bmax = marker.bins
    repeat = marker.repeat_unit or 1
    tolerance = repeat * config.bin_tolerance
    qc_flags = []

    # Step 1: Generate bins
    extended_min = bmin - repeat * config.bin_extend
    extended_max = bmax + repeat * config.bin_extend
    bins = np.arange(extended_min, extended_max + 1, repeat)

    # Step 2: Assign peaks to closest bin if within tolerance
    bin_to_peaks: Dict[float, List[Peak]] = defaultdict(list)
    for peak in peaks:
        closest_bin = min(bins, key=lambda b: abs(peak.position - b))
        if abs(peak.position - closest_bin) <= tolerance:
            bin_to_peaks[closest_bin].append(peak)

    if not bin_to_peaks:
        return [], [f"No peaks assigned to bins for {marker.marker}."]

    # Step 3: Within each bin, keep strongest peak
    assigned_peaks = []
    for bin_pos, bin_peaks in bin_to_peaks.items():
        if not bin_peaks:
            continue
        # Exclude peaks that are <1bp from another stronger peak
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
            strongest = max(filtered, key=lambda p: p.intensity)
            # Snap the peak position to the bin it's assigned to
            assigned_peaks.append(Peak(
                position=bin_pos,  # << use the bin, not the raw peak
                intensity=strongest.intensity,
                saturated=strongest.saturated,
                note=strongest.note + f" (snapped from {strongest.position:.2f})"
            ))

    if not assigned_peaks:
        return [], [f"All peaks filtered due to proximity or weakness for {marker.marker}."]

    # Step 4: Intensity threshold (relative to max in retained peaks)
    max_intensity = max(p.intensity for p in assigned_peaks)
    threshold = config.relative_peak_threshold * max_intensity
    final_peaks = [p for p in assigned_peaks if p.intensity >= threshold]
    n_removed = len(assigned_peaks) - len(final_peaks)

    if n_removed > 0:
        qc_flags.append(
            f"{n_removed} peak(s) below {int(config.relative_peak_threshold * 100)}% of retained max intensity were removed."
        )

    return final_peaks, qc_flags