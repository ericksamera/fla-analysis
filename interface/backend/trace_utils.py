# interface/backend/trace_utils.py

import numpy as np

def find_peak_bounds(x, y, center_x, threshold, max_span=2.0, flat_tolerance=3):
    """
    Enhanced version to handle saturated plateaus by detecting flatness
    or distance-based cutoffs.
    """
    center_idx = np.argmin(np.abs(x - center_x))
    center_val = x[center_idx]

    # Left search
    left = center_idx
    flat_count = 0
    while left > 0:
        if y[left] <= threshold:
            break
        if abs(x[left] - center_val) > max_span:
            break
        if abs(y[left] - y[left - 1]) < 1e-2:
            flat_count += 1
            if flat_count >= flat_tolerance:
                break
        else:
            flat_count = 0
        left -= 1

    # Right search
    right = center_idx
    flat_count = 0
    while right < len(y) - 1:
        if y[right] <= threshold:
            break
        if abs(x[right] - center_val) > max_span:
            break
        if abs(y[right] - y[right + 1]) < 1e-2:
            flat_count += 1
            if flat_count >= flat_tolerance:
                break
        else:
            flat_count = 0
        right += 1

    return left, right

def find_suppressed_regions(smap, trace, suppressed_peaks, baseline_percentile=5):
    """
    Given an electropherogram trace and a list of suppressed peaks,
    return index ranges where the trace should be grayed out.

    Each suppressed peak is a tuple of (position, height).
    The gray region extends left/right until the trace returns to baseline.
    """
    regions = []
    baseline = np.percentile(trace, baseline_percentile)

    for pos, height in suppressed_peaks:
        if height <= 0:
            continue  # skip invalid data

        # Set a more robust threshold to detect base return
        threshold = max(baseline + 50, 0.05 * height)
        left, right = find_peak_bounds(smap, trace, pos, threshold)

        # Expand region slightly for visibility
        left = max(0, left - 1)
        right = min(len(trace) - 1, right + 1)

        regions.append((left, right))
    return regions
