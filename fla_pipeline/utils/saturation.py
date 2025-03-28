# fla_pipeline/utils/saturation.py

from scipy.optimize import curve_fit
import numpy as np

# Base Gaussian function
def gaussian(x, A, mu, sigma):
    return A * np.exp(-((x - mu) ** 2) / (2 * sigma ** 2))

# Fit a Gaussian to given data
def fit_gaussian_peak(x_data: np.ndarray, y_data: np.ndarray):
    if len(x_data) < 3:
        return None
    A0 = np.max(y_data)
    mu0 = x_data[np.argmax(y_data)]
    sigma0 = np.std(x_data) or 1.0
    try:
        popt, _ = curve_fit(gaussian, x_data, y_data, p0=[A0, mu0, sigma0])
        return popt  # A, mu, sigma
    except Exception:
        return None

# Robust saturation correction
def correct_if_saturated(idx: int, smap: np.ndarray, signal: np.ndarray, flank: int = 5, threshold_drop: float = 0.02):
    peak_val = signal[idx]
    threshold = peak_val * (1 - threshold_drop)

    # 1. Identify plateau region
    left = idx
    while left > 0 and signal[left] >= threshold:
        left -= 1
    right = idx
    while right < len(signal) and signal[right] >= threshold:
        right += 1

    # 2. Expand region around plateau for full-Gaussian fit
    fit_start = max(0, left - flank)
    fit_end = min(len(signal), right + flank)
    x_fit = smap[fit_start:fit_end]
    y_fit = signal[fit_start:fit_end]

    if len(x_fit) < 3:
        return peak_val, smap[idx], 0.0

    # 3. Fit full Gaussian
    fit = fit_gaussian_peak(x_fit, y_fit)
    if fit is not None:
        A, mu, sigma = fit
        fwhm = sigma * 2.355

        # 4. Accept only reasonable fits
        in_plateau_bounds = smap[left] - 0.5 <= mu <= smap[right] + 0.5
        if in_plateau_bounds and fwhm <= 2.5:
            A = min(A, peak_val * 1.2)  # limit correction boost
            return A, mu, fwhm

    # Fallback
    return float(peak_val), float(smap[idx]), 0.0
