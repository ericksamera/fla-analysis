import numpy as np
from scipy.optimize import curve_fit

def gaussian(x, A, mu, sigma):
    return A * np.exp(-((x - mu) ** 2) / (2 * sigma**2))

def fit_gaussian_peak(x_data: np.ndarray, y_data: np.ndarray):
    if len(x_data) < 3:
        return None
    A0 = np.max(y_data)
    mu0 = x_data[np.argmax(y_data)]
    sigma0 = np.std(x_data) or 1.0

    try:
        popt, _ = curve_fit(gaussian, x_data, y_data, p0=[A0, mu0, sigma0])
        return popt  # A, mu, sigma
    except Exception as e:
        #logger.debug(f"Gaussian fit failed: {e}")
        return None

def correct_if_saturated(idx: int, smap: np.ndarray, signal: np.ndarray, window: int = 5):
    start = max(0, idx - window)
    end = min(len(signal), idx + window + 1)
    x = smap[start:end]
    y = signal[start:end]
    fit = fit_gaussian_peak(x, y)
    if fit is not None:
        A, _, _ = fit
        return A
    return signal[idx]
