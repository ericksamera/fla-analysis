#!/usr/bin/env python3
__description__ =\
"""
core/data_utils.py

This module provides common utility functions for FLA processing,
including the standard Gaussian function, Gaussian fitting utilities,
and a function to aggregate quality control (QC) penalty scores.
"""
__author__ = "Erick Samera"
__version__ = "1.0.0"
__comments__ = "Works"

import numpy as np
from scipy.optimize import curve_fit
from typing import Optional, Tuple, List
import logging

logger = logging.getLogger(__name__)


def estimate_saturated_peak_intensity(window: int, data_np: np.ndarray, smap: np.ndarray, peak_idx: int) -> float:
    """
    For a peak suspected of saturation, fit a Gaussian in a local window around the peak
    to estimate the true amplitude. Return the fitted amplitude if it's higher.
    """
    # Prepare local window around the peak
    start = max(peak_idx - window, 0)
    end = min(peak_idx + window + 1, len(data_np))

    x_data = np.array(smap[start:end], dtype=float)
    y_data = np.array(data_np[start:end], dtype=float)

    fit_params = fit_gaussian_to_xy_data(x_data, y_data)
    if fit_params is not None:
        A, mu, sigma = fit_params
        return max(A, data_np[peak_idx])  # Return whichever is larger
    else:
        return data_np[peak_idx]

def gaussian(x: np.ndarray, A: float, mu: float, sigma: float) -> np.ndarray:
    """
    Standard Gaussian function.

    Parameters:
        x (np.ndarray): Input array of x values.
        A (float): Amplitude of the Gaussian.
        mu (float): Mean (center) of the Gaussian.
        sigma (float): Standard deviation (width) of the Gaussian.

    Returns:
        np.ndarray: The evaluated Gaussian function.
    """
    return A * np.exp(-((x - mu) ** 2) / (2 * sigma ** 2))

def fit_gaussian_to_xy_data(x_data: np.ndarray, y_data: np.ndarray) -> Optional[Tuple[float, float, float]]:
    """
    Fits a Gaussian to the provided x and y data.

    Parameters:
        x_data (np.ndarray): Array of x values.
        y_data (np.ndarray): Array of y values.

    Returns:
        Optional[Tuple[float, float, float]]: A tuple (A, mu, sigma) representing the 
        fitted Gaussian parameters, or None if the fit fails or there is insufficient data.
    """
    if len(x_data) < 3:
        return None

    A0 = np.max(y_data)
    mu0 = x_data[np.argmax(y_data)]
    sigma0 = np.std(x_data) if np.std(x_data) > 0 else 1.0

    try:
        popt, _ = curve_fit(gaussian, x_data, y_data, p0=[A0, mu0, sigma0])
        return popt[0], popt[1], popt[2]
    except Exception as e:
        logger.debug(f"Gaussian fit failed: {e}")
        return None

def aggregate_qc_penalties(penalties: List[float]) -> float:
    """
    Aggregates quality control penalty values by summing them up.

    Parameters:
        penalties (List[float]): A list of penalty scores.

    Returns:
        float: The sum of all penalty scores.
    """
    return sum(penalties)