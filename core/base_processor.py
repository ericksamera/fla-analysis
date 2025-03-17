#!/usr/bin/env python3
__description__ =\
"""
core/base_processor.py

This module provides a base class for FLA processing. It encapsulates core business logic
for loading FSA data, detecting peaks, and performing common computations. This module is
completely decoupled from any UI or presentation logic, making it easier to test,
debug, and reuse in different contexts.
"""
__author__ = "Erick Samera"
__version__ = "1.0.0"
__comments__ = "Works"

import logging
import numpy as np
from dataclasses import dataclass
from typing import List, Dict, Optional, Any
from Bio import SeqIO
from scipy.signal import find_peaks

from core.data_utils import estimate_saturated_peak_intensity

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


@dataclass
class BaseConfig:
    """
    Base configuration parameters for FLA processing.
    These parameters are shared across different processor implementations.
    """
    bin_tolerance: int = 2
    min_peak_height: int = 50
    min_peak_position: int = 50
    relative_peak_threshold: float = 0.1
    ratio_weight: float = 1.0
    stutter_weight: float = 1.0
    intensity_weight: float = 3.0
    highest_peak_ratio_weight: float = 1.0
    qc_confidence_scale: float = 2.0
    stutter_ratio_threshold: float = 0.15
    stutter_tolerance: float = 1.0
    instrument_saturation_value: float = 30000.0
    global_het_ratio_lower_bound: float = 0.35
    gaussian_fit_window: int = 5
    ploidy: int = 2


@dataclass
class Peak:
    """
    Represents a detected peak in the FSA data.
    """
    position: float
    intensity: float
    saturated: bool = False
    note: str = ""


class BaseProcessor:
    """
    Abstract base class for FLA processors.
    Provides common functionality for:
      - Loading FSA data
      - Detecting peaks
      - Running the overall analysis pipeline

    Subclasses should override or extend methods as needed for specialized behavior.
    """

    def __init__(self, config: BaseConfig):
        """
        Initialize the processor with the given configuration.

        :param config: BaseConfig instance with processing parameters.
        """
        self.config = config
        self.fsa_data: Optional[Dict[str, Any]] = None
        self.detected_peaks: Dict[str, List[Peak]] = {}

    def load_fsa(self, file_path) -> bool:
        """Loads the .fsa file into memory, extracting SMAP and channel data arrays."""
        try:
            fsa_object = SeqIO.read(file_path, 'abi')
            smap = fsa_object.annotations.get('abif_raw', {}).get('SMap2', [])
            channels = {
                "6-FAM": fsa_object.annotations['abif_raw'].get('DATA9', []),
                "VIC":   fsa_object.annotations['abif_raw'].get('DATA10', []),
                "NED":   fsa_object.annotations['abif_raw'].get('DATA11', []),
                "PET":   fsa_object.annotations['abif_raw'].get('DATA12', []),
                "LIZ":   fsa_object.annotations['abif_raw'].get('DATA205', []),
            }
            self.fsa_data = {
                "name": fsa_object.name,
                "smap": smap,
                "channels": channels
            }
            return True
        except Exception as e:
            logger.error(f"Error loading FSA file: {e}")
            return False

    def detect_peaks(self) -> None:
        """
        Detects peaks in each channel using SciPy's find_peaks, with saturation handling.
        Stores them in self.detected_peaks.
        """
        if self.fsa_data is None:
            logger.error("FSA data not loaded.")
            return

        smap = np.array(self.fsa_data["smap"], dtype=float)
        self.detected_peaks = {}

        for channel, data in self.fsa_data["channels"].items():
            data_np = np.array(data, dtype=float)
            indices = find_peaks(data_np)[0]

            peaks: List[Peak] = []
            for i in indices:
                raw_intensity = data_np[i]
                if raw_intensity >= self.config.min_peak_height:
                    pos = smap[i] if i < len(smap) else None
                    if pos is not None and pos >= self.config.min_peak_position:
                        saturated = (raw_intensity >= self.config.instrument_saturation_value)
                        if saturated:
                            corrected_intensity = estimate_saturated_peak_intensity(
                                self.config.gaussian_fit_window,
                                data_np,
                                smap,
                                i)
                            peak = Peak(position=float(pos),
                                        intensity=float(corrected_intensity),
                                        saturated=True,
                                        note="Saturated peak; corrected via Gaussian.")
                        else:
                            peak = Peak(position=float(pos), intensity=float(raw_intensity))
                        peaks.append(peak)

            peaks.sort(key=lambda p: p.intensity, reverse=True)
            self.detected_peaks[channel] = peaks

    def process_all(self, file_path: str) -> Optional[Dict[str, Any]]:
        """
        High-level processing routine: loads FSA data, detects peaks,
        and returns the processed results.

        :param file_path: Path to the FSA file.
        :return: Dictionary containing FSA data and detected peaks, or None if processing fails.
        """
        if not self.load_fsa(file_path):
            logger.error("Failed to load FSA file.")
            return None

        self.detect_peaks()

        # Additional processing (e.g., marker analysis) can be performed here in subclasses.
        results = {
            "fsa_data": self.fsa_data,
            "detected_peaks": self.detected_peaks,
        }
        logger.info("Processing completed successfully.")
        return results