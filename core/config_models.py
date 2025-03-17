#!/usr/bin/env python3
__description__ =\
"""
core/config_models.py

This module defines configuration data models for FLA processing.
It includes:
  - BaseConfig: Common parameters shared across processing modes.
  - FLAConfig: Configuration for diploid FLA processing (inherits from BaseConfig).
  - FLAConfigPoly: Configuration for polyploid FLA processing with a different default ploidy.
  - MarkerConfig: Model for marker settings (channel, expected allele range, repeat unit, etc.)
"""
__author__ = "Erick Samera"
__version__ = "1.0.0"
__comments__ = "Works"

from dataclasses import dataclass
from typing import Tuple, Optional

@dataclass
class BaseConfig:
    """
    Base configuration parameters for FLA processing.
    These parameters are common to both diploid and polyploid processing.
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
    global_het_ratio_lower_bound: float = 0.3
    gaussian_fit_window: int = 5
    ploidy: int = 2

@dataclass
class FLAConfig(BaseConfig):
    """
    Configuration parameters for diploid FLA processing.
    Inherits all parameters from BaseConfig.
    """
    pass


@dataclass
class MarkerConfig:
    """
    Marker configuration for FLA processing.
    
    Attributes:
        marker (str): Name of the marker.
        channel (str): Dye channel associated with the marker.
        repeat_unit (Optional[int]): Expected repeat unit for the marker.
        bins (Tuple[int, int]): Tuple (min, max) representing the expected allele size range.
        stutter_ratio_threshold (Optional[float]): Threshold for stutter filtering.
        stutter_tolerance (Optional[float]): Tolerance for stutter filtering.
        het_ratio_lower_bound (Optional[float]): Lower bound for heterozygosity ratio.
        stutter_direction (str): Expected stutter direction; options: "left", "right", or "both".
    """
    marker: str
    channel: str
    repeat_unit: Optional[int]
    bins: Tuple[int, int]
    stutter_ratio_threshold: Optional[float] = None
    stutter_tolerance: Optional[float] = None
    het_ratio_lower_bound: Optional[float] = None
    stutter_direction: str = "both"
