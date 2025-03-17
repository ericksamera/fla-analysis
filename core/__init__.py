#!/usr/bin/env python3
"""
core/__init__.py

This package contains the core business logic for FLA processing.
It includes:
  - Base processing functionality (BaseProcessor)
  - Diploid and polyploid processing modules
  - Configuration models and utility functions

All heavy computation and file processing happens here, separated from the UI layer.
"""

from .base_processor import BaseProcessor
from .config_models import BaseConfig, FLAConfig, MarkerConfig
from .data_utils import gaussian, fit_gaussian_to_xy_data, aggregate_qc_penalties
from .fla_processor import FLAProcessor
from .fla_processor_poly import FLAPolyProcessor

__all__ = [
    "BaseProcessor",
    "BaseConfig",
    "FLAConfig",
    "MarkerConfig",
    "gaussian",
    "fit_gaussian_to_xy_data",
    "aggregate_qc_penalties",
    "FLAProcessor",
    "FLAPolyProcessor"
]