#!/usr/bin/env python3
"""
core/fla_processor_poly.py

This module implements the polyploid FLA processor by extending the BaseProcessor.
It focuses on advanced filtering and allele assignment for polyploid analysis using a
fraction-based approach, bin filtering, stutter filtering, and optional allele adjustment.
"""
__author__ = "Erick Samera"
__version__ = "1.0.0"
__comments__ = "Works"

import logging
import math
import numpy as np
from dataclasses import dataclass
from typing import List, Tuple, Optional, Dict

from core.base_processor import BaseProcessor, BaseConfig, Peak

logger = logging.getLogger(__name__)

# ------------------------------------------------------------------
# Data Classes for Polyploid Configuration and Marker Definitions
# ------------------------------------------------------------------

@dataclass
class FLAConfigPoly(BaseConfig):
    """
    Configuration parameters for polyploid FLA processing.
    Inherits common parameters from BaseConfig.
    """
    ploidy: int = 4  # Default polyploid level (e.g., tetraploid)

@dataclass
class MarkerConfig:
    """
    Marker configuration for processing.
    Each marker defines a channel, expected bin range, repeat unit, and optional stutter parameters.
    """
    marker: str
    channel: str
    repeat_unit: Optional[int]
    bins: Tuple[int, int]
    stutter_ratio_threshold: Optional[float] = None
    stutter_tolerance: Optional[float] = None
    het_ratio_lower_bound: Optional[float] = None
    stutter_direction: str = "both"  # Options: "left", "right", "both"

# ------------------------------------------------------------------
# Polyploid FLA Processor Class
# ------------------------------------------------------------------

class FLAPolyProcessor(BaseProcessor):
    """
    Polyploid FLA Processor that extends BaseProcessor.
    
    This class inherits FSA loading and peak detection from BaseProcessor.
    It then applies marker-specific processing to:
      - Bin peaks based on marker configurations (using a tolerance that may depend on repeat_unit)
      - Optionally filter out stutter peaks
      - Perform fraction-based allele calling (distributing copies according to intensity)
      - Adjust allele calls to the nearest repeat boundary
      - Compute QC flags and a genotype confidence score
    """
    
    def __init__(self, config: FLAConfigPoly, marker_config: Optional[List[MarkerConfig]] = None):
        super().__init__(config)
        self.config: FLAConfigPoly = config
        self.marker_config: List[MarkerConfig] = marker_config if marker_config is not None else self.default_marker_config()
        self.marker_results: Dict[str, dict] = {}

    def default_marker_config(self) -> List[MarkerConfig]:
        """
        Provides default marker configurations.
        Modify these defaults as needed.
        """
        return [
            MarkerConfig(marker="Dox09", channel="6-FAM", repeat_unit=3, bins=(181, 233), stutter_direction="left"),
            MarkerConfig(marker="Dox10", channel="NED",   repeat_unit=3, bins=(211, 234), stutter_direction="right"),
            MarkerConfig(marker="Dox11", channel="VIC",   repeat_unit=3, bins=(215, 300), stutter_direction="right"),
        ]
    
    def _bin_peaks(self, raw_peaks: List[Peak], b_range: Tuple[int, int], repeat_unit: Optional[int]) -> Tuple[List[Peak], List[str]]:
        """
        Filters peaks based on a bin range plus a tolerance.
        Tolerance is determined as repeat_unit.
        Also applies relative intensity filtering.
        
        Returns a tuple (filtered peaks, QC flags).
        """
        bin_min, bin_max = b_range
        tol = repeat_unit * self.config.bin_tolerance

        current_bin_peaks: List[Peak] = []
        qc_flags: List[str] = []

        for peak in raw_peaks:
            if (bin_min - tol) <= peak.position <= (bin_max + tol):
                if not (bin_min <= peak.position <= bin_max):
                    qc_flags.append(f"Peak at {peak.position:.2f} is out-of-bin but within ±{tol} bp of bin {b_range}.")
                current_bin_peaks.append(peak)

        if current_bin_peaks:
            max_intensity = max(p.intensity for p in current_bin_peaks)
            rel_thresh = self.config.relative_peak_threshold
            filtered_peaks = [p for p in current_bin_peaks if p.intensity >= rel_thresh * max_intensity]
            if len(filtered_peaks) < len(current_bin_peaks):
                qc_flags.append(f"Peaks below {rel_thresh*100:.0f}% of max intensity in bin {b_range} were filtered out.")
            current_bin_peaks = filtered_peaks

        return current_bin_peaks, qc_flags

    @staticmethod
    def _filter_stutter_peaks(peaks: List[Peak],
                              repeat_unit: float, 
                              stutter_ratio_threshold: float,
                              tolerance: float, 
                              stutter_direction: str) -> List[Peak]:
        """
        Removes likely stutter peaks based on the expected repeat_unit, a ratio threshold,
        and a tolerance value. The stutter_direction parameter determines which side(s) are considered.
        """
        if not peaks:
            return peaks

        sorted_peaks = sorted(peaks, key=lambda p: p.intensity, reverse=True)
        accepted_peaks: List[Peak] = []
        removed_indices = set()

        for i, main_peak in enumerate(sorted_peaks):
            if i in removed_indices:
                continue
            accepted_peaks.append(main_peak)
            for j in range(i+1, len(sorted_peaks)):
                if j in removed_indices:
                    continue
                candidate = sorted_peaks[j]
                stutter_offset = main_peak.position - candidate.position
                if abs(stutter_offset - repeat_unit) <= tolerance:
                    if (stutter_direction in ["left", "both"] and stutter_offset > 0) or \
                       (stutter_direction in ["right", "both"] and stutter_offset < 0):
                        if candidate.intensity <= stutter_ratio_threshold * main_peak.intensity:
                            removed_indices.add(j)
        final_peaks = [p for i, p in enumerate(sorted_peaks) if i not in removed_indices]
        return sorted(final_peaks, key=lambda p: p.intensity, reverse=True)

    def _call_alleles_polyploid(self, current_bin_peaks: List[Peak]) -> List[int]:
        """
        Uses a fraction-based method to assign allele copies.
        Each peak's intensity contributes fractionally to a copy count scaled to config.ploidy.
        Returns a sorted list of allele positions.
        """
        if not current_bin_peaks:
            return []

        total_intensity = sum(p.intensity for p in current_bin_peaks)
        if total_intensity <= 0:
            return []

        float_counts = [ (p.intensity / total_intensity) * self.config.ploidy for p in current_bin_peaks ]
        floor_counts = [int(math.floor(fc)) for fc in float_counts]
        sum_floor = sum(floor_counts)
        leftover = self.config.ploidy - sum_floor

        remainders = [fc - int(math.floor(fc)) for fc in float_counts]
        idx_sorted_by_rema = sorted(range(len(current_bin_peaks)), key=lambda i: remainders[i], reverse=True)
        final_counts = floor_counts[:]
        if leftover > 0:
            for idx in idx_sorted_by_rema[:leftover]:
                final_counts[idx] += 1
        elif leftover < 0:
            leftover_to_remove = abs(leftover)
            idx_sorted_asc = sorted(range(len(current_bin_peaks)), key=lambda i: remainders[i])
            for idx in idx_sorted_asc:
                if leftover_to_remove <= 0:
                    break
                if final_counts[idx] > 0:
                    final_counts[idx] -= 1
                    leftover_to_remove -= 1

        allele_positions = []
        for count, peak in zip(final_counts, current_bin_peaks):
            allele_positions.extend([round(peak.position)] * count)
        allele_positions.sort()
        return allele_positions

    def _adjust_alleles(self,
                        allele_positions: List[int],
                        b_range: Tuple[int, int],
                        repeat_unit: int) -> Tuple[List[int], List[str], float]:
        """
        Adjusts each allele call to the nearest repeat boundary.
        Returns the adjusted allele calls, any QC flags, and a penalty score based on deviation.
        """
        bin_min, bin_max = b_range
        adjusted_calls: List[int] = []
        qc_flags: List[str] = []
        score = 0.0

        for allele in allele_positions:
            expected_val = bin_min + round((allele - bin_min) / repeat_unit) * repeat_unit
            if abs(allele - expected_val) > self.config.bin_tolerance:
                qc_flags.append(f"Allele call {allele} in bin {b_range} adjusted to {expected_val} (excess dev).")
                score += abs(allele - expected_val) * self.config.intensity_weight
            adjusted_calls.append(expected_val)
        adjusted_calls.sort()
        return adjusted_calls, qc_flags, score

    def process_markers(self) -> None:
        """
        Processes markers by:
          1) Sorting detected peaks and binning them based on each marker's expected range.
          2) Optionally applying stutter filtering.
          3) Assigning allele copies fractionally.
          4) Adjusting allele calls to repeat boundaries.
          5) Computing QC flags and genotype confidence.
          6) Storing the results in self.marker_results.
        """
        self.marker_results = {}

        for marker in self.marker_config:
            marker_name = marker.marker
            channel = marker.channel
            b_range = marker.bins
            repeat_unit = marker.repeat_unit

            qc_flags: List[str] = []
            summary_qc_score = 0.0

            # Sort peaks by position in the given channel.
            raw_peaks = sorted(self.detected_peaks.get(channel, []), key=lambda p: p.position)
            
            # Bin the peaks according to the marker's range.
            current_bin_peaks, qc_flags_bin = self._bin_peaks(raw_peaks, b_range, repeat_unit)
            qc_flags.extend(qc_flags_bin)

            # Optionally filter stutter peaks if parameters are configured.
            if (repeat_unit is not None and repeat_unit > 0 and
                marker.stutter_ratio_threshold is not None and
                marker.stutter_tolerance is not None):
                current_bin_peaks = self._filter_stutter_peaks(
                    current_bin_peaks,
                    repeat_unit,
                    marker.stutter_ratio_threshold,
                    marker.stutter_tolerance,
                    marker.stutter_direction
                )

            # Perform fraction-based allele assignment.
            allele_positions = self._call_alleles_polyploid(current_bin_peaks)

            # Adjust allele calls to repeat boundaries if applicable.
            parsed_alleles = []
            if allele_positions and repeat_unit is not None and repeat_unit > 0:
                adj_alleles, qc_flags_adj, penalty_adj = self._adjust_alleles(allele_positions, b_range, repeat_unit)
                qc_flags.extend(qc_flags_adj)
                summary_qc_score += penalty_adj
                parsed_alleles = adj_alleles
            else:
                parsed_alleles = allele_positions

            # Compute a quick genotype confidence based on the QC penalty.
            genotype_confidence = float(np.exp(-self.config.qc_confidence_scale * summary_qc_score))

            binned_peaks_info = [vars(p) for p in current_bin_peaks]
            self.marker_results[marker_name] = {
                "channel": channel,
                "binned_peaks": binned_peaks_info,
                "allele_calls": allele_positions,
                "parsed_allele_calls": parsed_alleles,
                "QC_flags": qc_flags,
                "QC_score": summary_qc_score,
                "genotype_confidence": genotype_confidence
            }
            logger.info(f"Marker {marker_name} processed: alleles assigned {parsed_alleles}.")

    def process_all(self, file_path: str) -> bool:
        """
        High-level entry point:
          1) Loads the FSA file.
          2) Detects peaks (using BaseProcessor functionality).
          3) Processes markers for polyploid allele calling.
        
        Returns True if processing succeeds; otherwise, False.
        """
        if self.load_fsa(file_path):
            self.detect_peaks()
            self.process_markers()
            return True
        return False

# ------------------------------------------------------------------
# Testing / Command-Line Interface
# ------------------------------------------------------------------

if __name__ == "__main__":
    import argparse
    import json
    import os

    parser = argparse.ArgumentParser(description="Process an .fsa file for polyploid FLA with advanced filtering.")
    parser.add_argument("file", type=str, help="Path to the .fsa file")
    args = parser.parse_args()

    if not os.path.exists(args.file):
        print(f"Error: File {args.file} does not exist.")
        exit(1)

    # Create an instance of the polyploid processor using default configuration and marker settings.
    config = FLAConfigPoly()  # Use default polyploid configuration parameters.
    processor = FLAPolyProcessor(config)
    
    if processor.process_all(args.file):
        output = {
            "file": processor.fsa_data["name"],
            "fsa_data": {
                "smap": processor.fsa_data["smap"],
                "channels": { ch: list(map(float, arr)) for ch, arr in processor.fsa_data["channels"].items() }
            },
            "detected_peaks": {
                channel: [ { "position": peak.position, "intensity": peak.intensity, "saturated": peak.saturated, "note": peak.note } for peak in peaks ]
                for channel, peaks in processor.detected_peaks.items()
            },
            "marker_results": processor.marker_results
        }
        print(json.dumps(output, indent=4))
    else:
        print("Processing failed.")
