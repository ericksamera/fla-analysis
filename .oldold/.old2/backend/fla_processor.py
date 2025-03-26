#!/usr/bin/env python3
__description__ =\
"""
core/fla_processor.py

Diploid FLA Processor that inherits from BaseProcessor (see base_processor.py).
Implements:
  - Repeat-based peak binning
  - Ratio-based 2-allele calls
  - Optional 'snap' to the nearest repeat multiple
  - Extended QC flags
"""
__author__ = "Erick Samera"
__version__ = "1.0.0"
__comments__ = "Works"

import math
import numpy as np
import logging
from dataclasses import dataclass
from typing import List, Tuple, Optional, Dict, Any

# Import from your updated base processor and config
from core.base_processor import BaseProcessor, BaseConfig, Peak

logger = logging.getLogger(__name__)

@dataclass
class FLAConfig(BaseConfig):
    """
    Configuration parameters specifically for diploid FLA processing.
    Inherits from BaseConfig, so it already has attributes like:
      - bin_tolerance
      - min_peak_height
      - ratio_weight
      - etc.
    """
    # You can add diploid-specific parameters here if desired
    pass

@dataclass
class MarkerConfig:
    """
    Marker configuration for FLA processing.
    
    Attributes:
        marker (str): Name of the marker.
        channel (str): Dye channel associated with the marker (e.g. '6-FAM', 'VIC', 'NED').
        repeat_unit (Optional[int]): The expected repeat unit length (e.g. 3 for trinucleotide).
        bins (Tuple[int, int]): The (min, max) range of valid allele sizes for this marker.
        stutter_ratio_threshold (Optional[float]): (Not heavily used here, but you can incorporate if needed.)
        stutter_tolerance (Optional[float]): Tolerance for stutter filtering (if used).
        het_ratio_lower_bound (Optional[float]): Lower bound for ratio to consider 2 peaks heterozygous.
        stutter_direction (str): "left", "right", or "both"; for stutter filtering logic if desired.
    """
    marker: str
    channel: str
    repeat_unit: Optional[int]
    bins: Tuple[int, int]
    stutter_ratio_threshold: Optional[float] = None
    stutter_tolerance: Optional[float] = None
    het_ratio_lower_bound: Optional[float] = None
    stutter_direction: str = "both"

def convert_marker_dicts_to_configs(marker_dict: Dict[str, Dict]) -> List[MarkerConfig]:
    """
    Converts a dictionary of marker definitions into a list of MarkerConfig objects.

    Returns a list of MarkerConfig objects, each built from these entries.
    """
    marker_configs = []
    for marker_name, fields in marker_dict.items():
        config_kwargs = dict(fields)
        config_kwargs["marker"] = marker_name
        marker_config = MarkerConfig(**config_kwargs)
        marker_configs.append(marker_config)
    return marker_configs

class FLAProcessor(BaseProcessor):
    """
    Diploid FLA Processor:

      1) Inherits FSA loading & peak detection from BaseProcessor.
      2) Bins peaks for each marker using repeat-based tolerance if available.
      3) Calls up to 2 alleles (diploid) based on intensity ratio checks.
      4) Optionally aligns final calls to the nearest repeat multiple.
      5) Produces a genotype confidence from a QC penalty score.
    """

    def __init__(self, config: FLAConfig, marker_config: List[dict]):
        """
        :param config: FLAConfig instance with diploid processing parameters.
        :param marker_config: List of MarkerConfig objects describing each marker.
        """
        super().__init__(config)
        self.config: FLAConfig = config
        print(marker_config)
        self.marker_config = convert_marker_dicts_to_configs(marker_config)
        self.marker_results: Dict[str, dict] = {}

    # -------------------- Marker Processing Methods -------------------- #

    def _bin_peaks(
        self,
        raw_peaks: List[Peak],
        b_range: Tuple[int, int],
        repeat_unit: Optional[int]
    ) -> Tuple[List[Peak], List[str]]:
        """
        Filters peaks by a bin range plus a tolerance region based on repeat_unit if given.
        Also applies relative intensity filtering in that bin.

        Returns:
            (filtered_peaks, qc_flags)
        """
        bin_min, bin_max = b_range

        # If repeat_unit known, use repeat_unit*2 as bin expansion; else use config.bin_tolerance
        tol = repeat_unit * self.config.bin_tolerance

        qc_flags = []
        in_bin_peaks: List[Peak] = []

        # First pass: keep any peak that falls within [bin_min - tol, bin_max + tol]
        for pk in raw_peaks:
            if (bin_min - tol) <= pk.position <= (bin_max + tol):
                if not (bin_min <= pk.position <= bin_max):
                    qc_flags.append(
                        f"Peak at {pk.position:.2f} is out-of-bin but within ±{tol} of bin {b_range}."
                    )
                in_bin_peaks.append(pk)

        # Relative intensity filter
        if in_bin_peaks:
            max_intensity = max(p.intensity for p in in_bin_peaks)
            rel_thresh = self.config.relative_peak_threshold
            filtered = [p for p in in_bin_peaks if p.intensity >= rel_thresh * max_intensity]
            if len(filtered) < len(in_bin_peaks):
                qc_flags.append(
                    f"Peaks below {rel_thresh*100:.1f}% of max intensity in bin {b_range} were removed."
                )
            in_bin_peaks = filtered

        return in_bin_peaks, qc_flags

    def _dynamic_ratio_threshold(self, peak_distance: float, repeat_unit: Optional[int]) -> float:
        """
        Decide how high the ratio must be to consider two peaks heterozygous, 
        based on distance in bp vs. the known repeat_unit.
        """
        if not repeat_unit or repeat_unit <= 0:
            return self.config.global_het_ratio_lower_bound

        repeats_away = round(peak_distance / repeat_unit)
        # Example rules:
        if repeats_away <= 1:
            return 0.40  # If they differ by ~1 repeat, we need at least 0.40
        elif repeats_away == 2:
            return 0.30
        else:
            return 0.20

    def _call_alleles_no_combine(
        self,
        b_range: Tuple[int, int],
        current_bin_peaks: List[Peak],
        marker_name: str,
        repeat_unit: Optional[int],
        baseline_ratio: float
    ) -> Tuple[List[Tuple[int, int]], List[dict], List[str], float, Optional[float]]:
        """
        Chooses up to two peaks (diploid). If ratio is too low, we treat it as a single-allele call.

        Returns:
            - allele_calls: list of (allele1, allele2) ints
            - binned_peaks_info: details about the peaks chosen
            - qc_flags
            - qc_penalty
            - heterozygous_fraction (None if homozygous)
        """
        allele_calls: List[Tuple[int, int]] = []
        binned_peaks_info: List[dict] = []
        qc_flags: List[str] = []
        qc_penalties: List[float] = []
        het_fraction: Optional[float] = None

        if not current_bin_peaks:
            qc_flags.append(f"No peaks detected in bin {b_range} for marker {marker_name}.")
            return allele_calls, binned_peaks_info, qc_flags, 0.0, het_fraction

        # Sort peaks by intensity desc
        sorted_peaks = sorted(current_bin_peaks, key=lambda p: p.intensity, reverse=True)

        if len(sorted_peaks) == 1:
            # Single peak => homozygous
            main_peak = sorted_peaks[0]
            pair = (round(main_peak.position), round(main_peak.position))
            allele_calls.append(pair)
            qc_flags.append(f"Only one peak in bin {b_range} for {marker_name}; calling homozygous.")
            qc_penalties.append(self.config.ratio_weight * 0.5)  # a mild penalty
            binned_peaks_info.append({
                "allele_range": b_range,
                "peaks": [vars(main_peak)],
                "note": "Single-peak -> homozygous"
            })
        else:
            # Evaluate top 2
            major = sorted_peaks[0]
            minor = sorted_peaks[1]
            ratio = minor.intensity / major.intensity
            distance = abs(major.position - minor.position)

            # Combine baseline + distance-based threshold
            dynamic_threshold = self._dynamic_ratio_threshold(distance, repeat_unit)
            effective_min_ratio = max(dynamic_threshold, baseline_ratio)

            if effective_min_ratio <= ratio <= 1.0:
                # Heterozygous
                a_call = (round(major.position), round(minor.position))
                allele_calls.append(tuple(sorted(a_call)))
                het_fraction = ratio
                qc_penalties.append(self.config.ratio_weight * abs(ratio - 0.5))
                binned_peaks_info.append({
                    "allele_range": b_range,
                    "peaks": [vars(major), vars(minor)],
                    "ratio": ratio,
                    "distance": distance,
                    "dynamic_threshold": dynamic_threshold,
                })
            else:
                # Too imbalanced => treat as single allele
                qc_flags.append(
                    f"Peaks in bin {b_range} for {marker_name} have ratio {ratio:.2f}; calling homozygous."
                )
                pair = (round(major.position), round(major.position))
                allele_calls.append(pair)
                qc_penalties.append(self.config.ratio_weight * abs(ratio - 0.5))
                binned_peaks_info.append({
                    "allele_range": b_range,
                    "peaks": [vars(major), vars(minor)],
                    "ratio": ratio,
                    "note": "Imbalanced -> homozygous"
                })

        total_penalty = sum(qc_penalties)
        return allele_calls, binned_peaks_info, qc_flags, total_penalty, het_fraction

    def _adjust_alleles(
        self,
        allele_call: Tuple[int, int],
        b_range: Tuple[int, int],
        repeat_unit: int
    ) -> Tuple[Tuple[int, int], List[str], float]:
        """
        Snaps the allele calls to the nearest valid multiple of `repeat_unit` within the bin.

        Returns:
          - new_call (allele1, allele2)
          - qc_flags
          - shift penalty
        """
        bin_min, bin_max = b_range
        qc_flags: List[str] = []
        shift_penalty = 0.0
        snapped: List[int] = []

        for allele in allele_call:
            expected = bin_min + round((allele - bin_min) / repeat_unit) * repeat_unit
            if abs(allele - expected) > self.config.bin_tolerance:
                qc_flags.append(
                    f"Allele {allele} adjusted to {expected} (off by {allele - expected:.1f} bp)."
                )
                shift_penalty += abs(allele - expected) * self.config.intensity_weight
            snapped.append(expected)

        snapped.sort()
        return (snapped[0], snapped[1]), qc_flags, shift_penalty

    def process_markers(self) -> None:
        """
        Loops through each MarkerConfig, bins peaks, calls diploid alleles,
        and stores results in self.marker_results.
        """
        self.marker_results = {}

        # If no data/peaks loaded, do nothing
        if not self.fsa_data or not self.detected_peaks:
            logger.error("No FSA data or peaks to process.")
            return

        for marker in self.marker_config:
            mname = marker.marker
            channel = marker.channel
            b_range = marker.bins
            repeat_unit = marker.repeat_unit

            qc_flags: List[str] = []
            summary_score = 0.0
            heterozygous_fraction = None

            # Retrieve detected peaks for this channel
            if channel not in self.detected_peaks:
                logger.warning(f"Channel {channel} not present for marker {mname}.")
                continue

            raw_peaks = sorted(self.detected_peaks[channel], key=lambda p: p.position)
            # Bin the peaks
            binned_peaks, bin_flags = self._bin_peaks(raw_peaks, b_range, repeat_unit)
            qc_flags.extend(bin_flags)

            # baseline ratio threshold
            base_ratio = marker.het_ratio_lower_bound or self.config.global_het_ratio_lower_bound

            # call alleles
            calls, binned_info, call_flags, penalty, het_frac = self._call_alleles_no_combine(
                b_range, binned_peaks, mname, repeat_unit, base_ratio
            )
            qc_flags.extend(call_flags)
            summary_score += penalty
            if het_frac is not None:
                heterozygous_fraction = het_frac

            # adjust calls to nearest repeat boundary if desired
            parsed_calls = []
            if calls and repeat_unit and repeat_unit > 0:
                # Typically only 1 call in 'calls', so we adjust the last one
                final_call = calls[-1]
                adjusted_call, adj_flags, adj_penalty = self._adjust_alleles(
                    final_call, b_range, repeat_unit
                )
                qc_flags.extend(adj_flags)
                summary_score += adj_penalty
                parsed_calls.append(adjusted_call)
            elif calls:
                # If no repeat_unit, just store them as-is
                parsed_calls.append(calls[-1])

            # genotype confidence
            genotype_conf = float(np.exp(-self.config.qc_confidence_scale * summary_score))

            self.marker_results[mname] = {
                "channel": channel,
                "binned_peaks": binned_info,
                "allele_calls": calls,
                "parsed_allele_calls": parsed_calls,
                "heterozygous_fraction": heterozygous_fraction,
                "QC_flags": qc_flags,
                "QC_score": summary_score,
                "genotype_confidence": genotype_conf,
            }
            logger.info(f"Marker {mname} => {self.marker_results[mname]}")

    # -------------------- Main Pipeline -------------------- #

    def process_all(self, file_path: str) -> Optional[Dict[str, Any]]:
        """
        High-level method that:
          1) Loads the FSA data
          2) Detects peaks (via BaseProcessor)
          3) Processes markers (diploid logic)

        Returns a dictionary with the results (fsa_data, detected_peaks, marker_results)
        or None on failure.
        """
        # Step 1: load
        if not super().load_fsa(file_path):
            logger.error("Failed to load FSA file.")
            return None

        # Step 2: detect
        super().detect_peaks()

        # Step 3: per-marker analysis
        self.process_markers()

        # Prepare final output
        results = {
            "fsa_data": self.fsa_data,
            "detected_peaks": self.detected_peaks,
            "marker_results": self.marker_results,
        }
        logger.info("Diploid FLA processing complete.")
        return results


# -------------------- Command-Line Interface for Testing -------------------- #
if __name__ == "__main__":
    import argparse
    import json
    import os

    parser = argparse.ArgumentParser(description="Diploid FLA Processor (No peak combining, ratio-based calls).")
    parser.add_argument("file", type=str, help="Path to the .fsa file")
    args = parser.parse_args()

    if not os.path.exists(args.file):
        print(f"Error: File {args.file} does not exist.")
        exit(1)

    # Example marker configurations
    example_markers = [
        MarkerConfig(marker="Marker1", channel="6-FAM", repeat_unit=3, bins=(180, 240)),
        MarkerConfig(marker="Marker2", channel="VIC",   repeat_unit=4, bins=(200, 260)),
        MarkerConfig(marker="Marker3", channel="NED",   repeat_unit=3, bins=(150, 210)),
    ]

    # Create an FLAConfig and pass it + the marker config
    config = FLAConfig()
    processor = FLAProcessor(config, example_markers)

    result = processor.process_all(args.file)
    if result:
        # Show JSON output
        print(json.dumps(result, indent=4))
    else:
        print("Processing failed.")
