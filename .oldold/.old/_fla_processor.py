#!/usr/bin/env python3
"""
FLA Processor Module (Modified - No Peak Combining, Extended Raw Output, Saturation Handling)
--------------------------------------------------------------------------------------------
This module implements the FLAProcessor class to process fragment length analysis (FLA)
data from .fsa files. In this version:

1. Peaks are not combined via clustering.
2. We incorporate a local Gaussian fit for peaks that may be saturated at the instrument’s detection limit,
   potentially extending beyond the recorded max.
3. Distance-based ratio thresholds remain for determining heterozygous calls.

Author: Erick Samera (Refactored, with new distance-based ratio thresholds and saturation handling)
Version: 1.2.0-modified
"""

import numpy as np
import math
from dataclasses import dataclass
from typing import List, Tuple, Optional, Dict
from Bio import SeqIO
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
import logging

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# --- Data Classes for Configuration and Peak Data ---

@dataclass
class Peak:
    position: float
    intensity: float
    saturated: bool = False  # Flag if peak was instrument-saturated
    note: str = ""           # Additional annotations if needed

@dataclass
class FLAConfig:
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
class MarkerConfig:
    marker: str
    channel: str
    repeat_unit: Optional[int]
    bins: Tuple[int, int]
    stutter_ratio_threshold: Optional[float] = None
    stutter_tolerance: Optional[float] = None
    het_ratio_lower_bound: Optional[float] = None
    stutter_direction: str = "both"  # "left" for Dox09, "right" for others


# --- Helper Functions for Gaussian Fitting and QC ---

def gaussian(x: np.ndarray, A: float, mu: float, sigma: float) -> np.ndarray:
    """A standard Gaussian function."""
    return A * np.exp(-((x - mu) ** 2) / (2 * sigma ** 2))

def fit_gaussian_to_xy_data(x_data: np.ndarray, y_data: np.ndarray) -> Optional[Tuple[float, float, float]]:
    """
    Fits a Gaussian to (x_data, y_data) and returns (amplitude, mu, sigma).
    Returns None if the fit fails or data is insufficient.
    """
    if len(x_data) < 3:
        return None

    A0 = np.max(y_data)
    mu0 = x_data[np.argmax(y_data)]
    sigma0 = np.std(x_data) if np.std(x_data) > 0 else 1.0

    try:
        popt, _ = curve_fit(gaussian, x_data, y_data, p0=[A0, mu0, sigma0])
        # popt is [A, mu, sigma]
        return popt[0], popt[1], popt[2]
    except Exception as e:
        logger.debug(f"Gaussian fit failed: {e}")
        return None

def aggregate_qc_penalties(penalties: List[float]) -> float:
    """Sums up penalty values to produce a final QC penalty."""
    return sum(penalties)

# --- Main FLAProcessor Class (No Peak Combining, Distanced-Based Ratio Logic) ---

class FLAProcessor:
    def __init__(self, file_path: str, 
                 config: Optional[FLAConfig] = None, 
                 marker_config: Optional[List[MarkerConfig]] = None):
        self.file_path = file_path
        self.config = config if config is not None else FLAConfig()
        self.marker_config = marker_config if marker_config is not None else self.default_marker_config()
        self.fsa_data = None
        self.detected_peaks: Dict[str, List[Peak]] = {}
        self.marker_results: Dict[str, dict] = {}

    def default_marker_config(self) -> List[MarkerConfig]:
        return [
            MarkerConfig(marker="Dox09", channel="6-FAM", repeat_unit=3, bins=(181, 233), stutter_direction="left"),
            MarkerConfig(marker="Dox10", channel="NED",   repeat_unit=3, bins=(211, 234), stutter_direction="right"),
            MarkerConfig(marker="Dox11", channel="VIC",   repeat_unit=3, bins=(215, 300), stutter_direction="right"),
        ]

    def load_fsa(self) -> bool:
        """Loads the .fsa file into memory, extracting SMAP and channel data arrays."""
        try:
            fsa_object = SeqIO.read(self.file_path, 'abi')
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

    def _estimate_saturated_peak_intensity(self, data_np: np.ndarray, smap: np.ndarray, peak_idx: int) -> float:
        """
        For a peak suspected of saturation, fit a Gaussian in a local window around the peak
        to estimate the true amplitude. Return the fitted amplitude if it's higher.
        """
        # Prepare local window around the peak
        window = self.config.gaussian_fit_window
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

    def detect_peaks(self) -> None:
        """Detects peaks in each channel using SciPy's find_peaks, with saturation handling."""
        if self.fsa_data is None:
            logger.error("FSA data not loaded.")
            return

        smap = np.array(self.fsa_data["smap"], dtype=float)
        self.detected_peaks = {}

        for channel, data in self.fsa_data["channels"].items():
            data_np = np.array(data, dtype=float)
            # Basic peak detection
            indices = find_peaks(data_np)[0]

            peaks: List[Peak] = []
            for i in indices:
                raw_intensity = data_np[i]
                # Check if this peak is above min threshold and within position range
                if raw_intensity >= self.config.min_peak_height:
                    pos = smap[i] if i < len(smap) else None
                    if pos is not None and pos >= self.config.min_peak_position:
                        # Check if at or near saturation
                        saturated = False
                        if raw_intensity >= self.config.instrument_saturation_value:
                            saturated = True
                            # Attempt local Gaussian fit
                            corrected_intensity = self._estimate_saturated_peak_intensity(data_np, smap, i)
                            peak = Peak(position=float(pos),
                                        intensity=float(corrected_intensity),
                                        saturated=True,
                                        note="Saturated peak; corrected via Gaussian.")
                        else:
                            peak = Peak(position=float(pos),
                                        intensity=float(raw_intensity),
                                        saturated=False)
                        peaks.append(peak)

            # Sort peaks by intensity in descending order
            peaks.sort(key=lambda p: p.intensity, reverse=True)
            self.detected_peaks[channel] = peaks

    def _bin_peaks(self,
                   raw_peaks: List[Peak],
                   b_range: Tuple[int, int],
                   repeat_unit: Optional[int]) -> Tuple[List[Peak], List[str]]:
        """
        Filters peaks by a bin range plus a tolerance region based on repeat_unit,
        and applies a relative intensity threshold. Peaks that fall outside bin_min..bin_max
        but within this tolerance are flagged yet still included.
        """
        bin_min, bin_max = b_range

        # Use repeat_unit to set tolerance if available, else fall back to config.bin_tolerance
        if repeat_unit is not None:
            tol = repeat_unit * 2
        else:
            tol = self.config.bin_tolerance

        current_bin_peaks: List[Peak] = []
        qc_flags: List[str] = []

        # Collect all peaks that lie within (bin_min - tol, bin_max + tol)
        for peak in raw_peaks:
            if (bin_min - tol) <= peak.position <= (bin_max + tol):
                # If it's outside the strict bin, but within tol, just flag it
                if not (bin_min <= peak.position <= bin_max):
                    qc_flags.append(
                        f"Peak at {peak.position:.2f} is out-of-bin but within ±{tol} bp of bin {b_range}."
                    )
                current_bin_peaks.append(peak)

        # Apply relative-intensity filtering if we have any peaks
        if current_bin_peaks:
            max_intensity = max(p.intensity for p in current_bin_peaks)
            rel_thresh = self.config.relative_peak_threshold
            filtered_peaks = [p for p in current_bin_peaks if p.intensity >= rel_thresh * max_intensity]
            if len(filtered_peaks) < len(current_bin_peaks):
                qc_flags.append(
                    f"Peaks below {rel_thresh*100:.0f}% of max intensity in bin {b_range} were filtered out."
                )
            current_bin_peaks = filtered_peaks

        return current_bin_peaks, qc_flags

    @staticmethod
    def _filter_stutter_peaks(peaks: List[Peak],
                              repeat_unit: float, 
                              stutter_ratio_threshold: float,
                              tolerance: float, 
                              stutter_direction: str) -> List[Peak]:
        """
        Removes likely stutter peaks, factoring in the stutter direction.
        - "left": stutter is smaller than the main peak
        - "right": stutter is larger
        - "both": remove stutter on both sides
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

                # "Left" stutter => candidate is to the left of main peak
                if abs(stutter_offset - repeat_unit) <= tolerance:
                    if (stutter_direction in ["left", "both"] and stutter_offset > 0):
                        if candidate.intensity <= stutter_ratio_threshold * main_peak.intensity:
                            removed_indices.add(j)
                    elif (stutter_direction in ["right", "both"] and stutter_offset < 0):
                        if candidate.intensity <= stutter_ratio_threshold * main_peak.intensity:
                            removed_indices.add(j)

        return accepted_peaks

    def _compute_weighted_average(self, peaks: List[Peak]) -> float:
        """
        Computes an intensity-weighted average of peak positions.
        """
        total_intensity = sum(p.intensity for p in peaks)
        if total_intensity <= 0:
            return peaks[0].position
        return sum(p.position * p.intensity for p in peaks) / total_intensity

    def _dynamic_ratio_threshold(self, peak_distance: float, repeat_unit: Optional[int]) -> float:
        """
        Decide how high the ratio must be for calling heterozygosity, based on how many repeat
        units away the second peak is from the first. If only 1 repeat away, require a near 1:1 ratio.
        Larger distances allow a lower ratio.
        """
        if repeat_unit is None:
            # If no repeat unit, just use a baseline
            return self.config.global_het_ratio_lower_bound

        repeats_away = round(peak_distance / repeat_unit)

        # Example logic
        if repeats_away <= 1:
            return 0.4
        elif repeats_away == 2:
            return 0.3
        else:
            return 0.2

    def _call_alleles_no_combine(
        self,
        b_range: Tuple[int, int],
        current_bin_peaks: List[Peak],
        marker_name: str,
        repeat_unit: Optional[int],
        baseline_ratio: float
    ) -> Tuple[List[Tuple[int, int]], List[dict], List[str], float, Optional[float]]:
        """
        Processes peaks without combining them into clusters.
        Selects top two peaks by intensity for allele calling,
        with distance-based thresholds to decide heterozygosity.
        """
        allele_calls: List[Tuple[int, int]] = []
        binned_peaks_info: List[dict] = []
        qc_flags: List[str] = []
        qc_penalties: List[float] = []
        heterozygous_fraction: Optional[float] = None

        if not current_bin_peaks:
            qc_flags.append(f"No peak detected in bin {b_range} for marker {marker_name}")
            return allele_calls, binned_peaks_info, qc_flags, 0.0, heterozygous_fraction

        # Sort by intensity descending
        sorted_peaks = sorted(current_bin_peaks, key=lambda p: p.intensity, reverse=True)

        if len(sorted_peaks) == 1:
            # Only one peak => homozygous
            main = sorted_peaks[0]
            allele_call = (round(main.position), round(main.position))
            allele_calls.append(allele_call)
            qc_flags.append(f"Only one peak detected in bin {b_range} for marker {marker_name}.")
            qc_penalties.append(self.config.ratio_weight * 0.5)
            binned_peaks_info.append({
                "allele_range": b_range,
                "peaks": [vars(main)],
                "note": "Single peak"
            })
        else:
            # Take top two
            major = sorted_peaks[0]
            minor = sorted_peaks[1]
            ratio = minor.intensity / major.intensity

            # Distance-based threshold
            peak_distance = abs(major.position - minor.position)
            dynamic_threshold = self._dynamic_ratio_threshold(peak_distance, repeat_unit)
            effective_het_lower_bound = dynamic_threshold

            if effective_het_lower_bound <= ratio <= 1.0:
                # Heterozygous
                allele_call = (round(major.position), round(minor.position))
                allele_calls.append(tuple(sorted(allele_call)))
                heterozygous_fraction = ratio
                qc_penalties.append(self.config.ratio_weight * abs(ratio - 0.5))
                binned_peaks_info.append({
                    "allele_range": b_range,
                    "peaks": [vars(major), vars(minor)],
                    "ratio": ratio,
                    "repeat_unit": repeat_unit,
                    "peak_distance": peak_distance,
                    "distance_threshold": dynamic_threshold,
                })
            else:
                # Imbalanced ratio => treat as single allele call
                qc_flags.append(
                    f"Peaks in bin {b_range} for {marker_name} have imbalanced ratio "
                    f"({ratio:.2f}) with {peak_distance:.1f}bp separation."
                )
                allele_call = (round(major.position), round(major.position))
                allele_calls.append(allele_call)
                qc_penalties.append(self.config.ratio_weight * abs(ratio - 0.5))
                binned_peaks_info.append({
                    "allele_range": b_range,
                    "peaks": [vars(major), vars(minor)],
                    "ratio": ratio,
                    "note": "Ambiguous ratio",
                    "repeat_unit": repeat_unit,
                    "peak_distance": peak_distance,
                    "distance_threshold": dynamic_threshold
                })

        overall_penalty = aggregate_qc_penalties(qc_penalties)
        return allele_calls, binned_peaks_info, qc_flags, overall_penalty, heterozygous_fraction

    def _adjust_alleles(self,
                        allele_call: Tuple[int, int],
                        b_range: Tuple[int, int],
                        repeat_unit: int) -> Tuple[Tuple[int, int], List[str], float]:
        """
        Adjusts the allele calls so they align with expected repeat units.
        """
        bin_min, _ = b_range
        adjusted_call = []
        qc_flags = []
        score = 0.0

        for allele in allele_call:
            expected_val = bin_min + round((allele - bin_min) / repeat_unit) * repeat_unit
            if abs(allele - expected_val) > self.config.bin_tolerance:
                qc_flags.append(
                    f"Allele call {allele} in bin {b_range} deviates from expected structure; "
                    f"adjusted to {expected_val}."
                )
                score += abs(allele - expected_val) * self.config.intensity_weight
            adjusted_call.append(expected_val)

        return tuple(sorted(adjusted_call)), qc_flags, score

    def process_markers(self) -> None:
        """
        Main method to bin peaks by marker, filter stutters, call alleles,
        adjust allele calls to the nearest repeat unit, and compute QC.
        """
        self.marker_results = {}

        for marker in self.marker_config:
            marker_name = marker.marker
            channel = marker.channel
            b_range = marker.bins
            repeat_unit = marker.repeat_unit

            qc_flags: List[str] = []
            summary_qc_score = 0.0
            heterozygous_fraction: Optional[float] = None

            # Sort raw peaks by position
            raw_peaks = sorted(self.detected_peaks.get(channel, []), key=lambda p: p.position)
            current_bin_peaks, qc_flags_bin = self._bin_peaks(raw_peaks, b_range, repeat_unit)
            qc_flags.extend(qc_flags_bin)

            # Marker-specific stutter filtering if parameters are available
            if repeat_unit is not None and \
               marker.stutter_ratio_threshold is not None and \
               marker.stutter_tolerance is not None:
                current_bin_peaks = self._filter_stutter_peaks(
                    current_bin_peaks,
                    repeat_unit,
                    marker.stutter_ratio_threshold,
                    marker.stutter_tolerance,
                    marker.stutter_direction
                )

            # Baseline ratio threshold for heterozygosity
            baseline_ratio = (marker.het_ratio_lower_bound
                              if marker.het_ratio_lower_bound is not None
                              else self.config.global_het_ratio_lower_bound)

            # Call alleles
            allele_calls, binned_peaks, qc_flags_call, score_call, heterozygous_fraction = \
                self._call_alleles_no_combine(
                    b_range, current_bin_peaks, marker_name, repeat_unit, baseline_ratio
                )
            qc_flags.extend(qc_flags_call)
            summary_qc_score += score_call

            # Adjust allele calls to nearest repeat boundary if needed
            parsed_allele_calls = []
            if repeat_unit is not None and allele_calls:
                raw_call = allele_calls[-1]  # last call made
                adjusted_call, qc_flags_adj, score_adj = self._adjust_alleles(raw_call, b_range, repeat_unit)
                qc_flags.extend(qc_flags_adj)
                summary_qc_score += score_adj
                parsed_allele_calls.append(adjusted_call)
            elif allele_calls:
                parsed_allele_calls.append(allele_calls[-1])

            # Compute genotype confidence
            genotype_confidence = float(np.exp(-self.config.qc_confidence_scale * summary_qc_score))

            self.marker_results[marker_name] = {
                "channel": channel,
                "binned_peaks": binned_peaks,
                "allele_calls": allele_calls,
                "parsed_allele_calls": parsed_allele_calls,
                "heterozygous_fraction": heterozygous_fraction,
                "QC_flags": qc_flags,
                "QC_score": summary_qc_score,
                "genotype_confidence": genotype_confidence
            }

    def process_all(self) -> bool:
        """
        High-level entry point to load data, detect peaks (with saturation handling),
        then process markers.
        """
        if self.load_fsa():
            self.detect_peaks()
            self.process_markers()
            return True
        else:
            return False


# -------------------- Testing Section --------------------
if __name__ == "__main__":
    import argparse
    import json
    import os

    parser = argparse.ArgumentParser(description="Process an .fsa file for fragment length analysis with saturation handling.")
    parser.add_argument("file", type=str, help="Path to the .fsa file")
    args = parser.parse_args()

    if not os.path.exists(args.file):
        print(f"Error: File {args.file} does not exist.")
        exit(1)

    processor = FLAProcessor(args.file)
    success = processor.process_all()
    if success:
        output = {
            "file": processor.fsa_data["name"],
            "fsa_data": {
                "smap": processor.fsa_data["smap"],
                "channels": {
                    ch: list(map(float, arr)) for ch, arr in processor.fsa_data["channels"].items()
                }
            },
            "detected_peaks": {
                channel: [
                    {
                        "position": peak.position,
                        "intensity": peak.intensity,
                        "saturated": peak.saturated,
                        "note": peak.note
                    }
                    for peak in peaks
                ]
                for channel, peaks in processor.detected_peaks.items()
            },
            "marker_results": processor.marker_results
        }
        print(json.dumps(output, indent=4))
    else:
        print("Processing failed.")
