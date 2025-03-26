#!/usr/bin/env python3
"""
FLA Processor Module (Polyploid + Advanced Filtering Support)
------------------------------------------------------------------
This module implements a fragment length analysis (FLA) pipeline that:
 1) Loads .fsa files and detects peaks (including a saturation check
    and local Gaussian fitting for high-intensity peaks).
 2) Performs bin filtering to keep only peaks in the marker’s expected range.
 3) Performs stutter filtering based on stutter_ratio_threshold and tolerance.
 4) Uses a fraction-based approach to assign up to `ploidy` allele copies
    in each marker range. (E.g., in a tetraploid, if one peak has ~50% of
    the total intensity, it’s assigned two copies, etc.)
 5) Optionally adjusts final allele calls to repeat boundaries if needed.
 6) Collects QC flags, genotype confidence, and logs all results.

Author: Erick Samera (Refactored, with polyploid fraction-based allele assignment)
Version: 1.3.0-polyploid
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

# ------------------------------------------------------------------
# Data Classes for Configuration and Peak Data
# ------------------------------------------------------------------

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
    ploidy: int = 2  # For polyploid calling, can be 2 (diploid), 3 (triploid), 4 (tetraploid), etc.

@dataclass
class MarkerConfig:
    marker: str
    channel: str
    repeat_unit: Optional[int]
    bins: Tuple[int, int]
    stutter_ratio_threshold: Optional[float] = None
    stutter_tolerance: Optional[float] = None
    het_ratio_lower_bound: Optional[float] = None
    stutter_direction: str = "both"  # "left", "right", or "both"


# ------------------------------------------------------------------
# Helper Functions for Gaussian Fitting and QC
# ------------------------------------------------------------------

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

# ------------------------------------------------------------------
# Main FLAProcessor Class (Polyploid + Filter Logic)
# ------------------------------------------------------------------

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

    # ------------------------------------------------------------------
    # Saturation & Peak Detection
    # ------------------------------------------------------------------

    def _estimate_saturated_peak_intensity(self, data_np: np.ndarray, smap: np.ndarray, peak_idx: int) -> float:
        """
        For a peak suspected of saturation, fit a Gaussian in a local window around the peak
        to estimate the true amplitude. Return the fitted amplitude if it's higher.
        """
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
                # Check threshold
                if raw_intensity >= self.config.min_peak_height:
                    pos = smap[i] if i < len(smap) else None
                    if pos is not None and pos >= self.config.min_peak_position:
                        saturated = (raw_intensity >= self.config.instrument_saturation_value)
                        if saturated:
                            corrected_intensity = self._estimate_saturated_peak_intensity(data_np, smap, i)
                            peak = Peak(position=float(pos),
                                        intensity=float(corrected_intensity),
                                        saturated=True,
                                        note="Saturated peak; corrected via Gaussian.")
                        else:
                            peak = Peak(position=float(pos), intensity=float(raw_intensity))
                        peaks.append(peak)

            # Sort peaks by intensity descending
            peaks.sort(key=lambda p: p.intensity, reverse=True)
            self.detected_peaks[channel] = peaks

    # ------------------------------------------------------------------
    # Binning / Stutter Filtering
    # ------------------------------------------------------------------

    def _bin_peaks(self,
                   raw_peaks: List[Peak],
                   b_range: Tuple[int, int],
                   repeat_unit: Optional[int]) -> Tuple[List[Peak], List[str]]:
        """
        Filters peaks by a bin range plus a tolerance region based on repeat_unit,
        and applies a relative intensity threshold. 
        """
        bin_min, bin_max = b_range

        # Use repeat_unit to set a bigger tolerance if we want (e.g. ±2 repeats)
        if repeat_unit is not None and repeat_unit > 0:
            tol = repeat_unit * 2
        else:
            tol = self.config.bin_tolerance

        current_bin_peaks: List[Peak] = []
        qc_flags: List[str] = []

        # Keep all peaks within bin ± tol
        for peak in raw_peaks:
            if (bin_min - tol) <= peak.position <= (bin_max + tol):
                # If strictly out-of-bin but within tol => flag it
                if not (bin_min <= peak.position <= bin_max):
                    qc_flags.append(
                        f"Peak at {peak.position:.2f} is out-of-bin but within ±{tol} bp of bin {b_range}."
                    )
                current_bin_peaks.append(peak)

        # Apply relative-intensity filtering
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

                # "Left" stutter => candidate is smaller (left on bp axis)
                # "Right" stutter => candidate is bigger (right on bp axis)
                if abs(stutter_offset - repeat_unit) <= tolerance:
                    # stutter_offset>0 => candidate is to the left
                    if (stutter_direction in ["left", "both"] and stutter_offset > 0):
                        if candidate.intensity <= stutter_ratio_threshold * main_peak.intensity:
                            removed_indices.add(j)
                    # stutter_offset<0 => candidate is to the right
                    elif (stutter_direction in ["right", "both"] and stutter_offset < 0):
                        if candidate.intensity <= stutter_ratio_threshold * main_peak.intensity:
                            removed_indices.add(j)

        # Reconstruct the final list from accepted_peaks
        final_peaks = [p for i, p in enumerate(sorted_peaks) if i not in removed_indices]
        return sorted(final_peaks, key=lambda p: p.intensity, reverse=True)

    # ------------------------------------------------------------------
    # Polyploid Allele Calling (Fraction-Based)
    # ------------------------------------------------------------------

    def _call_alleles_polyploid(self, current_bin_peaks: List[Peak]) -> List[int]:
        """
        Fraction-based approach to distribute up to `config.ploidy` copies across the set of peaks.
        E.g. if ploidy=4 and a peak is ~50% of total intensity => 2 copies, etc.
        """
        if not current_bin_peaks:
            return []

        total_intensity = sum(p.intensity for p in current_bin_peaks)
        if total_intensity <= 0:
            return []

        # fraction for each peak
        float_counts = []
        for p in current_bin_peaks:
            frac = p.intensity / total_intensity
            # Multiply fraction by ploidy => "ideal" copy count
            float_counts.append(frac * self.config.ploidy)

        # Floor them, then distribute leftover
        floor_counts = [int(math.floor(fc)) for fc in float_counts]
        sum_floor = sum(floor_counts)
        leftover = self.config.ploidy - sum_floor

        # Sort peaks by fractional remainder descending
        remainders = [(fc - math.floor(fc)) for fc in float_counts]
        idx_sorted_by_rema = sorted(range(len(current_bin_peaks)), key=lambda i: remainders[i], reverse=True)

        final_counts = floor_counts[:]
        if leftover > 0:
            # Give leftover copies to largest remainders
            for idx in idx_sorted_by_rema[:leftover]:
                final_counts[idx] += 1
        elif leftover < 0:
            # We assigned too many copies => remove from smallest remainders first
            leftover_to_remove = abs(leftover)
            idx_sorted_asc = sorted(range(len(current_bin_peaks)), key=lambda i: remainders[i], reverse=False)
            for idx in idx_sorted_asc:
                if leftover_to_remove <= 0:
                    break
                if final_counts[idx] > 0:
                    final_counts[idx] -= 1
                    leftover_to_remove -= 1

        # Construct final allele set => each peak's position repeated final_counts[i] times
        allele_positions = []
        for count, peak in zip(final_counts, current_bin_peaks):
            allele_positions.extend([round(peak.position)] * count)
        # Sort ascending by position
        allele_positions.sort()

        return allele_positions

    def _adjust_alleles(
        self,
        allele_positions: List[int],
        b_range: Tuple[int, int],
        repeat_unit: int
    ) -> Tuple[List[int], List[str], float]:
        """
        Aligns each called allele to the nearest repeat boundary, if it deviates significantly.
        Returns the adjusted alleles, plus any QC flags + penalty score.
        """
        bin_min, bin_max = b_range
        adjusted_calls: List[int] = []
        qc_flags: List[str] = []
        score = 0.0

        for allele in allele_positions:
            # Round to nearest repeat boundary if we know repeat_unit
            # e.g. nearest multiple of repeat_unit
            #   offset = (allele - bin_min)
            #   # approximate # repeats
            #   n_repeats = round(offset / repeat_unit)
            #   new_value = bin_min + n_repeats * repeat_unit
            expected_val = bin_min + round((allele - bin_min) / repeat_unit) * repeat_unit
            # If difference is > bin_tolerance => penalize
            if abs(allele - expected_val) > self.config.bin_tolerance:
                qc_flags.append(
                    f"Allele call {allele} in bin {b_range} adjusted to {expected_val} (excess dev)."
                )
                score += abs(allele - expected_val) * self.config.intensity_weight

            adjusted_calls.append(expected_val)

        # Sort ascending
        adjusted_calls.sort()
        return adjusted_calls, qc_flags, score

    # ------------------------------------------------------------------
    # Putting it All Together in process_markers()
    # ------------------------------------------------------------------

    def process_markers(self) -> None:
        """
        1) Sort peaks by position.
        2) Bin them according to each MarkerConfig (± tolerance).
        3) Optionally do stutter filtering.
        4) Do fraction-based copy assignment up to ploidy.
        5) Adjust calls to nearest repeat boundary if requested.
        6) Compute QC + genotype confidence.
        7) Store results in self.marker_results[marker_name].
        """
        self.marker_results = {}

        for marker in self.marker_config:
            marker_name = marker.marker
            channel = marker.channel
            b_range = marker.bins
            repeat_unit = marker.repeat_unit

            qc_flags: List[str] = []
            summary_qc_score = 0.0

            # 1) Sort raw peaks by position
            raw_peaks = sorted(self.detected_peaks.get(channel, []), key=lambda p: p.position)
            
            # 2) Bin peaks
            current_bin_peaks, qc_flags_bin = self._bin_peaks(raw_peaks, b_range, repeat_unit)
            qc_flags.extend(qc_flags_bin)

            # 3) Stutter filtering if user configured it
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

            # 4) Fraction-based allele assignment
            allele_positions = self._call_alleles_polyploid(current_bin_peaks)

            # 5) If there's a repeat_unit, optionally align to boundary
            parsed_alleles = []
            if allele_positions and repeat_unit is not None and repeat_unit > 0:
                adj_alleles, qc_flags_adj, penalty_adj = self._adjust_alleles(
                    allele_positions, b_range, repeat_unit
                )
                qc_flags.extend(qc_flags_adj)
                summary_qc_score += penalty_adj
                parsed_alleles = adj_alleles
            else:
                parsed_alleles = allele_positions

            # 6) Quick QC confidence
            # You can add your own ratio or stutter-based penalties if you want
            genotype_confidence = float(np.exp(-self.config.qc_confidence_scale * summary_qc_score))

            # 7) Store results
            # For debugging: store the final binned peaks as well
            binned_peaks_info = [vars(p) for p in current_bin_peaks]
            self.marker_results[marker_name] = {
                "channel": channel,
                "binned_peaks": binned_peaks_info,
                "allele_calls": allele_positions,      # Raw fraction-based calls (unadjusted)
                "parsed_allele_calls": parsed_alleles, # Possibly adjusted
                "QC_flags": qc_flags,
                "QC_score": summary_qc_score,
                "genotype_confidence": genotype_confidence
            }

    def process_all(self) -> bool:
        """
        High-level entry point: load data, detect peaks, process markers (including bin, stutter, polyploid calls).
        """
        if self.load_fsa():
            self.detect_peaks()
            self.process_markers()
            return True
        return False


# -------------------- Testing Section --------------------
if __name__ == "__main__":
    import argparse
    import json
    import os

    parser = argparse.ArgumentParser(description="Process an .fsa file for polyploid FLA, with advanced filtering.")
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
