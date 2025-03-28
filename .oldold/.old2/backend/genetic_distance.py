#!/usr/bin/env python3
__description__ =\
"""
core/genetic_distance.py

This module defines a modular framework for computing genetic distance matrices from genotype data.
It includes an abstract base class for distance metrics and two implementations:
    - BruvoDistanceMetric (for Bruvo’s genetic distance)
    - NeiDistanceMetric (a simple Nei-like distance based on shared alleles)
The GeneticDistanceCalculator class accepts genotype data (as a pandas DataFrame) along with a list
of MarkerConfig objects and a chosen distance metric. It computes a distance matrix between samples,
averaging over markers with data in common.
"""
__author__ = "Erick Samera"
__version__ = "1.2.1"
__comments__ = ""

import pandas as pd
import numpy as np
import itertools
from collections import Counter
from typing import List, Dict, Tuple, Any, Optional
from abc import ABC, abstractmethod

from core.config_models import MarkerConfig

# --------------------- Utility Functions ---------------------

def perform_pcoa(distance_matrix: pd.DataFrame, n_components: int = 2) -> Tuple[np.ndarray, np.ndarray]:
    """
    Perform Principal Coordinates Analysis (PCoA) on the distance matrix and
    return the coordinates and explained variance ratios.
    """
    D = distance_matrix.values.astype(float)
    n = D.shape[0]
    if n < n_components:
        raise ValueError(f"Not enough samples for {n_components} components. Only {n} available.")
    D2 = D ** 2
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J.dot(D2).dot(J)
    eigvals, eigvecs = np.linalg.eigh(B)
    idx = np.argsort(eigvals)[::-1]
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]
    pos_idx = eigvals > 0
    coords = eigvecs[:, pos_idx] * np.sqrt(eigvals[pos_idx])
    if coords.shape[1] < n_components:
        coords = np.hstack([coords, np.zeros((n, n_components - coords.shape[1]))])
    explained_variance_ratio = eigvals / eigvals.sum()
    return coords[:, :n_components], explained_variance_ratio

def flatten_if_nested(thing: Any) -> List[float]:
    """
    Recursively flattens nested lists/tuples into a one-dimensional list of floats.
    """
    if isinstance(thing, (float, int)):
        return [float(thing)]
    if isinstance(thing, (list, tuple)):
        flattened = []
        for elem in thing:
            flattened.extend(flatten_if_nested(elem))
        return flattened
    try:
        return [float(thing)]
    except Exception:
        return []

def to_repeat_count(allele_size: float, repeat_unit: float) -> int:
    """
    Convert an allele size to its approximate repeat count using the given repeat unit.
    """
    if repeat_unit <= 0:
        repeat_unit = 1e-5
    return int(round(allele_size / repeat_unit))

def ensure_unique_sample_names(df: pd.DataFrame) -> pd.DataFrame:
    """
    Ensure that the "Sample Name" column in the DataFrame is unique.
    
    If a "Filename" column is present:
      1. For each row, extract the base sample name as the first token of "Sample Name" (split by "_").
      2. Create a unique key by combining the base sample name with the filename.
         This key ensures that rows from the same file are grouped together.
      3. Then, for all unique keys, if multiple keys share the same base name 
         (coming from different filenames), append a counter (e.g., base, base_2, ...).
         
    If no "Filename" column exists, the original duplicate-handling behavior is used.
    """
    if "Filename" in df.columns:
        # 1. Extract base sample names from the "Sample Name" column.
        df["Base Sample"] = df["Sample Name"].astype(str).str.strip().apply(lambda x: x.split("_")[0])
        # 2. Create a unique key that combines the base sample and the filename.
        df["Unique Key"] = df.apply(lambda row: f"{row['Base Sample']}|{row['Filename'].strip()}", axis=1)
        
        # 3. Build a mapping from each base sample to the list of unique keys that have that base.
        unique_keys = df["Unique Key"].unique()
        mapping = {}
        for key in unique_keys:
            base = key.split("|")[0]
            mapping.setdefault(base, []).append(key)
        
        # 4. Assign final sample names: if only one unique key exists for a base name, use the base;
        #    if multiple exist (from different filenames), append a counter.
        final_names = {}
        for base, keys in mapping.items():
            if len(keys) == 1:
                final_names[keys[0]] = base
            else:
                for i, key in enumerate(keys, start=1):
                    final_names[key] = base if i == 1 else f"{base}_{i}"
        
        # 5. Replace "Sample Name" with the final sample names.
        df["Sample Name"] = df["Unique Key"].map(final_names)
        
        # Optionally, drop the temporary columns.
        df.drop(columns=["Base Sample", "Unique Key"], inplace=True)
        return df

    # Fallback: if no "Filename" column exists, use the previous duplicate-handling method.
    if "Sample Name" in df.columns:
        df["Sample Name"] = df["Sample Name"].astype(str).str.strip()
        counts = df["Sample Name"].value_counts()
        duplicates = counts[counts > 1].index.tolist()
        if not duplicates:
            return df
        sample_name_map = {}
        for sample in duplicates:
            count = 1
            for idx in df[df["Sample Name"] == sample].index:
                new_name = sample if count == 1 else f"{sample}_{count}"
                sample_name_map[idx] = new_name
                count += 1
        for idx, new_name in sample_name_map.items():
            df.at[idx, "Sample Name"] = new_name
        return df

    return df

# --------------------- Distance Metric Base Class ---------------------

class GeneticDistanceMetric(ABC):
    """
    Abstract base class for genetic distance metrics.
    """
    @abstractmethod
    def compute_distance(self, geno1: Any, geno2: Any, marker: str, marker_repeat: int) -> float:
        """
        Computes the genetic distance between two genotype calls.
        
        Parameters:
            geno1, geno2: The genotype data (typically lists or nested lists of allele sizes).
            marker: The marker name (if needed for context).
            marker_repeat: The repeat unit for the marker.
            
        Returns:
            A float between 0 (identical) and 1 (maximally different).
        """
        pass

# --------------------- Bruvo's Distance Implementation ---------------------

class BruvoDistanceMetric(GeneticDistanceMetric):
    """
    Computes Bruvo's genetic distance.
    """
    def compute_distance(self, geno1: Any, geno2: Any, marker: str, marker_repeat: int) -> float:
        if not geno1 or not geno2:
            return 1.0
        # Convert allele sizes to repeat counts after flattening the input.
        geno1_counts = [to_repeat_count(a, marker_repeat) for a in flatten_if_nested(geno1)]
        geno2_counts = [to_repeat_count(a, marker_repeat) for a in flatten_if_nested(geno2)]
        n1, n2 = len(geno1_counts), len(geno2_counts)
        if n1 == 0 or n2 == 0:
            return 1.0

        def _min_perm_distance(gcounts1: List[int], gcounts2: List[int]) -> float:
            best = np.inf
            for perm in itertools.permutations(gcounts2):
                distances = []
                for a_val, b_val in zip(gcounts1, perm):
                    x = abs(a_val - b_val)
                    d = 1 - 2 ** (-x)
                    distances.append(d)
                avg = sum(distances) / len(distances) if distances else 1.0
                best = min(best, avg)
            return best if best != np.inf else 1.0

        if n1 != n2:
            # Allow genome addition/loss: expand the smaller list
            if n1 < n2:
                best_score = np.inf
                for exp in itertools.product(geno1_counts, repeat=(n2 - n1)):
                    full_geno1 = list(geno1_counts) + list(exp)
                    score = _min_perm_distance(full_geno1, geno2_counts)
                    best_score = min(best_score, score)
                return best_score
            else:
                best_score = np.inf
                for exp in itertools.product(geno2_counts, repeat=(n1 - n2)):
                    full_geno2 = list(geno2_counts) + list(exp)
                    score = _min_perm_distance(geno1_counts, full_geno2)
                    best_score = min(best_score, score)
                return best_score
        else:
            return _min_perm_distance(geno1_counts, geno2_counts)

# --------------------- Nei's Distance Implementation ---------------------

class NeiDistanceMetric(GeneticDistanceMetric):
    """
    Computes a simple Nei-like genetic distance based on the proportion of shared alleles.
    """
    def compute_distance(self, geno1: Any, geno2: Any, marker: str, marker_repeat: int = 1) -> float:
        alleles1 = flatten_if_nested(geno1)
        alleles2 = flatten_if_nested(geno2)
        if not alleles1 or not alleles2:
            return 1.0
        count1 = Counter(alleles1)
        count2 = Counter(alleles2)
        shared = 0
        for allele, cnt in count1.items():
            shared += min(cnt, count2.get(allele, 0))
        total = max(len(alleles1), len(alleles2))
        similarity = shared / total
        if similarity <= 0:
            return 1.0
        return -np.log(similarity)

# --------------------- Genetic Distance Calculator ---------------------

class GeneticDistanceCalculator:
    """
    Computes a genetic distance matrix from genotype data.
    
    Parameters:
        df (pd.DataFrame): Genotype data with at least the following columns:
            - "Sample Name": unique identifier for each sample.
            - "Marker": marker name.
            - "Genotype": genotype call (e.g., a list of allele sizes).
        marker_configs (List[MarkerConfig]): List of marker configurations.
        metric: Either a string ("bruvo" or "nei") or an instance of GeneticDistanceMetric.
    """
    def __init__(self, df: pd.DataFrame, marker_configs: List[MarkerConfig],
                 metric: Optional[Any] = "bruvo"):
        self.df = ensure_unique_sample_names(df.copy())
        self.marker_configs = marker_configs
        # Create a mapping for marker repeat units.
        self.marker_repeat = {mc.marker: (mc.repeat_unit if mc.repeat_unit is not None else 1)
                              for mc in marker_configs}
        # Select the distance metric.
        if isinstance(metric, str):
            metric = metric.lower()
            if metric == "bruvo":
                self.metric = BruvoDistanceMetric()
            elif metric == "nei":
                self.metric = NeiDistanceMetric()
            else:
                raise ValueError(f"Unsupported metric: {metric}")
        elif isinstance(metric, GeneticDistanceMetric):
            self.metric = metric
        else:
            raise ValueError("Metric must be either a string or an instance of GeneticDistanceMetric.")
        # Parse genotype data.
        self.genotypes, self.marker_modes = self._parse_and_impute_genotypes(self.df)
        self.samples = sorted(self.genotypes.keys())

    def _parse_and_impute_genotypes(self, df: pd.DataFrame) -> Tuple[Dict[str, Dict[str, Any]], Dict[str, List[float]]]:
        """
        Parse genotype data into a dictionary mapping sample to marker calls.
        Impute missing calls using the most common (mode) call for each marker.
        
        Returns:
            genotypes: {sample: {marker: genotype}}
            marker_modes: {marker: mode genotype}
        """
        marker_modes = {}
        genotypes = {}
        if df.empty:
            raise ValueError("Input DataFrame is empty.")
        for marker, group in df.groupby("Marker"):
            genotype_list = []
            for _, row in group.iterrows():
                parsed = row["Genotype"]
                if parsed:
                    genotype_list.append(parsed)
            if genotype_list:
                mode_tuple = Counter(tuple(g) for g in genotype_list if g).most_common(1)[0][0]
                marker_modes[marker] = list(mode_tuple)
        for _, row in df.iterrows():
            sample = row["Sample Name"]
            marker = row["Marker"]
            parsed = row["Genotype"]
            if not parsed:
                parsed = marker_modes.get(marker, [])
            genotypes.setdefault(sample, {})[marker] = parsed
        if not genotypes:
            raise ValueError("No valid genotype data found.")
        return genotypes, marker_modes

    def compute_distance_matrix(self) -> pd.DataFrame:
        """
        Compute the genetic distance matrix among all samples.
        
        For each pair of samples, the distance is computed as the average over markers that
        have genotype data for both samples.
        
        Returns:
            A symmetric pandas DataFrame where rows and columns correspond to sample names.
        """
        dist_matrix = pd.DataFrame(np.nan, index=self.samples, columns=self.samples)
        for i, s1 in enumerate(self.samples):
            for j, s2 in enumerate(self.samples):
                if j < i:
                    continue
                common_markers = set(self.genotypes[s1].keys()) & set(self.genotypes[s2].keys())
                if not common_markers:
                    d = 1.0
                else:
                    d_list = []
                    for marker in common_markers:
                        geno1 = self.genotypes[s1][marker]
                        geno2 = self.genotypes[s2][marker]
                        r_unit = self.marker_repeat.get(marker, 1)
                        d_marker = self.metric.compute_distance(geno1, geno2, marker, r_unit)
                        d_list.append(d_marker)
                    d = sum(d_list) / len(d_list) if d_list else 1.0
                dist_matrix.at[s1, s2] = d
                dist_matrix.at[s2, s1] = d
        np.fill_diagonal(dist_matrix.values, 0)
        return dist_matrix

# --------------------- Command-Line Interface for Testing ---------------------

if __name__ == "__main__":
    import argparse
    import sys
    import ast

    parser = argparse.ArgumentParser(description="Compute a genetic distance matrix from genotype data.")
    parser.add_argument("datafile", type=str, help="Path to a CSV file containing genotype data.")
    parser.add_argument("--metric", type=str, default="bruvo", choices=["bruvo", "nei"],
                        help="The genetic distance metric to use (default: bruvo).")
    args = parser.parse_args()

    try:
        df = pd.read_csv(args.datafile)
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        sys.exit(1)

    # The CSV file should have columns: "Sample Name", "Marker", and "Genotype".
    # The "Genotype" column is expected to be a string representation of a list (e.g., "[100, 102]").
    df["Genotype"] = df["Genotype"].apply(lambda x: ast.literal_eval(x) if pd.notnull(x) else [])

    # For demonstration, create default MarkerConfig objects for each unique marker.
    unique_markers = df["Marker"].unique()
    marker_configs = []
    for marker in unique_markers:
        # Assign a default repeat_unit of 1 and a dummy bin range.
        marker_configs.append(MarkerConfig(marker=marker, channel="", repeat_unit=1, bins=(0, 100)))

    try:
        calculator = GeneticDistanceCalculator(df, marker_configs, metric=args.metric)
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)

    matrix = calculator.compute_distance_matrix()
    print("Genetic Distance Matrix:")
    print(matrix)
