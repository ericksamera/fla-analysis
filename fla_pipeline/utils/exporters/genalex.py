# fla_pipeline/exporters/genalex.py

import io
import csv
import pandas as pd
import ast
from collections import defaultdict

from fla_pipeline.models.genotype import GenotypeResult

from .base import ExporterBase

class GenalexExporter(ExporterBase):
    def __init__(self, metadata=None):
        self.metadata = metadata or {}

    def export(self, df: pd.DataFrame) -> str:


        def extract_alleles(g):
            if isinstance(g, GenotypeResult):
                return g.alleles
            if isinstance(g, dict) and "alleles" in g:
                return g["alleles"]
            if isinstance(g, str):
                if "/" in g:
                    return g.split("/")[:2]
                try:
                    parsed = ast.literal_eval(g)
                    if isinstance(parsed, dict) and "alleles" in parsed:
                        return parsed["alleles"]
                    if isinstance(parsed, (list, tuple)):
                        return list(parsed)
                except Exception:
                    return []
            return []

        # Map population label to each sample
        df["_Pop"] = df["Sample Name"].map(
            lambda s: self.metadata.get(s, {}).get("Population", "Unknown")
        )

        # Sort by population and sample name for block grouping
        df = df.sort_values(by=["_Pop", "Sample Name"]).copy()

        # Identify markers and populations
        markers = sorted(df["Marker"].unique())
        all_samples = df["Sample Name"].unique()
        pop_to_samples = defaultdict(list)
        for sample, pop in zip(df["Sample Name"], df["_Pop"]):
            pop_to_samples[pop].append(sample)

        ordered_pop_labels = sorted(pop_to_samples.keys())
        ordered_samples = [s for pop in ordered_pop_labels for s in pop_to_samples[pop]]
        pop_sizes = [len(pop_to_samples[pop]) for pop in ordered_pop_labels]

        # Row 1: metadata
        row1 = [len(markers), len(ordered_samples), len(ordered_pop_labels)]

        # Row 2: title + population block sizes
        row2 = ["GenAlEx Codominant Export"] + [""] * (2 * len(markers) - 1) + pop_sizes

        # Row 3: headers
        header = ["Sample Name", "Pop"] + [m for marker in markers for m in (marker, "")]

        # Data rows
        data_rows = []
        for sample in ordered_samples:
            sample_df = df[df["Sample Name"] == sample]
            pop_label = sample_df["_Pop"].iloc[0]
            row = [sample, pop_label]

            for marker in markers:
                match = sample_df[sample_df["Marker"] == marker]
                if not match.empty:
                    g = match["Genotype"].iloc[0]
                    alleles = extract_alleles(g)
                    row.extend([str(a) for a in alleles[:2]] + [""] * (2 - len(alleles)))
                else:
                    row.extend(["", ""])
            data_rows.append(row)

        # Write CSV
        output = io.StringIO()
        writer = csv.writer(output)
        writer.writerow(row1)
        writer.writerow(row2)
        writer.writerow(header)
        for row in data_rows:
            writer.writerow([v if v != "" else "" for v in row])
        return output.getvalue()



    def filename(self) -> str:
        return "genalex.csv"

    def _assign_population_codes(self, samples):
        """
        Assigns numeric codes to populations based on metadata["Population"]
        """
        pop_labels = {
            self.metadata.get(s, {}).get("Population", "Unknown")
            for s in samples
        }
        sorted_pops = sorted(pop_labels)
        return {pop: i + 1 for i, pop in enumerate(sorted_pops)}
