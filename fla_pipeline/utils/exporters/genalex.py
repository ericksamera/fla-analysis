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

        if "sample_uid" not in df.columns:
            raise ValueError("Genotype results must include sample_uid column.")

        def get_population(uid):
            meta = self.metadata.get(uid)
            if meta and isinstance(meta, dict):
                return meta.get("Population", "Unknown")
            print(f"⚠️ No metadata match for UID: {uid}")
            return "Unknown"

        df["_Pop"] = df["sample_uid"].map(get_population)

        # Sort and prep
        df = df.sort_values(by=["_Pop", "Sample Name"]).copy()
        markers = sorted(df["Marker"].unique())
        unique_uids = df["sample_uid"].unique()

        pop_to_uids = defaultdict(list)
        for uid in unique_uids:
            pop = get_population(uid)
            pop_to_uids[pop].append(uid)

        ordered_pops = sorted(pop_to_uids.keys())
        ordered_uids = [uid for pop in ordered_pops for uid in pop_to_uids[pop]]
        pop_sizes = [len(pop_to_uids[pop]) for pop in ordered_pops]

        # Headers
        row1 = [len(markers), len(ordered_uids), len(ordered_pops)] + pop_sizes
        row2 = ["GenAlEx Codominant Export"] + ["", ""] + [pop for pop in ordered_pops] 
        header = ["Sample Name", "Pop"] + [m for marker in markers for m in (marker, "")]

        # Data rows
        data_rows = []
        for uid in ordered_uids:
            sample_df = df[df["sample_uid"] == uid]
            if sample_df.empty:
                continue

            sample_name = sample_df["Sample Name"].iloc[0]
            pop = sample_df["_Pop"].iloc[0]

            row = [sample_name, pop]
            for marker in markers:
                match = sample_df[sample_df["Marker"] == marker]
                if not match.empty:
                    g = match["Genotype"].iloc[0]
                    alleles = extract_alleles(g)
                    if not alleles:
                        row.extend(["0", "0"])
                    else:
                        row.extend([str(a) for a in alleles[:2]] + ["0"] * (2 - len(alleles)))
                else:
                    row.extend(["0", "0"])
            data_rows.append(row)

        # Write CSV
        output = io.StringIO()
        writer = csv.writer(output)
        writer.writerow(row1)
        writer.writerow(row2)
        writer.writerow(header)
        writer.writerows(data_rows)

        return output.getvalue()


    def filename(self) -> str:
        return "genalex.csv"
