import pandas as pd
from collections import defaultdict

def calculate_fst_per_marker(genotype_df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculates per-marker FST using the classic Nei 1973 definition:
    FST = (HT - HS) / HT
    where:
      - HT is total expected heterozygosity (global allele frequencies)
      - HS is mean expected heterozygosity within subpopulations
    """
    results = []
    grouped = genotype_df.groupby("Marker")

    for marker, marker_df in grouped:
        pop_groups = marker_df.groupby("_Pop")
        allele_freqs_by_pop = {}
        hs_list = []

        for pop, group in pop_groups:
            alleles = []
            for g in group["Genotype"]:
                parts = str(g).split("/")
                if len(parts) == 2:
                    alleles.extend(parts)

            if not alleles:
                continue

            total = len(alleles)
            counts = defaultdict(int)
            for a in alleles:
                counts[a] += 1
            freqs = {a: c / total for a, c in counts.items()}
            allele_freqs_by_pop[pop] = freqs

            hs = 1 - sum(f ** 2 for f in freqs.values())
            hs_list.append(hs)

        if not allele_freqs_by_pop or not hs_list:
            continue

        # Calculate total (global) frequencies
        global_counts = defaultdict(int)
        total_alleles = 0
        for freqs in allele_freqs_by_pop.values():
            for a, f in freqs.items():
                count = f * 2  # Each individual has two alleles
                global_counts[a] += count
                total_alleles += count

        if total_alleles == 0:
            continue

        global_freqs = {a: c / total_alleles for a, c in global_counts.items()}
        ht = 1 - sum(f ** 2 for f in global_freqs.values())
        hs = sum(hs_list) / len(hs_list)

        if ht > 0:
            fst = (ht - hs) / ht
        else:
            fst = 0.0

        results.append({"Marker": marker, "FST": fst, "HT": ht, "HS": hs})

    return pd.DataFrame(results).sort_values("FST", ascending=False)

