import pandas as pd
from collections import Counter
from scipy.stats import chisquare

def compute_hwe_stats(genotype_df: pd.DataFrame, marker_configs: dict) -> list[dict]:
    """
    Computes per-marker allele frequencies and HWE test results.
    Returns a list of dictionaries with marker, frequencies, and p-value (if available).
    """
    results = []

    for marker in sorted(marker_configs.keys()):
        sub_df = genotype_df[genotype_df["Marker"] == marker]
        alleles = []

        for g in sub_df["Genotype"]:
            parts = str(g).split("/")
            if len(parts) == 2:
                alleles.extend(parts)

        if not alleles:
            continue

        allele_counts = Counter(alleles)
        total_alleles = sum(allele_counts.values())
        freqs = {allele: count / total_alleles for allele, count in allele_counts.items()}
        sorted_alleles = sorted(freqs.items(), key=lambda x: -x[1])

        entry = {
            "marker": marker,
            "allele_frequencies": {a: round(f, 4) for a, f in sorted_alleles},
            "p_value": None,
            "hwe_tested": False,
            "hwe_message": ""
        }

        if len(allele_counts) == 2:
            alleles_sorted = sorted(allele_counts)
            p = freqs[alleles_sorted[0]]
            q = freqs[alleles_sorted[1]]
            expected = {
                f"{alleles_sorted[0]}/{alleles_sorted[0]}": p**2 * len(sub_df),
                f"{alleles_sorted[1]}/{alleles_sorted[1]}": q**2 * len(sub_df),
                f"{alleles_sorted[0]}/{alleles_sorted[1]}": 2*p*q * len(sub_df),
            }

            observed = Counter(sub_df["Genotype"])
            observed_counts = [observed.get(g, 0) for g in expected]
            expected_counts = list(expected.values())

            if all(e >= 1 for e in expected_counts):
                chi2, pval = chisquare(f_obs=observed_counts, f_exp=expected_counts)
                entry["p_value"] = round(pval, 4)
                entry["hwe_tested"] = True
            else:
                entry["hwe_message"] = "HWE test skipped (expected counts too low)."
        else:
            entry["hwe_message"] = "HWE test only available for biallelic markers."

        results.append(entry)

    return results


def display_hwe_results(hwe_results: list[dict]):
    import streamlit as st
    st.subheader("Hardy-Weinberg Equilibrium & Allele Summary")
    with st.expander("Per-marker allele frequencies and HWE test", expanded=False):
        for result in hwe_results:
            st.markdown(f"**{result['marker']}**")
            st.write("Allele frequencies:", {a: f"{f:.2f}" for a, f in result["allele_frequencies"].items()})

            if result["hwe_tested"]:
                st.write(f"Hardy-Weinberg p-value: `{result['p_value']}`")
            elif result["hwe_message"]:
                st.write(result["hwe_message"])