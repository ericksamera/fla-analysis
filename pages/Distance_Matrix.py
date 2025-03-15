import streamlit as st
import pandas as pd
import numpy as np
import itertools
import ast
import plotly.express as px
import plotly.figure_factory as ff
import plotly.graph_objects as go
from collections import Counter

st.title("Hierarchical Clustering & PCoA (Plotly)")

# --- UI for selecting clustering method ---
method = st.selectbox(
    "Select linkage method for hierarchical clustering",
    ["single", "complete", "average", "weighted", "centroid", "median", "ward"]
)

###############################################################################
# 1) Helper to parse genotype strings robustly
###############################################################################
def parse_genotype_string(genotype_str):
    """
    Attempt multiple ways to parse a genotype string like:
      "228, 231" or "228 231" or even "[228,231]" or "(228, 231)"
    Always returns a simple list of floats, e.g. [228.0, 231.0].
    """
    # If it's already a list/tuple in Python, just return it directly
    if not isinstance(genotype_str, str):
        # e.g. if it's already a list [228,231], return as float
        return flatten_if_nested(genotype_str)

    # 1) Try literal_eval (handles "[228,231]" or "(228,231)" or similar).
    try:
        obj = ast.literal_eval(genotype_str)
        return flatten_if_nested(obj)
    except:
        pass

    # 2) Try splitting by comma
    try:
        return [float(x.strip()) for x in genotype_str.split(",") if x.strip()]
    except:
        pass

    # 3) Try splitting by whitespace
    try:
        return [float(x.strip()) for x in genotype_str.split() if x.strip()]
    except:
        pass

    # If all fail, return an empty list
    return []

def flatten_if_nested(thing):
    """
    If 'thing' is a nested list or tuple, flatten it into a 1D list of floats.
    E.g. ( (228, 231), ) -> [228, 231]
         [(204.0, 204.0)] -> [204.0, 204.0]
    """
    if isinstance(thing, (float, int)):
        return [float(thing)]

    if isinstance(thing, (list, tuple)):
        flattened = []
        for elem in thing:
            if isinstance(elem, (list, tuple)):
                flattened.extend(flatten_if_nested(elem))
            else:
                flattened.append(float(elem))
        return flattened

    # If we still get something weird, wrap it as single float if possible
    try:
        return [float(thing)]
    except:
        return []

###############################################################################
# 2) Convert raw allele sizes to repeat counts (rounding), with tiny offset if needed
###############################################################################
def to_repeat_count(allele_size, repeat_unit):
    if repeat_unit <= 0:
        repeat_unit = 1e-5
    return int(round(allele_size / repeat_unit))

###############################################################################
# 3) Bruvo’s distance implementation that expands smaller genotype if requested
###############################################################################
def bruvo_distance_between_genotypes(geno1_raw, geno2_raw, repeat_unit,
                                     genome_add=True, genome_loss=True):
    """
    Convert each allele to 'repeat counts' then compute Bruvo’s distance:
        d = 1 - 2^(-|difference in repeat counts|)
    If genome_add or genome_loss are True, we attempt expansions for partial ploidy.
    """

    if not geno1_raw or not geno2_raw:
        return 1.0  # If either is empty, default distance is 1

    # Convert to repeat counts
    geno1_counts = [to_repeat_count(a, repeat_unit) for a in geno1_raw]
    geno2_counts = [to_repeat_count(a, repeat_unit) for a in geno2_raw]

    n1, n2 = len(geno1_counts), len(geno2_counts)
    if n1 == 0 or n2 == 0:
        return 1.0

    if genome_add or genome_loss:
        # Expand the smaller set in a simplistic "genome addition" approach
        if n1 < n2:
            expansions = itertools.product(geno1_counts, repeat=(n2 - n1))
            best_score = np.inf
            for exp in expansions:
                g1_full = list(geno1_counts) + list(exp)
                score = _min_perm_distance(g1_full, geno2_counts)
                if score < best_score:
                    best_score = score
            return best_score

        elif n2 < n1:
            expansions = itertools.product(geno2_counts, repeat=(n1 - n2))
            best_score = np.inf
            for exp in expansions:
                g2_full = list(geno2_counts) + list(exp)
                score = _min_perm_distance(geno1_counts, g2_full)
                if score < best_score:
                    best_score = score
            return best_score

        else:
            # equal length
            return _min_perm_distance(geno1_counts, geno2_counts)
    else:
        # Just do the old approach: multiply the smaller set to match size
        max_len = max(n1, n2)
        if n1 < max_len:
            geno1_counts *= max_len
        if n2 < max_len:
            geno2_counts *= max_len
        return _min_perm_distance(geno1_counts, geno2_counts)

def _min_perm_distance(gcounts1, gcounts2):
    """Compute the minimal average Bruvo distance across all permutations of gcounts2."""
    best = np.inf
    for perm in itertools.permutations(gcounts2):
        distances = []
        for a_val, b_val in zip(gcounts1, perm):
            x = abs(a_val - b_val)
            # Bruvo's formula:  d = 1 - 2^(-x)
            d = 1 - 2 ** (-x)
            distances.append(d)
        avg = sum(distances)/len(distances) if distances else 1.0
        if avg < best:
            best = avg
    return best if best != np.inf else 1.0

###############################################################################
# 4) Main script with improved parsing
###############################################################################
if ("genotype_results_df" in st.session_state
    and st.session_state.genotype_results_df is not None
    and not st.session_state.genotype_results_df.empty):

    df = st.session_state.genotype_results_df.copy()
    marker_counts = df.groupby("Marker").size()
    valid_markers = marker_counts[marker_counts >= 2].index.tolist()

    if valid_markers:
        df_valid = df[df["Marker"].isin(valid_markers)]
        st.write("Markers with genotype data for at least two individuals:")
        st.dataframe(df_valid)
    else:
        st.info("No markers have genotype data for at least two individuals.")

    # Compute mode genotype per marker for imputation
    marker_modes = {}
    for marker, group in df_valid.groupby("Marker"):
        genotypes_list = []
        for _, row in group.iterrows():
            genotype_raw = row["Genotype"]
            parsed = parse_genotype_string(genotype_raw)
            if parsed:
                genotypes_list.append(tuple(parsed))  # store as tuple for counting
        if genotypes_list:
            mode_tuple = Counter(genotypes_list).most_common(1)[0][0]
            marker_modes[marker] = list(mode_tuple)

    # Build dictionary: sample -> marker -> genotype
    genotypes = {}
    for _, row in df_valid.iterrows():
        sample = row["Sample Name"]
        marker = row["Marker"]
        genotype_raw = row["Genotype"]

        parsed = parse_genotype_string(genotype_raw)
        if not parsed:
            # fallback: mode genotype if available
            parsed = marker_modes.get(marker, [])

        genotypes.setdefault(sample, {})[marker] = parsed

    # Map from marker name to repeat unit
    marker_repeat = {}
    if "marker_list" in st.session_state and st.session_state.marker_list:
        for m in st.session_state.marker_list:
            marker_repeat[m["marker"]] = m.get("repeat_unit", 1)

    # Build Bruvo distance matrix
    samples = sorted(genotypes.keys())
    distance_matrix = pd.DataFrame(np.nan, index=samples, columns=samples)

    for i, s1 in enumerate(samples):
        for j, s2 in enumerate(samples):
            if j < i:
                continue
            common_markers = set(genotypes[s1].keys()).intersection(genotypes[s2].keys())
            if common_markers:
                d_list = []
                for marker in common_markers:
                    r_unit = marker_repeat.get(marker, 1)
                    geno1 = genotypes[s1][marker]
                    geno2 = genotypes[s2][marker]
                    d = bruvo_distance_between_genotypes(
                        geno1, geno2, r_unit,
                        genome_add=True,   # turn off if you don’t need expansions
                        genome_loss=True
                    )
                    d_list.append(d)

                if d_list:
                    avg_d = sum(d_list) / len(d_list)
                    distance_matrix.at[s1, s2] = avg_d
                    distance_matrix.at[s2, s1] = avg_d

    # Fill missing values with row-wise average, then symmetrize
    dm_filled = distance_matrix.copy()
    for i in range(len(dm_filled)):
        for j in range(len(dm_filled)):
            if i == j:
                dm_filled.iat[i, j] = 0.0
            elif pd.isnull(dm_filled.iat[i, j]):
                row_vals = dm_filled.iloc[i, :].dropna()
                dm_filled.iat[i, j] = row_vals.mean() if len(row_vals) > 0 else 0.0
    dm_filled = (dm_filled + dm_filled.T) / 2  # ensure symmetry

    # Output the matrix
    tsv_output = dm_filled.to_csv(sep="\t", index=True)
    st.write("Bruvo's distance matrix (tab-delimited):")
    st.text_area("Distance Matrix", tsv_output, height=300)

    # Perform PCoA
    D = dm_filled.values.astype(float)
    n = D.shape[0]
    D2 = D**2
    J = np.eye(n) - np.ones((n, n))/n
    B = -0.5 * J.dot(D2).dot(J)
    eigvals, eigvecs = np.linalg.eigh(B)
    idx = np.argsort(eigvals)[::-1]
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]
    # keep only positive eigenvalues
    pos_idx = eigvals > 0
    coords = eigvecs[:, :2] * np.sqrt(eigvals[:2])

    pcoa_df = pd.DataFrame(coords, columns=["PCoA1", "PCoA2"])
    pcoa_df["Sample"] = samples
    fig_pcoa = px.scatter(
        pcoa_df, x="PCoA1", y="PCoA2", text="Sample",
        title="PCoA of Bruvo's Distance Matrix",
        labels={"PCoA1": "PCoA 1", "PCoA2": "PCoA 2"}
    )
    fig_pcoa.update_traces(textposition="top center")
    st.plotly_chart(fig_pcoa)

    # Create Plotly Dendrogram
    fig = ff.create_dendrogram(dm_filled.values, labels=samples, orientation='left')
    fig.update_layout(title=f"Hierarchical Clustering Dendrogram ({method.capitalize()} Linkage)")
    st.plotly_chart(fig)

else:
    st.info("No genotype results available.")
