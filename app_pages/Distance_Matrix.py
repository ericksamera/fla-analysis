import streamlit as st
import pandas as pd
import numpy as np
import itertools
import ast
import plotly.express as px
import plotly.figure_factory as ff
import plotly.graph_objects as go
from collections import Counter

if st.session_state.config_changed:
    st.warning("Configuration has changed. Please re-analyze your data to apply the new settings.")

st.title("Hierarchical Clustering & PCoA (Plotly)")
method = st.selectbox(
    "Select linkage method for hierarchical clustering",
    ["single", "complete", "average", "weighted", "centroid", "median", "ward"]
)

###############################################################################
# Utility Functions (Same as Before)
###############################################################################
def parse_genotype_string(genotype_str):
    """Parse a genotype string robustly into a list of floats."""
    if not isinstance(genotype_str, str):
        return flatten_if_nested(genotype_str)
    try:
        obj = ast.literal_eval(genotype_str)  # e.g. "[228,231]", "(228, 231)"
        return flatten_if_nested(obj)
    except:
        pass
    try:
        return [float(x.strip()) for x in genotype_str.split(",") if x.strip()]
    except:
        pass
    try:
        return [float(x.strip()) for x in genotype_str.split() if x.strip()]
    except:
        pass
    return []

def flatten_if_nested(thing):
    """Flatten nested list/tuple elements into a 1D list of floats."""
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
    try:
        return [float(thing)]
    except:
        return []

def to_repeat_count(allele_size, repeat_unit):
    """Convert allele size to its approximate repeat count."""
    if repeat_unit <= 0:
        repeat_unit = 1e-5
    return int(round(allele_size / repeat_unit))

def _min_perm_distance(gcounts1, gcounts2):
    """Compute minimal average Bruvo distance across permutations of gcounts2."""
    best = np.inf
    for perm in itertools.permutations(gcounts2):
        distances = []
        for a_val, b_val in zip(gcounts1, perm):
            x = abs(a_val - b_val)
            d = 1 - 2 ** (-x)  # Bruvo's distance formula
            distances.append(d)
        avg = sum(distances)/len(distances) if distances else 1.0
        best = min(best, avg)
    return best if best != np.inf else 1.0

def bruvo_distance_between_genotypes(geno1_raw, geno2_raw, repeat_unit,
                                     genome_add=True, genome_loss=True):
    """Bruvo’s distance metric for two genotypes."""
    if not geno1_raw or not geno2_raw:
        return 1.0
    geno1_counts = [to_repeat_count(a, repeat_unit) for a in geno1_raw]
    geno2_counts = [to_repeat_count(a, repeat_unit) for a in geno2_raw]
    n1, n2 = len(geno1_counts), len(geno2_counts)
    if n1 == 0 or n2 == 0:
        return 1.0

    if genome_add or genome_loss:
        if n1 < n2:
            expansions = itertools.product(geno1_counts, repeat=(n2 - n1))
            best_score = np.inf
            for exp in expansions:
                g1_full = list(geno1_counts) + list(exp)
                score = _min_perm_distance(g1_full, geno2_counts)
                best_score = min(best_score, score)
            return best_score
        elif n2 < n1:
            expansions = itertools.product(geno2_counts, repeat=(n1 - n2))
            best_score = np.inf
            for exp in expansions:
                g2_full = list(geno2_counts) + list(exp)
                score = _min_perm_distance(geno1_counts, g2_full)
                best_score = min(best_score, score)
            return best_score
        else:
            # equal length
            return _min_perm_distance(geno1_counts, geno2_counts)
    else:
        # no genome add/loss
        max_len = max(n1, n2)
        if n1 < max_len:
            geno1_counts *= max_len
        if n2 < max_len:
            geno2_counts *= max_len
        return _min_perm_distance(geno1_counts, geno2_counts)

###############################################################################
# Class-Based Approach
###############################################################################
class GenotypeAnalyzer:
    def __init__(self, df_valid, marker_repeat):
        """
        df_valid: A filtered DataFrame containing valid markers for each sample.
        marker_repeat: Dict of {marker -> repeat_unit}, e.g. {"M1": 2, "M2": 4}.
        """
        self.df_valid = df_valid
        self.marker_repeat = marker_repeat
        
        # Parse + Impute
        self.genotypes, self.marker_modes = self._parse_and_impute_genotypes(df_valid)
        self.samples = sorted(self.genotypes.keys())

    def _parse_and_impute_genotypes(self, df_valid):
        """Create {sample -> {marker -> genotype_list}} and {marker -> mode_genotype_list}."""
        marker_modes = {}
        for marker, group in df_valid.groupby("Marker"):
            genotypes_list = []
            for _, row in group.iterrows():
                parsed = parse_genotype_string(row["Genotype"])
                if parsed:
                    genotypes_list.append(tuple(parsed))
            if genotypes_list:
                mode_tuple = Counter(genotypes_list).most_common(1)[0][0]
                marker_modes[marker] = list(mode_tuple)

        genotypes = {}
        for _, row in df_valid.iterrows():
            sample = row["Sample Name"]
            marker = row["Marker"]
            parsed = parse_genotype_string(row["Genotype"])
            if not parsed:
                parsed = marker_modes.get(marker, [])
            genotypes.setdefault(sample, {})[marker] = parsed
        return genotypes, marker_modes

    def compute_distance_matrix(self, distance_func=bruvo_distance_between_genotypes):
        """Build a distance matrix using a chosen distance function."""
        dist_matrix = pd.DataFrame(np.nan, index=self.samples, columns=self.samples)
        
        for i, s1 in enumerate(self.samples):
            for j, s2 in enumerate(self.samples):
                if j < i:
                    continue
                common_markers = set(self.genotypes[s1].keys()) & set(self.genotypes[s2].keys())
                if common_markers:
                    d_list = []
                    for marker in common_markers:
                        r_unit = self.marker_repeat.get(marker, 1)
                        geno1 = self.genotypes[s1][marker]
                        geno2 = self.genotypes[s2][marker]
                        d_list.append(distance_func(geno1, geno2, r_unit))
                    avg_d = sum(d_list) / len(d_list) if d_list else 1.0
                    dist_matrix.at[s1, s2] = avg_d
                    dist_matrix.at[s2, s1] = avg_d

        # Fill diagonals w/ 0, symmetrize missing
        dm_filled = dist_matrix.copy()
        for i in range(len(dm_filled)):
            for j in range(len(dm_filled)):
                if i == j:
                    dm_filled.iat[i, j] = 0.0
                elif pd.isnull(dm_filled.iat[i, j]):
                    row_vals = dm_filled.iloc[i, :].dropna()
                    dm_filled.iat[i, j] = row_vals.mean() if len(row_vals) > 0 else 0.0
        dm_filled = (dm_filled + dm_filled.T) / 2
        return dm_filled

    def perform_pcoa(self, distance_matrix, n_components=2):
        """Perform basic PCoA on a distance matrix."""
        D = distance_matrix.values.astype(float)
        n = D.shape[0]
        D2 = D ** 2
        J = np.eye(n) - np.ones((n, n)) / n
        B = -0.5 * J.dot(D2).dot(J)
        eigvals, eigvecs = np.linalg.eigh(B)
        idx = np.argsort(eigvals)[::-1]
        eigvals = eigvals[idx]
        eigvecs = eigvecs[:, idx]
        pos_idx = eigvals > 0
        coords = eigvecs[:, pos_idx][:, :n_components] * np.sqrt(eigvals[pos_idx][:n_components])
        return coords

    def create_dendrogram(self, distance_matrix, linkage_method="single"):
        """Create a Plotly dendrogram from the distance matrix."""
        fig = ff.create_dendrogram(distance_matrix.values, labels=self.samples, orientation='left')
        fig.update_layout(title=f"Hierarchical Clustering Dendrogram ({linkage_method.capitalize()} Linkage)")
        return fig


###############################################################################
# 7) Main Streamlit Logic
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
        st.stop()

    # Build marker_repeat from the session_state if available
    marker_repeat = {}
    if "marker_list" in st.session_state and st.session_state.marker_list:
        for m in st.session_state.marker_list:
            marker_repeat[m["marker"]] = m.get("repeat_unit", 1)

    # Create the analyzer
    analyzer = GenotypeAnalyzer(df_valid, marker_repeat)

    # Compute the distance matrix (using Bruvo by default)
    dm_filled = analyzer.compute_distance_matrix()

    # Show distance matrix
    with st.expander("Raw distance matrix", expanded=False):
        tsv_output = dm_filled.to_csv(sep="\t", index=True)
        st.write("Bruvo's distance matrix (tab-delimited):")
        st.text_area("Distance Matrix", tsv_output, height=300)

    # Perform PCoA
    coords = analyzer.perform_pcoa(dm_filled, n_components=2)
    pcoa_df = pd.DataFrame(coords, columns=["PCoA1", "PCoA2"])
    pcoa_df["Sample"] = analyzer.samples

    fig_pcoa = px.scatter(
        pcoa_df, x="PCoA1", y="PCoA2", text="Sample",
        title="PCoA of Bruvo's Distance Matrix",
        labels={"PCoA1": "PCoA 1", "PCoA2": "PCoA 2"}
    )
    fig_pcoa.update_traces(textposition="top center")
    st.plotly_chart(fig_pcoa)

    # Dendrogram
    fig_dend = analyzer.create_dendrogram(dm_filled, linkage_method=method)
    st.plotly_chart(fig_dend)

else:
    st.info("No genotype results available.")
