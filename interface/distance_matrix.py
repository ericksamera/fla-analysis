#!/usr/bin/env python3
"""
interface/distance_matrix.py

This module implements a page for visualizing hierarchical clustering and PCoA 
using Plotly. Genotype results are stored in st.session_state. The genetic 
distance matrix is now computed using a modular approach by importing 
the GeneticDistanceCalculator from the core package.
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.figure_factory as ff
from dataclasses import dataclass
import io
import csv
import scipy.cluster.hierarchy as sch

# --------------------------------------------------------------------------
# Additional: Rst Distance Metric for microsatellite analysis
# --------------------------------------------------------------------------
class RstDistanceMetric:
    """
    Computes a simple Rst-like genetic distance for microsatellite data.
    This implementation calculates the average squared difference between allele sizes,
    normalized by the squared range of all allele sizes.
    """
    def compute_distance(self, geno1, geno2, marker, marker_repeat):
        def flatten(geno):
            if isinstance(geno, (list, tuple)):
                flat = []
                for g in geno:
                    if isinstance(g, (list, tuple)):
                        flat.extend(g)
                    else:
                        flat.append(g)
                return flat
            return [geno]
        alleles1 = flatten(geno1)
        alleles2 = flatten(geno2)
        if not alleles1 or not alleles2:
            return 1.0
        try:
            alleles1 = [float(a) for a in alleles1]
            alleles2 = [float(a) for a in alleles2]
        except Exception:
            return 1.0
        all_alleles = alleles1 + alleles2
        range_sq = (max(all_alleles) - min(all_alleles)) ** 2
        if range_sq == 0:
            range_sq = 1
        squared_diffs = []
        for a in alleles1:
            for b in alleles2:
                squared_diffs.append((a - b) ** 2)
        avg_sq_diff = sum(squared_diffs) / len(squared_diffs)
        rst = avg_sq_diff / range_sq
        return min(max(rst, 0), 1)

# --------------------------------------------------------------------------
# Imports from Core Modules
# --------------------------------------------------------------------------
from core.genetic_distance import GeneticDistanceCalculator, ensure_unique_sample_names, perform_pcoa
from core.config_models import MarkerConfig

# --------------------------------------------------------------------------
# UI Title and Metric Selection
# --------------------------------------------------------------------------
st.title("Hierarchical Clustering & PCoA (Plotly)")

# Allow the user to choose the genetic distance metric.
selected_metric = st.radio("Select Genetic Distance Metric", ["bruvo", "nei", "rst"], horizontal=True)

@dataclass
class PCOAResults:
    """
    Stores PCoA coordinates and explained variance ratios.
    """
    coords: np.ndarray
    explained_variance_ratio: np.ndarray

@st.dialog("SSR Analysis", width="large")
def compute_distance_matrix_dialog(df_valid, marker_configs, distance_method, metric):
    with st.status("Performing analysis...", expanded=True) as status:
        with st.spinner(f"Generating distance matrix with `method='{distance_method}'` ...", show_time=True):
            # If "rst" is selected, pass an instance of RstDistanceMetric; otherwise use string.
            metric_param = RstDistanceMetric() if metric == "rst" else metric
            calculator = GeneticDistanceCalculator(df_valid, marker_configs, metric=metric_param)
            distance_matrix: pd.DataFrame = GeneticDistanceCalculator.compute_distance_matrix(calculator)
        st.write(f"Generated `distance matrix` with `method='{distance_method}'`!")
        
        st.session_state.calculator_samples = calculator.samples
        st.session_state.distance_matrix = distance_matrix

        with st.spinner("Performing `Principle Coordinate Analysis (PCoA)` ...", show_time=True):
            pcoa_results = PCOAResults(*perform_pcoa(distance_matrix, n_components=2))
            st.session_state.pcoa_results = pcoa_results
            st.write("Performed `Principle Coordinate Analysis (PCoA)` !")

        # --- Generate Scree and PCoA Plots ---
        pos_explained = st.session_state.pcoa_results.explained_variance_ratio[
            st.session_state.pcoa_results.explained_variance_ratio > 0
        ]
        scree_df = pd.DataFrame({
            "Principal Coordinate": [f"PCoA {i+1}" for i in range(len(pos_explained))],
            "Explained Variance": pos_explained
        })

        with st.spinner("Generating figures for PCoA ...", show_time=True):
            pcoa_df = pd.DataFrame(st.session_state.pcoa_results.coords, columns=["PCoA1", "PCoA2"])
            pcoa_df["Sample"] = st.session_state.calculator_samples
            fig_pcoa = px.scatter(
                pcoa_df, x="PCoA1", y="PCoA2", text="Sample",
                title=f"PCoA of {distance_method.capitalize()}'s Distance Matrix",
                labels={"PCoA1": "PCoA 1", "PCoA2": "PCoA 2"}
            )
            fig_pcoa.update_traces(textposition="top center")
            st.session_state.fig_pcoa = fig_pcoa

            fig_scree = px.bar(
                scree_df, x="Principal Coordinate", y="Explained Variance",
                title="Scree Plot of PCoA",
                labels={"Explained Variance": "Proportion of Variance Explained"}
            )
            st.session_state.fig_scree = fig_scree
            st.write("Generated `scree plot` and `PCoA plot` !")

        with st.spinner(f"Generating `dendrogram` with `method=''{distance_method}'` ...", show_time=True):
            fig_dend = create_dendrogram(distance_matrix, labels=st.session_state.calculator_samples, linkage_method=metric)
            st.session_state.fig_dend = fig_dend
            st.write(f"Generated `dendrogram` with `method='{distance_method}'` !")

        status.update(label="Analysis complete!", state="complete", expanded=False)

# --------------------------------------------------------------------------
# GenAlEx Conversion Functions
# --------------------------------------------------------------------------
def convert_to_genalex_csv(df):
    """
    Converts a genotype DataFrame into a GenAlEx-formatted CSV string.
    Assumes the input DataFrame has at least:
      - "Sample Name": sample identifier
      - "Marker": locus name
      - "Genotype": a tuple/list (or evaluable string) containing allele values
      - "Population": (optional) population assignment (defaults to "Unknown")
    
    Format:
      - First row: metadata (blank, num_loci, num_samples, etc.)
      - Second row: header (Sample No., Pop, then markers)
      - Subsequent rows: one per sample.
    """
    unique_markers = sorted(df["Marker"].unique())
    unique_samples = df["Sample Name"].unique()
    optional_title = "TITLE"
    csv_like = []
    params_row = f"{len(unique_markers)},{len(unique_samples)},1,{len(unique_samples)}"
    a2_row = f"{optional_title},,,POP1,"
    header = ["SAMPLE", "POP"]
    for marker in unique_markers:
        header += [marker, ""]
    header_row = ','.join(header)
    csv_like += [params_row.split(',')]
    csv_like += [a2_row.split(',')]
    csv_like += [header_row.split(',')]
    for sample in unique_samples:
        sample_df = df[df["Sample Name"] == sample]
        row = [sample, "1"]
        for marker in unique_markers:
            row_data = sample_df[sample_df["Marker"] == marker]
            if not row_data.empty:
                if not row_data["Genotype"].iloc[0]:
                    row.extend(["", ""])
                else:
                    if not row_data["Genotype"].iloc[0][0]:
                        row.extend(["", ""])
                    else:
                        genotype = row_data["Genotype"].iloc[0][0]
                        if isinstance(genotype, (list, tuple)) and len(genotype) >= 2:
                            row.extend([*genotype])
        csv_like.append(row)
    output = io.StringIO()
    writer = csv.writer(output, delimiter=",")
    blank_header_row = ["\u200b"*(i+1) for i in range(1 + (2 + (2 * len(unique_markers))))]
    writer.writerow(blank_header_row)
    for row in csv_like:
        writer.writerow([i if i else 0 for i in row])
    csv_text = output.getvalue()
    output.close()
    return csv_text

def get_genalex_df_from_csv(csv_text):
    return pd.read_csv(io.StringIO(csv_text), sep=",")

# --------------------------------------------------------------------------
# Dendrogram Helper Function
# --------------------------------------------------------------------------
def create_dendrogram(distance_matrix, labels, linkage_method="single"):
    fig = ff.create_dendrogram(distance_matrix.values, labels=labels, linkagefun=lambda x: sch.linkage(x, "average"), orientation='left')
    fig.update_layout(title="Hierarchical Clustering Dendrogram (UPGMA) Linkage)")
    return fig

# --------------------------------------------------------------------------
# Main Streamlit Logic
# --------------------------------------------------------------------------
if not st.session_state.genotype_results_df.empty:
    df = st.session_state.genotype_results_df.copy()
    
    # Filter out markers with genotype data in fewer than two individuals.
    marker_counts = df.groupby("Marker").size()
    valid_markers = marker_counts[marker_counts >= 2].index.tolist()
    
    if valid_markers:
        df_valid = df[df["Marker"].isin(valid_markers)]
        if len(set(df_valid["Filename"])) != len(set(df_valid["Sample Name"])):
            st.warning("Sample naming not unique! Appended a number to make sample names unique. Consider renaming samples in the `Experiment Setup`.")
        df_valid = ensure_unique_sample_names(df_valid)
        
        # Let the user select the data format for display.
        data_format = st.selectbox("Select Data Format", ["Standard Table", "GenAlEx Format"])
        if data_format == "Standard Table":
            st.dataframe(df_valid, use_container_width=True, column_order=(column for column in df_valid.columns if column != 'Filename'))
        elif data_format == "GenAlEx Format":
            genalex_csv_text = convert_to_genalex_csv(df_valid)
            genalex_df = get_genalex_df_from_csv(genalex_csv_text)
            st.dataframe(genalex_df, use_container_width=True, hide_index=True)
            # Remove metadata row for download
            download_csv_text = "\n".join(genalex_csv_text.splitlines()[1:])
            st.download_button(
                label="Download GenAlEx File", 
                data=download_csv_text, 
                file_name="genalex.csv", 
                mime="text/csv"
            )
    else:
        st.info("No markers have genotype data for at least two individuals.")
        st.stop()
    
    # Create MarkerConfig objects from the session marker list.
    marker_list = st.session_state.get("marker_list", [])
    if marker_list:
        marker_configs = [MarkerConfig(**m) for m in marker_list]
    else:
        st.error("No marker configuration found in session state.")
        st.stop()
    
    # Compute the distance matrix, PCoA, and dendrogram if not already generated.
    if not st.session_state.pcoa_results:
        # Here we use "bruvo" as default; you could extend to use the selected_metric.
        compute_distance_matrix_dialog(df_valid, marker_configs, distance_method=selected_metric, metric=selected_metric)
    
    with st.expander("Raw Distance Matrix", expanded=False):
        raw_matrix_delimiter = st.radio(
            "Table Delimiter",
            [",", "\\t"],
            captions=["(.csv)", "(.tsv)"],
            format_func=lambda x: f'`{x}`',
            horizontal=True
        )
        matrix_display_option = st.radio(
            "Matrix Display",
            ["Full Matrix", "Lower Triangle", "Upper Triangle"],
            horizontal=True
        )
        dm_display = st.session_state.distance_matrix.copy()
        if matrix_display_option == "Lower Triangle":
            dm_display = dm_display.where(np.tril(np.ones(dm_display.shape), k=0).astype(bool))
        elif matrix_display_option == "Upper Triangle":
            dm_display = dm_display.where(np.triu(np.ones(dm_display.shape), k=0).astype(bool))
        st.code(dm_display.to_csv(sep={',': ',', '\\t': '\t'}.get(raw_matrix_delimiter)),
                line_numbers=True, height=300)
    
    col_pcoa, col_scree = st.columns([4, 2])
    with col_pcoa:
        st.plotly_chart(st.session_state.fig_pcoa)
    with col_scree:
        st.plotly_chart(st.session_state.fig_scree)
    
    if st.session_state.distance_matrix.isnull().values.any():
        st.error("Error: Distance matrix contains NaN values. Check input genotypes.")
        st.write(st.session_state.distance_matrix)
    st.plotly_chart(st.session_state.fig_dend)
else:
    st.info("No genotype results available.")
