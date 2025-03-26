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
import pandas as pd

import scipy.cluster.hierarchy as sch

@dataclass
class PCOAResults:
    """
    Base configuration parameters for FLA processing.
    These parameters are shared across different processor implementations.
    """
    coords: np.ndarray
    explained_variance_ratio: np.ndarray

# Import the modular distance calculator and MarkerConfig from the core package.
from core.genetic_distance import GeneticDistanceCalculator, ensure_unique_sample_names, perform_pcoa
from core.config_models import MarkerConfig

st.title("Hierarchical Clustering & PCoA (Plotly)")

@st.dialog("SSR Analysis", width="large")
def compute_distance_matrix_dialog(df_valid, marker_configs, distance_method, metric):
    with st.status("Performing analysis...", expanded=True) as status:
        with st.spinner(f"Generating distance matrix with `method='{distance_method}'` ...", show_time=True):
            calculator = GeneticDistanceCalculator(df_valid, marker_configs, metric=distance_method)
            distance_matrix: pd.DataFrame = GeneticDistanceCalculator.compute_distance_matrix(calculator)
        st.write(f"Generated `distance matrix` with `method='{distance_method}'`!")
        
        st.session_state.calculator_samples = calculator.samples
        st.session_state.distance_matrix = distance_matrix

        with st.spinner("Performing `Principle Coordinate Analysis (PCoA)` ...", show_time=True):
            pcoa_results = PCOAResults(*perform_pcoa(distance_matrix, n_components=2))
            st.session_state.pcoa_results = pcoa_results
            st.write("Performed `Principle Coordinate Analysis (PCoA)` !")

        # --- Generate Plots and store them in session_state ---
        # Scree Plot Data
        pos_explained = st.session_state.pcoa_results.explained_variance_ratio[
            st.session_state.pcoa_results.explained_variance_ratio > 0
        ]
        scree_df = pd.DataFrame({
            "Principal Coordinate": [f"PCoA {i+1}" for i in range(len(pos_explained))],
            "Explained Variance": pos_explained
        })

        # PCoA Scatter Plot
        with st.spinner("Generating figures for PCoA ...", show_time=True):
            pcoa_df = pd.DataFrame(st.session_state.pcoa_results.coords, columns=["PCoA1", "PCoA2"])
            pcoa_df["Sample"] = st.session_state.calculator_samples
            fig_pcoa = px.scatter(
                pcoa_df, x="PCoA1", y="PCoA2", text="Sample",
                title="PCoA of Bruvo's Distance Matrix",
                labels={"PCoA1": "PCoA 1", "PCoA2": "PCoA 2"}
            )
            fig_pcoa.update_traces(textposition="top center")
            st.session_state.fig_pcoa = fig_pcoa
            

            # Scree Bar Plot
            fig_scree = px.bar(
                scree_df, x="Principal Coordinate", y="Explained Variance",
                title="Scree Plot of PCoA",
                labels={"Explained Variance": "Proportion of Variance Explained"}
            )
            st.session_state.fig_scree = fig_scree
            st.write("Generated `scree plot` and `PCoA plot` !")

        # Dendrogram Plot
        with st.spinner(f"Generating `dendrogram` with `method=''{distance_method}'` ...", show_time=True):
            fig_dend = create_dendrogram(distance_matrix, labels=st.session_state.calculator_samples, linkage_method=metric)
            st.session_state.fig_dend = fig_dend
            st.write(f"Generated `dendrogram` with `method='{distance_method}'` !")

        status.update(label="Analysis complete!", state="complete", expanded=False)

# -----------------------------------------------------------------------------
# Helper functions for visualization (PCoA and dendrogram)
# -----------------------------------------------------------------------------
def convert_to_genalex_csv(df):
    """
    Converts a genotype DataFrame into a GenAlEx-formatted CSV string.
    Assumes the input DataFrame df has at least the following columns:
      - "Sample Name" : sample identifier
      - "Marker"      : locus name
      - "Genotype"    : a tuple/list (or string that can be evaluated) containing two allele values
      - "Population"  : (optional) population assignment; if missing, defaults to "Unknown"
    
    The GenAlEx CSV format will be:
      - First row: metadata row (blank, num_samples, num_loci, then blanks)
      - Second row: header row (Sample No., Pop, then for each locus: Marker_1, Marker_2)
      - Subsequent rows: one row per sample.
    """
    # Get sorted unique markers
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
                if not row_data["Genotype"].iloc[0]: row.extend(["", ""])
                else:
                    if not row_data["Genotype"].iloc[0][0]: row.extend(["", ""])
                    else:
                        genotype = row_data["Genotype"].iloc[0][0]
                        if isinstance(genotype, (list, tuple)) and len(genotype) >= 2:
                            row.extend([*genotype])
        csv_like.append(row)

    
    # Write to a CSV string using tab as delimiter (or comma if you prefer)
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
    # Read CSV with tab delimiter and use the second row as header
    return pd.read_csv(io.StringIO(csv_text), sep=",")

def create_dendrogram(distance_matrix, labels, linkage_method="single"):
    """
    Create a Plotly dendrogram from the distance matrix with the given labels.
    """
    fig = ff.create_dendrogram(distance_matrix.values, labels=labels, linkagefun=lambda x: sch.linkage(x, "average"), orientation='left')
    fig.update_layout(title=f"Hierarchical Clustering Dendrogram (UPGMA) Linkage)")
    return fig

# -----------------------------------------------------------------------------
# Main Streamlit Logic
# -----------------------------------------------------------------------------
if not st.session_state.genotype_results_df.empty:
    df = st.session_state.genotype_results_df.copy()

    # Filter out markers that don't have genotype data for at least two individuals.
    marker_counts = df.groupby("Marker").size()
    valid_markers = marker_counts[marker_counts >= 2].index.tolist()

    if valid_markers:
        df_valid = df[df["Marker"].isin(valid_markers)]
        if len(set(list(df_valid["Filename"]))) != len(set(list(df_valid["Sample Name"]))):
            st.warning("Sample naming not unique! Appended a number to make sample names unique. Consider renaming samples in the `Experiment Setup`.")
        df_valid = ensure_unique_sample_names(df_valid)

        data_format = st.selectbox("Select Data Format", ["Standard Table", "GenAlEx Format"])

        if data_format == "Standard Table":
            st.dataframe(df_valid, use_container_width=True, column_order=(column for column in df_valid.columns if column != 'Filename'))
        elif data_format == "GenAlEx Format":
            # Convert the original DataFrame to a GenAlEx-formatted CSV string
            genalex_csv_text = convert_to_genalex_csv(df_valid)
            # Parse the CSV string back into a DataFrame (ignoring the metadata row if desired)
            genalex_df = get_genalex_df_from_csv(genalex_csv_text)
            st.dataframe(genalex_df, use_container_width=True, hide_index=True)
            
            # Remove the first row from the CSV text for download
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

    # Create MarkerConfig objects from the marker list stored in session_state.
    marker_list = st.session_state.get("marker_list", [])
    if marker_list:
        marker_configs = [MarkerConfig(**m) for m in marker_list]
    else:
        st.error("No marker configuration found in session state.")
        st.stop()
    
    # Guard clause: Only compute if PCoA hasn't been generated yet.
    if not st.session_state.pcoa_results:
        compute_distance_matrix_dialog(df_valid, marker_configs, distance_method="bruvo", metric="bruvo")


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
