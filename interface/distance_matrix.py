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
        with st.expander("Genotypes Table", expanded=False):
            st.dataframe(df_valid, use_container_width=True,
                         column_order=(column for column in df_valid.columns if column != 'Filename'))
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
