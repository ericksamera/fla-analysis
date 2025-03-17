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
import ast
import plotly.express as px
import plotly.figure_factory as ff

# Import the modular distance calculator and MarkerConfig from the core package.
from core.genetic_distance import GeneticDistanceCalculator, ensure_unique_sample_names, perform_pcoa
from core.config_models import MarkerConfig

st.title("Hierarchical Clustering & PCoA (Plotly)")
method = "average"

# -----------------------------------------------------------------------------
# Helper functions for visualization (PCoA and dendrogram)
# -----------------------------------------------------------------------------

def create_dendrogram(distance_matrix, labels, linkage_method="single"):
    """
    Create a Plotly dendrogram from the distance matrix with the given labels.
    """
    fig = ff.create_dendrogram(distance_matrix.values, labels=labels, orientation='left')
    fig.update_layout(title=f"Hierarchical Clustering Dendrogram ({linkage_method.capitalize()} Linkage)")
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
        if len(set(list(df_valid["Sample Name"]))) != len(list(df_valid["Sample Name"])):
            st.warning("Renaming.")
        df_valid = ensure_unique_sample_names(df_valid)
        with st.expander("Genotypes Table", expanded=False):
            st.dataframe(df_valid, use_container_width=True, column_order=(column for column in df_valid.columns if column != 'Filename'))
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

    # Instantiate the GeneticDistanceCalculator using the imported module.
    calculator = GeneticDistanceCalculator(df_valid, marker_configs, metric="bruvo")
    distance_matrix = calculator.compute_distance_matrix()

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

        dm_display = distance_matrix.copy()
        if matrix_display_option == "Lower Triangle":
            dm_display = dm_display.where(np.tril(np.ones(dm_display.shape), k=0).astype(bool))
        elif matrix_display_option == "Upper Triangle":
            dm_display = dm_display.where(np.triu(np.ones(dm_display.shape), k=0).astype(bool))

        st.code(dm_display.to_csv(sep={',': ',', '\\t': '\t'}.get(raw_matrix_delimiter)), line_numbers=True, height=300)

    coords, explained_variance_ratio = perform_pcoa(distance_matrix, n_components=2)

    pos_explained = explained_variance_ratio[explained_variance_ratio > 0]

    scree_df = pd.DataFrame({
        "Principal Coordinate": [f"PCoA {i+1}" for i in range(len(pos_explained))],
        "Explained Variance": pos_explained
    })

    pcoa_df = pd.DataFrame(coords, columns=["PCoA1", "PCoA2"])
    pcoa_df["Sample"] = calculator.samples

    fig_pcoa = px.scatter(
        pcoa_df, x="PCoA1", y="PCoA2", text="Sample",
        title="PCoA of Bruvo's Distance Matrix",
        labels={"PCoA1": "PCoA 1", "PCoA2": "PCoA 2"}
    )
    fig_pcoa.update_traces(textposition="top center")

    fig_scree = px.bar(
        scree_df, x="Principal Coordinate", y="Explained Variance",
        title="Scree Plot of PCoA",
        labels={"Explained Variance": "Proportion of Variance Explained"}
    )

    col_pcoa, col_scree = st.columns([4, 2])
    with col_pcoa:
        st.plotly_chart(fig_pcoa)
    with col_scree:
        st.plotly_chart(fig_scree)

    if distance_matrix.isnull().values.any():
        st.error("Error: Distance matrix contains NaN values. Check input genotypes.")
        st.write(distance_matrix)
    fig_dend = create_dendrogram(distance_matrix, labels=calculator.samples, linkage_method=method)
    st.plotly_chart(fig_dend)
else:
    st.info("No genotype results available.")
