#!/usr/bin/env python3
__description__ =\
"""
streamlit_app.py

Purpose:
  - Serves as the main Streamlit entry point for the FLA viewer.
  - Manages navigation and session state.
  - Uses the new Pages and Navigation system.
"""
__author__ = "Erick Samera"
__version__ = "1.2.1"
__comments__ = ""

import streamlit as st
import pandas as pd

# --------------------------------------------------
# Initialize Session State
# --------------------------------------------------
session_defaults = {
    "marker_list": [],
    "genotype_results_df": pd.DataFrame(),
    "detected_peaks_df": None,
    "PROCESSED_FLA": {},
    "uploaded_files": [],
    "uploaded_files_id_counter": 0,
    "config_changed": False,
    "ploidy": 2,
    "config_json": {}
}

for key, value in session_defaults.items():
    if key not in st.session_state:
        st.session_state[key] = value

# --------------------------------------------------
# Configure Streamlit App Layout
# --------------------------------------------------
st.set_page_config(
    page_title="abi-sauce | FLA-viewer",
    page_icon=":rainbow:",
    layout="wide",
    initial_sidebar_state="expanded"
)

# --------------------------------------------------
# Navigation Setup Using Pages
# --------------------------------------------------
def main():
    custom_pages = {"Analysis Tools": []}

    # Base pages
    custom_pages["Analysis Tools"].append(
        st.Page("interface/home.py", title="Home", icon=":material/info:")
    )
    custom_pages["Analysis Tools"].append(
        st.Page("interface/file_processor.py", title="Experiment Setup", icon=":material/file_present:")
    )

    # Conditional pages
    if st.session_state.PROCESSED_FLA:
        custom_pages["Analysis Tools"].append(
            st.Page("interface/peak_visualizer.py", title="Peak Viewer", icon=":material/insights:")
        )

    if st.session_state.marker_list and not st.session_state.genotype_results_df.empty:
        custom_pages["Analysis Tools"].append(
            st.Page("interface/distance_matrix.py", title="Analysis", icon=":material/analytics:")
        )

    # Initialize navigation
    pg = st.navigation(custom_pages)
    pg.run()


    st.divider()
    st.markdown("Check out my GitHub with the link below for some of my other projects.")
    st.caption(f'[@{__author__}](https://github.com/ericksamera) | v{__version__} | {__comments__}')

    with st.sidebar:
        if st.session_state.ploidy != 2 and st.session_state.marker_list:
            st.error(f"Running in polyploid mode ({st.session_state.ploidy}n)! Very experimental!", icon=":material/warning:")

        if st.session_state.config_changed:
            st.warning("Configuration has changed. Please re-analyze your data to apply the new settings.")

# --------------------------------------------------
if __name__ == "__main__":
    main()
