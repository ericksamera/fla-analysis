#!/usr/bin/env python3
__description__ =\
"""
Purpose: Streamlit wrapper for fla-viewer-advanced using preprocessed JSON files.
"""
__author__ = "Erick Samera"
__version__ = "1.2.1"
__comments__ = "stable enough; altered zoom"
# --------------------------------------------------
import streamlit as st
import pandas as pd
# --------------------------------------------------


if "marker_list" not in st.session_state:
    st.session_state.marker_list = []
if "genotype_results_df" not in st.session_state:
    st.session_state.genotype_results_df = None
if "detected_peaks_df" not in st.session_state:
    st.session_state.detected_peaks_df = None
if "PROCESSED_FLA" not in st.session_state:
    st.session_state.PROCESSED_FLA = {}
if "uploaded_files" not in st.session_state:
    st.session_state.uploaded_files = []
if "uploaded_files_id_counter" not in st.session_state:
    st.session_state.uploaded_files_id_counter = 0
if "config_changed" not in st.session_state:
    st.session_state.config_changed = False
if "marker_list" not in st.session_state:
    st.session_state.marker_list = []
if "genotype_results_df" not in st.session_state:
    st.session_state.genotype_results_df = None
if "detected_peaks_df" not in st.session_state:
    st.session_state.detected_peaks_df = None
if "PROCESSED_FLA" not in st.session_state:
    st.session_state.PROCESSED_FLA = {}
if "uploaded_files" not in st.session_state:
    st.session_state.uploaded_files = []
if "uploaded_files_id_counter" not in st.session_state:
    st.session_state.uploaded_files_id_counter = 0

if "ploidy" not in st.session_state:
    st.session_state.ploidy = 2

def main():
    st.set_page_config(
        page_title="abi-sauce | FLA-viewer",
        page_icon=':rainbow:',
        layout='wide',
        initial_sidebar_state='expanded')
    
    if "PROCESSED_FLA" not in st.session_state:
        st.session_state.PROCESSED_FLA = {}

    custom_pages = {"Analysis Tools": []}

    custom_pages["Analysis Tools"].append(
        st.Page("app_pages/home.py", title="Home", icon=":material/info:")
    )
    custom_pages["Analysis Tools"].append(
        st.Page("app_pages/File_Processor.py", title="Experiment Setup", icon=":material/file_present:")
    )

    if st.session_state.PROCESSED_FLA:
        custom_pages["Analysis Tools"].append(
            st.Page("app_pages/Peak_Visualizer.py", title="Peak Visualizer", icon=":material/insights:")
        )

        if st.session_state.marker_list and isinstance(st.session_state.genotype_results_df, pd.DataFrame):
            custom_pages["Analysis Tools"].append(
                st.Page("app_pages/Distance_Matrix.py", title="Distance Matrix", icon=":material/analytics:")
            )

    pg = st.navigation(custom_pages)
    pg.run()

    st.divider()
    st.markdown('Check out my GitHub with the link below for some of my other projects.')
    st.caption(f'[@{__author__}](https://github.com/ericksamera) | v{__version__} | {__comments__}')

    with st.sidebar:

        if (st.session_state.ploidy != 2) and st.session_state.marker_list:
            st.error(f"Running in polyploid mode ({st.session_state.ploidy}n)! Very experimental!", icon=":material/warning:")

        if st.session_state.config_changed:
            st.warning("Configuration has changed. Please re-analyze your data to apply the new settings.")
    
# --------------------------------------------------
if __name__ == "__main__":
    main()