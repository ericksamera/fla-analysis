# interface/backend/session.py

import streamlit as st
import pandas as pd

def initialize_session_state():
    defaults = {
        "samples": {},                     # {sample_id: Sample}
        "marker_list": [],                 # list of marker dicts
        "config": None,                    # GlobalConfig instance
        "uploaded_files": [],              # uploaded file metadata
        "config_changed": False,           # UI toggle when markers change
        "selected_sample": None,           # key into st.session_state.samples
        "ploidy": 2,                       # global ploidy setting
        "genotype_results_df": pd.DataFrame(),  # merged output
    }
    for key, value in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = value
