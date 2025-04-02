# interface/ui_global_config.py

import streamlit as st
from fla_pipeline.config.global_config import GlobalConfig

def instantiate_global_config():
    if "config" not in st.session_state or st.session_state.config is None:
        st.session_state.config = GlobalConfig()


def global_config_ui():
    with st.expander("Advanced options"):
        cfg = st.session_state.config
        
        col1, col2 = st.columns([3, 5])

        with col1:
            ploidy_labels = {
                1: "Haploid (n)",
                2: "Diploid (2n)",
            }
            default_label = ploidy_labels.get(cfg.ploidy, "Polyploid (> 2n)")
            options = list(ploidy_labels.values()) + ["Polyploid (> 2n)"]

            selected_label = st.radio("Ploidy", options=options, index=options.index(default_label), horizontal=True)

        with col2:
            if selected_label == "Polyploid (> 2n)":
                cfg.ploidy = st.slider("Polyploidy (> 2n)", min_value=3, max_value=8, value=max(cfg.ploidy, 3))
            else:
                cfg.ploidy = {v: k for k, v in ploidy_labels.items()}[selected_label]

        cfg.min_peak_height = st.number_input("Minimum Peak Height", value=1000 if not cfg.min_peak_height else cfg.min_peak_height)

        st.session_state.config = cfg