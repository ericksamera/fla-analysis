# interface/ui_global_config.py

import streamlit as st
from fla_pipeline.config.global_config import GlobalConfig

def global_config_ui():
    if "config" not in st.session_state or st.session_state.config is None:
        st.session_state.config = GlobalConfig()

    cfg = st.session_state.config
    col1, col2, col3 = st.columns(3)

    with col1:
        cfg.ploidy = st.selectbox("Ploidy", options=[1, 2], index=1 if cfg.ploidy == 2 else 0)

    st.session_state.config = cfg
