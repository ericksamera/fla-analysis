# streamlit_app.py

import streamlit as st
from interface.backend.session import initialize_session_state
from interface.backend.session_io import session_export_button, session_import_button


initialize_session_state()

st.set_page_config(
    page_title="abi-sauce | FLA viewer",
    page_icon="🌈",
    layout="wide",
    initial_sidebar_state="expanded",
)

def main():
    initialize_session_state()

    with st.sidebar:
        col_export, col_import  = st.columns(2, border=False)
        with col_import:
            session_export_button()
        with col_export:
            session_import_button()

    custom_pages = {"Analysis Tools": []}
    custom_pages["Analysis Tools"].append(
        st.Page("interface/home.py", title="Home", icon=":material/info:")
    )

    custom_pages["Analysis Tools"].append(
        st.Page("interface/file_processor.py", title="Experiment Setup", icon=":material/file_present:")
    )

    if st.session_state.samples:
        custom_pages["Analysis Tools"].append(
            st.Page("interface/sample_metadata.py", title="Sample Metadata Editor", icon=":material/view_list:")
        )

        custom_pages["Analysis Tools"].append(
            st.Page("interface/peak_visualizer.py", title="Peak Viewer", icon=":material/insights:")
        )

        if not st.session_state.genotype_results_df.empty:
            custom_pages["Analysis Tools"].append(
                st.Page("interface/population_analysis.py", title="Analysis", icon=":material/analytics:")
            )

    # Register pages
    page = st.navigation(custom_pages)
    page.run()

    st.divider()
    st.caption("abi-sauce | Developed by Erick Samera")

if __name__ == "__main__":
    main()