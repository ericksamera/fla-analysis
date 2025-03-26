# streamlit_app.py

import streamlit as st
from interface.backend.session import initialize_session_state

# -------------------------------------
# 🔧 Session Initialization
# -------------------------------------
initialize_session_state()

# -------------------------------------
# ⚙️ Streamlit App Setup
# -------------------------------------
st.set_page_config(
    page_title="abi-sauce | FLA viewer",
    page_icon="🌈",
    layout="wide",
    initial_sidebar_state="expanded",
)

# -------------------------------------
# 📚 Custom Navigation (No pages/ dir)
# -------------------------------------
def main():
    initialize_session_state()

    custom_pages = {"FLA Analysis": []}

    custom_pages["FLA Analysis"].append(
        st.Page("interface/home.py", title="Home", icon="🏠")
    )

    custom_pages["FLA Analysis"].append(
        st.Page("interface/file_processor.py", title="Upload + Setup", icon="📁")
    )

    if st.session_state.samples:
        custom_pages["FLA Analysis"].append(
            st.Page("interface/peak_visualizer.py", title="Peak Viewer", icon="📈")
        )

    # if not st.session_state.genotype_results_df.empty:
    #     custom_pages["FLA Analysis"].append(
    #         st.Page("interface/distance_matrix.py", title="Analysis", icon="📊")
    #     )

    # Register pages
    page = st.navigation(custom_pages)
    page.run()

    st.divider()
    st.caption("abi-sauce | Developed by Erick Samera")

if __name__ == "__main__":
    main()
