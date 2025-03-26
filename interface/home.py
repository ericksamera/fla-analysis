# interface/home.py

import streamlit as st

def run():
    st.title("🌈 Welcome to `abi-sauce`")
    st.markdown(
        """
        This app is designed for visualizing and analyzing **Fragment Length Analysis (FLA)** data.

        **Features:**
        - Upload `.fsa` trace files
        - Set marker bin ranges and repeat units
        - Automatically call diploid genotypes with QC
        - Visualize peaks and genotype results
        - Perform clustering & population analysis (PCoA, dendrograms)

        **Next step:** Go to the **Upload + Setup** page to begin.
        """
    )

run()