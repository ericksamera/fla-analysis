# interface/population_analysis.py

import streamlit as st
import pandas as pd
from fla_pipeline.analysis.config import DistanceConfig
from fla_pipeline.analysis.distance import DistanceCalculator
from fla_pipeline.config import MarkerConfig
from fla_pipeline.models.genotype import GenotypeResult
from interface.plotting.population_plots import plot_pcoa

def run():
    st.title("📊 Population Structure Analysis")

    if "genotype_results_df" not in st.session_state or st.session_state.genotype_results_df.empty:
        st.info("No genotype results available. Please process samples first.")
        return

    samples = st.session_state.samples
    marker_list = st.session_state.marker_list
    genotype_matrix = {
        sid: sample.marker_results
        for sid, sample in samples.items()
        if sample.marker_results
    }
    marker_configs = {
        m["marker"]: MarkerConfig(**m) for m in marker_list
    }
    metadata = {
        sid: sample.metadata for sid, sample in samples.items()
    }

    # Advanced settings
    with st.expander("⚙️ Advanced Settings", expanded=False):
        metric = st.selectbox("Distance metric", options=["bruvo"], index=0)
        min_conf = st.slider("Min genotype confidence", 0.0, 1.0, 0.8, step=0.05)
        min_samples = st.slider("Min sample count per marker", 1, 10, 2)
        impute = st.checkbox("Enable imputation for missing genotypes", value=False)
        debug_mode = st.checkbox("Enable debug logging", value=False)

    # Recalc button
    if st.button("🚀 Recalculate Distance Matrix"):
        cfg = DistanceConfig(
            metric=metric,
            min_confidence=min_conf,
            min_sample_count=min_samples,
            impute_missing=impute,
            debug_mode=debug_mode
        )
        calc = DistanceCalculator(genotype_matrix, marker_configs, sample_metadata=metadata)
        matrix = calc.compute_matrix(cfg)
        st.session_state["distance_matrix"] = matrix
        st.session_state["distance_debug_log"] = calc.get_debug_log()
        st.session_state["included_markers"] = calc.get_included_markers()

    if "distance_matrix" not in st.session_state:
        st.stop()

    matrix = st.session_state["distance_matrix"]
    st.subheader("🔢 Distance Matrix")
    st.dataframe(matrix, height=300, use_container_width=True)

    # PCoA
    calc = DistanceCalculator(genotype_matrix, marker_configs, sample_metadata=metadata)
    coords_array, explained = calc.perform_pcoa(matrix)

    coords_df = pd.DataFrame(coords_array, columns=["PCoA1", "PCoA2"])
    coords_df.index = matrix.index

    # --- Color-by metadata ---
    meta_keys = sorted({"Population", "Sample Name"} | {k for sample in samples.values() for k in sample.metadata.keys()})
    color_by = st.selectbox("Color PCoA by metadata field", options=["None"] + meta_keys, index=0)
    color_field = color_by if color_by != "None" else None

    fig = plot_pcoa(coords_df, metadata=metadata, color_by=color_field)
    st.plotly_chart(fig, use_container_width=True)

    # Marker summary
    st.subheader("🧬 Included Markers")
    with st.expander("Marker filter results", expanded=False):
        included = st.session_state.get("included_markers", [])
        excluded = sorted(set(marker_configs.keys()) - set(included))
        st.write(f"**Included markers ({len(included)}):** {', '.join(included)}")
        st.write(f"**Excluded markers ({len(excluded)}):** {', '.join(excluded)}")

    # Debug log
    if debug_mode:
        with st.expander("📜 Debug Log"):
            st.text_area("Distance Log", st.session_state.get("distance_debug_log", ""), height=300)

run()
