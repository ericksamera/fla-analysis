import streamlit as st
import pandas as pd

from fla_pipeline.analysis.config import DistanceConfig
from fla_pipeline.analysis.distance import DistanceCalculator
from fla_pipeline.analysis.fst_analysis import calculate_fst_per_marker
from fla_pipeline.analysis.hwe_analysis import compute_hwe_stats, display_hwe_results
from fla_pipeline.utils.table_builders import build_genotype_results_df
from fla_pipeline.utils.exporters import get_exporter
from interface.plotting.population_plots import plot_pcoa


def export_section(samples, genotype_df):
    format_option = st.selectbox("Choose export format", options=["GenAlEx"])
    metadata_by_uid = {
        sample.sample_uid: {
            "Sample Name": sample.metadata.get("Sample Name", sample.sample_id),
            "Population": sample.metadata.get("Population", "Unknown"),
        } for sample in samples.values()
    }
    exporter = get_exporter(format_option, metadata=metadata_by_uid)
    csv_data = exporter.export(genotype_df)
    st.download_button("Download Export", data=csv_data, file_name=exporter.filename(), mime="text/csv")


def fst_section(genotype_df, marker_configs):
    fst_df = calculate_fst_per_marker(genotype_df)
    with st.expander('FST stuff'):
        st.subheader("FST Per Marker")
        st.dataframe(fst_df, use_container_width=True)
        st.markdown("**Top divergent markers:**")
        st.write(fst_df[fst_df['FST'] > 0.15])


def distance_settings():
    with st.expander("Advanced Settings", expanded=False):
        metric = st.selectbox("Distance metric", options=["bruvo", "nei"], index=0)
        min_conf = st.slider("Min genotype confidence", 0.0, 1.0, 0.8, step=0.05)
        min_samples = st.slider("Min sample count per marker", 1, 10, 2)
        impute = st.checkbox("Enable imputation for missing genotypes", value=False)
        debug_mode = st.checkbox("Enable debug logging", value=False)
    return metric, min_conf, min_samples, impute, debug_mode


def run():
    st.title("Population Structure Analysis")
    samples = st.session_state.samples

    if not samples or not any(s.marker_results for s in samples.values()):
        st.info("No genotype results available. Please process samples first.")
        return

    genotype_df = build_genotype_results_df(samples)
    marker_list = st.session_state.marker_list
    genotype_matrix = {sid: s.marker_results for sid, s in samples.items() if s.marker_results}
    marker_configs = {m.marker: m for m in marker_list}
    metadata = {sid: s.metadata for sid, s in samples.items()}

    genotype_df["_Pop"] = genotype_df["sample_uid"].map(
        lambda uid: samples[uid].metadata.get("Population", "Unknown")
        if uid in samples else "Unknown"
    )


    hwe_results = compute_hwe_stats(genotype_df, marker_configs)
    
    

    metric, min_conf, min_samples, impute, debug_mode = distance_settings()

    if st.button("Recalculate Distance Matrix"):
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

    with st.expander("Distance Matrix"):
        matrix = st.session_state["distance_matrix"]
        st.subheader("Distance Matrix")
        st.dataframe(matrix, height=300, use_container_width=True)

        export_section(samples, genotype_df)

    calc = DistanceCalculator(genotype_matrix, marker_configs, sample_metadata=metadata)
    coords_array, explained = calc.perform_pcoa(matrix)
    coords_df = pd.DataFrame(coords_array, columns=["PCoA1", "PCoA2"])
    coords_df.index = matrix.index

    meta_keys = sorted({"Population", "Sample Name"} | {k for s in samples.values() for k in s.metadata.keys()})
    color_by = st.selectbox("Color PCoA by metadata field", options=["None"] + meta_keys, index=0)
    color_field = color_by if color_by != "None" else None

    fig = plot_pcoa(coords_df, metadata=metadata, color_by=color_field)
    st.plotly_chart(fig, use_container_width=True)

    st.subheader("🧬 Included Markers")
    with st.expander("Marker filter results", expanded=False):
        included = st.session_state.get("included_markers", [])
        excluded = sorted(set(marker_configs.keys()) - set(included))
        st.write(f"**Included markers ({len(included)}):** {', '.join(included)}")
        st.write(f"**Excluded markers ({len(excluded)}):** {', '.join(excluded)}")

    display_hwe_results(hwe_results)
    fst_section(genotype_df, marker_configs)

    if debug_mode:
        with st.expander("📜 Debug Log"):
            st.text_area("Distance Log", st.session_state.get("distance_debug_log", ""), height=300)

run()