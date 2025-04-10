# interface/ui_file_processor.py

from interface.ui_upload import upload_file_ui
from interface.ui_marker_config import marker_config_ui
from interface.ui_global_config import instantiate_global_config, global_config_ui

import streamlit as st
import pandas as pd
from fla_pipeline.models.genotype import GenotypeResult
from fla_pipeline.config.marker_config import MarkerConfig
from fla_pipeline.pipeline import run_pipeline


def run():
    st.title("Upload FSA Files & Configure Markers")

    upload_col, config_col = st.columns([1, 1])

    with upload_col:
        st.subheader("File Manager")
        upload_file_ui()

    with config_col:
        st.subheader("Marker Configuration")
        marker_config_ui()

    instantiate_global_config()

    if st.button("Run Genotype Calling", type="primary", use_container_width=True, disabled=not (st.session_state.samples and st.session_state.marker_list)):
        with st.spinner("Processing all uploaded samples..."):

            all_calls = []
            cfg = st.session_state.config
            marker_cfgs = []
            for i, m in enumerate(st.session_state.marker_list):
                if isinstance(m, MarkerConfig):
                    marker_cfgs.append(m)
                elif isinstance(m, dict):
                    try:
                        marker_cfgs.append(MarkerConfig(**m))
                    except Exception as e:
                        st.warning(f"Skipping malformed marker config at index {i}: {e}")
                else:
                    msg = f"Unknown marker config type at index {i}: {type(m)}"
                    st.warning(msg)
                    raise RuntimeError(msg)

            for sid, sample in st.session_state.samples.items():
                result = run_pipeline(sample.file_path, cfg, marker_cfgs, sample=sample)
                sample.marker_results = {
                    k: v if isinstance(v, GenotypeResult) else GenotypeResult(**v)
                    for k, v in result["marker_results"].items()
                }

                sample.metadata["max_liz_intensity"] = result.get("max_liz_intensity", 0.0)

                for marker, geno in sample.marker_results.items():
                    all_calls.append({
                        "sample_uid": sample.sample_uid,
                        "sample_id": sid,
                        "Sample Name": sample.metadata.get("Sample Name", sid),
                        "Marker": marker,
                        "Genotype": "/".join([str(i) for i in geno.alleles]),
                        "Confidence": geno.confidence,
                        "QC Flags": "; ".join(geno.qc_flags),
                    })

            st.session_state.genotype_results_df = pd.DataFrame(all_calls)
            st.success(f"Processed {len(st.session_state.samples)} samples.")
            st.rerun()
    global_config_ui()

run()