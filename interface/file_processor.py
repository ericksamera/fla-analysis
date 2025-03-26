# interface/file_processor.py

import streamlit as st
from fla_pipeline.models.sample import Sample
from fla_pipeline.config import MarkerConfig, GlobalConfig
from fla_pipeline.core.fsa_loader import load_fsa
from fla_pipeline.core.peak_detector import detect_peaks
import os
import tempfile
import uuid
import json

from fla_pipeline.pipeline import run_pipeline
from fla_pipeline.models.genotype import GenotypeResult
import pandas as pd


def run():
    st.title("📁 Upload FSA Files & Configure Markers")

    # --- File Upload ---
    uploaded = st.file_uploader(
        "Upload one or more `.fsa` or `.ab1` files", 
        type=["fsa", "ab1"], 
        accept_multiple_files=True
    )

    if uploaded:
        st.session_state.uploaded_files.clear()
        for f in uploaded:
            temp_path = save_temp_file(f)
            sample_id = os.path.splitext(f.name)[0]
            sample = Sample(sample_id=sample_id, file_path=temp_path)
            sample.fsa_data = load_fsa(temp_path)
            sample.peaks = detect_peaks(
                sample.fsa_data["smap"],
                sample.fsa_data["channels"],
                config=GlobalConfig()
            )
            st.session_state.samples[sample_id] = sample
            st.session_state.uploaded_files.append(sample_id)

        st.success(f"Uploaded and parsed {len(uploaded)} file(s).")
        st.rerun()

    # --- Marker Config ---
    st.subheader("🧬 Marker Configuration")
    with st.expander("Load marker config from JSON"):
        marker_file = st.file_uploader("Upload marker config (JSON)", type="json")
        if marker_file:
            marker_data = json.load(marker_file)
            try:
                st.session_state.marker_list = [MarkerConfig(**entry).dict() for entry in marker_data]
                st.success("Loaded marker config.")
                st.session_state.config_changed = True
            except Exception as e:
                st.error(f"Error parsing marker config: {e}")

    # --- Marker Table Preview ---
    if st.session_state.marker_list:
        st.dataframe(st.session_state.marker_list, use_container_width=True)
    else:
        st.info("No marker config loaded yet.")


def save_temp_file(uploaded_file):
    """Save uploaded file to a temp file and return the path."""
    suffix = "." + uploaded_file.name.split(".")[-1]
    temp_file = tempfile.NamedTemporaryFile(delete=False, suffix=suffix)
    temp_file.write(uploaded_file.getvalue())
    temp_file.close()
    return temp_file.name


# --- Process All Button ---
if st.session_state.samples and st.session_state.marker_list:
    st.subheader("⚙️ Process All Samples")

    if st.button("🚀 Run Genotype Calling"):
        with st.spinner("Processing all uploaded samples..."):

            all_calls = []
            global_cfg = GlobalConfig(ploidy=st.session_state.ploidy)

            marker_cfgs = [MarkerConfig(**m) for m in st.session_state.marker_list]

            for sample_id, sample in st.session_state.samples.items():
                result = run_pipeline(sample.file_path, global_cfg, marker_cfgs)
                sample.marker_results = {
                    k: GenotypeResult(**v) for k, v in result["marker_results"].items()
                }

                # Flatten into one-row-per-marker format for DataFrame
                for marker, geno in sample.marker_results.items():
                    all_calls.append({
                        "Sample Name": sample_id,
                        "Marker": marker,
                        "Alleles": geno.alleles,
                        "Confidence": geno.confidence,
                        "QC Flags": "; ".join(geno.qc_flags),
                    })

            st.session_state.genotype_results_df = pd.DataFrame(all_calls)
            st.success(f"Processed {len(st.session_state.samples)} samples.")
            st.rerun()

run()