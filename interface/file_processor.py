# interface/file_processor.py

import streamlit as st
import os
import tempfile
import json
import pandas as pd

from fla_pipeline.models.sample import Sample
from fla_pipeline.models.genotype import GenotypeResult
from fla_pipeline.config import MarkerConfig, GlobalConfig
from fla_pipeline.core.fsa_loader import load_fsa
from fla_pipeline.core.peak_detector import detect_peaks
from fla_pipeline.pipeline import run_pipeline

@st.dialog("Marker Overrides", width="large")
def edit_marker_overrides_dialog(index: int):
    marker = st.session_state.marker_list[index]
    overrides = marker.setdefault("overrides", {})

@st.dialog("Upload Files", width="large")
def show_upload_dialog():
    uploaded = st.file_uploader(
        "Upload .fsa or .ab1 files",
        type=["fsa", "ab1"],
        accept_multiple_files=True,
    )

    if uploaded:
        updated = False
        for f in uploaded:
            temp_path = save_temp_file(f)
            if not any(existing["name"] == f.name for existing in st.session_state.uploaded_files):
                updated = True
                st.session_state.uploaded_files.append({
                    "name": f.name,
                    "sample_name": os.path.splitext(f.name)[0],
                    "path": temp_path,
                })

        if updated:
            st.success("Files uploaded successfully.")
            st.rerun()


def run():
    st.title("📁 Upload FSA Files & Configure Markers")

    # --- Upload Section ---
    with st.container():
        st.subheader("📤 Upload Files")

        upload_col, table_col = st.columns([1, 2])

        with upload_col:
            if st.button("📤 Upload Files", use_container_width=True):
                show_upload_dialog()


        with table_col:
            show_uploaded_files_table()

    # --- Marker Config ---
    with st.container():
        st.subheader("🧬 Marker Configuration")

        with st.expander("Load marker config from JSON"):
            marker_file = st.file_uploader("Upload marker config (JSON)", type="json")
            if marker_file:
                marker_data = json.load(marker_file)
                try:
                    st.session_state.marker_list = [MarkerConfig(**entry).dict() for entry in marker_data]
                    st.success("Loaded marker config.")
                    st.session_state.config_changed = True
                    st.rerun()
                except Exception as e:
                    st.error(f"Error parsing marker config: {e}")

        st.markdown("### Edit Markers")
        for i, marker in enumerate(st.session_state.marker_list):
            cols = st.columns([2, 2, 1, 1, 1, 2], vertical_alignment="bottom")
            marker["marker"] = cols[0].text_input("Marker", value=marker["marker"], key=f"marker_{i}")
            marker["channel"] = cols[1].selectbox("Channel", ["6-FAM", "VIC", "NED", "PET"],
                                                   index=["6-FAM", "VIC", "NED", "PET"].index(marker["channel"]),
                                                   key=f"channel_{i}")
            marker["repeat_unit"] = cols[2].number_input("Repeat", 1, 10, value=marker.get("repeat_unit", 1), key=f"repeat_{i}")
            marker["bins"] = (
                cols[3].number_input("Min", 30, 500, value=marker["bins"][0], key=f"bmin_{i}"),
                cols[4].number_input("Max", 30, 500, value=marker["bins"][1], key=f"bmax_{i}")
            )
            with cols[5]:
                if st.button("Options", key=f"edit_overrides_{i}", use_container_width=True):
                    st.session_state.override_index = i

        st.markdown("### Actions")
        action_cols = st.columns(3)
        if action_cols[0].button("➕ Add Marker"):
            st.session_state.marker_list.append({
                "marker": f"Marker{len(st.session_state.marker_list)+1}",
                "channel": "6-FAM",
                "repeat_unit": 3,
                "bins": [180, 240],
                "overrides": {}
            })
            st.rerun()

        if action_cols[1].button("💾 Export Marker Config"):
            config_json = json.dumps(st.session_state.marker_list, indent=4)
            st.download_button(
                label="Download JSON",
                data=config_json,
                file_name="marker_config.json",
                mime="application/json"
            )

    # --- Processing Section ---
    if st.session_state.samples and st.session_state.marker_list:
        with st.container():
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
                        sample.metadata["max_liz_intensity"] = result.get("max_liz_intensity", 0.0)

                        for marker, geno in sample.marker_results.items():
                            all_calls.append({
                                "Sample Name": sample_id,
                                "Marker": marker,
                                "Genotype": "/".join([str(i) for i in geno.alleles]),
                                "Confidence": geno.confidence,
                                "QC Flags": "; ".join(geno.qc_flags),
                            })

                    st.session_state.genotype_results_df = pd.DataFrame(all_calls)
                    st.success(f"Processed {len(st.session_state.samples)} samples.")
                    st.rerun()


def save_temp_file(uploaded_file):
    suffix = "." + uploaded_file.name.split(".")[-1]
    temp_file = tempfile.NamedTemporaryFile(delete=False, suffix=suffix)
    temp_file.write(uploaded_file.getvalue())
    temp_file.close()
    return temp_file.name


def show_uploaded_files_table():
    st.markdown("### 📑 Uploaded Files")

    if not st.session_state.uploaded_files:
        st.info("No files uploaded yet.")
        return

    df = pd.DataFrame([
        {
            "Filename": f["name"],
            "Sample Name": f["sample_name"],
            "Delete?": False
        } for f in st.session_state.uploaded_files
    ])

    edited = st.data_editor(df, use_container_width=True, num_rows="fixed")

    # Apply edits
    for i, row in edited.iterrows():
        st.session_state.uploaded_files[i]["sample_name"] = row["Sample Name"]

    # Handle deletions
    to_delete = edited[edited["Delete?"]].index.tolist()
    if to_delete:
        st.session_state.uploaded_files = [
            f for i, f in enumerate(st.session_state.uploaded_files) if i not in to_delete
        ]
        for i in to_delete:
            sid = df.iloc[i]["Sample Name"]
            st.session_state.samples.pop(sid, None)
        st.rerun()

    # Parse uploaded files if not yet processed
    for f in st.session_state.uploaded_files:
        sid = f["sample_name"]
        if sid not in st.session_state.samples:
            try:
                sample = Sample(sample_id=sid, file_path=f["path"])
                sample.fsa_data = load_fsa(f["path"])
                peaks, max_liz, suppressed = detect_peaks(
                    sample.fsa_data["smap"],
                    sample.fsa_data["channels"],
                    config=GlobalConfig()
                )
                sample.peaks = peaks
                sample.suppressed_peaks = suppressed
                sample.metadata["max_liz_intensity"] = max_liz
                st.session_state.samples[sid] = sample
                st.rerun()  # Ensure fresh uploads are available immediately
            except Exception as e:
                st.error(f"Error parsing {f['name']}: {e}")



run()