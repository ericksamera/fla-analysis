# interface/file_processor.py

import streamlit as st
import json
import os
import tempfile
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

    st.markdown(f"### Editing Overrides for `{marker['marker']}`")

    overrides["relative_peak_threshold"] = st.slider(
        "Relative Peak Threshold",
        min_value=0.01, max_value=1.0,
        value=overrides.get("relative_peak_threshold", st.session_state.config.relative_peak_threshold),
        key=f"override_rel_thresh_dialog_{index}"
    )

    overrides["stutter_ratio_threshold"] = st.slider(
        "Stutter Ratio Threshold",
        min_value=0.01, max_value=1.0,
        value=overrides.get("stutter_ratio_threshold", 0.15),
        key=f"override_stutter_ratio_dialog_{index}"
    )

    overrides["stutter_tolerance"] = st.slider(
        "Stutter Tolerance",
        min_value=0.1, max_value=5.0,
        value=overrides.get("stutter_tolerance", 1.0),
        key=f"override_stutter_tolerance_dialog_{index}"
    )

    direction = overrides.get("stutter_direction", "both")
    overrides["stutter_direction"] = st.selectbox(
        "Stutter Direction",
        options=["left", "both", "right"],
        index=["left", "both", "right"].index(direction),
        key=f"override_stutter_dir_dialog_{index}"
    )

    if st.button("✅ Save and Close"):
        st.session_state.marker_list[index]["overrides"] = overrides
        st.session_state.config_changed = True
        st.rerun()



def save_temp_file(uploaded_file):
    suffix = "." + uploaded_file.name.split(".")[-1]
    temp_file = tempfile.NamedTemporaryFile(delete=False, suffix=suffix)
    temp_file.write(uploaded_file.getvalue())
    temp_file.close()
    return temp_file.name


def ensure_config_initialized():
    if "config" not in st.session_state or not isinstance(st.session_state.config, GlobalConfig):
        st.session_state.config = GlobalConfig()


def global_config_ui():
    cfg: GlobalConfig = st.session_state.config
    with st.expander("⚙️ Advanced Global Settings", expanded=False):
        cols = st.columns(3)
        cfg.ploidy = cols[0].number_input("Ploidy", min_value=1, max_value=8, value=cfg.ploidy, key="ploidy")
        cfg.min_peak_height = cols[1].number_input("Min Peak Height", value=cfg.min_peak_height, key="min_peak_height")
        cfg.min_peak_position = cols[2].number_input("Min Peak Position", value=cfg.min_peak_position, key="min_peak_pos")

        cfg.relative_peak_threshold = st.slider("Relative Peak Threshold", 0.01, 1.0, value=cfg.relative_peak_threshold)
        cfg.intensity_weight = st.slider("Intensity Weight", 0.1, 5.0, value=cfg.intensity_weight)
        cfg.ratio_weight = st.slider("Ratio Weight", 0.1, 5.0, value=cfg.ratio_weight)
        cfg.stutter_weight = st.slider("Stutter Weight", 0.1, 5.0, value=cfg.stutter_weight)
        cfg.qc_confidence_scale = st.slider("QC Confidence Scale", 0.1, 2.0, value=cfg.qc_confidence_scale)
        cfg.diploid_strategy = st.selectbox("Diploid Strategy", ["probabilistic", "strict_stutter", "lenient_ratio"], index=["probabilistic", "strict_stutter", "lenient_ratio"].index(cfg.diploid_strategy))


def marker_config_ui():
    st.subheader("🧬 Marker Configuration")

    with st.expander("Load marker config from JSON"):
        marker_file = st.file_uploader("Upload marker config (JSON)", type="json", key="marker_json")
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
        marker["channel"] = cols[1].selectbox("Channel", ["6-FAM", "VIC", "NED", "PET"], index=["6-FAM", "VIC", "NED", "PET",].index(marker["channel"]), key=f"channel_{i}")
        marker["repeat_unit"] = cols[2].number_input("Repeat", 1, 10, value=marker.get("repeat_unit", 1), key=f"repeat_{i}")
        marker["bins"] = (
            cols[3].number_input("Min", 30, 500, value=marker["bins"][0], key=f"bmin_{i}"),
            cols[4].number_input("Max", 30, 500, value=marker["bins"][1], key=f"bmax_{i}")
        )
        with cols[5]:
            if st.button("Options", key=f"edit_overrides_{i}", use_container_width=True):
                edit_marker_overrides_dialog(i)

        # Optional overrides UI
        # if f"show_overrides_{i}" not in st.session_state:
        #     st.session_state[f"show_overrides_{i}"] = False
        # with st.expander("⚙️ Marker Overrides", expanded=False):
        #     marker.setdefault("overrides", {})
        #     marker["overrides"]["relative_peak_threshold"] = st.slider(
        #         f"[{marker['marker']}] Rel. Threshold", 0.01, 1.0,
        #         value=marker["overrides"].get("relative_peak_threshold", st.session_state.config.relative_peak_threshold),
        #         key=f"override_rel_thresh_{i}"
        #     )

    st.markdown("### Actions")
    cols = st.columns(3)
    if cols[0].button("➕ Add Marker"):
        st.session_state.marker_list.append({
            "marker": f"Marker{len(st.session_state.marker_list)+1}",
            "channel": "6-FAM",
            "repeat_unit": 3,
            "bins": [180, 240],
            "overrides": {}
        })
        st.rerun()

    if cols[1].button("💾 Export Marker Config"):
        config_json = json.dumps(st.session_state.marker_list, indent=4)
        st.download_button(
            label="Download JSON",
            data=config_json,
            file_name="marker_config.json",
            mime="application/json"
        )


def run():
    ensure_config_initialized()
    st.title("📁 Upload FSA Files & Configure Markers")

    # --- Layout ---
    upload_col, config_col = st.columns([1, 1])

    # --- File Upload ---
    with upload_col:
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
                try:
                    sample.fsa_data = load_fsa(temp_path)
                    peak_dict, max_liz, suppressed_peaks = detect_peaks(
                        sample.fsa_data["smap"],
                        sample.fsa_data["channels"],
                        config=st.session_state.config
                    )
                    sample.suppressed_peaks = suppressed_peaks
                    sample.peaks = peak_dict
                    sample.metadata["max_liz_intensity"] = max_liz
                    st.session_state.samples[sample_id] = sample
                    st.session_state.uploaded_files.append(sample_id)
                except Exception as e:
                    st.error(f"Error processing {f.name}: {e}")

            st.success(f"Uploaded and parsed {len(uploaded)} file(s).")
            st.rerun()

    # --- Global Config Editor ---
    with config_col:
        marker_config_ui()

    # --- Marker Config ---
    

    # --- Process All Button ---
    global_config_ui()  
    if st.session_state.samples and st.session_state.marker_list:
        st.subheader("⚙️ Process All Samples")
        if st.button("🚀 Run Genotype Calling", use_container_width=True):
            with st.spinner("Processing all uploaded samples..."):
                all_calls = []
                global_cfg = st.session_state.config
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

run()
