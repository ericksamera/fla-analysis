import streamlit as st
import json
import os
import numpy as np
import pandas as pd
from _fla_processor import FLAProcessor, MarkerConfig, FLAConfig

# Initialize persistent session state keys if they don't already exist


def mark_config_changed():
    st.session_state.config_changed = True

st.title("FLA Processor Interface")

# Let user upload both .fsa or .json files
uploaded_files = st.file_uploader(
    "Upload .fsa or .json files",
    type=["fsa", "json"],
    accept_multiple_files=True
)

# Add new uploaded files to session state (assigning each a unique id)
if uploaded_files:
    for file in uploaded_files:
        if not any(f["name"] == file.name for f in st.session_state.uploaded_files):
            new_id = st.session_state.uploaded_files_id_counter
            st.session_state.uploaded_files_id_counter += 1
            st.session_state.uploaded_files.append({"id": new_id, "name": file.name, "file": file})

# --- Configurable Parameters for FLA Processor ---
with st.expander("Advanced Options"):
    bin_tolerance = st.number_input("Bin Tolerance", min_value=0, max_value=10, value=2, on_change=mark_config_changed)
    min_peak_height = st.number_input("Minimum Peak Height", min_value=10, max_value=1000, value=1000, on_change=mark_config_changed)
    min_peak_position = st.number_input("Minimum Peak Position", min_value=10, max_value=1000, value=50, on_change=mark_config_changed)
    relative_peak_threshold = st.slider("Relative Peak Threshold", min_value=0.01, max_value=1.0, value=0.1, on_change=mark_config_changed)

    ratio_weight = st.slider("Ratio Weight", min_value=0.1, max_value=5.0, value=1.0, on_change=mark_config_changed)
    stutter_weight = st.slider("Stutter Weight", min_value=0.1, max_value=5.0, value=1.0, on_change=mark_config_changed)
    intensity_weight = st.number_input("Intensity Weight", min_value=0.1, max_value=5.0, value=3.0, on_change=mark_config_changed)
    highest_peak_ratio_weight = st.slider("Highest Peak Ratio Weight", min_value=0.1, max_value=5.0, value=1.0, on_change=mark_config_changed)

    qc_confidence_scale = st.slider("QC Confidence Scale", min_value=0.1, max_value=5.0, value=2.0, on_change=mark_config_changed)
    stutter_ratio_threshold = st.slider("Stutter Ratio Threshold", min_value=0.01, max_value=1.0, value=0.15, on_change=mark_config_changed)
    stutter_tolerance = st.slider("Stutter Tolerance", min_value=0.1, max_value=5.0, value=1.0, on_change=mark_config_changed)

    global_het_ratio_lower_bound = st.slider("Heterozygous Ratio Lower Bound", min_value=0.1, max_value=1.0, value=0.35, on_change=mark_config_changed)
    instrument_saturation_value = st.number_input("Instrument Saturation Value", min_value=1000, max_value=50000, value=30000, on_change=mark_config_changed)
    gaussian_fit_window = st.slider("Gaussian Fit Window", min_value=3, max_value=10, value=5, on_change=mark_config_changed)

    fla_config = FLAConfig(
        bin_tolerance=bin_tolerance,
        min_peak_height=min_peak_height,
        min_peak_position=min_peak_position,
        relative_peak_threshold=relative_peak_threshold,
        ratio_weight=ratio_weight,
        stutter_weight=stutter_weight,
        intensity_weight=intensity_weight,
        highest_peak_ratio_weight=highest_peak_ratio_weight,
        qc_confidence_scale=qc_confidence_scale,
        stutter_ratio_threshold=stutter_ratio_threshold,
        stutter_tolerance=stutter_tolerance,
        global_het_ratio_lower_bound=global_het_ratio_lower_bound,
        instrument_saturation_value=instrument_saturation_value,
        gaussian_fit_window=gaussian_fit_window
    )

# --- Marker Configuration ---
st.header("Marker Configuration")
config_file = st.file_uploader("Upload a JSON config file (optional)", type=["json"])

if config_file:
    config_data = json.load(config_file)
    st.session_state.marker_list = config_data
    st.success("Configuration loaded successfully!")

# Marker editing with persistence
st.write("**Current Marker Configuration:**")
for i, marker in enumerate(st.session_state.marker_list):
    marker_cols = st.columns([1, 1, 2, 1])
    with marker_cols[0]:
        marker["marker"] = st.text_input("Marker Name", value=marker["marker"], key=f"marker_name_{i}")
    with marker_cols[1]:
        marker["channel"] = st.selectbox(
            "Dye Channel", ["6-FAM", "VIC", "NED", "PET", "LIZ"],
            index=["6-FAM", "VIC", "NED", "PET", "LIZ"].index(marker["channel"]),
            key=f"dye_channel_{i}"
        )
    with marker_cols[2]:
        repeat_col, allele_min_col, allele_max_col = st.columns(3)
        with repeat_col:
            marker["repeat_unit"] = st.number_input("Repeat", min_value=1, max_value=10, value=marker["repeat_unit"], key=f"repeat_unit_{i}")
        with allele_min_col:
            allele_min = st.number_input("Min", min_value=50, max_value=500, value=marker["bins"][0], key=f"allele_min_{i}")
        with allele_max_col:
            allele_max = st.number_input("Max", min_value=50, max_value=500, value=marker["bins"][1], key=f"allele_max_{i}")
    with marker_cols[3]:
        st.markdown(" ")
        st.markdown(" ")
        with st.popover("Options"):
            marker["bins"] = [allele_min, allele_max]
            marker["stutter_ratio_threshold"] = st.number_input("Stutter Ratio Threshold", min_value=0.01, max_value=1.0, value=marker["stutter_ratio_threshold"], key=f"stutter_ratio_{i}", on_change=mark_config_changed)
            marker["stutter_tolerance"] = st.number_input("Stutter Tolerance", min_value=0.1, max_value=5.0, value=marker["stutter_tolerance"], key=f"stutter_tolerance_{i}", on_change=mark_config_changed)
            marker["stutter_direction"] = st.segmented_control("Stutter Direction", selection_mode="single", options=["left", "both", "right"], default=marker["stutter_direction"], key=f"stutter_direction_{i}", on_change=mark_config_changed)
            if st.button(f"Remove Marker ({marker['marker']})", key=f"remove_marker_{i}",):
                st.session_state.marker_list.pop(i)
                st.rerun()

if st.button("Add Marker"):
    new_marker = {
        "marker": f"Marker {len(st.session_state.marker_list) + 1}",
        "channel": "6-FAM",
        "bins": [180, 240],
        "repeat_unit": 3,
        "stutter_ratio_threshold": 0.15,
        "stutter_tolerance": 1.0,
        "stutter_direction": "both"
    }
    st.session_state.marker_list.append(new_marker)
    mark_config_changed()
    st.rerun()

if st.button("Save Configuration"):
    config_json = json.dumps(st.session_state.marker_list, indent=4)
    st.download_button(label="Download Configuration", data=config_json, file_name="marker_config.json", mime="application/json")

def convert_to_viewer_format(raw_dict: dict) -> dict:
    color_mapping = {
        "6-FAM": "blue",
        "VIC": "green",
        "NED": "black",
        "PET": "red",
        "LIZ": "orange"
    }
    sample_name = raw_dict.get("file", "Unknown")
    fsa_data = raw_dict.get("fsa_data", {})
    smap = fsa_data.get("smap", [])
    channels_peaks = fsa_data.get("channels", {})
    detected_peaks = raw_dict.get("detected_peaks", {})
    peaks_x = {ch: [p["position"] for p in peaks] for ch, peaks in detected_peaks.items()}
    peaks_y = {ch: [p["intensity"] for p in peaks] for ch, peaks in detected_peaks.items()}
    return {
        "name": sample_name,
        "channels_peaks": channels_peaks,
        "colors": color_mapping,
        "smap": smap,
        "detected_peaks_x": peaks_x,
        "detected_peaks_y": peaks_y,
        "marker_results": raw_dict.get("marker_results", {})
    }

# Process/Load Files using session state's uploaded_files
if st.session_state.uploaded_files and st.button("Process / Load Files"):
    st.session_state.config_changed = False
    results = []
    # If markers exist, use marker analysis
    if len(st.session_state.marker_list) > 0:
        marker_configs = [MarkerConfig(**marker) for marker in st.session_state.marker_list]
        for file_dict in st.session_state.uploaded_files:
            uploaded_file = file_dict["file"]
            file_name = file_dict["name"]
            file_ext = os.path.splitext(file_name)[1].lower()
            if file_ext == ".fsa":
                file_path = os.path.join("/tmp", file_name)
                with open(file_path, "wb") as f:
                    f.write(uploaded_file.getbuffer())
                processor = FLAProcessor(file_path, config=fla_config, marker_config=marker_configs)
                success = processor.process_all()
                if success:
                    for marker, result in processor.marker_results.items():
                        results.append({
                            "Sample Name": file_name,
                            "Marker": marker,
                            "Genotype": result["parsed_allele_calls"],
                            "Genotype Confidence": f"{result['genotype_confidence']:.3f}",
                            "QC Flags": ", ".join(result["QC_flags"]) if result["QC_flags"] else "None"
                        })
                    raw_dict = {
                        "file": processor.fsa_data["name"],
                        "fsa_data": {
                            "smap": processor.fsa_data["smap"],
                            "channels": {ch: list(map(float, arr)) for ch, arr in processor.fsa_data["channels"].items()}
                        },
                        "detected_peaks": {
                            ch: [{"position": p.position, "intensity": p.intensity, "saturated": p.saturated, "note": p.note} for p in peaks]
                            for ch, peaks in processor.detected_peaks.items()
                        },
                        "marker_results": processor.marker_results
                    }
                    processed_for_viewer = convert_to_viewer_format(raw_dict)
                    st.session_state.PROCESSED_FLA[processed_for_viewer["name"]] = processed_for_viewer
                else:
                    st.error(f"Processing failed for {file_name}.")
            elif file_ext == ".json":
                st.write(f"**Loading JSON**: {file_name}")
                try:
                    loaded_data = json.load(uploaded_file)
                    processed_for_viewer = convert_to_viewer_format(loaded_data)
                    st.session_state.PROCESSED_FLA[processed_for_viewer["name"]] = processed_for_viewer
                    st.success(f"JSON loaded into session as '{processed_for_viewer['name']}'.")
                except Exception as exc:
                    st.error(f"Error loading JSON {file_name}: {exc}")
            else:
                st.warning(f"Skipping unknown file type: {file_name}")
        df = pd.DataFrame(results)
        st.dataframe(df)
        st.session_state.genotype_results_df = df
        if df.empty:
            st.info("No results to display.")
    else:
        # No markers: process file, detect peaks, and record each valid peak
        for file_dict in st.session_state.uploaded_files:
            uploaded_file = file_dict["file"]
            file_name = file_dict["name"]
            file_ext = os.path.splitext(file_name)[1].lower()
            if file_ext == ".fsa":
                file_path = os.path.join("/tmp", file_name)
                with open(file_path, "wb") as f:
                    f.write(uploaded_file.getbuffer())
                processor = FLAProcessor(file_path, config=fla_config, marker_config=[])
                if processor.load_fsa():
                    processor.detect_peaks()
                    raw_dict = {
                        "file": processor.fsa_data["name"],
                        "fsa_data": {
                            "smap": processor.fsa_data["smap"],
                            "channels": {ch: list(map(float, arr)) for ch, arr in processor.fsa_data["channels"].items()}
                        },
                        "detected_peaks": {
                            ch: [{"position": p.position, "intensity": p.intensity, "saturated": p.saturated, "note": p.note} for p in peaks]
                            for ch, peaks in processor.detected_peaks.items()
                        },
                        "marker_results": {}
                    }
                    processed_for_viewer = convert_to_viewer_format(raw_dict)
                    st.session_state.PROCESSED_FLA[processed_for_viewer["name"]] = processed_for_viewer
                    for channel, peaks in processor.detected_peaks.items():
                        if channel == "LIZ":
                            continue
                        if peaks:
                            max_intensity = max(p.intensity for p in peaks)
                            threshold = max(0.1 * max_intensity, min_peak_height)
                            for peak in peaks:
                                if peak.position >= 100 and peak.intensity >= threshold:
                                    results.append({
                                        "Sample Name": file_name,
                                        "Marker": channel,
                                        "Position": peak.position,
                                        "Intensity": peak.intensity
                                    })
                else:
                    st.error(f"Loading failed for {file_name}.")
            elif file_ext == ".json":
                st.write(f"**Loading JSON**: {file_name}")
                try:
                    loaded_data = json.load(uploaded_file)
                    processed_for_viewer = convert_to_viewer_format(loaded_data)
                    st.session_state.PROCESSED_FLA[processed_for_viewer["name"]] = processed_for_viewer
                    st.success(f"JSON loaded into session as '{processed_for_viewer['name']}'.")
                except Exception as exc:
                    st.error(f"Error loading JSON {file_name}: {exc}")
            else:
                st.warning(f"Skipping unknown file type: {file_name}")
        df = pd.DataFrame(results)
        st.dataframe(df)
        st.session_state.detected_peaks_df = df

        if df.empty:
            st.info("No results to display.")
    st.rerun()
if st.session_state.config_changed:
    st.warning("Configuration has changed. Please re-analyze your data to apply the new settings.")
