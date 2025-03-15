import streamlit as st
import json
import os
import numpy as np
import pandas as pd
from _fla_processor import FLAProcessor, MarkerConfig, FLAConfig

st.title("FLA Processor Interface")

# Let user upload both .fsa or .json
uploaded_files = st.file_uploader(
    "Upload .fsa or .json files",
    type=["fsa", "json"],
    accept_multiple_files=True
)

# --- Configurable Parameters for FLA Processor ---
with st.expander("Advanced Options"):
    bin_tolerance = st.number_input("Bin Tolerance", min_value=0, max_value=10, value=2)
    min_peak_height = st.number_input("Minimum Peak Height", min_value=10, max_value=1000, value=50)
    min_peak_position = st.number_input("Minimum Peak Position", min_value=10, max_value=1000, value=50)
    relative_peak_threshold = st.slider("Relative Peak Threshold", min_value=0.01, max_value=1.0, value=0.1)

    ratio_weight = st.slider("Ratio Weight", min_value=0.1, max_value=5.0, value=1.0)
    stutter_weight = st.slider("Stutter Weight", min_value=0.1, max_value=5.0, value=1.0)
    intensity_weight = st.slider("Intensity Weight", min_value=0.1, max_value=5.0, value=3.0)
    highest_peak_ratio_weight = st.slider("Highest Peak Ratio Weight", min_value=0.1, max_value=5.0, value=1.0)

    qc_confidence_scale = st.slider("QC Confidence Scale", min_value=0.1, max_value=5.0, value=2.0)
    stutter_ratio_threshold = st.slider("Stutter Ratio Threshold", min_value=0.01, max_value=1.0, value=0.15)
    stutter_tolerance = st.slider("Stutter Tolerance", min_value=0.1, max_value=5.0, value=1.0)

    global_het_ratio_lower_bound = st.slider("Heterozygous Ratio Lower Bound", min_value=0.1, max_value=1.0, value=0.35)

    instrument_saturation_value = st.number_input("Instrument Saturation Value", min_value=1000, max_value=50000, value=30000)
    gaussian_fit_window = st.slider("Gaussian Fit Window", min_value=3, max_value=10, value=5)

    # Create a configuration object dynamically
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

# Initialize or load marker list
if config_file:
    config_data = json.load(config_file)
    st.session_state.marker_list = config_data
    st.success("Configuration loaded successfully!")

if "marker_list" not in st.session_state:
    st.session_state.marker_list = []

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

# Display existing markers
for i, marker in enumerate(st.session_state.marker_list):
    marker_col, dye_col, numerics_col, adv_col = st.columns([1, 1, 2, 1])

    with marker_col:
        marker["marker"] = st.text_input("Marker Name", value=marker["marker"], key=f"marker_name_{i}")

    with dye_col:
        marker["channel"] = st.selectbox(
            "Dye Channel", ["6-FAM", "VIC", "NED", "PET", "LIZ"],
            index=["6-FAM", "VIC", "NED", "PET", "LIZ"].index(marker["channel"]),
            key=f"dye_channel_{i}"
        )

    with numerics_col:
        repeat_col, allele_min_col, allele_max_col = st.columns(3)
        with repeat_col:
            marker["repeat_unit"] = st.number_input(
                "Repeat", min_value=1, max_value=10, value=marker["repeat_unit"], key=f"repeat_unit_{i}"
            )
        with allele_min_col:
            allele_min = st.number_input(
                "Min", min_value=50, max_value=500, value=marker["bins"][0], key=f"allele_min_{i}"
            )
        with allele_max_col:
            allele_max = st.number_input(
                "Max", min_value=50, max_value=500, value=marker["bins"][1], key=f"allele_max_{i}"
            )

    with adv_col:
        st.caption(" ")
        st.text(" ")
        with st.popover("Options"):
            marker["bins"] = [allele_min, allele_max]
            marker["stutter_ratio_threshold"] = st.number_input(
                "Stutter Ratio Threshold",
                min_value=0.01,
                max_value=1.0,
                value=marker["stutter_ratio_threshold"],
                key=f"stutter_ratio_{i}"
            )
            marker["stutter_tolerance"] = st.number_input(
                "Stutter Tolerance",
                min_value=0.1,
                max_value=5.0,
                value=marker["stutter_tolerance"],
                key=f"stutter_tolerance_{i}"
            )
            marker["stutter_direction"] = st.segmented_control(
                "Stutter Direction",
                selection_mode="single",
                options=["left", "both", "right"],
                default=marker["stutter_direction"],
                key=f"stutter_direction_{i}"
            )

            if st.button(f"Remove Marker ({marker['marker']})", key=f"remove_marker_{i}"):
                st.session_state.marker_list.pop(i)
                st.rerun()

# Save Configuration
if st.button("Save Configuration"):
    config_json = json.dumps(st.session_state.marker_list, indent=4)
    st.download_button(
        label="Download Configuration",
        data=config_json,
        file_name="marker_config.json",
        mime="application/json"
    )

# Prepare a dictionary in session to store processed data
if "PROCESSED_FLA" not in st.session_state:
    st.session_state.PROCESSED_FLA = {}

def convert_to_viewer_format(raw_dict: dict) -> dict:
    """
    Convert the raw dictionary from FLAProcessor or user JSON 
    into the structure your viewer page expects:
    
    {
      "name": <string>,
      "channels_peaks": {...},
      "colors": {...},
      "smap": [...],
      "detected_peaks_x": {...},
      "detected_peaks_y": {...},
      "marker_results": {...}    # optional
    }
    """
    color_mapping = {
        "6-FAM": "blue",
        "VIC":   "green",
        "NED":   "black",
        "PET":   "red",
        "LIZ":   "orange"
    }

    # The main identifier: if "file" doesn't exist, fallback
    sample_name = raw_dict.get("file", "Unknown")

    # Pull the SMAP and channels
    fsa_data = raw_dict.get("fsa_data", {})
    smap = fsa_data.get("smap", [])
    channels_peaks = fsa_data.get("channels", {})
    
    # Build x/y from 'detected_peaks'
    detected_peaks = raw_dict.get("detected_peaks", {})
    peaks_x = {}
    peaks_y = {}
    for ch, peaks in detected_peaks.items():
        px = [p["position"] for p in peaks]
        py = [p["intensity"] for p in peaks]
        peaks_x[ch] = px
        peaks_y[ch] = py

    viewer_dict = {
        "name": sample_name,
        "channels_peaks": channels_peaks,
        "colors": color_mapping,
        "smap": smap,
        "detected_peaks_x": peaks_x,
        "detected_peaks_y": peaks_y,
        # If you'd like to keep marker results for advanced analysis in the viewer:
        "marker_results": raw_dict.get("marker_results", {})
    }
    return viewer_dict

# Button to process or load files
if uploaded_files and st.button("Process / Load Files"):
    results = []
    # Build MarkerConfigs from the session marker_list
    marker_configs = [MarkerConfig(**marker) for marker in st.session_state.marker_list]

    for uploaded_file in uploaded_files:
        file_name = uploaded_file.name
        file_ext = os.path.splitext(file_name)[1].lower()

        # If it's an FSA, run the FLAProcessor
        if file_ext == ".fsa":
            st.write(f"**Processing FSA**: {file_name}")
            file_path = os.path.join("/tmp", file_name)
            with open(file_path, "wb") as f:
                f.write(uploaded_file.getbuffer())

            processor = FLAProcessor(file_path, config=fla_config, marker_config=marker_configs)
            success = processor.process_all()

            if success:
                # Collect genotype results for display
                for marker, result in processor.marker_results.items():
                    results.append({
                        "Sample Name": file_name,
                        "Marker": marker,
                        "Genotype": result["parsed_allele_calls"],
                        "Genotype Confidence": f"{result['genotype_confidence']:.3f}",
                        "QC Flags": ", ".join(result["QC_flags"]) if result["QC_flags"] else "None"
                    })

                # Build raw JSON-like dict from the processor
                raw_dict = {
                    "file": processor.fsa_data["name"],
                    "fsa_data": {
                        "smap": processor.fsa_data["smap"],
                        "channels": {
                            ch: list(map(float, arr))
                            for ch, arr in processor.fsa_data["channels"].items()
                        }
                    },
                    "detected_peaks": {
                        ch: [
                            {
                                "position": p.position,
                                "intensity": p.intensity,
                                "saturated": p.saturated,
                                "note": p.note
                            }
                            for p in peaks
                        ]
                        for ch, peaks in processor.detected_peaks.items()
                    },
                    "marker_results": processor.marker_results
                }

                # Convert to viewer format
                processed_for_viewer = convert_to_viewer_format(raw_dict)
                st.session_state.PROCESSED_FLA[processed_for_viewer["name"]] = processed_for_viewer

            else:
                st.error(f"Processing failed for {file_name}. Check the file format and try again.")

        # If it's JSON, parse, then convert to viewer format
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

    # Display genotype results for any newly processed FSA
    if results:
        st.subheader("Genotype Results Table")
        df = pd.DataFrame(results)
        st.dataframe(df)

st.info(
    "All processed or loaded traces are now in `st.session_state.PROCESSED_FLA`.\n\n"
    "Each entry is a 'viewer dict' with keys: 'name', 'channels_peaks', 'colors', 'smap',\n"
    "'detected_peaks_x', 'detected_peaks_y', and optional 'marker_results'.\n\n"
    "Switch to your other Viewer page to see detailed electropherograms."
)