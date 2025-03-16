import streamlit as st
import json
import os
import numpy as np
import pandas as pd
from dataclasses import dataclass
from typing import Optional, Tuple, List, Dict

###############################################################################
# 1) Temporary local data classes for storing user config (FLAConfig, MarkerConfig)
###############################################################################
# We define these locally so we can gather user input (including ploidy) BEFORE 
# deciding which module (_fla_processor or _fla_process_with_ploidy) to import.
###############################################################################


@dataclass
class FLAConfigManager:
    """
    Manages advanced FLA processing parameters and ensures they persist across sessions.
    """
    bin_tolerance: int = 2
    min_peak_height: int = 1000
    min_peak_position: int = 50
    relative_peak_threshold: float = 0.1
    ratio_weight: float = 1.0
    stutter_weight: float = 1.0
    intensity_weight: float = 3.0
    highest_peak_ratio_weight: float = 1.0
    qc_confidence_scale: float = 2.0
    stutter_ratio_threshold: float = 0.15
    stutter_tolerance: float = 1.0
    global_het_ratio_lower_bound: float = 0.35
    gaussian_fit_window: int = 5
    ploidy: int = 2  # Default is diploid

    @staticmethod
    def initialize():
        """Ensures an instance of FLAConfigManager exists in session state."""
        if "fla_config" not in st.session_state:
            st.session_state.fla_config = FLAConfigManager()

    @staticmethod
    def update_from_ui():
        """Updates the FLAConfigManager instance with user-defined values."""
        config = st.session_state.fla_config
        with st.expander("Advanced Options"):
            config.ploidy = st.number_input("Ploidy (1-8)", min_value=1, max_value=8, value=config.ploidy)
            config.bin_tolerance = st.number_input("Bin Tolerance", min_value=0, max_value=10, value=config.bin_tolerance)
            config.min_peak_height = st.number_input("Minimum Peak Height", min_value=10, max_value=1000, value=config.min_peak_height)
            config.min_peak_position = st.number_input("Minimum Peak Position", min_value=10, max_value=1000, value=config.min_peak_position)
            config.relative_peak_threshold = st.slider("Relative Peak Threshold", min_value=0.01, max_value=1.0, value=config.relative_peak_threshold)
            config.ratio_weight = st.slider("Ratio Weight", min_value=0.1, max_value=5.0, value=config.ratio_weight)
            config.stutter_weight = st.slider("Stutter Weight", min_value=0.1, max_value=5.0, value=config.stutter_weight)
            config.intensity_weight = st.number_input("Intensity Weight", min_value=0.1, max_value=5.0, value=config.intensity_weight)
            config.highest_peak_ratio_weight = st.slider("Highest Peak Ratio Weight", min_value=0.1, max_value=5.0, value=config.highest_peak_ratio_weight)
            config.qc_confidence_scale = st.slider("QC Confidence Scale", min_value=0.1, max_value=5.0, value=config.qc_confidence_scale)
            config.stutter_ratio_threshold = st.slider("Stutter Ratio Threshold", min_value=0.01, max_value=1.0, value=config.stutter_ratio_threshold)
            config.stutter_tolerance = st.slider("Stutter Tolerance", min_value=0.1, max_value=5.0, value=config.stutter_tolerance)
            config.global_het_ratio_lower_bound = st.slider("Heterozygous Ratio Lower Bound", min_value=0.1, max_value=1.0, value=config.global_het_ratio_lower_bound)
            config.gaussian_fit_window = st.slider("Gaussian Fit Window", min_value=3, max_value=10, value=config.gaussian_fit_window)

    @staticmethod
    def get():
        """Returns the current configuration object from session state."""
        return st.session_state.fla_config

@dataclass
class FLAConfigLocal:
    bin_tolerance: int
    min_peak_height: int
    min_peak_position: int
    relative_peak_threshold: float
    ratio_weight: float
    stutter_weight: float
    intensity_weight: float
    highest_peak_ratio_weight: float
    qc_confidence_scale: float
    stutter_ratio_threshold: float
    stutter_tolerance: float
    global_het_ratio_lower_bound: float
    gaussian_fit_window: int
    ploidy: int

@dataclass
class MarkerConfigLocal:
    marker: str
    channel: str
    repeat_unit: int
    bins: Tuple[int, int]
    stutter_ratio_threshold: float
    stutter_tolerance: float
    stutter_direction: str = "both"

###############################################################################
# 2) Utility / Helper functions that do not depend on the chosen module
###############################################################################
def mark_config_changed():
    """
    Sets a session state flag indicating that configuration has changed
    and that re-analysis may be needed.
    """
    st.session_state.config_changed = True

def add_toast_message(msg: str, icon: str = None):
    """
    Appends a toast message with an optional icon to session_state,
    to be displayed after the next rerun.
    """
    st.session_state.toast_messages.append((msg, icon))

def convert_to_viewer_format(raw_dict: dict) -> dict:
    """
    Convert the dictionary of raw data (produced by either processor or a loaded JSON)
    into the format required by the viewer application.
    """
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

###############################################################################
# 3) Dialog / Upload Handling
###############################################################################
@st.dialog("Upload Files", width="large")
def upload_files_dialog():
    """
    A Streamlit dialog used for uploading .fsa or .json files.
    Files are stored in st.session_state.uploaded_files.
    """
    uploaded_files = st.file_uploader(
        "Upload .fsa or .json files",
        type=["fsa", "json"],
        accept_multiple_files=True
    )
    if uploaded_files:
        for file in uploaded_files:
            if not any(f["name"] == file.name for f in st.session_state.uploaded_files):
                new_id = st.session_state.uploaded_files_id_counter
                st.session_state.uploaded_files_id_counter += 1
                st.session_state.uploaded_files.append({
                    "name": file.name,
                    "sample_name": file.name,
                    "file": file
                })
        st.success("Files uploaded successfully!")
        st.rerun()

def display_uploaded_files_table():
    """
    Displays a table with the uploaded files, allowing rename or delete.
    """
    if st.session_state.uploaded_files:
        st.write("### Uploaded Files")

        df_uploaded_files = pd.DataFrame([
            {"Filename": f["name"], "Sample Name": f["sample_name"], "delete": False}
            for f in st.session_state.uploaded_files
        ])
        
        edited_df = st.data_editor(
            df_uploaded_files,
            use_container_width=True,
            num_rows="fixed",
            disabled=["Filename"],
            column_config={"delete": "Delete?", "Sample Name": "Sample Name"}
        )

        # Update session state with edited sample names
        for file_dict in st.session_state.uploaded_files:
            if file_dict["name"] in edited_df["Filename"].values:
                file_dict["sample_name"] = edited_df.loc[
                    edited_df["Filename"] == file_dict["name"], "Sample Name"
                ].values[0]

        # Process deletions
        if "delete" in edited_df.columns:
            to_delete = edited_df[edited_df["delete"] == True]["Filename"].tolist()
            if to_delete:
                mark_config_changed()
                st.session_state.uploaded_files = [
                    f for f in st.session_state.uploaded_files if f["name"] not in to_delete
                ]
                st.rerun()

###############################################################################
# 4) Advanced Options (including Ploidy) => returns FLAConfigLocal
###############################################################################
def display_advanced_options() -> FLAConfigLocal:
    """
    Renders the Advanced Options expander for configuring the analysis parameters,
    including the ploidy selection in the UI. Returns a local config object that
    we can use to decide which module to import.
    """
    with st.expander("Advanced Options"):
        ploidy_chosen = st.number_input(
            "Ploidy (1-8)",
            min_value=1, max_value=8,
            value=2,
            help="Choose '2' for diploid, or a higher number for polyploid analysis.",
            on_change=mark_config_changed
        )

        bin_tolerance = st.number_input(
            "Bin Tolerance", min_value=0, max_value=10, value=2, on_change=mark_config_changed
        )
        min_peak_height = st.number_input(
            "Minimum Peak Height", min_value=10, max_value=1000, value=1000, on_change=mark_config_changed
        )
        min_peak_position = st.number_input(
            "Minimum Peak Position", min_value=10, max_value=1000, value=50, on_change=mark_config_changed
        )
        relative_peak_threshold = st.slider(
            "Relative Peak Threshold", min_value=0.01, max_value=1.0, value=0.1,
            on_change=mark_config_changed
        )

        ratio_weight = st.slider(
            "Ratio Weight", min_value=0.1, max_value=5.0, value=1.0, on_change=mark_config_changed
        )
        stutter_weight = st.slider(
            "Stutter Weight", min_value=0.1, max_value=5.0, value=1.0, on_change=mark_config_changed
        )
        intensity_weight = st.number_input(
            "Intensity Weight", min_value=0.1, max_value=5.0, value=3.0, on_change=mark_config_changed
        )
        highest_peak_ratio_weight = st.slider(
            "Highest Peak Ratio Weight", min_value=0.1, max_value=5.0, value=1.0,
            on_change=mark_config_changed
        )

        qc_confidence_scale = st.slider(
            "QC Confidence Scale", min_value=0.1, max_value=5.0, value=2.0,
            on_change=mark_config_changed
        )
        stutter_ratio_threshold = st.slider(
            "Stutter Ratio Threshold", min_value=0.01, max_value=1.0, value=0.15,
            on_change=mark_config_changed
        )
        stutter_tolerance = st.slider(
            "Stutter Tolerance", min_value=0.1, max_value=5.0, value=1.0,
            on_change=mark_config_changed
        )

        global_het_ratio_lower_bound = st.slider(
            "Heterozygous Ratio Lower Bound", min_value=0.1, max_value=1.0, value=0.35,
            on_change=mark_config_changed
        )
        gaussian_fit_window = st.slider(
            "Gaussian Fit Window", min_value=3, max_value=10, value=5,
            on_change=mark_config_changed
        )

    return FLAConfigLocal(
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
        gaussian_fit_window=gaussian_fit_window,
        ploidy=ploidy_chosen
    )

###############################################################################
# 5) Marker Configuration
###############################################################################
def display_marker_configuration():
    """
    Renders the UI for marker configuration: upload config, add/edit/remove markers,
    and optionally save the configuration.
    """
    st.header("Marker Configuration")
    config_file = st.file_uploader("Upload a JSON config file (optional)", type=["json"])

    if config_file:
        config_data = json.load(config_file)
        st.session_state.marker_list = config_data
        st.success("Configuration loaded successfully!")

    st.write("**Current Marker Configuration:**")
    for i, marker in enumerate(st.session_state.marker_list):
        marker_cols = st.columns([1, 1, 2, 1])
        with marker_cols[0]:
            marker["marker"] = st.text_input(
                "Marker Name", value=marker["marker"], key=f"marker_name_{i}"
            )
        with marker_cols[1]:
            marker["channel"] = st.selectbox(
                "Dye Channel",
                ["6-FAM", "VIC", "NED", "PET", "LIZ"],
                index=["6-FAM", "VIC", "NED", "PET", "LIZ"].index(marker["channel"]),
                key=f"dye_channel_{i}"
            )
        with marker_cols[2]:
            repeat_col, allele_min_col, allele_max_col = st.columns(3)
            with repeat_col:
                marker["repeat_unit"] = st.number_input(
                    "Repeat", min_value=1, max_value=10,
                    value=marker["repeat_unit"], key=f"repeat_unit_{i}"
                )
            with allele_min_col:
                allele_min = st.number_input(
                    "Min", min_value=50, max_value=500,
                    value=marker["bins"][0], key=f"allele_min_{i}"
                )
            with allele_max_col:
                allele_max = st.number_input(
                    "Max", min_value=50, max_value=500,
                    value=marker["bins"][1], key=f"allele_max_{i}"
                )
        with marker_cols[3]:
            st.markdown(" ")
            st.markdown(" ")
            with st.popover("Options"):
                marker["bins"] = [allele_min, allele_max]
                marker["stutter_ratio_threshold"] = st.number_input(
                    "Stutter Ratio Threshold",
                    min_value=0.01, max_value=1.0,
                    value=marker["stutter_ratio_threshold"],
                    key=f"stutter_ratio_{i}",
                    on_change=mark_config_changed
                )
                marker["stutter_tolerance"] = st.number_input(
                    "Stutter Tolerance",
                    min_value=0.1, max_value=5.0,
                    value=marker["stutter_tolerance"],
                    key=f"stutter_tolerance_{i}",
                    on_change=mark_config_changed
                )
                marker["stutter_direction"] = st.segmented_control(
                    "Stutter Direction",
                    selection_mode="single",
                    options=["left", "both", "right"],
                    default=marker["stutter_direction"],
                    key=f"stutter_direction_{i}",
                    on_change=mark_config_changed
                )

                if st.button(f"Remove Marker ({marker['marker']})", key=f"remove_marker_{i}"):
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
        st.download_button(
            label="Download Configuration",
            data=config_json,
            file_name="marker_config.json",
            mime="application/json"
        )

###############################################################################
# 6) Dynamic Import based on ploidy (inside a helper function)
###############################################################################
def import_processor(ploidy: int):
    """
    Dynamically import either _fla_processor or _fla_process_with_ploidy
    based on the user's ploidy setting.
    Returns a triple: (FLAProcessor, MarkerConfig, FLAConfig).
    """
    st.session_state.ploidy = ploidy
    if ploidy == 2:
        from _fla_processor import FLAProcessor, MarkerConfig, FLAConfig
        return FLAProcessor, MarkerConfig, FLAConfig
    else:
        from _fla_process_with_ploidy import FLAProcessor, MarkerConfig, FLAConfig
        return FLAProcessor, MarkerConfig, FLAConfig

###############################################################################
# 7) File Processing Helpers
###############################################################################
def process_single_fsa(file_path: str,
                       fla_config,
                       marker_configs,
                       has_markers: bool,
                       FLAProcessorClass) -> dict:
    """
    Load & process a single FSA file. Uses the FLAProcessor class dynamically 
    provided by import_processor().
    """
    processor = FLAProcessorClass(file_path, config=fla_config, marker_config=marker_configs)

    # Attempt to load
    if not processor.load_fsa():
        return {}

    # Check if SMAP is broken
    if not processor.fsa_data.get("smap"):
        return {}

    if has_markers:
        success = processor.process_all()
        if not success:
            return {}
    else:
        processor.detect_peaks()

    # Build the raw dict
    raw_dict = {
        "file": processor.fsa_data["name"],
        "fsa_data": {
            "smap": processor.fsa_data["smap"],
            "channels": {
                ch: list(map(float, arr)) for ch, arr in processor.fsa_data["channels"].items()
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
    return raw_dict

def process_single_json(file_obj) -> dict:
    """Load a single JSON from the given file object."""
    try:
        return json.load(file_obj)
    except Exception as exc:
        st.error(f"Error loading JSON: {exc}")
        return {}

def process_or_load_files(fla_config):
    """
    Orchestrates the entire file processing. 
    1) Dynamically import the correct FLAProcessor based on ploidy in fla_config_local
    2) Build the final real FLAConfig for that module
    3) Process each uploaded file, storing results in st.session_state
    """
    st.session_state.config_changed = False

    # 1) Dynamically import
    FLAProcessorClass, MarkerConfigReal, FLAConfigReal = import_processor(fla_config.ploidy)

    # 2) Construct the FLAConfig instance for the chosen module
    real_fla_config = FLAConfigReal(**vars(fla_config))

    # 3) Build marker configs for the chosen module
    has_markers = len(st.session_state.marker_list) > 0
    if has_markers:
        marker_configs = []
        for m in st.session_state.marker_list:
            mc = MarkerConfigReal(
                marker=m["marker"],
                channel=m["channel"],
                repeat_unit=m["repeat_unit"],
                bins=tuple(m["bins"]),
                stutter_ratio_threshold=m["stutter_ratio_threshold"],
                stutter_tolerance=m["stutter_tolerance"],
                stutter_direction=m.get("stutter_direction", "both")
            )
            marker_configs.append(mc)
    else:
        marker_configs = []

    # 4) Process each uploaded file
    results = []
    files_to_remove = []

    for file_dict in st.session_state.uploaded_files:
        uploaded_file = file_dict["file"]
        file_name = file_dict["name"]
        file_ext = os.path.splitext(file_name)[1].lower()
        file_path = os.path.join("/tmp", file_name)

        # Write local copy
        with open(file_path, "wb") as f:
            f.write(uploaded_file.getbuffer())

        if file_ext == ".fsa":
            raw_dict = process_single_fsa(file_path, real_fla_config, marker_configs,
                                          has_markers, FLAProcessorClass)
            if not raw_dict:
                add_toast_message(f"{file_name} failed to load!", ':material/error:')
                files_to_remove.append(file_name)
                continue

            # Gather genotype data if markers are used
            if has_markers and "marker_results" in raw_dict:
                for marker_name, result in raw_dict["marker_results"].items():
                    results.append({
                        "Sample Name": file_name,
                        "Marker": marker_name,
                        "Genotype": result["parsed_allele_calls"],
                        "Genotype Confidence": f"{result['genotype_confidence']:.3f}",
                        "QC Flags": ", ".join(result["QC_flags"]) if result["QC_flags"] else "None"
                    })

            # Convert for viewer
            processed_for_viewer = convert_to_viewer_format(raw_dict)
            st.session_state.PROCESSED_FLA[processed_for_viewer["name"]] = processed_for_viewer

        elif file_ext == ".json":
            st.write(f"**Loading JSON**: {file_name}")
            loaded_data = process_single_json(uploaded_file)
            if not loaded_data:
                add_toast_message(f"Removing '{file_name}' (invalid JSON).")
                files_to_remove.append(file_name)
                continue

            processed_for_viewer = convert_to_viewer_format(loaded_data)
            st.session_state.PROCESSED_FLA[processed_for_viewer["name"]] = processed_for_viewer
            st.success(f"JSON loaded into session as '{processed_for_viewer['name']}'.")

            if has_markers and "marker_results" in loaded_data:
                for marker_name, result in loaded_data["marker_results"].items():
                    results.append({
                        "Sample Name": file_name,
                        "Marker": marker_name,
                        "Genotype": result["parsed_allele_calls"],
                        "Genotype Confidence": f"{result['genotype_confidence']:.3f}",
                        "QC Flags": ", ".join(result["QC_flags"]) if result["QC_flags"] else "None"
                    })
        else:
            st.warning(f"Skipping unknown file type: {file_name}")

    # Remove any failing files
    if files_to_remove:
        st.session_state.uploaded_files = [
            f for f in st.session_state.uploaded_files if f["name"] not in files_to_remove
        ]

    # Store results
    if has_markers:
        st.session_state.genotype_results_df = pd.DataFrame(results)
    else:
        if results:
            st.session_state.detected_peaks_df = pd.DataFrame(results)

    # Re-run if we removed any
    if files_to_remove:
        st.rerun()

###############################################################################
# 8) Main Script
###############################################################################
def main():
    """
    Main entry point for the Streamlit app: Title, file uploads, advanced options
    (including ploidy), marker configuration, and final processing of files.
    """
    # Toast messages from previous run

    FLAConfigManager.initialize()
    
    if "toast_messages" in st.session_state and st.session_state.toast_messages:
        for msg, icon in st.session_state.toast_messages:
            st.toast(msg, icon=icon)
        st.session_state.toast_messages.clear()

    st.header("FLA Processor Interface")
    tab1, tab2 = st.tabs(["Files", "Markers"])

    # Tab 1 => file management
    with tab1:
        upload_files_col, process_files_col = st.columns(2)
        display_uploaded_files_table()

    FLAConfigManager.update_from_ui()
    fla_config = FLAConfigManager.get()

    # Tab 2 => marker configuration
    with tab2:
        display_marker_configuration()

    # Back to tab1 => upload & analyze
    with tab1:
        with upload_files_col:
            if st.button("Upload Files", use_container_width=True, icon=":material/upload_file:"):
                upload_files_dialog()

        with process_files_col:
            if st.button(
                "Analyze Files",
                disabled=not bool(st.session_state.uploaded_files),
                use_container_width=True,
                type="primary"
            ):
                process_or_load_files(fla_config)

###############################################################################
# 9) Initialize session state + run main
###############################################################################
if "toast_messages" not in st.session_state:
    st.session_state.toast_messages = []

# -----------------------------------------------------------------------------
# Run it
# -----------------------------------------------------------------------------
main()
