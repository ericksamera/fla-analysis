#!/usr/bin/env python3
__description__ =\
"""
interface/file_processor.py

Implements a page in the Streamlit UI that handles:
 - File uploads (.fsa, .json)
 - Marker configuration
 - Advanced settings (diploid vs. polyploid, thresholds, etc.)
 - Processing (calling into core.fla_processor or core.fla_processor_poly)
"""
__author__ = "Erick Samera"
__version__ = "1.0.0"
__comments__ = "Works"

import streamlit as st
import json
import os
import pandas as pd

# Import your core modules:
from core.fla_processor import FLAProcessor, FLAConfig
from core.fla_processor_poly import FLAPolyProcessor, FLAConfigPoly
from core.config_models import MarkerConfig

###############################################################################
# 1) Session-State Helpers / Initialization
###############################################################################

def initialize_session_state():
    # Toast messages for success/error notifications
    if "toast_messages" not in st.session_state:
        st.session_state.toast_messages = []

def mark_config_changed():
    """Flag that config changed so we can track re-analysis or UI refresh needs."""
    st.session_state.config_changed = True

def add_toast_message(msg: str, icon: str = None):
    """
    Appends a toast message with an optional icon to session_state,
    to be displayed after the next rerun.
    """
    st.session_state.toast_messages.append((msg, icon))

###############################################################################
# 2) UI Components: Upload Dialog, File Table, Marker Config, Advanced Options
###############################################################################

@st.dialog("Delete all markers?", width="small")
def clear_markers_dialog():
    """
    A Streamlit dialog to confirm the deletion of uploaded .fsa or .json files.
    Files are stored in st.session_state.uploaded_files.
    """

    st.text("This cannot be undone.")
    st.text("Are you ready to face the consequences?")

    col_yes, col_no = st.columns(2)

    with col_yes:
        if st.button("Yes", use_container_width=True):
            st.session_state.config_json = {}
            st.session_state.marker_list = []
            st.rerun()
    
    with col_no:
        if st.button("Cancel", type="primary", use_container_width=True):
            st.rerun()

@st.dialog("Delete all files?", width="small")
def clear_files_dialog():
    """
    A Streamlit dialog to confirm the deletion of uploaded .fsa or .json files.
    Files are stored in st.session_state.uploaded_files.
    """

    st.text("This cannot be undone.")
    st.text("Are you ready to face the consequences?")

    col_yes, col_no = st.columns(2)

    file_manager_defaults = {
        "genotype_results_df": pd.DataFrame(),
        "detected_peaks_df": pd.DataFrame(),
        "PROCESSED_FLA": {},
        "uploaded_files": [],
        "uploaded_files_id_counter": 0,
    }

    with col_yes:
        if st.button("Yes", use_container_width=True):
            for key, value in file_manager_defaults.items():
                st.session_state[key] = value
            st.rerun()
    
    with col_no:
        if st.button("Cancel", type="primary", use_container_width=True):
            st.rerun()

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
                    "sample_name": file.name.split('_')[0],
                    "file": file
                })
        st.success("Files uploaded successfully!")
        st.rerun()

def display_uploaded_files_table():
    """
    Displays a table with the uploaded files, allowing rename or delete.
    """
    if not st.session_state.uploaded_files:
        st.info("Upload some files to get started!")
    else:
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
                for fname in to_delete:
                    add_toast_message(f"Deleting '{fname}' from upload list.", ":material/delete:")
                st.session_state.uploaded_files = [
                    f for f in st.session_state.uploaded_files if f["name"] not in to_delete
                ]
                st.rerun()

def display_marker_configuration():
    """
    Renders the UI for marker configuration: upload config, add/edit/remove markers,
    and optionally save the configuration.
    """
    st.write("#### Marker Configuration")
    config_file = st.file_uploader("Upload a JSON config file (optional)", type=["json"])

    col_right, col_save, col_dl = st.columns([4, 2, 2])

    if config_file:
        try:
            config_data = json.load(config_file)
            st.session_state.marker_list = config_data
            st.success("Configuration loaded successfully!")
        except Exception as exc:
            st.error(f"Error reading marker config: {exc}")

    for i, marker in enumerate(st.session_state.marker_list):
        marker_cols = st.columns([2, 2, 1, 1, 1, 2])  # Adjusted column widths

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
            marker["repeat_unit"] = st.number_input(
                "Repeat", min_value=1, max_value=10,
                value=marker["repeat_unit"], key=f"repeat_unit_{i}"
            )

        with marker_cols[3]:
            allele_min = st.number_input(
                "Min", min_value=50, max_value=500,
                value=marker["bins"][0], key=f"allele_min_{i}"
            )

        with marker_cols[4]:
            allele_max = st.number_input(
                "Max", min_value=50, max_value=500,
                value=marker["bins"][1], key=f"allele_max_{i}"
            )

        with marker_cols[5]:
            st.text("")
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
                if st.button(f"Remove ({marker['marker']})", key=f"remove_marker_{i}"):
                    removed_name = marker["marker"]
                    st.session_state.marker_list.pop(i)
                    add_toast_message(f"Removed marker '{removed_name}' from configuration.", ":material/delete:")
                    st.rerun()

    if st.session_state.marker_list:
        col_add_marker, col_del_marker = st.columns([6, 1])
        with col_add_marker:
            if st.button("Add Marker", icon=":material/add_circle_outline:", use_container_width=True):
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
                add_toast_message(f"Added new marker '{new_marker['marker']}'.")
                st.rerun()
        with col_del_marker:
            st.button("", icon=":material/delete_sweep:", use_container_width=True, type="primary", on_click=clear_markers_dialog)
    else:
        if st.button("Add Marker", icon=":material/add_circle_outline:", use_container_width=True):
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
            add_toast_message(f"Added new marker '{new_marker['marker']}'.")
            st.rerun()

    with col_save:
        if st.button("", icon=":material/save:", use_container_width=True):
            st.session_state.config_json = st.session_state.marker_list

    with col_right:
        marker_config_filename = st.text_input("filename", value="marker_config.json", label_visibility="collapsed", disabled=not bool(st.session_state.config_json))       

    with col_dl:
        st.download_button(
            label="",
            data=json.dumps(st.session_state.config_json, indent=4),
            file_name=marker_config_filename,
            mime="application/json",
            use_container_width=True,
            disabled=not bool(st.session_state.config_json),
            icon=":material/file_download:" if st.session_state.config_json else ":material/file_download_off:"
        )


def display_advanced_options():
    """
    Renders advanced analysis options, including ploidy.
    Saves the advanced config values to st.session_state so that they persist between reruns.
    Returns a tuple: (is_polyploid, chosen_ploidy, advanced_config_dict)
    """
    
    
    # Initialize defaults if not already set
    if "advanced_config" not in st.session_state:
        st.session_state.advanced_config = {
            "bin_tolerance": 2,
            "min_peak_height": 1000,
            "min_peak_position": 50,
            "relative_peak_threshold": 0.1,
            "stutter_ratio_threshold": 0.15,
            "stutter_tolerance": 1.0,
            "ploidy": 2
        }

    with st.expander("Additional Thresholds"):
        col_ploidy_opt, col_ploidy_n = st.columns([2, 8])
        with col_ploidy_opt:
            ploidy_mode = st.segmented_control("Analysis Mode", ["Diploid", "Polyploid"], selection_mode="single", default="Diploid", on_change=mark_config_changed)
        is_polyploid = (ploidy_mode == "Polyploid")

        if is_polyploid:
            with col_ploidy_n:
                st.session_state.ploidy = st.slider(
                    "Select Ploidy Level",
                    min_value=3,
                    max_value=8,
                    value=st.session_state.ploidy if st.session_state.ploidy != 2 else 3,
                    on_change=mark_config_changed
                )
        else:
            st.session_state.ploidy = 2
        bin_tolerance = st.slider(
            "Bin Tolerance", 0, 10, st.session_state.advanced_config["bin_tolerance"],
            on_change=mark_config_changed)
        min_peak_height = st.slider(
            "Minimum Peak Height",
            10, 5000, st.session_state.advanced_config["min_peak_height"],
            on_change=mark_config_changed)
        min_peak_position = st.slider(
            "Minimum Peak Position",
            10, 1000, st.session_state.advanced_config["min_peak_position"],
            on_change=mark_config_changed)
        relative_peak_threshold = st.slider(
            "Relative Peak Threshold",
            0.0, 1.0, st.session_state.advanced_config["relative_peak_threshold"],
            on_change=mark_config_changed)
        stutter_ratio_threshold = st.slider(
            "Default Stutter Ratio",
            0.0, 1.0, st.session_state.advanced_config["stutter_ratio_threshold"],
            on_change=mark_config_changed)
        stutter_tolerance = st.slider(
            "Default Stutter Tolerance",
            0.0, 5.0, st.session_state.advanced_config["stutter_tolerance"],
            on_change=mark_config_changed)

        # Save the updated values into session state
        st.session_state.advanced_config = {
            "bin_tolerance": bin_tolerance,
            "min_peak_height": min_peak_height,
            "min_peak_position": min_peak_position,
            "relative_peak_threshold": relative_peak_threshold,
            "stutter_ratio_threshold": stutter_ratio_threshold,
            "stutter_tolerance": stutter_tolerance
        }
        # Also store bin_tolerance separately if other parts of the code rely on it
        st.session_state.bin_tolerance = bin_tolerance

    return is_polyploid, st.session_state.ploidy, st.session_state.advanced_config


###############################################################################
# 3) Processor Selection & File Processing
###############################################################################

def choose_processor_and_config(is_polyploid: bool, ploidy: int, adv_cfg: dict):
    """
    Based on user input, build the correct config and processor class from the core modules.
    """
    if not is_polyploid:
        # Diploid
        config = FLAConfig(
            bin_tolerance=adv_cfg["bin_tolerance"],
            min_peak_height=adv_cfg["min_peak_height"],
            min_peak_position=adv_cfg["min_peak_position"],
            relative_peak_threshold=adv_cfg["relative_peak_threshold"],
            stutter_ratio_threshold=adv_cfg["stutter_ratio_threshold"],
            stutter_tolerance=adv_cfg["stutter_tolerance"],
            ploidy=2
        )
        processor_class = FLAProcessor
    else:
        # Polyploid
        config = FLAConfigPoly(
            bin_tolerance=adv_cfg["bin_tolerance"],
            min_peak_height=adv_cfg["min_peak_height"],
            min_peak_position=adv_cfg["min_peak_position"],
            relative_peak_threshold=adv_cfg["relative_peak_threshold"],
            stutter_ratio_threshold=adv_cfg["stutter_ratio_threshold"],
            stutter_tolerance=adv_cfg["stutter_tolerance"],
            ploidy=ploidy
        )
        processor_class = FLAPolyProcessor
    
    return processor_class, config

def build_marker_config_list():
    """
    Converts st.session_state.marker_list (plain dicts) into a list or dict
    of MarkerConfig objects expected by the processors.
    """
    markers_dict = {}
    for m in st.session_state.marker_list:
        mc = MarkerConfig(
            marker=m["marker"],
            channel=m["channel"],
            repeat_unit=m["repeat_unit"],
            bins=tuple(m["bins"]),
            stutter_ratio_threshold=m["stutter_ratio_threshold"],
            stutter_tolerance=m["stutter_tolerance"],
            stutter_direction=m["stutter_direction"]
        )
        # store it under the marker name as the key
        markers_dict[m["marker"]] = {
            "channel": mc.channel,
            "bins": mc.bins,
            "repeat_unit": mc.repeat_unit,
            "stutter_ratio_threshold": mc.stutter_ratio_threshold,
            "stutter_tolerance": mc.stutter_tolerance,
            "stutter_direction": mc.stutter_direction
        }
    return markers_dict

def build_marker_config_list_poly():
    """
    For polyploid usage, we might return a simple list of MarkerConfig objects.
    """
    marker_list = []
    for m in st.session_state.marker_list:
        mc = MarkerConfig(
            marker=m["marker"],
            channel=m["channel"],
            repeat_unit=m["repeat_unit"],
            bins=tuple(m["bins"]),
            stutter_ratio_threshold=m["stutter_ratio_threshold"],
            stutter_tolerance=m["stutter_tolerance"],
            stutter_direction=m["stutter_direction"]
        )
        marker_list.append(mc)
    return marker_list

def process_single_fsa(file_path: str, processor_class, config, marker_config, is_poly: bool):
    """
    Given a local file_path to an .fsa, constructs an instance of
    the chosen processor and runs it.
    """
    if is_poly:
        # Polyploid usage
        proc = processor_class(config, marker_config)  # marker_config is a list
        success = proc.process_all(file_path)
        if not success:
            return None
        results = {
            "fsa_data": proc.fsa_data,
            "detected_peaks": proc.detected_peaks,
            "marker_results": proc.marker_results
        }
        return results
    else:
        # Diploid usage
        proc = processor_class(config, marker_config)  # marker_config is a dict
        results = proc.process_all(file_path)
        return results

def process_single_json(file_obj):
    """Simply parse the JSON to mimic the same structure as an FSA-processed result."""
    try:
        return json.load(file_obj)
    except Exception as exc:
        return None  # We'll handle toasts outside

def process_all_files(processor_class, config, is_poly: bool):
    """
    Loop over st.session_state.uploaded_files and process each one,
    removing any that fail, storing results in st.session_state.
    """
    if not st.session_state.uploaded_files:
        st.warning("No files uploaded to process.")
        return

    if is_poly:
        marker_config_list = build_marker_config_list_poly()
    else:
        marker_config_list = build_marker_config_list()

    genotype_rows = []
    files_to_remove = []

    for file_dict in st.session_state.uploaded_files:
        file_obj = file_dict["file"]
        file_name = file_dict["name"]
        _, ext = os.path.splitext(file_name.lower())
        local_path = os.path.join("/tmp", file_name)

        # Write locally
        with open(local_path, "wb") as f:
            f.write(file_obj.getbuffer())

        if ext == ".fsa":
            results = process_single_fsa(local_path, processor_class, config, marker_config_list, is_poly)
            if not results["fsa_data"]["smap"]:
                files_to_remove.append(file_name)
                continue

            if "marker_results" in results:
                for mk_name, mk_data in results["marker_results"].items():
                    genotype = mk_data.get("parsed_allele_calls") or mk_data.get("alleles")
                    conf = mk_data.get("genotype_confidence", None)
                    flags = mk_data.get("QC_flags", [])
                    genotype_rows.append({
                        "Filename": file_dict["name"],
                        "Sample Name": file_dict["sample_name"],
                        "Marker": mk_name,
                        "Genotype": genotype,
                        "Genotype Confidence": f"{conf:.3f}" if conf is not None else "n/a",
                        "QC Flags": ", ".join(flags) if flags else "None"
                    })

            st.session_state.PROCESSED_FLA[file_name] = results

        elif ext == ".json":
            loaded_data = process_single_json(file_obj)
            if not loaded_data:
                add_toast_message(f"Removing '{file_name}' (invalid JSON).", ":material/error:")
                files_to_remove.append(file_name)
                continue

            st.session_state.PROCESSED_FLA[file_name] = loaded_data
            st.success(f"JSON loaded: {file_name}")
            if "marker_results" in loaded_data:
                for mk_name, mk_data in loaded_data["marker_results"].items():
                    genotype = mk_data.get("parsed_allele_calls") or mk_data.get("alleles")
                    conf = mk_data.get("genotype_confidence", None)
                    flags = mk_data.get("QC_flags", [])
                    genotype_rows.append({
                        "Filename": file_dict["name"],
                        "Sample Name": file_dict["sample_name"],
                        "Marker": mk_name,
                        "Genotype": genotype,
                        "Genotype Confidence": f"{conf:.3f}" if conf is not None else "n/a",
                        "QC Flags": ", ".join(flags) if flags else "None"
                    })
        else:
            pass
        for fname in files_to_remove:
            add_toast_message(f"Processing '{fname}' failed! Removed from uploaded files list.", ":material/error:")
        st.session_state.uploaded_files = [
            f for f in st.session_state.uploaded_files if f["name"] not in files_to_remove
        ]

    st.session_state.genotype_results_df = pd.DataFrame(genotype_rows)

###############################################################################
# 4) Main Page Entry Function
###############################################################################

def run():
    """
    Main Streamlit page function for uploading & processing files.
    """
    initialize_session_state()

    if st.session_state.toast_messages:
        for msg, icon in st.session_state.toast_messages:
            st.toast(msg, icon=icon)
        st.session_state.toast_messages.clear()

    st.title("Experiment Setup")
    st.write("Upload `.fsa` files to run diploid/polyploid analysis, or load `.json` results.")

    col_left, col_right = st.columns([2,2])
    with col_left:
        st.write("#### File Manager")
        if not st.session_state.uploaded_files:
            st.button("Upload Files", use_container_width=True, icon=":material/upload_file:", on_click=upload_files_dialog)
        else:
            upload_col, delete_all_col = st.columns(2)
            with upload_col:
                st.button("Upload Files", use_container_width=True, icon=":material/upload_file:", on_click=upload_files_dialog)
            with delete_all_col:
                st.button("Delete All Files", type="primary", use_container_width=True, icon=":material/delete_sweep:", on_click=clear_files_dialog)

        display_uploaded_files_table()

    with col_right:
        display_marker_configuration()

    # 4b) Marker config in its own section
    is_poly, chosen_ploidy, adv_cfg = display_advanced_options()

    # 4c) Process button
    st.write("---")
    if st.button("Process / Analyze All Files", type="primary", use_container_width=True):
        processor_class, config_obj = choose_processor_and_config(is_poly, chosen_ploidy, adv_cfg)
        process_all_files(processor_class, config_obj, is_poly)
        st.session_state.config_changed = False
        st.rerun()

    if not st.session_state.genotype_results_df.empty:
        st.subheader("Genotype Results")
        st.dataframe(st.session_state.genotype_results_df, use_container_width=True, column_order=(column for column in st.session_state.genotype_results_df.columns if column != 'Filename'))

run()