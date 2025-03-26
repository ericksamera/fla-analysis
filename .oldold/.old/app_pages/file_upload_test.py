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

# --- File Upload Modal ---
@st.dialog("Upload Files", width="large")
def upload_files_dialog():
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
                st.session_state.uploaded_files.append({"name": file.name, "sample_name": file.name, "file": file})
        st.success("Files uploaded successfully!")
        st.rerun()

if st.button("Upload Files"):
    upload_files_dialog()

# Display uploaded files with delete option and sample name input
if st.session_state.uploaded_files:
    st.write("### Uploaded Files")

    # Create DataFrame with UI-friendly column names
    df_uploaded_files = pd.DataFrame([
        {"Filename": f["name"], "Sample Name": f["sample_name"], "delete": False}  # Exclude 'file' object
        for f in st.session_state.uploaded_files
    ])
    
    # Display editable table
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
            file_dict["sample_name"] = edited_df.loc[edited_df["Filename"] == file_dict["name"], "Sample Name"].values[0]

    # Process deletions
    if "delete" in edited_df.columns:
        to_delete = edited_df[edited_df["delete"] == True]["Filename"].tolist()
        if to_delete:
            mark_config_changed()
            st.session_state.uploaded_files = [f for f in st.session_state.uploaded_files if f["name"] not in to_delete]
            st.rerun()
