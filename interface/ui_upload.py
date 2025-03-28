# interface/ui_upload.py

import os
import tempfile
import streamlit as st
import pandas as pd
from fla_pipeline.models.sample import Sample
from fla_pipeline.core.fsa_loader import load_fsa
from fla_pipeline.core.peak_detector import detect_peaks


def save_temp_file(uploaded_file):
    suffix = "." + uploaded_file.name.split(".")[-1]
    temp_file = tempfile.NamedTemporaryFile(delete=False, suffix=suffix)
    temp_file.write(uploaded_file.getvalue())
    temp_file.close()
    return temp_file.name

@st.dialog("Delete all files (and associated samples)?", width="small")
def clear_files_dialog():
    """
    A Streamlit dialog to confirm the deletion of uploaded .fsa or .json files.
    Files are stored in st.session_state.uploaded_files.
    """

    st.text("This cannot be undone.")
    st.text("Are you ready to face the consequences?")

    col_yes, col_no = st.columns(2)

    with col_yes:
        if st.button("Yes", use_container_width=True):
            st.session_state.uploaded_files = []
            st.session_state.samples = {}
            st.rerun()

    with col_no:
        if st.button("Cancel", type="primary", use_container_width=True):
            st.rerun()


@st.dialog("Upload Files", width="small")
def file_upload_dialog():
    uploaded = st.file_uploader(
        "Upload one or more `.fsa` files",
        type=["fsa"],
        accept_multiple_files=True,
    )

    if uploaded:
        updated = False
        for f in uploaded:
            temp_path = save_temp_file(f)
            filename_root = os.path.splitext(f.name)[0]
            sample_id = filename_root
            sample_name = filename_root.split("_")[0]

            if not any(existing["name"] == f.name for existing in st.session_state.uploaded_files):
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
                    sample.metadata["Sample Name"] = sample_name

                    st.session_state.samples[sample_id] = sample
                    st.session_state.uploaded_files.append({
                        "name": f.name,
                        "sample_name": sample_name,
                        "path": temp_path
                    })
                    updated = True
                except Exception as e:
                    st.error(f"Error processing {f.name}: {e}")

        if updated:
            st.success(f"Uploaded and parsed {len(uploaded)} new file(s).")
            st.rerun()


def show_uploaded_files_table():
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

    edited = st.data_editor(
        df,
        use_container_width=True,
        num_rows="fixed",
        column_config={"Filename": st.column_config.Column(disabled=True)}
    )

    for i, row in edited.iterrows():
        new_name = row["Sample Name"]
        st.session_state.uploaded_files[i]["sample_name"] = new_name
        sid = os.path.splitext(st.session_state.uploaded_files[i]["name"])[0]
        sample_obj = st.session_state.samples.get(sid)
        if sample_obj:
            sample_obj.metadata["Sample Name"] = new_name

    to_delete = edited[edited["Delete?"]].index.tolist()
    if to_delete:
        st.session_state.uploaded_files = [
            f for i, f in enumerate(st.session_state.uploaded_files) if i not in to_delete
        ]
        for i in to_delete:
            filename = df.iloc[i]["Filename"]
            sample_id = os.path.splitext(filename)[0]
            st.session_state.samples.pop(sample_id, None)
        st.rerun()


def upload_file_ui():
    if st.session_state.uploaded_files:
        col_upload, col_clear = st.columns([6, 1])
        with col_upload:
            if st.button("Upload Files", icon=":material/upload:", use_container_width=True):
                file_upload_dialog()
        with col_clear:
            if st.button("", icon=":material/delete_sweep:", use_container_width=True, type="primary", key="clear_uploaded_files"):
                clear_files_dialog()
    else:
        if st.button("Upload Files", icon=":material/upload:", use_container_width=True):
            file_upload_dialog()

    show_uploaded_files_table()