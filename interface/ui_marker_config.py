# interface/ui_marker_config.py

import streamlit as st
import json
from fla_pipeline.config.marker_config import MarkerConfig

@st.dialog("Edit Marker Overrides", width="large")
def edit_marker_overrides_dialog(marker_index: int):
    marker = st.session_state.marker_list[marker_index]

    col_save, col_delete = st.columns([6, 2])

    with col_save:
        st.button("Save", use_container_width=True)

    with col_delete:
        if st.button("Delete Marker", use_container_width=True, type="primary", key=f"delete_marker_{marker_index}"):
            del st.session_state.marker_list[marker_index]
            st.rerun()


@st.dialog("Delete all markers?", width="small")
def clear_markers_dialog():
    st.text("This cannot be undone.")
    st.text("Are you ready to face the consequences?")

    col_yes, col_no = st.columns(2)

    with col_yes:
        if st.button("Yes", use_container_width=True):
            st.session_state.marker_list = []
            st.session_state.config_changed = True
            st.rerun()

    with col_no:
        if st.button("Cancel", type="primary", use_container_width=True):
            st.rerun()

@st.dialog("Upload Marker Config", width="large")
def markers_upload_dialog():
    marker_file = st.file_uploader(
        "Upload a marker config (.json)",
        type=["json"],
        accept_multiple_files=False,
    )
    if marker_file:
        try:
            marker_data = json.loads(marker_file.read().decode("utf-8"))
            st.session_state.marker_list = [MarkerConfig(**entry) for entry in marker_data]
            st.success("Loaded marker config.")
            st.rerun()
        except Exception as e:
            st.error(f"Error parsing marker config: {e}")


def marker_config_ui():
    col_upload, col_download, col_filename = st.columns([1, 1, 4])

    with col_upload:
        if st.button("", use_container_width=True, icon=":material/upload_file:", key="upload_markers"):
            markers_upload_dialog()

    with col_filename:
        marker_config_filename = st.text_input("Filename", value="marker_config.json", label_visibility="collapsed")

    with col_download:
        st.download_button(
            label="",
            data=json.dumps([m.model_dump() for m in st.session_state.marker_list]),
            file_name=marker_config_filename,
            mime="application/json",
            use_container_width=True,
            key="download_markers",
            icon=":material/file_download:" if bool(st.session_state.marker_list) else ":material/file_download_off:",
            disabled=not bool(st.session_state.marker_list)
        )

    with st.container(height=300 if st.session_state.marker_list else 0, border=False):
        for i, marker in enumerate(st.session_state.marker_list):
            cols = st.columns([2, 2, 1, 1, 1, 2], vertical_alignment="bottom")

            marker.marker = cols[0].text_input("Marker", value=marker.marker, key=f"marker_{i}")
            marker.channel = cols[1].selectbox("Channel", ["6-FAM", "VIC", "NED", "PET"],
                                               index=["6-FAM", "VIC", "NED", "PET"].index(marker.channel),
                                               key=f"channel_{i}")
            marker.repeat_unit = cols[2].number_input("Repeat", 1, 10, value=marker.repeat_unit or 1, key=f"repeat_{i}")
            min_bin = cols[3].number_input("Min", 30, 500, value=marker.bins[0], key=f"bmin_{i}")
            max_bin = cols[4].number_input("Max", 30, 500, value=marker.bins[1], key=f"bmax_{i}")
            marker.bins = (min_bin, max_bin)

            with cols[5]:
                if st.button("", key=f"edit_overrides_{i}", icon=":material/settings:", use_container_width=True):
                    edit_marker_overrides_dialog(i)

    # Add/Delete controls
    if st.session_state.marker_list:
        col_add_marker, col_del_marker = st.columns([6, 1])
        with col_add_marker:
            if st.button("Add Marker", icon=":material/add_circle_outline:", use_container_width=True):
                st.session_state.marker_list.append(MarkerConfig(
                    marker=f"Marker{len(st.session_state.marker_list)+1}",
                    channel="6-FAM",
                    repeat_unit=3,
                    bins=(180, 240),
                    overrides={}
                ))
                st.rerun()
        with col_del_marker:
            if st.button("", icon=":material/delete_sweep:", use_container_width=True, type="primary", key="clear_markers"):
                clear_markers_dialog()
    else:
        if st.button("Add Marker", icon=":material/add_circle_outline:", use_container_width=True):
            st.session_state.marker_list.append(MarkerConfig(
                marker="Marker1",
                channel="6-FAM",
                repeat_unit=3,
                bins=(180, 240),
                overrides={}
            ))
            st.rerun()
