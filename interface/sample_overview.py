#!/usr/bin/env python3
__description__ =\
"""
interface/sample_overview.py

"""
__author__ = "Erick Samera"
__version__ = "1.0.0"
__comments__ = ""

import streamlit as st
import pandas as pd

def mark_config_changed():
    st.session_state.config_changed = True

def display_uploaded_files_table():
    if not st.session_state.uploaded_files:
        st.info("Upload some files to get started!")
    else:
        # Build a DataFrame that includes a default Population Group ("Default")
        df_uploaded_files = pd.DataFrame([
            {
                "Filename": f["name"],
                "Sample Name": f.get("sample_name", f["name"].split('_')[0]),
                "Population Group": f.get("population_group", "Default"),
                "Delete": False
            }
            for f in st.session_state.uploaded_files
        ])

        # Create a data editor with the columns:
        edited_df = st.data_editor(
            df_uploaded_files,
            use_container_width=True,
            num_rows="fixed",
            disabled=["Filename"],
            column_config={
                "Delete": "Delete?",
                "Sample Name": "Sample Name",
                "Population Group": "Population Group"
            },
            on_change=mark_config_changed,
        )

        # Update session state based on the edited table:
        for file_dict in st.session_state.uploaded_files:
            matching = edited_df.loc[edited_df["Filename"] == file_dict["name"]]
            if not matching.empty:
                file_dict["sample_name"] = matching["Sample Name"].iloc[0]
                file_dict["population_group"] = matching["Population Group"].iloc[0]

        # Process any deletions:
        to_delete = edited_df[edited_df["Delete"] == True]["Filename"].tolist()
        if to_delete:
            st.session_state.uploaded_files = [
                f for f in st.session_state.uploaded_files if f["name"] not in to_delete
            ]
            st.rerun()
df_uploaded_files = pd.DataFrame([{"Groups": ''}],)
st.data_editor(df_uploaded_files, num_rows='dynamic')

display_uploaded_files_table()
print(st.session_state.uploaded_files)