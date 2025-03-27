import streamlit as st
import pandas as pd

def run():

    def add_population_group():
        new_group = st.session_state.get("new_group_name", "").strip()
        if not new_group:
            return

        group_list = st.session_state.setdefault("population_groups", [])
        if new_group not in group_list:
            group_list.append(new_group)
            st.toast(f"Added group: {new_group}", icon="✨")

        st.session_state["new_group_name"] = ""

    st.title("Sample Metadata Editor")

    if not st.session_state.samples:
        st.warning("No samples found. Please upload files first.")
        return

    with st.expander("Extra options, maybe", expanded=True):
        st.session_state.setdefault("population_groups", ["Pop1", "Pop2", "Unknown"])

        left, right = st.columns([3, 2])

        with left:
            selection = st.pills(
                label="Population Groups",
                options=st.session_state["population_groups"],
                selection_mode="single",
                key="group_pills"
            )
            if selection:
                st.session_state["population_groups"].remove(selection)
                st.toast(f"Removed group: {selection}", icon="🗑️")
                st.rerun()
            st.text_input(
                "New group name",
                key="new_group_name",
                on_change=add_population_group,
                label_visibility="collapsed",
                placeholder="Type name + Enter to add"
            )

        with right:
            pass



    # Aggregate sample metadata
    rows = []
    all_populations = set()
    all_keys = {"Sample Name", "Population"}

    for sample in st.session_state.samples.values():
        sid = sample.sample_id
        meta = sample.metadata.copy()

        meta.setdefault("Sample Name", sid.split("_")[0])
        meta.setdefault("Population", "Unknown")

        all_populations.add(meta["Population"])

        rows.append({
            "Filename": sid + ".fsa",
            **meta
        })

    df = pd.DataFrame(rows)
    column_order = ["Filename", "Sample Name", "Population"]

    # --- Dynamic dropdown for population field ---
    pop_options = st.session_state.get("population_groups", []) or ["Unknown"]
    column_config = {
        "Population": st.column_config.SelectboxColumn(
            label="Population",
            help="Assign population group",
            options=pop_options,
            required=True
        )
    }



    edited = st.data_editor(
        df,
        column_order=column_order,
        column_config=column_config,
        disabled=["Filename"],
        num_rows="fixed",
        use_container_width=True,
        key="metadata_editor"
    )

    if st.button("💾 Save Metadata Changes"):
        for _, row in edited.iterrows():
            sid = row["Filename"].replace(".fsa", "")
            if sid not in st.session_state.samples:
                continue
            st.session_state.samples[sid].metadata = {
                k: row[k] for k in ["Sample Name", "Population"]
                if pd.notna(row[k])
            }
        st.success("Metadata updated.")

run()
