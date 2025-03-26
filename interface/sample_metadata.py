import streamlit as st
import pandas as pd

def run():
    st.title("🧾 Sample Metadata Editor")

    if not st.session_state.samples:
        st.warning("No samples found. Please upload files first.")
        return

    with st.expander("📋 Manage Population Groups", expanded=True):
        new_pop = st.text_input("Add new population group", key="new_pop_group")
        if st.button("➕ Add Group"):
            if "population_groups" not in st.session_state:
                st.session_state["population_groups"] = []
            if new_pop and new_pop not in st.session_state["population_groups"]:
                st.session_state["population_groups"].append(new_pop)
                st.success(f"Added group: {new_pop}")

        if st.session_state.get("population_groups"):
            to_remove = st.multiselect("Remove groups", st.session_state["population_groups"])
            if st.button("🗑️ Remove Selected"):
                st.session_state["population_groups"] = [
                    g for g in st.session_state["population_groups"] if g not in to_remove
                ]


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
