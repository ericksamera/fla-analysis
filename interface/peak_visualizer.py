# interface/peak_visualizer.py

import streamlit as st
from fla_pipeline.models.sample import Sample
from interface.plotting.channel_traces import make_total_trace_figure, make_per_channel_figure

class PeakVisualizerApp:
    def __init__(self):
        self.min_height = getattr(st.session_state.config, "min_peak_height", 0)

    def render(self):

        if not st.session_state.samples:
            st.warning("No samples found. Please upload files first.")
            return

        sample_ids = list(st.session_state.samples.keys())
        selected = st.selectbox("Select Sample", sample_ids)
        st.session_state.selected_sample = selected
        self.sample: Sample = st.session_state.samples[selected]

        st.markdown(f"**Sample:** `{self.sample.sample_id}`")

        self._init_sidebar()
        self._plot_total_trace()
        self._plot_per_channel()

    def _init_sidebar(self):
        if not self.sample.marker_results:
            return

        results = self.sample.marker_results
        icon_map = {
            'were removed': ':material/filter_alt:',
            'but within': ':material/auto_fix_normal:',
            'homozygous': ':material/mode:',
        }

        self._display_genotyping()

        if len(results) <= 5:
            tabs = st.sidebar.tabs(results.keys())
            for marker, tab in zip(results.keys(), tabs):
                with tab:
                    self._display_qc(marker, results[marker].qc_flags, icon_map)
        else:
            marker_list = list(results.keys())
            marker_selected = st.sidebar.selectbox("Select Marker", marker_list)
            self._display_qc(marker_selected, results[marker_selected].qc_flags, icon_map)

    def _display_qc(self, marker: str, flags: list[str], icon_map: dict):
        st.markdown(f"**{marker}**")
        if not flags:
            st.success("No QC flags!", icon=":material/check_circle_outline:")
            return
        for flag in flags:
            icon = next((icon for substr, icon in icon_map.items() if substr in flag), ":material/info:")
            st.warning(flag, icon=icon)

    def _display_genotyping(self):
        with st.sidebar.expander("📋 Genotype Summary", expanded=True):
            table_data = []
            marker_configs = {m["marker"]: m for m in st.session_state.marker_list}

            for marker_name, result in self.sample.marker_results.items():
                mcfg = marker_configs.get(marker_name, {})
                table_data.append({
                    "Marker": marker_name,
                    "Dye": mcfg.get("channel", "-"),
                    "Repeat": mcfg.get("repeat_unit", "-"),
                    "Genotype": "/".join(map(str, result.alleles)) if result.alleles else "-",
                    "Conf.": round(result.confidence, 3)
                })

            st.dataframe(table_data, use_container_width=True, hide_index=True)

    def _plot_total_trace(self):
        fig = make_total_trace_figure(
            self.sample,
            st.session_state.marker_list,
            st.session_state.config.model_dump() if hasattr(st.session_state.config, "model_dump") else {}
        )
        st.plotly_chart(fig, use_container_width=True)

    def _plot_per_channel(self):
        channels = self.sample.fsa_data["channels"]

        with st.expander("Per-channel traces"):
            for ch in channels:
                st.code(ch, language=None)
                fig = make_per_channel_figure(
                    ch,
                    self.sample,
                    st.session_state.marker_list,
                    st.session_state.config.model_dump() if hasattr(st.session_state.config, "model_dump") else {}
                )
                st.plotly_chart(fig, use_container_width=True)


# Streamlit entrypoint
def run():
    PeakVisualizerApp().render()

run()