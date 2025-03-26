# interface/peak_visualizer.py

import streamlit as st
import plotly.graph_objects as go
from fla_pipeline.models.sample import Sample

class PeakVisualizerApp:
    def __init__(self):
        self.color_map = {
            "6-FAM": "blue",
            "VIC": "green",
            "NED": "black",
            "PET": "red",
            "LIZ": "orange"
        }
        self.min_height = getattr(st.session_state.config, "min_peak_height", 0)

    def render(self):
        st.title("📈 Peak Visualizer")

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

        st.sidebar.header("🧪 QC Flags")
        results = self.sample.marker_results
        icon_map = {
            'were removed': ':material/filter_alt:',
            'but within': ':material/auto_fix_normal:',
            'homozygous': ':material/mode:',
        }

        if len(results) <= 5:
            tabs = st.sidebar.tabs(results.keys())
            for marker, tab in zip(results.keys(), tabs):
                with tab:
                    self._display_qc(marker, results[marker].qc_flags, icon_map)
        else:
            marker_list = list(results.keys())
            marker_selected = st.sidebar.selectbox("Select Marker", marker_list)
            self._display_qc(marker_selected, results[marker_selected].qc_flags, icon_map)

        # --- Genotype Summary Table ---
        with st.sidebar.expander("📋 Genotype Summary", expanded=True):
            table_data = []

            # Create quick lookup by marker name
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

    def _display_qc(self, marker: str, flags: list[str], icon_map: dict):
        st.markdown(f"**{marker}**")
        if not flags:
            st.success("No QC flags!", icon=":material/check_circle_outline:")
            return
        for flag in flags:
            icon = next((icon for substr, icon in icon_map.items() if substr in flag), ":material/info:")
            st.warning(flag, icon=icon)

    def _plot_total_trace(self):
        smap = self.sample.fsa_data["smap"]
        channels = self.sample.fsa_data["channels"]
        peaks = self.sample.peaks

        fig = go.Figure()

        # Bin and tolerance overlays
        if st.session_state.marker_list:
            for marker in st.session_state.marker_list:
                channel = marker.get("channel")
                repeat = marker.get("repeat_unit", 1)
                tol = st.session_state.get("bin_tolerance", 2) * repeat
                bmin, bmax = marker.get("bins", [0, 0])
                color = self._get_color(channel)

                fig.add_vrect(x0=bmin, x1=bmax, fillcolor=color, opacity=0.05, layer="below", line_width=0)
                fig.add_vrect(x0=bmin - tol, x1=bmax + tol, fillcolor=color, opacity=0.05, layer="below", line_width=0)

        # Channel traces and peaks
        for ch, trace in channels.items():
            color = self._get_color(ch)
            fig.add_trace(go.Scatter(x=smap, y=trace, mode="lines", name=f"{ch} trace", line=dict(color=color), hoverinfo="skip"))

            if ch != "LIZ" and ch in peaks:
                fig.add_trace(go.Scatter(
                    x=[p.position for p in peaks[ch]],
                    y=[p.intensity for p in peaks[ch]],
                    mode="markers",
                    name=f"{ch} peaks",
                    marker=dict(color=color, size=6),
                    showlegend=False,
                    hovertemplate="Size: %{x} bp<br>Height: %{y}"
                ))

        fig.update_layout(
            title="Electropherogram Trace",
            xaxis_title="Size (bp)",
            yaxis_title="Intensity",
            margin=dict(t=30, b=30),
            height=400,
            legend_title="Channels",
            dragmode="zoom"
        )
        st.plotly_chart(fig, use_container_width=True)

    def _plot_per_channel(self):
        st.subheader("🧬 Per-Channel Views")

        smap = self.sample.fsa_data["smap"]
        channels = self.sample.fsa_data["channels"]
        peaks = self.sample.peaks

        for ch, trace in channels.items():
            color = self._get_color(ch)
            ch_peaks = peaks.get(ch, [])

            fig = go.Figure()
            fig.add_trace(go.Scatter(
                x=smap, y=trace, mode="lines", name=f"{ch} trace",
                line=dict(color=color), hoverinfo="skip"
            ))


            if st.session_state.marker_list:
                for marker in st.session_state.marker_list:
                    channel = marker.get("channel")
                    if not channel == ch: continue
                    repeat = marker.get("repeat_unit", 1)
                    tol = st.session_state.get("bin_tolerance", 2) * repeat
                    bmin, bmax = marker.get("bins", [0, 0])
                    color = self._get_color(channel)

                    fig.add_vrect(x0=bmin, x1=bmax, fillcolor=color, opacity=0.05, layer="below", line_width=0)
                    fig.add_vrect(x0=bmin - tol, x1=bmax + tol, fillcolor=color, opacity=0.05, layer="below", line_width=0)

            visible_peaks = ch_peaks if ch == "LIZ" else [p for p in ch_peaks if p.intensity >= self.min_height]

            if visible_peaks:
                # Annotate each peak with any marker it falls within
                annotations = []
                for p in visible_peaks:
                    matched = []
                    for marker in st.session_state.marker_list:
                        if marker["channel"] != ch:
                            continue
                        bmin, bmax = marker.get("bins", [0, 0])
                        repeat = marker.get("repeat_unit", 1)
                        tol = st.session_state.get("bin_tolerance", 2) * repeat
                        if (bmin - tol) <= p.position <= (bmax + tol):
                            matched.append(marker["marker"])
                    marker_note = ", ".join(matched) if matched else "—"
                    annotations.append(f"Size: {p.position} bp<br>Height: {p.intensity}<br>Marker: {marker_note}")

                fig.add_trace(go.Scatter(
                    x=[p.position for p in visible_peaks],
                    y=[p.intensity for p in visible_peaks],
                    mode="markers",
                    name=f"({ch})",
                    marker=dict(color=color, size=6),
                    showlegend=False,
                    hoverinfo="text",
                    text=annotations
                ))

            marker_bins = [
                (m.get("bins", [0, 0]), m.get("repeat_unit", 1))
                for m in st.session_state.marker_list if m.get("channel") == ch
            ]

            if marker_bins:
                min_start = min(bmin for (bmin, _), _ in marker_bins)
                max_end = max(bmax for (_, bmax), _ in marker_bins)
                max_repeat = max(rpt for _, rpt in marker_bins)
                pad = 10 * max_repeat
                x_range = [max(0, min_start - pad), max_end + pad]

            elif ch_peaks:
                # Fallback: use actual detected peaks
                peak_positions = [p.position for p in ch_peaks]
                pad = 10  # or make configurable
                x_range = [max(0, min(peak_positions) - pad), max(peak_positions) + pad]

            else:
                # Fallback: full smap
                x_range = [min(smap), max(smap)]

            fig.update_layout(
                xaxis_title="Size (bp)",
                yaxis_title="Intensity",
                height=200,
                margin=dict(t=5, b=5),
                xaxis_range=x_range
            )

            st.plotly_chart(fig, use_container_width=True)

    def _get_color(self, channel: str) -> str:
        return self.color_map.get(channel, "gray")


# Streamlit entrypoint
def run():
    PeakVisualizerApp().render()

run()