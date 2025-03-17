#!/usr/bin/env python3
"""
interface/peak_visualizer.py

Purpose: Streamlit wrapper that can either upload or process .fsa files using FLAProcessor,
         or upload preprocessed .json files, and then view them in an interactive trace viewer.
"""
__author__ = "Erick Samera"
__version__ = "1.0.0"
__comments__ = "Works"

import streamlit as st
import plotly.graph_objects as go
import pandas as pd

if "PROCESSED_FLA" not in st.session_state:
    st.session_state.PROCESSED_FLA = {}
if "SELECTED_TRACE" not in st.session_state:
    st.session_state.SELECTED_TRACE = None
if "marker_list" not in st.session_state:
    st.session_state.marker_list = []

class App:
    def __init__(self):
        self.title = "FLA-viewer"
        self.emoji = ":rainbow:"
        st.markdown(
            """
            <style>
            [data-testid="stSidebar"][aria-expanded="true"]{
                min-width: 450px;
                max-width: 450px;
            }
            </style>
            """,
            unsafe_allow_html=True,
        )

    def _init_page(self) -> None:
        """
        Main entry point: sets up the sidebar and the plot window.
        """
        with st.sidebar:
            self._init_sidebar()
        self._update_plot_window()

    def _init_sidebar(self) -> None:
        """
        If processed data is present, show a trace selector and genotype summary.
        Otherwise, display file upload controls to either process a raw FSA file
        (if available) or load a preprocessed JSON file.
        """
        if st.session_state.PROCESSED_FLA:
            st.divider()
            with st.expander("**TRACE FILES:**", expanded=True):
                # Use the dictionary keys as the list of trace names.
                trace_names = list(st.session_state.PROCESSED_FLA.keys())
                st.session_state.SELECTED_TRACE = st.selectbox(
                    "Select trace file to view:",
                    options=trace_names,
                    index=0 if trace_names else 0
                )
            # Display genotype results if available.
            selected_trace = st.session_state.SELECTED_TRACE
            if selected_trace in st.session_state.PROCESSED_FLA:
                trace_data = st.session_state.PROCESSED_FLA[selected_trace]
                if "marker_results" in trace_data and trace_data["marker_results"]:
                    genotype_data = []
                    for marker, result in trace_data["marker_results"].items():
                        genotype_data.append({
                            "Marker": marker,
                            "Channel": result.get("channel", ""),
                            "Genotype": ", ".join(map(str, result.get("parsed_allele_calls", []))),
                            "Confidence": f"{result.get('genotype_confidence', 0):.3f}"
                        })
                    if genotype_data:
                        st.divider()
                        st.markdown("### **Predicted Genotypes**")
                        df_genotypes = pd.DataFrame(genotype_data)
                        st.dataframe(df_genotypes, use_container_width=True, hide_index=True)

    def _hex_to_rgb(self, hex_color: str):
        """
        Converts a hex color string (e.g., "#0000FF") to an RGB tuple (0, 0, 255).
        """
        hex_color = hex_color.lstrip("#")
        return tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))

    def _convert_for_plotting(self, data: dict) -> dict:
        color_mapping = {
            "6-FAM": "blue",
            "VIC": "green",
            "NED": "black",
            "PET": "red",
            "LIZ": "orange"
        }
        processed_dict = {}
        processed_dict["name"] = data.get("file", "Unknown")
        fsa_data = data.get("fsa_data", {})
        processed_dict["smap"] = fsa_data.get("smap", [])
        # Always use the raw channels data from fsa_data:
        processed_dict["channels_peaks"] = fsa_data.get("channels", {})
        
        # Convert detected peaks.
        detected_peaks = data.get("detected_peaks", {})
        peaks_x = {}
        peaks_y = {}
        for channel, peaks in detected_peaks.items():
            if isinstance(peaks, list) and peaks:
                if isinstance(peaks[0], dict):
                    peaks_x[channel] = [p.get("position", 0) for p in peaks]
                    peaks_y[channel] = [p.get("intensity", 0) for p in peaks]
                else:
                    # Assume it's a Peak object; convert using getattr.
                    peaks_x[channel] = [getattr(p, "position", 0) for p in peaks]
                    peaks_y[channel] = [getattr(p, "intensity", 0) for p in peaks]
            else:
                peaks_x[channel] = []
                peaks_y[channel] = []
        processed_dict["detected_peaks_x"] = peaks_x
        processed_dict["detected_peaks_y"] = peaks_y
        processed_dict["colors"] = color_mapping

        if "marker_results" in data:
            processed_dict["marker_results"] = data["marker_results"]

        return processed_dict


    def _plot_total_trace(self, trace_dict: dict) -> go.Figure:
        """
        Plots each channel as a line plus its detected peaks.
        """
        smap = trace_dict.get("smap", [])
        channels_data = trace_dict.get("channels_peaks", {})
        peaks_x = trace_dict.get("detected_peaks_x", {})
        peaks_y = trace_dict.get("detected_peaks_y", {})
        colors = trace_dict.get("colors", {})

        fig = go.Figure()

        all_intensities = []
        for vals in channels_data.values():
            all_intensities.extend(vals)
        max_y = max(all_intensities) if all_intensities else 100
        max_y *= 1.15

        max_x = max(smap) if smap else 100
        max_x += 20

        if "marker_list" in st.session_state and st.session_state.marker_list:
            for marker in st.session_state.marker_list:
                bin_range = marker.get("bins")
                channel = marker.get("channel")
                if bin_range and channel in colors:
                    fig.add_vrect(
                        x0=bin_range[0], x1=bin_range[1],
                        fillcolor=colors[channel], opacity=0.05,
                        layer="below", line_width=0
                    )
                    fig.add_vrect(
                        x0=bin_range[0]-(st.session_state.bin_tolerance*marker.get("repeat_unit")),
                        x1=bin_range[1]+(st.session_state.bin_tolerance*marker.get("repeat_unit")),
                        fillcolor=colors[channel], opacity=0.05,
                        layer="below", line_width=0
                    )

        for ch_name, intensity_list in channels_data.items():
            fig.add_trace(go.Scatter(
                x=smap,
                y=intensity_list,
                mode="lines",
                name=ch_name,
                fill="tozeroy",
                hoverinfo="skip",
                marker=dict(color=colors.get(ch_name, "gray"))
            ))
            if ch_name == "LIZ":
                continue
            fig.add_trace(go.Scatter(
                x=peaks_x.get(ch_name, []),
                y=peaks_y.get(ch_name, []),
                mode="markers",
                name=f"{ch_name}",
                showlegend=False,
                hoverinfo="x+y+text+name",
                hovertemplate="size (bp): %{x}<br>height: %{y}<br>",
                marker=dict(color=colors.get(ch_name, "gray"))
            ))

        fig.update_layout(
            xaxis=dict(range=[0, max_x], title="Size (bp)"),
            yaxis=dict(range=[0, max_y], title="Intensity"),
            margin=dict(l=0, r=0, t=0, b=0),
            legend_title_text="Channels",
            dragmode="zoom"
        )
        return fig

    def _plot_per_channel_traces(self, trace_dict: dict) -> None:
        """
        Creates an expander with individual channel sub-plots and optional peak tables.
        """
        smap = trace_dict.get("smap", [])
        channels_data = trace_dict.get("channels_peaks", {})
        peaks_x = trace_dict.get("detected_peaks_x", {})
        peaks_y = trace_dict.get("detected_peaks_y", {})
        colors = trace_dict.get("colors", {})

        with st.expander("Individual channels"):
            st.session_state.FULL_GENOTYPES = st.checkbox("Show all detected peaks", value=False)
            st.session_state.SHOW_INDIV_TABLES = st.checkbox("Show table with individual channels", value=False)

            for ch_name, ch_color in colors.items():
                x_vals = peaks_x.get(ch_name, [])
                y_vals = peaks_y.get(ch_name, [])
                valid_data = [(x, y) for x, y in zip(x_vals, y_vals) if x > 50]

                if not valid_data and not channels_data.get(ch_name):
                    continue

                fig = go.Figure()

                if "marker_list" in st.session_state and st.session_state.marker_list:
                    for marker in st.session_state.marker_list:
                        if marker.get("channel") != ch_name: continue
                        bin_range = marker.get("bins")
                        if bin_range:
                            fig.add_vrect(
                                x0=bin_range[0], x1=bin_range[1],
                                fillcolor=ch_color, opacity=0.05,
                                layer="below", line_width=0
                            )
                            fig.add_vrect(
                                x0=bin_range[0]-(st.session_state.bin_tolerance*marker.get("repeat_unit")),
                                x1=bin_range[1]+(st.session_state.bin_tolerance*marker.get("repeat_unit")),
                                fillcolor=ch_color, opacity=0.05,
                                layer="below", line_width=0
                            )

                fig.add_trace(go.Scatter(
                    x=smap,
                    y=channels_data.get(ch_name, []),
                    mode="lines",
                    name=ch_name,
                    fill="tozeroy",
                    hoverinfo="skip",
                    marker=dict(color=ch_color)
                ))

                if valid_data:
                    dx = [p[0] for p in valid_data]
                    dy = [p[1] for p in valid_data]
                    fig.add_trace(go.Scatter(
                        x=dx,
                        y=dy,
                        mode="markers",
                        name=f"{ch_name}",
                        showlegend=False,
                        hovertemplate="size (bp): %{x}<br>height: %{y}<br>",
                        marker=dict(color=ch_color, size=3)
                    ))

                max_y = max([y for _, y in valid_data]) * 1.25 if valid_data else 100
                max_x = max([x for x, _ in valid_data]) * 1.05 if valid_data else 100

                fig.update_layout(
                    xaxis=dict(range=[0, max_x], title="Size (bp)"),
                    yaxis=dict(range=[0, max_y], title="Intensity"),
                    height=200,
                    margin=dict(l=0, r=0, t=0, b=0),
                )

                st.markdown(f"### {ch_name}")
                st.plotly_chart(fig, use_container_width=True)

                if st.session_state.SHOW_INDIV_TABLES and valid_data:
                    table_rows = [{"size (bp)": x, "height": y} for x, y in valid_data]
                    avg_h = sum(row["height"] for row in table_rows) / len(table_rows)
                    if not st.session_state.FULL_GENOTYPES:
                        table_rows = [row for row in table_rows if row["height"] > avg_h]
                    st.dataframe(table_rows, use_container_width=True)
                st.divider()

    def _update_plot_window(self) -> None:
        """
        Shows the main trace and per-channel details for the selected trace.
        """
        if not st.session_state.SELECTED_TRACE:
            st.info("No trace selected. Please upload a file to visualize peaks.")
            return

        st.header(f"Viewing: {st.session_state.SELECTED_TRACE}")
        trace_key = st.session_state.SELECTED_TRACE
        trace_dict = st.session_state.PROCESSED_FLA.get(trace_key)
        if not trace_dict:
            st.error("Selected trace not found in session.")
            return

        # Ensure trace is in correct format
        if "smap" not in trace_dict:
            trace_dict = self._convert_for_plotting(trace_dict)

        fig = self._plot_total_trace(trace_dict)
        st.plotly_chart(fig, use_container_width=True)
        self._plot_per_channel_traces(trace_dict)

def run():
    """
    Entry point for the peak visualizer.
    """
    app = App()
    app._init_page()

run()
