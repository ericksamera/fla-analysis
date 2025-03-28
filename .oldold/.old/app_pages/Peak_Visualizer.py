#!/usr/bin/env python3
"""
Purpose: Streamlit wrapper that can either upload or process .fsa files using FLAProcessor,
         or upload preprocessed .json files, and then view them in an interactive trace viewer.

Author: Erick Samera
Version: 1.2.1
Comments: stable enough; altered zoom
"""
import streamlit as st
import plotly.graph_objects as go
import pandas as pd


try:
    from _fla_processor import FLAProcessor, FLAConfig, MarkerConfig
    HAVE_PROCESSOR = True
except ImportError:
    st.warning("Could not import _fla_processor. Processing .fsa files will not work.")
    HAVE_PROCESSOR = False

class App:
    def __init__(self):
        self.title = "FLA-viewer"
        self.emoji = ':rainbow:'
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
        Main entry point: sets up the sidebar, checks session_state, shows either
        the upload form or the viewer depending on what data we have.
        """
        with st.sidebar: self._init_sidebar()
        self._update_plot_window()

    def _init_sidebar(self) -> None:
        if st.session_state.PROCESSED_FLA:
            st.divider()
            with st.expander('**TRACE FILES:**', expanded=True):
                # Sort the processed list by name
                trace_list = list(st.session_state.PROCESSED_FLA.values())
                trace_names = [t['name'] for t in trace_list]
                st.session_state.SELECTED_TRACE = st.selectbox(
                    'Select trace file to view:',
                    options=trace_names
                )

            # Display genotype results if available
            selected_trace = st.session_state.SELECTED_TRACE
            if selected_trace in st.session_state.PROCESSED_FLA:
                trace_data = st.session_state.PROCESSED_FLA[selected_trace]
                if "marker_results" in trace_data:
                    genotype_data = []
                    for marker, result in trace_data["marker_results"].items():
                        genotype_data.append({
                            "Marker": marker,
                            "Channel": result["channel"],
                            "Genotype": ", ".join(map(str, result["parsed_allele_calls"])),
                            "Confidence": f"{result['genotype_confidence']:.3f}"
                        })
                    
                    if genotype_data:
                        st.divider()
                        st.markdown("### **Predicted Genotypes**")
                        df_genotypes = pd.DataFrame(genotype_data)
                        st.dataframe(df_genotypes, use_container_width=True)


    def _hex_to_rgb(self, hex_color: str):
        """
        Converts a hex color string (e.g., "#0000FF") to an RGB tuple (0, 0, 255).
        """
        hex_color = hex_color.lstrip("#")
        return tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))


    def _convert_for_plotting(self, data: dict) -> dict:
        """
        Converts raw data from FLAProcessor (or user-supplied JSON) to the structure
        needed by the viewer: 'smap', 'channels_peaks', 'detected_peaks_x', etc.
        """
        color_mapping = {
            '6-FAM': "blue",
            'VIC':   "green",
            'NED':   "black",
            'PET':   "red",
            'LIZ':   "orange"
        }
        processed_dict = {}
        processed_dict['name'] = data.get("file", "Unknown")
        fsa_data = data.get("fsa_data", {})
        processed_dict['smap'] = fsa_data.get("smap", [])
        processed_dict['channels_peaks'] = fsa_data.get("channels", {})

        # Build x,y for detected peaks
        detected_peaks = data.get("detected_peaks", {})
        peaks_x = {}
        peaks_y = {}
        for channel, peaks in detected_peaks.items():
            peaks_x[channel] = [p["position"] for p in peaks]
            peaks_y[channel] = [p["intensity"] for p in peaks]

        processed_dict['detected_peaks_x'] = peaks_x
        processed_dict['detected_peaks_y'] = peaks_y
        processed_dict['colors'] = color_mapping

        return processed_dict

    def _plot_total_trace(self, trace_dict: dict) -> go.Figure:
        """
        Similar to your original "plot_total_trace" approach:
        Plots each channel as a line + its detected peaks.
        """
        smap = trace_dict.get('smap', [])
        channels_data = trace_dict.get('channels_peaks', {})
        peaks_x = trace_dict.get('detected_peaks_x', {})
        peaks_y = trace_dict.get('detected_peaks_y', {})
        colors = trace_dict.get('colors', {})

        fig = go.Figure()

        # Find a suitable max for scaling
        all_ints = []
        for chvals in channels_data.values():
            all_ints.extend(chvals)
        max_y = max(all_ints) if all_ints else 100
        max_y *= 1.15

        max_x = max(smap) if smap else 100
        max_x += 20

        if "marker_list" in st.session_state:
            for marker in st.session_state.marker_list:
                bin_range = marker["bins"]
                channel = marker["channel"]

                if bin_range and channel in colors:
                    fig.add_vrect(
                        x0=bin_range[0], x1=bin_range[1],
                        fillcolor=colors[channel], opacity=0.1,
                        layer="below", line_width=0
                    )


        for ch_name, intensity_list in channels_data.items():
            fig.add_trace(go.Scatter(
                x=smap,
                y=intensity_list,
                mode='lines',
                name=ch_name,
                fill='tozeroy',
                hoverinfo='skip',
                marker=dict(color=colors.get(ch_name, "gray"))
            ))
            if ch_name == "LIZ": continue
            fig.add_trace(go.Scatter(
                x=peaks_x.get(ch_name, []),
                y=peaks_y.get(ch_name, []),
                mode='markers',
                name=f"{ch_name}",
                showlegend=False,
                hoverinfo='x+y+text+name',
                hovertemplate = "size (bp): %{x}<br>" + "height: %{y}<br>" if ch_name not in ('LIZ') else "",
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
        Creates an expander with individual channel sub-plots
        and optional peak tables.
        """
        smap = trace_dict.get('smap', [])
        channels_data = trace_dict.get('channels_peaks', {})
        peaks_x = trace_dict.get('detected_peaks_x', {})
        peaks_y = trace_dict.get('detected_peaks_y', {})
        colors = trace_dict.get('colors', {})

        with st.expander("Individual channels"):
            # The checkboxes from your snippet
            st.session_state.FULL_GENOTYPES = st.checkbox("Show all detected peaks", value=False)
            st.session_state.SHOW_INDIV_TABLES = st.checkbox("Show table with individual channels", value=False)

            for ch_name, ch_color in colors.items():
                x_vals = peaks_x.get(ch_name, [])
                y_vals = peaks_y.get(ch_name, [])
                # Filter out very low position if you prefer
                valid_data = [(x, y) for x, y in zip(x_vals, y_vals) if x > 50]

                if not valid_data and not channels_data.get(ch_name):
                    continue  # No data, skip channel

                fig = go.Figure()


                if "marker_list" in st.session_state:
                    for marker in st.session_state.marker_list:
                        bin_range = marker["bins"]
                        channel = marker["channel"]
                        if channel != ch_name: continue

                        if bin_range and channel in colors:
                            fig.add_vrect(
                                x0=bin_range[0], x1=bin_range[1],
                                fillcolor=colors[channel], opacity=0.1,
                                layer="below", line_width=0
                            )

                # Add the entire channel trace
                fig.add_trace(go.Scatter(
                    x=smap,
                    y=channels_data.get(ch_name, []),
                    mode='lines',
                    name=ch_name,
                    fill='tozeroy',
                    hoverinfo='skip',
                    marker=dict(color=ch_color)
                ))

                # Add the peak markers
                if valid_data:
                    dx = [p[0] for p in valid_data]
                    dy = [p[1] for p in valid_data]
                    fig.add_trace(go.Scatter(
                        x=dx,
                        y=dy,
                        mode='markers',
                        name=f"{ch_name}",
                        showlegend=False,
                        hovertemplate = "size (bp): %{x}<br>" + "height: %{y}<br>" if ch_name not in ('LIZ') else "",
                        marker=dict(color=ch_color, size=3)
                    ))
                # Scale
                max_y = max(y_vals) * 1.25 if y_vals else 100
                max_x = max(x_vals) * 1.05 if x_vals else 100

                fig.update_layout(
                    xaxis=dict(range=[0, max_x], title="Size (bp)"),
                    yaxis=dict(range=[0, max_y], title="Intensity"),
                    height=200,
                    margin=dict(l=0, r=0, t=0, b=0),
                )

                st.markdown(f"### {ch_name}")
                st.plotly_chart(fig, use_container_width=True)

                # If user wants tables
                if st.session_state.SHOW_INDIV_TABLES and valid_data:
                    table_rows = [{"size (bp)": x, "height": y} for x, y in valid_data]
                    avg_h = sum(row["height"] for row in table_rows) / len(table_rows)
                    # If not showing "all" peaks, filter above average
                    if not st.session_state.FULL_GENOTYPES:
                        table_rows = [row for row in table_rows if row["height"] > avg_h]

                    st.dataframe(table_rows, use_container_width=True)
                st.divider()

    def _update_plot_window(self) -> None:
        """
        Show the main trace and optional per-channel details
        for whichever trace the user selected in the sidebar.
        """
        st.header(f"Viewing: {st.session_state.SELECTED_TRACE}")
        trace_key = st.session_state.SELECTED_TRACE
        trace_dict = None
        # Find the matching dictionary
        for k, v in st.session_state.PROCESSED_FLA.items():
            if v['name'] == trace_key:
                trace_dict = v
                break
        if not trace_dict:
            st.error("Selected trace not found in session.")
            return

        # Plot the total trace
        fig = self._plot_total_trace(trace_dict)
        st.plotly_chart(fig, use_container_width=True)
        # Per-channel expansions
        self._plot_per_channel_traces(trace_dict)

streamlit_app = App()
streamlit_app._init_page()
