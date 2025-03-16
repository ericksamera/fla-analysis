#!/usr/bin/env python3
"""
Purpose: Streamlit wrapper that can either upload or process .fsa files using FLAProcessor,
         or upload preprocessed .json files, and then view them in an interactive trace viewer.

Author: Erick Samera
Version: 1.2.1
Comments: stable enough; altered zoom
"""
import streamlit as st
import json
import os
import plotly.graph_objects as go


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
        # st.set_page_config(
        #     page_title=f"abi-sauce | {self.title}",
        #     page_icon=self.emoji,
        #     layout='wide',
        #     initial_sidebar_state='expanded'
        # )
        # A bit of styling for the sidebar
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
        self._init_sidebar()

        if "PROCESSED_FLA" not in st.session_state or not st.session_state.PROCESSED_FLA:
            # If we have no data, let user upload FSA or JSON
            self._init_file_uploader()
        else:
            # Otherwise, let user view the data
            self._update_plot_window()

    def _init_sidebar(self) -> None:
        with st.sidebar:
            st.title(f"{self.emoji} abi-sauce | {self.title}")
            st.markdown('This script can process `.fsa` files (if FLAProcessor is available) or parse `.json` files, then display the resulting traces.')
            
            # If there's processed data in session, let user pick a trace
            if "PROCESSED_FLA" in st.session_state and st.session_state.PROCESSED_FLA:
                st.divider()
                with st.expander('**TRACE FILES:**', expanded=True):
                    # Sort the processed list by name
                    trace_list = list(st.session_state.PROCESSED_FLA.values())
                    # Let user pick which trace to view
                    trace_names = [t['name'] for t in trace_list]
                    st.session_state.SELECTED_TRACE = st.selectbox(
                        'Select trace file to view:',
                        options=trace_names
                    )

                # Button to reset session
                st.button('Reset & Upload New', type='primary', on_click=self._reset_state, use_container_width=True)

            st.divider()
            with st.expander('MORE INFO'):
                st.markdown(
                    'Fragment length analysis (FLA) determines the size of DNA fragments, '
                    'useful for genotyping microsatellites.\n\n'
                    'Capillary electrophoresis with fluorescently tagged fragments provides '
                    'high resolution data that is processed with specialized software.'
                )
                st.caption(f'[Erick Samera](https://github.com/ericksamera) | v1.2.1 | stable enough; altered zoom')

    def _init_file_uploader(self) -> None:
        """
        If we have no data in session, let user upload .fsa or .json files.
        """
        st.header("Upload .fsa or preprocessed .json files to view")
        uploaded_files = st.file_uploader(
            "Select files",
            type=["fsa", "json"],
            accept_multiple_files=True
        )
        if uploaded_files:
            if st.button("Process/Load files"):
                self._upload_files(uploaded_files)

    def _upload_files(self, _st_uploaded_files):
        """
        Processes .fsa files with FLAProcessor (if available), or loads JSON into session.
        """
        if "PROCESSED_FLA" not in st.session_state:
            st.session_state.PROCESSED_FLA = {}

        # If you want a default config or marker
        if HAVE_PROCESSOR:
            fla_config = FLAConfig()
            marker_list = [MarkerConfig(marker="DefaultMarker", channel="6-FAM", repeat_unit=3, bins=(150, 250))]
        else:
            fla_config, marker_list = None, []

        for up_file in _st_uploaded_files:
            file_name = up_file.name
            file_ext = os.path.splitext(file_name)[1].lower()

            if file_ext == ".fsa":
                if not HAVE_PROCESSOR:
                    st.error("Cannot process .fsa because FLAProcessor isn't available.")
                    continue
                st.write(f"**Processing FSA**: {file_name}")
                path = os.path.join("/tmp", file_name)
                with open(path, "wb") as f:
                    f.write(up_file.getbuffer())

                processor = FLAProcessor(path, config=fla_config, marker_config=marker_list)
                success = processor.process_all()
                if success:
                    # Build JSON data
                    json_out = {
                        "file": processor.fsa_data["name"],
                        "fsa_data": {
                            "smap": processor.fsa_data["smap"],
                            "channels": {
                                ch: list(map(float, arr))
                                for ch, arr in processor.fsa_data["channels"].items()
                            }
                        },
                        "detected_peaks": {
                            ch: [
                                {
                                    "position": p.position,
                                    "intensity": p.intensity,
                                    "saturated": p.saturated,
                                    "note": p.note
                                }
                                for p in peaks
                            ]
                            for ch, peaks in processor.detected_peaks.items()
                        },
                        "marker_results": processor.marker_results
                    }
                    # Convert to "viewer" format
                    st.session_state.PROCESSED_FLA[processor.fsa_data["name"]] = self._convert_for_plotting(json_out)
                else:
                    st.error(f"Could not process FSA: {file_name}")

            elif file_ext == ".json":
                st.write(f"**Parsing JSON**: {file_name}")
                try:
                    data = json.load(up_file)
                    # Convert to "viewer" format
                    key_name = data.get("file", file_name)
                    st.session_state.PROCESSED_FLA[key_name] = self._convert_for_plotting(data)
                except Exception as exc:
                    st.error(f"JSON parse error for {file_name}: {exc}")
            else:
                st.warning(f"Skipping unknown file type: {file_name}")

        # Once done, if everything worked, user can pick from the sidebar

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
        max_y *= 1.05

        max_x = max(smap) if smap else 100
        max_x += 20

        for ch_name, intensity_list in channels_data.items():
            fig.add_trace(go.Scatter(
                x=smap,
                y=intensity_list,
                mode='lines',
                name=ch_name,
                fill='tozeroy',
                marker=dict(color=colors.get(ch_name, "gray"))
            ))
            # Now add any peaks
            fig.add_trace(go.Scatter(
                x=peaks_x.get(ch_name, []),
                y=peaks_y.get(ch_name, []),
                mode='markers',
                name=f"{ch_name} peaks",
                showlegend=False,
                marker=dict(color=colors.get(ch_name, "gray"))
            ))

        fig.update_layout(
            xaxis=dict(range=[0, max_x], title="Size (bp)"),
            yaxis=dict(range=[0, max_y], title="Intensity"),
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
                valid_data = [(x, y) for x, y in zip(x_vals, y_vals) if x > 25]

                if not valid_data and not channels_data.get(ch_name):
                    continue  # No data, skip channel

                fig = go.Figure()
                # Add the entire channel trace
                fig.add_trace(go.Scatter(
                    x=smap,
                    y=channels_data.get(ch_name, []),
                    mode='lines',
                    name=ch_name,
                    fill='tozeroy',
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
                        name=f"{ch_name} peaks",
                        showlegend=False,
                        marker=dict(color=ch_color, size=3)
                    ))
                # Scale
                max_y = max(y_vals) * 1.25 if y_vals else 100
                max_x = max(x_vals) * 1.05 if x_vals else 100

                fig.update_layout(
                    xaxis=dict(range=[0, max_x], title="Size (bp)"),
                    yaxis=dict(range=[0, max_y], title="Intensity"),
                    height=200,
                    margin=dict(l=50, r=50, t=30, b=30)
                )

                st.subheader(f"Channel: {ch_name}")
                st.plotly_chart(fig, use_container_width=True)

                # If user wants tables
                if st.session_state.SHOW_INDIV_TABLES and valid_data:
                    table_rows = [{"size (bp)": x, "height": y} for x, y in valid_data]
                    avg_h = sum(row["height"] for row in table_rows) / len(table_rows)
                    # If not showing "all" peaks, filter above average
                    if not st.session_state.FULL_GENOTYPES:
                        table_rows = [row for row in table_rows if row["height"] > avg_h]

                    st.dataframe(table_rows, use_container_width=True)

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

    def _reset_state(self) -> None:
        """
        Clears everything in st.session_state so user can start fresh.
        """
        for key in list(st.session_state.keys()):
            del st.session_state[key]

streamlit_app = App()
streamlit_app._init_page()
