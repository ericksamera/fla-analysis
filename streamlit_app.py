#!/usr/bin/env python3
__description__ =\
"""
Purpose: Streamlit wrapper for fla-viewer-advanced using preprocessed JSON files.
"""
__author__ = "Erick Samera"
__version__ = "1.2.1"
__comments__ = "stable enough; altered zoom"
# --------------------------------------------------
import streamlit as st
import json
import plotly.graph_objects as go
import csv
from pathlib import Path

# --------------------------------------------------
class App:
    def __init__(self):
        """
        Initialize the Streamlit page configuration.
        """
        self.title = "FLA-viewer"
        self.emoji = ':rainbow:'
        st.set_page_config(
            page_title=f"abi-sauce | {self.title}",
            page_icon=self.emoji,
            layout='wide',
            initial_sidebar_state='expanded'
        )
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
        Instantiate the main page.
        """
        self._init_sidebar()
        # Use JSON uploader now instead of FSA uploader.
        if 'UPLOADED_JSON' not in st.session_state: 
            self._init_file_uploader()
        else:
            self._update_plot_window()
        return None
    def _init_sidebar(self) -> None:
        """
        Instantiate the sidebar.
        """
        with st.sidebar:
            st.title(f"{self.emoji} abi-sauce | {self.title}")
            st.markdown('This script processes pre-parsed FLA JSON files with extended raw data.')
            if 'UPLOADED_JSON' in st.session_state:
                st.divider()
                if 'FAILED_FILES' in st.session_state:
                    with st.expander('**⚠ WARNINGS:**', expanded=True):
                        for failed_file_name in st.session_state.FAILED_FILES:
                            st.error(f"{failed_file_name} failed to load!", icon="⚠")
                with st.expander('**TRACE FILES:**', expanded=True):
                    st.session_state.SORTED_LIST = list(st.session_state.PROCESSED_FLA.values())
                    st.session_state.SELECTED_TRACE = st.selectbox(
                        'Select trace file to view:', 
                        options=[trace_object['name'] for trace_object in st.session_state.SORTED_LIST]
                    )
                st.button('Reset & Upload New', type='primary', on_click=self._reset_state, use_container_width=True)
            st.divider()
            with st.expander('MORE INFO'):
                st.markdown(
                    'Fragment length analysis (FLA) is a technique for determining '
                    'the size of DNA fragments which is useful for genotyping microsatellites.'
                )
                st.markdown(
                    'Capillary electrophoresis with fluorescently tagged fragments provides '
                    'high resolution data that is processed with specialized software.'
                )
                st.caption(f'[@{__author__}](https://github.com/ericksamera)\t|\tv{__version__}\t|\t{__comments__}')
    def _init_file_uploader(self) -> None:
        """
        Initialize file upload handling for preprocessed JSON files.
        """
        st.header('Select preprocessed FLA JSON files')
        uploaded_files: list = st.file_uploader(
            'Upload JSON files produced by the FLA processing pipeline.',
            type=['json'],
            accept_multiple_files=True
        )
        with st.form("upload form", clear_on_submit=True):
            submit_button = st.form_submit_button("Load files", on_click=self._upload_files, args=(uploaded_files,))
            if submit_button and not uploaded_files:
                st.error('Select some JSON files!')
    def _upload_files(self, _st_uploaded_files: list) -> None:
        """
        Handle file upload and transform JSON content into the structure used for plotting.
        """
        if not _st_uploaded_files:
            return None
        st.session_state.UPLOADED_JSON = _st_uploaded_files
        st.session_state.PROCESSED_FLA = self._process_files(st.session_state.UPLOADED_JSON)
    def _process_files(self, _st_uploaded_files) -> dict:
        """
        Process uploaded JSON files and convert them into a dictionary with keys used for plotting.
        Expected JSON structure:
        {
            "file": <name>,
            "fsa_data": {"smap": [...], "channels": {channel: [intensity values]}},
            "detected_peaks": {channel: [{"position": float, "intensity": float}, ...]},
            "marker_results": {...}
        }
        This function builds additional keys:
         - "channels_peaks": same as fsa_data["channels"]
         - "detected_peaks_x": {channel: [position, ...]}
         - "detected_peaks_y": {channel: [intensity, ...]}
         - "colors": a static color mapping.
        """
        processed_files = {}
        # Define a static color mapping.
        color_mapping = {
            '6-FAM': "blue",
            'VIC': "green",
            'NED': "black",
            'PET': "red",
            'LIZ': "orange"
        }
        for file in sorted(_st_uploaded_files, key=lambda x: x.name):
            try:
                data = json.load(file)
            except Exception as e:
                if 'FAILED_FILES' not in st.session_state:
                    st.session_state.FAILED_FILES = []
                st.session_state.FAILED_FILES.append(file.name)
                continue
            # Build a new dictionary with expected keys.
            processed_dict = {}
            processed_dict['name'] = data.get("file", "Unknown")
            fsa_data = data.get("fsa_data", {})
            processed_dict['smap'] = fsa_data.get("smap", [])
            # Use "channels" from fsa_data as channels_peaks.
            processed_dict['channels_peaks'] = fsa_data.get("channels", {})
            # Transform detected_peaks into separate x and y lists.
            detected_peaks = data.get("detected_peaks", {})
            detected_peaks_x = {}
            detected_peaks_y = {}
            for channel, peaks in detected_peaks.items():
                detected_peaks_x[channel] = [peak.get("position") for peak in peaks]
                detected_peaks_y[channel] = [peak.get("intensity") for peak in peaks]
            processed_dict['detected_peaks_x'] = detected_peaks_x
            processed_dict['detected_peaks_y'] = detected_peaks_y
            processed_dict['colors'] = color_mapping
            processed_files[processed_dict['name']] = processed_dict
        return processed_files
    def _plot_total_trace(self, _trace_dict):
        """
        Plot the overall trace using the smap, channels, and detected peaks.
        """
        all_x_vals = []
        for dye_color in _trace_dict.get('channels_peaks', {}):
            all_x_vals += _trace_dict.get('smap', [])
        if not all_x_vals:
            max_x_val = 100
        else:
            max_x_val = max(all_x_vals)
        all_y_vals = []
        for dye_color in _trace_dict.get('channels_peaks', {}):
            all_y_vals += _trace_dict['channels_peaks'][dye_color]
        if not all_y_vals:
            max_y_val = 100
        else:
            max_y_val = max(all_y_vals) * 1.05

        max_heights_with_pad = {key: (max(_trace_dict['detected_peaks_y'][key]) * 1.05 if _trace_dict['detected_peaks_y'].get(key) else max_y_val)
                                for key in _trace_dict.get('channels_peaks', {})}

        fig = go.Figure()
        for dye_color, values in _trace_dict.get('channels_peaks', {}).items():
            fig.add_trace(
                go.Scatter(
                    mode='lines',
                    x=_trace_dict.get('smap', []),
                    y=values,
                    hoverinfo='skip',
                    marker=dict(size=0.8, color=_trace_dict['colors'].get(dye_color, "gray")),
                    name=dye_color,
                    fill='tozeroy',
                    legendgroup=dye_color
                )
            )
            fig.add_trace(
                go.Scatter(
                    mode='markers',
                    x=_trace_dict['detected_peaks_x'].get(dye_color, []),
                    y=_trace_dict['detected_peaks_y'].get(dye_color, []),
                    hoverinfo='x+y+text+name' if dye_color not in ('LIZ',) else 'skip',
                    hovertemplate="size (bp): %{x}<br>height: %{y}<br>",
                    marker=dict(size=0.8, color=_trace_dict['colors'].get(dye_color, "gray")),
                    name=dye_color,
                    showlegend=False,
                    legendgroup=dye_color
                )
            )
        fig.update_layout(
            legend_title_text='Channels',
            dragmode="zoom",
            xaxis=dict(range=(0, max_x_val+20)),
            yaxis=dict(range=(0, max_y_val))
        )
        buttons_list = [dict(args=[{"yaxis": dict(autorange=True)}], label="Scale to tallest", method="relayout")]
        for dye_color in _trace_dict.get('channels_peaks', {}):
            if dye_color not in max_heights_with_pad: continue
            buttons_list.append(dict(args=[{"yaxis": dict(autorange=False, range=(0, max_heights_with_pad[dye_color]))}],
                                       label=f"Scale to {dye_color}", method="relayout"))
        fig.update_layout(
            updatemenus=[dict(
                type="buttons",
                buttons=buttons_list,
                active=0,
                showactive=True,
            )]
        )
        return fig
    def _per_channel_bar_plot(self, _trace_dict: dict, _dye_name: str, _average_height: float):
        """
        Plot per-channel bar plot.
        """
        fig = go.Figure()
        fig.add_trace(go.Scattergl(
            mode='lines',
            x=_trace_dict.get('smap', []),
            y=_trace_dict['channels_peaks'].get(_dye_name, []),
            name=_dye_name,
            hoverinfo='skip',
            fill='tozeroy',
            marker=dict(color=_trace_dict['colors'].get(_dye_name, "gray"))
        ))
        fig.add_trace(go.Scattergl(
            mode='markers',
            x=_trace_dict['detected_peaks_x'].get(_dye_name, []),
            y=_trace_dict['detected_peaks_y'].get(_dye_name, []),
            hovertemplate="size (bp): %{x}<br>height: %{y}<br>",
            marker=dict(size=0.5, color=_trace_dict['colors'].get(_dye_name, "gray")),
            name=_dye_name,
            showlegend=False
        ))
        fig.update_layout(
            margin=dict(l=50, r=150, t=0, b=0),
            height=100,
            legend_title_text='Channels',
            xaxis=dict(range=(0, max(_trace_dict['detected_peaks_x'].get(_dye_name, [0]))*1.05)),
            yaxis=dict(range=(0, int(_average_height*1.25))),
            modebar=dict(remove=["zoom", "pan", "select", "lasso", "toImage"])
        )
        return fig
    def _plot_per_channel_traces(self, _trace_dict) -> None:
        """
        Plot individual channel traces and optionally display tables of predicted peaks.
        """
        with st.expander('Individual channels'):
            st.session_state.FULL_GENOTYPES = st.checkbox("Show all detected peaks", value=False)
            st.session_state.SHOW_INDIV_TABLES = st.checkbox("Show table with individual channel peaks", value=False)
            for dye_color in _trace_dict.get('colors', {}):
                averages = []
                full_genotype = []
                for x_val, y_val in zip(_trace_dict['detected_peaks_x'].get(dye_color, []),
                                          _trace_dict['detected_peaks_y'].get(dye_color, [])):
                    if x_val < 25:
                        continue
                    averages.append(y_val)
                    full_genotype.append({'size (bp)': x_val, 'height': y_val})
                if averages:
                    average_height = max(averages)
                    st.plotly_chart(
                        self._per_channel_bar_plot(_trace_dict, dye_color, average_height),
                        use_container_width=True
                    )
                if full_genotype:
                    avg_height = sum(geno['height'] for geno in full_genotype) / len(full_genotype)
                    filtered = [geno for geno in full_genotype if geno['height'] > avg_height]
                    if st.session_state.SHOW_INDIV_TABLES:
                        if dye_color not in ('LIZ',):
                            st.dataframe(full_genotype if st.session_state.FULL_GENOTYPES else filtered, use_container_width=True)
        return None
    def _update_plot_window(self) -> None:
        """
        Update the main plot window based on the selected JSON trace.
        """
        st.header(f"{st.session_state.SELECTED_TRACE}")
        trace_key = st.session_state.SELECTED_TRACE
        trace_dict = st.session_state.PROCESSED_FLA.get(trace_key)
        if not trace_dict:
            st.error("Selected trace not found.")
            return
        st.plotly_chart(
            self._plot_total_trace(trace_dict),
            use_container_width=True
        )
        self._plot_per_channel_traces(trace_dict)
        return None
    def _reset_state(self) -> None:
        """
        Reset session state.
        """
        for key in list(st.session_state.keys()):
            del st.session_state[key]

# --------------------------------------------------
if __name__=="__main__":
    streamlit_app = App()
    streamlit_app._init_page()