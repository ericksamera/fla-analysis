# interface/plotting/channel_traces.py

import plotly.graph_objects as go
from fla_pipeline.models.sample import Sample
from interface.backend.trace_utils import find_suppressed_regions

COLOR_MAP = {
    "6-FAM": "blue",
    "VIC": "green",
    "NED": "black",
    "PET": "red",
    "LIZ": "orange"
}

def make_total_trace_figure(sample: Sample, marker_list: list, config: dict) -> go.Figure:
    smap = sample.fsa_data["smap"]
    channels = sample.fsa_data["channels"]
    peaks = sample.peaks

    fig = go.Figure()

    # Bin and tolerance overlays
    if marker_list:
        y_max = _get_max_visible_peak_in_regions(peaks, marker_list, config)
        if y_max > 0:
            fig.update_yaxes(range=[0, y_max * 1.1])
        for marker in marker_list:
            channel = marker.channel
            repeat = marker.repeat_unit or 1
            bmin, bmax = marker.bins
            tol = config.get("bin_tolerance", 2) * repeat
            color = COLOR_MAP.get(channel, "gray")

            fig.add_vrect(x0=bmin, x1=bmax, fillcolor=color, opacity=0.05, layer="below", line_width=0)
            fig.add_vrect(x0=bmin - tol, x1=bmax + tol, fillcolor=color, opacity=0.05, layer="below", line_width=0)

    # Channel traces and peaks
    for ch, trace in channels.items():
        color = COLOR_MAP.get(ch, "gray")
        fig.add_trace(go.Scatter(x=smap, y=trace, mode="lines", name=f"{ch} trace", line=dict(color=color), hoverinfo="skip"))

        if ch != "LIZ" and ch in peaks:
            fig.add_trace(go.Scatter(
                x=[p.position for p in peaks[ch]],
                y=[p.intensity for p in peaks[ch]],
                mode="markers",
                name=f"{ch}",
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
    return fig

def make_per_channel_figure(ch: str, sample: Sample, marker_list: list, config: dict) -> go.Figure:
    smap = sample.fsa_data["smap"]
    trace = sample.fsa_data["channels"][ch]
    peaks = sample.peaks.get(ch, [])
    min_height = config.get("min_peak_height", 0)
    suppressed = sample.suppressed_peaks.get(ch, [])

    fig = go.Figure()
    gray_regions = find_suppressed_regions(smap, trace, suppressed)

    # Trace segments
    last_idx = 0
    for i, (left, right) in enumerate(gray_regions):
        if left > last_idx:
            fig.add_trace(go.Scatter(
                x=smap[last_idx:left], y=trace[last_idx:left],
                mode="lines", line=dict(color=COLOR_MAP.get(ch, "gray")),
                hoverinfo="skip", showlegend=(i == 0), name=f"{ch}"
            ))
        fig.add_trace(go.Scatter(
            x=smap[left:right], y=trace[left:right],
            mode="lines", line=dict(color="gray", dash="dot"), hoverinfo="skip",
            showlegend=False
        ))
        last_idx = right

    if last_idx < len(smap):
        fig.add_trace(go.Scatter(
            x=smap[last_idx:], y=trace[last_idx:],
            mode="lines", line=dict(color=COLOR_MAP.get(ch, "gray")), hoverinfo="skip",
            showlegend=False
        ))

    if marker_list and ch != "LIZ":
        y_max = _get_max_visible_peak_in_regions(sample.peaks, marker_list, config, target_channel=ch)
        if y_max > 0:
            fig.update_yaxes(range=[0, y_max * 1.1])
        for marker in marker_list:
            if marker.channel != ch:
                continue
            repeat = marker.repeat_unit
            tol = config.get("bin_tolerance", 2) * repeat
            bmin, bmax = marker.bins
            color = COLOR_MAP.get(ch, "gray")
            fig.add_vrect(x0=bmin, x1=bmax, fillcolor=color, opacity=0.05, layer="below", line_width=0)
            fig.add_vrect(x0=bmin - tol, x1=bmax + tol, fillcolor=color, opacity=0.05, layer="below", line_width=0)

    visible_peaks = peaks if ch == "LIZ" else [p for p in peaks if p.intensity >= min_height]
    if visible_peaks:
        annotations = []
        for p in visible_peaks:
            matched = []
            for m in marker_list:
                if m.channel != ch:
                    continue
                bmin, bmax = m.bins
                repeat = m.repeat_unit
                tol = config.get("bin_tolerance", 2) * repeat
                if (bmin - tol) <= p.position <= (bmax + tol):
                    matched.append(m.marker)
            marker_note = ", ".join(matched) if matched else "—"
            annotations.append(f"Size: {p.position} bp<br>Height: {p.intensity}<br>Marker: {marker_note}")

        fig.add_trace(go.Scatter(
            x=[p.position for p in visible_peaks],
            y=[p.intensity for p in visible_peaks],
            mode="markers", name=f"({ch})",
            marker=dict(color=COLOR_MAP.get(ch, "gray"), size=6),
            showlegend=False, hoverinfo="text", text=annotations
        ))

    # Determine x-axis range
    marker_bins = [
        (m.bins, m.repeat_unit)
        for m in marker_list if m.channel == ch
    ]
    if marker_bins:
        min_start = min(bmin for (bmin, _), _ in marker_bins)
        max_end = max(bmax for (_, bmax), _ in marker_bins)
        max_repeat = max(rpt for _, rpt in marker_bins)
        pad = 10 * max_repeat
        x_range = [max(0, min_start - pad), max_end + pad]
    elif peaks:
        peak_positions = [p.position for p in peaks]
        pad = 10
        x_range = [max(0, min(peak_positions) - pad), max(peak_positions) + pad]
    else:
        x_range = [min(smap), max(smap)]

    fig.update_layout(
        xaxis_title="Size (bp)",
        yaxis_title="Intensity",
        height=200,
        margin=dict(t=5, b=5),
        xaxis_range=x_range
    )
    return fig

def _get_max_visible_peak_in_regions(peaks_by_channel, marker_list, config, target_channel=None):
    tol_lookup = {
        m.marker: config.get("bin_tolerance", 2) * (m.repeat_unit or 1)
        for m in marker_list
    }

    max_intensity = 0
    for marker in marker_list:
        ch = marker.channel
        if target_channel and ch != target_channel:
            continue
        if ch not in peaks_by_channel:
            continue

        tol = tol_lookup[marker.marker]
        bmin, bmax = marker.bins

        visible_peaks = [
            p for p in peaks_by_channel[ch]
            if not getattr(p, "suppressed", False) and (bmin - tol) <= p.position <= (bmax + tol)
        ]

        local_max = max((p.intensity for p in visible_peaks), default=0)
        if local_max > max_intensity:
            max_intensity = local_max

    return max_intensity