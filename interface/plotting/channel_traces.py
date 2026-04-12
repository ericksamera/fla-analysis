# interface/plotting/channel_traces.py

import numpy as np

import plotly.graph_objects as go
from fla_pipeline.models.sample import Sample
from interface.backend.trace_utils import find_suppressed_regions

COLOR_MAP = {
    "6-FAM": "blue",
    "VIC": "green",
    "NED": "black",
    "PET": "red",
    "LIZ": "orange",
}

GENESCAN_1200_LIZ_SIZES = [
    20, 30, 40, 60, 80, 100, 114, 120, 140, 160, 180, 200, 214, 220, 240, 250,
    260, 280, 300, 314, 320, 340, 360, 380, 400, 414, 420, 440, 460, 480, 500,
    514, 520, 540, 560, 580, 600, 614, 620, 640, 660, 680, 700, 714, 720, 740,
    760, 780, 800, 820, 840, 850, 860, 880, 900, 920, 940, 960, 980, 1000, 1020,
    1040, 1060, 1080, 1100, 1120, 1160, 1200,
]
GENESCAN_1200_LIZ_MAJOR_LABELS = [20] + list(range(100, 1300, 100))


def _match_expected_liz_size(position, expected_sizes, tolerance=2.0):
    nearest = min(expected_sizes, key=lambda size: abs(size - position))
    delta = round(position - nearest, 2)
    return nearest, delta, abs(delta) <= tolerance


def _position_key(value):
    return round(float(value), 2)


def _get_prominent_liz_peaks(peaks, min_relative_intensity=0.05):
    if not peaks:
        return []

    max_intensity = max(getattr(p, "intensity", 0) for p in peaks)
    min_intensity = max(100, max_intensity * min_relative_intensity)
    prominent = [p for p in peaks if getattr(p, "intensity", 0) >= min_intensity]
    prominent = sorted(prominent, key=lambda peak: peak.position)
    if len(prominent) > len(GENESCAN_1200_LIZ_SIZES):
        prominent = sorted(
            sorted(prominent, key=lambda peak: peak.intensity, reverse=True)[
                : len(GENESCAN_1200_LIZ_SIZES)
            ],
            key=lambda peak: peak.position,
        )
    return prominent


def _detect_liz_standard(peaks, tolerance=2.0):
    prominent = _get_prominent_liz_peaks(peaks)
    if not prominent:
        return None

    observed = [p.position for p in prominent]
    matched = [
        size
        for size in GENESCAN_1200_LIZ_SIZES
        if any(abs(size - pos) <= tolerance for pos in observed)
    ]

    high_range_matches = [size for size in matched if size >= 600]
    if max(observed, default=0) < 1000:
        return None
    if len(matched) < 12 or len(high_range_matches) < 5:
        return None

    return {
        "name": "GeneScan 1200 LIZ",
        "expected_sizes": GENESCAN_1200_LIZ_SIZES,
        "matched_sizes": matched,
        "major_labels": [size for size in GENESCAN_1200_LIZ_MAJOR_LABELS if size <= 1200],
    }


def _get_nearest_expected_index(expected_logs, target_log, lower_idx, upper_idx):
    candidates = []
    insertion_idx = int(np.searchsorted(expected_logs, target_log))
    for idx in (insertion_idx - 1, insertion_idx):
        if lower_idx <= idx <= upper_idx:
            candidates.append(idx)

    if not candidates:
        return lower_idx

    return min(candidates, key=lambda idx: abs(expected_logs[idx] - target_log))


def _fit_scan_index_liz_standard(peaks):
    prominent = _get_prominent_liz_peaks(peaks)
    if len(prominent) < 20:
        return None

    observed_positions = np.array([p.position for p in prominent], dtype=float)
    if observed_positions[-1] - observed_positions[0] < 1000:
        return None

    expected_sizes = np.array(GENESCAN_1200_LIZ_SIZES, dtype=float)
    assigned_idx = np.round(
        np.linspace(0, len(expected_sizes) - 1, len(observed_positions))
    ).astype(int)
    assigned_sizes = expected_sizes[assigned_idx]

    expected_logs = np.log(expected_sizes)
    coeff = np.polyfit(np.log(assigned_sizes), observed_positions, 1)
    slope, intercept = coeff
    if slope <= 0:
        return None

    predicted_logs = (observed_positions - intercept) / slope
    refined_idx = []
    prev_idx = -1
    total_expected = len(expected_sizes)
    total_observed = len(observed_positions)
    for obs_i, predicted_log in enumerate(predicted_logs):
        lower_idx = prev_idx + 1
        upper_idx = total_expected - (total_observed - obs_i)
        idx = _get_nearest_expected_index(
            expected_logs, predicted_log, lower_idx, upper_idx
        )
        refined_idx.append(idx)
        prev_idx = idx

    assigned_sizes = expected_sizes[refined_idx]
    coeff = np.polyfit(np.log(assigned_sizes), observed_positions, 1)
    slope, intercept = coeff
    fitted_positions = slope * np.log(assigned_sizes) + intercept
    ss_res = float(np.sum((observed_positions - fitted_positions) ** 2))
    ss_tot = float(np.sum((observed_positions - observed_positions.mean()) ** 2))
    fit_r2 = 1.0 - (ss_res / ss_tot) if ss_tot else 1.0
    if slope <= 0 or fit_r2 < 0.95:
        return None

    assignment_by_position = {
        _position_key(peak.position): int(size)
        for peak, size in zip(prominent, assigned_sizes)
    }
    major_label_positions = {
        size: float(slope * np.log(size) + intercept)
        for size in GENESCAN_1200_LIZ_MAJOR_LABELS
        if size <= 1200
    }

    return {
        "name": "GeneScan 1200 LIZ",
        "expected_sizes": GENESCAN_1200_LIZ_SIZES,
        "assignment_by_position": assignment_by_position,
        "major_label_positions": major_label_positions,
        "fit_r2": fit_r2,
        "peak_count": len(prominent),
    }


def _add_liz_standard_overlay(fig, trace, standard_info):
    if not standard_info or len(trace) == 0:
        return

    y_max = float(np.max(trace))
    if y_max <= 0:
        return

    for size in standard_info["expected_sizes"]:
        fig.add_vline(
            x=size,
            line_width=1,
            line_dash="dot",
            line_color="rgba(255, 165, 0, 0.18)",
            layer="below",
        )

    for size in standard_info["major_labels"]:
        fig.add_annotation(
            x=size,
            y=y_max * 1.04,
            text=str(size),
            showarrow=False,
            yanchor="bottom",
            font=dict(size=9, color="orange"),
        )

    fig.add_annotation(
        x=0.99,
        y=0.98,
        xref="paper",
        yref="paper",
        text="Detected size standard: GeneScan 1200 LIZ",
        showarrow=False,
        xanchor="right",
        yanchor="top",
        bordercolor="orange",
        borderwidth=1,
        bgcolor="rgba(255,255,255,0.7)",
        font=dict(size=11),
    )

    fig.update_yaxes(range=[0, y_max * 1.12])


def _get_trace_x(trace, smap):
    if len(smap) == len(trace) and len(smap) > 0:
        return smap, "Size (bp)", True
    return np.arange(len(trace), dtype=float), "Scan index", False


def _peak_hover_title(has_size_map):
    return "Size" if has_size_map else "Scan index"


def _peak_hover_suffix(has_size_map):
    return " bp" if has_size_map else ""


def _get_scan_index_liz_note(peak, standard_info):
    assigned_size = standard_info["assignment_by_position"].get(
        _position_key(peak.position)
    )
    if assigned_size is None:
        return "off-ladder"
    return f"putative {assigned_size} bp"


def _add_scan_index_liz_overlay(fig, trace, standard_info):
    if not standard_info or len(trace) == 0:
        return

    y_max = float(np.max(trace))
    if y_max <= 0:
        return

    trace_end = max(len(trace) - 1, 0)
    for size, x_pos in standard_info["major_label_positions"].items():
        if not 0 <= x_pos <= trace_end:
            continue
        fig.add_vline(
            x=x_pos,
            line_width=1,
            line_dash="dot",
            line_color="rgba(255, 165, 0, 0.18)",
            layer="below",
        )
        fig.add_annotation(
            x=x_pos,
            y=y_max * 1.04,
            text=str(size),
            showarrow=False,
            yanchor="bottom",
            font=dict(size=9, color="orange"),
        )

    fig.add_annotation(
        x=0.99,
        y=0.98,
        xref="paper",
        yref="paper",
        text=(
            "Detected size standard: GeneScan 1200 LIZ "
            f"(scan-index fit, R²={standard_info['fit_r2']:.3f})"
        ),
        showarrow=False,
        xanchor="right",
        yanchor="top",
        bordercolor="orange",
        borderwidth=1,
        bgcolor="rgba(255,255,255,0.7)",
        font=dict(size=11),
    )

    fig.update_yaxes(range=[0, y_max * 1.12])


def make_total_trace_figure(
    sample: Sample, marker_list: list, config: dict
) -> go.Figure:
    smap = sample.fsa_data["smap"]
    channels = sample.fsa_data["channels"]
    peaks = sample.peaks or {}

    fig = go.Figure()
    first_trace = next(iter(channels.values()), np.array([]))
    _, xaxis_title, has_size_map = _get_trace_x(first_trace, smap)
    liz_scan_standard = (
        _fit_scan_index_liz_standard(peaks.get("LIZ", [])) if not has_size_map else None
    )

    # Bin and tolerance overlays
    if marker_list and has_size_map:
        y_max = _get_max_visible_peak_in_regions(peaks, marker_list, config)
        if y_max > 0:
            fig.update_yaxes(range=[0, y_max * 1.1])
        for marker in marker_list:
            channel = marker.channel
            repeat = marker.repeat_unit or 1
            bmin, bmax = marker.bins
            tol = config.get("bin_tolerance", 2) * repeat
            color = COLOR_MAP.get(channel, "gray")

            fig.add_vrect(
                x0=bmin,
                x1=bmax,
                fillcolor=color,
                opacity=0.05,
                layer="below",
                line_width=0,
            )
            fig.add_vrect(
                x0=bmin - tol,
                x1=bmax + tol,
                fillcolor=color,
                opacity=0.05,
                layer="below",
                line_width=0,
            )

    # Channel traces and peaks
    for ch, trace in channels.items():
        color = COLOR_MAP.get(ch, "gray")
        x, _, _ = _get_trace_x(trace, smap)
        fig.add_trace(
            go.Scatter(
                x=x,
                y=trace,
                mode="lines",
                name=f"{ch} trace",
                line=dict(color=color),
                hoverinfo="skip",
            )
        )

        show_detected_peaks = ch in peaks and (
            (ch != "LIZ" and has_size_map) or (ch == "LIZ" and not has_size_map)
        )
        if show_detected_peaks:
            peak_label = _peak_hover_title(has_size_map)
            peak_suffix = _peak_hover_suffix(has_size_map)
            peak_kwargs = {
                "x": [p.position for p in peaks[ch]],
                "y": [p.intensity for p in peaks[ch]],
                "mode": "markers",
                "name": f"{ch}",
                "marker": dict(color=color, size=6),
                "showlegend": False,
            }
            if ch == "LIZ" and not has_size_map and liz_scan_standard:
                peak_kwargs["hoverinfo"] = "text"
                peak_kwargs["text"] = [
                    f"Scan index: {p.position}<br>"
                    f"Height: {p.intensity}<br>"
                    f"1200 LIZ: {_get_scan_index_liz_note(p, liz_scan_standard)}"
                    for p in peaks[ch]
                ]
            else:
                peak_kwargs["hovertemplate"] = (
                    f"{peak_label}: %{{x}}{peak_suffix}<br>Height: %{{y}}"
                )
            fig.add_trace(go.Scatter(**peak_kwargs))

    if liz_scan_standard and "LIZ" in channels:
        _add_scan_index_liz_overlay(fig, channels["LIZ"], liz_scan_standard)

    fig.update_layout(
        title="Electropherogram Trace",
        xaxis_title=xaxis_title,
        yaxis_title="Intensity",
        margin=dict(t=30, b=30),
        height=400,
        legend_title="Channels",
        dragmode="zoom",
    )
    return fig


def make_per_channel_figure(
    ch: str, sample: Sample, marker_list: list, config: dict
) -> go.Figure:
    smap = sample.fsa_data["smap"]
    trace = sample.fsa_data["channels"][ch]
    peaks = (sample.peaks or {}).get(ch, [])
    min_height = config.get("min_peak_height", 0)
    suppressed = sample.suppressed_peaks.get(ch, [])

    fig = go.Figure()
    x, xaxis_title, has_size_map = _get_trace_x(trace, smap)
    gray_regions = (
        find_suppressed_regions(smap, trace, suppressed) if has_size_map else []
    )

    # Trace segments
    last_idx = 0
    for i, (left, right) in enumerate(gray_regions):
        if left > last_idx:
            fig.add_trace(
                go.Scatter(
                    x=x[last_idx:left],
                    y=trace[last_idx:left],
                    mode="lines",
                    line=dict(color=COLOR_MAP.get(ch, "gray")),
                    hoverinfo="skip",
                    showlegend=(i == 0),
                    name=f"{ch}",
                )
            )
        fig.add_trace(
            go.Scatter(
                x=x[left:right],
                y=trace[left:right],
                mode="lines",
                line=dict(color="gray", dash="dot"),
                hoverinfo="skip",
                showlegend=False,
            )
        )
        last_idx = right

    if last_idx < len(trace):
        fig.add_trace(
            go.Scatter(
                x=x[last_idx:],
                y=trace[last_idx:],
                mode="lines",
                line=dict(color=COLOR_MAP.get(ch, "gray")),
                hoverinfo="skip",
                showlegend=False,
            )
        )

    if marker_list and ch != "LIZ" and has_size_map:
        y_max = _get_max_visible_peak_in_regions(
            sample.peaks or {}, marker_list, config, target_channel=ch
        )
        if y_max > 0:
            fig.update_yaxes(range=[0, y_max * 1.1])
        for marker in marker_list:
            if marker.channel != ch:
                continue
            repeat = marker.repeat_unit
            tol = config.get("bin_tolerance", 2) * repeat
            bmin, bmax = marker.bins
            color = COLOR_MAP.get(ch, "gray")
            fig.add_vrect(
                x0=bmin,
                x1=bmax,
                fillcolor=color,
                opacity=0.05,
                layer="below",
                line_width=0,
            )
            fig.add_vrect(
                x0=bmin - tol,
                x1=bmax + tol,
                fillcolor=color,
                opacity=0.05,
                layer="below",
                line_width=0,
            )

    visible_peaks = (
        peaks if ch == "LIZ" else [p for p in peaks if p.intensity >= min_height]
    )
    liz_standard = (
        _detect_liz_standard(visible_peaks) if ch == "LIZ" and has_size_map else None
    )
    scan_index_liz_standard = (
        _fit_scan_index_liz_standard(visible_peaks)
        if ch == "LIZ" and not has_size_map
        else None
    )
    if visible_peaks:
        annotations = []
        peak_label = _peak_hover_title(has_size_map)
        peak_suffix = _peak_hover_suffix(has_size_map)
        for p in visible_peaks:
            if liz_standard:
                expected_size, delta, is_match = _match_expected_liz_size(
                    p.position, liz_standard["expected_sizes"]
                )
                expected_note = (
                    f"{expected_size} bp (Δ={delta:+.2f})"
                    if is_match
                    else f"off-ladder; nearest {expected_size} bp (Δ={delta:+.2f})"
                )
                annotations.append(
                    f"{peak_label}: {p.position}{peak_suffix}<br>"
                    f"Height: {p.intensity}<br>"
                    f"1200 LIZ: {expected_note}"
                )
                continue

            if scan_index_liz_standard:
                annotations.append(
                    f"{peak_label}: {p.position}{peak_suffix}<br>"
                    f"Height: {p.intensity}<br>"
                    f"1200 LIZ: {_get_scan_index_liz_note(p, scan_index_liz_standard)}"
                )
                continue

            matched = []
            if has_size_map:
                for m in marker_list:
                    if m.channel != ch:
                        continue
                    bmin, bmax = m.bins
                    repeat = m.repeat_unit
                    tol = config.get("bin_tolerance", 2) * repeat
                    if (bmin - tol) <= p.position <= (bmax + tol):
                        matched.append(m.marker)
            marker_note = ", ".join(matched) if matched else "—"
            if has_size_map:
                annotations.append(
                    f"{peak_label}: {p.position}{peak_suffix}<br>"
                    f"Height: {p.intensity}<br>"
                    f"Marker: {marker_note}"
                )
            else:
                annotations.append(
                    f"{peak_label}: {p.position}{peak_suffix}<br>"
                    f"Height: {p.intensity}"
                )

        fig.add_trace(
            go.Scatter(
                x=[p.position for p in visible_peaks],
                y=[p.intensity for p in visible_peaks],
                mode="markers",
                name=f"({ch})",
                marker=dict(color=COLOR_MAP.get(ch, "gray"), size=6),
                showlegend=False,
                hoverinfo="text",
                text=annotations,
            )
        )

    # Determine x-axis range
    marker_bins = (
        [(m.bins, m.repeat_unit) for m in marker_list if m.channel == ch]
        if has_size_map
        else []
    )
    if marker_bins:
        min_start = min(bmin for (bmin, _), _ in marker_bins)
        max_end = max(bmax for (_, bmax), _ in marker_bins)
        max_repeat = max(rpt for _, rpt in marker_bins)
        pad = 10 * max_repeat
        x_range = [max(0, min_start - pad), max_end + pad]
    elif peaks and has_size_map:
        peak_positions = [p.position for p in peaks]
        pad = 10
        x_range = [max(0, min(peak_positions) - pad), max(peak_positions) + pad]
    elif len(x):
        x_range = [float(x[0]), float(x[-1])]
    else:
        x_range = [0.0, 1.0]

    if liz_standard:
        _add_liz_standard_overlay(fig, trace, liz_standard)
    elif scan_index_liz_standard:
        _add_scan_index_liz_overlay(fig, trace, scan_index_liz_standard)

    fig.update_layout(
        xaxis_title=xaxis_title,
        yaxis_title="Intensity",
        height=200,
        margin=dict(t=5, b=5),
        xaxis_range=x_range,
    )
    return fig


def _get_max_visible_peak_in_regions(
    peaks_by_channel, marker_list, config, target_channel=None
):
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
            p
            for p in peaks_by_channel[ch]
            if not getattr(p, "suppressed", False)
            and (bmin - tol) <= p.position <= (bmax + tol)
        ]

        local_max = max((p.intensity for p in visible_peaks), default=0)
        if local_max > max_intensity:
            max_intensity = local_max

    return max_intensity
