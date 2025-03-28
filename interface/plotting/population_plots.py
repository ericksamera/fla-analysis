# interface/plotting/population_plots.py

import pandas as pd
import plotly.express as px
from typing import Optional

def plot_pcoa(
    coords: pd.DataFrame,
    metadata: dict,
    color_by: Optional[str] = None,
    title: str = "PCoA of Genetic Distances"
):
    df = coords.copy()
    df["SampleID"] = df.index

    # Get display name from metadata
    df["Label"] = df["SampleID"].apply(
        lambda sid: metadata.get(sid, {}).get("Sample Name", sid)
    )

    # Add coloring by metadata field
    if color_by:
        df[color_by] = df["SampleID"].apply(
            lambda sid: metadata.get(sid, {}).get(color_by, "—")
        )
        color_arg = color_by
    else:
        color_arg = None

    fig = px.scatter(
        df,
        x="PCoA1",
        y="PCoA2",
        color=color_arg,
        text="Label",
        hover_data=["SampleID", "Label"] + ([color_by] if color_by else []),
        title=title,
        labels={"PCoA1": "PCoA 1", "PCoA2": "PCoA 2"},
    )

    fig.update_traces(
        textposition="top center",
        marker=dict(size=10, line=dict(width=0.5, color="DarkSlateGrey"))
    )

    fig.update_layout(
        height=550,
        margin=dict(t=40, b=40, l=40, r=40),
        hovermode="closest"
    )

    return fig
