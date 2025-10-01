from __future__ import annotations

from typing import List, Optional

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go


def histogram(df: pd.DataFrame, column: str):
    fig = px.histogram(df, x=column, nbins=30, title=f"Histograma: {column}")
    fig.update_layout(margin=dict(l=20, r=20, t=40, b=20))
    return fig


def scatter(df: pd.DataFrame, x: str, y: str, color: Optional[str] = None):
    fig = px.scatter(df, x=x, y=y, color=color, title=f"Dispersión: {x} vs {y}")
    fig.update_layout(margin=dict(l=20, r=20, t=40, b=20))
    return fig


def corr_heatmap(corr: pd.DataFrame):
    fig = go.Figure(
        data=go.Heatmap(
            z=corr.values,
            x=corr.columns,
            y=corr.index,
            colorscale="RdBu",
            zmin=-1,
            zmax=1,
            colorbar=dict(title="r")
        )
    )
    fig.update_layout(title="Matriz de correlación", xaxis_nticks=36, margin=dict(l=60, r=20, t=40, b=60))
    return fig

