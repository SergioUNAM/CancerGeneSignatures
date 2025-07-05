import plotly.graph_objects as go
import numpy as np
import pandas as pd

SETTINGS = {
    "color_thresholds": [1, 2],
    "colors": ["#9fccff", "#FFFBCA", "#FFCCE1"],
}


def _color_scale(series: pd.Series) -> list:
    return np.select(
        [series > SETTINGS["color_thresholds"][1], series > SETTINGS["color_thresholds"][0]],
        SETTINGS["colors"][::-1][:2],
        SETTINGS["colors"][0],
    )


def tabla_fold_change(df: pd.DataFrame) -> go.Figure:
    """Genera una tabla interactiva con los resultados de Fold Change."""
    df_display = df.copy()
    df_display["fold_change"] = df_display["fold_change"].round(2)
    colores = _color_scale(df_display["fold_change"])
    fig = go.Figure(
        data=[go.Table(
            header=dict(values=list(df_display.columns)),
            cells=dict(values=[df_display[c] for c in df_display.columns], fill_color=[colores]),
        )]
    )
    return fig
