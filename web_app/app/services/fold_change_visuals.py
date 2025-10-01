from __future__ import annotations

from typing import Iterable, Sequence

import numpy as np
import pandas as pd
import plotly.graph_objects as go

_TABLE_SETTINGS = {
    "thresholds": (1.0, 2.0),
    "colors": ("#9fccff", "#FFFBCA", "#FFCCE1"),
    "decimals": 2,
    "header_fill": "#2C3E50",
    "header_font_color": "white",
    "header_font_size": 14,
    "cell_font_size": 12,
    "char_width": 7,
    "min_column_width": 60,
    "max_column_width": 220,
}


def _ensure_columns(df: pd.DataFrame, required: Sequence[str]) -> None:
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Columnas requeridas ausentes: {', '.join(missing)}")


def _calculate_column_widths(df: pd.DataFrame) -> Iterable[int]:
    widths = []
    for column in df.columns:
        header_len = len(str(column))
        values = df[column].astype(str)
        data_len = int(values.map(len).max()) if not values.empty else 0
        width = max(header_len, data_len) * _TABLE_SETTINGS["char_width"]
        width = int(np.clip(width, _TABLE_SETTINGS["min_column_width"], _TABLE_SETTINGS["max_column_width"]))
        widths.append(width)
    return widths


def _color_scale(series: pd.Series) -> np.ndarray:
    thresholds = _TABLE_SETTINGS["thresholds"]
    colors = _TABLE_SETTINGS["colors"]
    return np.select(
        [series >= thresholds[1], series >= thresholds[0]],
        [colors[2], colors[1]],
        default=colors[0],
    )


def build_fold_change_table(
    consolidated: pd.DataFrame,
    *,
    reference_gene: str,
    decimals: int | None = None,
) -> go.Figure:
    """Genera una tabla comparativa de fold change por método."""
    use_decimals = decimals if decimals is not None else _TABLE_SETTINGS["decimals"]
    cols = ["target", "fold_change_promedio", "fold_change_gen_ref"]
    _ensure_columns(consolidated, cols)
    data = consolidated.loc[:, cols].copy()
    data = data.sort_values("target").reset_index(drop=True)
    data["fold_change_promedio"] = data["fold_change_promedio"].round(use_decimals)
    data["fold_change_gen_ref"] = data["fold_change_gen_ref"].round(use_decimals)

    display_df = data.astype({"fold_change_promedio": str, "fold_change_gen_ref": str})
    color_prom = _color_scale(data["fold_change_promedio"])
    color_ref = _color_scale(data["fold_change_gen_ref"])

    fig = go.Figure(data=[go.Table(
        columnwidth=list(_calculate_column_widths(display_df)),
        header=dict(
            values=[f"<b>{name}</b>" for name in display_df.columns],
            fill_color=_TABLE_SETTINGS["header_fill"],
            font=dict(color=_TABLE_SETTINGS["header_font_color"], size=_TABLE_SETTINGS["header_font_size"]),
            align="center",
        ),
        cells=dict(
            values=[display_df[col] for col in display_df.columns],
            fill_color=[["#EEEEEE"] * len(display_df), color_prom, color_ref],
            align="center",
            font=dict(size=_TABLE_SETTINGS["cell_font_size"]),
        ),
    )])

    fig.update_layout(
        title=dict(
            text=f"Comparación de fold change (gen de referencia: {reference_gene})",
            x=0.5,
            xanchor="center",
        ),
        margin=dict(t=80, b=20, l=20, r=20),
        height=320 + len(display_df) * 6,
    )
    return fig


def build_fold_change_chart(consolidated: pd.DataFrame) -> go.Figure:
    """Genera gráfico combinado ΔΔCT vs fold change para ambos métodos."""
    required = [
        "target",
        "delta_delta_ct_promedio",
        "delta_delta_ct_gen_ref",
        "fold_change_promedio",
        "fold_change_gen_ref",
    ]
    _ensure_columns(consolidated, required)
    data = consolidated.sort_values("target").reset_index(drop=True)

    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=data["target"],
        y=data["delta_delta_ct_promedio"],
        name="ΔΔCT (Promedio)",
        marker_color="#1f77b4",
        opacity=0.8,
        yaxis="y",
    ))
    fig.add_trace(go.Bar(
        x=data["target"],
        y=data["delta_delta_ct_gen_ref"],
        name="ΔΔCT (Gen ref)",
        marker_color="#ff7f0e",
        opacity=0.8,
        yaxis="y",
    ))
    fig.add_trace(go.Scatter(
        x=data["target"],
        y=data["fold_change_promedio"],
        name="Fold change (Promedio)",
        mode="markers+lines",
        marker=dict(color="#2ca02c", size=9, symbol="diamond"),
        line=dict(color="#2ca02c", width=2, dash="dot"),
        yaxis="y2",
    ))
    fig.add_trace(go.Scatter(
        x=data["target"],
        y=data["fold_change_gen_ref"],
        name="Fold change (Gen ref)",
        mode="markers+lines",
        marker=dict(color="#d62728", size=9, symbol="diamond"),
        line=dict(color="#d62728", width=2, dash="dot"),
        yaxis="y2",
    ))

    fig.update_layout(
        title=dict(text="Comparativo ΔΔCT vs Fold change", x=0.5, xanchor="center"),
        template="plotly_white",
        barmode="group",
        yaxis=dict(
            title="ΔΔCT",
            showgrid=True,
            gridcolor="lightgray",
        ),
        yaxis2=dict(
            title="Fold change (log)",
            overlaying="y",
            side="right",
            type="log",
        ),
        legend=dict(orientation="h", x=1, xanchor="right", y=1.1, bgcolor="rgba(255,255,255,0.8)"),
        height=640,
        margin=dict(t=80, b=120, l=40, r=60),
    )

    fig.add_hline(y=0.0, line_width=1, line_dash="dot", line_color="#666666")
    fig.add_hline(y=1.0, line_width=1, line_dash="dot", line_color="#d62728", yref="y2")
    return fig


__all__ = [
    "build_fold_change_table",
    "build_fold_change_chart",
]
