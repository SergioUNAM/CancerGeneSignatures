from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go


def build_fc_detail_figure(
    consolidated: pd.DataFrame,
    df_expr: pd.DataFrame,
    *,
    exclude_stable: bool,
    reference_gene: str,
    scale: str = "log",
) -> go.Figure:
    """Genera la figura comparativa de ΔΔCt y Fold Change."""

    scale = (scale or "log").lower()
    y2_type = "log" if scale == "log" else "linear"

    try:
        if exclude_stable:
            targets_plot = df_expr.loc[df_expr["nivel_expresion"] != "estable", "target"]
            subset = consolidated[consolidated["target"].isin(targets_plot)]
        else:
            subset = consolidated
    except Exception:
        subset = consolidated

    x_vals = subset.get("target", pd.Series(dtype=str))
    ddct_mean = subset.get("delta_delta_ct_promedio", pd.Series(dtype=float))
    ddct_ref = subset.get("delta_delta_ct_gen_ref", pd.Series(dtype=float))
    fc_mean = subset.get("fold_change_promedio", pd.Series(dtype=float))
    fc_ref = subset.get("fold_change_gen_ref", pd.Series(dtype=float))

    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            x=x_vals,
            y=ddct_mean,
            name="ΔΔCT (Promedios)",
            marker_color="#1f77b4",
            opacity=0.85,
            yaxis="y",
            hovertemplate="Gen=%{x}<br>ΔΔCT Promedios=%{y:.3f}<extra></extra>",
        )
    )
    fig.add_trace(
        go.Bar(
            x=x_vals,
            y=ddct_ref,
            name="ΔΔCT (Gen Ref)",
            marker_color="#ff7f0e",
            opacity=0.85,
            yaxis="y",
            hovertemplate="Gen=%{x}<br>ΔΔCT GenRef=%{y:.3f}<extra></extra>",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=x_vals,
            y=fc_mean,
            name="Fold Change (Promedios)",
            mode="markers+lines",
            marker=dict(color="#2ca02c", size=8, symbol="diamond"),
            line=dict(color="#2ca02c", width=2, dash="dot"),
            yaxis="y2",
            hovertemplate="Gen=%{x}<br>FC Promedios=%{y:.3f}<extra></extra>",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=x_vals,
            y=fc_ref,
            name="Fold Change (Gen Ref)",
            mode="markers+lines",
            marker=dict(color="#d62728", size=8, symbol="diamond"),
            line=dict(color="#d62728", width=2, dash="dot"),
            yaxis="y2",
            hovertemplate="Gen=%{x}<br>FC GenRef=%{y:.3f}<extra></extra>",
        )
    )

    try:
        ref_mask = x_vals == reference_gene
        fig.add_trace(
            go.Scatter(
                x=x_vals[ref_mask],
                y=fc_ref[ref_mask],
                mode="markers",
                name="Ref gene",
                marker=dict(color="black", size=12, symbol="star"),
                yaxis="y2",
                hovertemplate="Gen de referencia=%{x}<br>FC GenRef=%{y:.3f}<extra></extra>",
            )
        )
    except Exception:
        pass

    fig.update_layout(
        title=dict(text="Análisis comparativo de métodos de cálculo", x=0.5),
        template="plotly_white",
        barmode="group",
        yaxis=dict(title="ΔΔCT", showgrid=True, gridcolor="lightgray"),
        yaxis2=dict(
            title=f"Fold Change ({scale})",
            overlaying="y",
            side="right",
            type=y2_type,
            showgrid=False,
        ),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        height=600,
        margin=dict(b=80, t=80, l=60, r=60),
    )
    return fig


def build_expression_distribution(
    df_expr: pd.DataFrame,
    *,
    exclude_stable: bool,
) -> go.Figure:
    """Genera la gráfica de barras con la distribución de niveles de expresión."""

    order_levels = ["estable", "subexpresado", "sobreexpresado"]
    if exclude_stable:
        df_plot = df_expr[df_expr["nivel_expresion"] != "estable"]
    else:
        df_plot = df_expr
    counts = df_plot["nivel_expresion"].value_counts().reindex(order_levels, fill_value=0)
    return px.bar(
        x=counts.index,
        y=counts.values,
        labels={"x": "Nivel de expresión", "y": "Frecuencia"},
        title="Distribución de niveles de expresión",
    )


def build_expression_treemap(df_expr: pd.DataFrame) -> Optional[go.Figure]:
    """Treemap de genes por nivel con tamaño proporcional a |log2FC|."""

    if df_expr.empty:
        return None

    treedata = df_expr.copy()
    try:
        treedata["log2fc"] = np.log2(treedata["fold_change"].clip(lower=1e-12))
    except Exception:
        treedata["log2fc"] = 0.0

    return px.treemap(
        treedata,
        path=["nivel_expresion", "target"],
        values=treedata["log2fc"].abs(),
        color="log2fc",
        color_continuous_scale="RdBu",
        color_continuous_midpoint=0.0,
        title="Genes por nivel (tamaño = |log2FC|, color = log2FC)",
    )


__all__ = [
    "build_fc_detail_figure",
    "build_expression_distribution",
    "build_expression_treemap",
]
