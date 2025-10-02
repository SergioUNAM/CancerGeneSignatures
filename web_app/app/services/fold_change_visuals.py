from __future__ import annotations

from typing import Iterable, Sequence, Tuple

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


# ===================== Nuevas comparativas de Fold Change =====================

def _log2fc_columns(consolidated: pd.DataFrame) -> pd.DataFrame:
    required = [
        "target",
        "delta_delta_ct_promedio",
        "delta_delta_ct_gen_ref",
        "fold_change_promedio",
        "fold_change_gen_ref",
    ]
    _ensure_columns(consolidated, required)
    df = consolidated.copy()
    # log2FC = -ΔΔCt (consistente con fold_change = 2**(-ΔΔCt))
    df["log2fc_promedio"] = -df["delta_delta_ct_promedio"].astype(float)
    df["log2fc_gen_ref"] = -df["delta_delta_ct_gen_ref"].astype(float)
    return df


def build_fc_methods_bars(
    consolidated: pd.DataFrame,
    *,
    top_n: int = 30,
    order_by: str = "max_abs_log2fc",
) -> go.Figure:
    """Barras agrupadas de log2FC por gen para ambos métodos.

    - top_n: limita a los genes con mayor |log2FC| (máximo entre métodos).
    - order_by: actualmente soporta "max_abs_log2fc".
    """
    df = _log2fc_columns(consolidated)
    score = df[["log2fc_promedio", "log2fc_gen_ref"]].abs().max(axis=1)
    df = df.assign(_score=score)
    df = df.sort_values("_score", ascending=False).head(int(top_n))
    df = df.sort_values("log2fc_promedio", ascending=True)  # orden estable para lectura

    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            x=df["log2fc_promedio"],
            y=df["target"],
            name="log2FC (Promedio)",
            orientation="h",
            marker_color="#1f77b4",
            opacity=0.85,
        )
    )
    fig.add_trace(
        go.Bar(
            x=df["log2fc_gen_ref"],
            y=df["target"],
            name="log2FC (Gen ref)",
            orientation="h",
            marker_color="#ff7f0e",
            opacity=0.85,
        )
    )

    fig.update_layout(
        barmode="overlay",
        title=dict(text="Comparativa de log2FC por gen (Top-N)", x=0.5, xanchor="center"),
        template="plotly_white",
        xaxis_title="log2FC",
        yaxis_title="Gen",
        height=max(500, 18 * len(df)),
        margin=dict(t=60, b=40, l=80, r=20),
        legend=dict(orientation="h", x=1, xanchor="right", y=1.05),
    )
    fig.add_vline(x=0.0, line_width=1, line_dash="dot", line_color="#777")
    return fig


def build_fc_methods_scatter(fc_table: pd.DataFrame) -> go.Figure:
    """Gráfica comparativa única: log2FC avanzada vs promedio (color = gen de referencia)."""

    required = [
        "target",
        "log2fc_advanced",
        "log2fc_promedio",
        "log2fc_gen_ref",
    ]
    _ensure_columns(fc_table, required)

    df = fc_table.loc[:, required].copy()
    df.dropna(subset=["log2fc_advanced", "log2fc_promedio"], inplace=True)

    diff_prom = (df["log2fc_advanced"] - df["log2fc_promedio"]).abs()
    diff_ref = (df["log2fc_advanced"] - df["log2fc_gen_ref"]).abs()
    size = (diff_prom.combine(diff_ref, max) * 12).fillna(6.0) + 6.0

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=df["log2fc_advanced"],
            y=df["log2fc_promedio"],
            mode="markers",
            name="Genes",
            marker=dict(
                size=size,
                sizemode="diameter",
                color=df["log2fc_gen_ref"],
                colorscale="RdBu",
                showscale=True,
                colorbar=dict(title="log2FC (Gen ref)"),
                line=dict(width=0.5, color="#444"),
            ),
            text=df["target"],
            hovertemplate=(
                "<b>%{text}</b><br>log2FC Avanzada: %{x:.2f}<br>"
                "log2FC Promedio: %{y:.2f}<br>log2FC Gen ref: %{marker.color:.2f}"
                "<extra></extra>"
            ),
        )
    )

    all_vals = pd.concat([df["log2fc_advanced"], df["log2fc_promedio"]])
    if not all_vals.empty and all_vals.notna().any():
        lo = float(np.nanquantile(all_vals, 0.02))
        hi = float(np.nanquantile(all_vals, 0.98))
        fig.add_trace(
            go.Scatter(
                x=[lo, hi],
                y=[lo, hi],
                mode="lines",
                name="y = x",
                line=dict(color="#666", dash="dot"),
            )
        )

    fig.update_layout(
        title=dict(text="Comparativa única de log2FC (Avanzada vs Promedio)", x=0.5, xanchor="center"),
        template="plotly_white",
        xaxis_title="log2FC (Avanzada)",
        yaxis_title="log2FC (Promedio)",
        height=560,
        margin=dict(t=60, b=60, l=60, r=20),
        legend=dict(orientation="h", x=1, xanchor="right", y=1.05),
    )
    return fig


def summarize_fc_methods(
    fc_table: pd.DataFrame,
    *,
    thresholds: Tuple[float, float] = (1.5, 2.0),
    top_k_discrepancies: int = 10,
) -> dict:
    """Resumen cuantitativo de concordancia y discrepancias entre los tres métodos."""

    required = [
        "target",
        "log2fc_advanced",
        "log2fc_promedio",
        "log2fc_gen_ref",
    ]
    _ensure_columns(fc_table, required)

    df = fc_table.loc[:, required].dropna(subset=["log2fc_advanced", "log2fc_promedio", "log2fc_gen_ref"], how="all")
    if df.empty:
        return {
            "pairwise_metrics": {},
            "counts": {},
            "top_discrepancies": [],
        }

    metrics = {}
    pairs = {
        ("advanced", "promedio"): (df["log2fc_advanced"], df["log2fc_promedio"]),
        ("advanced", "gen_ref"): (df["log2fc_advanced"], df["log2fc_gen_ref"]),
        ("promedio", "gen_ref"): (df["log2fc_promedio"], df["log2fc_gen_ref"]),
    }

    for (a_name, b_name), (series_a, series_b) in pairs.items():
        mask = series_a.notna() & series_b.notna()
        if mask.sum() == 0:
            metrics[(a_name, b_name)] = {"pearson_r": np.nan, "rmse": np.nan, "mae": np.nan, "n": 0}
            continue
        a = series_a[mask].astype(float)
        b = series_b[mask].astype(float)
        pearson_r = float(np.corrcoef(a, b)[0, 1]) if len(a) > 1 else np.nan
        rmse = float(np.sqrt(np.nanmean((a - b) ** 2)))
        mae = float(np.nanmean(np.abs(a - b)))
        metrics[(a_name, b_name)] = {
            "pearson_r": pearson_r,
            "rmse": rmse,
            "mae": mae,
            "n": int(len(a)),
        }

    counts = {}
    log2_cols = {
        "advanced": df["log2fc_advanced"],
        "promedio": df["log2fc_promedio"],
        "gen_ref": df["log2fc_gen_ref"],
    }
    for thr in thresholds:
        log2_thr = float(np.log2(thr))
        entry = {}
        masks = {}
        for name, series in log2_cols.items():
            mask = series.abs() >= log2_thr
            entry[name] = int(mask.sum())
            masks[name] = mask
        # intersecciones
        all_mask = np.logical_and.reduce([masks[name] for name in masks])
        any_mask = np.logical_or.reduce([masks[name] for name in masks])
        entry["todos"] = int(all_mask.sum())
        entry["alguno"] = int(any_mask.sum())
        counts[f">={thr}x"] = entry

    diff_prom = (df["log2fc_advanced"] - df["log2fc_promedio"]).abs().fillna(0)
    diff_ref = (df["log2fc_advanced"] - df["log2fc_gen_ref"]).abs().fillna(0)
    diff_max = diff_prom.combine(diff_ref, max)
    order = np.argsort(-diff_max.values)[: int(top_k_discrepancies)]
    top_disc = [
        {
            "gene": str(df["target"].iloc[i]),
            "log2fc_advanced": float(df["log2fc_advanced"].iloc[i]),
            "log2fc_promedio": float(df["log2fc_promedio"].iloc[i]),
            "log2fc_gen_ref": float(df["log2fc_gen_ref"].iloc[i]),
            "max_abs_diff": float(diff_max.iloc[i]),
        }
        for i in order if i < len(df)
    ]

    return {
        "pairwise_metrics": metrics,
        "counts": counts,
        "top_discrepancies": top_disc,
    }


def _classify_fold_change(series: pd.Series) -> pd.Categorical:
    bins = [-float("inf"), 1.0, 2.0, float("inf")]
    labels = ["subexpresado", "estable", "sobreexpresado"]
    return pd.cut(series.astype(float), bins=bins, labels=labels, right=False)


def build_expression_classification_table(fc_table: pd.DataFrame) -> pd.DataFrame:
    required = [
        "target",
        "fold_change_advanced",
        "fold_change_promedio",
        "fold_change_gen_ref",
    ]
    _ensure_columns(fc_table, required)

    table = fc_table.loc[:, required].copy()
    table.rename(columns={
        "fold_change_advanced": "fc_advanced",
        "fold_change_promedio": "fc_promedio",
        "fold_change_gen_ref": "fc_gen_ref",
    }, inplace=True)

    table["clasificacion_advanced"] = _classify_fold_change(table["fc_advanced"])
    table["clasificacion_promedio"] = _classify_fold_change(table["fc_promedio"])
    table["clasificacion_gen_ref"] = _classify_fold_change(table["fc_gen_ref"])

    output_columns = [
        "target",
        "clasificacion_advanced",
        "clasificacion_promedio",
        "clasificacion_gen_ref",
        "fc_advanced",
        "fc_promedio",
        "fc_gen_ref",
    ]

    log2_cols = {"log2fc_advanced", "log2fc_promedio", "log2fc_gen_ref"}
    if log2_cols.issubset(fc_table.columns):
        table["delta_log2fc_adv_prom"] = fc_table["log2fc_advanced"] - fc_table["log2fc_promedio"]
        table["delta_log2fc_adv_ref"] = fc_table["log2fc_advanced"] - fc_table["log2fc_gen_ref"]
        table["delta_log2fc_prom_ref"] = fc_table["log2fc_promedio"] - fc_table["log2fc_gen_ref"]
        table["max_abs_delta_log2fc"] = (
            table[["delta_log2fc_adv_prom", "delta_log2fc_adv_ref", "delta_log2fc_prom_ref"]]
            .abs()
            .max(axis=1)
        )
        output_columns.extend([
            "delta_log2fc_adv_prom",
            "delta_log2fc_adv_ref",
            "delta_log2fc_prom_ref",
            "max_abs_delta_log2fc",
        ])

    return table[output_columns].copy()


def build_classification_summary_chart(classification_table: pd.DataFrame) -> go.Figure:
    categories = ["subexpresado", "estable", "sobreexpresado"]
    method_columns = {
        "clasificacion_advanced": "Avanzada",
        "clasificacion_promedio": "Promedio",
        "clasificacion_gen_ref": "Gen ref",
    }

    counts = {
        label: classification_table[col].value_counts().reindex(categories, fill_value=0)
        for col, label in method_columns.items()
    }
    counts_df = pd.DataFrame(counts, index=categories)

    fig = go.Figure()
    for method in counts_df.columns:
        fig.add_trace(
            go.Bar(
                x=counts_df.index,
                y=counts_df[method],
                name=method,
                text=counts_df[method],
                textposition="outside",
            )
        )

    fig.update_layout(
        title=dict(text="Distribución de niveles de expresión por método", x=0.5, xanchor="center"),
        template="plotly_white",
        barmode="group",
        xaxis_title="Clasificación de expresión",
        yaxis_title="Número de genes",
        height=420,
        margin=dict(t=60, b=60, l=40, r=20),
        legend=dict(orientation="h", x=1, xanchor="right", y=1.1),
    )
    return fig


def build_three_way_venn(
    sets_by_method: dict[str, Iterable[str]],
    labels: Tuple[str, str, str],
) -> go.Figure:
    """Venn diagram (3 sets) highlighting overlaps between methods."""

    a_label, b_label, c_label = labels
    set_a = set(map(str, sets_by_method.get("advanced", [])))
    set_b = set(map(str, sets_by_method.get("global_mean", [])))
    set_c = set(map(str, sets_by_method.get("refgene", [])))

    only_a = len(set_a - set_b - set_c)
    only_b = len(set_b - set_a - set_c)
    only_c = len(set_c - set_a - set_b)
    ab = len((set_a & set_b) - set_c)
    ac = len((set_a & set_c) - set_b)
    bc = len((set_b & set_c) - set_a)
    abc = len(set_a & set_b & set_c)

    fig = go.Figure()
    circles = [
        dict(type="circle", x0=-1.1, y0=-0.9, x1=0.1, y1=1.1, line=dict(color="#1f77b4", width=3), opacity=0.25),
        dict(type="circle", x0=-0.1, y0=-0.9, x1=1.1, y1=1.1, line=dict(color="#ff7f0e", width=3), opacity=0.25),
        dict(type="circle", x0=-0.6, y0=-0.1, x1=0.6, y1=1.9, line=dict(color="#2ca02c", width=3), opacity=0.25),
    ]
    fig.update_layout(shapes=circles)

    annotations = [
        dict(x=-0.85, y=-0.2, text=str(only_a), showarrow=False),
        dict(x=0.85, y=-0.2, text=str(only_b), showarrow=False),
        dict(x=0.0, y=1.4, text=str(only_c), showarrow=False),
        dict(x=0.0, y=-0.4, text=str(ab), showarrow=False),
        dict(x=-0.35, y=0.75, text=str(ac), showarrow=False),
        dict(x=0.35, y=0.75, text=str(bc), showarrow=False),
        dict(x=0.0, y=0.25, text=str(abc), font=dict(size=18, color="#111"), showarrow=False),
        dict(x=-0.9, y=1.2, text=a_label, showarrow=False, font=dict(color="#1f77b4")),
        dict(x=0.9, y=1.2, text=b_label, showarrow=False, font=dict(color="#ff7f0e")),
        dict(x=0.0, y=1.9, text=c_label, showarrow=False, font=dict(color="#2ca02c")),
    ]
    fig.update_layout(
        annotations=annotations,
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        showlegend=False,
        margin=dict(t=40, b=40, l=40, r=40),
        height=500,
        width=600,
        plot_bgcolor="white",
        paper_bgcolor="white",
        title=dict(text="Genes compartidos entre métodos", x=0.5, xanchor="center"),
    )
    return fig


__all__ += [
    "build_fc_methods_bars",
    "build_fc_methods_scatter",
    "summarize_fc_methods",
    "build_expression_classification_table",
    "build_classification_summary_chart",
    "build_three_way_venn",
]
