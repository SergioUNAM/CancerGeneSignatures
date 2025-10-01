from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.figure_factory as ff
from plotly.subplots import make_subplots
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist


def _zscore_rows(mat: pd.DataFrame) -> pd.DataFrame:
    arr = mat.values.astype(float)
    means = np.nanmean(arr, axis=1, keepdims=True)
    stds = np.nanstd(arr, axis=1, keepdims=True)
    stds[stds == 0] = 1.0
    z = (arr - means) / stds
    return pd.DataFrame(z, index=mat.index, columns=mat.columns)


def _order_by_dendrogram(mat: pd.DataFrame, metric: str = "euclidean", method: str = "average") -> tuple[pd.Index, pd.Index]:
    # Relleno simple para distancias (media de fila)
    filled = mat.copy()
    if filled.isna().values.any():
        row_means = filled.mean(axis=1)
        filled = filled.apply(lambda row: row.fillna(row_means[row.name]))

    row_dist = pdist(filled.values, metric=metric)
    col_dist = pdist(filled.values.T, metric=metric)
    row_linkage = linkage(row_dist, method=method)
    col_linkage = linkage(col_dist, method=method)
    row_order = mat.index[leaves_list(row_linkage)]
    col_order = mat.columns[leaves_list(col_linkage)]
    return row_order, col_order


def build_dendrogram_heatmap(
    matrix: pd.DataFrame,
    *,
    title: str = "",
    zscore_by_gene: bool = True,
    distance: str = "euclidean",
    linkage_method: str = "average",
) -> go.Figure:
    """Crea un heatmap con dendrogramas (filas y columnas) usando Plotly.

    - matrix: DataFrame index=genes, columns=tests, values=log2_rel_expr
    - zscore_by_gene: estandariza por fila para comparabilidad visual
    """
    if matrix is None or matrix.empty:
        return go.Figure()

    mat = matrix.copy()
    if zscore_by_gene:
        mat = _zscore_rows(mat)

    # Si no hay suficientes filas/columnas para clustering, devolver solo heatmap
    if mat.shape[0] < 2 or mat.shape[1] < 2:
        return go.Figure(
            data=[go.Heatmap(
                z=mat.values,
                x=mat.columns.tolist(),
                y=mat.index.tolist(),
                colorscale=[[0.0, "#2c7bb6"], [0.5, "#ffffb2"], [1.0, "#d7191c"]],
                zmid=0,
            )],
            layout=dict(title=title)
        )

    row_order, col_order = _order_by_dendrogram(mat, metric=distance, method=linkage_method)
    mat_ord = mat.loc[row_order, col_order]

    # Dendrogramas
    dendro_top = ff.create_dendrogram(mat_ord.values.T, orientation="bottom", labels=mat_ord.columns.tolist())
    dendro_side = ff.create_dendrogram(mat_ord.values, orientation="right", labels=mat_ord.index.tolist())

    for i in range(len(dendro_top["data"])):
        dendro_top["data"][i]["marker"]["color"] = "#606060"
    for i in range(len(dendro_side["data"])):
        dendro_side["data"][i]["marker"]["color"] = "#606060"

    fig = make_subplots(
        rows=2,
        cols=2,
        row_heights=[0.2, 0.8],
        column_widths=[0.2, 0.8],
        specs=[[{}, {}], [{}, {}]],
        horizontal_spacing=0.02,
        vertical_spacing=0.02,
    )

    for trace in dendro_top["data"]:
        fig.add_trace(trace, row=1, col=2)
    for trace in dendro_side["data"]:
        fig.add_trace(trace, row=2, col=1)

    heat = go.Heatmap(
        z=mat_ord.values,
        x=mat_ord.columns.tolist(),
        y=mat_ord.index.tolist(),
        colorscale=[[0.0, "#2c7bb6"], [0.5, "#ffffb2"], [1.0, "#d7191c"]],
        coloraxis="coloraxis",
        zmid=0,
    )
    fig.add_trace(heat, row=2, col=2)

    fig.update_layout(
        title=title,
        coloraxis=dict(colorscale=[[0.0, "#2c7bb6"], [0.5, "#ffffb2"], [1.0, "#d7191c"]]),
        showlegend=False,
        margin=dict(t=60, l=0, r=0, b=0),
    )

    # Ocultar ticks innecesarios en dendrogramas
    fig.update_xaxes(visible=False, row=1, col=2)
    fig.update_yaxes(visible=False, row=2, col=1)
    # Ejes de heatmap con ticks visibles
    fig.update_xaxes(tickangle=45, row=2, col=2)
    fig.update_yaxes(autorange="reversed", row=2, col=2)

    return fig


__all__ = [
    "build_dendrogram_heatmap",
]
