from __future__ import annotations

import math
import pandas as pd
import plotly.graph_objects as go


def fc_comparison_table(df: pd.DataFrame) -> go.Figure:
    """
    Genera una tabla comparativa de fold change por método (promedio vs gen de referencia).
    Espera columnas: ['target', 'fold_change_promedio', 'fold_change_gen_ref'] en `df`.
    """
    cols = ['target', 'fold_change_promedio', 'fold_change_gen_ref']
    data = df[cols].copy()
    data = data.sort_values(by='target')
    data['fold_change_promedio'] = data['fold_change_promedio'].round(2)
    data['fold_change_gen_ref'] = data['fold_change_gen_ref'].round(2)

    header_values = ['Gen (target)', 'FC (promedio)', 'FC (gen de referencia)']
    cell_values = [
        data['target'].tolist(),
        data['fold_change_promedio'].tolist(),
        data['fold_change_gen_ref'].tolist(),
    ]

    fig = go.Figure(data=[go.Table(
        header=dict(
            values=header_values,
            fill_color='#2C3E50',
            font_color='white',
            align='center',
            height=36,
        ),
        cells=dict(
            values=cell_values,
            align=['left', 'center', 'center'],
            height=30,
        )
    )])

    fig.update_layout(
        title='Comparativa de fold change',
        margin=dict(l=10, r=10, t=40, b=10)
    )
    return fig

