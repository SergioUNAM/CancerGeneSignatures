from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple

import numpy as np
import pandas as pd


@dataclass
class FoldChangeResult:
    consolidated: pd.DataFrame
    by_means: pd.DataFrame
    by_refgene: pd.DataFrame
    reference_gene: str


def _stats_by_group(df: pd.DataFrame, group_name: str) -> pd.DataFrame:
    stats = (
        df.groupby('target', as_index=False)
        .agg(ct_promedio=('ct', 'mean'), ct_std=('ct', 'std'))
        .rename(columns={
            'ct_promedio': f'ct_promedio_{group_name}',
            'ct_std': f'ct_std_{group_name}',
        })
    )
    return stats


def compute_fold_change(controles: pd.DataFrame, muestras: pd.DataFrame) -> FoldChangeResult:
    """
    Implementa los dos métodos descritos en el notebook:
    - Promedios globales por grupo
    - Gen de referencia (más estable por menor std promedio)
    Recibe DataFrames largos con columnas: ['test','target','ct'] para controles y muestras.
    """
    stats_controles = _stats_by_group(controles, 'controles')
    stats_muestras = _stats_by_group(muestras, 'muestras')

    # Promedios globales por grupo
    stats_controles['promedio_general_controles'] = stats_controles['ct_promedio_controles'].mean()
    stats_muestras['promedio_general_muestras'] = stats_muestras['ct_promedio_muestras'].mean()

    # ΔCt respecto al promedio global
    stats_controles['delta_ct_promedio_controles'] = (
        stats_controles['ct_promedio_controles'] - stats_controles['promedio_general_controles']
    )
    stats_muestras['delta_ct_promedio_muestras'] = (
        stats_muestras['ct_promedio_muestras'] - stats_muestras['promedio_general_muestras']
    )

    # Consolidado inicial
    df_consolidado = stats_controles.merge(
        stats_muestras, on='target', how='outer', suffixes=('_controles', '_muestras')
    )

    # Gen de referencia: menor std promedio
    df_consolidado['stability'] = (
        df_consolidado['ct_std_controles'] + df_consolidado['ct_std_muestras']
    ) / 2.0
    reference_gene = df_consolidado.loc[df_consolidado['stability'].idxmin(), 'target']

    # ΔCt respecto al gen de referencia
    ref_c = stats_controles.query('target == @reference_gene')['ct_promedio_controles'].values[0]
    ref_m = stats_muestras.query('target == @reference_gene')['ct_promedio_muestras'].values[0]

    stats_controles['delta_ct_gen_ref_controles'] = stats_controles['ct_promedio_controles'] - ref_c
    stats_muestras['delta_ct_gen_ref_muestras'] = stats_muestras['ct_promedio_muestras'] - ref_m

    # Consolidado final (con ambas métricas)
    df_consolidado = stats_controles.merge(
        stats_muestras, on='target', how='outer', suffixes=('_controles', '_muestras')
    )

    # ΔΔCt y Fold Change (promedio global)
    df_consolidado['delta_delta_ct_promedio'] = (
        df_consolidado['delta_ct_promedio_muestras'] - df_consolidado['delta_ct_promedio_controles']
    )
    df_consolidado['fold_change_promedio'] = 2 ** (-df_consolidado['delta_delta_ct_promedio'])

    # ΔΔCt y Fold Change (gen de referencia)
    df_consolidado['delta_delta_ct_gen_ref'] = (
        df_consolidado['delta_ct_gen_ref_muestras'] - df_consolidado['delta_ct_gen_ref_controles']
    )
    df_consolidado['fold_change_gen_ref'] = 2 ** (-df_consolidado['delta_delta_ct_gen_ref'])

    cols_promedio = ['target'] + [c for c in df_consolidado.columns if 'promedio' in c]
    cols_genref = ['target'] + [c for c in df_consolidado.columns if 'gen_ref' in c]
    by_means = df_consolidado[cols_promedio].rename(columns={'fold_change_promedio': 'fold_change'}).copy()
    by_ref = df_consolidado[cols_genref].rename(columns={'fold_change_gen_ref': 'fold_change'}).copy()

    return FoldChangeResult(
        consolidated=df_consolidado,
        by_means=by_means,
        by_refgene=by_ref,
        reference_gene=str(reference_gene),
    )

