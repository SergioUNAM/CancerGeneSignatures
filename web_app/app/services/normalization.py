from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import pandas as pd

from app.core.reference_normalization import (
    AdvancedNormalizationError,
    AdvancedNormalizationResult,
    run_advanced_normalization,
)


@dataclass(frozen=True)
class AdvancedNormalizationParams:
    alpha: float = 0.05
    n_candidates: int = 20
    k_refs: int = 2
    bootstrap_iter: int = 200
    permutation_iter: int = 100
    random_seed: Optional[int] = 123


def build_total_dataframe(controles: pd.DataFrame, muestras: pd.DataFrame) -> pd.DataFrame:
    """Une controles y muestras agregando la etiqueta 'tipo'."""
    ctrl = controles.copy()
    samp = muestras.copy()
    ctrl["tipo"] = "Control"
    samp["tipo"] = "Muestra"
    cols = sorted(set(ctrl.columns).union(samp.columns))
    ctrl = ctrl.reindex(columns=cols)
    samp = samp.reindex(columns=cols)
    combined = pd.concat([ctrl, samp], ignore_index=True)
    return combined


def execute_advanced_normalization(
    controles: pd.DataFrame,
    muestras: pd.DataFrame,
    params: Optional[AdvancedNormalizationParams] = None,
) -> AdvancedNormalizationResult:
    params = params or AdvancedNormalizationParams()
    df_total = build_total_dataframe(controles, muestras)
    result = run_advanced_normalization(
        df_total,
        alpha=params.alpha,
        n_candidates=params.n_candidates,
        k_refs=params.k_refs,
        bootstrap_iter=params.bootstrap_iter,
        permutation_iter=params.permutation_iter,
        random_seed=params.random_seed,
    )
    return result


__all__ = [
    "AdvancedNormalizationParams",
    "AdvancedNormalizationError",
    "AdvancedNormalizationResult",
    "build_total_dataframe",
    "execute_advanced_normalization",
]
