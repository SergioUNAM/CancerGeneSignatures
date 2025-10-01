from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List, Optional, Sequence

import pandas as pd

from app.core.imputation import procesar_ct_column
from app.core.fold_change import compute_fold_change


class FoldChangePreparationError(Exception):
    """Error al preparar datos para el cálculo de fold change."""


@dataclass(frozen=True)
class ImputationOutput:
    controles: pd.DataFrame
    muestras: pd.DataFrame
    summary: pd.DataFrame
    message: Optional[str]
    policy: str
    remaining_missing_ctrl: int
    remaining_missing_samples: int

    @property
    def has_missing(self) -> bool:
        return bool(self.remaining_missing_ctrl or self.remaining_missing_samples)


@dataclass(frozen=True)
class QualityMetrics:
    ctrl_targets: set[str]
    sample_targets: set[str]
    common_targets: set[str]
    ctrl_nan_ratio: float
    sample_nan_ratio: float
    ctrl_counts: pd.Series
    sample_counts: pd.Series

    def eligible_genes(self, minimum: int) -> List[str]:
        if minimum <= 1:
            return sorted(self.common_targets)
        return sorted(
            t for t in self.common_targets
            if self.ctrl_counts.get(t, 0) >= minimum and self.sample_counts.get(t, 0) >= minimum
        )


@dataclass(frozen=True)
class FoldChangeResult:
    computation: object
    expression_table: pd.DataFrame


def apply_undetermined_policy(
    controles: pd.DataFrame,
    muestras: pd.DataFrame,
    *,
    policy: str,
    fixed_value: Optional[float] = None,
    group_columns: Optional[Sequence[str]] = None,
) -> ImputationOutput:
    """Aplica la política de 'Undetermined/ND' y devuelve dataframes preparados."""

    policy = (policy or "nan").lower()
    ctrl = controles.copy()
    samp = muestras.copy()
    group_by = list(group_columns) if group_columns else None

    kwargs = {"columna_ct": "ct"}
    if group_by:
        kwargs["by"] = group_by

    message: Optional[str] = None
    if policy in {"ctmax", "value"}:
        try:
            if policy == "value":
                kwargs["max_ct"] = float(fixed_value) if fixed_value is not None else 40.0
            ctrl, samp = procesar_ct_column(ctrl, samp, **kwargs)
        except ValueError as exc:  # mantener semántica anterior
            raise FoldChangePreparationError(str(exc)) from exc
        message = ctrl.attrs.get("imputation_info")
    elif policy == "nan":
        message = "La política actual mantiene los valores 'Undetermined' como NaN."
    else:
        raise FoldChangePreparationError(f"Política de undetermined no soportada: {policy}")

    summary = _build_imputation_summary(ctrl, samp)
    remaining_ctrl = int(ctrl["ct"].isna().sum()) if "ct" in ctrl.columns else 0
    remaining_samp = int(samp["ct"].isna().sum()) if "ct" in samp.columns else 0

    return ImputationOutput(
        controles=ctrl,
        muestras=samp,
        summary=summary,
        message=message,
        policy=policy,
        remaining_missing_ctrl=remaining_ctrl,
        remaining_missing_samples=remaining_samp,
    )


def compute_quality_metrics(controles: pd.DataFrame, muestras: pd.DataFrame) -> QualityMetrics:
    """Calcula métricas previas al fold change (genes, NaN, conteos)."""

    ctrl_targets = set(controles.get("target", pd.Series(dtype=str)).dropna().astype(str))
    sample_targets = set(muestras.get("target", pd.Series(dtype=str)).dropna().astype(str))
    common_target = ctrl_targets.intersection(sample_targets)

    ctrl_nan_ratio = float(controles["ct"].isna().mean()) if "ct" in controles.columns and not controles.empty else 1.0
    sample_nan_ratio = float(muestras["ct"].isna().mean()) if "ct" in muestras.columns and not muestras.empty else 1.0

    ctrl_counts = controles.dropna(subset=["ct"]).groupby("target")["ct"].size() if not controles.empty else pd.Series(dtype=int)
    sample_counts = muestras.dropna(subset=["ct"]).groupby("target")["ct"].size() if not muestras.empty else pd.Series(dtype=int)

    return QualityMetrics(
        ctrl_targets=ctrl_targets,
        sample_targets=sample_targets,
        common_targets=common_target,
        ctrl_nan_ratio=ctrl_nan_ratio,
        sample_nan_ratio=sample_nan_ratio,
        ctrl_counts=ctrl_counts,
        sample_counts=sample_counts,
    )


def compute_fold_change_with_expression(
    controles: pd.DataFrame,
    muestras: pd.DataFrame,
    *,
    fold_source: str = "promedios",
    precomputed: Optional[object] = None,
) -> FoldChangeResult:
    """Calcula (o reutiliza) el fold change y devuelve tabla consolidada para visualizaciones."""

    fc = precomputed or compute_fold_change(controles, muestras)
    use_col = "fold_change_promedio" if fold_source == "promedios" else "fold_change_gen_ref"
    df_expr = fc.consolidated[["target", use_col]].rename(columns={use_col: "fold_change"}).copy()
    df_expr["nivel_expresion"] = pd.cut(
        df_expr["fold_change"],
        bins=[-float("inf"), 1.0, 2.0, float("inf")],
        labels=["subexpresado", "estable", "sobreexpresado"],
        right=False,
    )
    return FoldChangeResult(computation=fc, expression_table=df_expr)


def _build_imputation_summary(ctrl: pd.DataFrame, samp: pd.DataFrame) -> pd.DataFrame:
    if "ct_imputed" not in ctrl.columns and "ct_imputed" not in samp.columns:
        return pd.DataFrame()

    ctrl_imp = ctrl.loc[ctrl.get("ct_imputed", False)].copy()
    ctrl_imp["grupo"] = "control"
    samp_imp = samp.loc[samp.get("ct_imputed", False)].copy()
    samp_imp["grupo"] = "muestra"
    return pd.concat([ctrl_imp, samp_imp], ignore_index=True, sort=False)


__all__ = [
    "apply_undetermined_policy",
    "compute_quality_metrics",
    "compute_fold_change_with_expression",
    "FoldChangePreparationError",
    "ImputationOutput",
    "QualityMetrics",
    "FoldChangeResult",
]
