from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List, Optional, Sequence, Tuple

import pandas as pd

from app.core.io import LoadResult
from app.core.qpcr import (
    classify_by_prefixes as _classify_by_prefixes,
    classify_by_suffixes as _classify_by_suffixes,
    melt_wide_to_long,
)
from app.core.cleaning import drop_machine_controls

try:
    from app.core.qpcr import classify_by_regex as _classify_by_regex  # type: ignore
except Exception:  # pragma: no cover - fallback para instalaciones antiguas
    def _classify_by_regex(df_long, ctrl_pattern: str, sample_pattern: str):
        t = df_long.copy()
        t["test_str"] = t["test"].astype(str)
        ctrl = (
            t[t["test_str"].str.contains(ctrl_pattern, regex=True, na=False)]
            if ctrl_pattern else t.iloc[0:0]
        )
        samp = (
            t[t["test_str"].str.contains(sample_pattern, regex=True, na=False)]
            if sample_pattern else t.iloc[0:0]
        )
        return (
            ctrl.drop(columns=["test_str"]) if not ctrl.empty else ctrl,
            samp.drop(columns=["test_str"]) if not samp.empty else samp,
        )


@dataclass(frozen=True)
class ExtractionSummary:
    """Resumen básico derivado de los datos qPCR cargados."""

    sample_names: List[str]
    genes: List[str]
    wells: List[str]


@dataclass(frozen=True)
class ClassificationResult:
    """Resultado de una clasificación de pruebas en controles vs muestras."""

    controles: pd.DataFrame
    muestras: pd.DataFrame

    def is_valid(self) -> bool:
        return not self.controles.empty and not self.muestras.empty


def build_long_table(
    load_result: LoadResult,
    *,
    drop_controls: bool = True,
    control_markers: Sequence[str] = ("PPC", "RTC"),
) -> Tuple[pd.DataFrame, Optional[str]]:
    """Convierte el dataframe ancho en largo y filtra controles si es posible."""

    long_df = melt_wide_to_long(load_result.df)
    warning: Optional[str] = None
    if drop_controls:
        try:
            long_df = drop_machine_controls(long_df, column="target", controls=list(control_markers))
        except Exception as exc:
            warning = str(exc)
    return long_df, warning


def summarize_extraction(load_result: LoadResult) -> ExtractionSummary:
    """Construye el resumen presentado en la UI (tests, genes, pozos)."""

    sample_names: List[str] = []
    if isinstance(load_result.meta, dict):
        meta_tests = load_result.meta.get("sample_names")
        if isinstance(meta_tests, Iterable):
            sample_names = [str(v).strip() for v in meta_tests if str(v).strip()]

    if not sample_names:
        sample_names = [
            str(c).strip()
            for c in load_result.df.columns
            if str(c).strip() and c not in ("Well", "Target Name")
        ]

    genes_series = load_result.df.get("Target Name")
    genes = (
        [str(g).strip() for g in genes_series.dropna().unique().tolist() if str(g).strip()]
        if genes_series is not None
        else []
    )

    wells_series = load_result.df.get("Well")
    wells = (
        [str(w).strip() for w in wells_series.dropna().unique().tolist() if str(w).strip()]
        if wells_series is not None
        else []
    )

    return ExtractionSummary(sample_names=sample_names, genes=genes, wells=wells)


def classify_tests_by_prefixes(
    long_df: pd.DataFrame,
    control_prefixes: Sequence[str],
    sample_prefixes: Sequence[str],
) -> ClassificationResult:
    ctrl = _normalize_affixes(control_prefixes)
    samp = _normalize_affixes(sample_prefixes)
    ctrl_df, samp_df = _classify_by_prefixes(long_df, ctrl, samp)
    return ClassificationResult(controles=ctrl_df.copy(), muestras=samp_df.copy())


def classify_tests_by_suffixes(
    long_df: pd.DataFrame,
    control_suffixes: Sequence[str],
    sample_suffixes: Sequence[str],
) -> ClassificationResult:
    ctrl = _normalize_affixes(control_suffixes)
    samp = _normalize_affixes(sample_suffixes)
    ctrl_df, samp_df = _classify_by_suffixes(long_df, ctrl, samp)
    return ClassificationResult(controles=ctrl_df.copy(), muestras=samp_df.copy())


def classify_tests_by_regex(
    long_df: pd.DataFrame,
    control_pattern: str,
    sample_pattern: str,
) -> ClassificationResult:
    ctrl_df, samp_df = _classify_by_regex(long_df, control_pattern, sample_pattern)
    return ClassificationResult(controles=ctrl_df.copy(), muestras=samp_df.copy())


def classify_tests_by_selection(
    long_df: pd.DataFrame,
    selected_controls: Sequence[str],
    selected_samples: Sequence[str],
) -> ClassificationResult:
    controls = set(str(v) for v in selected_controls if str(v).strip())
    samples = set(str(v) for v in selected_samples if str(v).strip())
    ctrl_df = long_df[long_df["test"].astype(str).isin(controls)] if controls else long_df.iloc[0:0]
    samp_df = long_df[long_df["test"].astype(str).isin(samples)] if samples else long_df.iloc[0:0]
    return ClassificationResult(controles=ctrl_df.copy(), muestras=samp_df.copy())


def detect_collisions(ctrl_df: pd.DataFrame, samp_df: pd.DataFrame) -> List[str]:
    if ctrl_df.empty or samp_df.empty:
        return []
    ctrl = set(ctrl_df["test"].astype(str))
    samp = set(samp_df["test"].astype(str))
    return sorted(ctrl.intersection(samp))


def apply_collision_strategy(
    ctrl_df: pd.DataFrame,
    samp_df: pd.DataFrame,
    collisions: Sequence[str],
    strategy: str,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if not collisions:
        return ctrl_df, samp_df

    collisions_set = set(collisions)
    strategy = (strategy or "").lower()
    if strategy == "priorizar muestras":
        ctrl_df = ctrl_df[~ctrl_df["test"].astype(str).isin(collisions_set)]
    elif strategy == "excluir colisiones":
        ctrl_df = ctrl_df[~ctrl_df["test"].astype(str).isin(collisions_set)]
        samp_df = samp_df[~samp_df["test"].astype(str).isin(collisions_set)]
    else:  # por defecto priorizar controles
        samp_df = samp_df[~samp_df["test"].astype(str).isin(collisions_set)]
    return ctrl_df, samp_df


def _normalize_affixes(values: Sequence[str]) -> List[str]:
    norm = [str(v).strip() for v in values if str(v).strip()]
    return list(dict.fromkeys(norm))  # preserva orden sin duplicados


__all__ = [
    "ExtractionSummary",
    "build_long_table",
    "summarize_extraction",
    "ClassificationResult",
    "classify_tests_by_prefixes",
    "classify_tests_by_suffixes",
    "classify_tests_by_regex",
    "classify_tests_by_selection",
    "detect_collisions",
    "apply_collision_strategy",
]
