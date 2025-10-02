from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import pandas as pd

from app.core.reference_normalization import (
    AdvancedNormalizationError,
    AdvancedNormalizationResult,
    run_advanced_normalization,
)
from app.core.fold_change import compute_fold_change
from app.core.reference_normalization import build_heatmap_matrix


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


# ===================== Utilidades de comparación de métodos ======================

def _normalize_by_global_mean(df_total: pd.DataFrame) -> pd.DataFrame:
    """Normaliza por promedio global por test: log2_rel_expr = -(Ct - mean(Ct)_test)."""
    out = df_total.copy()
    mean_ct = out.groupby("test")["ct"].mean().rename("ref_ct")
    out = out.merge(mean_ct, on="test", how="left")
    out["delta_ct"] = out["ct"] - out["ref_ct"]
    out["log2_rel_expr"] = -out["delta_ct"]
    return out


def _normalize_by_single_reference(df_total: pd.DataFrame, reference_gene: str) -> Tuple[pd.DataFrame, list[str]]:
    """Normaliza usando un único gen de referencia por test.

    - Calcula Ct_ref por test como la media de Ct del gen de referencia (si existe en dicho test).
    - Si falta en algún test, deja NaN y reporta advertencias.
    """
    out = df_total.copy()
    ref_rows = out[out["target"].astype(str) == str(reference_gene)]
    warn: list[str] = []
    if ref_rows.empty:
        warn.append("El gen de referencia no está presente en el DataFrame combinado.")
        out["ref_ct"] = pd.NA
    else:
        ref_ct = ref_rows.groupby("test")["ct"].mean().rename("ref_ct")
        out = out.merge(ref_ct, on="test", how="left")
        missing = out.loc[out["ref_ct"].isna(), "test"].nunique()
        total = out["test"].nunique()
        if missing:
            warn.append(
                f"Faltan Ct de referencia en {missing}/{total} tests; se omiten en el cálculo de expresión."
            )
    out["delta_ct"] = out["ct"] - out["ref_ct"]
    out["log2_rel_expr"] = -out["delta_ct"]
    return out, warn


def build_method_matrices(
    controles: pd.DataFrame,
    muestras: pd.DataFrame,
    adv_result: AdvancedNormalizationResult,
    *,
    genes_for_heatmap: Optional[Sequence[str]] = None,
) -> Tuple[Dict[str, pd.DataFrame], str, Dict[str, List[str]]]:
    """Construye matrices (genes x tests) de log2_rel_expr para tres métodos:
    - "advanced": referencias múltiples elegidas por el pipeline avanzado
    - "global_mean": promedio global por test
    - "refgene": un único gen de referencia básico (más estable) por test

    Devuelve: (dict de matrices, nombre_gen_ref_basico, advertencias_por_método)
    """
    # Base combinada
    df_total = build_total_dataframe(controles, muestras)

    # 1) Avanzado — ya tenemos df_norm; reutilizar matriz generada o reconstruir según genes solicitados
    adv_norm = adv_result.df_norm.copy()
    adv_refs = list(adv_result.reference_result.references)

    if genes_for_heatmap is None:
        # si no se especifica, usa la lista del resultado (df_heatmap ya contiene orden base)
        adv_mat = adv_result.df_heatmap.copy()
        heat_genes = adv_mat.index.tolist()
    else:
        heat_genes = list(dict.fromkeys(str(g) for g in genes_for_heatmap))
        adv_mat = build_heatmap_matrix(adv_norm, heat_genes)

    heat_genes = list(dict.fromkeys(heat_genes))
    if not heat_genes:
        heat_genes = list(dict.fromkeys(adv_norm["target"].dropna().astype(str).tolist()))
        adv_mat = build_heatmap_matrix(adv_norm, heat_genes)
    all_tests = sorted(df_total["test"].dropna().astype(str).unique().tolist())

    def _align_matrix(mat: Optional[pd.DataFrame]) -> pd.DataFrame:
        out = (mat.copy() if isinstance(mat, pd.DataFrame) else pd.DataFrame())
        if heat_genes:
            out = out.reindex(heat_genes)
        if all_tests:
            out = out.reindex(columns=all_tests)
        return out

    adv_mat = _align_matrix(adv_mat)

    warnings: Dict[str, List[str]] = {"advanced": [], "global_mean": [], "refgene": []}

    # 2) Global mean por test (manejar fallos sin romper la interfaz)
    gm_mat: pd.DataFrame = pd.DataFrame()
    try:
        gm_norm = _normalize_by_global_mean(df_total)
        gm_mat = build_heatmap_matrix(gm_norm, heat_genes)
    except Exception as err:
        warnings["global_mean"].append(
            f"No se pudo calcular el método de promedio global: {err}"
        )
    gm_mat = _align_matrix(gm_mat)

    # 3) Gen de referencia básico — obtener gen más estable del método básico
    basic_ref_gene = adv_refs[0] if adv_refs else ""
    ref_mat: pd.DataFrame = pd.DataFrame()
    try:
        basic_fc = compute_fold_change(controles, muestras)
        basic_ref_gene = basic_fc.reference_gene
        ref_norm, ref_warn = _normalize_by_single_reference(df_total, basic_ref_gene)
        if ref_warn:
            warnings["refgene"].extend(ref_warn)
        ref_mat = build_heatmap_matrix(ref_norm, heat_genes)
    except Exception as err:
        warnings["refgene"].append(
            f"No se pudo calcular el método del gen de referencia básico: {err}"
        )
    ref_mat = _align_matrix(ref_mat)

    return {"advanced": adv_mat, "global_mean": gm_mat, "refgene": ref_mat}, basic_ref_gene, warnings


__all__ = [
    "AdvancedNormalizationParams",
    "AdvancedNormalizationError",
    "AdvancedNormalizationResult",
    "build_total_dataframe",
    "execute_advanced_normalization",
    "build_method_matrices",
]
