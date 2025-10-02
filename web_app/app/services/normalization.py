from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Sequence, Tuple

import pandas as pd
import numpy as np

from app.core.reference_normalization import (
    AdvancedNormalizationError,
    AdvancedNormalizationResult,
    DifferentialExpressionResult,
    run_advanced_normalization,
)
from app.core.fold_change import compute_fold_change, FoldChangeResult
from app.core.reference_normalization import build_heatmap_matrix, evaluate_differential_expression


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
) -> Tuple[Dict[str, pd.DataFrame], str]:
    """Construye matrices (genes x tests) de log2_rel_expr para tres métodos:
    - "advanced": referencias múltiples elegidas por el pipeline avanzado
    - "global_mean": promedio global por test
    - "refgene": un único gen de referencia básico (más estable) por test

    Devuelve: (dict de matrices, nombre_gen_ref_basico)
    """
    # Base combinada
    df_total = build_total_dataframe(controles, muestras)

    # 1) Avanzado — ya tenemos df_norm; reutilizar matriz generada o reconstruir según genes solicitados
    adv_norm = adv_result.df_norm.copy()
    refgene_adv = ", ".join(adv_result.reference_result.references)
    if genes_for_heatmap is None:
        # si no se especifica, usa la lista del resultado (df_heatmap ya contiene orden base)
        adv_mat = adv_result.df_heatmap.copy()
        heat_genes = adv_mat.index.tolist()
    else:
        heat_genes = list(genes_for_heatmap)
        adv_mat = build_heatmap_matrix(adv_norm, heat_genes)

    # 2) Global mean por test
    gm_norm = _normalize_by_global_mean(df_total)
    gm_mat = build_heatmap_matrix(gm_norm.rename(columns={"log2_rel_expr": "log2_rel_expr"}), heat_genes)

    # 3) Gen de referencia básico — obtener gen más estable del método básico
    basic_fc = compute_fold_change(controles, muestras)
    basic_ref_gene = basic_fc.reference_gene
    ref_norm, _ = _normalize_by_single_reference(df_total, basic_ref_gene)
    ref_mat = build_heatmap_matrix(ref_norm.rename(columns={"log2_rel_expr": "log2_rel_expr"}), heat_genes)

    return {"advanced": adv_mat, "global_mean": gm_mat, "refgene": ref_mat}, basic_ref_gene


__all__ = [
    "AdvancedNormalizationParams",
    "AdvancedNormalizationError",
    "AdvancedNormalizationResult",
    "build_total_dataframe",
    "execute_advanced_normalization",
    "build_method_matrices",
]

# ===================== Selección por método y persistencia ======================

def compute_basic_normalizations(
    controles: pd.DataFrame,
    muestras: pd.DataFrame,
    adv_result: AdvancedNormalizationResult,
) -> tuple[pd.DataFrame, pd.DataFrame, str]:
    """Devuelve (df_norm_global_mean, df_norm_refgene, basic_ref_gene).

    - "global_mean": normalización por promedio global por test.
    - "refgene": normalización por gen de referencia básico (elegido por std mínima).
    """
    df_total = build_total_dataframe(controles, muestras)

    # Global mean
    gm_norm = _normalize_by_global_mean(df_total)

    # Single reference (básico) tomado de compute_fold_change
    basic_fc = compute_fold_change(controles, muestras)
    basic_ref_gene = basic_fc.reference_gene
    ref_norm, _ = _normalize_by_single_reference(df_total, basic_ref_gene)

    return gm_norm, ref_norm, basic_ref_gene


def detect_significant_genes_by_method(
    controles: pd.DataFrame,
    muestras: pd.DataFrame,
    adv_result: AdvancedNormalizationResult,
    *,
    alpha: float = 0.05,
) -> tuple[dict[str, list[str]], dict[str, pd.DataFrame], str]:
    """Calcula genes diferencialmente expresados (q < alpha) para tres métodos.

    Devuelve:
    - dict "significant" por método: {"advanced", "global_mean", "refgene"} -> lista de genes
    - dict "tables" por método: tablas completas con estadísticas y q
    - basic_ref_gene: nombre del gen de referencia básico utilizado
    """
    # Avanzada: ya viene en adv_result
    adv_sig = list(adv_result.differential.significant_genes)
    adv_table = adv_result.differential.stats.copy()

    gm_norm, ref_norm, basic_ref_gene = compute_basic_normalizations(controles, muestras, adv_result)

    # Evaluación de significancia reaprovechando la función general
    gm_diff: DifferentialExpressionResult = evaluate_differential_expression(gm_norm, alpha=alpha)
    ref_diff: DifferentialExpressionResult = evaluate_differential_expression(ref_norm, alpha=alpha)

    significant = {
        "advanced": adv_sig,
        "global_mean": list(gm_diff.significant_genes),
        "refgene": list(ref_diff.significant_genes),
    }
    tables = {
        "advanced": adv_table,
        "global_mean": gm_diff.stats,
        "refgene": ref_diff.stats,
    }
    return significant, tables, basic_ref_gene


def build_heatmaps_by_method(
    controles: pd.DataFrame,
    muestras: pd.DataFrame,
    adv_result: AdvancedNormalizationResult,
    *,
    alpha: float = 0.05,
    top_n_fallback: int = 20,
) -> tuple[dict[str, pd.DataFrame], dict[str, list[str]], str]:
    """Construye una matriz heatmap por método, usando los genes seleccionados por cada método.

    Si algún método no produce genes significativos, usa los Top-N por menor q.
    Devuelve: (matrices, gene_lists, basic_ref_gene)
    """
    sig, tables, basic_ref_gene = detect_significant_genes_by_method(
        controles, muestras, adv_result, alpha=alpha
    )

    df_total = build_total_dataframe(controles, muestras)

    # Normalizados por método
    adv_norm = adv_result.df_norm.copy()
    gm_norm = _normalize_by_global_mean(df_total)
    ref_norm, _ = _normalize_by_single_reference(df_total, basic_ref_gene)

    matrices: dict[str, pd.DataFrame] = {}
    gene_lists: dict[str, list[str]] = {}

    for method, df_norm in (
        ("advanced", adv_norm),
        ("global_mean", gm_norm),
        ("refgene", ref_norm),
    ):
        genes = list(sig.get(method) or [])
        if not genes:
            table = tables.get(method)
            if table is not None and not table.empty:
                genes = table.sort_values(["q", "t_abs"], ascending=[True, False])["gene"].head(top_n_fallback).tolist()
        gene_lists[method] = genes
        matrices[method] = build_heatmap_matrix(df_norm, genes)

    return matrices, gene_lists, basic_ref_gene


def save_gene_sets(
    gene_lists: dict[str, list[str]],
    *,
    output_dir: str = "resultados/genes_significativos",
    basic_ref_gene: str | None = None,
) -> dict[str, str]:
    """Guarda los conjuntos de genes por método como CSV (una columna 'gene').

    Devuelve dict método -> ruta del archivo.
    """
    import os

    os.makedirs(output_dir, exist_ok=True)
    paths: dict[str, str] = {}
    for method, genes in gene_lists.items():
        fname = f"genes_significativos_{method}.csv"
        if method == "refgene" and basic_ref_gene:
            fname = f"genes_significativos_refgene_{basic_ref_gene}.csv"
        path = os.path.join(output_dir, fname)
        pd.DataFrame({"gene": genes}).to_csv(path, index=False)
        paths[method] = path
    return paths


def build_expression_datasets(
    controles: pd.DataFrame,
    muestras: pd.DataFrame,
    adv_result: AdvancedNormalizationResult,
    *,
    include_reference_ct: bool = True,
) -> tuple[dict[str, pd.DataFrame], str]:
    """Devuelve DataFrames normalizados (formato largo) para cada método."""

    gm_norm, ref_norm, basic_ref_gene = compute_basic_normalizations(controles, muestras, adv_result)

    datasets: dict[str, pd.DataFrame] = {
        "advanced": adv_result.df_norm.copy(),
        "global_mean": gm_norm.copy(),
        "refgene": ref_norm.copy(),
    }

    if not include_reference_ct:
        for df in datasets.values():
            for col in ("ref_ct",):
                if col in df.columns:
                    df.drop(columns=col, inplace=True)

    return datasets, basic_ref_gene


def build_combined_ddct_fc_tables(
    adv_result: AdvancedNormalizationResult,
    fold_change_result: FoldChangeResult,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Construye tablas consolidadas (ΔΔCt y Fold Change) integrando los tres métodos."""

    consolidated = fold_change_result.consolidated.copy()
    # Log2FC para métodos básicos
    consolidated["log2fc_promedio"] = -consolidated["delta_delta_ct_promedio"].astype(float)
    consolidated["log2fc_gen_ref"] = -consolidated["delta_delta_ct_gen_ref"].astype(float)

    adv_norm = adv_result.df_norm.copy()
    ctrl_mask = adv_norm["tipo"] == "Control"
    ctrl_means = adv_norm.loc[ctrl_mask].groupby("target")["log2_rel_expr"].mean()
    case_means = adv_norm.loc[~ctrl_mask].groupby("target")["log2_rel_expr"].mean()

    all_targets = sorted(set(consolidated["target"].astype(str)) | set(ctrl_means.index.astype(str)) | set(case_means.index.astype(str)))
    ctrl_means = ctrl_means.reindex(all_targets)
    case_means = case_means.reindex(all_targets)

    adv_log2fc = (case_means - ctrl_means).reindex(all_targets)
    adv_ddct = (-adv_log2fc).reindex(all_targets)
    adv_fc = np.power(2.0, -adv_ddct)

    advanced_df = pd.DataFrame(
        {
            "target": all_targets,
            "delta_delta_ct_advanced": adv_ddct.to_numpy(dtype=float, na_value=np.nan),
            "fold_change_advanced": adv_fc.to_numpy(dtype=float, na_value=np.nan),
            "log2fc_advanced": adv_log2fc.to_numpy(dtype=float, na_value=np.nan),
        }
    )

    merged = consolidated.reset_index(drop=True).merge(advanced_df, on="target", how="outer")

    ddct_table = merged[[
        "target",
        "delta_delta_ct_advanced",
        "delta_delta_ct_promedio",
        "delta_delta_ct_gen_ref",
    ]].copy()

    fc_table = merged[[
        "target",
        "fold_change_advanced",
        "fold_change_promedio",
        "fold_change_gen_ref",
        "log2fc_advanced",
        "log2fc_promedio",
        "log2fc_gen_ref",
    ]].copy()

    return ddct_table, fc_table


__all__ += [
    "compute_basic_normalizations",
    "detect_significant_genes_by_method",
    "build_heatmaps_by_method",
    "save_gene_sets",
    "build_expression_datasets",
    "build_combined_ddct_fc_tables",
]
