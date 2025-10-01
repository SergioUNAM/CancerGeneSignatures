from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from typing import Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind


@dataclass(frozen=True)
class ReferenceSelectionResult:
    """Information about the chosen reference set and search space."""

    references: Tuple[str, ...]
    score: float
    candidate_scores: pd.DataFrame
    coverage_warnings: List[str]


@dataclass(frozen=True)
class DifferentialExpressionResult:
    """Differential expression statistics after normalization."""

    stats: pd.DataFrame
    alpha: float
    significant_genes: List[str]
    bootstrap_freq: pd.Series
    false_positive_rate: float


@dataclass(frozen=True)
class AdvancedNormalizationResult:
    """Aggregate output of the advanced normalization pipeline."""

    df_base: pd.DataFrame
    df_norm: pd.DataFrame
    df_heatmap: pd.DataFrame
    reference_result: ReferenceSelectionResult
    differential: DifferentialExpressionResult

    def summary(self) -> dict:
        diff = self.differential
        ref = self.reference_result
        return {
            "refs_elegidas": list(ref.references),
            "score_estabilidad": round(ref.score, 4),
            f"n_genes_signif_q<{diff.alpha:.2f}": int(len(diff.significant_genes)),
            "fpr_empirica_aprox": (round(diff.false_positive_rate, 4)
                                    if not np.isnan(diff.false_positive_rate)
                                    else np.nan),
        }


class AdvancedNormalizationError(Exception):
    """Error wrapper for issues during the advanced normalization pipeline."""


REQUIRED_COLUMNS = {"test", "target", "ct", "tipo"}


def _validate_input(df: pd.DataFrame) -> pd.DataFrame:
    if not REQUIRED_COLUMNS.issubset(df.columns):
        faltan = REQUIRED_COLUMNS.difference(df.columns)
        raise AdvancedNormalizationError(
            f"Faltan columnas obligatorias para la normalización avanzada: {sorted(faltan)}"
        )
    out = df.copy()
    out["test"] = out["test"].astype(str)
    out["target"] = out["target"].astype(str)
    out["tipo"] = out["tipo"].astype(str)
    out["ct"] = pd.to_numeric(out["ct"], errors="coerce")
    if out["ct"].isna().all():
        raise AdvancedNormalizationError("Todos los valores Ct son NaN después de la conversión numérica.")
    return out


def stability_score_for_refs(
    df: pd.DataFrame,
    refs: Sequence[str],
    *,
    weight_intra: float = 0.7,
    weight_inter: float = 0.3,
) -> float:
    """Compute stability score: weighted intra-group SD and inter-group mean difference."""

    refs = tuple(dict.fromkeys(refs))  # unique sequence while preserving order
    if not refs:
        return float("inf")

    subset = df[df["target"].isin(refs)].copy()
    if subset.empty:
        return float("inf")

    # Mean Ct per reference gene and test
    per_gene = (
        subset.groupby(["test", "target"], as_index=False)["ct"]
        .mean()
        .pivot(index="test", columns="target", values="ct")
    )

    # Average Ct across references (row-wise mean)
    ref_ct = per_gene.mean(axis=1)
    counts = per_gene.notna().sum(axis=1)
    valid_ref_ct = ref_ct[counts > 0]
    if valid_ref_ct.empty:
        return float("inf")

    tmp = df.merge(valid_ref_ct.rename("ref_ct"), on="test", how="left")
    ref_rows = tmp[tmp["target"].isin(refs)].copy()
    if ref_rows["ref_ct"].isna().all():
        return float("inf")

    # Intra-group variability (average of SD over groups)
    by_tipo_test = (
        ref_rows.groupby(["tipo", "test"])["ref_ct"].mean()
        .groupby(level=0)
    )
    s_intra = by_tipo_test.std(ddof=1)
    s_intra_mean = float(s_intra.mean()) if s_intra.notna().any() else float("inf")

    # Inter-group difference (max difference of group means)
    m_by_grp = ref_rows.groupby("tipo")["ref_ct"].mean()
    inter_diff = float(m_by_grp.max() - m_by_grp.min()) if len(m_by_grp) >= 2 else float("inf")

    return weight_intra * s_intra_mean + weight_inter * inter_diff


def _compute_reference_ct(df: pd.DataFrame, refs: Sequence[str]) -> Tuple[pd.Series, List[str]]:
    subset = df[df["target"].isin(refs)].copy()
    if subset.empty:
        return pd.Series(dtype=float), ["No se encontraron Ct para las referencias seleccionadas."]

    per_gene = (
        subset.groupby(["test", "target"], as_index=False)["ct"]
        .mean()
        .pivot(index="test", columns="target", values="ct")
    )
    ref_ct = per_gene.mean(axis=1)
    counts = per_gene.notna().sum(axis=1)
    missing = counts[counts < len(refs)]
    warnings: List[str] = []
    if not missing.empty:
        total = len(counts)
        warnings.append(
            "Algunas muestras carecen de todos los genes de referencia ("
            f"{len(missing)}/{total}). Se calculará el promedio con los Ct disponibles."
        )
    return ref_ct, warnings


def normalize_with_references(df: pd.DataFrame, refs: Sequence[str]) -> Tuple[pd.DataFrame, List[str]]:
    ref_ct, warnings = _compute_reference_ct(df, refs)
    if ref_ct.empty:
        raise AdvancedNormalizationError("No hay Ct válidos para las referencias seleccionadas.")
    out = df.merge(ref_ct.rename("ref_ct"), on="test", how="left")
    out["delta_ct"] = out["ct"] - out["ref_ct"]
    out["log2_rel_expr"] = -out["delta_ct"]
    return out, warnings


def evaluate_differential_expression(
    df_norm: pd.DataFrame,
    *,
    alpha: float = 0.05,
) -> DifferentialExpressionResult:
    rows: List[tuple] = []
    for gene, sub in df_norm.groupby("target", sort=False):
        ctrl = sub.loc[sub["tipo"] == "Control", "log2_rel_expr"].dropna()
        case = sub.loc[sub["tipo"] != "Control", "log2_rel_expr"].dropna()
        if len(ctrl) > 1 and len(case) > 1:
            t_stat, p_val = ttest_ind(ctrl, case, equal_var=False)
            pooled_sd = np.sqrt((ctrl.var(ddof=1) + case.var(ddof=1)) / 2.0)
            cohen_d = (case.mean() - ctrl.mean()) / pooled_sd if pooled_sd > 0 else np.nan
            rows.append(
                (
                    gene,
                    float(abs(t_stat)),
                    float(p_val),
                    float(cohen_d),
                    float(ctrl.mean()),
                    float(case.mean()),
                    float(ctrl.std(ddof=1) if len(ctrl) > 1 else np.nan),
                    float(case.std(ddof=1) if len(case) > 1 else np.nan),
                )
            )
    columns = [
        "gene",
        "t_abs",
        "p",
        "cohen_d",
        "mean_ctrl",
        "mean_case",
        "sd_ctrl",
        "sd_case",
    ]
    stats = pd.DataFrame(rows, columns=columns)
    if not stats.empty:
        stats = stats.sort_values("p")
        stats["q"] = _fdr_bh(stats["p"])  # type: ignore[arg-type]
        stats = stats.sort_values(["q", "t_abs"], ascending=[True, False]).reset_index(drop=True)
    else:
        stats["q"] = pd.Series(dtype=float)
    significant = stats.loc[stats["q"] < alpha, "gene"].tolist() if not stats.empty else []
    return DifferentialExpressionResult(
        stats=stats,
        alpha=alpha,
        significant_genes=significant,
        bootstrap_freq=pd.Series(dtype=float, name="bootstrap_freq"),
        false_positive_rate=float("nan"),
    )


def bootstrap_significance(
    df_base: pd.DataFrame,
    refs: Sequence[str],
    *,
    iterations: int = 300,
    alpha: float = 0.05,
    random_state: Optional[np.random.Generator] = None,
) -> pd.Series:
    if iterations <= 0:
        return pd.Series(dtype=float, name="bootstrap_freq")
    rng = random_state or np.random.default_rng()
    tests = df_base["test"].unique()
    genes = df_base["target"].unique()
    counts = dict.fromkeys(genes, 0)

    for _ in range(iterations):
        sampled_tests = rng.choice(tests, size=len(tests), replace=True)
        sampled_rows = [df_base[df_base["test"] == t] for t in sampled_tests]
        sub = pd.concat(sampled_rows, ignore_index=True)
        sub_norm, _ = normalize_with_references(sub, refs)
        diff = evaluate_differential_expression(sub_norm, alpha=alpha)
        sig = set(diff.significant_genes)
        for g in sig:
            counts[g] = counts.get(g, 0) + 1

    freq = pd.Series(counts, name="bootstrap_freq", dtype=float) / float(iterations)
    return freq.sort_values(ascending=False)


def permutation_false_positive_rate(
    df_base: pd.DataFrame,
    refs: Sequence[str],
    *,
    iterations: int = 200,
    alpha: float = 0.05,
    random_state: Optional[np.random.Generator] = None,
) -> float:
    if iterations <= 0:
        return float("nan")
    rng = random_state or np.random.default_rng()
    tests = df_base["test"].unique()
    if len(tests) == 0:
        return float("nan")
    fprs: List[float] = []
    tipo_por_test = df_base.groupby("test")["tipo"].agg(lambda x: x.iloc[0])

    for _ in range(iterations):
        permuted = pd.Series(rng.permutation(tipo_por_test.values), index=tipo_por_test.index)
        sub = df_base.copy()
        sub["tipo_perm"] = sub["test"].map(permuted)
        sub_norm, _ = normalize_with_references(sub, refs)

        pvals: List[float] = []
        for _, gsub in sub_norm.groupby("target", sort=False):
            ctrl = gsub.loc[gsub["tipo_perm"] == "Control", "log2_rel_expr"].dropna()
            case = gsub.loc[gsub["tipo_perm"] != "Control", "log2_rel_expr"].dropna()
            if len(ctrl) > 1 and len(case) > 1:
                _, p_val = ttest_ind(ctrl, case, equal_var=False)
                pvals.append(float(p_val))
        if pvals:
            q_perm = _fdr_bh(pvals)
            fprs.append(float((q_perm < alpha).mean()))

    return float(np.nanmean(fprs)) if fprs else float("nan")


def select_reference_set(
    df_base: pd.DataFrame,
    *,
    n_candidates: int = 20,
    k_refs: int = 2,
) -> ReferenceSelectionResult:
    sd_por_gen = (
        df_base.groupby("target")["ct"].std()
        .replace([np.inf, -np.inf], np.nan)
        .dropna()
        .sort_values()
    )
    candidatos = sd_por_gen.index.tolist()[:n_candidates]
    if len(candidatos) < k_refs:
        raise AdvancedNormalizationError(
            f"No hay suficientes candidatos ({len(candidatos)}) para K={k_refs}."
        )

    results = []
    mejor_refs: Optional[Tuple[str, ...]] = None
    mejor_score = float("inf")
    for refs in combinations(candidatos, k_refs):
        score = stability_score_for_refs(df_base, refs)
        results.append({"refs": refs, "score": score})
        if score < mejor_score:
            mejor_score = score
            mejor_refs = tuple(refs)

    scores_df = pd.DataFrame(results).sort_values("score").reset_index(drop=True)
    assert mejor_refs is not None
    return ReferenceSelectionResult(
        references=mejor_refs,
        score=float(mejor_score),
        candidate_scores=scores_df,
        coverage_warnings=[],
    )


def build_heatmap_matrix(df_norm: pd.DataFrame, genes_for_heatmap: Sequence[str]) -> pd.DataFrame:
    if not genes_for_heatmap:
        return pd.DataFrame()
    pivot = df_norm[df_norm["target"].isin(genes_for_heatmap)].pivot_table(
        index="target",
        columns="test",
        values="log2_rel_expr",
        aggfunc="mean",
    )
    return pivot.sort_index()


def run_advanced_normalization(
    df_total: pd.DataFrame,
    *,
    alpha: float = 0.05,
    n_candidates: int = 20,
    k_refs: int = 2,
    bootstrap_iter: int = 300,
    permutation_iter: int = 200,
    random_seed: Optional[int] = 123,
) -> AdvancedNormalizationResult:
    df_base = _validate_input(df_total)

    ref_selection = select_reference_set(df_base, n_candidates=n_candidates, k_refs=k_refs)
    df_norm, warn_cov = normalize_with_references(df_base, ref_selection.references)
    ref_selection.coverage_warnings.extend(warn_cov)

    diff = evaluate_differential_expression(df_norm, alpha=alpha)

    rng = np.random.default_rng(random_seed) if random_seed is not None else None
    bootstrap_freq = bootstrap_significance(
        df_base,
        ref_selection.references,
        iterations=bootstrap_iter,
        alpha=alpha,
        random_state=rng,
    )

    stats = diff.stats.copy()
    if not stats.empty:
        stats = stats.merge(
            bootstrap_freq.rename("bootstrap_freq"),
            left_on="gene",
            right_index=True,
            how="left",
        ).fillna({"bootstrap_freq": 0.0})
    else:
        stats["bootstrap_freq"] = pd.Series(dtype=float)

    sig_genes = stats.loc[stats["q"] < diff.alpha, "gene"].tolist() if not stats.empty else diff.significant_genes

    fpr_emp = permutation_false_positive_rate(
        df_base,
        ref_selection.references,
        iterations=permutation_iter,
        alpha=alpha,
        random_state=rng,
    )

    diff = DifferentialExpressionResult(
        stats=stats,
        alpha=diff.alpha,
        significant_genes=sig_genes,
        bootstrap_freq=bootstrap_freq,
        false_positive_rate=fpr_emp,
    )

    heat_genes = diff.significant_genes if diff.significant_genes else diff.stats.head(20)["gene"].tolist()
    df_heat = build_heatmap_matrix(df_norm, heat_genes)

    return AdvancedNormalizationResult(
        df_base=df_base,
        df_norm=df_norm,
        df_heatmap=df_heat,
        reference_result=ref_selection,
        differential=diff,
    )


__all__ = [
    "AdvancedNormalizationError",
    "AdvancedNormalizationResult",
    "DifferentialExpressionResult",
    "ReferenceSelectionResult",
    "build_heatmap_matrix",
    "evaluate_differential_expression",
    "normalize_with_references",
    "run_advanced_normalization",
    "select_reference_set",
    "stability_score_for_refs",
]
def _fdr_bh(pvals: Sequence[float]) -> np.ndarray:
    """Benjamini-Hochberg FDR correction."""

    pvals_arr = np.asarray(pvals, dtype=float)
    if pvals_arr.size == 0:
        return np.array([], dtype=float)
    order = np.argsort(pvals_arr)
    ranked = np.empty_like(pvals_arr)
    n = float(pvals_arr.size)
    adjusted = (n / (np.arange(pvals_arr.size, dtype=float) + 1.0)) * pvals_arr[order]
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    ranked[order] = adjusted
    return ranked
