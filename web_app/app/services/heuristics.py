from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from src.core.bibliography import (
    interpret_gene_relations,
    summarize_relations_by_gene,
)


@dataclass(frozen=True)
class SankeyData:
    nodes: List[str]
    sources: List[int]
    targets: List[int]
    values: List[float]


def compute_heuristic_summary(focused_df: pd.DataFrame) -> pd.DataFrame:
    """Devuelve el resumen heurístico por gen a partir de la bibliografía clasificada."""

    if focused_df is None or focused_df.empty:
        return pd.DataFrame()
    interpreted = interpret_gene_relations(focused_df)
    if interpreted is None or interpreted.empty:
        return pd.DataFrame()
    summary = summarize_relations_by_gene(interpreted)
    return summary if isinstance(summary, pd.DataFrame) else pd.DataFrame()


def build_heatmap_data(summary: pd.DataFrame, top: int = 30) -> Optional[pd.DataFrame]:
    """Construye la matriz normalizada para el heatmap genes × funciones."""

    if summary is None or summary.empty:
        return None
    func_cols = [
        c for c in summary.columns
        if c.endswith("_score") and not c.startswith(("upregulated", "downregulated", "prognosis_"))
    ]
    if not func_cols:
        return None
    sm = summary.copy()
    sm["_total_score"] = sm[func_cols].sum(axis=1)
    top_sm = sm.sort_values("_total_score", ascending=False).head(top)
    data = top_sm[func_cols]
    norm = (data - data.min()) / (data.max() - data.min() + 1e-9)
    norm.index = top_sm["Gene"].tolist()
    norm.columns = [c.replace("_score", "") for c in func_cols]
    return norm


def build_sankey_data(summary: pd.DataFrame) -> Optional[SankeyData]:
    """Prepara los nodos y enlaces para el Sankey de expresión → pronóstico."""

    if summary is None or summary.empty:
        return None
    bad = summary.get("prognosis_bad_score")
    good = summary.get("prognosis_good_score")
    if bad is None or good is None:
        return None
    if not hasattr(bad, "values") or not hasattr(good, "values"):
        return None
    dom = np.where(
        bad.values > good.values,
        "prognosis_bad",
        np.where(good.values > 0, "prognosis_good", "prognosis_uncertain"),
    )
    expr_effect = summary.get("net_expression_effect")
    mid = expr_effect.astype(str).tolist() if isinstance(expr_effect, pd.Series) else ["uncertain"] * len(summary)
    left = summary["Gene"].astype(str).tolist()
    right = dom.tolist()
    nodes = list(dict.fromkeys(left + mid + right))
    idx = {n: i for i, n in enumerate(nodes)}
    sources = [idx[l] for l in left] + [idx[m] for m in mid]
    targets = [idx[m] for m in mid] + [idx[r] for r in right]
    values = [1.0] * len(sources)
    return SankeyData(nodes=nodes, sources=sources, targets=targets, values=values)


def merge_with_expression(summary: pd.DataFrame, df_expr: Optional[pd.DataFrame], exclude_stable: bool) -> pd.DataFrame:
    """Fusiona el resumen heurístico con niveles de expresión."""

    if summary is None or summary.empty or df_expr is None or df_expr.empty:
        return pd.DataFrame()
    lev = df_expr[["target", "nivel_expresion", "fold_change"]].drop_duplicates("target").rename(columns={"target": "Gene"})
    merged = pd.merge(summary, lev, on="Gene", how="left")
    if exclude_stable:
        merged = merged[merged["nivel_expresion"] != "estable"]
    return merged


def build_function_long(merged: pd.DataFrame) -> pd.DataFrame:
    """Devuelve los datos en formato largo (Gene, función, score, nivel, fold_change)."""

    if merged is None or merged.empty:
        return pd.DataFrame()
    func_score_cols = [
        c for c in merged.columns
        if c.endswith("_score") and not c.startswith(("upregulated", "downregulated", "prognosis_"))
    ]
    if not func_score_cols:
        return pd.DataFrame()
    records: List[Dict[str, object]] = []
    for _, row in merged.iterrows():
        gene = row.get("Gene")
        lvl = row.get("nivel_expresion")
        fc = row.get("fold_change")
        for fscore in func_score_cols:
            score = row.get(fscore)
            if score and float(score) > 0:
                records.append({
                    "Gene": gene,
                    "funcion": fscore.replace("_score", ""),
                    "score": float(score),
                    "flag": 1 if float(score) >= 1.2 else 0,
                    "nivel_expresion": lvl,
                    "fold_change": fc,
                })
    return pd.DataFrame(records)


def compute_function_counts(long_df: pd.DataFrame) -> pd.DataFrame:
    """Calcula la agregación función × nivel con sumatoria de scores."""

    if long_df is None or long_df.empty:
        return pd.DataFrame()
    agg = (
        long_df.groupby(["nivel_expresion", "funcion"])
        ["score"].sum()
        .reset_index()
        .rename(columns={"score": "total_score"})
    )
    return agg


def map_functions_to_hallmarks(
    summary: pd.DataFrame,
    signatures_df: Optional[pd.DataFrame],
    cancer_type: str,
) -> pd.DataFrame:
    """Relaciona funciones heurísticas con Hallmarks a partir de las firmas generadas."""

    if summary is None or summary.empty:
        return pd.DataFrame()
    if signatures_df is None or signatures_df.empty:
        return pd.DataFrame()
    base = signatures_df
    if "cancer_type" in base.columns:
        base = base[base["cancer_type"] == cancer_type]
    records: List[Dict[str, object]] = []
    for _, row in base.iterrows():
        for col in row.index:
            if col.startswith("hallmark_") and col.endswith("_genes"):
                hallmark = col[len("hallmark_"):-len("_genes")]
                genes = row[col]
                if isinstance(genes, list):
                    for gene in genes:
                        records.append({"Gene": gene, "hallmark": hallmark})
    if not records:
        return pd.DataFrame()
    hall_df = pd.DataFrame(records).drop_duplicates()
    return hall_df


__all__ = [
    "SankeyData",
    "compute_heuristic_summary",
    "build_heatmap_data",
    "build_sankey_data",
    "merge_with_expression",
    "build_function_long",
    "compute_function_counts",
    "map_functions_to_hallmarks",
]
