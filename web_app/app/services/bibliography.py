from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import pandas as pd

from src.core.bibliography import search_pubmed_by_genes


@dataclass(frozen=True)
class PubMedRequest:
    context: str
    max_per_gene: int
    email: str
    api_key: Optional[str] = None


def fetch_pubmed_articles(genes_df: pd.DataFrame, request: PubMedRequest) -> pd.DataFrame:
    """Consulta PubMed para un conjunto de genes, devolviendo el DataFrame resultante."""

    if genes_df.empty:
        return pd.DataFrame()

    return search_pubmed_by_genes(
        genes_df,
        symbol_col="target",
        ensembl_col="ensembl_id",
        selected_context=request.context,
        max_per_gene=int(request.max_per_gene),
        progress=None,
        logger=None,
        email=request.email,
        api_key=request.api_key,
    )


def merge_expression_levels(pubmed_df: pd.DataFrame, genes_df: pd.DataFrame) -> pd.DataFrame:
    """Anexa la columna 'nivel_expresion' a los resultados de PubMed si está disponible."""

    if pubmed_df.empty:
        return pubmed_df

    merge_cols = ["target", "Gene"]
    if "Gene" in pubmed_df.columns:
        left_key = "Gene"
    elif "target" in pubmed_df.columns:
        left_key = "target"
    else:
        return pubmed_df

    genes_clean = genes_df.drop_duplicates(subset=["target"]).copy()
    genes_clean.rename(columns={"target": left_key}, inplace=True)
    return pd.merge(pubmed_df, genes_clean[[left_key, "nivel_expresion"]], on=left_key, how="left")


__all__ = [
    "PubMedRequest",
    "fetch_pubmed_articles",
    "merge_expression_levels",
]
