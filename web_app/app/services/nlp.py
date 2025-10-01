from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional

import pandas as pd

try:
    from app.core.bibliography import filter_bibliography_by_cancer  # type: ignore
except Exception:  # pragma: no cover - dependencia opcional
    filter_bibliography_by_cancer = None  # type: ignore

from app.integrations.google_nlp import GoogleNLPClient, aggregate_insights


@dataclass(frozen=True)
class NLPCorpusInfo:
    filtered: pd.DataFrame
    gene_filtered: pd.DataFrame
    title_column: str
    abstract_column: Optional[str]
    gene_column: Optional[str]
    total_all: int
    total_gene: int
    total_filtered: int


def prepare_corpus(
    df: pd.DataFrame,
    *,
    cancer_label: str,
    gene: Optional[str],
    apply_filter: bool,
) -> NLPCorpusInfo:
    """Normaliza el DataFrame de bibliografía y calcula métricas básicas."""

    if df is None or df.empty:
        raise ValueError("El DataFrame de bibliografía está vacío.")

    data = df.copy()
    title_col = _select_title_column(data)
    if title_col is None:
        raise ValueError("El CSV debe incluir una columna 'Title' o 'article_title'.")
    abstract_col = _select_abstract_column(data)

    text_col = data[title_col].astype(str) + " " + (data[abstract_col].astype(str) if abstract_col else "")
    data["_text"] = text_col
    total_all = int((data["_text"].str.strip() != "").sum())

    gene_col = _select_gene_column(data)
    gene_filtered = data
    if gene and gene_col:
        mask = gene_filtered[gene_col].astype(str).str.strip().str.lower() == gene.strip().lower()
        gene_filtered = gene_filtered[mask]
    total_gene = int((gene_filtered["_text"].str.strip() != "").sum())

    filtered = gene_filtered
    if apply_filter and filter_bibliography_by_cancer is not None:
        try:
            filtered = filter_bibliography_by_cancer(filtered, cancer_label)
        except Exception:
            pass
    total_filtered = int((filtered["_text"].str.strip() != "").sum())

    gene_filtered = gene_filtered.drop(columns=["_text"]) if "_text" in gene_filtered else gene_filtered
    filtered = filtered.drop(columns=["_text"]) if "_text" in filtered else filtered

    return NLPCorpusInfo(
        filtered=filtered,
        gene_filtered=gene_filtered,
        title_column=title_col,
        abstract_column=abstract_col,
        gene_column=gene_col,
        total_all=total_all,
        total_gene=total_gene,
        total_filtered=total_filtered,
    )


def extract_texts(
    df: pd.DataFrame,
    *,
    title_column: str,
    abstract_column: Optional[str],
    limit: int,
) -> List[str]:
    """Construye la lista de textos combinando título y abstract."""

    if df.empty:
        return []
    texts = df[title_column].astype(str)
    if abstract_column:
        texts = texts + " " + df[abstract_column].astype(str)
    result = [s.strip() for s in texts.tolist() if isinstance(s, str) and s.strip()]
    return result[: max(1, int(limit))]


def analyze_texts(
    texts: List[str],
    *,
    api_key: Optional[str],
    language: Optional[str],
) -> dict:
    """Ejecuta Google NLP sobre los textos preparados y retorna agregaciones."""

    if not texts:
        return {}
    client = GoogleNLPClient(api_key=api_key or None, default_language=language, sleep_between=0.0)
    return aggregate_insights(
        texts,
        client,
        do_entities=True,
        do_entity_sentiment=True,
        do_sentiment=True,
        do_categories=True,
        language=language,
        max_chars_per_doc=8000,
    )


def _select_title_column(df: pd.DataFrame) -> Optional[str]:
    for col in ("Title", "article_title"):
        if col in df.columns:
            return col
    return None


def _select_abstract_column(df: pd.DataFrame) -> Optional[str]:
    for col in ("Abstract", "abstract"):
        if col in df.columns:
            return col
    return None


def _select_gene_column(df: pd.DataFrame) -> Optional[str]:
    for col in ("Gene", "gene", "Symbol", "symbol", "target"):
        if col in df.columns:
            return col
    return None


__all__ = [
    "NLPCorpusInfo",
    "prepare_corpus",
    "extract_texts",
    "analyze_texts",
]
