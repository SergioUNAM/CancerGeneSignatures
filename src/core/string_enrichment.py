from __future__ import annotations

import os
from typing import Iterable, List, Optional

import pandas as pd
import requests


# Default to explicit version used in notebooks for reproducibility
STRING_ENRICH_URL = "https://version-12-0.string-db.org/api/json/enrichment"


def _join_identifiers(genes: Iterable[str]) -> str:
    # STRING expects identifiers separated by %0d (URL-encoded CR)
    cleaned = [str(g).strip() for g in genes if isinstance(g, str) and str(g).strip()]
    return "%0d".join(cleaned)


def run_string_enrichment(
    genes: Iterable[str],
    species: int = 9606,
    sources: Optional[List[str]] = None,
    caller_identity: Optional[str] = None,
    timeout: int = 20,
    base_url: Optional[str] = None,
) -> pd.DataFrame:
    """
    Query STRING's enrichment API for a list of gene symbols.

    Parameters
    - genes: iterable of gene symbols (HGNC for human).
    - species: NCBI Taxon ID (9606 = human).
    - sources: Optional list of sources to include, e.g. ["GO", "KEGG"].
    - caller_identity: Optional caller name for STRING logs.
    - timeout: request timeout in seconds.

    Returns
    - DataFrame with enrichment results. Empty if none or on failure.
    """

    ids = _join_identifiers(genes)
    if not ids:
        return pd.DataFrame()

    params = {
        "identifiers": ids,
        "species": int(species),
        "caller_identity": caller_identity or os.getenv("CGS_CALLER", "CancerGeneSignatures"),
    }
    # Note: STRING supports "source" to filter categories, e.g., GO:BP, GO:MF, GO:CC, KEGG, Reactome
    if sources:
        # Accept common shorthands and raw values
        # Examples: ["GO", "KEGG"], ["GO:BP"], ["KEGG"]
        params["source"] = "%0d".join(sources)

    url = base_url or STRING_ENRICH_URL
    data = None
    try:
        # Prefer POST per notebook example
        resp = requests.post(url, data=params, timeout=timeout)
        resp.raise_for_status()
        data = resp.json()
    except Exception:
        # Fallback to GET if POST fails
        try:
            resp = requests.get(url, params=params, timeout=timeout)
            resp.raise_for_status()
            data = resp.json()
        except Exception:
            return pd.DataFrame()
    if not isinstance(data, list) or not data:
        return pd.DataFrame()

    # Normalize JSON to a tidy DataFrame
    df = pd.json_normalize(data, sep="_")

    # Expected fields (defensive selection)
    cols_map = {
        "category": "category",
        "description": "term",
        "term": "term_id",
        "p_value": "p_value",
        "fdr": "fdr",
        "inputGenes": "input_genes",
        "number_of_genes": "term_genes",
        "preferredNames": "preferred_names",
    }
    out = pd.DataFrame()
    for src, dst in cols_map.items():
        if src in df.columns:
            out[dst] = df[src]

    # Split preferred names for readability
    if "preferred_names" in out.columns:
        out["preferred_names"] = out["preferred_names"].apply(
            lambda v: ", ".join(v) if isinstance(v, list) else (v or "")
        )

    # Sort by FDR then p-value
    if "fdr" in out.columns:
        out = out.sort_values(["fdr", "p_value"], na_position="last")

    return out.reset_index(drop=True)


def filter_enrichment(
    df: pd.DataFrame,
    include_categories: Optional[List[str]] = None,
    max_fdr: float = 0.05,
    min_term_genes: int = 2,
    top_n: Optional[int] = 25,
) -> pd.DataFrame:
    """Apply simple, useful filters to an enrichment DataFrame."""
    if df is None or df.empty:
        return pd.DataFrame()
    filt = df.copy()
    if include_categories:
        # Match if the category starts with any desired prefix (e.g., GO, GO:BP, KEGG)
        incl = tuple(include_categories)
        filt = filt[filt["category"].astype(str).str.upper().str.startswith(tuple(c.upper() for c in incl))]
    if "fdr" in filt.columns:
        filt = filt[filt["fdr"] <= float(max_fdr)]
    if "term_genes" in filt.columns:
        filt = filt[filt["term_genes"] >= int(min_term_genes)]
    if top_n is not None and top_n > 0:
        filt = filt.head(int(top_n))
    return filt.reset_index(drop=True)


def enrich_by_levels(
    df_expr: pd.DataFrame,
    symbol_col: str = "target",
    level_col: str = "nivel_expresion",
    levels: Optional[List[str]] = None,
    species: int = 9606,
    sources: Optional[List[str]] = None,
    caller_identity: Optional[str] = None,
    timeout: int = 20,
) -> dict:
    """
    Run STRING enrichment per expression level and return a dict with:
      - combined: concatenated DataFrame with a level column
      - by_level: mapping level -> DataFrame
    """
    if df_expr is None or df_expr.empty:
        return {"combined": pd.DataFrame(), "by_level": {}}

    levels = levels or sorted(df_expr[level_col].dropna().astype(str).unique().tolist())
    by_level = {}
    frames = []
    for lvl in levels:
        sub = df_expr[df_expr[level_col] == lvl]
        genes = sub[symbol_col].dropna().astype(str).unique().tolist()
        if not genes:
            by_level[lvl] = pd.DataFrame()
            continue
        df_enr = run_string_enrichment(
            genes,
            species=species,
            sources=sources,
            caller_identity=caller_identity,
            timeout=timeout,
        )
        if not df_enr.empty:
            df_enr = df_enr.copy()
            df_enr[level_col] = lvl
        by_level[lvl] = df_enr
        frames.append(df_enr)
    combined = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    return {"combined": combined, "by_level": by_level}


def dfs_to_excel_bytes(dfs: List[pd.DataFrame], sheet_names: List[str]) -> bytes:
    """Export multiple DataFrames to an Excel bytes object."""
    import io
    buf = io.BytesIO()
    with pd.ExcelWriter(buf, engine="openpyxl") as writer:
        for df, name in zip(dfs, sheet_names):
            (df if isinstance(df, pd.DataFrame) else pd.DataFrame()).to_excel(writer, index=False, sheet_name=name[:31] or "Sheet")
    buf.seek(0)
    return buf.getvalue()
