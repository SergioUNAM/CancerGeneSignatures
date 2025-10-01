from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache
from typing import List, Tuple

import requests
import pandas as pd


BASE_URL = "https://rest.ensembl.org"
HEADERS = {
    "Content-Type": "application/json",
    "Accept": "application/json",
    "User-Agent": "CancerGeneSignatures/0.1 (streamlit app)"
}


def _with_content_type(url: str) -> str:
    return url + ("&" if "?" in url else "?") + "content-type=application/json"


def _get(session: requests.Session, url: str, timeout: float = 10.0):
    try:
        return session.get(_with_content_type(url), headers=HEADERS, timeout=timeout)
    except Exception:
        return None


def _try_xrefs_symbol(session: requests.Session, gene: str) -> Tuple[str, str]:
    r = _get(session, f"{BASE_URL}/xrefs/symbol/homo_sapiens/{gene}?content-type=application/json")
    if r and r.ok:
        data = r.json() or []
        # Preferir entradas de tipo 'gene'
        for item in data:
            if item.get('type') == 'gene':
                ensembl_id = item.get('id') or ''
                desc = item.get('description') or ''
                return ensembl_id, desc
        # Fallback a la primera entrada
        if data:
            ensembl_id = data[0].get('id') or ''
            desc = data[0].get('description') or ''
            return ensembl_id, desc
    return '', ''


def _try_lookup_id_desc(session: requests.Session, ensembl_id: str) -> str:
    if not ensembl_id:
        return ''
    r = _get(session, f"{BASE_URL}/lookup/id/{ensembl_id}?expand=1")
    if r and r.ok:
        data = r.json() or {}
        return data.get('description') or ''
    return ''


def _try_lookup_symbol(session: requests.Session, gene: str) -> Tuple[str, str]:
    r = _get(session, f"{BASE_URL}/lookup/symbol/homo_sapiens/{gene}?expand=1")
    if r and r.ok:
        data = r.json() or {}
        ensembl_id = data.get('id') or ''
        desc = data.get('description') or ''
        return ensembl_id, desc
    return '', ''


def _try_mygene_summary(session: requests.Session, gene: str) -> str:
    try:
        r = session.get(
            f"https://mygene.info/v3/query?q={gene}&species=human&size=1",
            headers={"Accept": "application/json", "User-Agent": HEADERS["User-Agent"]},
            timeout=10.0,
        )
        if r.ok:
            hits = (r.json() or {}).get('hits', [])
            if hits:
                return hits[0].get('summary') or ''
    except Exception:
        pass
    return ''


@lru_cache(maxsize=2048)
def _fetch_ensembl_id(gene: str) -> Tuple[str, str]:
    """
    Busca (ensembl_id, description) robustamente:
    1) xrefs/symbol (prefer type=gene)
    2) si falta descripción, lookup/id/{id}
    3) si falta id/desc, lookup/symbol
    4) si aún falta desc, resumen desde MyGene.info
    """
    s = requests.Session()
    gene_str = str(gene).strip()

    ensembl_id, desc = _try_xrefs_symbol(s, gene_str)
    # Siempre intenta completar/refrescar descripción con lookup/id si ya hay ID
    if ensembl_id:
        d_lookup = _try_lookup_id_desc(s, ensembl_id)
        desc = d_lookup or desc
    if (not ensembl_id) or (not desc):
        e2, d2 = _try_lookup_symbol(s, gene_str)
        ensembl_id = ensembl_id or e2
        desc = desc or d2
    if ensembl_id and not desc:
        desc = _try_mygene_summary(s, gene_str)

    return (ensembl_id if ensembl_id else 'Not found', desc if desc else 'No description')


def add_ensembl_info_batch(df: pd.DataFrame, symbol_col: str = 'target', max_workers: int = 6) -> pd.DataFrame:
    symbols: List[str] = df[symbol_col].astype(str).tolist()
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        results = list(ex.map(_fetch_ensembl_id, symbols))
    ids, descriptions = zip(*results) if results else ([], [])
    out = df.copy()
    out['ensembl_id'] = ids
    out['description'] = descriptions
    return out
