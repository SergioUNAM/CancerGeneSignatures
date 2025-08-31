from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from typing import List, Tuple

import requests
import pandas as pd


def _fetch_ensembl_id(gene: str) -> Tuple[str, str]:
    """Busca un ID en Ensembl para un símbolo de gen humano. Devuelve (ensembl_id, description)."""
    base = "https://rest.ensembl.org"
    headers = {"Content-Type": "application/json"}
    try:
        r = requests.get(f"{base}/xrefs/symbol/homo_sapiens/{gene}?content-type=application/json", headers=headers, timeout=10)
        if r.ok:
            data = r.json()
            if data:
                ensembl_id = data[0].get('id', 'Not found')
                desc = data[0].get('description', '')
                return ensembl_id, desc or 'No description'
    except Exception:
        pass
    return 'Not found', 'No description'


def add_ensembl_info_batch(df: pd.DataFrame, symbol_col: str = 'target', max_workers: int = 6) -> pd.DataFrame:
    symbols: List[str] = df[symbol_col].astype(str).tolist()
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        results = list(ex.map(_fetch_ensembl_id, symbols))
    ids, descriptions = zip(*results)
    out = df.copy()
    out['ensembl_id'] = ids
    out['description'] = descriptions
    return out

