from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional

import numpy as np
import pandas as pd


@dataclass
class PreprocessConfig:
    id_columns: Optional[List[str]] = None  # columns to keep untouched
    numeric_only: bool = True
    normalization: str = "none"  # one of: none|zscore|minmax


def select_numeric(df: pd.DataFrame, id_columns: Optional[Iterable[str]] = None) -> pd.DataFrame:
    ids = list(id_columns or [])
    num = df.select_dtypes(include=[np.number]).copy()
    if ids:
        for c in ids:
            if c in df.columns and c not in num.columns:
                num[c] = df[c]
        # reorder: ids first
        cols = [c for c in ids if c in num.columns] + [c for c in num.columns if c not in ids]
        num = num[cols]
    return num


def normalize(df: pd.DataFrame, method: str) -> pd.DataFrame:
    if method == "none":
        return df
    out = df.copy()
    numeric_cols = out.select_dtypes(include=[np.number]).columns
    if method == "zscore":
        out[numeric_cols] = (out[numeric_cols] - out[numeric_cols].mean()) / (out[numeric_cols].std(ddof=0) + 1e-12)
    elif method == "minmax":
        mins = out[numeric_cols].min()
        maxs = out[numeric_cols].max()
        out[numeric_cols] = (out[numeric_cols] - mins) / (maxs - mins + 1e-12)
    else:
        raise ValueError(f"Unknown normalization: {method}")
    return out


def preprocess(df: pd.DataFrame, cfg: PreprocessConfig) -> pd.DataFrame:
    data = select_numeric(df, cfg.id_columns) if cfg.numeric_only else df.copy()
    data = normalize(data, cfg.normalization)
    return data

