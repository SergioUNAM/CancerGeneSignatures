from __future__ import annotations

from dataclasses import dataclass
from typing import Dict

import numpy as np
import pandas as pd


@dataclass
class AnalysisResult:
    summary: pd.DataFrame
    correlation: pd.DataFrame


def summarize(df: pd.DataFrame) -> pd.DataFrame:
    desc = df.describe(include=[np.number]).T
    # Add NA counts
    desc["na_count"] = df.isna().sum()
    return desc


def correlation(df: pd.DataFrame, method: str = "pearson") -> pd.DataFrame:
    return df.select_dtypes(include=[np.number]).corr(method=method)


def run_basic_analysis(df: pd.DataFrame) -> AnalysisResult:
    return AnalysisResult(summary=summarize(df), correlation=correlation(df))

