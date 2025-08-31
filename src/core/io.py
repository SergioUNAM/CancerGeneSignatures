from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple, Union

import pandas as pd


@dataclass
class LoadResult:
    df: pd.DataFrame
    source_name: str
    sheet_name: Optional[str] = None


def is_excel(filename: str) -> bool:
    return filename.lower().endswith((".xlsx", ".xlsm", ".xls"))


def is_csv(filename: str) -> bool:
    return filename.lower().endswith((".csv", ".tsv"))


def list_excel_sheets(file: Union[str, Path, "IO"]):
    """Return sheet names if file is an Excel workbook; else empty list."""
    try:
        xls = pd.ExcelFile(file)
        return list(xls.sheet_names)
    except Exception:
        return []


def load_table(
    file: Union[str, Path, "IO"],
    file_name: Optional[str] = None,
    sheet_name: Optional[str] = None,
    csv_sep: str = ",",
    csv_decimal: str = ".",
) -> LoadResult:
    """
    Load CSV/TSV or Excel into a DataFrame, with optional sheet selection.
    - file: path or file-like object
    - file_name: used to infer type when file-like
    """
    source_name = file_name or (str(file) if not hasattr(file, "name") else getattr(file, "name", "uploaded"))

    if is_excel(source_name):
        df = pd.read_excel(file, sheet_name=sheet_name)  # requires openpyxl in most cases
        return LoadResult(df=df, source_name=source_name, sheet_name=sheet_name)

    if is_csv(source_name):
        sep = "\t" if source_name.lower().endswith(".tsv") else csv_sep
        df = pd.read_csv(file, sep=sep, decimal=csv_decimal)
        return LoadResult(df=df, source_name=source_name, sheet_name=None)

    # Try flexible fallback: attempt CSV then Excel
    try:
        df = pd.read_csv(file, sep=csv_sep, decimal=csv_decimal)
        return LoadResult(df=df, source_name=source_name, sheet_name=None)
    except Exception:
        df = pd.read_excel(file, sheet_name=sheet_name)
        return LoadResult(df=df, source_name=source_name, sheet_name=sheet_name)

