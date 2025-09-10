from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple, Union

import numpy as np
import pandas as pd


@dataclass
class LoadResult:
    df: pd.DataFrame
    source_name: str
    sheet_name: Optional[str] = None
    meta: Optional[dict] = None


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


# --- qPCR helpers ---

def _cell_str(x) -> str:
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return ""
    return str(x).strip()


def parse_qpcr_wide(
    file: Union[str, Path, "IO"],
    sheet_name: Optional[str] = None,
    undetermined_policy: str = "nan",  # one of: nan|ctmax|value
    undetermined_value: float = 40.0,
    header_mode: str = "auto",  # one of: auto|coords
    header_row_idx: Optional[int] = None,  # 0-based; default 3 (A4/B4) when coords
    well_col_idx: Optional[int] = None,    # 0-based; default 0 (A)
    target_col_idx: Optional[int] = None,  # 0-based; default 1 (B)
) -> LoadResult:
    """
    Parse qPCR-like Excel where headers include a row with sample codes above a row of
    column labels like: [Well, Target Name, CT, CT, CT, ...]. Returns a wide DataFrame
    with columns: Well, Target Name, <sample1>, <sample2>, ...
    """
    raw = pd.read_excel(file, sheet_name=sheet_name, header=None, dtype=object)

    # Determine header by mode
    header_row = None
    well_col = None
    target_col = None

    def _detect_header_auto() -> tuple[int, int, int]:
        for ridx in range(len(raw)):
            row_vals = [(_cell_str(v).lower()) for v in raw.iloc[ridx].tolist()]
            if "well" in row_vals and any("target" in v for v in row_vals):
                h_row = ridx
                w_col = row_vals.index("well")
                t_col = next(i for i, v in enumerate(row_vals) if "target" in v)
                return h_row, w_col, t_col
        raise ValueError("No se encontró una fila de encabezado con 'Well' y 'Target'.")

    def _detect_header_coords() -> tuple[int, int, int]:
        h_row = 3 if header_row_idx is None else int(header_row_idx)
        w_col = 0 if well_col_idx is None else int(well_col_idx)
        t_col = 1 if target_col_idx is None else int(target_col_idx)
        return h_row, w_col, t_col

    if header_mode == "coords":
        header_row, well_col, target_col = _detect_header_coords()
    else:
        header_row, well_col, target_col = _detect_header_auto()

    # Determine sample columns: those after target_col that have any non-empty value below
    data_start = header_row + 1
    candidate_cols = list(range(target_col + 1, raw.shape[1]))
    sample_cols = [c for c in candidate_cols if raw.iloc[data_start:, c].notna().any()]

    # Heurística para encontrar la fila con nombres de pruebas/muestras:
    # En muchas plantillas, la fila justo arriba de la cabecera contiene 'CT'.
    # Los nombres de prueba suelen estar 2 filas arriba de la cabecera.
    def is_ct_row(ridx: int) -> bool:
        if ridx is None or ridx < 0:
            return False
        vals = [(_cell_str(v).lower()) for v in raw.iloc[ridx, sample_cols].tolist()]
        nonempty = [v for v in vals if v]
        if not nonempty:
            return False
        ct_ratio = sum(v == "ct" for v in vals) / max(1, len(vals))
        return ct_ratio >= 0.4  # mayoría son 'CT'

    candidate_rows = [header_row - 1, header_row - 2, header_row - 3]
    sample_name_row = None

    # Elegir la fila con más patrones de "código de muestra" (con guion),
    # y en empate, con más dígitos y más no vacíos distintos de 'ct'.
    eligible: list[tuple[int, int, int, int]] = []  # (ridx, hyphen_count, digit_count, good_count)
    for ridx in candidate_rows:
        if ridx is None or ridx < 0:
            continue
        if is_ct_row(ridx):
            continue
        vals = [(_cell_str(v)) for v in raw.iloc[ridx, sample_cols].tolist()]
        good = [v for v in vals if v and v.lower() != "ct"]
        if len(good) >= max(1, int(0.3 * len(sample_cols))):
            hyphen_count = sum('-' in v for v in good)
            digit_count = sum(any(ch.isdigit() for ch in v) for v in good)
            eligible.append((ridx, hyphen_count, digit_count, len(good)))

    if eligible:
        # Seleccionar con mayor hyphen_count, luego digit_count, luego good_count
        eligible.sort(key=lambda t: (t[1], t[2], t[3]), reverse=True)
        sample_name_row = eligible[0][0]
    if sample_name_row is None and header_row > 1 and not is_ct_row(header_row - 2):
        sample_name_row = header_row - 2

    sample_names: List[str] = []
    for idx, c in enumerate(sample_cols, start=1):
        name = _cell_str(raw.iat[sample_name_row, c]) if sample_name_row is not None else ""
        if not name or name.lower() == "ct":
            name = f"Sample_{idx}"
        sample_names.append(name)

    # Build clean DataFrame
    cols = [well_col, target_col] + sample_cols
    df = raw.iloc[data_start:, cols].copy()
    df.columns = ["Well", "Target Name"] + sample_names
    # Drop fully empty rows (e.g., tail)
    df = df[~(df["Well"].isna() & df["Target Name"].isna())]
    # Normalize strings
    df["Well"] = df["Well"].apply(_cell_str)
    df["Target Name"] = df["Target Name"].apply(_cell_str)

    # Handle 'Undetermined' and cast to numeric for sample columns
    sample_cols_names = sample_names
    undet_tokens = {"undetermined", "undet", "nd", "neg"}
    for c in sample_cols_names:
        s = df[c].apply(_cell_str)
        mask_undet = s.str.lower().isin(undet_tokens)
        # Coerce numeric
        num = pd.to_numeric(df[c], errors="coerce")
        # Apply policy
        if undetermined_policy == "nan":
            num = num.mask(mask_undet, np.nan)
        elif undetermined_policy == "ctmax":
            # Use max observed per column (excluding undet/NaN) or fallback undetermined_value
            col_max = pd.to_numeric(s, errors="coerce").max()
            fill_val = float(col_max) if pd.notna(col_max) else float(undetermined_value)
            num = num.mask(mask_undet, fill_val)
        elif undetermined_policy == "value":
            num = num.mask(mask_undet, float(undetermined_value))
        else:
            raise ValueError("undetermined_policy debe ser uno de: nan|ctmax|value")
        df[c] = num

    meta = {
        "header_row": int(header_row),
        "sample_name_row": int(sample_name_row) if sample_name_row is not None else None,
        "sample_names": sample_names,
    }
    return LoadResult(df=df, source_name=getattr(file, "name", "uploaded"), sheet_name=sheet_name, meta=meta)
