from __future__ import annotations

from typing import Optional, Sequence, Tuple

import pandas as pd

_SENTINEL_MISSING = "__CGS_MISSING__"
_GLOBAL_KEY = "__GLOBAL__"


def _ensure_column(df: pd.DataFrame, column: str) -> None:
    if column not in df.columns:
        raise KeyError(f"La columna '{column}' no se encuentra en el DataFrame.")


def _clean_ct_column(
    df: pd.DataFrame,
    column: str,
    *,
    inplace: bool,
    undetermined_tokens: Sequence[str],
) -> pd.DataFrame:
    data = df if inplace else df.copy()
    _ensure_column(data, column)
    token_map = {str(t).strip().lower(): pd.NA for t in undetermined_tokens}
    series = (
        data[column]
        .astype("string")
        .str.strip()
        .str.lower()
        .replace(token_map)
    )
    data[column] = pd.to_numeric(series, errors="coerce")
    return data


def _build_group_key(df: pd.DataFrame, keys: Sequence[str]) -> pd.Series:
    if not keys:
        return pd.Series([_GLOBAL_KEY] * len(df), index=df.index)
    joined = None
    for key in keys:
        if key not in df.columns:
            raise KeyError(f"La columna de agrupación '{key}' no se encuentra en el DataFrame.")
        part = df[key].astype("string").fillna(_SENTINEL_MISSING)
        joined = part if joined is None else joined.str.cat(part, sep="||")
    assert joined is not None  # for mypy; guarded by not keys check
    return joined


def procesar_ct_column(
    df1: pd.DataFrame,
    df2: pd.DataFrame,
    columna_ct: str = "ct",
    *,
    max_ct: Optional[float] = None,
    by: Optional[Sequence[str]] = None,
    undetermined_tokens: Sequence[str] = (
        "undetermined",
        "undet",
        "nd",
        "na",
    ),
    inplace: bool = False,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Imputa valores Ct ausentes o no determinados.

    Por defecto, limpia valores string como "Undetermined" y los convierte en NaN. Luego
    imputa usando:
    - ``max_ct`` si se provee (límite técnico).
    - Máximo observado por grupos definidos en ``by`` si no se entrega ``max_ct``.
      Si aún quedan NaN, se aplica un fallback al máximo global.
    - Máximo global si no se proporciona ``by``.

    Se añade ``<columna_ct>_imputed`` para marcar filas imputadas y se registra información
    en ``DataFrame.attrs['imputation_info']``.
    """

    group_cols: Sequence[str] = tuple(by) if by else ()

    a = _clean_ct_column(df1, columna_ct, inplace=inplace, undetermined_tokens=undetermined_tokens)
    b = _clean_ct_column(df2, columna_ct, inplace=inplace, undetermined_tokens=undetermined_tokens)

    a_flag = a[columna_ct].isna()
    b_flag = b[columna_ct].isna()

    combined = pd.concat([a[columna_ct], b[columna_ct]], ignore_index=True)
    if max_ct is not None:
        fill_value = float(max_ct)
        a.loc[a_flag, columna_ct] = fill_value
        b.loc[b_flag, columna_ct] = fill_value
        used_info = f"Imputación por límite técnico max_ct={fill_value}"
    else:
        if group_cols:
            a_keys = _build_group_key(a, group_cols)
            b_keys = _build_group_key(b, group_cols)
            comb_keys = pd.concat([a_keys, b_keys], axis=0)
            grp_max = combined.groupby(comb_keys).max()
            grp_dict = grp_max.to_dict()
            a_group_vals = a_keys.map(grp_dict)
            b_group_vals = b_keys.map(grp_dict)
            a.loc[a_flag, columna_ct] = a_group_vals[a_flag]
            b.loc[b_flag, columna_ct] = b_group_vals[b_flag]
            used_info = f"Imputación por máximo observado estratificado por {list(group_cols)}"
        else:
            used_info = "Imputación por máximo observado global"

        combined_after = pd.concat([a[columna_ct], b[columna_ct]], ignore_index=True)
        remaining_nan = combined_after.isna()
        if remaining_nan.any():
            fallback_max = combined_after.max()
            if pd.isna(fallback_max):
                raise ValueError(
                    "No hay valores numéricos válidos de Ct tras la imputación; "
                    "proporciona 'max_ct' o revisa los datos."
                )
            a_missing = a[columna_ct].isna()
            b_missing = b[columna_ct].isna()
            a.loc[a_missing, columna_ct] = fallback_max
            b.loc[b_missing, columna_ct] = fallback_max
            used_info += f" + fallback global={fallback_max:.3f}"

    a[f"{columna_ct}_imputed"] = a_flag
    b[f"{columna_ct}_imputed"] = b_flag

    for frame, name in ((a, "df1"), (b, "df2")):
        if not pd.api.types.is_float_dtype(frame[columna_ct]):
            frame[columna_ct] = frame[columna_ct].astype(float)
        out_of_range = frame.loc[
            (frame[columna_ct] < 0) | (frame[columna_ct] > 50),
            columna_ct,
        ]
        if not out_of_range.empty:
            raise ValueError(f"{name}: Ct fuera de rango plausible: {out_of_range.unique()}")

        frame.attrs["imputation_info"] = used_info

    return a, b


def ctmax(
    df1: pd.DataFrame,
    df2: pd.DataFrame,
    columna_ct: str = "ct",
    *,
    max_ct: Optional[float] = None,
    by: Optional[Sequence[str]] = None,
    undetermined_tokens: Sequence[str] = (
        "undetermined",
        "undet",
        "nd",
        "na",
    ),
    inplace: bool = False,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Alias for procesar_ct_column, maintained for compatibility.
    Preserves the exact behavior of procesar_ct_column, including group-wise maximum and global fallback.
    """
    return procesar_ct_column(
        df1,
        df2,
        columna_ct=columna_ct,
        max_ct=max_ct,
        by=by,
        undetermined_tokens=undetermined_tokens,
        inplace=inplace,
    )
