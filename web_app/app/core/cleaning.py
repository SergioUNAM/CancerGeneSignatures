from __future__ import annotations

from typing import Iterable, List

import pandas as pd


def drop_machine_controls(
    df_long: pd.DataFrame,
    column: str = "target",
    controls: Iterable[str] | None = None,
) -> pd.DataFrame:
    """
    Elimina filas cuyo `column` inicia con alguno de los prefijos de controles de máquina.

    - df_long: DataFrame en formato largo con columnas al menos: `target` y opcionalmente `test`, `ct`.
    - column: columna sobre la que se aplica el filtro (por defecto `target`).
    - controls: prefijos de controles de máquina. Si no se proveen, usa una lista por defecto.

    Devuelve un nuevo DataFrame filtrado.
    """
    if column not in df_long.columns:
        raise ValueError(f"La columna '{column}' no existe en el DataFrame.")

    default_controls: List[str] = [
        "BUBL",         # Burbujas o similares
        "NTC",          # No Template Control
        "NTC\t",
        "NTC\n",
        "CC_NEGATIVO",  # Control negativo
        "NEGATIVO",
        "CONTROL_NEGATIVO",
        "CONTROL_POSITIVO",
        "POS",          # Positivo
        "NEG",          # Negativo
        "CONTROL",
    ]

    # Unificar: usar lista por defecto + extras provistos (p. ej., PPC/RTC)
    prefixes = (default_controls + list(controls)) if controls is not None else default_controls
    series = df_long[column].astype(str)
    mask_keep = ~series.str.startswith(tuple(prefixes))
    return df_long.loc[mask_keep].copy()
