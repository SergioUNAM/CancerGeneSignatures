from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Tuple, List
from collections import Counter
import re

import pandas as pd


@dataclass
class QpcrTemplateConfig:
    contexto_biologico: Optional[str] = None
    metodo: Optional[str] = None
    tipo_cancer: Optional[str] = None
    prefijo_controles: Optional[str] = None
    prefijo_muestras: Optional[str] = None


def try_extract_template_config(file, sheet_name: Optional[str] = None) -> QpcrTemplateConfig:
    """
    Intenta extraer configuración de la parte superior de la plantilla Excel (estilo UIMEO).
    Supone que los valores están en celdas B1..B5 (índices 0..4, col 1). Es tolerante a errores.
    """
    try:
        head = pd.read_excel(file, sheet_name=sheet_name, header=None, nrows=8)
    except Exception:
        return QpcrTemplateConfig()

    def at(r, c):
        try:
            v = head.iat[r, c]
            return None if pd.isna(v) else str(v)
        except Exception:
            return None

    return QpcrTemplateConfig(
        tipo_cancer=at(0, 1),
        metodo=at(1, 1),
        contexto_biologico=at(2, 1),
        prefijo_controles=(at(3, 1) or ""),
        prefijo_muestras=(at(4, 1) or ""),
    )


def melt_wide_to_long(df_wide: pd.DataFrame) -> pd.DataFrame:
    """
    Convierte tabla ancha qPCR (Well, Target Name, <muestras>...) a formato largo:
    columnas: test (muestra), target, ct
    """
    id_vars = [c for c in ["Well", "Target Name"] if c in df_wide.columns]
    value_vars = [c for c in df_wide.columns if c not in id_vars]
    long = df_wide.melt(id_vars=id_vars, value_vars=value_vars, var_name="test", value_name="ct")
    # target
    if "Target Name" in long.columns:
        long = long.rename(columns={"Target Name": "target"})
    # Tipar ct
    long["ct"] = pd.to_numeric(long["ct"], errors="coerce")
    return long


def classify_tests(df_long: pd.DataFrame, ctrl_prefix: str, sample_prefix: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Separa df_long en controles y muestras por prefijo en columna 'test' (case-insensitive)."""
    t = df_long.copy()
    t["test_up"] = t["test"].astype(str).str.upper()
    ctrl = ctrl_prefix.strip().upper() if ctrl_prefix else ""
    samp = sample_prefix.strip().upper() if sample_prefix else ""

    controles_df = t[t["test_up"].str.startswith(ctrl)] if ctrl else t.iloc[0:0]
    muestras_df = t[t["test_up"].str.startswith(samp)] if samp else t.iloc[0:0]
    # Limpiar columnas auxiliares
    controles_df = controles_df.drop(columns=["test_up"]) if not controles_df.empty else controles_df
    muestras_df = muestras_df.drop(columns=["test_up"]) if not muestras_df.empty else muestras_df
    return controles_df, muestras_df


def suggest_name_affixes(sample_names: List[str], top_n: int = 10) -> Dict[str, List[Tuple[str, int]]]:
    """
    Sugiere prefijos y sufijos comunes a partir de los nombres de prueba.
    - Prefijo: parte antes del primer '-' (o hasta primer '_'); como alternativa, letras iniciales + dígitos iniciales.
    - Sufijo: dígitos finales o parte después del último '-'.
    Devuelve listas ordenadas por frecuencia descendente.
    """
    prefixes: Counter = Counter()
    suffixes: Counter = Counter()
    for name in sample_names:
        n = str(name).strip()
        if not n:
            continue
        # prefijos candidatos
        token_dash = n.split('-')[0]
        token_us = n.split('_')[0]
        m_lead = re.match(r'^([A-Za-z0-9]+)', n)
        for cand in {token_dash, token_us, m_lead.group(1) if m_lead else ''}:
            if cand:
                prefixes[cand] += 1
        # sufijos candidatos
        m_tailnum = re.search(r'(\d+)$', n)
        token_last_dash = n.split('-')[-1] if '-' in n else ''
        if m_tailnum:
            suffixes[m_tailnum.group(1)] += 1
        if token_last_dash:
            suffixes[token_last_dash] += 1
    pref_list = prefixes.most_common(top_n)
    suff_list = suffixes.most_common(top_n)
    return {"prefixes": pref_list, "suffixes": suff_list}


def classify_by_prefixes(df_long: pd.DataFrame, ctrl_prefixes: List[str], sample_prefixes: List[str]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    t = df_long.copy()
    t["test_str"] = t["test"].astype(str)
    ctrl = t[t["test_str"].str.startswith(tuple(ctrl_prefixes))] if ctrl_prefixes else t.iloc[0:0]
    samp = t[t["test_str"].str.startswith(tuple(sample_prefixes))] if sample_prefixes else t.iloc[0:0]
    return ctrl.drop(columns=["test_str"]) if not ctrl.empty else ctrl, samp.drop(columns=["test_str"]) if not samp.empty else samp


def classify_by_suffixes(df_long: pd.DataFrame, ctrl_suffixes: List[str], sample_suffixes: List[str]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    t = df_long.copy()
    t["test_str"] = t["test"].astype(str)
    ctrl = t[t["test_str"].str.endswith(tuple(ctrl_suffixes))] if ctrl_suffixes else t.iloc[0:0]
    samp = t[t["test_str"].str.endswith(tuple(sample_suffixes))] if sample_suffixes else t.iloc[0:0]
    return ctrl.drop(columns=["test_str"]) if not ctrl.empty else ctrl, samp.drop(columns=["test_str"]) if not samp.empty else samp


def classify_by_regex(df_long: pd.DataFrame, ctrl_pattern: str, sample_pattern: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Clasifica usando expresiones regulares sobre la columna 'test'.
    Si un patrón está vacío, devuelve DataFrame vacío para ese grupo.
    """
    t = df_long.copy()
    t["test_str"] = t["test"].astype(str)
    ctrl = (
        t[t["test_str"].str.contains(ctrl_pattern, regex=True, na=False)]
        if ctrl_pattern else t.iloc[0:0]
    )
    samp = (
        t[t["test_str"].str.contains(sample_pattern, regex=True, na=False)]
        if sample_pattern else t.iloc[0:0]
    )
    return ctrl.drop(columns=["test_str"]) if not ctrl.empty else ctrl, samp.drop(columns=["test_str"]) if not samp.empty else samp
