# -*- coding: utf-8 -*-
from __future__ import annotations

import os
import streamlit as st
import pandas as pd
from typing import Optional

# Ensure project root is importable para acceder al paquete `app`
import sys
from pathlib import Path
_PROJ_ROOT = Path(__file__).resolve().parents[2]
if str(_PROJ_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJ_ROOT))

from app.core.string_enrichment import enrich_by_levels, filter_enrichment, dfs_to_excel_bytes

st.set_page_config(page_title="CGS — Enriquecimiento STRING", page_icon="🧬", layout="wide")
st.title("Enriquecimiento funcional (STRING)")
st.caption("Analiza términos enriquecidos por nivel de expresión. Usa los resultados de qPCR de la página principal.")

df_expr: Optional[pd.DataFrame] = st.session_state.get('df_expr')
if df_expr is None or df_expr.empty:
    st.info("Primero ejecuta el análisis qPCR en la página principal para generar la tabla de expresión categorizada.")
    st.stop()

# Configuración
order_levels = ['subexpresado', 'estable', 'sobreexpresado']
levels = st.multiselect("Niveles a considerar", options=order_levels, default=['subexpresado', 'sobreexpresado'])
cat_options = ["GO", "GO:BP", "GO:MF", "GO:CC", "KEGG", "Reactome"]
sel_cats = st.multiselect("Categorías (fuentes)", options=cat_options, default=["GO", "KEGG"])
c1, c2, c3, c4 = st.columns(4)
with c1:
    max_fdr = st.number_input("FDR máx.", value=0.05, min_value=0.0, max_value=1.0, step=0.01, format="%.2f")
with c2:
    min_size = st.number_input("Mín. genes por término", value=3, min_value=1, max_value=100, step=1)
with c3:
    top_n = st.number_input("Top N", value=25, min_value=1, max_value=200, step=1)
with c4:
    species = st.selectbox("Especie (NCBI Taxon)", options=[9606, 10090, 10116], index=0)

run = st.button("Ejecutar enriquecimiento (STRING)", type="primary")

if run:
    with st.spinner("Consultando STRING por nivel…"):
        res = enrich_by_levels(
            df_expr,
            symbol_col='target',
            level_col='nivel_expresion',
            levels=levels or order_levels,
            species=int(species),
            sources=sel_cats,
        )
    combined = res.get('combined')
    if combined is None or combined.empty:
        st.info("Sin resultados de enriquecimiento con los parámetros actuales.")
    else:
        st.success(f"Resultados obtenidos: {len(combined)} filas totales")
        # Filtrar por criterios
        filtered = filter_enrichment(combined, include_categories=sel_cats, max_fdr=float(max_fdr), min_term_genes=int(min_size), top_n=int(top_n))
        st.dataframe(filtered)
        # Exportar a Excel por nivel
        by_level = res.get('by_level') or {}
        sheets = []
        names = []
        for lvl, df in by_level.items():
            f = filter_enrichment(df, include_categories=sel_cats, max_fdr=float(max_fdr), min_term_genes=int(min_size), top_n=int(top_n))
            sheets.append(f)
            names.append(str(lvl))
        xls = dfs_to_excel_bytes(sheets, names)
        st.download_button(
            label="Descargar resultados por nivel (XLSX)",
            data=xls,
            file_name="string_enrichment.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            use_container_width=True,
        )
