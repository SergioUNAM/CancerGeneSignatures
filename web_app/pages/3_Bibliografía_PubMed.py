# -*- coding: utf-8 -*-
from __future__ import annotations

import os
import streamlit as st
import pandas as pd
from pathlib import Path
import sys

# Ensure project root on path
_PROJ_ROOT = Path(__file__).resolve().parents[2]
if str(_PROJ_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJ_ROOT))

from src.core.bibliography import search_pubmed_by_genes, classify_bibliography, aggregate_counts_by_level_and_cancer

st.set_page_config(page_title="CGS — Bibliografía PubMed", page_icon="📚", layout="wide")
st.title("Bibliografía (PubMed)")
st.caption("Busca artículos por gen y contexto. Usa la tabla de expresión de la página principal.")

df_expr: pd.DataFrame | None = st.session_state.get('df_expr')
context_sel_label = st.session_state.get('context_sel_label', 'TEM')
if df_expr is None or df_expr.empty:
    st.info("Primero ejecuta el análisis qPCR en la página principal para generar la tabla de expresión.")
    st.stop()

env_email = os.getenv("NCBI_EMAIL", "")
env_key = os.getenv("NCBI_API_KEY", "")
c1, c2 = st.columns(2)
with c1:
    ncbi_email_input = st.text_input("NCBI Email (obligatorio)", value=env_email, placeholder="tu_email@dominio.com")
with c2:
    ncbi_api_key_input = st.text_input("NCBI API Key (opcional)", value=env_key, placeholder="api_key")
max_per_gene = st.number_input("Máximo de artículos por gen", value=100, min_value=10, max_value=300, step=10)
run = st.button("Buscar en PubMed", disabled=not bool(ncbi_email_input.strip()))

if run:
    # Preparar genes (sin necesitar Ensembl ID aquí)
    tmp = df_expr[['target', 'nivel_expresion']].drop_duplicates('target').copy()
    tmp['ensembl_id'] = ''
    genes_df = tmp[['target', 'ensembl_id', 'nivel_expresion']]

    os.environ["NCBI_EMAIL"] = ncbi_email_input.strip()
    if ncbi_api_key_input.strip():
        os.environ["NCBI_API_KEY"] = ncbi_api_key_input.strip()

    prog = st.progress(0)
    status = st.empty()

    def _on_progress(i: int, total: int, gene: str) -> None:
        pct = int((i / max(1, total)) * 100)
        prog.progress(pct)
        status.info(f"Procesando {i}/{total}: {gene}")

    with st.spinner("Consultando PubMed por gen…"):
        bib = search_pubmed_by_genes(
            genes_df.rename(columns={'target': 'target'}),
            symbol_col='target',
            ensembl_col='ensembl_id',
            selected_context=context_sel_label,
            max_per_gene=int(max_per_gene),
            progress=_on_progress,
            logger=None,
        )
    prog.progress(100)
    status.success("Consulta PubMed finalizada")
    if bib is None or bib.empty:
        st.info("No se encontraron artículos para los parámetros actuales.")
    else:
        bib2 = pd.merge(
            bib,
            genes_df.rename(columns={'target': 'Gene'})[['Gene', 'nivel_expresion']],
            on='Gene', how='left'
        )
        st.dataframe(bib2.head(50))
        st.download_button(
            label="Descargar bibliografía (CSV)",
            data=bib2.to_csv(index=False),
            file_name="bibliografia_pubmed.csv",
            mime="text/csv",
            use_container_width=True,
        )

        st.markdown("### Clasificación por tipo de cáncer y contexto")
        classified = classify_bibliography(bib2)
        st.session_state["bibliografia_clasificada"] = classified.copy()
        st.dataframe(classified.head(50))
        st.download_button(
            label="Descargar bibliografía clasificada (CSV)",
            data=classified.to_csv(index=False),
            file_name="bibliografia_clasificada.csv",
            mime="text/csv",
            use_container_width=True,
        )

        agg = aggregate_counts_by_level_and_cancer(classified)
        if not agg.empty:
            import plotly.express as px
            order_levels = ['sobreexpresado', 'estable', 'subexpresado']
            agg['nivel_expresion'] = pd.Categorical(agg['nivel_expresion'], categories=order_levels, ordered=True)
            fig = px.bar(
                agg.sort_values(['nivel_expresion','n'], ascending=[True, False]),
                x='cancer_type', y='n', color='nivel_expresion', barmode='group',
                labels={'cancer_type': 'Tipo de cáncer', 'n': 'Número de artículos', 'nivel_expresion': 'Nivel'},
                title='Artículos por tipo de cáncer y nivel de expresión'
            )
            fig.update_layout(height=500, margin=dict(t=60, b=60))
            st.plotly_chart(fig, use_container_width=True)

