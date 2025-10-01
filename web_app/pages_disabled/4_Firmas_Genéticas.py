# -*- coding: utf-8 -*-
from __future__ import annotations

import streamlit as st
import pandas as pd
from pathlib import Path
import sys
from typing import Optional

# Ensure project root on path
_PROJ_ROOT = Path(__file__).resolve().parents[2]
if str(_PROJ_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJ_ROOT))

from app.core.signatures import create_signatures, HallmarkConfig

st.set_page_config(page_title="CGS — Firmas Genéticas", page_icon="🧬", layout="wide")
st.title("Firmas Genéticas y Hallmarks")
st.caption("Construye firmas ajustando parámetros; requiere bibliografía clasificada de la página PubMed.")

bib = st.session_state.get("bibliografia_clasificada")
if not isinstance(bib, pd.DataFrame) or bib.empty:
    st.info("Primero genera la bibliografía clasificada en la página 'Bibliografía (PubMed)'.")
    st.stop()

ctx = st.selectbox("Contexto biológico", options=sorted(bib['Context'].dropna().unique().tolist()) if 'Context' in bib.columns else ["General"])
hall_gmt = st.text_input("Ruta a Hallmark GMT (opcional)", value="")
back_gmt = st.text_input("Ruta a background GMT (opcional)", value="")

run = st.button("Construir firmas", type="primary")

if run:
    cfg = HallmarkConfig(hallmark_gmt=hall_gmt or None, background_gmt=back_gmt or None)
    with st.spinner("Calculando firmas…"):
        df_sigs = create_signatures(bib, contexto_biologico=ctx, hallmark_cfg=cfg)
    if df_sigs is None or df_sigs.empty:
        st.info("No se pudieron construir firmas con los parámetros actuales.")
    else:
        st.session_state['df_signatures'] = df_sigs.copy()
        st.success(f"Firmas generadas: {len(df_sigs)} registros")
        st.dataframe(df_sigs.head(50))
        st.download_button(
            label="Descargar firmas (CSV)",
            data=df_sigs.to_csv(index=False),
            file_name="firmas_geneticas.csv",
            mime="text/csv",
            use_container_width=True,
        )
