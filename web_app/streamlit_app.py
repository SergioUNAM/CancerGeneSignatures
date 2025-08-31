from __future__ import annotations

import io
import zipfile
from datetime import datetime
from pathlib import Path
from typing import Optional

import pandas as pd
import streamlit as st

from src.core.io import LoadResult, list_excel_sheets, load_table, parse_qpcr_wide
from src.core.qpcr import try_extract_template_config, melt_wide_to_long, classify_tests
from src.core.preprocessing import PreprocessConfig, preprocess
from src.core.analysis import run_basic_analysis
from src.core.plots import corr_heatmap, histogram, scatter


st.set_page_config(page_title="CancerGeneSignatures - Web App", layout="wide")

st.title("CancerGeneSignatures • Análisis interactivo")
st.write(
    "Sube un archivo Excel/CSV, elige parámetros y genera resultados con tablas y gráficas."
)

with st.sidebar:
    st.header("1) Datos de entrada")
    uploaded = st.file_uploader("Archivo (.xlsx, .csv, .tsv)", type=["xlsx", "xls", "csv", "tsv"])
    sheet: Optional[str] = None
    if uploaded is not None and uploaded.name.lower().endswith((".xlsx", "xls")):
        try:
            # Leer nombres de hojas (reabrir el buffer cada vez por seguridad)
            sheets = list_excel_sheets(uploaded)
            if sheets:
                sheet = st.selectbox("Hoja de Excel", options=sheets, index=0)
        except Exception as e:
            st.warning(f"No se pudieron listar hojas: {e}")

    st.header("2) Parámetros")
    normalization = st.selectbox("Normalización", options=["none", "zscore", "minmax"], index=0)
    data_format = st.radio(
        "Formato de datos",
        options=["Tabla simple", "qPCR (Well/Target Name + CT por muestra)"],
        index=1,
    )
    if data_format.startswith("qPCR"):
        st.caption("Se intentará detectar encabezados con Well/Target Name y nombres de muestras.")
        undet_policy = st.selectbox(
            "Tratamiento de 'Undetermined'",
            options=["nan", "ctmax", "value"],
            help="nan: deja como faltante; ctmax: reemplaza por el CT máximo observado de la muestra; value: usar un valor fijo.",
        )
        undet_value = st.number_input("Valor fijo para 'value'", value=40.0, step=0.5)
    else:
        undet_policy = "nan"
        undet_value = 40.0
    # Configuración de contexto biológico / cáncer
    st.header("3) Configuración (proyecto)")
    auto_config = False
    contexto_biologico = st.text_input("Contexto biológico", value="")
    tipo_cancer = st.text_input("Tipo de cáncer", value="")
    metodo = st.text_input("Método (descriptivo)", value="")

    st.caption("Prefijos para clasificar pruebas en controles/muestras (p. ej., 4GB, 3CG)")
    pref_ctrl = st.text_input("Prefijo controles", value="")
    pref_samp = st.text_input("Prefijo muestras", value="")

    id_cols_input = st.text_input(
        "Columnas ID (opcional, separadas por comas)",
        value="",
        help="Ej.: muestra, grupo"
    )
    id_columns = [c.strip() for c in id_cols_input.split(",") if c.strip()]
    run_btn = st.button("Ejecutar análisis", type="primary")


def build_results_zip(df_raw: pd.DataFrame, df_proc: pd.DataFrame, analysis,
                      cfg: Optional[dict] = None,
                      qpcr_long: Optional[pd.DataFrame] = None,
                      controles: Optional[pd.DataFrame] = None,
                      muestras: Optional[pd.DataFrame] = None):
    mem = io.BytesIO()
    time_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    with zipfile.ZipFile(mem, mode="w", compression=zipfile.ZIP_DEFLATED) as z:
        # Datos
        z.writestr(f"results_{time_id}/datos_raw.csv", df_raw.to_csv(index=False))
        z.writestr(f"results_{time_id}/datos_procesados.csv", df_proc.to_csv(index=False))
        z.writestr(f"results_{time_id}/resumen.csv", analysis.summary.to_csv())
        z.writestr(f"results_{time_id}/correlacion.csv", analysis.correlation.to_csv())
        # Config y qPCR
        if cfg is not None:
            import json
            z.writestr(f"results_{time_id}/config.json", json.dumps(cfg, ensure_ascii=False, indent=2))
        if qpcr_long is not None:
            z.writestr(f"results_{time_id}/qpcr_long.csv", qpcr_long.to_csv(index=False))
        if controles is not None and not controles.empty:
            z.writestr(f"results_{time_id}/controles.csv", controles.to_csv(index=False))
        if muestras is not None and not muestras.empty:
            z.writestr(f"results_{time_id}/muestras.csv", muestras.to_csv(index=False))
    mem.seek(0)
    return mem, f"resultados_{time_id}.zip"


df_loaded: Optional[LoadResult] = None
if uploaded is not None and (run_btn or st.session_state.get("auto_run", True)):
    try:
        if data_format.startswith("qPCR") and uploaded.name.lower().endswith((".xlsx", "xls")):
            df_loaded = parse_qpcr_wide(
                uploaded,
                sheet_name=sheet,
                undetermined_policy=undet_policy,
                undetermined_value=float(undet_value),
            )
            # Autocompletar columnas ID típicas del formato qPCR
            if not id_columns:
                id_columns = ["Well", "Target Name"]
            # Intento de autoconfig desde plantilla
            try:
                cfg = try_extract_template_config(uploaded, sheet_name=sheet)
                if cfg.contexto_biologico and not contexto_biologico:
                    contexto_biologico = cfg.contexto_biologico
                if cfg.tipo_cancer and not tipo_cancer:
                    tipo_cancer = cfg.tipo_cancer
                if cfg.metodo and not metodo:
                    metodo = cfg.metodo
                if cfg.prefijo_controles and not pref_ctrl:
                    pref_ctrl = cfg.prefijo_controles
                if cfg.prefijo_muestras and not pref_samp:
                    pref_samp = cfg.prefijo_muestras
            except Exception:
                pass
        else:
            df_loaded = load_table(uploaded, file_name=uploaded.name, sheet_name=sheet)
    except Exception as e:
        st.error(f"Error al cargar el archivo: {e}")

if df_loaded is not None:
    st.subheader("Vista previa de datos")
    st.caption(f"Fuente: {df_loaded.source_name} | Hoja: {df_loaded.sheet_name or '-'} | Forma: {df_loaded.df.shape}")
    st.dataframe(df_loaded.df.head(20), use_container_width=True)

    cfg = PreprocessConfig(id_columns=id_columns or None, numeric_only=True, normalization=normalization)
    try:
        df_proc = preprocess(df_loaded.df, cfg)
    except Exception as e:
        st.error(f"Error en preprocesamiento: {e}")
        st.stop()

    analysis = run_basic_analysis(df_proc)

    st.subheader("Configuración del proyecto")
    st.json({
        "contexto_biologico": contexto_biologico or None,
        "tipo_cancer": tipo_cancer or None,
        "metodo": metodo or None,
        "prefijos": {"controles": pref_ctrl or None, "muestras": pref_samp or None},
    })

    st.subheader("Resultados")
    st.markdown("- Resumen estadístico de variables numéricas\n- Matriz de correlación\n- Clasificación de pruebas (controles vs. muestras)")
    with st.expander("Resumen (describe)", expanded=True):
        st.dataframe(analysis.summary, use_container_width=True)
    with st.expander("Correlación (tabla)", expanded=False):
        st.dataframe(analysis.correlation, use_container_width=True)

    st.subheader("Visualizaciones")
    # Correlación
    st.plotly_chart(corr_heatmap(analysis.correlation), use_container_width=True)

    # Histograma rápido si hay columnas numéricas
    num_cols = df_proc.select_dtypes(include=["number"]).columns.tolist()
    if num_cols:
        sel_col = st.selectbox("Histograma de columna", options=num_cols)
        st.plotly_chart(histogram(df_proc, sel_col), use_container_width=True)

    # Clasificación qPCR si aplica
    long_df = None
    controles_df = None
    muestras_df = None
    if data_format.startswith("qPCR"):
        long_df = melt_wide_to_long(df_loaded.df)
        controles_df, muestras_df = classify_tests(long_df, pref_ctrl, pref_samp)
        with st.expander("Controles detectados", expanded=False):
            if controles_df is not None and not controles_df.empty:
                st.write(
                    f"{controles_df['test'].nunique()} pruebas • {controles_df['target'].nunique()} genes"
                )
                st.dataframe(controles_df.head(20))
            else:
                st.info("No se detectaron controles (verifica el prefijo).")
        with st.expander("Muestras detectadas", expanded=False):
            if muestras_df is not None and not muestras_df.empty:
                st.write(
                    f"{muestras_df['test'].nunique()} pruebas • {muestras_df['target'].nunique()} genes"
                )
                st.dataframe(muestras_df.head(20))
            else:
                st.info("No se detectaron muestras (verifica el prefijo).")

    # Descarga
    cfg_json = {
        "contexto_biologico": contexto_biologico or None,
        "tipo_cancer": tipo_cancer or None,
        "metodo": metodo or None,
        "prefijos": {"controles": pref_ctrl or None, "muestras": pref_samp or None},
        "archivo": df_loaded.source_name,
        "hoja": df_loaded.sheet_name,
        "normalizacion": normalization,
    }
    mem, name = build_results_zip(
        df_loaded.df,
        df_proc,
        analysis,
        cfg=cfg_json,
        qpcr_long=long_df,
        controles=controles_df,
        muestras=muestras_df,
    )
    st.download_button("Descargar resultados (ZIP)", data=mem, file_name=name, mime="application/zip")


st.markdown("---")
st.markdown(
    """
    Sugerencias de uso:
    1. Sube un archivo Excel (elige la hoja) o CSV/TSV.
    2. Opcionalmente indica columnas ID (no numéricas) separadas por coma.
    3. Elige la normalización y ejecuta el análisis.
    4. Revisa tablas y gráficas; descarga resultados en ZIP.
    """
)
