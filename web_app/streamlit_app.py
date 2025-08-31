from __future__ import annotations

import io
import zipfile
from datetime import datetime
from pathlib import Path
import json
from typing import Optional

import pandas as pd
import streamlit as st

from src.core.io import LoadResult, list_excel_sheets, load_table, parse_qpcr_wide
from src.core.qpcr import (
    try_extract_template_config,
    melt_wide_to_long,
    classify_tests,
    suggest_name_affixes,
    classify_by_prefixes,
    classify_by_suffixes,
)
from src.core.preprocessing import PreprocessConfig, preprocess
from src.core.analysis import run_basic_analysis
from src.core.plots import corr_heatmap, histogram, scatter


st.set_page_config(page_title="CancerGeneSignatures - Web App", layout="wide")

st.title("Análisis de datos de expresión diferencial por qPCR para identificar firmas génicas y redes de interacción en cánceres específicos")
st.write(
    "A partir de datos de qPCR: visualización y redes génicas, sustentadas en evidencia bibliográfica, para profundizar en la comprensión de las vías oncológicas"
)

# Cargar menú de configuración desde JSON
def load_menu() -> dict:
    default = {
        "version": 1,
        "menu": {
            "cancer_types": ["Breast Cancer", "Melanoma", "Colon Cancer"],
            "contexts": [
                {"key": "TEM", "label": "Cáncer y TEM"},
                {"key": "micro_rnas", "label": "Cáncer y micro RNAs"},
            ],
            "normalization_methods": [
                {"key": "reference_gene", "label": "gen de referencia"},
                {"key": "means", "label": "promedios"},
            ],
        },
    }
    cfg_path = Path(__file__).parent / "config" / "menu.json"
    try:
        with open(cfg_path, "r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        st.warning("No se pudo cargar web_app/config/menu.json; usando valores por defecto.")
        return default

MENU = load_menu()

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
    # Transformación numérica local (independiente del dominio)
    normalization = st.selectbox("Normalización numérica", options=["none", "zscore", "minmax"], index=0)
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
    # Configuración basada en JSON (no se lee del Excel)
    st.header("3) Configuración (desde JSON)")
    cancer_type = st.selectbox("Tipo de cáncer", options=MENU["menu"]["cancer_types"], index=0)
    context_labels = [c["label"] for c in MENU["menu"]["contexts"]]
    context_sel_label = st.selectbox("Contexto", options=context_labels, index=0)
    context_sel = next((c for c in MENU["menu"]["contexts"] if c["label"] == context_sel_label), MENU["menu"]["contexts"][0])
    norm_labels = [n["label"] for n in MENU["menu"]["normalization_methods"]]
    norm_sel_label = st.selectbox("Método de normalización (dominio)", options=norm_labels, index=0)
    norm_sel = next((n for n in MENU["menu"]["normalization_methods"] if n["label"] == norm_sel_label), MENU["menu"]["normalization_methods"][0])

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
        "context": {"key": context_sel.get("key"), "label": context_sel.get("label")},
        "cancer_type": cancer_type,
        "normalization_method": {"key": norm_sel.get("key"), "label": norm_sel.get("label")},
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
        # Sugerir prefijos/sufijos a partir de nombres de prueba detectados
        sample_names = (
            df_loaded.meta.get("sample_names")
            if (df_loaded is not None and df_loaded.meta)
            else [c for c in df_loaded.df.columns if c not in ("Well", "Target Name")]
        )
        aff = suggest_name_affixes(sample_names)

        st.markdown("### Sugerencias de nombres")
        col_a, col_b = st.columns(2)
        with col_a:
            st.caption("Posibles prefijos (frecuentes)")
            pref_opts = [f"{p} ({n})" for p, n in aff["prefixes"]]
            sel_ctrl_pref = st.multiselect("Prefijos controles", options=pref_opts, default=[])
            sel_samp_pref = st.multiselect("Prefijos muestras", options=pref_opts, default=[])
            def _clean(vals):
                return [v.split(' (')[0] for v in vals]
            ctrl_pref_vals = _clean(sel_ctrl_pref)
            samp_pref_vals = _clean(sel_samp_pref)
        with col_b:
            st.caption("Posibles sufijos (frecuentes)")
            suff_opts = [f"{s} ({n})" for s, n in aff["suffixes"]]
            sel_ctrl_suff = st.multiselect("Sufijos controles", options=suff_opts, default=[])
            sel_samp_suff = st.multiselect("Sufijos muestras", options=suff_opts, default=[])
            def _clean2(vals):
                return [v.split(' (')[0] for v in vals]
            ctrl_suff_vals = _clean2(sel_ctrl_suff)
            samp_suff_vals = _clean2(sel_samp_suff)

        # Combinar entradas manuales y sugeridas
        ctrl_prefixes = [p for p in [pref_ctrl] if p] + ctrl_pref_vals
        samp_prefixes = [p for p in [pref_samp] if p] + samp_pref_vals
        ctrl_suffixes = ctrl_suff_vals
        samp_suffixes = samp_suff_vals

        # Clasificar por prefijos y/o sufijos
        by_pref_ctrl, by_pref_samp = classify_by_prefixes(long_df, ctrl_prefixes, samp_prefixes)
        by_suf_ctrl, by_suf_samp = classify_by_suffixes(long_df, ctrl_suffixes, samp_suffixes)
        # Unir si se usan ambos métodos
        import pandas as _pd
        controles_df = _pd.concat([d for d in [by_pref_ctrl, by_suf_ctrl] if d is not None and not d.empty], ignore_index=True).drop_duplicates() if (not by_pref_ctrl.empty or not by_suf_ctrl.empty) else by_pref_ctrl
        muestras_df = _pd.concat([d for d in [by_pref_samp, by_suf_samp] if d is not None and not d.empty], ignore_index=True).drop_duplicates() if (not by_pref_samp.empty or not by_suf_samp.empty) else by_pref_samp
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
        "context": {"key": context_sel.get("key"), "label": context_sel.get("label")},
        "cancer_type": cancer_type,
        "normalization_method": {"key": norm_sel.get("key"), "label": norm_sel.get("label")},
        "prefijos": {"controles": pref_ctrl or None, "muestras": pref_samp or None},
        "archivo": df_loaded.source_name,
        "hoja": df_loaded.sheet_name,
        "normalizacion_numerica": normalization,
        "menu_version": MENU.get("version"),
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
