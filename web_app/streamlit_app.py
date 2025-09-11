# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# CancerGeneSignatures - Web App
# -----------------------------------------------------------------------------
# Esta app de Streamlit implementa el flujo del notebook de análisis qPCR:
#   1) Carga de Excel (.xlsx/.xls)
#   2) Selección de parámetros desde menu.json
#   3) Filtrado de controles de máquina (PPC, RTC)
#   4) Selección de controles/muestras (tests o prefijos sugeridos)
#   5) Cálculo de Fold Change (promedios vs gen de referencia)
#   6) Visualizaciones y descarga de resultados
# -----------------------------------------------------------------------------

from __future__ import annotations

import io, os, json, sys
from pathlib import Path
from typing import Optional
import logging

import pandas as pd
import numpy as np
import streamlit as st

# Ensure project root is importable so `src.*` works when running from web_app/
_PROJ_ROOT = Path(__file__).resolve().parents[1]
if str(_PROJ_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJ_ROOT))

# Importamos funciones propias del proyecto
from src.core.io import LoadResult, list_excel_sheets, parse_qpcr_wide

# Hot-reload core IO to pick up signature changes during dev sessions
try:
    import importlib
    import src.core.io as _cgs_io
    _cgs_io = importlib.reload(_cgs_io)
    LoadResult = _cgs_io.LoadResult  # type: ignore
    list_excel_sheets = _cgs_io.list_excel_sheets  # type: ignore
    parse_qpcr_wide = _cgs_io.parse_qpcr_wide  # type: ignore
except Exception:
    pass
from src.core.qpcr import (
    melt_wide_to_long,
    classify_tests,  # case-insensitive
    suggest_name_affixes,
)
from src.core.cleaning import drop_machine_controls
from src.core.fold_change import compute_fold_change
from src.core.tables import fc_comparison_table
from src.core.ensembl import add_ensembl_info_batch
from src.core.string_enrichment import (
    run_string_enrichment,
    filter_enrichment,
    enrich_by_levels,
    dfs_to_excel_bytes,
)
from src.core.bibliography import (
    search_pubmed_by_genes,
    classify_bibliography,
    aggregate_counts_by_level_and_cancer,
)
from src.core.signatures import (
    create_signatures,
)

# Hot-reload STRING enrichment helpers during dev
try:
    import importlib as _importlib
    import src.core.string_enrichment as _cgs_str
    _cgs_str = _importlib.reload(_cgs_str)
    run_string_enrichment = _cgs_str.run_string_enrichment  # type: ignore
    filter_enrichment = _cgs_str.filter_enrichment  # type: ignore
    enrich_by_levels = _cgs_str.enrich_by_levels  # type: ignore
    dfs_to_excel_bytes = _cgs_str.dfs_to_excel_bytes  # type: ignore
except Exception:
    pass

# -----------------------------------------------------------------------------
# Configuración de la página Streamlit
# -----------------------------------------------------------------------------
st.set_page_config(
    page_title="CancerGeneSignatures - Web App",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.title("Análisis de datos de expresión diferencial por qPCR")
st.write(
    "De datos de qPCR a redes génicas sustentadas en evidencia bibliográfica "
    "para profundizar en la comprensión de las vías oncológicas."
)

# -----------------------------------------------------------------------------
# Funciones auxiliares
# -----------------------------------------------------------------------------
def _safe_concat(*dfs: Optional[pd.DataFrame]) -> pd.DataFrame:
    """Une varios DataFrames ignorando los que estén vacíos o None."""
    parts = [d for d in dfs if d is not None and isinstance(d, pd.DataFrame) and not d.empty]
    return pd.concat(parts, ignore_index=True).drop_duplicates() if parts else pd.DataFrame()

# -----------------------------------------------------------------------------
# Logging
# Configure a basic logger once; use env CGS_LOGLEVEL to override
if "_log_configured" not in st.session_state:
    level_name = os.getenv("CGS_LOGLEVEL", "INFO").upper()
    level = getattr(logging, level_name, logging.INFO)
    logging.basicConfig(level=level, format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
    st.session_state["_log_configured"] = True
logger = logging.getLogger("cgs.web_app")

 
# -----------------------------------------------------------------------------
# Menú de configuración (JSON con fallback)
# -----------------------------------------------------------------------------
def _default_menu() -> dict:
    return {
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

@st.cache_data
def load_menu() -> dict:
    env_path = os.getenv("CGS_MENU_PATH")
    if env_path:
        try:
            with open(env_path, "r", encoding="utf-8") as f:
                return json.load(f)
        except Exception:
            st.warning("No se pudo cargar menú desde CGS_MENU_PATH; usando valores por defecto.")
            return _default_menu()
    cfg_path = Path(__file__).parent / "config" / "menu.json"
    try:
        with open(cfg_path, "r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        st.warning("No se encontró config/menu.json; usando valores por defecto.")
        return _default_menu()

MENU = load_menu()

 

# -----------------------------------------------------------------------------
# Sidebar: parámetros de entrada y configuración
# -----------------------------------------------------------------------------
with st.sidebar:
    st.header("1) Datos de entrada (Excel qPCR)")
    uploaded = st.file_uploader("Archivo (.xlsx, .xls)", type=["xlsx", "xls"])
    sheet: Optional[str] = None
    if uploaded is not None:
        uploaded.seek(0)  # Reinicia el buffer
        if uploaded.name.lower().endswith((".xlsx", ".xls")):
            try:
                sheets = list_excel_sheets(uploaded)
                if sheets:
                    sheet = st.selectbox("Hoja de Excel", options=sheets, index=0)
            except Exception as e:
                st.warning(f"No se pudieron listar hojas: {e}")

    st.header("2) Parámetros del estudio")
    cancer_type = st.selectbox("Tipo de cáncer", MENU["menu"]["cancer_types"], index=0)
    context_sel_label = st.selectbox("Contexto", [c["label"] for c in MENU["menu"]["contexts"]], index=0)
    norm_sel_label = st.selectbox("Método preferido", [n["label"] for n in MENU["menu"]["normalization_methods"]], index=1)
    run_btn = st.button("Procesar archivo", type="primary")

# -----------------------------------------------------------------------------
# Carga de archivo y preprocesamiento
# -----------------------------------------------------------------------------
df_loaded: Optional[LoadResult] = st.session_state.get("df_loaded")
if uploaded is not None and run_btn:
    try:
        uploaded.seek(0)
        # Primero intentamos con coordenadas fijas (A4/B4)
        try:
            df_loaded = parse_qpcr_wide(
                uploaded,
                sheet_name=sheet,
                header_mode="coords",
                header_row_idx=3,  # A4/B4
                well_col_idx=0,
                target_col_idx=1,
            )
        except TypeError:
            # Compatibilidad con versiones antiguas sin nuevos parámetros
            df_loaded = parse_qpcr_wide(uploaded, sheet_name=sheet)
        logger.info(f"Archivo cargado: name={df_loaded.source_name}, sheet={df_loaded.sheet_name}, shape={df_loaded.df.shape}")
        # Persistir en sesión para evitar perderlo al enviar el formulario
        st.session_state["df_loaded"] = df_loaded
    except Exception as e:
        st.warning(f"Encabezado A4/B4 no válido o firma previa detectada ({e}). Probando detección automática…")
        logger.warning("Fallo modo coords; intentando auto")
        try:
            uploaded.seek(0)
            try:
                df_loaded = parse_qpcr_wide(uploaded, sheet_name=sheet, header_mode="auto")
            except TypeError:
                df_loaded = parse_qpcr_wide(uploaded, sheet_name=sheet)
            st.session_state["df_loaded"] = df_loaded
            logger.info("Carga con modo auto exitosa")
        except Exception as e2:
            st.error(f"Error al cargar el archivo: {e2}")
            logger.exception("Fallo al cargar archivo (auto)")

# -----------------------------------------------------------------------------
# Vista previa, análisis y gráficas
# -----------------------------------------------------------------------------
if df_loaded is not None:
    logger.debug("Usando df_loaded desde sesión para renderizado")
    st.subheader("Vista previa de datos (qPCR)")
    st.caption(f"Archivo: {df_loaded.source_name} | Hoja: {df_loaded.sheet_name or '-'} | Forma: {df_loaded.df.shape}")
    st.dataframe(df_loaded.df.head(20))

    # Mostrar parámetros elegidos desde menú
    c1, c2, c3 = st.columns(3)
    with c1:
        st.info(f"Contexto: {context_sel_label}")
    with c2:
        st.info(f"Método preferido: {norm_sel_label}")
    with c3:
        st.info(f"Tipo cáncer: {cancer_type}")

    # Construir largo y filtrar controles de máquina por defecto
    long_df = melt_wide_to_long(df_loaded.df)
    try:
        long_df = drop_machine_controls(long_df, column="target", controls=["PPC", "RTC"])  # como en el notebook
    except Exception as e:
        st.warning(f"No se pudieron filtrar controles de máquina: {e}")
        logger.warning(f"No se filtraron controles de máquina: {e}")

    # Resultados de la extracción (análogos al notebook)
    st.subheader("Resultados de la extracción")
    try:
        # Nombres de pruebas (desde meta o columnas)
        sample_names = []
        if isinstance(df_loaded.meta, dict):
            sample_names = df_loaded.meta.get('sample_names') or []
        if not sample_names:
            sample_names = [c for c in df_loaded.df.columns if c not in ("Well", "Target Name")]
        # Genes y pozos desde el dataframe ancho
        genes_col = df_loaded.df.get('Target Name')
        genes = [str(g).strip() for g in genes_col.dropna().unique().tolist()] if genes_col is not None else []
        wells_col = df_loaded.df.get('Well')
        pozos = [str(w).strip() for w in wells_col.dropna().unique().tolist()] if wells_col is not None else []

        st.markdown("- Nombres de las pruebas realizadas")
        tests_filtered = sample_names[1:] if len(sample_names) > 1 else sample_names
        st.info(f"Total {len(tests_filtered)}: {', '.join(tests_filtered)}")
        st.markdown("- Genes objetivo analizados")
        genes_filtered = [g for g in genes if g]
        st.info(f"Total {len(genes_filtered)}: {', '.join(genes_filtered)}")
        st.markdown("- Pozos detectados")
        st.info(f"Total {len(pozos)}: {', '.join(pozos)}")
        logger.debug(f"Extracción -> tests={len(sample_names)}, genes={len(genes_filtered)}, pozos={len(pozos)}")
    except Exception as e:
        st.warning(f"No se pudo mostrar el resumen de extracción: {e}")
        logger.warning(f"Fallo en resumen de extracción: {e}")

    # Clasificación con mejor UX: prefijos sugeridos o selección manual
    st.subheader("Clasificación de controles y muestras")
    st.caption("Usa prefijos sugeridos/manuales o selecciona las pruebas directamente.")

    file_key = f"assign_{df_loaded.source_name}:{df_loaded.sheet_name}"
    state = st.session_state.setdefault(file_key, {})

    unique_tests = sorted(long_df['test'].astype(str).dropna().unique().tolist())

    tab_pref, tab_select = st.tabs(["Por prefijos", "Selección manual"])

    with tab_pref:
        sugg = suggest_name_affixes(unique_tests, top_n=10)
        top_prefixes = [p for p, _ in (sugg.get('prefixes') or [])]
        colp1, colp2 = st.columns(2)
        with colp1:
            ctrl_suggest = st.selectbox("Sugerencia (controles)", options=["(vacío)"] + top_prefixes, index=0)
            ctrl_default = "" if ctrl_suggest == "(vacío)" else ctrl_suggest
            ctrl_prefix = st.text_input("Prefijo controles", value=state.get('ctrl_prefix', ctrl_default))
        with colp2:
            samp_suggest = st.selectbox("Sugerencia (muestras)", options=["(vacío)"] + top_prefixes, index=0)
            samp_default = "" if samp_suggest == "(vacío)" else samp_suggest
            samp_prefix = st.text_input("Prefijo muestras", value=state.get('samp_prefix', samp_default))
        submitted_pref = st.button("Clasificar por prefijos", type="primary")

    with tab_select:
        colm1, colm2 = st.columns(2)
        with colm1:
            selected_ctrl = st.multiselect("Pruebas de controles", options=unique_tests, default=state.get('selected_ctrl', []))
        with colm2:
            selected_samp = st.multiselect("Pruebas de muestras", options=unique_tests, default=state.get('selected_samp', []))
        submitted_sel = st.button("Clasificar por selección", type="secondary")

    clear_cls = st.button("Limpiar clasificación")
    if clear_cls:
        for k in ('ctrl_prefix','samp_prefix','selected_ctrl','selected_samp','controles_df','muestras_df'):
            state.pop(k, None)
        logger.info("Clasificación limpiada por el usuario")

    controles_df = state.get('controles_df', pd.DataFrame())
    muestras_df = state.get('muestras_df', pd.DataFrame())

    if submitted_pref:
        if not ctrl_prefix and not samp_prefix:
            st.warning("Debes ingresar al menos un prefijo (controles o muestras)")
        pref_ctrl_df, pref_samp_df = classify_tests(long_df, ctrl_prefix, samp_prefix)
        state['ctrl_prefix'] = ctrl_prefix
        state['samp_prefix'] = samp_prefix
        state['controles_df'] = pref_ctrl_df
        state['muestras_df'] = pref_samp_df
        controles_df = pref_ctrl_df
        muestras_df = pref_samp_df
        logger.info(f"Clasificación por prefijos: ctrl='{ctrl_prefix}' -> {len(controles_df)} filas, samp='{samp_prefix}' -> {len(muestras_df)} filas")
        if controles_df.empty or muestras_df.empty:
            st.warning("Alguna de las categorías resultó vacía. Revisa los prefijos o usa 'Selección manual'.")

    if submitted_sel:
        # Validar colisiones
        inter = set(selected_ctrl).intersection(set(selected_samp))
        if inter:
            st.warning(f"Hay pruebas en ambas categorías: {', '.join(sorted(inter))}")
        sel_ctrl_df = long_df[long_df['test'].astype(str).isin(selected_ctrl)] if selected_ctrl else long_df.iloc[0:0]
        sel_samp_df = long_df[long_df['test'].astype(str).isin(selected_samp)] if selected_samp else long_df.iloc[0:0]
        state['selected_ctrl'] = selected_ctrl
        state['selected_samp'] = selected_samp
        state['controles_df'] = sel_ctrl_df
        state['muestras_df'] = sel_samp_df
        controles_df = sel_ctrl_df
        muestras_df = sel_samp_df
        logger.info(f"Clasificación por selección: ctrl={len(selected_ctrl)} pruebas -> {len(controles_df)} filas, samp={len(selected_samp)} pruebas -> {len(muestras_df)} filas")
        if controles_df.empty or muestras_df.empty:
            st.warning("Alguna de las categorías quedó vacía. Selecciona al menos una prueba en cada lado.")

    st.write("Controles clasificados:", len(controles_df))
    st.write("Muestras clasificadas:", len(muestras_df))
    # Resumen de clasificación (análogos al notebook)
    if not controles_df.empty or not muestras_df.empty:
        st.subheader("Resumen de clasificación")
        for tipo, df_tipo in [("Controles", controles_df), ("Muestras", muestras_df)]:
            if not df_tipo.empty:
                uniq = df_tipo['test'].astype(str).unique().tolist()
                st.success(f"{tipo}: {len(uniq)} pruebas → {', '.join(uniq)}")
            else:
                st.warning(f"No se detectaron {tipo.lower()} con los prefijos actuales.")

    extras = {}
    # Imputación estilo notebook: NaN -> valor máximo global de CT
    if not controles_df.empty and not muestras_df.empty:
        import pandas as pd
        v_max = pd.concat([controles_df['ct'], muestras_df['ct']]).max()
        # Mensaje informativo del valor máximo usado (estilo notebook)
        max_str = f"{v_max:.2f}" if pd.notna(v_max) else "NaN"
        st.info(f"Valor máximo usado para imputación de Ct: {max_str}")
        controles_df['ct'] = pd.to_numeric(controles_df['ct'], errors='coerce').fillna(v_max)
        muestras_df['ct'] = pd.to_numeric(muestras_df['ct'], errors='coerce').fillna(v_max)

        # Guardar CSV limpios
        extras['controles_limpios.csv'] = controles_df.to_csv(index=False)
        extras['muestras_limpias.csv'] = muestras_df.to_csv(index=False)

        # Calcular Fold Change (promedio y gen de referencia)
        try:
            fc = compute_fold_change(controles_df, muestras_df)
            logger.info(f"FC consolidado shape: {fc.consolidated.shape}")
            st.subheader("Fold Change y ΔΔCt")
            m1, m2 = st.columns(2)
            with m1:
                st.metric("Gen de referencia (automático)", fc.reference_gene)
            with m2:
                st.caption("Elegido por menor desviación estándar promedio")

            with st.expander("Tabla consolidada", expanded=False):
                st.dataframe(fc.consolidated)
            st.plotly_chart(fc_comparison_table(fc.consolidated), use_container_width=True)

            import plotly.graph_objects as go
            fig = go.Figure()
            x_vals = fc.consolidated['target']
            fig.add_trace(go.Bar(x=x_vals, y=fc.consolidated['delta_delta_ct_promedio'], name='ΔΔCT (Promedios)', marker_color='#1f77b4', opacity=0.85, yaxis='y'))
            fig.add_trace(go.Bar(x=x_vals, y=fc.consolidated['delta_delta_ct_gen_ref'], name='ΔΔCT (Gen Ref)', marker_color='#ff7f0e', opacity=0.85, yaxis='y'))
            fig.add_trace(go.Scatter(x=x_vals, y=fc.consolidated['fold_change_promedio'], name='Fold Change (Promedios)', mode='markers+lines', marker=dict(color='#2ca02c', size=8, symbol='diamond'), line=dict(color='#2ca02c', width=2, dash='dot'), yaxis='y2'))
            fig.add_trace(go.Scatter(x=x_vals, y=fc.consolidated['fold_change_gen_ref'], name='Fold Change (Gen Ref)', mode='markers+lines', marker=dict(color='#d62728', size=8, symbol='diamond'), line=dict(color='#d62728', width=2, dash='dot'), yaxis='y2'))
            fig.update_layout(
                title=dict(text='Análisis comparativo de métodos de cálculo', x=0.5),
                template='plotly_white', barmode='group',
                yaxis=dict(title='ΔΔCT', showgrid=True, gridcolor='lightgray'),
                yaxis2=dict(title='Fold Change (log)', overlaying='y', side='right', type='log', showgrid=False),
                legend=dict(orientation='h', yanchor='bottom', y=1.02, xanchor='right', x=1),
                height=600, margin=dict(b=80, t=80, l=60, r=60)
            )
            st.plotly_chart(fig, use_container_width=True)

            # Clasificación por nivel de expresión (por método preferido del menú)
            st.subheader("Clasificación por nivel de expresión")
            default_idx = 1 if norm_sel_label == 'gen de referencia' else 0
            fc_source = st.radio("Fuente de Fold Change", ["promedios", "gen de referencia"], horizontal=True, index=default_idx)
            use_col = 'fold_change_promedio' if fc_source == 'promedios' else 'fold_change_gen_ref'
            df_expr = fc.consolidated[['target', use_col]].rename(columns={use_col: 'fold_change'}).copy()
            df_expr['nivel_expresion'] = pd.cut(
                df_expr['fold_change'],
                bins=[-float('inf'), 1.0, 2.0, float('inf')],
                labels=['subexpresado', 'estable', 'sobreexpresado'],
                right=False,
            )
            cexp1, cexp2 = st.columns(2)
            with cexp1:
                st.dataframe(df_expr)
            with cexp2:
                import plotly.express as px
                order_levels = ['estable', 'subexpresado', 'sobreexpresado']
                counts = df_expr['nivel_expresion'].value_counts().reindex(order_levels, fill_value=0)
                bar = px.bar(x=counts.index, y=counts.values, labels={'x': 'Nivel de expresión', 'y': 'Frecuencia'}, title='Distribución de niveles de expresión')
                st.plotly_chart(bar, use_container_width=True)

            extras['fold_change_consolidado.csv'] = fc.consolidated.to_csv(index=False)
            extras['expresion_categorizada.csv'] = df_expr.to_csv(index=False)

            # Anotación Ensembl (IDs y descripciones) sobre los targets clasificados
            st.subheader("Anotación Ensembl (IDs y descripciones)")
            st.caption("Consulta Ensembl para cada gen (requiere conexión a internet). Incluye exploración interactiva.")

            try:
                df_to_annot = df_expr[['target', 'nivel_expresion', 'fold_change']].drop_duplicates(subset=['target']).reset_index(drop=True)
                logger.info(f"Ensembl: anotando {len(df_to_annot)} genes (max_workers=3)")
                with st.spinner("Consultando Ensembl…"):
                    ensembl_df = add_ensembl_info_batch(df_to_annot, symbol_col='target', max_workers=3)
                desc_series = ensembl_df['description'].fillna('').astype(str).str.strip()
                ensembl_df['has_desc'] = desc_series.ne('') & desc_series.ne('No description')
                extras['ensembl_anotado.csv'] = ensembl_df.to_csv(index=False)
                try:
                    encontrados = int((ensembl_df['ensembl_id'] != 'Not found').sum())
                    logger.info(f"Ensembl: completado. IDs encontrados {encontrados}/{len(ensembl_df)}")
                except Exception:
                    pass

                tab_resumen, tab_explorar, tab_enlaces = st.tabs(["Resumen", "Explorar", "Enlaces"])

                # Resumen: métricas y gráfico de calidad de anotación por nivel de expresión
                with tab_resumen:
                    total = len(ensembl_df)
                    encontrados = int((ensembl_df['ensembl_id'] != 'Not found').sum())
                    con_desc = int((ensembl_df['has_desc']).sum())
                    m1, m2, m3 = st.columns(3)
                    m1.metric("Genes anotados (totales)", total)
                    m2.metric("Con Ensembl ID", encontrados)
                    m3.metric("Con descripción", con_desc)

                    import plotly.express as px
                    counts = (
                        ensembl_df.assign(desc=lambda d: d['has_desc'].map({True: 'con_descripción', False: 'sin_descripción'}))
                        .groupby(['nivel_expresion', 'desc']).size().reset_index(name='n')
                    )
                    order_levels = ['estable', 'subexpresado', 'sobreexpresado']
                    counts['nivel_expresion'] = pd.Categorical(counts['nivel_expresion'], categories=order_levels, ordered=True)
                    bar = px.bar(
                        counts.sort_values(['nivel_expresion','desc']),
                        x='nivel_expresion', y='n', color='desc', barmode='stack',
                        labels={'nivel_expresion': 'Nivel de expresión', 'n': 'Número de genes', 'desc': 'Descripción'},
                        title='Cobertura de anotación por nivel de expresión'
                    )
                    st.plotly_chart(bar, use_container_width=True)

                # Explorar: filtros por gen, descripción y nivel; descarga del subconjunto
                with tab_explorar:
                    f1, f2 = st.columns(2)
                    with f1:
                        q_gene = st.text_input("Filtrar por gen (contiene)", "").strip().lower()
                    with f2:
                        q_desc = st.text_input("Filtrar por descripción (contiene)", "").strip().lower()
                    order_levels = ['estable', 'subexpresado', 'sobreexpresado']
                    sel_levels = st.multiselect("Niveles", order_levels, default=order_levels)

                    filt = ensembl_df.copy()
                    if q_gene:
                        filt = filt[filt['target'].astype(str).str.lower().str.contains(q_gene, na=False)]
                    if q_desc:
                        filt = filt[filt['description'].astype(str).str.lower().str.contains(q_desc, na=False)]
                    if sel_levels:
                        filt = filt[filt['nivel_expresion'].isin(sel_levels)]

                    cols_show = ['target', 'nivel_expresion', 'fold_change', 'ensembl_id', 'description']
                    st.dataframe(filt[cols_show])

                    st.download_button(
                        label="Descargar resultado filtrado (CSV)",
                        data=filt[cols_show].to_csv(index=False),
                        file_name="ensembl_filtrado.csv",
                        mime="text/csv",
                        use_container_width=True,
                    )

                # Enlaces: listado con links a Ensembl y descripciones legibles
                with tab_enlaces:
                    st.caption("Navega por enlaces directos a Ensembl (máx. 100 primeros)")
                    subset = ensembl_df.copy().head(100)
                    for _, row in subset.iterrows():
                        gene = str(row['target'])
                        eid = str(row['ensembl_id'])
                        desc = str(row['description'])
                        lvl = str(row['nivel_expresion'])
                        url = f"https://www.ensembl.org/Homo_sapiens/Gene/Summary?g={eid}" if eid and eid != 'Not found' else None
                        if url:
                            st.markdown(f"- [{gene}]({url}) · {lvl} — {desc}")
                        else:
                            st.markdown(f"- {gene} · {lvl} — {desc}")

            except Exception as e:
                st.warning(f"No se pudo anotar con Ensembl: {e}")
            
            # Enriquecimiento funcional (STRING)
            st.subheader("Enriquecimiento funcional (STRING)")
            st.caption("Analiza términos GO/KEGG/Reactome enriquecidos por nivel de expresión. Requiere conexión a internet.")

            # Selección de niveles a procesar por separado
            order_levels = ['subexpresado', 'estable', 'sobreexpresado']
            sel_levels = st.multiselect(
                "Niveles a considerar",
                options=order_levels,
                default=['subexpresado', 'sobreexpresado']
            )

            # Fuentes/categorías a consultar en STRING
            cat_options = ["GO", "GO:BP", "GO:MF", "GO:CC", "KEGG", "Reactome"]
            sel_cats = st.multiselect(
                "Categorías (fuentes)",
                options=cat_options,
                default=["GO", "KEGG"]
            )

            cconf1, cconf2, cconf3, cconf4 = st.columns(4)
            with cconf1:
                max_fdr = st.number_input("FDR máx.", value=0.05, min_value=0.0, max_value=1.0, step=0.01, format="%.2f")
            with cconf2:
                min_size = st.number_input("Mín. genes por término", value=3, min_value=1, max_value=100, step=1)
            with cconf3:
                top_n = st.number_input("Top N", value=25, min_value=1, max_value=200, step=1)
            with cconf4:
                species = st.selectbox("Especie (NCBI Taxon)", options=[9606], index=0)

            run_enrich = st.button("Ejecutar enriquecimiento (STRING)", type="primary")

            if run_enrich:
                try:
                    # Ejecutar enriquecimiento por nivel y concatenar
                    levels_to_use = sel_levels or order_levels
                    base = df_expr[df_expr['nivel_expresion'].isin(levels_to_use)].copy()
                    logger.info(f"STRING: niveles={levels_to_use}, fuentes={sel_cats}, genes_totales={base['target'].nunique()}")
                    with st.spinner("Consultando STRING por nivel…"):
                        enr_dict = enrich_by_levels(
                            base,
                            symbol_col='target',
                            level_col='nivel_expresion',
                            levels=levels_to_use,
                            species=int(species),
                            sources=sel_cats,
                            caller_identity="UIMEO",
                        )
                    combined = enr_dict.get('combined')
                    by_level = enr_dict.get('by_level', {})
                    try:
                        logger.info(f"STRING: filas combinadas={0 if combined is None else len(combined)}")
                    except Exception:
                        pass

                    if combined is None or combined.empty:
                        st.info("Sin resultados de enriquecimiento para los parámetros actuales.")
                    else:
                        # Filtros comunes
                        combined_f = filter_enrichment(
                            combined,
                            include_categories=sel_cats or None,
                            max_fdr=float(max_fdr),
                            min_term_genes=int(min_size),
                            top_n=int(top_n),
                        )
                        st.dataframe(combined_f)

                        # Pestañas por nivel + resumen
                        tab_res, *tabs_levels = st.tabs(["Resumen"] + [lvl.capitalize() for lvl in levels_to_use])

                        # Resumen: barras -log10(FDR) top N combinadas
                        with tab_res:
                            try:
                                import numpy as np
                                import plotly.express as px
                                plot_df = combined_f.copy()
                                if not plot_df.empty and 'fdr' in plot_df.columns:
                                    plot_df["neglog10_fdr"] = -np.log10(plot_df["fdr"].clip(lower=1e-300))
                                    plot_df["label"] = plot_df.apply(lambda r: f"{r.get('term','')} ({r.get('category','')})", axis=1)
                                    fig = px.bar(
                                        plot_df.sort_values(["nivel_expresion", "neglog10_fdr"], ascending=[True, True]),
                                        x="neglog10_fdr",
                                        y="label",
                                        color="nivel_expresion",
                                        orientation="h",
                                        labels={"neglog10_fdr": "-log10(FDR)", "label": "Término"},
                                        title="Enriquecimiento combinado por nivel"
                                    )
                                    st.plotly_chart(fig, use_container_width=True)
                            except Exception:
                                pass

                        # Una pestaña por nivel con su propio gráfico y tabla
                        for lvl, tab in zip(levels_to_use, tabs_levels):
                            with tab:
                                lvl_df = by_level.get(lvl, pd.DataFrame())
                                if lvl_df is None or lvl_df.empty:
                                    st.info("Sin términos enriquecidos para este nivel.")
                                    continue
                                lvl_df_f = filter_enrichment(
                                    lvl_df,
                                    include_categories=sel_cats or None,
                                    max_fdr=float(max_fdr),
                                    min_term_genes=int(min_size),
                                    top_n=int(top_n),
                                )
                                st.dataframe(lvl_df_f)
                                try:
                                    import numpy as np
                                    import plotly.express as px
                                    p = lvl_df_f.copy()
                                    if not p.empty and 'fdr' in p.columns:
                                        p["neglog10_fdr"] = -np.log10(p["fdr"].clip(lower=1e-300))
                                        p["label"] = p.apply(lambda r: f"{r.get('term','')} ({r.get('category','')})", axis=1)
                                        fig_l = px.bar(
                                            p.sort_values("neglog10_fdr", ascending=True),
                                            x="neglog10_fdr",
                                            y="label",
                                            orientation="h",
                                            labels={"neglog10_fdr": "-log10(FDR)", "label": "Término"},
                                            title=f"Enriquecimiento ({lvl})"
                                        )
                                        st.plotly_chart(fig_l, use_container_width=True)
                                except Exception:
                                    pass

                        # Descargas: CSV combinado y Excel con hojas por nivel
                        st.download_button(
                            label="Descargar enriquecimiento combinado (CSV)",
                            data=combined_f.to_csv(index=False),
                            file_name="string_enrichment_combined.csv",
                            mime="text/csv",
                            use_container_width=True,
                        )

                        try:
                            sheets = [combined_f] + [by_level.get(lvl, pd.DataFrame()) for lvl in levels_to_use]
                            names = ["combinado"] + [f"{lvl}" for lvl in levels_to_use]
                            xlsx_bytes = dfs_to_excel_bytes(sheets, names)
                            st.download_button(
                                label="Descargar enriquecimiento (Excel, por nivel)",
                                data=xlsx_bytes,
                                file_name="string_enrichment_by_level.xlsx",
                                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                                use_container_width=True,
                            )
                        except Exception:
                            pass
                except Exception as e:
                    st.warning(f"No se pudo ejecutar el enriquecimiento: {e}")
        except Exception as e:
            st.error(f"Error calculando Fold Change: {e}")

        # Bibliografía (PubMed)
        st.subheader("Bibliografía (PubMed)")
        st.caption(
            "Busca artículos por gen y contexto (TEM o micro RNAs) en PubMed. "
            "Puedes ingresar tu email/API key aquí para la demo."
        )
        # Entrada directa (demo): email y API key NCBI
        env_email = os.getenv("NCBI_EMAIL", "")
        env_key = os.getenv("NCBI_API_KEY", "")
        ccreds1, ccreds2 = st.columns([2, 2])
        with ccreds1:
            ncbi_email_input = st.text_input("NCBI Email (obligatorio)", value=env_email, placeholder="tu_email@dominio.com")
        with ccreds2:
            ncbi_api_key_input = st.text_input("NCBI API Key (opcional)", value=env_key, placeholder="api_key")

        max_per_gene = st.number_input("Máximo de artículos por gen", value=100, min_value=10, max_value=300, step=10)
        run_pubmed = st.button("Buscar en PubMed", disabled=not bool(ncbi_email_input.strip()))

        if run_pubmed:
            try:
                # Usar el DataFrame anotado si está disponible; si no, usar df_expr
                if 'ensembl_df' in locals() and isinstance(ensembl_df, pd.DataFrame) and not ensembl_df.empty:
                    genes_df = ensembl_df[['target', 'ensembl_id', 'nivel_expresion']].drop_duplicates('target')
                else:
                    tmp = df_expr[['target', 'nivel_expresion']].drop_duplicates('target')
                    tmp['ensembl_id'] = ''
                    genes_df = tmp[['target', 'ensembl_id', 'nivel_expresion']]

                # Para simplificar la demo: setear variables de entorno en tiempo de ejecución
                os.environ["NCBI_EMAIL"] = ncbi_email_input.strip()
                if ncbi_api_key_input.strip():
                    os.environ["NCBI_API_KEY"] = ncbi_api_key_input.strip()
                logger.info(f"PubMed: consultando {len(genes_df)} genes (contexto={context_sel_label}, max_per_gene={int(max_per_gene)})")
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
                        logger=logger,
                    )
                prog.progress(100)
                status.success("Consulta PubMed finalizada")
                if bib is None or bib.empty:
                    st.info("No se encontraron artículos para los parámetros actuales.")
                else:
                    # Merge nivel_expresion
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

                    # Clasificación
                    st.markdown("### Clasificación por tipo de cáncer y contexto")
                    try:
                        classified = classify_bibliography(bib2)
                        try:
                            st.session_state["bibliografia_clasificada"] = classified.copy()
                        except Exception:
                            pass
                        st.dataframe(classified.head(50))
                        st.download_button(
                            label="Descargar bibliografía clasificada (CSV)",
                            data=classified.to_csv(index=False),
                            file_name="bibliografia_clasificada.csv",
                            mime="text/csv",
                            use_container_width=True,
                        )

                        # Gráfica de barras: número de estudios por tipo de cáncer y nivel
                        agg = aggregate_counts_by_level_and_cancer(classified)
                        if not agg.empty:
                            import plotly.express as px
                            order_levels = ['sobreexpresado', 'estable', 'subexpresado']
                            agg['nivel_expresion'] = pd.Categorical(agg['nivel_expresion'], categories=order_levels, ordered=True)
                            figb = px.bar(
                                agg.sort_values(['cancer_type', 'nivel_expresion']),
                                x='cancer_type', y='count', color='nivel_expresion',
                                barmode='group',
                                title=f"Número de estudios por tipo de cáncer y nivel de expresión ({context_sel_label})",
                                labels={'cancer_type': 'Tipo de cáncer', 'count': 'Número de estudios', 'nivel_expresion': 'Nivel de expresión'},
                            )
                            figb.update_layout(xaxis_tickangle=45)
                            st.plotly_chart(figb, use_container_width=True)

                            # Dispersión (burbujas) alternativa
                            figs = px.scatter(
                                agg,
                                x='count', y='cancer_type', size='count', color='nivel_expresion',
                                title=f"Artículos por tipo de cáncer y nivel ({context_sel_label})",
                                labels={'count': 'Número de artículos', 'cancer_type': 'Tipo de cáncer'},
                            )
                            st.plotly_chart(figs, use_container_width=True)

                    except Exception as e:
                        st.warning(f"No se pudo clasificar la bibliografía: {e}")

            except Exception as e:
                st.warning(f"No se pudo ejecutar la búsqueda en PubMed: {e}")

        # Botones de descarga de resultados (si hay datos)
        if extras:
            st.subheader("Descargas")
            for fname, data in extras.items():
                st.download_button(
                    label=f"Descargar {fname}",
                    data=data,
                    file_name=fname,
                    mime="text/csv",
                    use_container_width=True,
                )

    # (Opcional) Puedes exportar manualmente desde cada tabla mostrada en pantalla.

    # -------------------------------------------------------------------------
    # Firmas genéticas (a partir de bibliografía clasificada)
    # -------------------------------------------------------------------------
    st.markdown("---")
    st.header("Firmas genéticas")
    st.caption("Genera firmas por tipo de cáncer y nivel, con enriquecimiento de Hallmarks (MSigDB). Requiere gseapy y acceso a los GMT locales.")

    # Intentar recuperar bibliografía clasificada de esta sesión
    bib_class = st.session_state.get("bibliografia_clasificada")
    if 'classified' in locals() and isinstance(classified, pd.DataFrame) and not classified.empty:
        bib_class = classified.copy()
        st.session_state["bibliografia_clasificada"] = bib_class

    if bib_class is None or bib_class.empty:
        st.info("Ejecuta primero la sección 'Bibliografía (PubMed)' y clasificación para habilitar firmas. O bien, sube un CSV clasificado.")
        uploaded_bib = st.file_uploader("Opcional: subir bibliografía clasificada (CSV)", type=["csv"])
        if uploaded_bib is not None:
            try:
                bib_class = pd.read_csv(uploaded_bib)
                st.session_state["bibliografia_clasificada"] = bib_class
                st.success("Bibliografía cargada correctamente.")
            except Exception as e:
                st.error(f"No se pudo leer el CSV: {e}")

    if bib_class is not None and not bib_class.empty:
        c1, c2, c3 = st.columns([2, 2, 2])
        with c1:
            hall_gmt = st.text_input("Ruta GMT Hallmarks", value=str(_PROJ_ROOT / "gen-sets_GSEA_MSigDB/gsea_hallmarks_formatted.gmt"))
        with c2:
            back_gmt = st.text_input("Ruta GMT background (opcional)", value=str(_PROJ_ROOT / "gen-sets_GSEA_MSigDB/C5- ontology gene sets.gmt"))
        with c3:
            ctx = st.selectbox("Contexto biológico", ["Cáncer y TEM", "Cáncer y micro RNAs"], index=0)

        run_sig = st.button("Generar firmas")

        if run_sig:
            try:
                from src.core.signatures import HallmarkConfig
                cfg = HallmarkConfig(hallmark_gmt=hall_gmt, background_gmt=back_gmt)
                with st.spinner("Calculando firmas (incluye enriquecimiento de Hallmarks)…"):
                    df_sigs = create_signatures(bib_class, contexto_biologico=ctx, hallmark_cfg=cfg)
                if df_sigs is None or df_sigs.empty:
                    st.info("No se generaron firmas para los datos disponibles.")
                else:
                    st.success(f"Firmas generadas: {len(df_sigs)} filas")
                    # Vista segura para Streamlit/Arrow: convertir columnas *_genes (listas) a string
                    df_sigs_display = df_sigs.copy()
                    try:
                        list_cols = [c for c in df_sigs_display.columns if c.startswith('hallmark_') and c.endswith('_genes')]
                        for c in list_cols:
                            df_sigs_display[c] = df_sigs_display[c].apply(lambda v: ", ".join(v) if isinstance(v, list) else (str(v) if pd.notna(v) else ""))
                    except Exception:
                        pass
                    st.dataframe(df_sigs_display)
                    # Descarga CSV
                    st.download_button(
                        label="Descargar firmas (CSV)",
                        data=df_sigs.to_csv(index=False),
                        file_name="firmas_geneticas.csv",
                        mime="text/csv",
                        use_container_width=True,
                    )

                    # Visualización rápida: Sunburst por tipo seleccionado
                    tipos = sorted(df_sigs['cancer_type'].dropna().astype(str).unique().tolist())
                    sel_tipo = st.selectbox("Tipo de cáncer para visualizar", tipos, index=0 if tipos else None)
                    if sel_tipo:
                        try:
                            import plotly.express as px
                            # construir registros planos: por nivel, gene y hallmark presentes
                            recs = []
                            for _, row in df_sigs[df_sigs['cancer_type'] == sel_tipo].iterrows():
                                nivel = row.get('nivel_expresion')
                                genes_firma = row.get('genes', []) if isinstance(row.get('genes'), list) else []
                                counts_firma = row.get('conteo_articulos_por_gene', []) if isinstance(row.get('conteo_articulos_por_gene'), list) else []
                                gene_to_count = dict(zip(genes_firma, counts_firma))
                                # detectar columnas hallmark
                                for col in row.index:
                                    if col.startswith('hallmark_') and col.endswith('_genes'):
                                        term = col[len('hallmark_'):-len('_genes')]
                                        genes_hm = row[col] if isinstance(row[col], list) else []
                                        pval_col = f"hallmark_{term}_pvalue"
                                        pval = row.get(pval_col, None)
                                        for g in genes_hm:
                                            recs.append({
                                                'nivel_expresion': nivel,
                                                'gene': g,
                                                'hallmark': term.replace('HALLMARK_', '').replace('_', ' '),
                                                'articles': gene_to_count.get(g, 0),
                                                'pvalue': pval if pd.notna(pval) else 1.0,
                                            })
                            if recs:
                                flat = pd.DataFrame(recs)
                                flat['log_p'] = -np.log10(flat['pvalue'].replace(0, 1e-300))
                                tabs = st.tabs([f"{lvl}" for lvl in sorted(flat['nivel_expresion'].dropna().unique().tolist())])
                                for lvl, tab in zip(sorted(flat['nivel_expresion'].dropna().unique().tolist()), tabs):
                                    with tab:
                                        d = flat[flat['nivel_expresion'] == lvl]
                                        if d.empty:
                                            st.info("Sin datos para este nivel.")
                                            continue
                                        # Sunburst gene -> hallmark
                                        fig = px.sunburst(
                                            d,
                                            path=['gene', 'hallmark'],
                                            values='articles',
                                            color='log_p',
                                            color_continuous_scale='RdBu_r',
                                            title=f"Firmas - {sel_tipo} - {lvl}",
                                        )
                                        st.plotly_chart(fig, use_container_width=True)
                        except Exception as e:
                            st.warning(f"No se pudo renderizar la visualización de firmas: {e}")
            except Exception as e:
                st.error(f"No se pudieron generar las firmas: {e}")


# -----------------------------------------------------------------------------
# Guía rápida
# -----------------------------------------------------------------------------
st.markdown("---")
st.markdown("""
**Cómo usar la aplicación**
1. Sube un archivo Excel (.xlsx/.xls) con formato qPCR y elige la hoja.
2. Selecciona parámetros del estudio desde el menú (contexto, tipo de cáncer, método preferido).
3. Ingresa manualmente los prefijos de controles y muestras y presiona "Clasificar y calcular".
4. Revisa Fold Change (promedios vs gen de referencia), gráficas y la clasificación por nivel de expresión.
5. Exporta desde los widgets (descarga en cada tabla/gráfico) según necesidad.
""")
