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
import streamlit as st

# Ensure project root is importable so `src.*` works when running from web_app/
_PROJ_ROOT = Path(__file__).resolve().parents[1]
if str(_PROJ_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJ_ROOT))

# Importamos funciones propias del proyecto
from src.core.io import LoadResult, list_excel_sheets, parse_qpcr_wide
from src.core.qpcr import (
    melt_wide_to_long,
    classify_tests,  # case-insensitive
    suggest_name_affixes,
)
from src.core.cleaning import drop_machine_controls
from src.core.fold_change import compute_fold_change
from src.core.tables import fc_comparison_table
from src.core.ensembl import add_ensembl_info_batch

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
        df_loaded = parse_qpcr_wide(uploaded, sheet_name=sheet)
        logger.info(f"Archivo cargado: name={df_loaded.source_name}, sheet={df_loaded.sheet_name}, shape={df_loaded.df.shape}")
        # Persistir en sesión para evitar perderlo al enviar el formulario
        st.session_state["df_loaded"] = df_loaded
    except Exception as e:
        st.error(f"Error al cargar el archivo: {e}")
        logger.exception("Fallo al cargar archivo")

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
                with st.spinner("Consultando Ensembl…"):
                    ensembl_df = add_ensembl_info_batch(df_to_annot, symbol_col='target', max_workers=3)
                desc_series = ensembl_df['description'].fillna('').astype(str).str.strip()
                ensembl_df['has_desc'] = desc_series.ne('') & desc_series.ne('No description')
                extras['ensembl_anotado.csv'] = ensembl_df.to_csv(index=False)

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
        except Exception as e:
            st.error(f"Error calculando Fold Change: {e}")

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
