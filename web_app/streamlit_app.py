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

import io, os, sys, re
from pathlib import Path
from typing import Any, Optional
import logging

import pandas as pd
import numpy as np
import streamlit as st

from app.config.loader import load_app_config, ConfigError

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

# -----------------------------------------------------------------------------
# Configuración centralizada de la aplicación
# -----------------------------------------------------------------------------
_config_warnings: list[str] = []
try:
    secrets_source = st.secrets  # type: ignore[attr-defined]
except Exception:
    secrets_source = None
try:
    APP_CONFIG = load_app_config(
        secrets=secrets_source if secrets_source is not None else None,
        warn=_config_warnings.append,
    )
except ConfigError as exc:
    st.error(f"Error cargando configuración: {exc}")
    st.stop()

MENU = APP_CONFIG.menu
SERVICES = APP_CONFIG.services
from src.core.qpcr import (
    suggest_name_affixes,
)
from src.core.tables import fc_comparison_table
from src.core.ensembl import add_ensembl_info_batch
from src.core.string_enrichment import dfs_to_excel_bytes
# Hot-reload bibliography helpers during dev (ensure new functions are visible)
try:
    import importlib as _importlib
    import src.core.bibliography as _cgs_bib
    _cgs_bib = _importlib.reload(_cgs_bib)
    search_pubmed_by_genes = _cgs_bib.search_pubmed_by_genes  # type: ignore
    classify_bibliography = _cgs_bib.classify_bibliography  # type: ignore
    aggregate_counts_by_level_and_cancer = _cgs_bib.aggregate_counts_by_level_and_cancer  # type: ignore
    filter_bibliography_by_cancer = getattr(_cgs_bib, 'filter_bibliography_by_cancer')  # type: ignore
    interpret_gene_relations = getattr(_cgs_bib, 'interpret_gene_relations')  # type: ignore
    summarize_relations_by_gene = getattr(_cgs_bib, 'summarize_relations_by_gene')  # type: ignore
except Exception:
    # Fallbacks (should not be needed if module is available)
    from src.core.bibliography import (
        search_pubmed_by_genes,
        classify_bibliography,
        aggregate_counts_by_level_and_cancer,
    )
from src.core.signatures import (
    create_signatures,
)
from src.core.visuals import (
    hallmarks_polar_chart,
    fingerprint_heatmap,
    clustered_expression_by_level,
)
from src.integrations.google_nlp import GoogleNLPClient, aggregate_insights
from app.state import AppSessionState
from app.services.qpcr import (
    build_long_table,
    summarize_extraction,
    classify_tests_by_prefixes,
    classify_tests_by_suffixes,
    classify_tests_by_regex,
    classify_tests_by_selection,
    detect_collisions,
    apply_collision_strategy,
)
from app.services.fold_change import (
    apply_undetermined_policy,
    compute_quality_metrics,
    compute_fold_change_with_expression,
    FoldChangePreparationError,
    QualityMetrics,
)
from app.services.visuals import (
    build_fc_detail_figure,
    build_expression_distribution,
    build_expression_treemap,
)
from app.services.bibliography import (
    PubMedRequest,
    fetch_pubmed_articles,
    merge_expression_levels,
)

# Cache de firmas para evitar recomputar al cambiar solo la selección visual
import hashlib
@st.cache_data(show_spinner=False)
def _compute_signatures_cached(bib_csv: str, ctx: str, hall_gmt: str, back_gmt: str, sel_tipo: str) -> pd.DataFrame:
    try:
        df_bib = pd.read_csv(io.StringIO(bib_csv))
    except Exception:
        return pd.DataFrame()
    from src.core.signatures import HallmarkConfig
    cfg = HallmarkConfig(hallmark_gmt=hall_gmt, background_gmt=back_gmt)
    return create_signatures(df_bib, contexto_biologico=ctx, hallmark_cfg=cfg)

# Cache de anotación Ensembl (evitar llamadas repetidas a red en cada rerun)
@st.cache_data(show_spinner=False, ttl=3600)
def _annotate_ensembl_cached(df_csv: str, max_workers: int = 3) -> pd.DataFrame:
    try:
        df_in = pd.read_csv(io.StringIO(df_csv))
    except Exception:
        return pd.DataFrame()
    return add_ensembl_info_batch(df_in, symbol_col='target', max_workers=max_workers)

@st.cache_data(show_spinner=False, ttl=3600)
def _string_enrichment_cached(
    df_expr_csv: str,
    levels: list[str],
    species: int,
    sources: list[str],
    max_fdr: float,
    min_size: int,
    top_n: int,
) -> dict:
    try:
        df_expr = pd.read_csv(io.StringIO(df_expr_csv))
    except Exception:
        return {"combined": pd.DataFrame(), "by_level": {}}
    result = perform_string_enrichment(
        df_expr,
        levels=levels,
        species=species,
        sources=sources,
        max_fdr=max_fdr,
        min_term_genes=min_size,
        top_n=top_n,
    )
    return {
        "combined": result.combined,
        "by_level": result.by_level,
    }

@st.cache_data(show_spinner=False, ttl=86400)
def _pubmed_cached(genes_df_csv: str, context: str, max_per_gene: int, email: str, api_key: str|None) -> pd.DataFrame:
    try:
        genes_df = pd.read_csv(io.StringIO(genes_df_csv))
    except Exception:
        return pd.DataFrame()
    request = PubMedRequest(
        context=context,
        max_per_gene=int(max_per_gene),
        email=email,
        api_key=api_key,
    )
    return fetch_pubmed_articles(genes_df, request)

# -----------------------------------------------------------------------------
# Configuración de la página Streamlit
# -----------------------------------------------------------------------------
st.set_page_config(
    page_title=APP_CONFIG.project_name,
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

for _cfg_warning in _config_warnings:
    st.warning(_cfg_warning)

app_state = AppSessionState.load()

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
# Configure a basic logger once con el nivel definido en AppConfig
if "_log_configured" not in st.session_state:
    level_name = APP_CONFIG.log_level.upper() if APP_CONFIG.log_level else "INFO"
    level = getattr(logging, level_name, logging.INFO)
    logging.basicConfig(level=level, format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
    st.session_state["_log_configured"] = True
logger = logging.getLogger("cgs.web_app")

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
    cancer_type = st.selectbox("Tipo de cáncer", MENU.cancer_types, index=0)
    context_options = [ctx.label for ctx in MENU.contexts]
    context_default_idx = context_options.index(app_state.context_label) if app_state.context_label in context_options else 0
    context_sel_label = st.selectbox("Contexto", context_options, index=context_default_idx)
    norm_options = [method.label for method in MENU.normalization_methods]
    norm_default_idx = 1 if len(norm_options) > 1 else 0
    norm_sel_label = st.selectbox("Método preferido", norm_options, index=norm_default_idx)

    st.header("3) Política 'Undetermined/ND'")
    und_options = ["nan", "ctmax", "value"]
    und_default_idx = und_options.index(app_state.undetermined.policy) if app_state.undetermined.policy in und_options else 0
    und_policy = st.selectbox(
        "Cómo tratar valores 'Undetermined' (ND)",
        options=und_options,
        index=und_default_idx,
        help="\n- nan: deja como NaN (se puede imputar después).\n- ctmax: usa el Ct máximo observado por columna.\n- value: usa un valor fijo (p. ej., 40).",
    )
    und_value = st.number_input(
        "Valor fijo para 'value'",
        value=float(app_state.undetermined.value),
        min_value=0.0,
        max_value=100.0,
        step=0.5,
    )

    # Preferencia global de gráficos: excluir genes 'estables'
    exclude_stable = st.checkbox(
        "Excluir genes 'estables' en gráficos",
        value=bool(app_state.exclude_stable),
        help="Cuando está activado, los gráficos omiten los genes clasificados como 'estables'.",
    )
    run_btn = st.button("Procesar archivo", type="primary")

app_state.undetermined.policy = und_policy
app_state.undetermined.value = float(und_value)
app_state.context_label = context_sel_label
app_state.exclude_stable = bool(exclude_stable)
app_state.persist()

# -----------------------------------------------------------------------------
# Carga de archivo y preprocesamiento
# -----------------------------------------------------------------------------
df_loaded: Optional[LoadResult] = app_state.df_loaded
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
                undetermined_policy=und_policy,
                undetermined_value=float(und_value),
            )
        except TypeError:
            # Compatibilidad con versiones antiguas sin nuevos parámetros
            df_loaded = parse_qpcr_wide(uploaded, sheet_name=sheet)
        logger.info(f"Archivo cargado: name={df_loaded.source_name}, sheet={df_loaded.sheet_name}, shape={df_loaded.df.shape}")
        # Persistir en sesión para evitar perderlo al enviar el formulario
        app_state.df_loaded = df_loaded
        app_state.persist()
    except Exception as e:
        st.warning(f"Encabezado A4/B4 no válido o firma previa detectada ({e}). Probando detección automática…")
        logger.warning("Fallo modo coords; intentando auto")
        try:
            uploaded.seek(0)
            try:
                df_loaded = parse_qpcr_wide(
                    uploaded,
                    sheet_name=sheet,
                    header_mode="auto",
                    undetermined_policy=und_policy,
                    undetermined_value=float(und_value),
                )
            except TypeError:
                df_loaded = parse_qpcr_wide(uploaded, sheet_name=sheet)
            app_state.df_loaded = df_loaded
            app_state.persist()
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

    # Construir dataframe largo y filtrar controles de máquina por defecto
    long_df, controls_warning = build_long_table(df_loaded)
    if controls_warning:
        st.warning(f"No se pudieron filtrar controles de máquina: {controls_warning}")
        logger.warning(f"No se filtraron controles de máquina: {controls_warning}")

    # Resultados de la extracción (análogos al notebook)
    st.subheader("Resultados de la extracción")
    try:
        extraction = summarize_extraction(df_loaded)
        st.markdown("- Nombres de las pruebas realizadas")
        st.info(f"Total {len(extraction.sample_names)}: {', '.join(extraction.sample_names)}")
        st.markdown("- Genes objetivo analizados")
        st.info(f"Total {len(extraction.genes)}: {', '.join(extraction.genes)}")
        st.markdown("- Pozos detectados")
        st.info(f"Total {len(extraction.wells)}: {', '.join(extraction.wells)}")
        logger.debug(
            "Extracción -> tests=%s, genes=%s, pozos=%s",
            len(extraction.sample_names),
            len(extraction.genes),
            len(extraction.wells),
        )
    except Exception as e:
        st.warning(f"No se pudo mostrar el resumen de extracción: {e}")
        logger.warning(f"Fallo en resumen de extracción: {e}")

    # Clasificación con mejor UX: prefijos, sufijos, regex o selección manual
    st.subheader("Clasificación de controles y muestras")
    st.caption("Clasifica por prefijos/sufijos/regex o selecciona directamente.")

    file_key = f"assign_{df_loaded.source_name}:{df_loaded.sheet_name}"
    state = st.session_state.setdefault(file_key, {})

    unique_tests = sorted(long_df['test'].astype(str).dropna().unique().tolist())

    tab_pref, tab_suff, tab_regex, tab_select = st.tabs(["Por prefijos", "Por sufijos", "Por regex", "Selección manual"])

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
        # Vista previa
        prev_ctrl = [t for t in unique_tests if ctrl_prefix and str(t).startswith(ctrl_prefix)]
        prev_samp = [t for t in unique_tests if samp_prefix and str(t).startswith(samp_prefix)]
        st.caption(f"Previa → Controles: {len(prev_ctrl)} | Muestras: {len(prev_samp)}")
        if prev_ctrl:
            st.caption("Controles: " + ", ".join(prev_ctrl[:20]) + (" …" if len(prev_ctrl) > 20 else ""))
        if prev_samp:
            st.caption("Muestras: " + ", ".join(prev_samp[:20]) + (" …" if len(prev_samp) > 20 else ""))
        submitted_pref = st.button("Clasificar por prefijos", type="primary")

    with tab_suff:
        sugg = suggest_name_affixes(unique_tests, top_n=10)
        top_suffixes = [s for s, _ in (sugg.get('suffixes') or [])]
        cols1, cols2 = st.columns(2)
        with cols1:
            ctrl_suff = st.text_input("Sufijos controles (coma-sep.)", value=state.get('ctrl_suff', ''))
        with cols2:
            samp_suff = st.text_input("Sufijos muestras (coma-sep.)", value=state.get('samp_suff', ''))
        st.caption("Sugerencias sufijos: " + ", ".join(top_suffixes))
        suff_ctrl_list = [s.strip() for s in ctrl_suff.split(',') if s.strip()]
        suff_samp_list = [s.strip() for s in samp_suff.split(',') if s.strip()]
        prev_ctrl = [t for t in unique_tests if any(str(t).endswith(s) for s in suff_ctrl_list)]
        prev_samp = [t for t in unique_tests if any(str(t).endswith(s) for s in suff_samp_list)]
        st.caption(f"Previa → Controles: {len(prev_ctrl)} | Muestras: {len(prev_samp)}")
        if prev_ctrl:
            st.caption("Controles: " + ", ".join(prev_ctrl[:20]) + (" …" if len(prev_ctrl) > 20 else ""))
        if prev_samp:
            st.caption("Muestras: " + ", ".join(prev_samp[:20]) + (" …" if len(prev_samp) > 20 else ""))
        submitted_suff = st.button("Clasificar por sufijos", type="primary")

    with tab_regex:
        colr1, colr2 = st.columns(2)
        with colr1:
            ctrl_re = st.text_input("Regex controles", value=state.get('ctrl_re', ''))
        with colr2:
            samp_re = st.text_input("Regex muestras", value=state.get('samp_re', ''))
        try:
            prev_ctrl = [t for t in unique_tests if (ctrl_re and re.search(ctrl_re, str(t)) is not None)]
        except Exception:
            prev_ctrl = []
        try:
            prev_samp = [t for t in unique_tests if (samp_re and re.search(samp_re, str(t)) is not None)]
        except Exception:
            prev_samp = []
        st.caption(f"Previa → Controles: {len(prev_ctrl)} | Muestras: {len(prev_samp)}")
        if prev_ctrl:
            st.caption("Controles: " + ", ".join(prev_ctrl[:20]) + (" …" if len(prev_ctrl) > 20 else ""))
        if prev_samp:
            st.caption("Muestras: " + ", ".join(prev_samp[:20]) + (" …" if len(prev_samp) > 20 else ""))
        submitted_regex = st.button("Clasificar por regex", type="primary")

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

    # Utilidad local para resolver colisiones
    def _resolve_collisions(ctrl_df: pd.DataFrame, samp_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
        collisions = detect_collisions(ctrl_df, samp_df)
        if not collisions:
            return ctrl_df, samp_df
        preview = ', '.join(collisions[:10]) + (' …' if len(collisions) > 10 else '')
        st.warning(f"Colisiones: {len(collisions)} pruebas aparecen en ambos grupos → {preview}")
        choice = st.radio(
            "Resolver colisiones",
            ["priorizar controles", "priorizar muestras", "excluir colisiones"],
            horizontal=True,
            index=0,
        )
        return apply_collision_strategy(ctrl_df, samp_df, collisions, choice)

    if submitted_pref:
        if not ctrl_prefix and not samp_prefix:
            st.warning("Debes ingresar al menos un prefijo (controles o muestras)")
        result = classify_tests_by_prefixes(long_df, [ctrl_prefix], [samp_prefix])
        pref_ctrl_df, pref_samp_df = _resolve_collisions(result.controles, result.muestras)
        state['ctrl_prefix'] = ctrl_prefix
        state['samp_prefix'] = samp_prefix
        state['controles_df'] = pref_ctrl_df
        state['muestras_df'] = pref_samp_df
        controles_df = pref_ctrl_df
        muestras_df = pref_samp_df
        logger.info(f"Clasificación por prefijos: ctrl='{ctrl_prefix}' -> {len(controles_df)} filas, samp='{samp_prefix}' -> {len(muestras_df)} filas")
        if controles_df.empty or muestras_df.empty:
            st.warning("Alguna de las categorías resultó vacía. Revisa los prefijos o usa 'Selección manual'.")

    if 'submitted_suff' in locals() and submitted_suff:
        result = classify_tests_by_suffixes(long_df, suff_ctrl_list, suff_samp_list)
        pref_ctrl_df, pref_samp_df = _resolve_collisions(result.controles, result.muestras)
        state['ctrl_suff'] = ctrl_suff
        state['samp_suff'] = samp_suff
        state['controles_df'] = pref_ctrl_df
        state['muestras_df'] = pref_samp_df
        controles_df = pref_ctrl_df
        muestras_df = pref_samp_df
        logger.info(f"Clasificación por sufijos: ctrl={suff_ctrl_list} -> {len(controles_df)} filas, samp={suff_samp_list} -> {len(muestras_df)} filas")
        if controles_df.empty or muestras_df.empty:
            st.warning("Alguna de las categorías resultó vacía. Revisa los sufijos o usa otra pestaña.")

    if 'submitted_regex' in locals() and submitted_regex:
        try:
            result = classify_tests_by_regex(long_df, ctrl_re, samp_re)
            pref_ctrl_df, pref_samp_df = _resolve_collisions(result.controles, result.muestras)
            state['ctrl_re'] = ctrl_re
            state['samp_re'] = samp_re
            state['controles_df'] = pref_ctrl_df
            state['muestras_df'] = pref_samp_df
            controles_df = pref_ctrl_df
            muestras_df = pref_samp_df
            logger.info(f"Clasificación por regex: ctrl='{ctrl_re}' -> {len(controles_df)} filas, samp='{samp_re}' -> {len(muestras_df)} filas")
            if controles_df.empty or muestras_df.empty:
                st.warning("Alguna de las categorías resultó vacía. Revisa los patrones o usa otra pestaña.")
        except re.error as ex:
            st.error(f"Regex inválido: {ex}")

    if submitted_sel:
        result = classify_tests_by_selection(long_df, selected_ctrl, selected_samp)
        collisions = detect_collisions(result.controles, result.muestras)
        if collisions:
            preview = ', '.join(collisions[:10]) + (' …' if len(collisions) > 10 else '')
            st.warning(f"Hay pruebas en ambas categorías: {preview}")
        sel_ctrl_df, sel_samp_df = result.controles, result.muestras
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

    extras: dict[str, str] = {}
    ok_for_fc = False
    ready_for_fc = False
    metrics: Optional[QualityMetrics] = None
    imputation_output = None

    if not controles_df.empty or not muestras_df.empty:
        st.subheader("Resumen de clasificación")
        ok_for_fc = True
        for tipo, df_tipo in [("Controles", controles_df), ("Muestras", muestras_df)]:
            if not df_tipo.empty:
                uniq = df_tipo['test'].astype(str).unique().tolist()
                st.success(f"{tipo}: {len(uniq)} pruebas → {', '.join(uniq)}")
            else:
                st.warning(f"No se detectaron {tipo.lower()} con los criterios actuales.")

        impute_by_cols = [c for c in ['plate_id', 'target'] if c in controles_df.columns and c in muestras_df.columns]
        try:
            imputation_output = apply_undetermined_policy(
                controles_df,
                muestras_df,
                policy=policy_applied,
                fixed_value=float(app_state.undetermined.value),
                group_columns=impute_by_cols or None,
            )
        except FoldChangePreparationError as exc:
            st.error(f"Imputación fallida: {exc}")
            controles_df = pd.DataFrame()
            muestras_df = pd.DataFrame()
            ok_for_fc = False
        else:
            controles_df = imputation_output.controles
            muestras_df = imputation_output.muestras
            if imputation_output.message:
                if imputation_output.policy == 'nan':
                    st.info(imputation_output.message)
                else:
                    st.success(imputation_output.message)

            if not imputation_output.summary.empty:
                counts = imputation_output.summary['grupo'].value_counts()
                st.metric("Total imputaciones", len(imputation_output.summary))
                st.caption(
                    f"Controles imputados: {int(counts.get('control', 0))} · "
                    f"Muestras imputadas: {int(counts.get('muestra', 0))}"
                )
                cols_show = [
                    c for c in ['grupo', 'test', 'target', 'plate_id', 'well', 'ct', 'ct_imputed']
                    if c in imputation_output.summary.columns
                ]
                st.dataframe(imputation_output.summary[cols_show])
                extras['ct_imputaciones.csv'] = imputation_output.summary.to_csv(index=False)
                st.download_button(
                    label="Descargar detalle de imputaciones",
                    data=imputation_output.summary.to_csv(index=False),
                    file_name="detalle_imputaciones.csv",
                    mime="text/csv",
                    use_container_width=True,
                )
            elif imputation_output and imputation_output.policy != 'nan':
                st.info("No se detectaron filas que requieran imputación con la política seleccionada.")

        if not controles_df.empty and not muestras_df.empty:
            metrics = compute_quality_metrics(controles_df, muestras_df)
            st.subheader("Calidad de datos (post-imputación)")
            mqa1, mqa2, mqa3, mqa4 = st.columns(4)
            with mqa1:
                st.metric("Genes (controles)", len(metrics.ctrl_targets))
            with mqa2:
                st.metric("Genes (muestras)", len(metrics.sample_targets))
            with mqa3:
                st.metric("Genes en común", len(metrics.common_targets))
            with mqa4:
                st.metric(
                    "NaN ct (ctrl/mues)",
                    f"{metrics.ctrl_nan_ratio:.0%} / {metrics.sample_nan_ratio:.0%}"
                )

            if not metrics.common_targets:
                st.error("No hay intersección de genes entre controles y muestras. Ajusta tu clasificación.")
                ok_for_fc = False
            if metrics.ctrl_nan_ratio >= 1.0 or metrics.sample_nan_ratio >= 1.0:
                st.error("Todos los valores Ct son NaN en algún grupo. Revisa la política de 'Undetermined' o tus datos.")
                ok_for_fc = False

            n_min = st.number_input(
                "Mínimo de réplicas por gen y grupo",
                min_value=1,
                max_value=10,
                value=1,
                step=1,
            )
            eligible = metrics.eligible_genes(int(n_min)) if metrics else []
            st.info(f"Genes que cumplen n≥{int(n_min)} en ambos grupos: {len(eligible)}")
            if not eligible:
                st.warning("Ningún gen cumple el mínimo de réplicas en ambos grupos. Considera reducir el umbral o revisar datos.")

            if imputation_output and imputation_output.has_missing:
                st.warning(
                    "Persisten valores Ct faltantes. Ajusta la política en la barra lateral y vuelve a procesar el archivo."
                )
                ready_for_fc = False
            else:
                ready_for_fc = ok_for_fc and bool(metrics.common_targets)

    ready_for_fc = ready_for_fc and not controles_df.empty and not muestras_df.empty

    if not controles_df.empty and not muestras_df.empty and ok_for_fc:
        controles_df = controles_df.copy()
        muestras_df = muestras_df.copy()
        controles_df.loc[:, 'ct'] = pd.to_numeric(controles_df['ct'], errors='coerce')
        muestras_df.loc[:, 'ct'] = pd.to_numeric(muestras_df['ct'], errors='coerce')

        st.subheader("Imputación de Ct (resumen)")
        policy_applied = app_state.undetermined.policy
        st.caption(f"Política seleccionada: {policy_applied}")

        if not controles_df.empty and not muestras_df.empty:
            extras['controles_limpios.csv'] = controles_df.to_csv(index=False)
            extras['muestras_limpias.csv'] = muestras_df.to_csv(index=False)

    if not ready_for_fc and not controles_df.empty and not muestras_df.empty:
        st.stop()

    if ready_for_fc:
        try:
            fc_result = compute_fold_change_with_expression(controles_df, muestras_df)
            fc = fc_result.computation
            df_expr = fc_result.expression_table
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

            exclude_stable = bool(app_state.exclude_stable)
            y2_scale = st.radio("Escala FC", ["log", "lineal"], horizontal=True, index=0)
            fc_figure = build_fc_detail_figure(
                fc.consolidated,
                df_expr,
                exclude_stable=exclude_stable,
                reference_gene=fc.reference_gene,
                scale=y2_scale,
            )
            st.plotly_chart(fc_figure, use_container_width=True)

            # Clasificación por nivel de expresión (por método preferido del menú)
            st.subheader("Clasificación por nivel de expresión")
            default_idx = 1 if norm_sel_label == 'gen de referencia' else 0
            fc_source = st.radio("Fuente de Fold Change", ["promedios", "gen de referencia"], horizontal=True, index=default_idx)
            expr_result = compute_fold_change_with_expression(
                controles_df,
                muestras_df,
                fold_source=fc_source,
                precomputed=fc,
            )
            df_expr = expr_result.expression_table
            cexp1, cexp2 = st.columns(2)
            with cexp1:
                st.dataframe(df_expr)
            with cexp2:
                dist_fig = build_expression_distribution(
                    df_expr,
                    exclude_stable=bool(app_state.exclude_stable),
                )
                st.plotly_chart(dist_fig, use_container_width=True)

            # Genes por nivel de expresión: lista/tablas y treemap
            st.markdown("### Genes por nivel de expresión")
            # Base para visualizaciones (respetar preferencia de exclusión de 'estables')
            df_expr_vis = df_expr[df_expr['nivel_expresion'] != 'estable'].copy() if bool(app_state.exclude_stable) else df_expr.copy()
            # Orden sugerido
            levels_vis = ['subexpresado', 'estable', 'sobreexpresado']
            levels_vis = [lvl for lvl in levels_vis if lvl in df_expr_vis['nivel_expresion'].astype(str).unique()]
            view_mode = st.radio("Vista", ["Lista por nivel", "Treemap"], horizontal=True, index=0)
            if view_mode == "Lista por nivel":
                tabs = st.tabs([lvl.capitalize() for lvl in levels_vis]) if levels_vis else []
                for lvl, tab in zip(levels_vis, tabs):
                    with tab:
                        sub = df_expr_vis[df_expr_vis['nivel_expresion'] == lvl].copy()
                        if sub.empty:
                            st.info("Sin genes para este nivel.")
                            continue
                        # Ordenar por |log2FC| descendente para resaltar extremos
                        try:
                            sub['log2fc'] = np.log2(sub['fold_change'].clip(lower=1e-12))
                            sub = sub.sort_values(sub['log2fc'].abs(), ascending=False)
                        except Exception:
                            pass
                        st.dataframe(sub[['target', 'fold_change', 'nivel_expresion']])
                        st.download_button(
                            label=f"Descargar genes ({lvl})",
                            data=sub[['target', 'fold_change', 'nivel_expresion']].to_csv(index=False),
                            file_name=f"genes_{lvl}.csv",
                            mime="text/csv",
                            use_container_width=True,
                        )
            else:
                treemap_fig = build_expression_treemap(df_expr_vis)
                if treemap_fig is None:
                    st.info("Sin genes para graficar con los filtros actuales.")
                else:
                    st.plotly_chart(treemap_fig, use_container_width=True)

            # Persistir para páginas
            st.session_state['fc_consolidated'] = fc.consolidated.copy()
            st.session_state['reference_gene'] = fc.reference_gene
            st.session_state['df_expr'] = df_expr.copy()
            extras['fold_change_consolidado.csv'] = fc.consolidated.to_csv(index=False)
            extras['expresion_categorizada.csv'] = df_expr.to_csv(index=False)

            # Anotación Ensembl (IDs y descripciones) sobre los targets clasificados
            st.subheader("Anotación Ensembl (IDs y descripciones)")
            st.caption("Consulta Ensembl para cada gen (requiere conexión a internet). Incluye exploración interactiva.")

            try:
                df_to_annot = df_expr[['target', 'nivel_expresion', 'fold_change']].drop_duplicates(subset=['target']).reset_index(drop=True)
                # Clave estable por lista de genes
                genes_key = ",".join(sorted(df_to_annot['target'].dropna().astype(str).unique().tolist()))
                ensembl_state_key = st.session_state.get('ensembl_key')
                ensembl_state_df = st.session_state.get('ensembl_df')

                cols = st.columns([1,1,2])
                with cols[0]:
                    need_compute = not (ensembl_state_key == genes_key and isinstance(ensembl_state_df, pd.DataFrame) and not ensembl_state_df.empty)
                    ens_workers = st.number_input("Hilos Ensembl", value=3, min_value=1, max_value=16, step=1, help="Límite de concurrencia para llamadas a Ensembl")
                    if need_compute:
                        if st.button("Anotar Ensembl", key="btn_annot_ens"):
                            with st.spinner("Consultando Ensembl…"):
                                ensembl_df = _annotate_ensembl_cached(df_to_annot[['target','nivel_expresion','fold_change']].sort_values('target').to_csv(index=False), max_workers=int(ens_workers))
                            st.session_state['ensembl_df'] = ensembl_df.copy()
                            st.session_state['ensembl_key'] = genes_key
                    else:
                        st.success("Usando anotación en caché para esta lista de genes.")
                        if st.button("Recalcular", key="btn_recalc_ens"):
                            with st.spinner("Consultando Ensembl…"):
                                ensembl_df = _annotate_ensembl_cached(df_to_annot[['target','nivel_expresion','fold_change']].sort_values('target').to_csv(index=False), max_workers=int(ens_workers))
                            st.session_state['ensembl_df'] = ensembl_df.copy()
                            st.session_state['ensembl_key'] = genes_key

                # Obtener df desde sesión si existe
                ensembl_df = st.session_state.get('ensembl_df')
                if isinstance(ensembl_df, pd.DataFrame) and not ensembl_df.empty:
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
                    if isinstance(ensembl_df, pd.DataFrame) and not ensembl_df.empty:
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
                    else:
                        st.info("Genera la anotación para ver enlaces.")

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
            # Respetar preferencia global de exclusión de 'estables'
            exclude_stable_pref = bool(app_state.exclude_stable)
            if exclude_stable_pref and 'estable' in sel_levels:
                sel_levels = [lvl for lvl in sel_levels if lvl != 'estable']
                st.caption("Preferencia activa: excluyendo 'estable' del enriquecimiento")

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
                species_default = int(getattr(SERVICES.string, "species", 9606) or 9606)
                species = st.number_input(
                    "Especie (NCBI Taxon)",
                    value=species_default,
                    min_value=1,
                    max_value=1_000_000,
                    step=1,
                )

            run_enrich = st.button("Ejecutar enriquecimiento (STRING)", type="primary")

            if run_enrich:
                try:
                    # Ejecutar enriquecimiento por nivel y concatenar (con caché)
                    levels_to_use = sel_levels or order_levels
                    if exclude_stable_pref:
                        levels_to_use = [lvl for lvl in levels_to_use if lvl != 'estable']
                    base = df_expr[df_expr['nivel_expresion'].isin(levels_to_use)].copy()
                    logger.info(f"STRING: niveles={levels_to_use}, fuentes={sel_cats}, genes_totales={base['target'].nunique()}")
                    with st.spinner("Consultando STRING por nivel…"):
                        res = _string_enrichment_cached(
                            base.to_csv(index=False),
                            levels=levels_to_use,
                            species=int(species),
                            sources=sel_cats,
                            max_fdr=float(max_fdr),
                            min_size=int(min_size),
                            top_n=int(top_n),
                        )
                    combined = res.get('combined')
                    by_level = res.get('by_level', {})
                    try:
                        logger.info(f"STRING: filas combinadas={0 if combined is None else len(combined)}")
                    except Exception:
                        pass

                    if combined is None or combined.empty:
                        st.info("Sin resultados de enriquecimiento para los parámetros actuales.")
                    else:
                        combined_f = combined
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
                                lvl_df_f = lvl_df
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
        pubmed_cfg = SERVICES.pubmed
        try:
            sec_email = st.secrets["NCBI_EMAIL"]
        except Exception:
            sec_email = ""
        try:
            sec_key = st.secrets["NCBI_API_KEY"]
        except Exception:
            sec_key = ""
        default_email = (sec_email or pubmed_cfg.default_email or "")
        default_api_key = (sec_key or (pubmed_cfg.api_key or ""))
        ccreds1, ccreds2 = st.columns([2, 2])
        with ccreds1:
            ncbi_email_input = st.text_input(
                "NCBI Email (obligatorio)",
                value=default_email,
                placeholder="tu_email@dominio.com",
            )
        with ccreds2:
            ncbi_api_key_input = st.text_input(
                "NCBI API Key (opcional)",
                value=default_api_key,
                placeholder="api_key",
            )

        max_default = int(pubmed_cfg.max_per_gene or 100)
        max_default = max(10, min(300, max_default))
        max_per_gene = st.number_input(
            "Máximo de artículos por gen",
            value=max_default,
            min_value=10,
            max_value=300,
            step=10,
        )
        run_pubmed = st.button("Buscar en PubMed", disabled=not bool(ncbi_email_input.strip()))
        st.caption("Las consultas a PubMed se cachean (24h) por combinación de genes, contexto y credenciales.")

        if run_pubmed:
            try:
                # Usar el DataFrame anotado si está disponible; si no, usar df_expr
                if 'ensembl_df' in locals() and isinstance(ensembl_df, pd.DataFrame) and not ensembl_df.empty:
                    genes_df = ensembl_df[['target', 'ensembl_id', 'nivel_expresion']].drop_duplicates('target')
                else:
                    tmp = df_expr[['target', 'nivel_expresion']].drop_duplicates('target')
                    tmp['ensembl_id'] = ''
                    genes_df = tmp[['target', 'ensembl_id', 'nivel_expresion']]
                # Respetar preferencia global de exclusión de 'estables'
                if bool(app_state.exclude_stable):
                    n_before = len(genes_df)
                    genes_df = genes_df[genes_df['nivel_expresion'] != 'estable']
                    n_after = len(genes_df)
                    if n_after < n_before:
                        st.caption(f"Preferencia activa: excluyendo 'estable' en PubMed ({n_before-n_after} genes menos)")

                # Usar credenciales explícitas (sin escribir en os.environ)
                email = ncbi_email_input.strip()
                api_key = ncbi_api_key_input.strip() or None
                logger.info(f"PubMed: consultando {len(genes_df)} genes (contexto={context_sel_label}, max_per_gene={int(max_per_gene)})")
                prog = st.progress(0)
                status = st.empty()

                def _on_progress(i: int, total: int, gene: str) -> None:
                    pct = int((i / max(1, total)) * 100)
                    prog.progress(pct)
                    status.info(f"Procesando {i}/{total}: {gene}")

                with st.spinner("Consultando PubMed por gen…"):
                    # Cache por combinación de parámetros
                    bib = _pubmed_cached(
                        genes_df.to_csv(index=False),
                        context_sel_label,
                        int(max_per_gene),
                        email,
                        api_key,
                    )
                prog.progress(100)
                status.success("Consulta PubMed finalizada")
                if bib is None or bib.empty:
                    st.info("No se encontraron artículos para los parámetros actuales.")
                else:
                    bib2 = merge_expression_levels(bib, genes_df[['target', 'nivel_expresion']])
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
                        # Enfocar al cáncer seleccionado en el panel
                        focused = filter_bibliography_by_cancer(classified, cancer_type)
                        try:
                            st.session_state["bibliografia_clasificada"] = focused.copy()
                        except Exception:
                            pass
                        st.dataframe(focused.head(50))
                        st.download_button(
                            label="Descargar bibliografía clasificada (CSV)",
                            data=focused.to_csv(index=False),
                            file_name="bibliografia_clasificada.csv",
                            mime="text/csv",
                            use_container_width=True,
                        )

                        # Gráfica de barras: número de estudios por tipo de cáncer y nivel
                        agg = aggregate_counts_by_level_and_cancer(focused)
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
with st.container():
    st.caption("Genera firmas por tipo de cáncer y nivel, con enriquecimiento de Hallmarks (MSigDB). Requiere gseapy y acceso a los GMT locales.")

    # Intentar recuperar bibliografía clasificada de esta sesión
    bib_class = st.session_state.get("bibliografia_clasificada")
    if 'classified' in locals() and isinstance(classified, pd.DataFrame) and not classified.empty:
        bib_class = classified.copy()
        st.session_state["bibliografia_clasificada"] = bib_class

    # Controles unificados (siempre visibles). El botón se desactiva si falta bibliografía.
    ready_bib = isinstance(bib_class, pd.DataFrame) and not bib_class.empty
    if not ready_bib:
        st.info("Ejecuta primero 'Bibliografía (PubMed)' y clasificación, o sube un CSV clasificado.")
    uploaded_bib = st.file_uploader("Opcional: subir bibliografía clasificada (CSV)", type=["csv"], key="upl_bib_class")
    if uploaded_bib is not None:
        try:
            bib_class = pd.read_csv(uploaded_bib)
            st.session_state["bibliografia_clasificada"] = bib_class
            ready_bib = True
            st.success("Bibliografía cargada correctamente.")
        except Exception as e:
            st.error(f"No se pudo leer el CSV: {e}")

    import os as _os
    def_hall = str(_PROJ_ROOT / "gen-sets_GSEA_MSigDB/gsea_hallmarks_formatted.gmt")
    def_back = str(_PROJ_ROOT / "gen-sets_GSEA_MSigDB/C5- ontology gene sets.gmt")
    hall_gmt = st.session_state.get("hall_gmt_path", def_hall)
    back_gmt = st.session_state.get("back_gmt_path", def_back)

    c1, c2, c3 = st.columns([2, 2, 2])
    with c1:
        hall_gmt = st.text_input("Ruta GMT Hallmarks", value=hall_gmt)
        st.session_state["hall_gmt_path"] = hall_gmt
        if not _os.path.exists(hall_gmt):
            upl = st.file_uploader("Subir GMT de Hallmarks", type=["gmt"], key="upl_hallmark_gmt")
            if upl is not None:
                try:
                    tmp_path = str(_PROJ_ROOT / "gen-sets_GSEA_MSigDB/_uploaded_hallmarks.gmt")
                    with open(tmp_path, "wb") as fh:
                        fh.write(upl.getbuffer())
                    hall_gmt = tmp_path
                    st.session_state["hall_gmt_path"] = hall_gmt
                    st.success("Hallmarks GMT cargado.")
                except Exception as e:
                    st.error(f"No se pudo guardar el GMT de Hallmarks: {e}")
    with c2:
        back_gmt = st.text_input("Ruta GMT background (opcional)", value=back_gmt)
        st.session_state["back_gmt_path"] = back_gmt
        if back_gmt and not _os.path.exists(back_gmt):
            upl_b = st.file_uploader("Subir GMT de background (opcional)", type=["gmt"], key="upl_back_gmt")
            if upl_b is not None:
                try:
                    tmp_path_b = str(_PROJ_ROOT / "gen-sets_GSEA_MSigDB/_uploaded_background.gmt")
                    with open(tmp_path_b, "wb") as fh:
                        fh.write(upl_b.getbuffer())
                    back_gmt = tmp_path_b
                    st.session_state["back_gmt_path"] = back_gmt
                    st.success("Background GMT cargado.")
                except Exception as e:
                    st.error(f"No se pudo guardar el GMT de background: {e}")
    with c3:
        ctx = st.selectbox("Contexto biológico", ["Cáncer y TEM", "Cáncer y micro RNAs"], index=0, key="sig_ctx")

    # Selección de tipo de cáncer (previa a generación)
    # Usar el tipo de cáncer definido en los parámetros del estudio
    sel_tipo_gen = cancer_type
    st.caption(f"Tipo de cáncer (parámetros del estudio): {sel_tipo_gen}")

    run_sig = st.button("Generar firmas", disabled=not ready_bib, key="btn_run_signatures")

    if run_sig:
        try:
            with st.spinner("Calculando firmas (incluye enriquecimiento de Hallmarks)…"):
                bib_csv = bib_class.to_csv(index=False)
                df_sigs = _compute_signatures_cached(bib_csv, ctx, hall_gmt, back_gmt, sel_tipo_gen or "")
            if df_sigs is None or df_sigs.empty:
                st.info("No se generaron firmas para los datos disponibles.")
            else:
                try:
                    st.session_state["df_signatures"] = df_sigs.copy()
                except Exception:
                    pass
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

                # Visualización se realiza debajo usando datos en sesión
        except Exception as e:
            st.error(f"No se pudieron generar las firmas: {e}")

        # Visualización: disponible siempre que existan firmas en sesión
        df_sigs_viz = st.session_state.get("df_signatures")
        if isinstance(df_sigs_viz, pd.DataFrame) and not df_sigs_viz.empty:
            st.subheader("Visualización de firmas")
            sel_tipo = sel_tipo_gen
            st.caption(f"Visualizando tipo de cáncer (parámetros del estudio): {sel_tipo}")
            if sel_tipo:
                try:
                    import plotly.express as px
                    recs = []
                    for _, row in df_sigs_viz[df_sigs_viz['cancer_type'] == sel_tipo].iterrows():
                        nivel = row.get('nivel_expresion')
                        genes_firma = row.get('genes', []) if isinstance(row.get('genes'), list) else []
                        counts_firma = row.get('conteo_articulos_por_gene', []) if isinstance(row.get('conteo_articulos_por_gene'), list) else []
                        gene_to_count = dict(zip(genes_firma, counts_firma))
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
                    if not recs:
                        st.info("No hay hallmarks asociados para este tipo de cáncer. Mostrando conteo de artículos por gen como alternativa.")
                        # Fallback: barras por nivel con conteo de artículos por gen
                        base = df_sigs_viz[df_sigs_viz['cancer_type'] == sel_tipo].copy()
                        # Expandir por nivel y pares (gen, conteo)
                        rows = []
                        for _, r in base.iterrows():
                            lvl = r.get('nivel_expresion')
                            genes = r.get('genes', []) if isinstance(r.get('genes'), list) else []
                            counts = r.get('conteo_articulos_por_gene', []) if isinstance(r.get('conteo_articulos_por_gene'), list) else []
                            for g, c in zip(genes, counts):
                                rows.append({'nivel_expresion': lvl, 'gene': g, 'articles': c})
                        if rows:
                            import plotly.express as px
                            flatc = pd.DataFrame(rows)
                            levels_sorted = sorted(flatc['nivel_expresion'].dropna().unique().tolist())
                            tabs = st.tabs([f"{lvl}" for lvl in levels_sorted])
                            for lvl, tab in zip(levels_sorted, tabs):
                                with tab:
                                    d = flatc[flatc['nivel_expresion'] == lvl]
                                    if d.empty:
                                        st.info("Sin datos para este nivel.")
                                        continue
                                    d = d.sort_values('articles', ascending=False).head(50)
                                    fig = px.bar(d, x='gene', y='articles', title=f"Artículos por gen — {sel_tipo} — {lvl}")
                                    fig.update_layout(height=600, margin=dict(t=60, b=100, l=20, r=20), xaxis_tickangle=45)
                                    st.plotly_chart(fig, use_container_width=True)
                        else:
                            st.info("No hay datos suficientes para este tipo de cáncer.")
                    else:
                        flat = pd.DataFrame(recs)
                        flat['log_p'] = -np.log10(flat['pvalue'].replace(0, 1e-300))
                        levels_sorted = sorted(flat['nivel_expresion'].dropna().unique().tolist())
                        tabs = st.tabs([f"{lvl}" for lvl in levels_sorted])
                        for lvl, tab in zip(levels_sorted, tabs):
                            with tab:
                                d = flat[flat['nivel_expresion'] == lvl]
                                if d.empty:
                                    st.info("Sin datos para este nivel.")
                                    continue
                                fig = px.sunburst(
                                    d,
                                    path=['gene', 'hallmark'],
                                    values='articles',
                                    color='log_p',
                                    color_continuous_scale='RdBu_r',
                                    title=f"Firmas - {sel_tipo} - {lvl}",
                                )
                                fig.update_layout(height=800, margin=dict(t=60, b=20, l=20, r=20))
                                st.plotly_chart(fig, use_container_width=True)
                except Exception as e:
                    st.warning(f"No se pudo renderizar la visualización de firmas: {e}")

            # Visualizaciones adicionales (polar, fingerprint, clustergrama)
            try:
                st.subheader("Visualizaciones adicionales")
                # Niveles disponibles en firmas para el tipo seleccionado
                niveles_disp = (
                    df_sigs_viz[df_sigs_viz['cancer_type'] == sel_tipo]['nivel_expresion']
                    .dropna().astype(str).unique().tolist()
                )
                niveles_disp = [n for n in ['Sobreexpresados', 'Subexpresados', 'estable'] if n in niveles_disp] or niveles_disp

                c1, c2 = st.columns(2)
                with c1:
                    with st.expander("Gráfico polar de Hallmarks", expanded=False):
                        if niveles_disp:
                            lvl = st.selectbox("Nivel", niveles_disp, key="viz_polar_lvl")
                            figp = hallmarks_polar_chart(df_sigs_viz, sel_tipo, lvl)
                            st.plotly_chart(figp, use_container_width=True)
                        else:
                            st.info("No hay niveles disponibles para el tipo seleccionado.")
                with c2:
                    with st.expander("Fingerprint (Genes × Hallmarks)", expanded=False):
                        if niveles_disp:
                            lvl2 = st.selectbox("Nivel", niveles_disp, key="viz_fp_lvl")
                            figf = fingerprint_heatmap(df_sigs_viz, sel_tipo, lvl2)
                            st.plotly_chart(figf, use_container_width=True)
                        else:
                            st.info("No hay niveles disponibles para el tipo seleccionado.")

                with st.expander("Clustergrama de expresión relativa", expanded=False):
                    try:
                        # Recuperar clasificación y datos de expresión si existen
                        if df_loaded is None:
                            st.info("Primero carga y clasifica un Excel de qPCR para ver el clustergrama.")
                        else:
                            file_key = f"assign_{df_loaded.source_name}:{df_loaded.sheet_name}"
                            state_local = st.session_state.get(file_key, {})
                            controles_df = state_local.get('controles_df', pd.DataFrame())
                            muestras_df = state_local.get('muestras_df', pd.DataFrame())
                            if (controles_df is None or controles_df.empty) or (muestras_df is None or muestras_df.empty):
                                st.info("Clasifica controles y muestras en la sección inicial para habilitar el clustergrama.")
                            else:
                                # Listas de genes desde las firmas para el tipo y niveles clásicos
                                lista_sobre, lista_sub = [], []
                                base = df_sigs_viz[df_sigs_viz['cancer_type'] == sel_tipo]
                                try:
                                    lista_sobre = base.loc[base['nivel_expresion'] == 'Sobreexpresados', 'genes'].iloc[0]
                                except Exception:
                                    lista_sobre = []
                                try:
                                    lista_sub = base.loc[base['nivel_expresion'] == 'Subexpresados', 'genes'].iloc[0]
                                except Exception:
                                    lista_sub = []
                                lista_sobre = lista_sobre if isinstance(lista_sobre, list) else ([] if pd.isna(lista_sobre) else [str(lista_sobre)])
                                lista_sub = lista_sub if isinstance(lista_sub, list) else ([] if pd.isna(lista_sub) else [str(lista_sub)])

                                metodo = st.radio("Método de normalización", ["gen de referencia", "promedios"], index=0, horizontal=True, key="viz_cluster_method")
                                ref_default = st.session_state.get('reference_gene', '')
                                gen_ref = st.text_input("Gen de referencia (si aplica)", value=str(ref_default) if ref_default else "", key="viz_cluster_ref")
                                figs = {}
                                try:
                                    figs = clustered_expression_by_level(controles_df, muestras_df, lista_sobre, lista_sub, metodo_elegido=metodo, gen_referencia=(gen_ref or None))
                                except Exception as ex:
                                    st.warning(f"No se pudo generar el clustergrama: {ex}")
                                    figs = {}
                                # Mostrar
                                if not figs:
                                    st.info("Sin datos suficientes para construir el clustergrama con los parámetros actuales.")
                                else:
                                    for k in ("Sobreexpresados", "Subexpresados"):
                                        if k in figs:
                                            st.plotly_chart(figs[k], use_container_width=True)
                    except Exception:
                        st.info("Clustergrama no disponible (faltan dependencias opcionales o datos).")
            except Exception:
                pass


# -------------------------------------------------------------------------
# Heurística (PubMed) y relación con niveles/firmas
# -------------------------------------------------------------------------
st.markdown("---")
st.header("Heurística bibliográfica y relación con niveles/firmas")

focused = st.session_state.get("bibliografia_clasificada")
if not (isinstance(focused, pd.DataFrame) and not focused.empty):
    st.info("Primero ejecuta la búsqueda en PubMed y la clasificación para generar la heurística.")
else:
    try:
        st.subheader("Interpretación heurística de relaciones (enfocada al cáncer seleccionado)")
        interp = interpret_gene_relations(focused)
        if interp is None or interp.empty:
            st.info("No se pudieron etiquetar relaciones con las heurísticas actuales.")
        else:
            summary = summarize_relations_by_gene(interp)
            st.dataframe(summary.sort_values('heuristic_summary').reset_index(drop=True))
            st.download_button(
                label="Descargar resumen heurístico (CSV)",
                data=summary.to_csv(index=False),
                file_name="resumen_heuristico_relaciones.csv",
                mime="text/csv",
                use_container_width=True,
            )

            # Heatmap genes × funciones (scores normalizados)
            try:
                import numpy as np
                import plotly.express as px
                func_cols = [
                    c for c in summary.columns
                    if c.endswith('_score') and not c.startswith(('upregulated','downregulated','prognosis_'))
                ]
                if func_cols:
                    sm = summary.copy()
                    sm['_total_score'] = sm[func_cols].sum(axis=1)
                    top = sm.sort_values('_total_score', ascending=False).head(30)
                    Z = (top[func_cols] - top[func_cols].min()) / (top[func_cols].max() - top[func_cols].min() + 1e-9)
                    fig_hm = px.imshow(
                        Z.values,
                        labels=dict(x='Función', y='Gen', color='Score norm.'),
                        x=[c.replace('_score','') for c in func_cols],
                        y=top['Gene'].tolist(),
                        aspect='auto',
                        title='Heatmap genes × funciones (scores normalizados)'
                    )
                    fig_hm.update_layout(margin=dict(l=60,r=20,t=40,b=40))
                    st.plotly_chart(fig_hm, use_container_width=True)
            except Exception:
                pass

            # Sankey: Gen → efecto expresión → pronóstico dominante
            try:
                import numpy as np
                import plotly.graph_objects as go
                bad = summary.get('prognosis_bad_score', 0)
                good = summary.get('prognosis_good_score', 0)
                if hasattr(bad, 'values') and hasattr(good, 'values'):
                    dom = np.where(bad.values > good.values, 'prognosis_bad', np.where(good.values > 0, 'prognosis_good', 'prognosis_uncertain'))
                    left = summary['Gene'].tolist()
                    midl = summary['net_expression_effect'].astype(str).tolist() if 'net_expression_effect' in summary.columns else ['uncertain'] * len(left)
                    right = dom.tolist()
                    nodes = list(dict.fromkeys(left + midl + right))
                    idx = {n:i for i,n in enumerate(nodes)}
                    links = []
                    links += [(idx[l], idx[m], 1) for l, m in zip(left, midl)]
                    links += [(idx[m], idx[r], 1) for m, r in zip(midl, right)]
                    link = dict(source=[s for s,_,_ in links], target=[t for _,t,_ in links], value=[v for _,_,v in links])
                    node = dict(label=nodes, pad=12, thickness=12)
                    fig_sk = go.Figure(go.Sankey(node=node, link=link))
                    fig_sk.update_layout(margin=dict(l=10,r=10,t=10,b=10), title='Flujo expresión → pronóstico')
                    st.plotly_chart(fig_sk, use_container_width=True)
            except Exception:
                pass

            # Relación heurística ↔ niveles de expresión y firmas
            try:
                import numpy as np
                import plotly.express as px
                df_levels = st.session_state.get('df_expr')
                if isinstance(df_levels, pd.DataFrame) and not df_levels.empty:
                    lev = df_levels[['target', 'nivel_expresion', 'fold_change']].drop_duplicates('target').rename(columns={'target':'Gene'})
                    merged = pd.merge(summary, lev, on='Gene', how='left')
                    if bool(app_state.exclude_stable):
                        merged = merged[merged['nivel_expresion'] != 'estable']
                    func_score_cols = [c for c in merged.columns if c.endswith('_score') and not c.startswith(('upregulated','downregulated','prognosis_'))]
                    func_flags = [c.replace('_score','') for c in func_score_cols]
                    long_rows = []
                    for fscore, fname in zip(func_score_cols, func_flags):
                        sub = merged[["Gene", "nivel_expresion", fscore]].copy()
                        sub['funcion'] = fname
                        sub['flag'] = (sub[fscore] >= 1.2).astype(int)
                        long_rows.append(sub.rename(columns={fscore: 'score'}))
                    long = pd.concat(long_rows, ignore_index=True) if long_rows else pd.DataFrame()
                    if not long.empty:
                        st.subheader("Niveles × funciones")
                        colm1, colm2, colm3 = st.columns([1,1,1])
                        with colm1:
                            agg_mode = st.selectbox("Métrica", ["conteo", "score"], index=1, key="agg_mode_heu")
                        with colm2:
                            top_funcs = st.number_input("Top funciones", min_value=3, max_value=20, value=8, step=1, key="top_funcs_heu")
                        with colm3:
                            norm_scores = st.checkbox("Normalizar por columna", value=True, key="norm_scores_heu")
                        if agg_mode == 'conteo':
                            pivot = long.groupby(['nivel_expresion','funcion'])['flag'].sum().unstack(fill_value=0)
                        else:
                            pivot = long.groupby(['nivel_expresion','funcion'])['score'].sum().unstack(fill_value=0)
                        totals = pivot.sum(axis=0).sort_values(ascending=False)
                        top_cols = totals.head(int(top_funcs)).index.tolist()
                        pivot_top = pivot[top_cols]
                        if norm_scores:
                            pivot_top = (pivot_top - pivot_top.min()) / (pivot_top.max() - pivot_top.min() + 1e-9)
                        fig_lv = px.imshow(
                            pivot_top.values,
                            labels=dict(x='Función', y='Nivel expresión', color='Valor'),
                            x=top_cols,
                            y=pivot_top.index.tolist(),
                            aspect='auto',
                            title=f"Niveles × Funciones ({agg_mode})"
                        )
                        st.plotly_chart(fig_lv, use_container_width=True)

                        # Treemap: Nivel → Función → Gen (score×|log2FC|)
                        trea = long.copy().merge(lev, on='Gene', how='left')
                        trea['log2fc_abs'] = np.log2(trea['fold_change'].clip(lower=1e-12)).abs()
                        trea['weight'] = trea['score'] * (1.0 + trea['log2fc_abs'])
                        fig_tree = px.treemap(
                            trea,
                            path=['nivel_expresion','funcion','Gene'],
                            values='weight',
                            color='score', color_continuous_scale='Viridis',
                            title='Treemap: Nivel → Función → Gen (peso = score×|log2FC|)'
                        )
                        st.plotly_chart(fig_tree, use_container_width=True)

                        # Función → Hallmarks si hay firmas
                        df_sigs_v = st.session_state.get('df_signatures')
                        if isinstance(df_sigs_v, pd.DataFrame) and not df_sigs_v.empty:
                            st.subheader("Funciones ↔ Hallmarks (firmas)")
                            try:
                                base_s = df_sigs_v[df_sigs_v.get('cancer_type') == cancer_type]
                                recs = []
                                for _, r in base_s.iterrows():
                                    for col in r.index:
                                        if col.startswith('hallmark_') and col.endswith('_genes'):
                                            term = col[len('hallmark_'):-len('_genes')]
                                            genes_hm = r[col] if isinstance(r[col], list) else []
                                            for g in genes_hm:
                                                recs.append({'Gene': g, 'hallmark': term})
                                if recs:
                                    map_hm = pd.DataFrame(recs).drop_duplicates()
                                    mf = pd.merge(long[['Gene','funcion','flag']], map_hm, on='Gene', how='inner')
                                    agg_hm = mf.groupby(['funcion','hallmark'])['flag'].sum().reset_index(name='n')
                                    import plotly.graph_objects as go
                                    funcs = agg_hm['funcion'].unique().tolist()
                                    terms = agg_hm['hallmark'].unique().tolist()
                                    nodes = funcs + terms
                                    idx = {n:i for i,n in enumerate(nodes)}
                                    src = [idx[f] for f in agg_hm['funcion']]
                                    tgt = [idx[t] for t in agg_hm['hallmark']]
                                    val = agg_hm['n'].astype(float).tolist()
                                    fig_fh = go.Figure(go.Sankey(node=dict(label=nodes, pad=12, thickness=12), link=dict(source=src, target=tgt, value=val)))
                                    fig_fh.update_layout(margin=dict(l=10,r=10,t=10,b=10), title='Funciones → Hallmarks (conteos)')
                                    st.plotly_chart(fig_fh, use_container_width=True)
                            except Exception:
                                pass
            except Exception:
                pass
    except Exception as e:
        st.warning(f"No se pudo ejecutar la heurística: {e}")

# -----------------------------------------------------------------------------
# Insights (Google NLP) sobre bibliografía filtrada por el cáncer seleccionado
# -----------------------------------------------------------------------------
st.markdown("---")
st.header("Insights de la literatura (Google NLP)")

google_cfg = SERVICES.google_nlp
try:
    sec_gkey = st.secrets["GOOGLE_NLP_API_KEY"]
except Exception:
    sec_gkey = ""
default_google_key = sec_gkey or (google_cfg.api_key or "")

ins_col1, ins_col2, ins_col3 = st.columns([2,1,1])
with ins_col1:
    g_api_key = st.text_input(
        "Google NLP API Key",
        value=default_google_key,
        type="password",
        help="Se recomienda usar st.secrets['GOOGLE_NLP_API_KEY'] o la variable de entorno GOOGLE_NLP_API_KEY.",
    )
with ins_col2:
    n_docs = st.number_input("Máx. artículos a analizar", min_value=5, max_value=200, value=30, step=5)
with ins_col3:
    lang_sel = st.selectbox("Idioma (opcional)", options=["auto", "es", "en"], index=0, help="Déjalo en 'auto' para detección automática.")

apply_cancer_filter = st.checkbox(
    "Filtrar por el tipo de cáncer seleccionado",
    value=True,
    help="Desactívalo si tu CSV incluye artículos de otros cánceres y quieres analizarlos igualmente.",
)
run_insights = st.button("Generar insights (Google NLP)")

uploaded_bib_ins = st.file_uploader(
    "Subir bibliografía (CSV) — columnas esperadas: Title, Abstract, Gene (opcional)",
    type=["csv"],
    key="upl_bib_nlp_csv",
    help="Si subes un archivo aquí, se usará para los insights en lugar de la bibliografía de la sesión.",
)

# Selección de gen (desde dataset o manual si no existe columna)
gene_selected: Optional[str] = None
gene_source_df = None
def _read_uploaded_csv(uploaded_file):
    if uploaded_file is None:
        return None
    # Intenta distintas estrategias para archivos regionales (coma/punto y coma) y BOM
    for attempt in range(3):
        try:
            uploaded_file.seek(0)
        except Exception:
            pass
        try:
            if attempt == 0:
                return pd.read_csv(uploaded_file)
            elif attempt == 1:
                return pd.read_csv(uploaded_file, sep=None, engine='python')
            else:
                return pd.read_csv(uploaded_file, sep=';', encoding='utf-8-sig')
        except Exception:
            continue
    return None

try:
    if uploaded_bib_ins is not None:
        gene_source_df = _read_uploaded_csv(uploaded_bib_ins)
    else:
        gene_source_df = st.session_state.get("bibliografia_clasificada")
except Exception:
    gene_source_df = None

if isinstance(gene_source_df, pd.DataFrame) and not gene_source_df.empty:
    # Buscar columnas de gen razonables
    gene_col = None
    for c in ["Gene", "gene", "Symbol", "symbol", "target"]:
        if c in gene_source_df.columns:
            gene_col = c
            break
    if gene_col:
        uniq_genes = (
            gene_source_df[gene_col].dropna().astype(str).str.strip().unique().tolist()
        )
        uniq_genes = sorted([g for g in uniq_genes if g])
        if uniq_genes:
            gene_selected = st.selectbox("Gen a analizar (desde archivo/datos)", options=uniq_genes, index=0)
    if not gene_selected:
        gene_selected = st.text_input("Gen a analizar (manual si no hay columna)", value="", placeholder="Ej.: TP53")
else:
    gene_selected = st.text_input("Gen a analizar", value="", placeholder="Ej.: TP53")

@st.cache_data(show_spinner=False, ttl=3600)
def _prep_nlp_texts_cached(df_csv: str, cancer_label: str, gene: Optional[str], limit: int, apply_filter: bool) -> list[str]:
    try:
        df = pd.read_csv(io.StringIO(df_csv))
    except Exception:
        return []
    t = df.copy()
    if apply_filter:
        try:
            from src.core.bibliography import filter_bibliography_by_cancer
            t = filter_bibliography_by_cancer(t, cancer_label)
        except Exception:
            pass
    # Filtrar por gen si aplica
    if gene:
        gcand = None
        for c in ["Gene", "gene", "Symbol", "symbol", "target"]:
            if c in t.columns:
                gcand = c
                break
        if gcand:
            t = t[t[gcand].astype(str).str.strip().str.lower() == str(gene).strip().lower()]
    # Tomar columnas de texto disponibles
    title_col = "Title" if "Title" in t.columns else ("article_title" if "article_title" in t.columns else None)
    abs_col = "Abstract" if "Abstract" in t.columns else None
    if not title_col:
        return []
    texts = (t[title_col].astype(str) + " " + (t[abs_col].astype(str) if abs_col else "")).tolist()
    texts = [s.strip() for s in texts if isinstance(s, str) and s.strip()]
    return texts[: max(1, int(limit))]

if run_insights:
    # Priorizar archivo subido; si no, usar bibliografía de sesión
    if uploaded_bib_ins is not None:
        # Reutiliza el DataFrame ya leído (evita puntero al final del buffer)
        if isinstance(gene_source_df, pd.DataFrame) and not gene_source_df.empty:
            bib_df = gene_source_df.copy()
        else:
            try:
                bib_df = _read_uploaded_csv(uploaded_bib_ins)
            except Exception as e:
                st.error(f"No se pudo leer el CSV subido: {e}")
                bib_df = None
        if bib_df is None:
            st.error("No se pudo leer el CSV subido (prueba con separador ';' o codificación UTF-8).")
    else:
        bib_df = st.session_state.get("bibliografia_clasificada")
    if not (isinstance(bib_df, pd.DataFrame) and not bib_df.empty):
        st.info("Primero genera o carga la bibliografía clasificada en la sección anterior.")
    elif not (g_api_key or google_cfg.api_key):
        st.warning("Falta la API Key de Google NLP. Ingresa una o configúrala en el entorno/secrets.")
    else:
        with st.spinner("Analizando textos con Google NLP…"):
            # Conteos informativos antes del análisis
            try:
                df_all = bib_df.copy()
                # Arma columna de texto combinada
                title_col = "Title" if "Title" in df_all.columns else ("article_title" if "article_title" in df_all.columns else None)
                abs_col = "Abstract" if "Abstract" in df_all.columns else None
                if title_col is None:
                    raise RuntimeError("El CSV debe incluir una columna 'Title' o 'article_title'.")
                df_all['_text'] = df_all[title_col].astype(str) + ' ' + (df_all[abs_col].astype(str) if abs_col else '')
                total_docs_all = int((df_all['_text'].astype(str).str.strip() != '').sum())
                # Filtrar por gen
                df_gene = df_all.copy()
                if gene_selected:
                    gcand = None
                    for c in ["Gene", "gene", "Symbol", "symbol", "target"]:
                        if c in df_gene.columns:
                            gcand = c
                            break
                    if gcand:
                        df_gene = df_gene[df_gene[gcand].astype(str).str.strip().str.lower() == str(gene_selected).strip().lower()]
                total_docs_gene = int((df_gene['_text'].astype(str).str.strip() != '').sum())
                # Filtrar por cáncer (si aplica)
                df_cancer_gene = df_gene.copy()
                if apply_cancer_filter:
                    try:
                        from src.core.bibliography import filter_bibliography_by_cancer
                        df_cancer_gene = filter_bibliography_by_cancer(df_cancer_gene, cancer_type)
                    except Exception:
                        pass
                total_docs_final = int(len(df_cancer_gene))
                st.caption(
                    f"Artículos: {total_docs_all} en total → {total_docs_gene} para el gen seleccionado"
                    + (f" → {total_docs_final} tras filtrar por '{cancer_type}'" if apply_cancer_filter else "")
                )
            except Exception:
                total_docs_final = None

            texts = _prep_nlp_texts_cached(
                bib_df.to_csv(index=False), cancer_type, gene_selected or None, int(n_docs), bool(apply_cancer_filter)
            )
            if not texts:
                if apply_cancer_filter:
                    st.info("No hay textos para analizar con el filtro de cáncer activo. Desactívalo para analizar todo tu CSV.")
                else:
                    st.info("No hay textos válidos en el CSV (requiere columnas Title/Abstract con contenido).")
            else:
                try:
                    glang = None if lang_sel == "auto" else lang_sel
                    client = GoogleNLPClient(api_key=(g_api_key or None), default_language=glang, sleep_between=0.0)
                    res = aggregate_insights(
                        texts,
                        client,
                        do_entities=True,
                        do_entity_sentiment=True,
                        do_sentiment=True,
                        do_categories=True,
                        language=glang,
                        max_chars_per_doc=8000,
                    )
                except Exception as e:
                    st.error(f"Error llamando Google NLP: {e}")
                    res = None

                if res:
                    # Breve párrafo resumen para el gen seleccionado (conclusiones)
                    try:
                        gname = (gene_selected or '').strip()
                        sent = res.get("sentiment", {}) or {}
                        cats = res.get("categories", []) or []
                        ents = res.get("entities", []) or []
                        entsent = res.get("entity_sentiment", []) or []
                        n_docs_used = int(sent.get('n', 0))
                        avg_s = float(sent.get('avg_score', 0.0))
                        avg_m = float(sent.get('avg_magnitude', 0.0))
                        # Top categorías
                        top_cats = [c.get('category','') for c in cats[:3]] if cats else []
                        # Entidades salientes
                        top_ent_names = [e.get('name','') for e in ents[:3]] if ents else []
                        # Entidades con polaridad marcada
                        pos_ents = sorted([e for e in entsent if float(e.get('avg_sentiment',0)) > 0.2], key=lambda x: (x.get('avg_sentiment',0), x.get('avg_magnitude',0)), reverse=True)[:3]
                        neg_ents = sorted([e for e in entsent if float(e.get('avg_sentiment',0)) < -0.2], key=lambda x: (abs(x.get('avg_sentiment',0)), x.get('avg_magnitude',0)), reverse=True)[:3]
                        pos_list = [p.get('name','') for p in pos_ents]
                        neg_list = [n.get('name','') for n in neg_ents]
                        # Construir narrativa breve orientada a conclusiones
                        partes = []
                        partes.append(
                            f"Se revisaron {n_docs_used} artículos" + (f" sobre {gname}" if gname else "") + (f" en {cancer_type}" if apply_cancer_filter else "") + "."
                        )
                        if top_cats:
                            partes.append(f"Los temas dominantes fueron: {', '.join(top_cats)}.")
                        if top_ent_names:
                            partes.append(f"Las entidades más relevantes incluyeron {', '.join(top_ent_names)}.")
                        if pos_list or neg_list:
                            if pos_list:
                                partes.append(f"Se observaron asociaciones positivas destacadas con {', '.join(pos_list)}.")
                            if neg_list:
                                partes.append(f"También emergieron señales negativas en torno a {', '.join(neg_list)}.")
                        partes.append(f"El tono global fue {('ligeramente negativo' if avg_s < -0.1 else 'ligeramente positivo' if avg_s > 0.1 else 'neutral')} (score {avg_s:.2f}, magnitud {avg_m:.2f}).")
                        resumen = " ".join(partes)
                        st.info(resumen)
                    except Exception:
                        pass

                    # Complemento con heurística propia si hay bibliografía disponible (refuerza conclusiones)
                    try:
                        df_sub = bib_df.copy()
                        # Filtrar por gen
                        gcand = None
                        for c in ["Gene", "gene", "Symbol", "symbol", "target"]:
                            if c in df_sub.columns:
                                gcand = c
                                break
                        if gcand and gname:
                            df_sub = df_sub[df_sub[gcand].astype(str).str.strip().str.lower() == gname.lower()]
                        # Filtrar por cáncer si aplica
                        if apply_cancer_filter:
                            try:
                                from src.core.bibliography import filter_bibliography_by_cancer
                                df_sub = filter_bibliography_by_cancer(df_sub, cancer_type)
                            except Exception:
                                pass
                        # Aplicar interpretación heurística
                        ih = interpret_gene_relations(df_sub)
                        summ = summarize_relations_by_gene(ih)
                        if isinstance(summ, pd.DataFrame) and not summ.empty:
                            # Toma fila del gen
                            row = summ.iloc[0]
                            effects = []
                            if float(row.get('upregulated_score', 0)) > float(row.get('downregulated_score', 0)):
                                effects.append('tendencia a sobreexpresión')
                            elif float(row.get('downregulated_score', 0)) > float(row.get('upregulated_score', 0)):
                                effects.append('tendencia a subexpresión')
                            prog = 'desfavorable' if float(row.get('prognosis_bad_score', 0)) > float(row.get('prognosis_good_score', 0)) else ('favorable' if float(row.get('prognosis_good_score', 0)) > 0 else 'incierta')
                            flags = [k for k in ['proliferation','invasion','migration','metastasis','drug_resistance','drug_sensitivity','apoptosis','emt_related'] if row.get(k, 0) > 0]
                            texto = (
                                f"Heurística complementaria: para {gname or 'el gen seleccionado'} se observa {', '.join(effects) if effects else 'un efecto de expresión incierto'}, "
                                f"y una señal pronóstica {prog}. Funciones implicadas: {', '.join(flags) if flags else 'no concluyentes'}."
                            )
                            st.caption(texto)
                    except Exception:
                        pass

                    # Resumen de sentimiento global
                    sent = res.get("sentiment", {}) or {}
                    s1, s2, s3 = st.columns(3)
                    with s1:
                        st.metric("Sentimiento promedio", f"{float(sent.get('avg_score', 0.0)):.2f}")
                    with s2:
                        st.metric("Magnitud promedio", f"{float(sent.get('avg_magnitude', 0.0)):.2f}")
                    with s3:
                        st.metric("Artículos analizados", int(sent.get('n', 0)))

                    # Entidades más salientes
                    ents = pd.DataFrame(res.get("entities") or [])
                    if not ents.empty:
                        st.subheader("Entidades más relevantes")
                        show = ents.rename(columns={"salience_sum": "saliencia_acum", "mentions": "menciones"})
                        st.dataframe(show.head(50))

                    # Entidades con mayor carga emocional
                    entsent = pd.DataFrame(res.get("entity_sentiment") or [])
                    if not entsent.empty:
                        st.subheader("Entidades con sentimiento (promedios)")
                        show2 = entsent.rename(columns={
                            "avg_sentiment": "sent_prom",
                            "avg_magnitude": "magn_prom",
                            "mentions": "menciones",
                        })
                        st.dataframe(show2.head(50))

                    # Categorías temáticas
                    cats = pd.DataFrame(res.get("categories") or [])
                    if not cats.empty:
                        st.subheader("Categorías (classifyText)")
                        try:
                            import plotly.express as px
                            figc = px.bar(cats.head(20), x="confidence_sum", y="category", orientation="h", title="Top categorías")
                            st.plotly_chart(figc, use_container_width=True)
                        except Exception:
                            pass
                        st.dataframe(cats)

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
