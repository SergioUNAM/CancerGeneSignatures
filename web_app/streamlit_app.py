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

import io, os, json, sys, re
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
    classify_by_prefixes,
    classify_by_suffixes,
)
# Compatibilidad: si la instalación local aún no expone classify_by_regex, usar fallback local
try:
    from src.core.qpcr import classify_by_regex  # type: ignore
except Exception:
    def classify_by_regex(df_long, ctrl_pattern: str, sample_pattern: str):
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
        return (ctrl.drop(columns=["test_str"]) if not ctrl.empty else ctrl,
                samp.drop(columns=["test_str"]) if not samp.empty else samp)
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
def _string_enrichment_cached(df_expr_csv: str, levels: list[str], species: int, sources: list[str], max_fdr: float, min_size: int, top_n: int) -> dict:
    try:
        df_expr = pd.read_csv(io.StringIO(df_expr_csv))
    except Exception:
        return {"combined": pd.DataFrame(), "by_level": {}}
    from src.core.string_enrichment import enrich_by_levels, filter_enrichment
    res = enrich_by_levels(df_expr, symbol_col='target', level_col='nivel_expresion', levels=levels, species=species, sources=sources)
    combined = res.get('combined')
    by_level = res.get('by_level') or {}
    # Aplicar filtros aquí para cachear resultado final
    out_by_level = {}
    if isinstance(by_level, dict):
        for lvl, d in by_level.items():
            out_by_level[lvl] = filter_enrichment(d, include_categories=sources, max_fdr=max_fdr, min_term_genes=min_size, top_n=top_n)
    out_combined = filter_enrichment(combined, include_categories=sources, max_fdr=max_fdr, min_term_genes=min_size, top_n=top_n)
    return {"combined": out_combined, "by_level": out_by_level}

@st.cache_data(show_spinner=False, ttl=86400)
def _pubmed_cached(genes_df_csv: str, context: str, max_per_gene: int, email: str, api_key: str|None) -> pd.DataFrame:
    try:
        genes_df = pd.read_csv(io.StringIO(genes_df_csv))
    except Exception:
        return pd.DataFrame()
    return search_pubmed_by_genes(
        genes_df.rename(columns={'target': 'target'}),
        symbol_col='target',
        ensembl_col='ensembl_id',
        selected_context=context,
        max_per_gene=int(max_per_gene),
        progress=None,
        logger=None,
        email=email,
        api_key=api_key,
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

    st.header("3) Política 'Undetermined/ND'")
    und_policy = st.selectbox(
        "Cómo tratar valores 'Undetermined' (ND)",
        options=["nan", "ctmax", "value"],
        index=0,
        help="\n- nan: deja como NaN (se puede imputar después).\n- ctmax: usa el Ct máximo observado por columna.\n- value: usa un valor fijo (p. ej., 40).",
    )
    und_value = st.number_input(
        "Valor fijo para 'value'",
        value=40.0,
        min_value=0.0,
        max_value=100.0,
        step=0.5,
    )
    # Persistir en sesión para usarlo durante todo el flujo
    st.session_state["und_policy"] = und_policy
    st.session_state["und_value"] = float(und_value)
    st.session_state["context_sel_label"] = context_sel_label

    # Preferencia global de gráficos: excluir genes 'estables'
    exclude_stable = st.checkbox(
        "Excluir genes 'estables' en gráficos",
        value=bool(st.session_state.get("exclude_stable", False)),
        help="Cuando está activado, los gráficos omiten los genes clasificados como 'estables'.",
    )
    st.session_state["exclude_stable"] = bool(exclude_stable)
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
                undetermined_policy=und_policy,
                undetermined_value=float(und_value),
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
                df_loaded = parse_qpcr_wide(
                    uploaded,
                    sheet_name=sheet,
                    header_mode="auto",
                    undetermined_policy=und_policy,
                    undetermined_value=float(und_value),
                )
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
        # Usar todos los nombres de prueba (sin omitir el primero) y filtrar vacíos
        tests_filtered = [s for s in sample_names if str(s).strip()]
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
        inter = set(ctrl_df['test'].astype(str)).intersection(set(samp_df['test'].astype(str)))
        if not inter:
            return ctrl_df, samp_df
        st.warning(f"Colisiones: {len(inter)} pruebas aparecen en ambos grupos → {', '.join(sorted(list(inter))[:10])}{' …' if len(inter)>10 else ''}")
        choice = st.radio("Resolver colisiones", ["priorizar controles", "priorizar muestras", "excluir colisiones"], horizontal=True, index=0)
        if choice == "priorizar controles":
            samp_df = samp_df[~samp_df['test'].astype(str).isin(inter)]
        elif choice == "priorizar muestras":
            ctrl_df = ctrl_df[~ctrl_df['test'].astype(str).isin(inter)]
        else:
            ctrl_df = ctrl_df[~ctrl_df['test'].astype(str).isin(inter)]
            samp_df = samp_df[~samp_df['test'].astype(str).isin(inter)]
        return ctrl_df, samp_df

    if submitted_pref:
        if not ctrl_prefix and not samp_prefix:
            st.warning("Debes ingresar al menos un prefijo (controles o muestras)")
        pref_ctrl_df, pref_samp_df = classify_by_prefixes(long_df, [ctrl_prefix] if ctrl_prefix else [], [samp_prefix] if samp_prefix else [])
        pref_ctrl_df, pref_samp_df = _resolve_collisions(pref_ctrl_df, pref_samp_df)
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
        pref_ctrl_df, pref_samp_df = classify_by_suffixes(long_df, suff_ctrl_list, suff_samp_list)
        pref_ctrl_df, pref_samp_df = _resolve_collisions(pref_ctrl_df, pref_samp_df)
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
            pref_ctrl_df, pref_samp_df = classify_by_regex(long_df, ctrl_re, samp_re)
            pref_ctrl_df, pref_samp_df = _resolve_collisions(pref_ctrl_df, pref_samp_df)
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

        # Panel de calidad antes de FC
        st.subheader("Calidad de datos (pre-FC)")
        ok_for_fc = True
        # Métricas básicas
        ctrl_targets = set(controles_df['target'].dropna().astype(str)) if not controles_df.empty else set()
        samp_targets = set(muestras_df['target'].dropna().astype(str)) if not muestras_df.empty else set()
        common_targets = ctrl_targets.intersection(samp_targets)
        ctrl_nan_ratio = float(controles_df['ct'].isna().mean()) if not controles_df.empty else 1.0
        samp_nan_ratio = float(muestras_df['ct'].isna().mean()) if not muestras_df.empty else 1.0

        mqa1, mqa2, mqa3, mqa4 = st.columns(4)
        with mqa1:
            st.metric("Genes (controles)", len(ctrl_targets))
        with mqa2:
            st.metric("Genes (muestras)", len(samp_targets))
        with mqa3:
            st.metric("Genes en común", len(common_targets))
        with mqa4:
            st.metric("NaN ct (ctrl/mues)", f"{ctrl_nan_ratio:.0%} / {samp_nan_ratio:.0%}")

        if len(common_targets) == 0:
            st.error("No hay intersección de genes entre controles y muestras. Ajusta tu clasificación.")
            ok_for_fc = False
        if ctrl_nan_ratio >= 1.0 or samp_nan_ratio >= 1.0:
            st.error("Todos los valores Ct son NaN en algún grupo. Revisa la política de 'Undetermined' o tus datos.")
            ok_for_fc = False

        # N mínimo por gen y grupo
        n_min = st.number_input("Mínimo de réplicas por gen y grupo", min_value=1, max_value=10, value=1, step=1)
        if common_targets:
            counts_ctrl = controles_df.dropna(subset=['ct']).groupby('target')['ct'].size()
            counts_samp = muestras_df.dropna(subset=['ct']).groupby('target')['ct'].size()
            eligible = [t for t in common_targets if counts_ctrl.get(t, 0) >= n_min and counts_samp.get(t, 0) >= n_min]
            st.info(f"Genes que cumplen n≥{n_min} en ambos grupos: {len(eligible)}")
            if len(eligible) == 0:
                st.warning("Ningún gen cumple el mínimo de réplicas en ambos grupos. Considera reducir el umbral o revisar datos.")

    extras = {}
    # Imputación estilo notebook: NaN -> valor máximo global de CT (solo si política es 'nan')
    if not controles_df.empty and not muestras_df.empty and ok_for_fc:
        import pandas as pd
        # Coerción a numérico
        controles_df = controles_df.copy()
        muestras_df = muestras_df.copy()
        controles_df.loc[:, 'ct'] = pd.to_numeric(controles_df['ct'], errors='coerce')
        muestras_df.loc[:, 'ct'] = pd.to_numeric(muestras_df['ct'], errors='coerce')

        und_policy = st.session_state.get('und_policy', 'nan')
        if und_policy == 'nan':
            v_max = pd.concat([controles_df['ct'], muestras_df['ct']]).max()
            max_str = f"{v_max:.2f}" if pd.notna(v_max) else "NaN"
            st.info(f"Imputación posterior: NaN de Ct → {max_str} (máximo global)")
            controles_df.loc[:, 'ct'] = controles_df['ct'].fillna(v_max)
            muestras_df.loc[:, 'ct'] = muestras_df['ct'].fillna(v_max)
        else:
            st.info("Valores 'Undetermined' ya tratados durante la carga (sin imputación adicional).")

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
            # Aplicar preferencia de exclusión de 'estables' solo para gráficos
            exclude_stable = bool(st.session_state.get('exclude_stable', False))
            try:
                # Usar df_expr ya calculado para determinar objetivos a mostrar
                targets_plot = (
                    df_expr.loc[df_expr['nivel_expresion'] != 'estable', 'target']
                    if exclude_stable else df_expr['target']
                )
                consolidated_plot = fc.consolidated[fc.consolidated['target'].isin(targets_plot)]
            except Exception:
                consolidated_plot = fc.consolidated
            x_vals = consolidated_plot['target']
            ddct_mean = consolidated_plot['delta_delta_ct_promedio']
            ddct_ref = consolidated_plot['delta_delta_ct_gen_ref']
            fc_mean = consolidated_plot['fold_change_promedio']
            fc_ref = consolidated_plot['fold_change_gen_ref']

            y2_scale = st.radio("Escala FC", ["log", "lineal"], horizontal=True, index=0)
            y2_type = 'log' if y2_scale == 'log' else 'linear'

            fig.add_trace(go.Bar(
                x=x_vals, y=ddct_mean, name='ΔΔCT (Promedios)', marker_color='#1f77b4', opacity=0.85, yaxis='y',
                hovertemplate='Gen=%{x}<br>ΔΔCT Promedios=%{y:.3f}<extra></extra>'
            ))
            fig.add_trace(go.Bar(
                x=x_vals, y=ddct_ref, name='ΔΔCT (Gen Ref)', marker_color='#ff7f0e', opacity=0.85, yaxis='y',
                hovertemplate='Gen=%{x}<br>ΔΔCT GenRef=%{y:.3f}<extra></extra>'
            ))
            fig.add_trace(go.Scatter(
                x=x_vals, y=fc_mean, name='Fold Change (Promedios)', mode='markers+lines',
                marker=dict(color='#2ca02c', size=8, symbol='diamond'),
                line=dict(color='#2ca02c', width=2, dash='dot'), yaxis='y2',
                hovertemplate='Gen=%{x}<br>FC Promedios=%{y:.3f}<extra></extra>'
            ))
            fig.add_trace(go.Scatter(
                x=x_vals, y=fc_ref, name='Fold Change (Gen Ref)', mode='markers+lines',
                marker=dict(color='#d62728', size=8, symbol='diamond'),
                line=dict(color='#d62728', width=2, dash='dot'), yaxis='y2',
                hovertemplate='Gen=%{x}<br>FC GenRef=%{y:.3f}<extra></extra>'
            ))

            # Resaltar el gen de referencia
            try:
                ref_gene = fc.reference_gene
                ref_mask = x_vals == ref_gene
                fig.add_trace(go.Scatter(
                    x=x_vals[ref_mask], y=fc_ref[ref_mask], mode='markers', name='Ref gene',
                    marker=dict(color='black', size=12, symbol='star'), yaxis='y2',
                    hovertemplate='Gen de referencia=%{x}<br>FC GenRef=%{y:.3f}<extra></extra>'
                ))
            except Exception:
                pass

            fig.update_layout(
                title=dict(text='Análisis comparativo de métodos de cálculo', x=0.5),
                template='plotly_white', barmode='group',
                yaxis=dict(title='ΔΔCT', showgrid=True, gridcolor='lightgray'),
                yaxis2=dict(title=f"Fold Change ({y2_scale})", overlaying='y', side='right', type=y2_type, showgrid=False),
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
                # Aplicar exclusión de 'estables' en el gráfico de distribución si corresponde
                exclude_stable = bool(st.session_state.get('exclude_stable', False))
                df_expr_plot = df_expr[df_expr['nivel_expresion'] != 'estable'] if exclude_stable else df_expr
                counts = df_expr_plot['nivel_expresion'].value_counts().reindex(order_levels, fill_value=0)
                bar = px.bar(x=counts.index, y=counts.values, labels={'x': 'Nivel de expresión', 'y': 'Frecuencia'}, title='Distribución de niveles de expresión')
                st.plotly_chart(bar, use_container_width=True)

            # Genes por nivel de expresión: lista/tablas y treemap
            st.markdown("### Genes por nivel de expresión")
            # Base para visualizaciones (respetar preferencia de exclusión de 'estables')
            df_expr_vis = df_expr[df_expr['nivel_expresion'] != 'estable'].copy() if bool(st.session_state.get('exclude_stable', False)) else df_expr.copy()
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
                try:
                    import plotly.express as px
                    treedata = df_expr_vis.copy()
                    if treedata.empty:
                        st.info("Sin genes para graficar con los filtros actuales.")
                    else:
                        treedata['log2fc'] = np.log2(treedata['fold_change'].clip(lower=1e-12))
                        fig_t = px.treemap(
                            treedata,
                            path=['nivel_expresion', 'target'],
                            values=treedata['log2fc'].abs(),
                            color='log2fc',
                            color_continuous_scale='RdBu',
                            color_continuous_midpoint=0.0,
                            title='Genes por nivel (tamaño = |log2FC|, color = log2FC)'
                        )
                        st.plotly_chart(fig_t, use_container_width=True)
                except Exception as _:
                    st.info("No se pudo generar el treemap con los datos actuales.")

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
            exclude_stable_pref = bool(st.session_state.get('exclude_stable', False))
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
                species = st.selectbox("Especie (NCBI Taxon)", options=[9606], index=0)

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
        # Credenciales: intentar st.secrets si existen; tolerar ausencia de secrets.toml
        try:
            sec_email = st.secrets["NCBI_EMAIL"]
        except Exception:
            sec_email = ""
        try:
            sec_key = st.secrets["NCBI_API_KEY"]
        except Exception:
            sec_key = ""
        env_email = sec_email or os.getenv("NCBI_EMAIL", "")
        env_key = sec_key or os.getenv("NCBI_API_KEY", "")
        ccreds1, ccreds2 = st.columns([2, 2])
        with ccreds1:
            ncbi_email_input = st.text_input("NCBI Email (obligatorio)", value=env_email, placeholder="tu_email@dominio.com")
        with ccreds2:
            ncbi_api_key_input = st.text_input("NCBI API Key (opcional)", value=env_key, placeholder="api_key")

        max_per_gene = st.number_input("Máximo de artículos por gen", value=100, min_value=10, max_value=300, step=10)
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
                if bool(st.session_state.get('exclude_stable', False)):
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
                    if bool(st.session_state.get('exclude_stable', False)):
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
