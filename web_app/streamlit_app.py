from __future__ import annotations

import io
import logging
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import streamlit as st

# Asegura que el paquete ``web_app`` sea importable al ejecutar con ``streamlit run``
_PROJ_ROOT = Path(__file__).resolve().parents[1]
if str(_PROJ_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJ_ROOT))

from app.config.loader import ConfigError, load_app_config
from app.core.ensembl import add_ensembl_info_batch
from app.core.io import LoadResult, list_excel_sheets, parse_qpcr_wide
from app.core.qpcr import suggest_name_affixes
from app.core.reference_normalization import DifferentialExpressionResult
from app.services.fold_change import (
    FoldChangePreparationError,
    QualityMetrics,
    apply_undetermined_policy,
    compute_fold_change_with_expression,
    compute_quality_metrics,
)
from app.services.normalization import (
    AdvancedNormalizationError,
    AdvancedNormalizationParams,
    AdvancedNormalizationResult,
    build_method_matrices,
    execute_advanced_normalization,
)
from app.services.fold_change_visuals import (
    build_fold_change_chart,
    build_fold_change_table,
)
from app.services.heatmap_visuals import build_dendrogram_heatmap
from app.services.qpcr import (
    build_long_table,
    classify_tests_by_prefixes,
    classify_tests_by_regex,
    classify_tests_by_selection,
    classify_tests_by_suffixes,
    detect_collisions,
    summarize_extraction,
)
from app.state import AppSessionState


@st.cache_data(show_spinner=False, ttl=3600)
def _annotate_ensembl_cached(df_csv: str, max_workers: int = 3) -> pd.DataFrame:
    try:
        df = pd.read_csv(io.StringIO(df_csv))
    except Exception:
        return pd.DataFrame()
    if df.empty:
        return pd.DataFrame()
    return add_ensembl_info_batch(df, symbol_col="target", max_workers=max_workers)


def _format_list_preview(items: List[str], limit: int = 20) -> str:
    if not items:
        return "(sin coincidencias)"
    head = items[:limit]
    suffix = " …" if len(items) > limit else ""
    return ", ".join(head) + suffix


def _store_classification(
    state: Dict[str, object],
    ctrl_df: pd.DataFrame,
    samp_df: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    ctrl_df = ctrl_df.copy()
    samp_df = samp_df.copy()
    collisions = detect_collisions(ctrl_df, samp_df)
    if collisions:
        preview = _format_list_preview(collisions, limit=10)
        st.warning(
            "Algunas pruebas aparecen tanto en controles como en muestras. "
            "Se priorizarán como controles: "
            f"{preview}"
        )
        samp_df = samp_df[~samp_df["test"].astype(str).isin(collisions)].copy()
    state["controles_df"] = ctrl_df
    state["muestras_df"] = samp_df
    return ctrl_df, samp_df


def _clear_classification_state(state: Dict[str, object]) -> None:
    for key in (
        "ctrl_prefix",
        "samp_prefix",
        "ctrl_suffix",
        "samp_suffix",
        "ctrl_regex",
        "samp_regex",
        "selected_ctrl",
        "selected_samp",
        "controles_df",
        "muestras_df",
    ):
        state.pop(key, None)
    # Limpiar entradas asociadas en session_state de Streamlit
    for ui_key in list(st.session_state.keys()):
        if any(ui_key.startswith(prefix) for prefix in (
            "ctrl_prefix_input::",
            "samp_prefix_input::",
            "ctrl_suffix_input::",
            "samp_suffix_input::",
            "pref_ctrl_suggest::",
            "pref_samp_suggest::",
            "suff_ctrl_suggest::",
            "suff_samp_suggest::",
        )):
            st.session_state.pop(ui_key, None)


def _render_classification_ui(long_df: pd.DataFrame, file_key: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    state = st.session_state.setdefault(file_key, {})
    unique_tests = sorted(long_df["test"].astype(str).dropna().unique().tolist())
    sugg = suggest_name_affixes(unique_tests, top_n=10)
    suggested_prefixes = [p for p, _ in (sugg.get("prefixes") or [])]
    suggested_suffixes = [s for s, _ in (sugg.get("suffixes") or [])]

    tab_pref, tab_suff, tab_regex, tab_manual = st.tabs(
        ["Por prefijos", "Por sufijos", "Por regex", "Selección manual"]
    )

    with tab_pref:
        col_pref_ctrl, col_pref_samp = st.columns(2)
        ctrl_suggest = ""
        samp_suggest = ""
        if suggested_prefixes:
            ctrl_options = ["(vacío)"] + suggested_prefixes
            samp_options = ["(vacío)"] + suggested_prefixes
            with col_pref_ctrl:
                ctrl_suggest = st.selectbox(
                    "Sugerencia controles",
                    options=ctrl_options,
                    index=0,
                    key=f"pref_ctrl_suggest::{file_key}",
                )
            with col_pref_samp:
                samp_suggest = st.selectbox(
                    "Sugerencia muestras",
                    options=samp_options,
                    index=0,
                    key=f"pref_samp_suggest::{file_key}",
                )
        with col_pref_ctrl:
            default_ctrl = state.get("ctrl_prefix") or ("" if ctrl_suggest in {"", "(vacío)"} else ctrl_suggest)
            ctrl_prefix = st.text_input(
                "Prefijo controles",
                value=default_ctrl,
                help="Introduce un prefijo común de los tests de control.",
                key=f"ctrl_prefix_input::{file_key}",
            )
        with col_pref_samp:
            default_samp = state.get("samp_prefix") or ("" if samp_suggest in {"", "(vacío)"} else samp_suggest)
            samp_prefix = st.text_input(
                "Prefijo muestras",
                value=default_samp,
                help="Introduce un prefijo común de los tests de muestra.",
                key=f"samp_prefix_input::{file_key}",
            )
        if suggested_prefixes:
            st.caption("Sugerencias detectadas: " + ", ".join(suggested_prefixes))
        prev_ctrl = [t for t in unique_tests if ctrl_prefix and str(t).startswith(ctrl_prefix)]
        prev_samp = [t for t in unique_tests if samp_prefix and str(t).startswith(samp_prefix)]
        st.caption(f"Controles previstos: {_format_list_preview(prev_ctrl)}")
        st.caption(f"Muestras previstas: {_format_list_preview(prev_samp)}")
        if st.button("Clasificar por prefijos", key=f"btn_pref_{file_key}"):
            if not ctrl_prefix and not samp_prefix:
                st.warning("Introduce al menos un prefijo para clasificar.")
            else:
                result = classify_tests_by_prefixes(
                    long_df,
                    [ctrl_prefix] if ctrl_prefix else [],
                    [samp_prefix] if samp_prefix else [],
                )
                ctrl_df, samp_df = _store_classification(state, result.controles, result.muestras)
                state["ctrl_prefix"] = ctrl_prefix
                state["samp_prefix"] = samp_prefix
                st.success(
                    f"Clasificación aplicada → controles: {len(ctrl_df)}, muestras: {len(samp_df)}"
                )

    with tab_suff:
        col_suff_ctrl, col_suff_samp = st.columns(2)
        suff_suggest_ctrl = ""
        suff_suggest_samp = ""
        if suggested_suffixes:
            suff_ctrl_options = ["(vacío)"] + suggested_suffixes
            suff_samp_options = ["(vacío)"] + suggested_suffixes
            with col_suff_ctrl:
                suff_suggest_ctrl = st.selectbox(
                    "Sugerencia controles",
                    options=suff_ctrl_options,
                    index=0,
                    key=f"suff_ctrl_suggest::{file_key}",
                )
            with col_suff_samp:
                suff_suggest_samp = st.selectbox(
                    "Sugerencia muestras",
                    options=suff_samp_options,
                    index=0,
                    key=f"suff_samp_suggest::{file_key}",
                )
        with col_suff_ctrl:
            default_ctrl_suff = state.get("ctrl_suffix") or ("" if suff_suggest_ctrl in {"", "(vacío)"} else suff_suggest_ctrl)
            ctrl_suffix = st.text_input(
                "Sufijos controles (separados por coma)",
                value=default_ctrl_suff,
                key=f"ctrl_suffix_input::{file_key}",
            )
        with col_suff_samp:
            default_samp_suff = state.get("samp_suffix") or ("" if suff_suggest_samp in {"", "(vacío)"} else suff_suggest_samp)
            samp_suffix = st.text_input(
                "Sufijos muestras (separados por coma)",
                value=default_samp_suff,
                key=f"samp_suffix_input::{file_key}",
            )
        if suggested_suffixes:
            st.caption("Sugerencias detectadas: " + ", ".join(suggested_suffixes))
        ctrl_suffixes = [s.strip() for s in ctrl_suffix.split(",") if s.strip()]
        samp_suffixes = [s.strip() for s in samp_suffix.split(",") if s.strip()]
        prev_ctrl = [t for t in unique_tests if any(str(t).endswith(s) for s in ctrl_suffixes)]
        prev_samp = [t for t in unique_tests if any(str(t).endswith(s) for s in samp_suffixes)]
        st.caption(f"Controles previstos: {_format_list_preview(prev_ctrl)}")
        st.caption(f"Muestras previstas: {_format_list_preview(prev_samp)}")
        if st.button("Clasificar por sufijos", key=f"btn_suff_{file_key}"):
            if not ctrl_suffixes and not samp_suffixes:
                st.warning("Define al menos un sufijo para clasificar.")
            else:
                result = classify_tests_by_suffixes(long_df, ctrl_suffixes, samp_suffixes)
                ctrl_df, samp_df = _store_classification(state, result.controles, result.muestras)
                state["ctrl_suffix"] = ctrl_suffix
                state["samp_suffix"] = samp_suffix
                st.success(
                    f"Clasificación aplicada → controles: {len(ctrl_df)}, muestras: {len(samp_df)}"
                )

    with tab_regex:
        col_regex_ctrl, col_regex_samp = st.columns(2)
        with col_regex_ctrl:
            ctrl_regex = st.text_input(
                "Regex controles",
                value=state.get("ctrl_regex", ""),
            )
        with col_regex_samp:
            samp_regex = st.text_input(
                "Regex muestras",
                value=state.get("samp_regex", ""),
            )
        try:
            ctrl_pattern = re.compile(ctrl_regex) if ctrl_regex else None
        except re.error:
            ctrl_pattern = None
            st.warning("Regex de controles inválido, actualízalo antes de clasificar.")
        try:
            samp_pattern = re.compile(samp_regex) if samp_regex else None
        except re.error:
            samp_pattern = None
            st.warning("Regex de muestras inválido, actualízalo antes de clasificar.")
        prev_ctrl = [t for t in unique_tests if ctrl_pattern and ctrl_pattern.search(str(t))]
        prev_samp = [t for t in unique_tests if samp_pattern and samp_pattern.search(str(t))]
        st.caption(f"Controles previstos: {_format_list_preview(prev_ctrl)}")
        st.caption(f"Muestras previstas: {_format_list_preview(prev_samp)}")
        if st.button("Clasificar por regex", key=f"btn_regex_{file_key}"):
            if not ctrl_regex and not samp_regex:
                st.warning("Introduce al menos un patrón regex para clasificar.")
            else:
                try:
                    result = classify_tests_by_regex(long_df, ctrl_regex, samp_regex)
                except re.error as exc:  # type: ignore[name-defined]
                    st.error(f"Regex inválido: {exc}")
                else:
                    ctrl_df, samp_df = _store_classification(state, result.controles, result.muestras)
                    state["ctrl_regex"] = ctrl_regex
                    state["samp_regex"] = samp_regex
                    st.success(
                        f"Clasificación aplicada → controles: {len(ctrl_df)}, muestras: {len(samp_df)}"
                    )

    with tab_manual:
        col_sel_ctrl, col_sel_samp = st.columns(2)
        with col_sel_ctrl:
            selected_ctrl = st.multiselect(
                "Pruebas controles",
                options=unique_tests,
                default=state.get("selected_ctrl", []),
            )
        with col_sel_samp:
            selected_samp = st.multiselect(
                "Pruebas muestras",
                options=unique_tests,
                default=state.get("selected_samp", []),
            )
        if st.button("Clasificar selección", key=f"btn_sel_{file_key}"):
            if not selected_ctrl and not selected_samp:
                st.warning("Selecciona al menos una prueba en cualquiera de los grupos.")
            else:
                result = classify_tests_by_selection(long_df, selected_ctrl, selected_samp)
                ctrl_df, samp_df = _store_classification(state, result.controles, result.muestras)
                state["selected_ctrl"] = selected_ctrl
                state["selected_samp"] = selected_samp
                st.success(
                    f"Clasificación aplicada → controles: {len(ctrl_df)}, muestras: {len(samp_df)}"
                )

    st.button(
        "Limpiar clasificación",
        key=f"btn_clear_{file_key}",
        on_click=_clear_classification_state,
        kwargs={"state": state},
    )

    controles_df = state.get("controles_df")
    muestras_df = state.get("muestras_df")
    if not isinstance(controles_df, pd.DataFrame):
        controles_df = pd.DataFrame()
    if not isinstance(muestras_df, pd.DataFrame):
        muestras_df = pd.DataFrame()
    return controles_df, muestras_df


def _quality_summary(metrics: QualityMetrics) -> None:
    cols = st.columns(4)
    cols[0].metric("Genes controles", len(metrics.ctrl_targets))
    cols[1].metric("Genes muestras", len(metrics.sample_targets))
    cols[2].metric("Genes comunes", len(metrics.common_targets))
    cols[3].metric(
        "NaN Ct (ctrl/mues)",
        f"{metrics.ctrl_nan_ratio:.0%} / {metrics.sample_nan_ratio:.0%}"
    )


def _build_expression_summary(
    diff: DifferentialExpressionResult,
    alpha: float,
) -> pd.DataFrame:
    stats = diff.stats.copy()
    if stats.empty:
        return pd.DataFrame(
            columns=[
                "target",
                "fold_change",
                "log2_fc",
                "mean_ctrl",
                "mean_case",
                "q",
                "bootstrap_freq",
                "nivel_expresion",
                "significativo",
            ]
        )
    df = stats.rename(columns={"gene": "target"}).copy()
    if "bootstrap_freq" not in df.columns:
        df["bootstrap_freq"] = diff.bootstrap_freq.reindex(df["target"]).fillna(0.0).values
    df["log2_fc"] = df["mean_case"] - df["mean_ctrl"]
    df["fold_change"] = np.power(2.0, df["log2_fc"])
    df["nivel_expresion"] = pd.cut(
        df["fold_change"],
        bins=[-np.inf, 1.0, 2.0, np.inf],
        labels=["subexpresado", "estable", "sobreexpresado"],
        right=False,
    )
    df["significativo"] = df["q"] < float(alpha)
    ordered = df.sort_values(["significativo", "q"], ascending=[False, True])
    columns = [
        "target",
        "fold_change",
        "log2_fc",
        "mean_ctrl",
        "mean_case",
        "q",
        "bootstrap_freq",
        "nivel_expresion",
        "significativo",
    ]
    return ordered[columns]


def _render_advanced_results(
    result: AdvancedNormalizationResult,
    params: AdvancedNormalizationParams,
) -> pd.DataFrame:
    st.subheader("Resultados de la normalización avanzada")
    summary = result.summary()
    if summary:
        metric_cols = st.columns(len(summary))
        for (label, value), col in zip(summary.items(), metric_cols):
            with col:
                display = ", ".join(map(str, value)) if isinstance(value, (list, tuple)) else value
                st.metric(label.replace("_", " "), display)

    ref_res = result.reference_result
    if not ref_res.candidate_scores.empty:
        candidates = ref_res.candidate_scores.copy()
        candidates["refs"] = candidates["refs"].apply(lambda ref: ", ".join(ref))
        st.markdown("**Ranking de combinaciones candidata**")
        st.dataframe(candidates[["refs", "score"]].head(30))
        st.download_button(
            "Descargar ranking (CSV)",
            candidates.to_csv(index=False),
            file_name="normalizacion_avanzada_candidatos.csv",
            use_container_width=True,
        )

    st.markdown("**Tabla normalizada (df_norm)**")
    st.dataframe(result.df_norm.head(100))
    st.download_button(
        "Descargar df_norm (CSV)",
        result.df_norm.to_csv(index=False),
        file_name="normalizacion_avanzada_df_norm.csv",
        use_container_width=True,
    )

    stats_df = result.differential.stats.copy()
    if not stats_df.empty:
        st.markdown("**Estadísticas diferenciales**")
        st.dataframe(stats_df.head(50))
        st.download_button(
            "Descargar estadísticas (CSV)",
            stats_df.to_csv(index=False),
            file_name="normalizacion_avanzada_stats.csv",
            use_container_width=True,
        )
    else:
        st.info("No se obtuvieron estadísticas significativas con los parámetros actuales.")

    heat_df = result.df_heatmap
    if isinstance(heat_df, pd.DataFrame) and not heat_df.empty:
        try:
            fig_heat = build_dendrogram_heatmap(heat_df, title="Avanzada — Heatmap con dendrogramas")
            st.plotly_chart(fig_heat, use_container_width=True)
        except Exception:
            st.dataframe(heat_df)

    # (comparativo 3 métodos se realiza en la vista principal, donde están controles/muestras)

    expr_summary = _build_expression_summary(result.differential, params.alpha)
    if expr_summary.empty:
        st.info("No se generó un resumen de expresión. Ajusta parámetros o revisa la clasificación.")
    else:
        st.markdown("**Resumen de expresión diferenciada**")
        st.dataframe(expr_summary.head(50))
        st.download_button(
            "Descargar resumen (CSV)",
            expr_summary.to_csv(index=False),
            file_name="normalizacion_avanzada_resumen_expresion.csv",
            use_container_width=True,
        )
        counts = expr_summary["nivel_expresion"].value_counts()
        st.success(
            "Resultados: "
            f"{counts.get('subexpresado', 0)} subexpresados · "
            f"{counts.get('estable', 0)} estables · "
            f"{counts.get('sobreexpresado', 0)} sobreexpresados"
        )
    return expr_summary


def _render_ensembl_section(expr_summary: pd.DataFrame) -> None:
    st.subheader("Anotación Ensembl (IDs génicos)")
    if expr_summary is None or expr_summary.empty:
        st.info("Ejecuta la normalización avanzada para generar un listado de genes.")
        return

    genes_df = (
        expr_summary[["target", "nivel_expresion", "fold_change"]]
        .drop_duplicates(subset=["target"])
        .sort_values("target")
        .reset_index(drop=True)
    )
    if genes_df.empty:
        st.info("No hay genes disponibles para anotar.")
        return

    genes_key = ",".join(genes_df["target"].astype(str).tolist())
    session_key = st.session_state.get("_ensembl_genes_key")
    cached_df = st.session_state.get("_ensembl_df")

    col_left, col_right = st.columns([1, 1])
    with col_left:
        workers = st.number_input(
            "Hilos simultáneos",
            min_value=1,
            max_value=16,
            value=3,
            step=1,
            help="Número máximo de hilos para consultas Ensembl.",
        )
        need_refresh = (
            session_key != genes_key
            or not isinstance(cached_df, pd.DataFrame)
            or cached_df.empty
        )
        if need_refresh:
            if st.button("Consultar Ensembl", key="btn_run_ensembl"):
                with st.spinner("Consultando Ensembl…"):
                    annotated = _annotate_ensembl_cached(
                        genes_df.to_csv(index=False),
                        max_workers=int(workers),
                    )
                st.session_state["_ensembl_df"] = annotated
                st.session_state["_ensembl_genes_key"] = genes_key
                cached_df = annotated
        else:
            st.success("Usando anotación almacenada para esta lista de genes.")
            if st.button("Recalcular anotación", key="btn_refresh_ensembl"):
                with st.spinner("Consultando Ensembl…"):
                    annotated = _annotate_ensembl_cached(
                        genes_df.to_csv(index=False),
                        max_workers=int(workers),
                    )
                st.session_state["_ensembl_df"] = annotated
                st.session_state["_ensembl_genes_key"] = genes_key
                cached_df = annotated

    if isinstance(cached_df, pd.DataFrame) and not cached_df.empty:
        annotated = cached_df.copy()
        annotated["has_desc"] = (
            annotated["description"].fillna("").astype(str).str.strip().ne("")
        )
        tot = len(annotated)
        with col_right:
            st.metric("Genes consultados", tot)
            st.metric("Con Ensembl ID", int((annotated["ensembl_id"] != "Not found").sum()))
            st.metric("Con descripción", int(annotated["has_desc"].sum()))

        st.dataframe(annotated)
        st.download_button(
            "Descargar anotaciones (CSV)",
            annotated.to_csv(index=False),
            file_name="ensembl_anotado.csv",
            use_container_width=True,
        )

        st.markdown("**Enlaces rápidos**")
        top = annotated.head(100)
        for _, row in top.iterrows():
            gene = str(row.get("target", ""))
            eid = str(row.get("ensembl_id", ""))
            desc = str(row.get("description", ""))
            lvl = str(row.get("nivel_expresion", ""))
            if eid and eid != "Not found":
                url = f"https://www.ensembl.org/Homo_sapiens/Gene/Summary?g={eid}"
                st.markdown(f"- [{gene}]({url}) · {lvl} — {desc}")
            else:
                st.markdown(f"- {gene} · {lvl} — {desc}")
    else:
        st.info("Ejecuta la consulta para obtener IDs Ensembl.")


def main() -> None:
    st.set_page_config(
        page_title="CancerGeneSignatures",
        layout="wide",
        initial_sidebar_state="expanded",
    )

    warnings: List[str] = []
    try:
        secrets_source = st.secrets  # type: ignore[attr-defined]
    except Exception:
        secrets_source = None

    try:
        app_config = load_app_config(
            secrets=secrets_source if secrets_source is not None else None,
            warn=warnings.append,
        )
    except ConfigError as exc:
        st.error(f"Error cargando configuración: {exc}")
        st.stop()
        return

    if "_log_configured" not in st.session_state:
        level_name = getattr(app_config, "log_level", "INFO")
        level = getattr(logging, str(level_name).upper(), logging.INFO)
        logging.basicConfig(
            level=level,
            format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        )
        st.session_state["_log_configured"] = True
    logger = logging.getLogger("cgs.web_app")

    for warning_msg in warnings:
        st.warning(warning_msg)

    app_state = AppSessionState.load()

    st.title("Flujo qPCR → Normalización avanzada → IDs Ensembl")
    st.write(
        "Carga tu hoja qPCR, clasifica controles y muestras, ejecuta la normalización avanzada "
        "y obtén IDs Ensembl con descripciones asociadas."
    )

    with st.sidebar:
        st.header("1) Datos de entrada")
        uploaded = st.file_uploader("Archivo Excel", type=["xlsx", "xls"])
        sheet: Optional[str] = None
        if uploaded is not None:
            uploaded.seek(0)
            try:
                sheets = list_excel_sheets(uploaded)
            except Exception as exc:
                st.warning(f"No se pudieron listar hojas: {exc}")
                sheets = []
            if sheets:
                sheet = st.selectbox("Hoja", options=sheets, index=0)

        st.header("2) Valores 'Undetermined/ND'")
        und_policy = st.selectbox(
            "Política",
            options=["nan", "ctmax", "value"],
            index=["nan", "ctmax", "value"].index(app_state.undetermined.policy)
            if app_state.undetermined.policy in {"nan", "ctmax", "value"}
            else 0,
        )
        und_value = st.number_input(
            "Valor fijo (si aplica)",
            value=float(app_state.undetermined.value),
            min_value=0.0,
            max_value=100.0,
            step=0.5,
        )

        st.header("3) Opciones del estudio")
        cancer_types = app_config.menu.cancer_types
        if not cancer_types:
            st.error("El menú de configuración no define tipos de cáncer disponibles.")
            cancer_types = ["Sin definir"]
        default_cancer = (
            cancer_types.index(app_state.cancer_type)
            if app_state.cancer_type in cancer_types
            else 0
        )
        cancer_selected = st.selectbox(
            "Tipo de cáncer",
            options=cancer_types,
            index=default_cancer,
        )

        context_options = [ctx.label for ctx in app_config.menu.contexts]
        default_context_idx = (
            context_options.index(app_state.context_label)
            if app_state.context_label in context_options
            else 0
        )
        context_selected = st.selectbox(
            "Contexto biológico",
            options=context_options,
            index=default_context_idx,
        )

        run_btn = st.button("Procesar archivo", type="primary")

    app_state.undetermined.policy = und_policy
    app_state.undetermined.value = float(und_value)
    app_state.cancer_type = str(cancer_selected)
    app_state.context_label = str(context_selected)
    app_state.persist()

    df_loaded: Optional[LoadResult] = app_state.df_loaded

    if uploaded is not None and run_btn:
        try:
            uploaded.seek(0)
            df_loaded = parse_qpcr_wide(
                uploaded,
                sheet_name=sheet,
                header_mode="coords",
                header_row_idx=3,
                well_col_idx=0,
                target_col_idx=1,
                undetermined_policy=und_policy,
                undetermined_value=float(und_value),
            )
        except Exception as exc:
            st.warning(f"Cabecera esperada no encontrada ({exc}); intentando detección automática…")
            uploaded.seek(0)
            try:
                df_loaded = parse_qpcr_wide(
                    uploaded,
                    sheet_name=sheet,
                    header_mode="auto",
                    undetermined_policy=und_policy,
                    undetermined_value=float(und_value),
                )
            except Exception as exc_auto:
                st.error(f"No se pudo leer el archivo: {exc_auto}")
                logger.exception("Fallo al parsear archivo", exc_info=exc_auto)
                df_loaded = None
        if df_loaded is not None:
            app_state.df_loaded = df_loaded
            app_state.persist()
            logger.info(
                "Archivo cargado: %s | hoja=%s | forma=%s",
                df_loaded.source_name,
                df_loaded.sheet_name,
                df_loaded.df.shape,
            )

    if df_loaded is None:
        st.info("Carga un archivo Excel y pulsa 'Procesar archivo' para comenzar.")
        return

    st.subheader("Vista previa de datos")
    st.caption(
        f"Archivo: {df_loaded.source_name} | Hoja: {df_loaded.sheet_name or '-'} | "
        f"Forma: {df_loaded.df.shape}"
    )
    st.dataframe(df_loaded.df.head(20))

    st.subheader("Resumen de extracción")
    try:
        extraction = summarize_extraction(df_loaded)
        st.info(
            f"Pruebas detectadas ({len(extraction.sample_names)}): "
            f"{', '.join(extraction.sample_names)}"
        )
        st.info(
            f"Genes objetivo ({len(extraction.genes)}): {', '.join(extraction.genes)}"
        )
        st.info(
            f"Pozos detectados ({len(extraction.wells)}): {', '.join(extraction.wells)}"
        )
    except Exception as exc:
        st.warning(f"No se pudo generar el resumen: {exc}")

    long_df, controls_warning = build_long_table(df_loaded)
    if controls_warning:
        st.warning(f"No se filtraron los controles de máquina: {controls_warning}")

    st.markdown("**Configuración actual del estudio**")
    col_conf1, col_conf2 = st.columns(2)
    with col_conf1:
        st.info(f"Tipo de cáncer: {app_state.cancer_type or app_config.menu.cancer_types[0]}")
    with col_conf2:
        st.info(f"Contexto biológico: {app_state.context_label or app_config.menu.contexts[0].label}")

    st.divider()
    st.subheader("Clasificación de controles vs muestras")
    file_key = f"assign::{df_loaded.source_name}:{df_loaded.sheet_name}"
    controles_df, muestras_df = _render_classification_ui(long_df, file_key)

    st.write(f"Controles clasificados: {len(controles_df)} filas")
    st.write(f"Muestras clasificadas: {len(muestras_df)} filas")

    if controles_df.empty or muestras_df.empty:
        st.info("Clasifica al menos una prueba como control y una como muestra para continuar.")
        return

    try:
        imputation = apply_undetermined_policy(
            controles_df,
            muestras_df,
            policy=und_policy,
            fixed_value=float(und_value),
            group_columns=[c for c in ["plate_id", "target"] if c in controles_df.columns and c in muestras_df.columns] or None,
        )
    except FoldChangePreparationError as exc:
        st.error(f"No se pudo aplicar la política 'Undetermined': {exc}")
        return

    controles_df = imputation.controles
    muestras_df = imputation.muestras
    if imputation.message:
        st.caption(imputation.message)
    if not imputation.summary.empty:
        st.markdown("**Detalle de imputaciones**")
        st.dataframe(imputation.summary.head(50))
        st.download_button(
            "Descargar imputaciones (CSV)",
            imputation.summary.to_csv(index=False),
            file_name="detalle_imputaciones.csv",
            use_container_width=True,
        )

    st.divider()
    st.subheader("Calidad de datos")
    metrics = compute_quality_metrics(controles_df, muestras_df)
    _quality_summary(metrics)
    if not metrics.common_targets:
        st.error("No hay genes en común entre controles y muestras. Revisa la clasificación.")
        return

    st.divider()
    st.subheader("ΔΔCT y fold change")
    try:
        fold_result = compute_fold_change_with_expression(controles_df, muestras_df)
    except ValueError as exc:
        st.error(f"No se pudo calcular el fold change: {exc}")
        return

    fc_core = fold_result.computation
    st.info(f"Gen de referencia seleccionado: {fc_core.reference_gene}")

    st.markdown("**Tabla consolidada de métricas**")
    st.dataframe(fc_core.consolidated.head(100))
    st.download_button(
        "Descargar consolidado (CSV)",
        fc_core.consolidated.to_csv(index=False),
        file_name="fold_change_consolidado.csv",
        use_container_width=True,
    )

    col_prom, col_ref = st.columns(2)
    with col_prom:
        st.markdown("**Método: Promedio global**")
        st.dataframe(fc_core.by_means.head(100))
        st.download_button(
            "Descargar promedios (CSV)",
            fc_core.by_means.to_csv(index=False),
            file_name="fold_change_promedios.csv",
            use_container_width=True,
        )
    with col_ref:
        st.markdown("**Método: Gen de referencia**")
        st.dataframe(fc_core.by_refgene.head(100))
        st.download_button(
            "Descargar gen de referencia (CSV)",
            fc_core.by_refgene.to_csv(index=False),
            file_name="fold_change_gen_referencia.csv",
            use_container_width=True,
        )

    st.markdown("**Comparativa visual**")
    table_fig = build_fold_change_table(fc_core.consolidated, reference_gene=fc_core.reference_gene)
    st.plotly_chart(table_fig, use_container_width=True)

    chart_fig = build_fold_change_chart(fc_core.consolidated)
    st.plotly_chart(chart_fig, use_container_width=True)

    st.markdown("**Clasificación por nivel de expresión**")
    st.dataframe(fold_result.expression_table.head(100))
    st.download_button(
        "Descargar clasificación (CSV)",
        fold_result.expression_table.to_csv(index=False),
        file_name="fold_change_categorizacion.csv",
        use_container_width=True,
    )

    st.divider()
    st.subheader("Normalización avanzada")
    adv_key = f"adv::{file_key}"
    adv_state = st.session_state.setdefault(adv_key, {})

    with st.expander("Parámetros", expanded=False):
        alpha = st.number_input(
            "Umbral FDR (q)",
            min_value=0.001,
            max_value=0.5,
            value=float(getattr(adv_state.get("params"), "alpha", 0.05)),
            step=0.005,
            format="%.3f",
            key=f"adv_alpha_{adv_key}",
        )
        col_a, col_b, col_c = st.columns(3)
        with col_a:
            n_candidates = st.slider(
                "Genes candidatos",
                min_value=4,
                max_value=30,
                value=int(getattr(adv_state.get("params"), "n_candidates", 20)),
            )
        with col_b:
            k_refs = st.selectbox(
                "Tamaño set referencia",
                options=[1, 2, 3],
                index=[1, 2, 3].index(int(getattr(adv_state.get("params"), "k_refs", 2))),
            )
        with col_c:
            random_seed = st.number_input(
                "Semilla aleatoria",
                value=int(getattr(adv_state.get("params"), "random_seed", 123) or 123),
                step=1,
            )
        col_boot, col_perm = st.columns(2)
        with col_boot:
            bootstrap_iter = st.number_input(
                "Iteraciones bootstrap",
                min_value=0,
                max_value=1000,
                value=int(getattr(adv_state.get("params"), "bootstrap_iter", 200)),
                step=50,
            )
        with col_perm:
            permutation_iter = st.number_input(
                "Permutaciones",
                min_value=0,
                max_value=600,
                value=int(getattr(adv_state.get("params"), "permutation_iter", 100)),
                step=20,
            )

    params = AdvancedNormalizationParams(
        alpha=float(alpha),
        n_candidates=int(n_candidates),
        k_refs=int(k_refs),
        bootstrap_iter=int(bootstrap_iter),
        permutation_iter=int(permutation_iter),
        random_seed=int(random_seed),
    )

    if st.button("Ejecutar normalización avanzada", key=f"btn_run_adv_{adv_key}"):
        try:
            result = execute_advanced_normalization(controles_df, muestras_df, params=params)
        except AdvancedNormalizationError as exc:
            st.error(f"La normalización avanzada falló: {exc}")
        else:
            adv_state["result"] = result
            adv_state["params"] = params
            st.success("Normalización avanzada completada.")

    result = adv_state.get("result")
    if not isinstance(result, AdvancedNormalizationResult):
        st.info("Configura los parámetros y ejecuta la normalización para ver resultados.")
        return

    expr_summary = _render_advanced_results(result, adv_state.get("params", params))
    st.session_state["_expression_summary"] = expr_summary

    st.divider()
    # Comparativo de métodos: avanzada vs promedio global vs gen de referencia básico
    st.subheader("Comparativo de expresión relativa (3 métodos)")
    try:
        # Selección de conjunto de genes y opciones de visualización
        genes_all = sorted(result.df_norm["target"].dropna().astype(str).unique().tolist())
        col_opts1, col_opts2 = st.columns([2, 1])
        with col_opts1:
            genes_mode = st.radio(
                "Genes para heatmap",
                options=["Significativos/Top-20", "Panel completo (todos)"],
                horizontal=True,
                help="Usa el panel completo para visualizar los 89 genes (o todos los disponibles).",
            )
        with col_opts2:
            zscore_rows = st.checkbox(
                "Estandarizar por gen (z-score)",
                value=True,
                help="Normaliza cada fila para resaltar patrones relativos entre tests.",
            )

        genes_for_heatmap = genes_all if genes_mode.endswith("(todos)") else None

        matrices, basic_ref_gene = build_method_matrices(
            controles_df,
            muestras_df,
            result,
            genes_for_heatmap=genes_for_heatmap,
        )

        refs_adv = ", ".join(result.reference_result.references)
        if len(result.reference_result.references) == 1 and refs_adv == str(basic_ref_gene):
            st.caption(
                "Las referencias avanzadas coinciden con el gen de referencia básico (K=1), "
                "por eso las visualizaciones pueden verse muy similares."
            )

        tab_adv, tab_gm, tab_ref = st.tabs(["Avanzada", "Promedio global", f"Gen ref básico ({basic_ref_gene})"])
        with tab_adv:
            fig_a = build_dendrogram_heatmap(
                matrices["advanced"],
                title="Avanzada — z-score por gen" if zscore_rows else "Avanzada",
                zscore_by_gene=zscore_rows,
            )
            st.plotly_chart(fig_a, use_container_width=True)
            st.download_button(
                "Descargar matriz avanzada (CSV)",
                matrices["advanced"].to_csv(),
                file_name="heatmap_avanzada.csv",
                use_container_width=True,
            )
        with tab_gm:
            fig_g = build_dendrogram_heatmap(
                matrices["global_mean"],
                title="Promedio global — z-score por gen" if zscore_rows else "Promedio global",
                zscore_by_gene=zscore_rows,
            )
            st.plotly_chart(fig_g, use_container_width=True)
            st.download_button(
                "Descargar matriz promedio global (CSV)",
                matrices["global_mean"].to_csv(),
                file_name="heatmap_promedio_global.csv",
                use_container_width=True,
            )
        with tab_ref:
            fig_r = build_dendrogram_heatmap(
                matrices["refgene"],
                title=f"Gen de referencia ({basic_ref_gene}) — z-score por gen" if zscore_rows else f"Gen de referencia ({basic_ref_gene})",
                zscore_by_gene=zscore_rows,
            )
            st.plotly_chart(fig_r, use_container_width=True)
            st.download_button(
                "Descargar matriz gen de referencia (CSV)",
                matrices["refgene"].to_csv(),
                file_name="heatmap_gen_referencia.csv",
                use_container_width=True,
            )
    except Exception as exc:
        st.info(f"No se pudo construir el comparativo de métodos: {exc}")

    st.divider()
    _render_ensembl_section(expr_summary)


if __name__ == "__main__":
    main()
