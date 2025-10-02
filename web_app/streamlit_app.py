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
    build_combined_ddct_fc_tables,
    build_expression_datasets,
    build_heatmaps_by_method,
    build_method_matrices,
    execute_advanced_normalization,
    save_gene_sets,
)
from app.services.fold_change_visuals import (
    build_fc_methods_scatter,
    build_expression_classification_table,
    build_classification_summary_chart,
    summarize_fc_methods,
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


def _notify_parameter_change(state_key: str, value: object, message: str) -> None:
    """Emite un mensaje informativo cuando un parámetro interactivo cambia."""

    feedback_state = st.session_state.setdefault("_param_feedback", {})
    if feedback_state.get(state_key) != value:
        feedback_state[state_key] = value
        st.info(message)


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
        st.dataframe(candidates[["refs", "score"]].head(30), width="stretch")
        st.download_button(
            "Descargar ranking (CSV)",
            candidates.to_csv(index=False),
            file_name="normalizacion_avanzada_candidatos.csv",
            width="stretch",
        )

    st.markdown("**Tabla normalizada (df_norm)**")
    st.dataframe(result.df_norm.head(100), width="stretch")
    st.download_button(
        "Descargar df_norm (CSV)",
        result.df_norm.to_csv(index=False),
        file_name="normalizacion_avanzada_df_norm.csv",
        width="stretch",
    )

    stats_df = result.differential.stats.copy()
    if not stats_df.empty:
        st.markdown("**Estadísticas diferenciales**")
        st.dataframe(stats_df.head(50), width="stretch")
        st.download_button(
            "Descargar estadísticas (CSV)",
            stats_df.to_csv(index=False),
            file_name="normalizacion_avanzada_stats.csv",
            width="stretch",
        )
    else:
        st.info("No se obtuvieron estadísticas significativas con los parámetros actuales.")

    heat_df = result.df_heatmap
    if isinstance(heat_df, pd.DataFrame) and not heat_df.empty:
        try:
            fig_heat = build_dendrogram_heatmap(heat_df, title="Avanzada — Heatmap con dendrogramas")
            st.plotly_chart(fig_heat, width="stretch")
        except Exception:
            st.dataframe(heat_df, width="stretch")

    # (comparativo 3 métodos se realiza en la vista principal, donde están controles/muestras)

    expr_summary = _build_expression_summary(result.differential, params.alpha)
    if expr_summary.empty:
        st.info("No se generó un resumen de expresión. Ajusta parámetros o revisa la clasificación.")
    else:
        st.markdown("**Resumen de expresión diferenciada**")
        st.dataframe(expr_summary.head(50), width="stretch")
        st.download_button(
            "Descargar resumen (CSV)",
            expr_summary.to_csv(index=False),
            file_name="normalizacion_avanzada_resumen_expresion.csv",
            width="stretch",
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

        st.dataframe(annotated, width="stretch")
        st.download_button(
            "Descargar anotaciones (CSV)",
            annotated.to_csv(index=False),
            file_name="ensembl_anotado.csv",
            width="stretch",
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
    st.markdown(
        "Identificamos primero qué pozos pertenecen a controles y a muestras, "
        "tal como se explica en la sección de métodos de muchos artículos de qPCR: "
        "sin una separación clara entre ambos grupos, cualquier métrica posterior "
        "es ambigua. Usa este bloque como una sección de *Materiales y métodos* "
        "interactiva: selecciona la estrategia que mejor refleje tu diseño experimental "
        "y valida en la tabla previa que las etiquetas se vean consistentes."
    )
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
        st.dataframe(imputation.summary.head(50), width="stretch")
        st.download_button(
            "Descargar imputaciones (CSV)",
            imputation.summary.to_csv(index=False),
            file_name="detalle_imputaciones.csv",
            width="stretch",
        )

    st.divider()
    st.subheader("Calidad de datos")
    st.markdown(
        "Antes de avanzar a cálculos reproducibles, documentamos la calidad del panel. "
        "Las métricas que ves aquí juegan el mismo papel que los controles de calidad en "
        "la sección de resultados de un paper: verifican cobertura de genes, presencia de "
        "valores faltantes y estabilidad básica. Si detectas un desequilibrio, vale la pena "
        "volver sobre la clasificación o incluso replantear los filtros iniciales."
    )
    metrics = compute_quality_metrics(controles_df, muestras_df)
    _quality_summary(metrics)
    if not metrics.common_targets:
        st.error("No hay genes en común entre controles y muestras. Revisa la clasificación.")
        return

    st.divider()
    st.subheader("Normalización avanzada")
    st.markdown(
        "Esta etapa replica la justificación metodológica de la normalización avanzada. "
        "Al ajustar los parámetros estás, en la práctica, explorando distintas hipótesis "
        "sobre qué subconjuntos de genes de referencia mejor capturan la estabilidad del experimento. "
        "Los mensajes de retroalimentación resumen qué combinación quedó activa, del mismo modo "
        "que un artículo detalla en texto los criterios de optimización adoptados."
    )
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

    if st.button("Ejecutar normalización avanzada", key=f"btn_run_adv_{adv_key}", width="stretch"):
        with st.status("Procesando normalización avanzada…", expanded=True) as status:
            status.write(
                f"FDR={params.alpha:.3f}, candidatos={params.n_candidates}, K={params.k_refs}, "
                f"bootstrap={params.bootstrap_iter}, permutaciones={params.permutation_iter}"
            )
            try:
                result = execute_advanced_normalization(controles_df, muestras_df, params=params)
            except AdvancedNormalizationError as exc:
                status.update(label=f"Error durante la normalización: {exc}", state="error")
            else:
                adv_state["result"] = result
                adv_state["params"] = params
                status.update(label="Normalización avanzada completada", state="complete")

    result = adv_state.get("result")
    if not isinstance(result, AdvancedNormalizationResult):
        st.info("Configura los parámetros y ejecuta la normalización avanzada para continuar con el pipeline.")
        st.stop()

    params_used = adv_state.get("params", params)
    _notify_parameter_change(
        f"adv_params::{adv_key}",
        (
            params_used.alpha,
            params_used.n_candidates,
            params_used.k_refs,
            params_used.bootstrap_iter,
            params_used.permutation_iter,
        ),
        f"Resultados vigentes con α={params_used.alpha:.3f}, candidatos={params_used.n_candidates} y K={params_used.k_refs}.",
    )

    expr_summary = _render_advanced_results(result, params_used)
    st.session_state["_expression_summary"] = expr_summary

    st.divider()
    st.subheader("Expresión relativa comparada (3 métodos)")
    st.markdown(
        "Con la normalización avanzada como base, contrastamos los tres enfoques disponibles. "
        "Piensa esta sección como el panel principal de figuras: cada heatmap es el equivalente "
        "a una figura de resultados que sustenta la discusión. Elegir si ver solo los genes "
        "significativos o el panel completo permite alternar entre una narrativa centrada en hallazgos "
        "y otra que incluye la totalidad del contexto experimental."
    )
    genes_all = sorted(result.df_norm["target"].dropna().astype(str).unique().tolist())
    if not genes_all:
        st.warning("No se encontraron genes para construir el comparativo de expresión.")
    else:
        col_opts1, col_opts2 = st.columns([2, 1])
        with col_opts1:
            genes_mode = st.radio(
                "Genes para heatmap",
                options=["Significativos/Top-20", "Panel completo (todos)"],
                horizontal=True,
                help="Elige si visualizar únicamente los genes destacados o todo el panel disponible.",
                key=f"heatmap_mode::{adv_key}",
            )
        with col_opts2:
            zscore_rows = st.checkbox(
                "Estandarizar por gen (z-score)",
                value=True,
                help="Normaliza cada fila para resaltar patrones relativos entre tests.",
                key=f"heatmap_zscore::{adv_key}",
            )

        _notify_parameter_change(
            f"heatmap_comparativo::{adv_key}",
            (genes_mode, zscore_rows),
            "Heatmaps recalculados con los parámetros seleccionados.",
        )

        genes_for_heatmap = genes_all if genes_mode.endswith("(todos)") else None
        matrices, basic_ref_gene = build_method_matrices(
            controles_df,
            muestras_df,
            result,
            genes_for_heatmap=genes_for_heatmap,
        )
        refs_adv = ", ".join(result.reference_result.references)
        st.caption(
            f"Referencias avanzadas: {refs_adv or 'sin definir'}. Gen de referencia básico: {basic_ref_gene}."
        )

        tab_adv, tab_gm, tab_ref = st.tabs(["Avanzada", "Promedio global", f"Gen ref básico ({basic_ref_gene})"])
        with tab_adv:
            fig_a = build_dendrogram_heatmap(
                matrices.get("advanced", pd.DataFrame()),
                title="Avanzada — z-score por gen" if zscore_rows else "Avanzada",
                zscore_by_gene=zscore_rows,
            )
            st.plotly_chart(fig_a, width="stretch")
            st.download_button(
                "Descargar matriz avanzada (CSV)",
                matrices.get("advanced", pd.DataFrame()).to_csv(),
                file_name="heatmap_avanzada.csv",
                width="stretch",
            )
        with tab_gm:
            fig_g = build_dendrogram_heatmap(
                matrices.get("global_mean", pd.DataFrame()),
                title="Promedio global — z-score por gen" if zscore_rows else "Promedio global",
                zscore_by_gene=zscore_rows,
            )
            st.plotly_chart(fig_g, width="stretch")
            st.download_button(
                "Descargar matriz promedio global (CSV)",
                matrices.get("global_mean", pd.DataFrame()).to_csv(),
                file_name="heatmap_promedio_global.csv",
                width="stretch",
            )
        with tab_ref:
            fig_r = build_dendrogram_heatmap(
                matrices.get("refgene", pd.DataFrame()),
                title=f"Gen de referencia ({basic_ref_gene}) — z-score por gen" if zscore_rows else f"Gen de referencia ({basic_ref_gene})",
                zscore_by_gene=zscore_rows,
            )
            st.plotly_chart(fig_r, width="stretch")
            st.download_button(
                "Descargar matriz gen de referencia (CSV)",
                matrices.get("refgene", pd.DataFrame()).to_csv(),
                file_name="heatmap_gen_referencia.csv",
                width="stretch",
            )

    st.divider()
    st.subheader("Datasets de niveles de expresión por método")
    st.markdown(
        "Aquí dejamos listos los conjuntos de datos que alimentarán los análisis posteriores. "
        "Seleccionar un método equivale a explicitar, como harías en un paper, cuál será el marco "
        "de referencia para la interpretación de resultados downstream. Cada descarga ofrece tanto "
        "el formato largo (útil para estadísticas) como la matriz genes×tests (propicia para visualizaciones)."
    )
    expression_datasets, basic_ref_gene_expr = build_expression_datasets(
        controles_df, muestras_df, result
    )
    method_labels = {
        "advanced": "Avanzada (referencias seleccionadas)",
        "global_mean": "Promedio global por test",
        "refgene": f"Gen de referencia básico ({basic_ref_gene_expr})",
    }
    default_method = st.session_state.get("selected_expression_method", "advanced")
    selected_method = st.radio(
        "Selecciona el dataset principal para las siguientes secciones",
        options=list(method_labels.keys()),
        index=list(method_labels.keys()).index(default_method) if default_method in method_labels else 0,
        format_func=lambda key: method_labels.get(key, key),
        key=f"expr_dataset::{adv_key}",
    )
    _notify_parameter_change(
        f"expr_dataset_choice::{adv_key}",
        selected_method,
        f"Trabajando con el dataset '{method_labels.get(selected_method, selected_method)}'.",
    )
    st.session_state["selected_expression_method"] = selected_method

    tabs_expr = st.tabs([method_labels[m] for m in method_labels])
    for method_key, tab in zip(method_labels.keys(), tabs_expr):
        with tab:
            df_method = expression_datasets.get(method_key, pd.DataFrame()).copy()
            if df_method.empty:
                st.warning("No hay datos disponibles para este método.")
                continue
            st.caption(
                f"Filas: {len(df_method)}, genes únicos: {df_method['target'].nunique()}, tests: {df_method['test'].nunique()}"
            )
            st.dataframe(df_method.head(100), width="stretch")
            pivot = (
                df_method.pivot_table(index="target", columns="test", values="log2_rel_expr", aggfunc="mean")
                .sort_index()
            )
            col_long, col_pivot = st.columns(2)
            with col_long:
                st.download_button(
                    "Descargar formato largo (CSV)",
                    df_method.to_csv(index=False),
                    file_name=f"expresion_{method_key}_largo.csv",
                    width="stretch",
                )
            with col_pivot:
                st.download_button(
                    "Descargar matriz genes×tests (CSV)",
                    pivot.to_csv(),
                    file_name=f"expresion_{method_key}_matriz.csv",
                    width="stretch",
                )

    st.divider()
    st.subheader("Genes diferencialmente expresados y heatmaps dedicados")
    st.markdown(
        "Esta sección funciona como un *supplementary table* interactivo: reportamos qué genes "
        "separan de manera consistente a los grupos según cada método. El mensaje de éxito resume "
        "los conteos, y la lista unificada opera como la base de la narrativa biológica que seguirá. "
        "Los heatmaps dedicados permiten inspeccionar visualmente la separación que justificará "
        "las conclusiones en etapas de anotación y enriquecimiento."
    )
    try:
        col_params1, col_params2, col_params3 = st.columns([1, 1, 1])
        with col_params1:
            alpha_sel = st.number_input(
                "Alpha FDR (q)",
                min_value=0.001,
                max_value=0.25,
                value=0.05,
                step=0.005,
                format="%.3f",
                key=f"genes_alpha::{adv_key}",
            )
        with col_params2:
            topn_fb = st.number_input(
                "Top-N fallback si no hay DE",
                min_value=5,
                max_value=100,
                value=20,
                step=5,
                key=f"genes_topn::{adv_key}",
            )
        with col_params3:
            zscore_rows_bymethod = st.checkbox(
                "Estandarizar por gen (z-score)",
                value=True,
                key=f"genes_zscore::{adv_key}",
            )

        _notify_parameter_change(
            f"genes_de_params::{adv_key}",
            (alpha_sel, topn_fb, zscore_rows_bymethod),
            "Recalculando genes diferencialmente expresados con los parámetros indicados… listo.",
        )

        matrices_by_method, gene_lists, basic_ref_gene = build_heatmaps_by_method(
            controles_df,
            muestras_df,
            result,
            alpha=float(alpha_sel),
            top_n_fallback=int(topn_fb),
        )

        counts_msg = (
            f"Avanzada: {len(gene_lists.get('advanced', []))} | "
            f"Promedio: {len(gene_lists.get('global_mean', []))} | "
            f"Gen ref ({basic_ref_gene}): {len(gene_lists.get('refgene', []))}"
        )
        st.success(f"Genes diferencialmente expresados detectados → {counts_msg}")

        col_dl1, col_dl2, col_dl3, col_save = st.columns([1, 1, 1, 1])
        with col_dl1:
            st.download_button(
                "Descargar genes DE (avanzada)",
                pd.DataFrame({"gene": gene_lists.get("advanced", [])}).to_csv(index=False),
                file_name="genes_significativos_avanzada.csv",
                width="stretch",
            )
        with col_dl2:
            st.download_button(
                "Descargar genes DE (promedio)",
                pd.DataFrame({"gene": gene_lists.get("global_mean", [])}).to_csv(index=False),
                file_name="genes_significativos_promedio.csv",
                width="stretch",
            )
        with col_dl3:
            st.download_button(
                f"Descargar genes DE (gen ref: {basic_ref_gene})",
                pd.DataFrame({"gene": gene_lists.get("refgene", [])}).to_csv(index=False),
                file_name=f"genes_significativos_refgene_{basic_ref_gene}.csv",
                width="stretch",
            )
        with col_save:
            do_save = st.checkbox("Guardar conjuntos en resultados/", value=False, key=f"genes_save::{adv_key}")
            if do_save:
                paths = save_gene_sets(gene_lists, basic_ref_gene=str(basic_ref_gene))
                st.caption("Guardado: " + ", ".join(f"{k}: {v}" for k, v in paths.items()))

        union_genes = sorted(set().union(*gene_lists.values())) if gene_lists else []
        st.session_state["genes_heatmap_union"] = union_genes
        with st.expander("Lista completa de genes a priorizar en etapas posteriores"):
            union_df = pd.DataFrame({"gene": union_genes})
            st.dataframe(union_df, width="stretch", hide_index=True)
            st.download_button(
                "Descargar lista unificada (CSV)",
                union_df.to_csv(index=False),
                file_name="genes_diferenciales_unificados.csv",
                width="stretch",
            )

        venn_fig = build_three_way_venn(
            gene_lists,
            labels=("Avanzada", "Promedio", f"Gen ref ({basic_ref_gene})"),
        )
        st.plotly_chart(venn_fig, width="stretch")
        st.caption(
            "El diagrama de Venn resume cuántos genes son exclusivos o compartidos por los métodos; "
            "puedes priorizar aquellos en la intersección triple para discusiones posteriores."
        )

        st.write("Heatmaps por método con sus propios genes DE")
        tab_a2, tab_g2, tab_r2 = st.tabs(["Avanzada (DE)", "Promedio (DE)", f"Gen ref (DE: {basic_ref_gene})"])
        with tab_a2:
            fig_a2 = build_dendrogram_heatmap(
                matrices_by_method.get("advanced", pd.DataFrame()),
                title="Avanzada — genes DE (z-score)" if zscore_rows_bymethod else "Avanzada — genes DE",
                zscore_by_gene=zscore_rows_bymethod,
            )
            st.plotly_chart(fig_a2, width="stretch")
        with tab_g2:
            fig_g2 = build_dendrogram_heatmap(
                matrices_by_method.get("global_mean", pd.DataFrame()),
                title="Promedio global — genes DE (z-score)" if zscore_rows_bymethod else "Promedio global — genes DE",
                zscore_by_gene=zscore_rows_bymethod,
            )
            st.plotly_chart(fig_g2, width="stretch")
        with tab_r2:
            fig_r2 = build_dendrogram_heatmap(
                matrices_by_method.get("refgene", pd.DataFrame()),
                title=f"Gen de referencia ({basic_ref_gene}) — genes DE (z-score)" if zscore_rows_bymethod else f"Gen de referencia ({basic_ref_gene}) — genes DE",
                zscore_by_gene=zscore_rows_bymethod,
            )
            st.plotly_chart(fig_r2, width="stretch")
    except Exception as exc:
        st.info(f"No se pudo calcular la selección por método: {exc}")

    st.divider()
    st.subheader("Tablas finales: ΔΔCt y Fold Change (3 métodos integrados)")
    st.markdown(
        "Finalmente, consolidamos las métricas clave en tablas comparables con las que se publican "
        "en la sección de resultados cuantitativos de los artículos. Separar ΔΔCt y Fold Change "
        "ayuda a distinguir la interpretación en escala logarítmica y lineal, respectivamente. "
        "El bloque de resumen cuantitativo aporta, en formato textual, los mismos indicadores "
        "(correlaciones, discrepancias) que suelen describirse en la discusión para sustentar "
        "la robustez de los hallazgos."
    )
    try:
        fold_result = compute_fold_change_with_expression(controles_df, muestras_df)
    except ValueError as exc:
        st.error(f"No se pudo calcular el fold change: {exc}")
        return

    fc_core = fold_result.computation
    ddct_table, fc_table = build_combined_ddct_fc_tables(result, fc_core)
    ddct_table = ddct_table.sort_values("target").reset_index(drop=True)
    fc_table = fc_table.sort_values("target").reset_index(drop=True)

    st.markdown("**Tabla ΔΔCt (tres métodos)**")
    st.dataframe(ddct_table, width="stretch")
    st.download_button(
        "Descargar tabla ΔΔCt (CSV)",
        ddct_table.to_csv(index=False),
        file_name="tabla_delta_delta_ct_3_metodos.csv",
        width="stretch",
    )

    ddct_help = {
        "target": "Nombre del gen evaluado.",
        "delta_delta_ct_advanced": "ΔΔCt calculado con el conjunto de referencias avanzadas (negativo → sobreexpresión en muestras).",
        "delta_delta_ct_promedio": "ΔΔCt usando el método de promedios globales por grupo.",
        "delta_delta_ct_gen_ref": "ΔΔCt empleando el gen de referencia básico seleccionado automáticamente.",
    }
    with st.expander("Guía de columnas ΔΔCt"):
        for col, desc in ddct_help.items():
            st.markdown(f"- `{col}`: {desc}")

    st.markdown("**Tabla Fold Change (tres métodos)**")
    fc_table = fc_table.assign(
        max_abs_log2fc=fc_table[["log2fc_advanced", "log2fc_promedio", "log2fc_gen_ref"]].abs().max(axis=1)
    )
    st.dataframe(fc_table, width="stretch")
    st.download_button(
        "Descargar tabla Fold Change (CSV)",
        fc_table.to_csv(index=False),
        file_name="tabla_fold_change_3_metodos.csv",
        width="stretch",
    )

    fc_help = {
        "target": "Nombre del gen evaluado.",
        "fold_change_advanced": "Fold change lineal derivado de la normalización avanzada.",
        "fold_change_promedio": "Fold change lineal usando promedios globales.",
        "fold_change_gen_ref": "Fold change lineal usando el gen de referencia básico.",
        "log2fc_advanced": "log2FC correspondiente al método avanzado (positivo → sobreexpresión).",
        "log2fc_promedio": "log2FC obtenido con el método de promedios globales.",
        "log2fc_gen_ref": "log2FC basado en el gen de referencia básico.",
        "max_abs_log2fc": "Mayor |log2FC| observado entre los tres métodos (útil para ordenamiento).",
    }
    with st.expander("Guía de columnas Fold Change"):
        for col, desc in fc_help.items():
            st.markdown(f"- `{col}`: {desc}")

    classification_table = build_expression_classification_table(fc_table)
    st.markdown("**Clasificación por niveles de expresión (3 métodos)**")
    st.markdown(
        "Cada gen se etiqueta como subexpresado, estable o sobreexpresado según `<1`, `1–2` o `≥2` en fold change, "
        "repitiendo el criterio para los tres enfoques. Las columnas `delta_log2fc` muestran cómo difieren "
        "los log2FC entre métodos; valores altos anticipan discrepancias que merecen discusión."
    )
    st.dataframe(classification_table.head(150), width="stretch")
    st.download_button(
        "Descargar clasificaciones (CSV)",
        classification_table.to_csv(index=False),
        file_name="clasificacion_niveles_expresion_3_metodos.csv",
        width="stretch",
    )
    class_chart = build_classification_summary_chart(classification_table)
    st.plotly_chart(class_chart, width="stretch")
    st.caption(
        "La gráfica resume la distribución de niveles de expresión por método, facilitando comparar "
        "si algún enfoque sesga hacia la sobreexpresión o subexpresión." 
    )

    total_genes = int(fc_table["target"].notna().sum())
    if total_genes > 0:
        max_slider = min(150, total_genes)
        min_slider = min(10, max_slider)
        top_for_chart = st.slider(
            "Top de genes por |log2FC| para la gráfica comparativa",
            min_value=1 if max_slider == 1 else min_slider,
            max_value=max_slider,
            value=min_slider,
            step=1,
            key=f"chart_top::{adv_key}",
        )
        _notify_parameter_change(
            f"chart_top_feedback::{adv_key}",
            top_for_chart,
            f"Gráfica actualizada con los {top_for_chart} genes más extremos.",
        )
        chart_df = fc_table.sort_values("max_abs_log2fc", ascending=False).head(int(top_for_chart))
        chart_fig = build_fc_methods_scatter(chart_df)
        st.plotly_chart(chart_fig, width="stretch")
        st.caption(
            "Color → log2FC del método de gen de referencia. Tamaño → discrepancia frente al método avanzado."
        )

        summary = summarize_fc_methods(fc_table)
        with st.expander("Resumen cuantitativo de concordancia"):
            pair_labels = {
                ("advanced", "promedio"): "Avanzada vs Promedio",
                ("advanced", "gen_ref"): "Avanzada vs Gen ref",
                ("promedio", "gen_ref"): "Promedio vs Gen ref",
            }
            for pair, metrics in summary.get("pairwise_metrics", {}).items():
                label = pair_labels.get(pair, f"{pair[0]} vs {pair[1]}")
                st.write(
                    f"{label}: r={metrics['pearson_r']:.3f} | RMSE={metrics['rmse']:.3f} | MAE={metrics['mae']:.3f} (N={metrics['n']})"
                )
            for label, counts in summary.get("counts", {}).items():
                st.write(
                    f"Umbral {label}: Avanzada={counts['advanced']}, Promedio={counts['promedio']}, Gen ref={counts['gen_ref']}, Todos={counts['todos']}, Alguno={counts['alguno']}"
                )
            disc = summary.get("top_discrepancies", [])
            if disc:
                disc_df = pd.DataFrame(disc)
                st.markdown("Genes con mayor discrepancia respecto al método avanzado")
                st.dataframe(disc_df, width="stretch", hide_index=True)
                st.download_button(
                    "Descargar discrepancias (CSV)",
                    disc_df.to_csv(index=False),
                    file_name="discrepancias_vs_avanzado.csv",
                    width="stretch",
                )

    st.session_state["fold_change_expression_table"] = classification_table

    st.divider()
    _render_ensembl_section(expr_summary)


if __name__ == "__main__":
    main()
