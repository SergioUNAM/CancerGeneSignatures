from __future__ import annotations

import io
import logging
import re
import sys
from pathlib import Path
from typing import List, Optional, Sequence

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
    build_three_way_venn,
    build_volcano_plot,
    summarize_fc_methods,
)
from app.services.heatmap_visuals import build_dendrogram_heatmap
from app.services.qpcr import build_long_table, summarize_extraction
from app.state import AppSessionState
from app.ui.components import (
    build_step_sequence,
    render_pipeline_progress,
    render_sidebar_progress,
)
from app.ui.sections import render_classification_section


@st.cache_data(show_spinner=False, ttl=3600)
def _annotate_ensembl_cached(df_csv: str, max_workers: int = 3) -> pd.DataFrame:
    try:
        df = pd.read_csv(io.StringIO(df_csv))
    except Exception:
        return pd.DataFrame()
    if df.empty:
        return pd.DataFrame()
    return add_ensembl_info_batch(df, symbol_col="target", max_workers=max_workers)


_PIPELINE_STEPS_INFO = (
    (
        "Carga de datos",
        "Sube tu Excel qPCR, elige la hoja y valida la vista previa.",
    ),
    (
        "Clasificación",
        "Separa controles y muestras para habilitar métricas de calidad.",
    ),
    (
        "Normalización",
        "Ejecuta la búsqueda de referencias y genera resultados reproducibles.",
    ),
    (
        "Resultados y exportación",
        "Descarga heatmaps, tablas y listas de genes diferencialmente expresados.",
    ),
    (
        "Anotación Ensembl",
        "Completa IDs y descripciones para tus genes prioritarios.",
    ),
)


def _compute_pipeline_completion(
    df_loaded: Optional[LoadResult],
) -> Sequence[bool]:
    """Evalúa qué etapas del pipeline están completas en función del estado actual."""

    has_data = df_loaded is not None
    classification_done = False
    normalization_done = False
    exports_ready = False
    ensembl_ready = False

    if has_data and df_loaded is not None:
        file_key = f"assign::{df_loaded.source_name}:{df_loaded.sheet_name}"
        assign_state = st.session_state.get(file_key, {})
        ctrl_df = assign_state.get("controles_df")
        samp_df = assign_state.get("muestras_df")
        classification_done = (
            isinstance(ctrl_df, pd.DataFrame)
            and isinstance(samp_df, pd.DataFrame)
            and not ctrl_df.empty
            and not samp_df.empty
        )

        adv_key = f"adv::{file_key}"
        adv_state = st.session_state.get(adv_key, {})
        normalization_done = isinstance(
            adv_state.get("result"),
            AdvancedNormalizationResult,
        )
        exports_ready = normalization_done

    ensembl_df = st.session_state.get("_ensembl_df")
    if isinstance(ensembl_df, pd.DataFrame) and not ensembl_df.empty:
        ensembl_ready = True

    return has_data, classification_done, normalization_done, exports_ready, ensembl_ready


def _pipeline_status_labels(completions: Sequence[bool]) -> Sequence[str]:
    statuses: List[str] = []
    active_assigned = False
    for completed in completions:
        if completed:
            statuses.append("complete")
        elif not active_assigned:
            statuses.append("active")
            active_assigned = True
        else:
            statuses.append("pending")
    if not statuses:
        return statuses
    if not active_assigned:
        # Todas las etapas se marcaron como completas.
        return statuses
    return statuses


def _build_pipeline_steps(
    df_loaded: Optional[LoadResult],
) -> tuple[Sequence[object], Sequence[bool]]:
    """Genera los pasos del pipeline y sus estados de completitud actuales."""

    completions = _compute_pipeline_completion(df_loaded)
    statuses = _pipeline_status_labels(completions)
    step_defs = [
        (title, description, status)
        for (title, description), status in zip(_PIPELINE_STEPS_INFO, statuses)
    ]
    return build_step_sequence(step_defs), completions


def _eq(formula: str) -> None:
    """Render LaTeX inline helper."""
    st.latex(formula)


def _notify_parameter_change(state_key: str, value: object, message: str) -> None:
    """Emite un mensaje informativo cuando un parámetro interactivo cambia."""

    feedback_state = st.session_state.setdefault("_param_feedback", {})
    if feedback_state.get(state_key) != value:
        feedback_state[state_key] = value
        st.info(message)


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
    df_loaded = app_state.df_loaded

    # Definimos llaves compartidas para los selectores ubicados en el cuerpo principal.
    UND_POLICY_KEY = "_und_policy"
    UND_VALUE_KEY = "_und_value"
    STUDY_CANCER_KEY = "_study_cancer"
    STUDY_CONTEXT_KEY = "_study_context"

    und_defaults = {"nan", "ctmax", "value"}
    policy_default = (
        app_state.undetermined.policy
        if app_state.undetermined.policy in und_defaults
        else "nan"
    )
    value_default = float(app_state.undetermined.value)

    if UND_POLICY_KEY not in st.session_state:
        st.session_state[UND_POLICY_KEY] = policy_default
    elif st.session_state[UND_POLICY_KEY] not in und_defaults:
        st.session_state[UND_POLICY_KEY] = policy_default

    if UND_VALUE_KEY not in st.session_state:
        st.session_state[UND_VALUE_KEY] = value_default

    cancer_types = app_config.menu.cancer_types
    if not cancer_types:
        st.error("El menú de configuración no define tipos de cáncer disponibles.")
        cancer_types = ["Sin definir"]
    cancer_default = (
        app_state.cancer_type
        if app_state.cancer_type in cancer_types
        else cancer_types[0]
    )
    if STUDY_CANCER_KEY not in st.session_state:
        st.session_state[STUDY_CANCER_KEY] = cancer_default
    elif st.session_state[STUDY_CANCER_KEY] not in cancer_types:
        st.session_state[STUDY_CANCER_KEY] = cancer_default

    context_options = [ctx.label for ctx in app_config.menu.contexts]
    if not context_options:
        context_options = ["Sin contexto"]
    context_default = (
        app_state.context_label
        if app_state.context_label in context_options
        else context_options[0]
    )
    if STUDY_CONTEXT_KEY not in st.session_state:
        st.session_state[STUDY_CONTEXT_KEY] = context_default
    elif st.session_state[STUDY_CONTEXT_KEY] not in context_options:
        st.session_state[STUDY_CONTEXT_KEY] = context_default

    st.title("Flujo qPCR → Normalización avanzada → Anotación Ensembl")
    st.caption(
        "Sigue las etapas guiadas para pasar de Ct crudos a resultados exportables con anotación genética."
    )

    progress_placeholder = st.container()

    with st.expander("Guía rápida del flujo", expanded=False):
        st.markdown(
            "**Ruta rápida**\n"
            "1) Carga el archivo Excel y revisa pozos problemáticos u outliers.\n"
            "2) Clasifica controles y muestras (prefijos, sufijos, regex o selección manual).\n"
            "3) Elige la política para valores indeterminados y revisa las métricas de calidad.\n"
            "4) Ejecuta la **normalización avanzada** y explora heatmaps, diagrama de Venn y volcano plot.\n"
            "5) Exporta datasets, listas de genes y las anotaciones Ensembl."
        )

    with st.expander("¿Cómo se imputan los valores 'Undetermined/ND'?", expanded=False):
        st.markdown(
            "**Política de imputación**\n"
            "• `nan`: conserva Ct no determinados; útil para visualizar huecos reales.\n"
            "• `ctmax`: sustituye por el Ct máximo observado del grupo; aproxima el límite de detección.\n"
            "• `value`: usa un Ct fijo (p. ej. 40) para penalizar indetectables de forma uniforme."
        )

    pipeline_steps: Sequence[object] = []

    with st.sidebar:
        status_placeholder = st.container()
        st.divider()
        st.subheader("Paso 1 · Datos de entrada")
        st.caption("Sube tu Excel qPCR y confirma su estructura antes de avanzar.")
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
        process_disabled = uploaded is None
        run_btn = st.button("Procesar archivo", type="primary", disabled=process_disabled)

        has_data = df_loaded is not None

        load_feedback: Optional[str] = None
        load_error: Optional[str] = None
        load_success = False
        current_und_policy = str(st.session_state.get(UND_POLICY_KEY, policy_default))
        current_und_value = float(st.session_state.get(UND_VALUE_KEY, value_default))

        if uploaded is not None and run_btn:
            candidate: Optional[LoadResult] = None
            uploaded.seek(0)
            try:
                candidate = parse_qpcr_wide(
                    uploaded,
                    sheet_name=sheet,
                    header_mode="coords",
                    header_row_idx=3,
                    well_col_idx=0,
                    target_col_idx=1,
                    undetermined_policy=current_und_policy,
                    undetermined_value=current_und_value,
                )
            except Exception as exc:
                load_feedback = (
                    f"Cabecera esperada no encontrada ({exc}); intentando detección automática…"
                )
                uploaded.seek(0)
                try:
                    candidate = parse_qpcr_wide(
                        uploaded,
                        sheet_name=sheet,
                        header_mode="auto",
                        undetermined_policy=current_und_policy,
                        undetermined_value=current_und_value,
                    )
                except Exception as exc_auto:
                    load_error = f"No se pudo leer el archivo: {exc_auto}"
                    logger.exception("Fallo al parsear archivo", exc_info=exc_auto)
                    candidate = None
                else:
                    df_loaded = candidate
                    load_success = True
            else:
                df_loaded = candidate
                load_success = True

            if load_success and df_loaded is not None:
                app_state.df_loaded = df_loaded
                logger.info(
                    "Archivo cargado: %s | hoja=%s | forma=%s",
                    df_loaded.source_name,
                    df_loaded.sheet_name,
                    df_loaded.df.shape,
                )

        has_data = df_loaded is not None

        if load_error:
            st.error(load_error)

        if has_data and df_loaded is not None:
            st.success(
                "Archivo activo: "
                f"{df_loaded.source_name} · {df_loaded.df.shape[0]} filas × {df_loaded.df.shape[1]} columnas."
            )
            if load_feedback:
                st.caption(load_feedback)
        else:
            st.info(
                "Carga un archivo Excel de qPCR y pulsa 'Procesar archivo' para habilitar el resto de pasos."
            )

        pipeline_steps, _ = _build_pipeline_steps(df_loaded)
        with status_placeholder:
            st.markdown("### Estado del flujo")
            render_sidebar_progress(pipeline_steps)

    # Recuperamos los valores vigentes tras posibles interacciones en el cuerpo principal.
    und_policy = str(st.session_state.get(UND_POLICY_KEY, policy_default))
    und_value = float(st.session_state.get(UND_VALUE_KEY, value_default))
    cancer_selected = str(st.session_state.get(STUDY_CANCER_KEY, cancer_default))
    context_selected = str(st.session_state.get(STUDY_CONTEXT_KEY, context_default))

    app_state.undetermined.policy = und_policy
    app_state.undetermined.value = float(und_value)
    app_state.cancer_type = str(cancer_selected)
    app_state.context_label = str(context_selected)
    app_state.df_loaded = df_loaded
    app_state.persist()

    with progress_placeholder:
        render_pipeline_progress(pipeline_steps)

    st.markdown(
        "La barra lateral queda enfocada en la carga de datos y el estado general; las decisiones analíticas "
        "se gestionan ahora en el cuerpo principal para ofrecer más contexto al configurarlas."
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

    st.subheader("Política de imputación de Ct indeterminados")
    st.caption(
        "Define cómo se reemplazan los valores 'Undetermined/ND' al volver a procesar el archivo. "
        "Los cambios requieren pulsar nuevamente **Procesar archivo** para aplicarse a la tabla base."
    )
    imp_cols = st.columns((2, 1))
    with imp_cols[0]:
        st.selectbox(
            "Política",
            options=["nan", "ctmax", "value"],
            key=UND_POLICY_KEY,
            disabled=df_loaded is None,
        )
    with imp_cols[1]:
        st.number_input(
            "Valor fijo (si aplica)",
            min_value=0.0,
            max_value=100.0,
            step=0.5,
            key=UND_VALUE_KEY,
            disabled=df_loaded is None,
        )
    st.caption(
        "• `nan`: conserva el Ct sin determinar; útil para visualizar huecos reales.\n"
        "• `ctmax`: sustituye por el Ct máximo válido observado en el grupo.\n"
        "• `value`: aplica un Ct constante (p. ej. 40) para penalizar indetectables."
        )

    st.divider()

    st.subheader("Etiquetas del estudio")
    st.caption("Selecciona cómo se documentará el contexto biológico del análisis.")
    study_cols = st.columns(2)
    with study_cols[0]:
        st.selectbox(
            "Tipo de cáncer",
            options=cancer_types,
            key=STUDY_CANCER_KEY,
            disabled=df_loaded is None,
        )
    with study_cols[1]:
        st.selectbox(
            "Contexto biológico",
            options=context_options,
            key=STUDY_CONTEXT_KEY,
            disabled=df_loaded is None,
        )

    # Actualizamos las variables locales tras posibles cambios de los widgets.
    und_policy = str(st.session_state.get(UND_POLICY_KEY, policy_default))
    und_value = float(st.session_state.get(UND_VALUE_KEY, value_default))
    cancer_selected = str(st.session_state.get(STUDY_CANCER_KEY, cancer_default))
    context_selected = str(st.session_state.get(STUDY_CONTEXT_KEY, context_default))

    st.markdown("**Configuración actual del estudio**")
    col_conf1, col_conf2 = st.columns(2)
    with col_conf1:
        st.info(f"Tipo de cáncer: {cancer_selected}")
    with col_conf2:
        st.info(f"Contexto biológico: {context_selected}")

    st.divider()

    long_df, controls_warning = build_long_table(df_loaded)
    if controls_warning:
        st.warning(f"No se filtraron los controles de máquina: {controls_warning}")

    st.subheader("Clasificación de controles vs muestras")
    st.markdown(
        "El primer paso es definir con precisión **qué pozos corresponden a controles y cuáles a muestras**. "
        "Sin esta separación, las métricas posteriores carecen de validez. Este bloque funciona como una sección "
        "de *Materiales y métodos* interactiva: seleccione la estrategia que refleje su diseño experimental y confirme "
        "en la vista previa que las etiquetas sean coherentes."
    )
    st.caption(
        "Revisión práctica: asegúrese de que cada regla capture solo los pozos esperados. Si un test aparece en ambos grupos, "
        "ajuste la regla o utilice la selección manual."
    )
    file_key = f"assign::{df_loaded.source_name}:{df_loaded.sheet_name}"
    controles_df, muestras_df = render_classification_section(long_df, file_key)

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
        "Previo al análisis, se documenta la **calidad del panel**. Estas métricas cubren cobertura génica, presencia de valores faltantes "
        "y estabilidad básica. Si se observa desbalance, conviene reetiquetar o revisar los filtros iniciales." 
    )
    st.caption(
        "Como referencia práctica, valores de NaN superiores al 10% suelen justificar repetir mediciones o ajustar la estrategia de imputación."
    )
    metrics = compute_quality_metrics(controles_df, muestras_df)
    _quality_summary(metrics)
    if not metrics.common_targets:
        st.error("No hay genes en común entre controles y muestras. Revisa la clasificación.")
        return

    st.divider()
    st.subheader("Normalización avanzada")
    st.markdown(
        "En esta etapa se identifica la **referencia más estable** y se aplican las transformaciones clásicas de qPCR. "
        "Los parámetros ajustables representan hipótesis alternativas sobre la estabilidad de subconjuntos de genes de referencia."
    )
    _eq(r"S = 0.7\,\sigma_{\text{intra}} + 0.3\,\lvert \mu_{\text{Control}} - \mu_{\text{Muestra}} \rvert")
    _eq(r"\Delta Ct_i = Ct_i - Ct_{\text{ref}}")
    _eq(r"\log_2(\mathrm{expr}_i) = -\,\Delta Ct_i,\quad \Delta\Delta Ct = \Delta Ct_{\text{muestra}} - \Delta Ct_{\text{control}},\quad \mathrm{FC} = 2^{-\Delta\Delta Ct}")
    st.markdown(
        "El módulo automatiza la búsqueda de referencia y aplica ΔCt/ΔΔCt de manera explícita y reproducible."
    )
    st.markdown(
        """
        En la puntuación \(S\) se privilegia la estabilidad intra-grupo (70 %) y se deja un 30 % para penalizar desplazamientos
        sistemáticos entre controles y muestras. Esta elección sigue la recomendación usual en textos de diseño experimental, donde se sugiere ponderar con
        mayor peso la varianza interna cuando el objetivo es fijar un patrón de referencia robusto, manteniendo un término de
        sesgo para evitar que la referencia arrastre diferencias reales. Si la variabilidad intra-grupo es pequeña y la brecha
        entre medias es nula, el score se minimiza, lo que garantiza que los genes elegidos se comporten como verdaderos
        “housekeeping”.
        """
    )
    st.caption(
        "Cuando el algoritmo avanzado termina seleccionando el mismo gen (K = 1) o un promedio casi idéntico al del método básico, "
        "los perfiles normalizados coinciden de forma natural: ambos restan el mismo Ct y obtienen log2FC equivalentes. "
        "Las diferencias emergen cuando la referencia óptima requiere combinar varios genes o cuando los controles muestran un sesgo "
        "que el método de promedios no corrige; en esos casos, las tablas comparativas y el resumen cuantitativo revelan "
        "discrepancias en log2FC mientras el básico permanece más conservador o más ruidoso."
    )
    adv_key = f"adv::{file_key}"
    adv_state = st.session_state.setdefault(adv_key, {})

    with st.expander("Parámetros", expanded=True):
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
                options=[1, 2, 3, 4, 5],
                index=[1, 2, 3, 4, 5].index(int(getattr(adv_state.get("params"), "k_refs", 2)))
                if int(getattr(adv_state.get("params"), "k_refs", 2)) in {1, 2, 3, 4, 5}
                else 1,
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
        st.markdown(
            """
            **Interpretación de los parámetros**

            - *Umbral FDR (q)*: nivel de corrección por pruebas múltiples. Valores pequeños (≤0.05) privilegian especificidad; subirlo (0.1) amplía la lista de genes pero aumenta el riesgo de falsos positivos.
            - *Genes candidatos*: número de genes con menor variabilidad que se exploran como referencias. Un N reducido da rapidez pero puede omitir combinaciones útiles; valores altos permiten descubrir referencias alternativas a costa de mayor tiempo de cómputo.
            - *Tamaño set referencia (K)*: cantidad de genes que se promedian. K=1 replica el método básico; K≥2 mitiga ruido y sesgos. El sistema garantiza que el número de candidatos sea al menos igual a K; si la cohorte es pequeña, reduzca K para evitar combinaciones imposibles.
            - *Semilla aleatoria*: fija el generador pseudoaleatorio para reproducir exactamente bootstrap y permutaciones. Cambiarla puede alterar levemente las frecuencias empíricas.
            - *Iteraciones bootstrap*: número de re-muestreos del conjunto de pruebas. Más iteraciones (≥300) estabilizan la frecuencia de genes significativos; pocas iteraciones agilizan el cálculo pero dan estimaciones más ruidosas.
            - *Permutaciones*: barajadas Control/Muestra usadas para estimar la tasa de falsos positivos. Incrementarlas (≥200) refina la estimación del FPR; valores bajos aceleran la corrida sacrificando precisión.
            """
        )

    params = AdvancedNormalizationParams(
        alpha=float(alpha),
        n_candidates=max(int(n_candidates), int(k_refs)),
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
        "Con la referencia elegida se contrastan **tres enfoques**. Cada heatmap actúa como figura de resultados: "
        "es posible mostrar únicamente los genes significativos o todo el panel, según se busque una narrativa focalizada "
        "o un contexto más amplio."
    )
    st.caption(
        "Columnas homogéneas sugieren consistencia entre muestras, mientras que bandas horizontales intensas reflejan genes con cambios robustos."
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
        "Se generan los conjuntos de datos que alimentarán los análisis posteriores. Seleccionar un método equivale a fijar el **marco de referencia** "
        "para interpretar resultados subsecuentes. Cada pestaña ofrece formato largo y matriz genes×tests."
    )
    st.caption(
        "Es recomendable conservar los tres para justificar comparaciones y la elección metodológica en la discusión."
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
        "Se listan los genes que distinguen consistentemente a los grupos según cada método. "
        "La lista unificada constituye la base para la interpretación biológica y la búsqueda bibliográfica."
    )
    st.caption(
        "La intersección triple en el diagrama de Venn suele señalar candidatos robustos; los exclusivos pueden evidenciar efectos dependientes del método de normalización."
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

        matrices_by_method, gene_lists, stats_by_method, basic_ref_gene = build_heatmaps_by_method(
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
        "Se presentan tablas de **ΔΔCt** y **Fold Change**, equivalentes a las incluidas en artículos científicos. "
        "Separar ΔΔCt (escala logarítmica) y FC (escala lineal) facilita la interpretación."
    )
    st.caption(
        "Signos invertidos entre métodos en ΔΔCt sugieren revisar la referencia. En FC, se recomienda priorizar genes con |log2FC| ≥ 1 y q < α."
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

    volcano_fig = build_volcano_plot(
        stats_by_method,
        labels={
            "advanced": "Avanzada",
            "global_mean": "Promedio",
            "refgene": f"Gen ref ({basic_ref_gene})",
        },
        q_threshold=float(alpha_sel),
        log2fc_threshold=np.log2(2.0),
    )
    st.plotly_chart(volcano_fig, width="stretch")
    st.caption(
        "El gráfico volcano resume magnitud (log2FC) y significancia (-log10 p). Los genes resaltados superan los umbrales de α y diferencia mínima."
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
            st.markdown(
                """
                Este panel resume cuán alineados están los tres métodos cuando miramos los log2FC por gen:

                * **r, RMSE, MAE** cuantifican la similitud entre pares de métodos. r≈1 implica ordenamientos casi idénticos;
                  RMSE y MAE bajos indican magnitudes homogéneas. N especifica cuántos genes entraron en el cálculo.
                * **Umbrales ≥1.5x, ≥2x** muestran cuántos genes superan esos fold change mínimos en cada método, cuántos coinciden todos
                  y cuántos al menos uno. Útil para detectar métodos más conservadores o liberales.
                * **Genes con mayor discrepancia** lista los casos donde la diferencia en log2FC respecto al método avanzado es más alta;
                  revisa estos genes antes de sacar conclusiones.
                """
            )
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
