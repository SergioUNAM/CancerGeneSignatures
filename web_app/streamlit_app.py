from __future__ import annotations

import io
import logging
import re
import sys
from dataclasses import asdict
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
    ExportRegistry,
    Highlight,
    build_step_sequence,
    render_highlight_pills,
    render_pipeline_progress,
    render_sidebar_progress,
    render_export_panel,
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


_ADVANCED_PRESETS = {
    "Rápido": AdvancedNormalizationParams(
        alpha=0.1,
        n_candidates=12,
        k_refs=1,
        bootstrap_iter=100,
        permutation_iter=60,
        random_seed=123,
    ),
    "Equilibrado": AdvancedNormalizationParams(
        alpha=0.05,
        n_candidates=20,
        k_refs=2,
        bootstrap_iter=200,
        permutation_iter=100,
        random_seed=123,
    ),
    "Exhaustivo": AdvancedNormalizationParams(
        alpha=0.03,
        n_candidates=28,
        k_refs=3,
        bootstrap_iter=400,
        permutation_iter=220,
        random_seed=123,
    ),
}

_ADVANCED_PRESET_NOTES = {
    "Rápido": "Optimiza tiempos de cálculo con un K simple (1) y menos permutaciones.",
    "Equilibrado": "Configura un balance entre sensibilidad y robustez, ideal para corridas estándar.",
    "Exhaustivo": "Amplía la búsqueda de referencias y las iteraciones para detectar efectos sutiles.",
}


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
    exports: ExportRegistry,
    export_prefix: str,
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
        exports.register_dataframe(
            key=f"{export_prefix}::advanced_ranking",
            section="Normalización avanzada",
            label="Ranking de combinaciones candidata",
            file_name="normalizacion_avanzada_candidatos.csv",
            dataframe=candidates,
            method_label="Avanzada",
            description="Score de combinaciones evaluadas durante la búsqueda de referencias.",
            parameters={
                "α": params.alpha,
                "K": params.k_refs,
                "Bootstrap": params.bootstrap_iter,
                "Permutaciones": params.permutation_iter,
            },
        )

    st.markdown("**Tabla normalizada (df_norm)**")
    st.dataframe(result.df_norm.head(100), width="stretch")
    exports.register_dataframe(
        key=f"{export_prefix}::advanced_df_norm",
        section="Normalización avanzada",
        label="Tabla normalizada (df_norm)",
        file_name="normalizacion_avanzada_df_norm.csv",
        dataframe=result.df_norm,
        method_label="Avanzada",
        description="Resultados de ΔCt ajustados con la combinación óptima de referencias.",
        parameters={
            "α": params.alpha,
            "Refs": ", ".join(ref_res.references),
        },
    )

    stats_df = result.differential.stats.copy()
    if not stats_df.empty:
        st.markdown("**Estadísticas diferenciales**")
        st.dataframe(stats_df.head(50), width="stretch")
        exports.register_dataframe(
            key=f"{export_prefix}::advanced_stats",
            section="Normalización avanzada",
            label="Estadísticas diferenciales",
            file_name="normalizacion_avanzada_stats.csv",
            dataframe=stats_df,
            method_label="Avanzada",
            description="Resultados completos de pruebas de significancia tras la normalización avanzada.",
            parameters={
                "α": params.alpha,
                "Refs": ", ".join(ref_res.references),
            },
        )
    else:
        st.info("No se obtuvieron estadísticas significativas con los parámetros actuales.")

    heat_df = result.df_heatmap
    if isinstance(heat_df, pd.DataFrame) and not heat_df.empty:
        try:
            fig_heat = build_dendrogram_heatmap(heat_df, title="Avanzada — Heatmap con dendrogramas")
            st.plotly_chart(
                fig_heat,
                width="stretch",
                key=f"classification_heatmap::{export_prefix}",
            )
        except Exception:
            st.dataframe(heat_df, width="stretch")

    # (comparativo 3 métodos se realiza en la vista principal, donde están controles/muestras)

    expr_summary = _build_expression_summary(result.differential, params.alpha)
    if expr_summary.empty:
        st.info("No se generó un resumen de expresión. Ajusta parámetros o revisa la clasificación.")
    else:
        st.markdown("**Resumen de expresión diferenciada**")
        st.dataframe(expr_summary.head(50), width="stretch")
        exports.register_dataframe(
            key=f"{export_prefix}::advanced_expr_summary",
            section="Normalización avanzada",
            label="Resumen de expresión diferenciada",
            file_name="normalizacion_avanzada_resumen_expresion.csv",
            dataframe=expr_summary,
            method_label="Avanzada",
            description="Síntesis de genes con ΔΔCt y fold change tras la normalización avanzada.",
            parameters={
                "α": params.alpha,
                "Top candidatos": params.n_candidates,
            },
        )
        counts = expr_summary["nivel_expresion"].value_counts()
        st.success(
            "Resultados: "
            f"{counts.get('subexpresado', 0)} subexpresados · "
            f"{counts.get('estable', 0)} estables · "
            f"{counts.get('sobreexpresado', 0)} sobreexpresados"
        )
    return expr_summary


def _prepare_heatmap_export(matrix: pd.DataFrame, *, zscore: bool) -> pd.DataFrame:
    if not isinstance(matrix, pd.DataFrame) or matrix.empty:
        return pd.DataFrame()
    export_df = matrix.copy()
    if zscore:
        arr = export_df.values.astype(float)
        means = np.nanmean(arr, axis=1, keepdims=True)
        stds = np.nanstd(arr, axis=1, keepdims=True)
        stds[stds == 0] = 1.0
        export_df = pd.DataFrame(arr - means, index=export_df.index, columns=export_df.columns)
        export_df = export_df.divide(stds, axis=0)
    return export_df


def _render_ensembl_section(
    expr_summary: pd.DataFrame,
    exports: ExportRegistry,
    export_prefix: str,
) -> None:
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
        exports.register_dataframe(
            key=f"{export_prefix}::ensembl",
            section="Anotación Ensembl",
            label="Listado anotado Ensembl",
            file_name="ensembl_anotado.csv",
            dataframe=annotated,
            method_label="Consulta API Ensembl",
            description="Tabla con IDs, descripciones y niveles de expresión para la lista priorizada de genes.",
            parameters={
                "Genes": tot,
                "Hilos": int(workers),
            },
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
    exports = ExportRegistry()

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

    st.subheader("Paso 1 · Carga de datos")
    st.caption(
        "Sube tu Excel qPCR, valida su estructura y procesa la tabla base antes de avanzar al siguiente paso."
    )

    load_feedback: Optional[str] = None
    load_error: Optional[str] = None
    load_success = False

    upload_cols = st.columns((2, 1))
    with upload_cols[0]:
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
    with upload_cols[1]:
        st.markdown(
            "**Consejos**\n"
            "- Prefiere archivos con cabeceras en filas 1–4 y pozos en la primera columna.\n"
            "- Si cambias la política de imputación deberás reprocesar el archivo para reflejar los ajustes."
        )

    current_und_policy = str(st.session_state.get(UND_POLICY_KEY, policy_default))
    current_und_value = float(st.session_state.get(UND_VALUE_KEY, value_default))

    process_disabled = uploaded is None
    run_btn = st.button("Procesar archivo", type="primary", disabled=process_disabled)

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

    with st.sidebar:
        st.markdown("### Estado del flujo")
        render_sidebar_progress(pipeline_steps)
        st.caption("Usa la vista principal para completar cada etapa y observarás aquí el progreso consolidado.")

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
        "La barra lateral funciona como monitor de progreso; la carga de archivos y las decisiones analíticas "
        "se gestionan dentro de los pasos guiados del cuerpo principal para ofrecer contexto inmediato."
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
        exports.register_dataframe(
            key=f"{file_key}::imputaciones",
            section="Preparación e imputación",
            label="Detalle de imputaciones",
            file_name="detalle_imputaciones.csv",
            dataframe=imputation.summary,
            description="Registro de Ct reemplazados al aplicar la política de valores indeterminados.",
            parameters={
                "Política": und_policy,
                "Valor fijo": und_value if und_policy == "value" else "N/A",
            },
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

    with st.expander("Configuración avanzada", expanded=True):
        stored_params = adv_state.get("params")
        if "form" not in adv_state:
            if isinstance(stored_params, AdvancedNormalizationParams):
                adv_state["form"] = asdict(stored_params)
            else:
                adv_state["form"] = asdict(_ADVANCED_PRESETS["Equilibrado"])

        preset_options = [*list(_ADVANCED_PRESETS.keys()), "Personalizado"]
        preset_key = f"adv_preset_choice_{adv_key}"
        current_preset = adv_state.get("preset", "Equilibrado")
        if current_preset not in preset_options:
            current_preset = "Personalizado"

        selected_preset = st.radio(
            "Estilo de configuración",
            options=preset_options,
            horizontal=True,
            index=preset_options.index(current_preset),
            key=preset_key,
            help="Aplica valores sugeridos para equilibrar rapidez, sensibilidad o exhaustividad.",
        )

        if selected_preset != adv_state.get("preset"):
            adv_state["preset"] = selected_preset
            if selected_preset in _ADVANCED_PRESETS:
                adv_state["form"] = asdict(_ADVANCED_PRESETS[selected_preset])

        form_values = adv_state.get("form", asdict(_ADVANCED_PRESETS["Equilibrado"]))
        note = _ADVANCED_PRESET_NOTES.get(selected_preset)
        if note:
            st.caption(note)

        st.divider()
        left_col, right_col = st.columns(2, gap="large")

        with left_col:
            st.markdown("##### Sensibilidad y descubrimiento")
            alpha = st.slider(
                "Umbral FDR (q)",
                min_value=0.001,
                max_value=0.5,
                value=float(form_values.get("alpha", 0.05)),
                step=0.001,
                help="Controla la tolerancia a falsos positivos tras corrección FDR.",
            )
            n_candidates = st.slider(
                "Genes candidatos",
                min_value=4,
                max_value=32,
                value=int(form_values.get("n_candidates", 20)),
                help="Cantidad de genes evaluados como posibles referencias.",
            )
            k_refs_value = int(form_values.get("k_refs", 2))
            k_refs_options = [1, 2, 3, 4, 5]
            k_refs_index = k_refs_options.index(k_refs_value) if k_refs_value in k_refs_options else 1
            k_refs = st.selectbox(
                "Tamaño del set de referencia (K)",
                options=k_refs_options,
                index=k_refs_index,
                help="Número de genes promedio usados como referencia agregada.",
            )

        with right_col:
            st.markdown("##### Robustez y reproducibilidad")
            bootstrap_iter = st.slider(
                "Iteraciones bootstrap",
                min_value=0,
                max_value=600,
                value=int(form_values.get("bootstrap_iter", 200)),
                step=20,
                help="Remuestreos internos para estabilizar la selección de referencias.",
            )
            permutation_iter = st.slider(
                "Permutaciones",
                min_value=0,
                max_value=400,
                value=int(form_values.get("permutation_iter", 100)),
                step=20,
                help="Intercambios Control/Muestra para estimar tasa de falsos positivos.",
            )
            random_seed = st.number_input(
                "Semilla aleatoria",
                value=int(form_values.get("random_seed", 123) or 123),
                step=1,
                help="Fija el generador pseudoaleatorio y permite reproducir resultados.",
            )

        n_candidates = max(int(n_candidates), int(k_refs))
        st.caption("La app asegura que el número de candidatos siempre sea ≥ K para construir combinaciones válidas.")
        adv_state["form"] = {
            "alpha": float(alpha),
            "n_candidates": n_candidates,
            "k_refs": int(k_refs),
            "bootstrap_iter": int(bootstrap_iter),
            "permutation_iter": int(permutation_iter),
            "random_seed": int(random_seed),
        }

        if selected_preset in _ADVANCED_PRESETS:
            preset_dict = asdict(_ADVANCED_PRESETS[selected_preset])
            if any(adv_state["form"][key] != preset_dict[key] for key in preset_dict):
                adv_state["preset"] = "Personalizado"
                st.session_state[preset_key] = "Personalizado"
                selected_preset = "Personalizado"
                st.caption("Has ajustado los valores: el preset se marcó como Personalizado.")

        with st.expander("¿Qué controla cada parámetro?", expanded=False):
            st.markdown(
                """
                - **Umbral FDR (q)**: reduce falsos positivos al precio de omitir señales débiles si es muy estricto.
                - **Genes candidatos**: explora más combinaciones para encontrar referencias estables.
                - **Tamaño del set (K)**: promedio de genes usados como ancla; K≥2 amortigua sesgos individuales.
                - **Bootstrap**: refina la frecuencia con la que aparece cada gen en la referencia.
                - **Permutaciones**: estima cuántos genes surgirían por azar al intercambiar etiquetas.
                - **Semilla**: asegura reproducibilidad cuando se comparten resultados con otras personas.
                """
            )

    params = AdvancedNormalizationParams(
        alpha=float(adv_state["form"]["alpha"]),
        n_candidates=int(adv_state["form"]["n_candidates"]),
        k_refs=int(adv_state["form"]["k_refs"]),
        bootstrap_iter=int(adv_state["form"]["bootstrap_iter"]),
        permutation_iter=int(adv_state["form"]["permutation_iter"]),
        random_seed=int(adv_state["form"]["random_seed"]),
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

    expr_summary = _render_advanced_results(result, params_used, exports, adv_key)
    st.session_state["_expression_summary"] = expr_summary

    st.divider()
    st.info(
        "Las descargas a partir de esta etapa se concentran al final en el panel de exportación consolidado."
    )
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

        heatmap_labels = {
            "advanced": "Matriz avanzada para heatmap",
            "global_mean": "Matriz promedio global para heatmap",
            "refgene": f"Matriz gen de referencia ({basic_ref_gene})",
        }
        for method_key, label in heatmap_labels.items():
            export_df = _prepare_heatmap_export(
                matrices.get(method_key, pd.DataFrame()),
                zscore=zscore_rows,
            )
            if not export_df.empty:
                exports.register_dataframe(
                    key=f"{adv_key}::heatmap::{method_key}",
                    section="Visualizaciones comparativas",
                    label=label,
                    file_name=f"heatmap_{method_key}.csv",
                    dataframe=export_df,
                    method_label={
                        "advanced": "Avanzada",
                        "global_mean": "Promedio global",
                        "refgene": f"Gen ref ({basic_ref_gene})",
                    }.get(method_key, method_key),
                    description="Matriz utilizada para construir el heatmap con dendrogramas.",
                    parameters={
                        "Genes": export_df.shape[0],
                        "Tests": export_df.shape[1],
                        "z-score": "Sí" if zscore_rows else "No",
                        "Selección": genes_mode,
                    },
                )

        tab_adv, tab_gm, tab_ref = st.tabs(["Avanzada", "Promedio global", f"Gen ref básico ({basic_ref_gene})"])
        with tab_adv:
            fig_a = build_dendrogram_heatmap(
                matrices.get("advanced", pd.DataFrame()),
                title="Avanzada — z-score por gen" if zscore_rows else "Avanzada",
                zscore_by_gene=zscore_rows,
            )
            st.plotly_chart(
                fig_a,
                width="stretch",
                key=f"heatmap_adv::{adv_key}::{genes_mode}::{zscore_rows}",
            )
        with tab_gm:
            fig_g = build_dendrogram_heatmap(
                matrices.get("global_mean", pd.DataFrame()),
                title="Promedio global — z-score por gen" if zscore_rows else "Promedio global",
                zscore_by_gene=zscore_rows,
            )
            st.plotly_chart(
                fig_g,
                width="stretch",
                key=f"heatmap_global::{adv_key}::{genes_mode}::{zscore_rows}",
            )
        with tab_ref:
            fig_r = build_dendrogram_heatmap(
                matrices.get("refgene", pd.DataFrame()),
                title=f"Gen de referencia ({basic_ref_gene}) — z-score por gen" if zscore_rows else f"Gen de referencia ({basic_ref_gene})",
                zscore_by_gene=zscore_rows,
            )
            st.plotly_chart(
                fig_r,
                width="stretch",
                key=f"heatmap_ref::{adv_key}::{genes_mode}::{zscore_rows}",
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
                exports.register_dataframe(
                    key=f"{adv_key}::expression::{method_key}::long",
                    section="Datasets de expresión",
                    label=f"Dataset en formato largo ({method_labels[method_key]})",
                    file_name=f"expresion_{method_key}_largo.csv",
                    dataframe=df_method,
                    method_label=method_labels[method_key],
                    description="Tabla long-form con log2(expr) relativa para cada test y gen.",
                    parameters={
                        "Genes": df_method["target"].nunique(),
                        "Tests": df_method["test"].nunique(),
                    },
                )
                st.caption("Disponible en el panel de exportación consolidado.")
            with col_pivot:
                exports.register_dataframe(
                    key=f"{adv_key}::expression::{method_key}::matrix",
                    section="Datasets de expresión",
                    label=f"Matriz genes×tests ({method_labels[method_key]})",
                    file_name=f"expresion_{method_key}_matriz.csv",
                    dataframe=pivot,
                    method_label=method_labels[method_key],
                    description="Matriz pivoteada para análisis estadístico adicional.",
                    parameters={
                        "Genes": pivot.shape[0],
                        "Tests": pivot.shape[1],
                    },
                    include_index=True,
                )
                st.caption("Disponible en el panel de exportación consolidado.")

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
        export_gene_labels = {
            "advanced": ("Genes DE — método avanzado", "genes_significativos_avanzada.csv"),
            "global_mean": ("Genes DE — promedio global", "genes_significativos_promedio.csv"),
            "refgene": (
                f"Genes DE — gen de referencia ({basic_ref_gene})",
                f"genes_significativos_refgene_{basic_ref_gene}.csv",
            ),
        }
        for method_key, (label, filename) in export_gene_labels.items():
            genes = gene_lists.get(method_key, [])
            if genes:
                exports.register_dataframe(
                    key=f"{adv_key}::genes::{method_key}",
                    section="Genes diferencialmente expresados",
                    label=label,
                    file_name=filename,
                    dataframe=pd.DataFrame({"gene": genes}),
                    method_label={
                        "advanced": "Avanzada",
                        "global_mean": "Promedio global",
                        "refgene": f"Gen ref ({basic_ref_gene})",
                    }.get(method_key, method_key),
                    description="Lista priorizada según el método de normalización seleccionado.",
                    parameters={"α": alpha_sel, "Top fallback": topn_fb},
                )
        with col_dl1:
            st.metric("Genes DE (avanzada)", len(gene_lists.get("advanced", [])))
            st.caption("Descarga disponible en el panel consolidado.")
        with col_dl2:
            st.metric("Genes DE (promedio)", len(gene_lists.get("global_mean", [])))
            st.caption("Descarga disponible en el panel consolidado.")
        with col_dl3:
            st.metric(
                f"Genes DE (gen ref {basic_ref_gene})",
                len(gene_lists.get("refgene", [])),
            )
            st.caption("Descarga disponible en el panel consolidado.")
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
            if not union_df.empty:
                exports.register_dataframe(
                    key=f"{adv_key}::genes::union",
                    section="Genes diferencialmente expresados",
                    label="Lista unificada de genes DE",
                    file_name="genes_diferenciales_unificados.csv",
                    dataframe=union_df,
                    description="Unión de genes significativos detectados por los tres métodos.",
                    parameters={"α": alpha_sel, "Top fallback": topn_fb},
                )
            st.caption("Descarga disponible en el panel consolidado.")

        venn_fig = build_three_way_venn(
            gene_lists,
            labels=("Avanzada", "Promedio", f"Gen ref ({basic_ref_gene})"),
        )
        st.plotly_chart(
            venn_fig,
            width="stretch",
            key=f"venn_diagram::{adv_key}::{selected_method}",
        )
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
            st.plotly_chart(
                fig_a2,
                width="stretch",
                key=f"heatmap_focus_adv::{adv_key}::{selected_method}",
            )
        with tab_g2:
            fig_g2 = build_dendrogram_heatmap(
                matrices_by_method.get("global_mean", pd.DataFrame()),
                title="Promedio global — genes DE (z-score)" if zscore_rows_bymethod else "Promedio global — genes DE",
                zscore_by_gene=zscore_rows_bymethod,
            )
            st.plotly_chart(
                fig_g2,
                width="stretch",
                key=f"heatmap_focus_global::{adv_key}::{selected_method}",
            )
        with tab_r2:
            fig_r2 = build_dendrogram_heatmap(
                matrices_by_method.get("refgene", pd.DataFrame()),
                title=f"Gen de referencia ({basic_ref_gene}) — genes DE (z-score)" if zscore_rows_bymethod else f"Gen de referencia ({basic_ref_gene}) — genes DE",
                zscore_by_gene=zscore_rows_bymethod,
            )
            st.plotly_chart(
                fig_r2,
                width="stretch",
                key=f"heatmap_focus_ref::{adv_key}::{selected_method}",
            )
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
    total_ddct = int(ddct_table["target"].notna().sum())
    ddct_neg = int((ddct_table["delta_delta_ct_advanced"] < -0.25).sum())
    ddct_pos = int((ddct_table["delta_delta_ct_advanced"] > 0.25).sum())
    ddct_neutral = int(
        (ddct_table["delta_delta_ct_advanced"].abs() <= 0.25).sum()
    )
    render_highlight_pills(
        [
            Highlight(
                label="Genes evaluados",
                value=total_ddct,
                help="Total de genes presentes tras la normalización consolidada.",
            ),
            Highlight(
                label="ΔΔCt < 0 (sobreexpresados)",
                value=ddct_neg,
                help="Genes con ΔΔCt ≤ -0.25 en el método avanzado → sobreexpresión en muestras.",
            ),
            Highlight(
                label="ΔΔCt estable (±0.25)",
                value=ddct_neutral,
                help="Genes con variación marginal según el método avanzado (|ΔΔCt| ≤ 0.25).",
            ),
            Highlight(
                label="ΔΔCt > 0 (subexpresados)",
                value=ddct_pos,
                help="Genes con ΔΔCt ≥ 0.25 en el método avanzado → subexpresión en muestras.",
            ),
        ],
        key=f"ddct-overview::{adv_key}",
    )
    with st.expander("Ver tabla ΔΔCt completa"):
        st.dataframe(ddct_table, width="stretch")
    exports.register_dataframe(
        key=f"{adv_key}::tables::ddct",
        section="Tablas finales ΔΔCt/FC",
        label="Tabla ΔΔCt (3 métodos)",
        file_name="tabla_delta_delta_ct_3_metodos.csv",
        dataframe=ddct_table,
        description="ΔΔCt comparativo para los métodos avanzado, promedio y gen de referencia.",
        parameters={
            "α": alpha_sel,
            "Método expresión": method_labels.get(selected_method, selected_method),
        },
    )
    st.caption("Descarga disponible en el panel consolidado.")

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
    high_abs = int((fc_table["max_abs_log2fc"].abs() >= 1).sum())
    high_adv = int((fc_table["fold_change_advanced"].abs() >= 2).sum())
    under_adv = int((fc_table["fold_change_advanced"] <= 0.5).sum())
    render_highlight_pills(
        [
            Highlight(
                label="|log2FC| ≥ 1 (cualquier método)",
                value=high_abs,
                help="Genes que alcanzan un cambio de al menos 2× en alguno de los tres métodos.",
            ),
            Highlight(
                label="Sobreexpresión ≥ 2× (Avanzada)",
                value=high_adv,
                help="Genes con fold change ≥ 2 bajo el método avanzado.",
            ),
            Highlight(
                label="Subexpresión ≤ 0.5× (Avanzada)",
                value=under_adv,
                help="Genes cuyo fold change avanzado sugiere al menos 50% de reducción.",
            ),
        ],
        key=f"fc-overview::{adv_key}",
    )
    with st.expander("Ver tabla Fold Change completa"):
        st.dataframe(fc_table, width="stretch")
    exports.register_dataframe(
        key=f"{adv_key}::tables::fold_change",
        section="Tablas finales ΔΔCt/FC",
        label="Tabla Fold Change (3 métodos)",
        file_name="tabla_fold_change_3_metodos.csv",
        dataframe=fc_table,
        description="Fold change y log2FC comparables entre los tres enfoques de normalización.",
        parameters={
            "α": alpha_sel,
            "Top fallback": topn_fb,
        },
    )
    st.caption("Descarga disponible en el panel consolidado.")

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
    class_counts = classification_table["clasificacion_advanced"].value_counts().reindex(
        ["sobreexpresado", "estable", "subexpresado"],
        fill_value=0,
    )
    render_highlight_pills(
        [
            Highlight(
                label="Sobreexpresados (Avanzada)",
                value=int(class_counts["sobreexpresado"]),
                help="Genes clasificados como sobreexpresados en el método avanzado.",
            ),
            Highlight(
                label="Estables (Avanzada)",
                value=int(class_counts["estable"]),
                help="Genes con comportamiento estable (0.5×–2×) según el método avanzado.",
            ),
            Highlight(
                label="Subexpresados (Avanzada)",
                value=int(class_counts["subexpresado"]),
                help="Genes clasificados como subexpresados en el método avanzado.",
            ),
        ],
        key=f"classification-overview::{adv_key}",
    )
    with st.expander("Ver tabla de clasificación por métodos"):
        st.dataframe(classification_table, width="stretch")
    exports.register_dataframe(
        key=f"{adv_key}::tables::classification",
        section="Tablas finales ΔΔCt/FC",
        label="Clasificación de niveles de expresión",
        file_name="clasificacion_niveles_expresion_3_metodos.csv",
        dataframe=classification_table,
        description="Categorización de genes (subexpresado/estable/sobreexpresado) para los tres métodos.",
        parameters={
            "α": alpha_sel,
            "Top fallback": topn_fb,
        },
    )
    st.caption("Descarga disponible en el panel consolidado.")
    class_chart = build_classification_summary_chart(classification_table)
    st.plotly_chart(
        class_chart,
        width="stretch",
        key=f"classification_outcome::{selected_method}",
    )
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
    st.plotly_chart(
        volcano_fig,
        width="stretch",
        key=f"volcano_plot::{selected_method}",
    )
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
        st.plotly_chart(
            chart_fig,
            width="stretch",
            key=f"export_chart::{adv_key}::{selected_method}::{top_for_chart}",
        )
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
                exports.register_dataframe(
                    key=f"{adv_key}::tables::discrepancias",
                    section="Tablas finales ΔΔCt/FC",
                    label="Discrepancias vs método avanzado",
                    file_name="discrepancias_vs_avanzado.csv",
                    dataframe=disc_df,
                    description="Genes cuyo log2FC difiere más respecto al método avanzado.",
                    parameters={
                        "α": alpha_sel,
                        "Top fallback": topn_fb,
                    },
                )
                st.caption("Descarga disponible en el panel consolidado.")

    st.session_state["fold_change_expression_table"] = classification_table

    st.divider()
    _render_ensembl_section(expr_summary, exports, adv_key)

    export_context = {
        "Archivo": f"{df_loaded.source_name} · {df_loaded.sheet_name}" if df_loaded else "Sin archivo",
        "Tipo de cáncer": cancer_selected,
        "Contexto biológico": context_selected,
        "Política ND": f"{und_policy} ({und_value})" if und_policy == "value" else und_policy,
        "Preset normalización": adv_state.get("preset", "Equilibrado"),
        "α FDR": f"{params_used.alpha:.3f}",
        "K refs": params_used.k_refs,
        "Dataset principal": method_labels.get(selected_method, selected_method),
    }
    render_export_panel(exports, study_context=export_context)


if __name__ == "__main__":
    main()
