"""Componentes UI para clasificar pruebas qPCR en controles y muestras."""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import md5
from html import escape
from typing import Dict, Iterable, List, Tuple

import pandas as pd
import streamlit as st
from streamlit.components.v1 import html as st_html
from app.ui.components import render_classification_preview

from app.services.qpcr import (
    PrefixGrouping,
    classify_tests_by_prefixes,
    detect_collisions,
    group_tests_by_initial_prefix,
)

__all__ = [
    "ClassificationResult",
    "format_list_preview",
    "store_classification",
    "clear_classification_state",
    "render_classification_section",
]


@dataclass
class ClassificationResult:
    """Agrupa los resultados de clasificación junto con metadatos de resumen."""

    controles: pd.DataFrame
    muestras: pd.DataFrame
    total_tests: int

    @property
    def coverage(self) -> float:
        """Porcentaje de pruebas asignadas a algún grupo."""

        if self.total_tests == 0:
            return 0.0
        assigned = len(
            set(self.controles["test"].astype(str)).union(self.muestras["test"].astype(str))
        )
        return assigned / float(self.total_tests)


def format_list_preview(items: Iterable[str], limit: int = 20) -> str:
    """Devuelve una representación amigable de una lista para tooltips o captions."""

    sequence = [str(item) for item in items if str(item)]
    if not sequence:
        return "(sin coincidencias)"
    head = sequence[:limit]
    suffix = " …" if len(sequence) > limit else ""
    return ", ".join(head) + suffix


def store_classification(
    state: Dict[str, object],
    ctrl_df: pd.DataFrame,
    samp_df: pd.DataFrame,
) -> ClassificationResult:
    """Guarda en el estado y resuelve colisiones entre controles y muestras."""

    ctrl_df = ctrl_df.copy()
    samp_df = samp_df.copy()
    collisions = detect_collisions(ctrl_df, samp_df)
    if collisions:
        preview = format_list_preview(collisions, limit=10)
        st.warning(
            "Algunas pruebas aparecen tanto en controles como en muestras. "
            "Se priorizarán como controles: "
            f"{preview}"
        )
        samp_df = samp_df[~samp_df["test"].astype(str).isin(collisions)].copy()

    state["controles_df"] = ctrl_df
    state["muestras_df"] = samp_df
    total_tests = len(pd.concat([ctrl_df["test"], samp_df["test"]]).astype(str).unique())
    return ClassificationResult(ctrl_df, samp_df, total_tests=total_tests)


def clear_classification_state(state: Dict[str, object]) -> None:
    """Elimina entradas relacionadas con clasificación del estado persistente."""

    for key in (
        "ctrl_prefix",
        "samp_prefix",
        "ctrl_prefixes",
        "samp_prefixes",
        "ctrl_unassigned",
        "samp_unassigned",
        "ctrl_suffix",
        "samp_suffix",
        "ctrl_regex",
        "samp_regex",
        "selected_ctrl",
        "selected_samp",
        "controles_df",
        "muestras_df",
        "_last_selection",
    ):
        state.pop(key, None)

    # Limpiar entradas asociadas en session_state de Streamlit
    for ui_key in list(st.session_state.keys()):
        if any(
            ui_key.startswith(prefix)
        for prefix in (
            "ctrl_prefix_input::",
            "samp_prefix_input::",
            "ctrl_suffix_input::",
            "samp_suffix_input::",
            "pref_ctrl_suggest::",
            "pref_samp_suggest::",
            "suff_ctrl_suggest::",
            "suff_samp_suggest::",
            "ctrl_prefixes::",
            "samp_prefixes::",
            "ctrl_unassigned::",
            "samp_unassigned::",
            "auto_apply::",
        )
        ):
            st.session_state.pop(ui_key, None)


def _resolve_manual_assignments(
    long_df: pd.DataFrame,
    ctrl_manual: Iterable[str],
    samp_manual: Iterable[str],
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    ctrl_df = long_df[long_df["test"].astype(str).isin(ctrl_manual)]
    samp_df = long_df[long_df["test"].astype(str).isin(samp_manual)]
    return ctrl_df, samp_df


def _render_sticky_preview(
    *,
    file_key: str,
    total_tests: int,
    ctrl_preview: List[str],
    samp_preview: List[str],
    assigned_tests: int,
    coverage: float,
    overlap_prefixes: Iterable[str],
    overlap_tests: Iterable[str],
    auto_apply: bool,
    manual_ctrl_count: int,
    manual_samp_count: int,
    applied_ctrl_count: int,
    applied_samp_count: int,
    prefix_count: int,
) -> None:
    """Renderiza el panel de vista previa sticky con el componente reusable."""

    render_classification_preview(
        file_key=file_key,
        total_tests=total_tests,
        ctrl_preview=ctrl_preview,
        samp_preview=samp_preview,
        assigned_tests=assigned_tests,
        coverage=coverage,
        overlap_prefixes=overlap_prefixes,
        overlap_tests=overlap_tests,
        auto_apply=auto_apply,
        manual_ctrl_count=manual_ctrl_count,
        manual_samp_count=manual_samp_count,
        applied_ctrl_count=applied_ctrl_count,
        applied_samp_count=applied_samp_count,
        prefix_count=prefix_count,
    )


def render_classification_section(long_df: pd.DataFrame, file_key: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Renderiza el bloque de clasificación de la UI y devuelve los dataframes resultantes."""

    state = st.session_state.setdefault(file_key, {})

    grouping: PrefixGrouping = group_tests_by_initial_prefix(long_df)
    prefix_groups = grouping.groups
    without_prefix = grouping.without_prefix
    prefix_options = sorted(prefix_groups.keys())
    label_map = {
        prefix: f"{prefix} · {len(prefix_groups[prefix])} tests"
        for prefix in prefix_options
    }

    unique_tests = sorted(long_df["test"].astype(str).dropna().unique().tolist())
    total_tests = len(unique_tests)

    if not prefix_options:
        st.warning("No se detectaron prefijos alfanuméricos en los nombres de prueba.")

    default_ctrl = [p for p in state.get("ctrl_prefixes", []) if p in prefix_options]
    default_samp = [p for p in state.get("samp_prefixes", []) if p in prefix_options]

    ctrl_key = f"ctrl_prefixes::{file_key}"
    samp_key = f"samp_prefixes::{file_key}"

    if ctrl_key in st.session_state:
        st.session_state[ctrl_key] = [p for p in st.session_state[ctrl_key] if p in prefix_options]
    else:
        st.session_state[ctrl_key] = default_ctrl

    if samp_key in st.session_state:
        st.session_state[samp_key] = [p for p in st.session_state[samp_key] if p in prefix_options]
    else:
        st.session_state[samp_key] = default_samp

    ctrl_unassigned_key = f"ctrl_unassigned::{file_key}"
    samp_unassigned_key = f"samp_unassigned::{file_key}"

    default_ctrl_unassigned = [t for t in state.get("ctrl_unassigned", []) if t in without_prefix]
    default_samp_unassigned = [t for t in state.get("samp_unassigned", []) if t in without_prefix]

    if ctrl_unassigned_key in st.session_state:
        st.session_state[ctrl_unassigned_key] = [
            t for t in st.session_state[ctrl_unassigned_key] if t in without_prefix
        ]
    else:
        st.session_state[ctrl_unassigned_key] = default_ctrl_unassigned

    if samp_unassigned_key in st.session_state:
        st.session_state[samp_unassigned_key] = [
            t for t in st.session_state[samp_unassigned_key] if t in without_prefix
        ]
    else:
        st.session_state[samp_unassigned_key] = default_samp_unassigned

    prefix_summary_rows = [
        {
            "Prefijo": prefix,
            "Tests": len(tests),
            "Ejemplos": format_list_preview(sorted(set(map(str, tests))), limit=5),
        }
        for prefix, tests in prefix_groups.items()
    ]
    prefix_summary = pd.DataFrame(prefix_summary_rows, columns=["Prefijo", "Tests", "Ejemplos"]).sort_values(
        "Prefijo"
    ) if prefix_summary_rows else pd.DataFrame(columns=["Prefijo", "Tests", "Ejemplos"])

    main_col, preview_col = st.columns([3, 2], gap="large")

    with main_col:
        st.markdown("#### 1. Explora los prefijos detectados")
        st.caption(
            "Revisa el catálogo de prefijos detectados con ejemplos representativos antes de asignarlos."
        )

        if not prefix_summary.empty:
            st.dataframe(
                prefix_summary,
                use_container_width=True,
                hide_index=True,
            )
        else:
            st.info("La tabla de prefijos aparecerá cuando se detecten patrones alfanuméricos.")

        st.markdown("#### 2. Selecciona reglas de clasificación")
        selector_cols = st.columns(2, gap="medium")
        with selector_cols[0]:
            ctrl_prefixes = st.multiselect(
                "Controles (prefijos)",
                options=prefix_options,
                default=default_ctrl,
                format_func=lambda value: label_map.get(value, value),
                key=ctrl_key,
            )
        with selector_cols[1]:
            samp_prefixes = st.multiselect(
                "Muestras (prefijos)",
                options=prefix_options,
                default=default_samp,
                format_func=lambda value: label_map.get(value, value),
                key=samp_key,
            )

        ctrl_unassigned: List[str] = []
        samp_unassigned: List[str] = []
        if without_prefix:
            st.markdown("#### 3. Ajustes manuales para tests sin prefijo")
            st.caption(
                "Asigne manualmente pozos sin prefijo detectado para conservarlos en el análisis."
            )
            manual_cols = st.columns(2, gap="medium")
            with manual_cols[0]:
                ctrl_unassigned = st.multiselect(
                    "Controles", options=without_prefix, default=default_ctrl_unassigned, key=ctrl_unassigned_key
                )
            with manual_cols[1]:
                samp_unassigned = st.multiselect(
                    "Muestras", options=without_prefix, default=default_samp_unassigned, key=samp_unassigned_key
                )
        else:
            st.caption("No se detectaron pozos fuera de las reglas de prefijos.")

        ctrl_preview = sorted(
            set(test for prefix in ctrl_prefixes for test in prefix_groups.get(prefix, []))
            .union(set(ctrl_unassigned))
        )
        samp_preview = sorted(
            set(test for prefix in samp_prefixes for test in prefix_groups.get(prefix, []))
            .union(set(samp_unassigned))
        )

        overlap_prefixes = set(ctrl_prefixes).intersection(samp_prefixes)
        overlap_tests = set(ctrl_preview).intersection(samp_preview)

        assigned_tests = len(set(ctrl_preview).union(samp_preview))
        coverage = (assigned_tests / total_tests) if total_tests else 0.0

        auto_apply_key = f"auto_apply::{file_key}"
        auto_apply_default = state.get("auto_apply", True)
        auto_apply_value = st.checkbox(
            "Aplicar automáticamente al cambiar la selección",
            value=st.session_state.get(auto_apply_key, auto_apply_default),
            key=auto_apply_key,
            help="Cuando está activo, la clasificación se actualiza en cuanto cambian los prefijos seleccionados.",
        )
        state["auto_apply"] = auto_apply_value

        selection_signature = (
            tuple(sorted(ctrl_prefixes)),
            tuple(sorted(samp_prefixes)),
            tuple(sorted(ctrl_unassigned)),
            tuple(sorted(samp_unassigned)),
        )

        if not auto_apply_value and selection_signature != tuple(state.get("_last_selection", ())) and (
            ctrl_preview or samp_preview
        ):
            st.info("Pulsa 'Aplicar clasificación' para actualizar la separación de controles y muestras.")

        def apply_classification() -> ClassificationResult:
            result = classify_tests_by_prefixes(long_df, ctrl_prefixes, samp_prefixes)
            ctrl_df = result.controles
            samp_df = result.muestras

            if ctrl_unassigned:
                extra_ctrl, _ = _resolve_manual_assignments(long_df, ctrl_unassigned, [])
                ctrl_df = (
                    pd.concat([ctrl_df, extra_ctrl], ignore_index=True)
                    .drop_duplicates()
                    .reset_index(drop=True)
                )
            if samp_unassigned:
                _, extra_samp = _resolve_manual_assignments(long_df, [], samp_unassigned)
                samp_df = (
                    pd.concat([samp_df, extra_samp], ignore_index=True)
                    .drop_duplicates()
                    .reset_index(drop=True)
                )

            classification = store_classification(state, ctrl_df, samp_df)
            state["ctrl_prefixes"] = ctrl_prefixes
            state["samp_prefixes"] = samp_prefixes
            state["ctrl_unassigned"] = ctrl_unassigned
            state["samp_unassigned"] = samp_unassigned
            state["_last_selection"] = selection_signature
            return classification

        classification_applied = False
        if auto_apply_value and (ctrl_preview or samp_preview):
            if selection_signature != tuple(state.get("_last_selection", ())):
                classification = apply_classification()
                st.success(
                    f"Clasificación aplicada → controles: {len(classification.controles)}, "
                    f"muestras: {len(classification.muestras)}"
                )
                classification_applied = True
        else:
            if st.button("Aplicar clasificación", key=f"btn_pref_{file_key}"):
                if not ctrl_preview or not samp_preview:
                    st.warning("Selecciona al menos un elemento en cada grupo para clasificar.")
                else:
                    classification = apply_classification()
                    st.success(
                        f"Clasificación aplicada → controles: {len(classification.controles)}, "
                        f"muestras: {len(classification.muestras)}"
                    )
                    classification_applied = True

        stored_ctrl_df = state.get("controles_df")
        stored_samp_df = state.get("muestras_df")
        if not isinstance(stored_ctrl_df, pd.DataFrame):
            stored_ctrl_df = pd.DataFrame()
        if not isinstance(stored_samp_df, pd.DataFrame):
            stored_samp_df = pd.DataFrame()

        if not stored_ctrl_df.empty or not stored_samp_df.empty:
            st.markdown("#### Vista previa de resultados aplicados")
            tabs = st.tabs(["Controles", "Muestras"])
            with tabs[0]:
                st.caption(f"{len(stored_ctrl_df)} filas")
                st.dataframe(stored_ctrl_df.head(50), use_container_width=True, height=260)
            with tabs[1]:
                st.caption(f"{len(stored_samp_df)} filas")
                st.dataframe(stored_samp_df.head(50), use_container_width=True, height=260)

        applied_ctrl_count = len(stored_ctrl_df)
        applied_samp_count = len(stored_samp_df)

    _render_sticky_preview(
        file_key=file_key,
        total_tests=total_tests,
        ctrl_preview=ctrl_preview,
        samp_preview=samp_preview,
        assigned_tests=assigned_tests,
        coverage=coverage,
        overlap_prefixes=overlap_prefixes,
        overlap_tests=overlap_tests,
        auto_apply=auto_apply_value,
        manual_ctrl_count=len(ctrl_unassigned),
        manual_samp_count=len(samp_unassigned),
        applied_ctrl_count=applied_ctrl_count,
        applied_samp_count=applied_samp_count,
        prefix_count=len(prefix_options),
    )

    if not classification_applied:
        ctrl_df = state.get("controles_df")
        samp_df = state.get("muestras_df")
    else:
        ctrl_df = state.get("controles_df")
        samp_df = state.get("muestras_df")

    if not isinstance(ctrl_df, pd.DataFrame):
        ctrl_df = pd.DataFrame()
    if not isinstance(samp_df, pd.DataFrame):
        samp_df = pd.DataFrame()

    return ctrl_df, samp_df
