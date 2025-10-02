"""Renderizado de un resumen visual del flujo principal."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence

import streamlit as st


@dataclass(frozen=True)
class PipelineStep:
    """Representa un paso del pipeline junto con su estado."""

    title: str
    description: str
    status: str  # ``complete`` | ``active`` | ``pending``


_STATUS_EMOJI = {
    "complete": "✅",
    "active": "🟡",
    "pending": "⚪️",
}


def _normalize_status(value: str) -> str:
    if value not in _STATUS_EMOJI:
        return "pending"
    return value


def render_pipeline_progress(steps: Sequence[PipelineStep]) -> None:
    """Dibuja los pasos del pipeline en columnas con un indicador de estado."""

    if not steps:
        return

    normalized = [PipelineStep(step.title, step.description, _normalize_status(step.status)) for step in steps]
    cols = st.columns(len(normalized))
    for idx, (col, step) in enumerate(zip(cols, normalized), start=1):
        emoji = _STATUS_EMOJI.get(step.status, _STATUS_EMOJI["pending"])
        with col:
            st.markdown(f"#### {emoji} {idx}. {step.title}")
            st.caption(step.description)


def build_step_sequence(step_defs: Iterable[tuple[str, str, str]]) -> Sequence[PipelineStep]:
    """Convierte una secuencia de tuplas en ``PipelineStep``s."""

    return [PipelineStep(title, description, status) for title, description, status in step_defs]
