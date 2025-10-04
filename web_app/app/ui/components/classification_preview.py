"""Classification preview board built on top of `insight_board`.

This replaces the ad-hoc HTML preview in the classification section with a
reusable component, keeping the light palette consistent across the app.
"""

from __future__ import annotations

from dataclasses import dataclass
from html import escape
from typing import Iterable, List, Optional

import re
import streamlit as st

from .insight_board import Badge, Card, render_insight_board

__all__ = ["render_classification_preview"]


def _slugify_id(s: str) -> str:
    safe = re.sub(r"[^a-zA-Z0-9_-]+", "-", str(s))
    if safe and safe[0].isdigit():
        safe = "id-" + safe
    return safe


def render_classification_preview(
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
    """Render sticky classification preview using the generic insight board.

    Mirrors the previous layout: badges row, two cards (Controles/Muestras),
    footer with application summary, and a banner for collision status.
    """

    mode_label = ("Auto" if auto_apply else "Manual") + " · aplicación"
    mode_tone = "success" if auto_apply else "warning"
    coverage_pct = f"{coverage:.0%}" if total_tests else "0%"
    manual_total = int(manual_ctrl_count + manual_samp_count)

    # Banner logic: prefer explicit collisions; otherwise show success
    overlaps_msgs: List[str] = []
    if overlap_prefixes:
        overlaps_msgs.append(
            "Los prefijos se solapan entre grupos: "
            + escape(", ".join(sorted({str(p) for p in overlap_prefixes})))
        )
    if overlap_tests:
        # Only show a short preview list to avoid very long banners
        sample_list = sorted({str(t) for t in overlap_tests})[:8]
        suffix = " …" if len(set(map(str, overlap_tests))) > 8 else ""
        overlaps_msgs.append(
            "Algunas pruebas coinciden en ambos grupos: "
            + escape(", ".join(sample_list) + suffix)
        )
    if overlaps_msgs:
        banner_text = " · ".join(overlaps_msgs)
        banner_tone = "danger"
    else:
        banner_text = "Sin colisiones detectadas entre prefijos ni tests."
        banner_tone = "success"

    footer_text: Optional[str] = None
    if (applied_ctrl_count or applied_samp_count):
        footer_text = (
            f"Clasificación aplicada → {applied_ctrl_count} controles · {applied_samp_count} muestras."
        )

    key = f"classification-preview::{file_key}"
    container_id = f"insight-board-{_slugify_id(key)}"

    # Sticky wrapper CSS targeting the board container id
    st.markdown(
        f"""
        <style>
        #{container_id} {{ position: sticky; top: 5.8rem; z-index: 1; }}
        </style>
        """,
        unsafe_allow_html=True,
    )

    render_insight_board(
        badges=[
            Badge(label=mode_label, tone=mode_tone),
            Badge(label=f"{assigned_tests} / {total_tests} tests asignados", tone="muted"),
            Badge(label=f"Cobertura {coverage_pct}", tone=("success" if coverage >= 0.999 else "muted")),
            Badge(label=f"Prefijos detectados: {int(prefix_count)}", tone="muted"),
            Badge(label=f"Asignaciones manuales: {manual_total}", tone=("warning" if manual_total else "muted")),
        ],
        left=Card(title="Controles", count=len(ctrl_preview), chips=ctrl_preview[:6]),
        right=Card(title="Muestras", count=len(samp_preview), chips=samp_preview[:6]),
        footer_text=footer_text,
        banner_text=banner_text,
        banner_tone=banner_tone,
        variant="light",
        key=key,
    )

