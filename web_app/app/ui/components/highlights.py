"""Small highlight pill components for summarizing results."""

from __future__ import annotations

from dataclasses import dataclass
from html import escape
from typing import Iterable, List, Optional

import streamlit as st

__all__ = ["Highlight", "render_highlight_pills"]


@dataclass
class Highlight:
    """Represents a pill metric with label, value and optional helper text."""

    label: str
    value: str | int | float
    help: Optional[str] = None

    def render_value(self) -> str:
        if isinstance(self.value, float):
            return f"{self.value:,.2f}".rstrip("0").rstrip(".")
        return f"{self.value:,}" if isinstance(self.value, int) else str(self.value)


def render_highlight_pills(highlights: Iterable[Highlight], *, key: str) -> None:
    """Render a responsive pill container with the given highlights."""

    items: List[Highlight] = [h for h in highlights if h is not None]
    if not items:
        return

    container_id = f"pill-container-{key}"
    pills_html = "".join(
        f"""
        <div class='pill' title='{escape(h.help or '')}'>
            <span class='value'>{escape(h.render_value())}</span>
            <span class='label'>{escape(h.label)}</span>
        </div>
        """
        for h in items
    )

    st.markdown(
        f"""
        <style>
            #{container_id} {{
                display: grid;
                gap: 0.6rem;
                grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
                margin: 0.8rem 0 1.1rem 0;
            }}
            #{container_id} .pill {{
                background: linear-gradient(135deg, rgba(30, 64, 175, 0.85), rgba(14, 116, 144, 0.82));
                border-radius: 14px;
                padding: 0.75rem 1rem;
                color: #f8fafc;
                box-shadow: 0 12px 30px rgba(15, 23, 42, 0.25);
                border: 1px solid rgba(148, 163, 184, 0.22);
                display: flex;
                flex-direction: column;
                gap: 0.3rem;
                min-height: 86px;
            }}
            #{container_id} .pill .value {{
                font-size: 1.4rem;
                font-weight: 700;
                letter-spacing: -0.03em;
            }}
            #{container_id} .pill .label {{
                font-size: 0.85rem;
                text-transform: uppercase;
                letter-spacing: 0.08em;
                color: rgba(241, 245, 249, 0.9);
            }}
            #{container_id} .pill:hover {{
                transform: translateY(-2px);
                transition: transform 120ms ease;
            }}
        </style>
        <div id='{container_id}'>
            {pills_html}
        </div>
        """,
        unsafe_allow_html=True,
    )
