"""Reusable summary/insight board components for compact previews.

The goal is to provide a compact header with badges and two large cards
showing counts and representative items (as chips), plus an optional
footer line and banner message. It is intentionally lightweight and
CSS-driven so it can be reused in multiple sections (ΔΔCt/FC, etc.).
"""

from __future__ import annotations

from dataclasses import dataclass
import re
from html import escape
from typing import Iterable, List, Optional, Sequence

import streamlit as st

__all__ = [
    "Badge",
    "Card",
    "render_insight_board",
]


@dataclass
class Badge:
    label: str
    tone: str = "muted"  # one of: primary, success, warning, danger, muted


@dataclass
class Card:
    title: str
    count: int
    chips: Sequence[str]
    # Máximo de chips a mostrar. None = sin límite
    chips_max: Optional[int] = 8
    # Lista de chips a resaltar (p. ej., controles de máquina)
    accent_chips: Optional[Sequence[str]] = None
    accent_tone: str = "warning"


def _tone_color(tone: str) -> tuple[str, str]:
    """Default (dark) tone colors kept for backward compatibility."""

    tone = (tone or "").lower()
    if tone == "primary":
        return ("#1f6feb", "#ffffff")
    if tone == "success":
        return ("#176b47", "#e8fff4")
    if tone == "warning":
        return ("#7c4a03", "#fff6e5")
    if tone == "danger":
        return ("#e11d1d", "#ffffff")
    return ("#2a2f3a", "#e6edf3")


def render_insight_board(
    *,
    badges: Iterable[Badge] | None = None,
    left: Card,
    right: Optional[Card] = None,
    footer_text: Optional[str] = None,
    banner_text: Optional[str] = None,
    banner_tone: str = "success",
    variant: str = "dark",  # 'dark' (default, app-wide) or 'light' (grays)
    key: str,
) -> None:
    """Render a 2-card summary board with badges and optional banner.

    - `badges`: small rounded labels shown above the cards.
    - `left`/`right`: two cards, each with a title, a large count, and chips.
    - `footer_text`: short line rendered below the cards area.
    - `banner_text`: optional full-width banner with success/warn/danger tone.
    - `key`: unique identifier to scope CSS and DOM ids.
    """

    def _slugify_id(s: str) -> str:
        # Replace any non [a-zA-Z0-9_-] with dash to keep CSS selectors valid
        safe = re.sub(r"[^a-zA-Z0-9_-]+", "-", s)
        if safe and safe[0].isdigit():
            safe = "id-" + safe
        return safe

    # Variant first, to drive tone mappings below
    v = (variant or "dark").lower()

    def tone_color(name: str) -> tuple[str, str]:
        name = (name or "").lower()
        if v == "light":
            # Light palette aligned with the classification preview board
            if name == "primary":
                return ("rgba(37, 99, 235, 0.18)", "#0b1220")  # blue tint
            if name == "success":
                return ("rgba(16, 185, 129, 0.18)", "#065f46")
            if name == "warning":
                return ("rgba(245, 158, 11, 0.18)", "#7c2d12")
            if name == "danger":
                # keep vivid red for explicit alerts/visual tests
                return ("#e11d1d", "#ffffff")
            # muted
            return ("rgba(148, 163, 184, 0.25)", "#0f172a")
        # default dark
        return _tone_color(name)

    container_id = f"insight-board-{_slugify_id(str(key))}"
    badge_html = ""
    items: List[Badge] = list(badges or [])
    for b in items:
        bg, fg = tone_color(b.tone)
        badge_html += (
            f"<span class='badge' style='background:{bg}; color:{fg}'>{escape(b.label)}</span>"
        )

    def card_html(card: Card) -> str:
        items = list(card.chips or [])
        max_n = card.chips_max if card.chips_max is not None else len(items)
        show = items[:max(0, max_n)]
        accent_set = {str(x).upper() for x in (card.accent_chips or [])}
        # Accent style per chip using the tone color
        abg, afg = tone_color(card.accent_tone)
        chips = "".join(
            (
                f"<span class='chip' style='background:{abg}; color:{afg}; border-color:rgba(255,255,255,0.18)'>"
                f"{escape(str(c))}</span>"
                if str(c).upper() in accent_set
                else f"<span class='chip'>{escape(str(c))}</span>"
            )
            for c in show
        )
        return f"""
        <div class='card'>
          <div class='card-title'>{escape(card.title)}</div>
          <div class='card-count'>{int(card.count)}</div>
          <div class='chip-list'>{chips}</div>
        </div>
        """

    banner_html = ""
    if banner_text:
        bg, fg = tone_color(banner_tone)
        banner_html = f"<div class='banner' style='background:{bg}; color:{fg}'>{escape(banner_text)}</div>"

    if v == "light":
        # Light, grayish palette with high contrast accents
        bg_grad = (
            "radial-gradient(1200px 400px at 0% 0%, rgba(148, 163, 184, 0.15), transparent), "
            "linear-gradient(180deg, rgba(245,247,250,0.98), rgba(235,238,243,0.92))"
        )
        card_bg = "rgba(248, 250, 252, 0.98)"  # slate-50
        border = "rgba(15, 23, 42, 0.10)"
        text = "#0f172a"  # slate-900
        title = "#334155"  # slate-700
        chip_bg = (
            "linear-gradient(180deg, rgba(37, 99, 235, 0.12), rgba(37, 99, 235, 0.08))"
        )
        chip_border = "rgba(59, 130, 246, 0.28)"
        shadow = "0 18px 38px rgba(15, 23, 42, 0.12)"
        banner_border = "rgba(148, 163, 184, 0.35)"
    else:
        # Default dark palette used by the rest of boards
        bg_grad = (
            "radial-gradient(1200px 400px at 0% 0%, rgba(30, 58, 138, 0.10), transparent), "
            "linear-gradient(180deg, rgba(2, 6, 23, 0.65), rgba(2, 6, 23, 0.45))"
        )
        card_bg = "rgba(30, 41, 59, 0.55)"
        border = "rgba(148, 163, 184, 0.22)"
        text = "#e6edf3"
        title = "#cbd5e1"
        chip_bg = "rgba(51, 65, 85, 0.7)"
        chip_border = border
        shadow = "0 18px 40px rgba(2,6,23,0.45)"
        banner_border = "rgba(255,255,255,0.12)"

    cols_css = "repeat(2, minmax(280px, 1fr))" if right is not None else "1fr"

    st.markdown(
        f"""
        <style>
        #{container_id} {{
            border-radius: 16px;
            border: 1px solid {border};
            background: {bg_grad};
            padding: 14px 16px 10px 16px;
            margin: 10px 0 12px 0;
            box-shadow:{shadow};
        }}
        #{container_id} .badges {{
            display:flex; flex-wrap:wrap; gap:8px; align-items:center; margin-bottom: 12px;
        }}
        #{container_id} .badge {{
            display:inline-block; padding: 6px 10px; border-radius: 999px; font-size: 12px; font-weight: 600;
            border: 1px solid {banner_border};
        }}
        #{container_id} .cards {{
            display:grid; grid-template-columns: {cols_css}; gap: 14px; margin-bottom: 10px;
        }}
        #{container_id} .card {{
            background: {card_bg};
            border: 1px solid {border};
            border-radius: 14px;
            padding: 12px 14px;
            color: {text};
        }}
        #{container_id} .card-title {{
            font-weight: 700; letter-spacing: 0.04em; color: {title}; text-transform: uppercase; font-size: 13px;
        }}
        #{container_id} .card-count {{
            font-size: 34px; font-weight: 800; margin-top: 6px; margin-bottom: 10px; letter-spacing: -0.02em;
        }}
        #{container_id} .chip-list {{ display:flex; flex-wrap:wrap; gap:8px; }}
        #{container_id} .chip {{
            display:inline-block; background: {chip_bg}; color:{text}; border:1px solid {chip_border};
            padding:4px 10px; border-radius: 999px; font-size: 12px;
        }}
        #{container_id} .footer {{ color:{title}; margin: 2px 2px 8px 2px; font-size: 13px; opacity:0.9; }}
        #{container_id} .banner {{
            margin-top: 10px; padding: 10px 12px; border-radius: 12px; font-weight: 600; border:1px solid {banner_border};
        }}
        @media (max-width: 900px) {{
            #{container_id} .cards {{ grid-template-columns: 1fr; }}
        }}
        </style>
        <div id="{container_id}">
          <div class='badges'>{badge_html}</div>
          <div class='cards'>
            {card_html(left)}
            {card_html(right) if right is not None else ''}
          </div>
          {f"<div class='footer'>{escape(footer_text)}</div>" if footer_text else ""}
          {banner_html}
        </div>
        """,
        unsafe_allow_html=True,
    )
