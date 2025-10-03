"""Generic card components for summaries that adapt to light/dark themes.

Components
----------
- InfoCard: dataclass for a simple titled card with text or chips.
- render_info_cards(cards, key, columns): renders a responsive grid of cards.

Design
------
We compute colors from Streamlit theme (`theme.base`) and expose a
rounded, soft-contrast style that works in light and dark.
"""

from __future__ import annotations

from dataclasses import dataclass
from html import escape
from typing import Iterable, List, Optional, Sequence

import re
import streamlit as st

__all__ = ["InfoCard", "render_info_cards"]


@dataclass
class InfoCard:
    title: str
    text: Optional[str] = None
    chips: Optional[Sequence[str]] = None
    emphasis: bool = False  # slightly stronger background


def _slug(s: str) -> str:
    s = re.sub(r"[^a-zA-Z0-9_-]+", "-", s)
    if s and s[0].isdigit():
        s = "id-" + s
    return s


def _theme_colors() -> dict[str, str]:
    base = (st.get_option("theme.base") or "light").lower()
    if base == "dark":
        # Dark style inspired by the classification board screenshot
        return dict(
            bg=(
                "radial-gradient(1200px 400px at 0% 0%, rgba(30, 58, 138, 0.10), transparent), "
                "linear-gradient(180deg, rgba(2,6,23,0.65), rgba(2,6,23,0.45))"
            ),
            card_bg="rgba(30, 41, 59, 0.55)",
            card_bg_emph="rgba(30, 41, 59, 0.70)",
            border="rgba(148, 163, 184, 0.22)",
            shadow="0 18px 40px rgba(2,6,23,0.45)",
            text="#e6edf3",
            title="#cbd5e1",
            chip_bg="linear-gradient(180deg, rgba(37, 99, 235, 0.18), rgba(37, 99, 235, 0.12))",
            chip_text="#e6edf3",
        )
    # light
    return dict(
        bg=(
            "radial-gradient(1200px 400px at 0% 0%, rgba(59, 130, 246, 0.10), transparent), "
            "linear-gradient(180deg, rgba(241,245,249,0.95), rgba(241,245,249,0.75))"
        ),
        card_bg="rgba(244, 248, 255, 0.96)",
        card_bg_emph="rgba(235, 244, 255, 0.98)",
        border="rgba(15, 23, 42, 0.08)",
        shadow="0 18px 40px rgba(15, 23, 42, 0.10)",
        text="#0f172a",
        title="#0b1220",
        chip_bg="linear-gradient(180deg, rgba(37, 99, 235, 0.12), rgba(37, 99, 235, 0.08))",
        chip_text="#0b1220",
    )


def render_info_cards(cards: Iterable[InfoCard], *, key: str, columns: int = 1) -> None:
    items: List[InfoCard] = [c for c in cards if c is not None]
    if not items:
        return
    cols = max(1, min(4, int(columns)))
    colors = _theme_colors()
    cid = f"cards-{_slug(str(key))}"

    def _card_html(c: InfoCard) -> str:
        chips_html = ""
        if c.chips:
            chips_html = "<div class='chips'>" + "".join(
                f"<span class='chip'>{escape(str(x))}</span>" for x in c.chips
            ) + "</div>"
        content_html = (
            f"<div class='text'>{escape(c.text)}</div>" if (c.text and not c.chips) else chips_html
        )
        bg = colors["card_bg_emph"] if c.emphasis else colors["card_bg"]
        return f"""
        <div class='card' style='background:{bg}; box-shadow:{colors['shadow']}; border-color:{colors['border']}; color:{colors['text']}'>
          <div class='title' style='color:{colors['title']}'>{escape(c.title)}</div>
          {content_html}
        </div>
        """

    html = "".join(_card_html(c) for c in items)
    st.markdown(
        f"""
        <style>
        #{cid} {{
            background:{colors['bg']}; padding: 12px 12px; border-radius: 22px; border:1px solid {colors['border']};
            box-shadow:{colors['shadow']};
        }}
        #{cid} .grid {{
            display:grid; gap: 12px; grid-template-columns: repeat({cols}, minmax(260px, 1fr));
        }}
        #{cid} .card {{
            border:1px solid; border-radius: 18px; padding: 14px 16px;
        }}
        #{cid} .title {{
            font-weight: 800; letter-spacing: .035em; text-transform: none; margin-bottom: 8px;
        }}
        #{cid} .text {{
            line-height: 1.45; font-size: 0.95rem;
        }}
        #{cid} .chips {{ display:flex; flex-wrap:wrap; gap:8px; }}
        #{cid} .chip {{
            display:inline-block; padding: 6px 10px; border-radius: 999px; background:{colors['chip_bg']}; color:{colors['chip_text']};
            border:1px solid {colors['border']}; font-size: 12px;
        }}
        @media (max-width: 960px) {{ #{cid} .grid {{ grid-template-columns: 1fr; }} }}
        </style>
        <div id='{cid}'><div class='grid'>{html}</div></div>
        """,
        unsafe_allow_html=True,
    )
