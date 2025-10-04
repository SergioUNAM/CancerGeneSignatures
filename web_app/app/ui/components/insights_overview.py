"""Composable overview for rendering multiple insight boards.

Provides a flexible API to assemble and render the existing compact
insight boards (extraction, ΔΔCt, FC) in a stack or tabbed layout.

Usage sketch (inside a Streamlit page):

    specs = build_default_insight_specs(data)
    render_insights_overview(specs, layout="stack", key_prefix="overview")

Where `data` is an `InsightsData` instance carrying the minimal inputs
needed to compute each section's `BoardSpec`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from html import escape
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import math
import re
import pandas as pd
import numpy as np
import streamlit as st

from .insight_board import Badge, Card, render_insight_board

__all__ = [
    "BoardSpec",
    "InsightsData",
    "build_extraction_specs",
    "build_ddct_specs",
    "build_ddct_comparative_spec",
    "build_fc_specs",
    "build_fc_comparative_spec",
    "build_concordance_spec_from_summary",
    "build_default_insight_specs",
    "render_insights_overview",
]


# ---------- Data contracts ----------


@dataclass
class BoardSpec:
    """Minimal render contract for a single insight board."""

    key: str
    badges: Sequence[Badge]
    left: Card
    right: Optional[Card] = None
    footer_text: Optional[str] = None
    banner_text: Optional[str] = None
    banner_tone: str = "success"
    variant: str = "dark"  # or "light"
    tab: Optional[str] = None  # Optional title when using tabs


@dataclass
class InsightsData:
    """Inputs required to generate board specs for overview sections."""

    # Extraction
    source_name: str
    sheet_name: Optional[str]
    tests: Sequence[str] = field(default_factory=list)
    wells: Sequence[str] = field(default_factory=list)
    genes: Sequence[str] = field(default_factory=list)
    control_markers: Sequence[str] = field(default_factory=lambda: ("PPC", "RTC"))

    # Context key for uniqueness across the page
    context_key: str = "default"

    # ΔΔCt table with columns: target, delta_delta_ct_*
    ddct_table: Optional[pd.DataFrame] = None
    # FC table with columns: target, fold_change_*, log2fc_*
    fc_table: Optional[pd.DataFrame] = None

    # Column mappings for flexibility
    ddct_cols: Dict[str, Tuple[str, str]] = field(
        default_factory=lambda: {
            "advanced": ("Avanzada", "delta_delta_ct_advanced"),
            "promedio": ("Promedio", "delta_delta_ct_promedio"),
            "gen_ref": ("Gen ref", "delta_delta_ct_gen_ref"),
        }
    )
    fc_cols: Dict[str, Tuple[str, str, str]] = field(
        default_factory=lambda: {
            "advanced": ("Avanzada", "fold_change_advanced", "log2fc_advanced"),
            "promedio": ("Promedio", "fold_change_promedio", "log2fc_promedio"),
            "gen_ref": ("Gen ref", "fold_change_gen_ref", "log2fc_gen_ref"),
        }
    )

    # Thresholds
    ddct_neutral: float = 0.25
    fc_over: float = 2.0
    fc_under: float = 0.5
    log2fc_abs: float = 1.0


# ---------- Helpers ----------


def _slug(s: str) -> str:
    s = re.sub(r"[^a-zA-Z0-9_-]+", "-", str(s))
    if s and s[0].isdigit():
        s = "id-" + s
    return s


def _sgn(x: Any) -> int:
    try:
        fx = float(x)
    except Exception:
        return 0
    if math.isfinite(fx) is False:
        return 0
    if abs(fx) < 1e-12:
        return 0
    return 1 if fx > 0 else -1


# ---------- Builders: Extraction ----------


def build_extraction_specs(data: InsightsData) -> List[BoardSpec]:
    tests = list(map(str, data.tests or []))
    wells = list(map(str, data.wells or []))
    genes = list(map(str, data.genes or []))

    # Board 1: Tests vs Wells
    b1 = BoardSpec(
        tab="Extracción",
        badges=[
            Badge(label="Extracción", tone="primary"),
            Badge(label=f"Archivo: {data.source_name}", tone="muted"),
            Badge(label=f"Hoja: {data.sheet_name or '-'}", tone="muted"),
        ],
        left=Card(title="Pruebas detectadas", count=len(tests), chips=tests),
        right=Card(title="Pozos detectados", count=len(wells), chips=wells),
        footer_text=f"Listas detectadas → {len(tests)} pruebas · {len(wells)} pozos.",
        variant="light",
        key=f"extraction-board-main::{data.source_name}:{data.sheet_name}",
    )

    # Board 2: Targets (genes + machine controls)
    controls_set = {str(x).upper() for x in (data.control_markers or [])}
    machine_controls = [g for g in genes if g.upper() in controls_set]
    gene_targets = [g for g in genes if g.upper() not in controls_set]
    combined = gene_targets + machine_controls
    b2 = BoardSpec(
        tab="Targets",
        badges=[
            Badge(label="Targets (genes + controles)", tone="success"),
            Badge(label=f"Total: {len(combined)}", tone="muted"),
            Badge(label=f"Controles: {len(machine_controls)}", tone="warning"),
        ],
        left=Card(
            title="Listado de targets",
            count=len(combined),
            chips=combined,
            chips_max=None,
            accent_chips=machine_controls,
            accent_tone="danger",
        ),
        right=None,
        footer_text=f"Genes: {len(gene_targets)} · Controles: {len(machine_controls)}.",
        variant="light",
        key=f"extraction-board-targets::{data.source_name}:{data.sheet_name}",
    )

    return [b1, b2]


# ---------- Builders: ΔΔCt ----------


def build_ddct_specs(data: InsightsData) -> List[BoardSpec]:
    df = data.ddct_table
    if df is None or df.empty:
        return []
    if "target" not in df.columns:
        return []
    res: List[BoardSpec] = []
    for m_key, (m_label, col) in data.ddct_cols.items():
        if col not in df.columns:
            continue
        series = df[col]
        n_eval = int(series.notna().sum())
        coverage = float(series.notna().mean())
        neg = int((series < -abs(data.ddct_neutral)).sum())
        pos = int((series > abs(data.ddct_neutral)).sum())
        neutral = int((series.abs() <= abs(data.ddct_neutral)).sum())
        over_list = df.sort_values(col, ascending=True)["target"].head(6).tolist()
        under_list = df.sort_values(col, ascending=False)["target"].head(6).tolist()
        res.append(
            BoardSpec(
                tab=f"ΔΔCt · {m_label}",
                badges=[
                    Badge(label=f"ΔΔCt · {m_label}", tone="primary"),
                    Badge(label=f"Cobertura {coverage:.0%}", tone="success"),
                    Badge(label=f"Genes evaluados: {n_eval}", tone="muted"),
                ],
                left=Card(title="Sobreexpresados (ΔΔCt < 0)", count=neg, chips=over_list),
                right=Card(title="Subexpresados (ΔΔCt > 0)", count=pos, chips=under_list),
                footer_text=(
                    f"Clasificación aplicada → {neg} sobreexpresados · {neutral} estables · {pos} subexpresados."
                ),
                variant="light",
                key=f"ddct-board::{data.context_key}::{m_key}",
            )
        )
    return res


def build_ddct_comparative_spec(data: InsightsData) -> List[BoardSpec]:
    df = data.ddct_table
    if df is None or df.empty:
        return []
    cols = [t[1] for t in data.ddct_cols.values() if t[1] in df.columns]
    if len(cols) < 2:
        return []
    _dd = df[cols].copy()
    coverage = float(_dd.notna().all(axis=1).mean())
    sgn_cols = {c: _dd[c].apply(_sgn) for c in cols}
    base = sgn_cols[cols[0]]
    mismatch_mask = np.zeros(len(base), dtype=bool)
    for c in cols[1:]:
        mismatch_mask |= (sgn_cols[c] != base) & (base != 0)
    ddct_mismatch = int(mismatch_mask.sum())
    ddct_match = int((((pd.DataFrame(sgn_cols).eq(base, axis=0)).all(axis=1)) & (base != 0)).sum())

    _dd2 = df[["target", *cols]].copy()
    _dd2 = _dd2.assign(spread=_dd2[cols].max(axis=1) - _dd2[cols].min(axis=1))
    disc_chips = _dd2.loc[mismatch_mask].sort_values("spread", ascending=False)["target"].head(6).tolist()
    stable_chips = _dd2.loc[(~mismatch_mask) & (base != 0)].sort_values("spread")["target"].head(6).tolist()

    return [
        BoardSpec(
            tab="ΔΔCt · Comparativo",
            badges=[
                Badge(label="Comparativo ΔΔCt", tone="primary"),
                Badge(label=f"Cobertura {coverage:.0%}", tone="success"),
                Badge(label=f"Discrepancias de signo: {ddct_mismatch}", tone=("warning" if ddct_mismatch else "success")),
            ],
            left=Card(title="Coinciden (signo)", count=ddct_match, chips=stable_chips),
            right=Card(title="Discrepan (signo)", count=ddct_mismatch, chips=disc_chips),
            banner_text=(
                "Sin discrepancias detectadas entre métodos en ΔΔCt." if ddct_mismatch == 0 else None
            ),
            banner_tone="success",
            variant="light",
            key=f"ddct-board::comparative::{data.context_key}",
        )
    ]


# ---------- Builders: FC ----------


def build_fc_specs(data: InsightsData) -> List[BoardSpec]:
    df = data.fc_table
    if df is None or df.empty:
        return []
    if "target" not in df.columns:
        return []
    res: List[BoardSpec] = []
    for m_key, (m_label, fc_col, l2_col) in data.fc_cols.items():
        if fc_col not in df.columns or l2_col not in df.columns:
            continue
        fc_s = df[fc_col]
        l2_s = df[l2_col]
        n_eval = int(fc_s.notna().sum())
        coverage = float(fc_s.notna().mean())
        over2 = int((fc_s >= float(data.fc_over)).sum())
        under05 = int((fc_s <= float(data.fc_under)).sum())
        above1 = int((l2_s.abs() >= float(data.log2fc_abs)).sum())
        top_over_fc = df.sort_values(fc_col, ascending=False)["target"].head(6).tolist()
        top_under_fc = df.sort_values(fc_col, ascending=True)["target"].head(6).tolist()
        res.append(
            BoardSpec(
                tab=f"FC · {m_label}",
                badges=[
                    Badge(label=f"Fold Change · {m_label}", tone="primary"),
                    Badge(label=f"Cobertura {coverage:.0%}", tone="success"),
                    Badge(label=f"Genes evaluados: {n_eval}", tone="muted"),
                ],
                left=Card(title=f"Sobreexpresión (≥ {data.fc_over:g}×)", count=over2, chips=top_over_fc),
                right=Card(title=f"Subexpresión (≤ {data.fc_under:g}×)", count=under05, chips=top_under_fc),
                footer_text=f"Umbral rápido: |log2FC| ≥ {data.log2fc_abs:g} → {above1} genes ({m_label}).",
                variant="light",
                key=f"fc-board::{data.context_key}::{m_key}",
            )
        )
    return res


def build_fc_comparative_spec(data: InsightsData) -> List[BoardSpec]:
    df = data.fc_table
    if df is None or df.empty:
        return []
    cols = [t[2] for t in data.fc_cols.values() if t[2] in df.columns]
    if len(cols) < 2:
        return []
    _fc = df[cols].copy()
    coverage = float(_fc.notna().all(axis=1).mean())
    sgn_cols = {c: _fc[c].apply(_sgn) for c in cols}
    base = sgn_cols[cols[0]]
    mismatch_mask = np.zeros(len(base), dtype=bool)
    for c in cols[1:]:
        mismatch_mask |= (sgn_cols[c] != base) & (base != 0)
    fc_mismatch = int(mismatch_mask.sum())
    fc_match = int((((pd.DataFrame(sgn_cols).eq(base, axis=0)).all(axis=1)) & (base != 0)).sum())

    _fc2 = df[["target", *cols]].copy()
    _fc2 = _fc2.assign(spread=_fc2[cols].max(axis=1) - _fc2[cols].min(axis=1))
    disc_chips_fc = _fc2.loc[mismatch_mask].sort_values("spread", ascending=False)["target"].head(6).tolist()
    stable_chips_fc = _fc2.loc[(~mismatch_mask) & (base != 0)].sort_values("spread")["target"].head(6).tolist()

    return [
        BoardSpec(
            tab="FC · Comparativo",
            badges=[
                Badge(label="Comparativo log2FC", tone="primary"),
                Badge(label=f"Cobertura {coverage:.0%}", tone="success"),
                Badge(label=f"Discrepancias de signo: {fc_mismatch}", tone=("warning" if fc_mismatch else "success")),
            ],
            left=Card(title="Coinciden (signo)", count=fc_match, chips=stable_chips_fc),
            right=Card(title="Discrepan (signo)", count=fc_mismatch, chips=disc_chips_fc),
            banner_text=(
                "Sin discrepancias detectadas entre métodos en log2FC." if fc_mismatch == 0 else None
            ),
            banner_tone="success",
            variant="light",
            key=f"fc-board::comparative::{data.context_key}",
        )
    ]


# ---------- Default assembly ----------


def build_default_insight_specs(
    data: InsightsData,
    include: Sequence[str] = (
        "extraction_main",
        "extraction_targets",
        "ddct_per_method",
        "ddct_comparative",
        "fc_per_method",
        "fc_comparative",
    ),
) -> List[BoardSpec]:
    specs: List[BoardSpec] = []
    inc = set(include or [])
    if {"extraction_main", "extraction_targets"} & inc:
        specs.extend(build_extraction_specs(data))
        if "extraction_main" not in inc:
            specs = [s for s in specs if not s.key.startswith("extraction-board-main::")]
        if "extraction_targets" not in inc:
            specs = [s for s in specs if not s.key.startswith("extraction-board-targets::")]
    if "ddct_per_method" in inc:
        specs.extend(build_ddct_specs(data))
    if "ddct_comparative" in inc:
        specs.extend(build_ddct_comparative_spec(data))
    if "fc_per_method" in inc:
        specs.extend(build_fc_specs(data))
    if "fc_comparative" in inc:
        specs.extend(build_fc_comparative_spec(data))
    return specs


# ---------- Concordance summary (log2FC) ----------


def build_concordance_spec_from_summary(
    summary: dict,
    *,
    context_key: str,
    key_prefix: str = "concordance",
) -> List[BoardSpec]:
    """Create a single insight board summarizing pairwise metrics and discrepancies.

    - summary: output of summarize_fc_methods(fc_table)
    - context_key: used to namespace the board key
    """

    if not isinstance(summary, dict):
        return []
    pairwise = summary.get("pairwise_metrics", {}) or {}
    counts = summary.get("counts", {}) or {}
    disc = summary.get("top_discrepancies", []) or []

    # Build pair chips: label with r and error metrics
    pair_labels = {
        ("advanced", "promedio"): "Avanzada vs Promedio",
        ("advanced", "gen_ref"): "Avanzada vs Gen ref",
        ("promedio", "gen_ref"): "Promedio vs Gen ref",
    }
    pair_chips: List[str] = []
    r_values: List[float] = []
    n_values: List[int] = []
    for pair, metrics in pairwise.items():
        label = pair_labels.get(tuple(pair), f"{pair}")
        try:
            r = float(metrics.get("pearson_r", float("nan")))
        except Exception:
            r = float("nan")
        rmse = metrics.get("rmse", float("nan"))
        mae = metrics.get("mae", float("nan"))
        n = int(metrics.get("n", 0))
        r_values.append(r if r == r else 0.0)  # handle NaN
        n_values.append(n)
        pair_chips.append(f"{label} · r={r:.2f} · RMSE={rmse:.2f} · MAE={mae:.2f} (N={n})")

    # Top discrepancy genes chips
    top_genes = [str(item.get("gene", "?")) for item in (disc or [])][:6]

    # Threshold summary for ≥2× and ≥1.5× if available
    thr2 = counts.get(">=2.0x") or counts.get(">=2x") or {}
    thr15 = counts.get(">=1.5x") or {}
    footer_bits = []
    if thr2:
        footer_bits.append(
            f"≥2× → Todos: {thr2.get('todos', 0)} · Alguno: {thr2.get('alguno', 0)}"
        )
    if thr15:
        footer_bits.append(
            f"≥1.5× → Todos: {thr15.get('todos', 0)} · Alguno: {thr15.get('alguno', 0)}"
        )
    footer_text = " | ".join(footer_bits) if footer_bits else None

    # Badges with quick glance metrics
    avg_r = (sum(r_values) / len(r_values)) if r_values else 0.0
    n_min = min(n_values) if n_values else 0

    spec = BoardSpec(
        tab="Concordancia",
        badges=[
            Badge(label="Concordancia log2FC", tone="primary"),
            Badge(label=f"r medio {avg_r:.2f}", tone=("success" if avg_r >= 0.9 else "muted")),
            Badge(label=f"N min: {n_min}", tone="muted"),
        ],
        left=Card(title="Pares (r, RMSE, MAE)", count=len(pair_chips), chips=pair_chips, chips_max=None),
        right=Card(title="Top discrepancias", count=len(top_genes), chips=top_genes),
        footer_text=footer_text,
        variant="light",
        key=f"{key_prefix}::{context_key}",
    )
    return [spec]


# ---------- Rendering ----------


def render_insights_overview(
    specs: Iterable[BoardSpec],
    *,
    layout: str = "stack",  # or "tabs"
    key_prefix: str = "overview",
) -> None:
    items: List[BoardSpec] = [s for s in specs if s is not None]
    if not items:
        return

    layout = (layout or "stack").lower()
    cid = f"insights-{_slug(key_prefix)}"

    if layout == "tabs":
        # Build tab labels; fall back to a readable label from the first badge or key
        tab_labels = []
        for s in items:
            if s.tab:
                tab_labels.append(str(s.tab))
            elif s.badges:
                tab_labels.append(str(getattr(s.badges[0], "label", "")) or s.key)
            else:
                tab_labels.append(s.key)
        tabs = st.tabs(tab_labels)
        for s, tab in zip(items, tabs):
            with tab:
                render_insight_board(
                    badges=s.badges,
                    left=s.left,
                    right=s.right,
                    footer_text=s.footer_text,
                    banner_text=s.banner_text,
                    banner_tone=s.banner_tone,
                    variant=s.variant,
                    key=s.key,
                )
        return

    # Default: stack
    st.markdown(
        f"""
        <style>
        #{cid} {{ display:flex; flex-direction:column; gap: 10px; }}
        </style>
        <div id='{cid}'></div>
        """,
        unsafe_allow_html=True,
    )
    for s in items:
        render_insight_board(
            badges=s.badges,
            left=s.left,
            right=s.right,
            footer_text=s.footer_text,
            banner_text=s.banner_text,
            banner_tone=s.banner_tone,
            variant=s.variant,
            key=s.key,
        )
