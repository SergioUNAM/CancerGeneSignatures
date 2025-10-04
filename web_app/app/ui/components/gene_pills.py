"""Integrated gene pills view for ΔΔCt/FC (3 methods).

Renders a filterable grid of "gene pills" showing direction and magnitude
of expression changes using log2FC across the three methods
(`advanced`, `promedio`, `gen_ref`).

Inputs: a consolidated `fc_table` with at least these columns:
 - `target`
 - `log2fc_advanced`, `log2fc_promedio`, `log2fc_gen_ref`
Optionally supports q-value columns: `q_advanced`, `q_promedio`, `q_gen_ref`.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple
from html import escape as _escape

import math
import re
import pandas as pd
import numpy as np
import streamlit as st
import plotly.graph_objects as go
from streamlit.components.v1 import html as st_html

__all__ = ["render_gene_pills_fc"]
__all__.append("render_gene_methods_view")


def _slug(s: str) -> str:
    s = re.sub(r"[^a-zA-Z0-9_-]+", "-", str(s))
    if s and s[0].isdigit():
        s = "id-" + s
    return s


_METHOD_META = {
    "advanced": dict(label="Avanzada", shape="diamond"),
    "promedio": dict(label="Promedio", shape="square"),
    "gen_ref": dict(label="Gen ref", shape="circle"),
}


def _sgn(x: float, eps: float = 0.25) -> int:
    try:
        v = float(x)
    except Exception:
        return 0
    if not math.isfinite(v):
        return 0
    if abs(v) <= eps:
        return 0
    return 1 if v > 0 else -1


def _extract_log2fc_columns(df: pd.DataFrame) -> Dict[str, str]:
    mapping = {}
    for key in ("advanced", "promedio", "gen_ref"):
        col = f"log2fc_{key}"
        if col in df.columns:
            mapping[key] = col
    return mapping


def _extract_q_columns(df: pd.DataFrame) -> Dict[str, str]:
    mapping = {}
    for key in ("advanced", "promedio", "gen_ref"):
        for candidate in (f"q_{key}", f"qvalue_{key}", f"padj_{key}"):
            if candidate in df.columns:
                mapping[key] = candidate
                break
    return mapping


def render_gene_pills_fc(
    fc_table: pd.DataFrame,
    *,
    key: str,
    max_cap: float = 3.0,  # caps bar length at |log2FC| == 3
) -> None:
    if not isinstance(fc_table, pd.DataFrame) or fc_table.empty or "target" not in fc_table.columns:
        st.info("No hay datos para mostrar el explorador de genes.")
        return

    l2_cols = _extract_log2fc_columns(fc_table)
    if not l2_cols:
        st.info("La tabla no contiene columnas log2FC necesarias.")
        return
    q_cols = _extract_q_columns(fc_table)

    cid = f"genes-{_slug(key)}"
    with st.container():
        st.markdown("### Explorador de genes (3 métodos integrados)")

        # Filtros básicos (condensados)
        fc1, fc2 = st.columns([3, 2], gap="small")
        with fc1:
            term = st.text_input("Buscar gen", value="", key=f"q::{key}")
        with fc2:
            thr = st.slider(
                "Umbral |log2FC|",
                min_value=0.0,
                max_value=float(max_cap),
                value=0.0,
                step=0.05,
                key=f"thr::{key}",
                help="Filtro por magnitud del cambio: |log2FC|≥umbral. Referencias: 1≈2×, 1.5≈2.8×, 2≈4×.",
            )

        # Filtros avanzados (plegados por defecto)
        with st.expander("Más filtros", expanded=False):
            ac1, ac2, ac3 = st.columns([1, 1, 2], gap="small")
            with ac1:
                dir_opt = st.selectbox(
                    "Dirección",
                    options=["Todos", "Sobreexpresados", "Subexpresados", "Neutros"],
                    key=f"dir::{key}",
                )
            with ac2:
                method_names = ["Consenso (mediana log2FC)"] + [
                    _METHOD_META[m]["label"] for m in l2_cols.keys()
                ]
                method_label = st.selectbox(
                    "Método",
                    options=method_names,
                    key=f"method::{key}",
                )
            with ac3:
                only_disc = st.checkbox(
                    "Ver solo discrepancias entre métodos",
                    value=False,
                    key=f"disc_only::{key}",
                )
            # Paginación (compacta)
            page_size = 60
            st.caption("Paginación automática (60 por página)")

        # Prepare data
        data = fc_table[["target", *l2_cols.values(), *q_cols.values()]].copy()
        # consensus/selected value
        if method_label.startswith("Consenso"):
            vals = data[list(l2_cols.values())].median(axis=1, skipna=True)
        else:
            # map back from label to key
            inv = {v["label"]: k for k, v in _METHOD_META.items()}
            sel_key = inv.get(method_label, "advanced")
            vals = data[l2_cols.get(sel_key)].astype(float)
        data = data.assign(_val=vals)

        # derive signs by method and consensus
        signs_by_method = {}
        for m, col in l2_cols.items():
            signs_by_method[m] = data[col].apply(_sgn)
        data["_sgn"] = data["_val"].apply(_sgn)
        # discrepancies mask: not all equal and at least one non-zero
        if len(signs_by_method) >= 2:
            _mvals = list(signs_by_method.values())
            base = _mvals[0]
            neq = False
            for s in _mvals[1:]:
                neq = neq | (s != base)
            disc_mask = neq & (base != 0)
        else:
            disc_mask = pd.Series(False, index=data.index)

        # coverage = fraction with all available methods notna
        cov = float(data[list(l2_cols.values())].notna().all(axis=1).mean())

        # Summary (línea compacta)
        over_n = int((data["_sgn"] > 0).sum())
        under_n = int((data["_sgn"] < 0).sum())
        neutral_n = int((data["_sgn"] == 0).sum())
        # Mostrar modo y umbral activos
        mode_txt = method_label if 'method_label' in locals() else 'Consenso (mediana log2FC)'
        st.caption(
            f"Vista: {mode_txt} · Umbral |log2FC| ≥ {thr:.2f}  —  "
            f"Sobreexpresados: {over_n} · Subexpresados: {under_n} · Neutros: {neutral_n} · Cobertura: {cov:.0%}"
        )

        # Apply filters
        mask = pd.Series(True, index=data.index)
        if dir_opt == "Sobreexpresados":
            mask &= data["_sgn"] > 0
        elif dir_opt == "Subexpresados":
            mask &= data["_sgn"] < 0
        elif dir_opt == "Neutros":
            mask &= data["_sgn"] == 0
        if thr > 0:
            mask &= data["_val"].abs() >= float(thr)
        if term:
            term_up = term.strip().upper()
            mask &= data["target"].astype(str).str.upper().str.contains(re.escape(term_up))
        if only_disc:
            mask &= disc_mask

        view = data.loc[mask, ["target", "_val"] + list(l2_cols.values())].copy()
        view = view.sort_values("target").reset_index(drop=True)

        total = len(view)
        pages = max(1, int(math.ceil(total / float(page_size))))
        page = st.slider(
            "Página",
            min_value=1,
            max_value=pages,
            value=1,
            step=1,
            key=f"page::{key}",
        )
        start = (page - 1) * int(page_size)
        end = start + int(page_size)
        show = view.iloc[start:end]

        # Render grid
        pills_html = []
        for _, row in show.iterrows():
            gene = str(row["target"]) if pd.notna(row["target"]) else ""
            val = float(row["_val"]) if pd.notna(row["_val"]) else 0.0
            width_pct = max(6.0, min(100.0, (abs(val) / float(max_cap)) * 100.0))
            # color mapping
            if _sgn(val) > 0:
                base = (46, 125, 50)   # #2e7d32
                fade = (165, 214, 167) # #a5d6a7
            elif _sgn(val) < 0:
                base = (94, 53, 177)   # #5e35b1
                fade = (209, 196, 233) # #d1c4e9
            else:
                base = (176, 176, 176) # neutral
                fade = (220, 220, 220)
            bar_bg = f"linear-gradient(90deg, rgba({base[0]},{base[1]},{base[2]},0.95), rgba({fade[0]},{fade[1]},{fade[2]},0.85))"

            # mini method icons
            icons_html = []
            for m, meta in _METHOD_META.items():
                if m not in l2_cols:
                    continue
                mval = row[l2_cols[m]]
                ms = _sgn(mval)
                if ms > 0:
                    color = "rgba(46,125,50,0.95)"
                elif ms < 0:
                    color = "rgba(94,53,177,0.95)"
                else:
                    color = "rgba(176,176,176,0.9)"
                if meta["shape"] == "square":
                    shape_html = f"<span class='sq' style='background:{color}'></span>"
                elif meta["shape"] == "circle":
                    shape_html = f"<span class='ci' style='background:{color}'></span>"
                else:  # diamond
                    shape_html = f"<span class='di' style='background:{color}'></span>"
                icons_html.append(shape_html)

            # tooltip
            tip_parts = [gene]
            for m, col in l2_cols.items():
                tip_parts.append(f"{_METHOD_META[m]['label']}: {row[col]:.2f}")
            title_attr = " | ".join(tip_parts)
            title_attr_esc = _escape(title_attr, quote=True)

            pills_html.append(
                f"""
                <div class='pill' title="{title_attr_esc}">
                  <div class='header'>
                    <span class='icons'>{''.join(icons_html)}</span>
                    <span class='gene'>{_escape(gene)}</span>
                  </div>
                  <div class='bar'>
                    <span class='fill' style='width:{width_pct:.1f}%; background:{bar_bg}'></span>
                  </div>
                </div>
                """
            )

        html = "".join(pills_html)
        # Approximate height: 60px per row, assume 4 columns minimum
        approx_rows = max(1, int(math.ceil(max(1, len(show)) / 6.0)))
        height = int(120 + approx_rows * 56)
        st_html(
            f"""
            <style>
            #{cid} .grid {{
                display:grid; grid-template-columns: repeat(auto-fill, minmax(180px, 1fr));
                gap: 8px; align-items:start;
            }}
            #{cid} .pill {{
                background: linear-gradient(180deg, rgba(241,245,249,0.96), rgba(241,245,249,0.80));
                border:1px solid rgba(15, 23, 42, 0.10); border-radius: 10px; padding: 6px 8px;
                box-shadow: 0 6px 14px rgba(15,23,42,0.07);
            }}
            #{cid} .header {{ display:flex; align-items:center; gap:6px; margin-bottom:2px; }}
            #{cid} .gene {{
                display:inline-block; padding:2px 8px; border-radius:999px; font-weight:600; font-size:11px;
                background: linear-gradient(180deg, rgba(203,213,225,0.45), rgba(203,213,225,0.28));
                border:1px solid rgba(148,163,184,0.22); color:#0b1220;
            }}
            #{cid} .bar {{ position: relative; height: 3px; background: rgba(148,163,184,0.25); border-radius: 4px; overflow:hidden; }}
            #{cid} .fill {{ position:absolute; left:0; top:0; bottom:0; border-radius: 4px; }}
            #{cid} .icons {{ display:inline-flex; gap:6px; align-items:center; }}
            #{cid} .sq {{ width:7px; height:7px; display:inline-block; border-radius:2px; }}
            #{cid} .ci {{ width:7px; height:7px; display:inline-block; border-radius:999px; }}
            #{cid} .di {{ width:7px; height:7px; display:inline-block; transform: rotate(45deg); border-radius:2px; }}
            #{cid} .legend {{
                margin-top:6px; display:flex; flex-wrap:wrap; gap:10px; font-size:12px; color:#334155;
                align-items:center;
            }}
            #{cid} .lg-chip {{ display:inline-flex; align-items:center; gap:6px; }}
            #{cid} .lg-swatch {{ width:12px; height:6px; border-radius:4px; display:inline-block; }}
            </style>
            <div id='{cid}'>
              <div class='grid'>{html}</div>
              <div class='legend'>
                <span class='lg-chip'><span class='sq' style='background:rgba(15,23,42,0.6); width:9px; height:9px;'></span> Promedio</span>
                <span class='lg-chip'><span class='di' style='background:rgba(15,23,42,0.6); width:9px; height:9px;'></span> Avanzada</span>
                <span class='lg-chip'><span class='ci' style='background:rgba(15,23,42,0.6); width:9px; height:9px;'></span> Gen ref</span>
                <span class='lg-chip'><span class='lg-swatch' style='background:linear-gradient(90deg, rgba(46,125,50,0.95), rgba(165,214,167,0.85))'></span> Sobreexpresión</span>
                <span class='lg-chip'><span class='lg-swatch' style='background:linear-gradient(90deg, rgba(94,53,177,0.95), rgba(209,196,233,0.85))'></span> Subexpresión</span>
                <span class='lg-chip'><span class='lg-swatch' style='background:rgba(176,176,176,0.9)'></span> Neutro</span>
              </div>
            </div>
            """,
            height=height,
        )

        # Mini lollipop for a selected gene (manual selection for now)
        if not show.empty:
            sel_gene = st.selectbox(
                "Detalle (mini‑lollipop) para gen",
                options=show["target"].astype(str).tolist(),
                key=f"gene_detail::{key}::{page}",
            )
            row_all = data.loc[data["target"].astype(str) == sel_gene].head(1)
            if not row_all.empty:
                xs = []
                ys = []
                texts = []
                colors = []
                for m in ("promedio", "advanced", "gen_ref"):
                    if m in l2_cols:
                        xs.append(_METHOD_META[m]["label"])
                        val = float(row_all[l2_cols[m]].iloc[0]) if pd.notna(row_all[l2_cols[m]].iloc[0]) else 0.0
                        ys.append(val)
                        texts.append(f"{_METHOD_META[m]['label']}: {val:.2f}")
                        colors.append("#2e7d32" if val > 0 else ("#5e35b1" if val < 0 else "#888"))
                fig = go.Figure()
                fig.add_trace(go.Scatter(x=xs, y=ys, mode="lines+markers", marker=dict(size=10, color=colors)))
                fig.update_layout(
                    title=dict(text=f"{sel_gene} · log2FC por método", x=0.5, xanchor="center"),
                    template="plotly_white",
                    yaxis_title="log2FC",
                    height=260,
                    margin=dict(t=50, b=30, l=40, r=20),
        )
                st.plotly_chart(fig, use_container_width=True, key=f"mini_lollipop::{key}::{sel_gene}")


def render_gene_methods_view(
    fc_table: pd.DataFrame,
    *,
    key: str,
    max_cap: float = 3.0,
    eps: float = 0.25,
) -> None:
    """Minimal viewer of per-gene expression across methods (no filters).

    Shows, for cada gen, tres micro‑barras (Promedio, Avanzada, Gen ref)
    codificando dirección y magnitud (|log2FC|) con un cap de `max_cap`.
    """

    if not isinstance(fc_table, pd.DataFrame) or fc_table.empty or "target" not in fc_table.columns:
        st.info("No hay datos para el visor simple.")
        return

    l2_cols = _extract_log2fc_columns(fc_table)
    if not l2_cols:
        st.info("Faltan columnas de log2FC para el visor simple.")
        return

    q_cols = _extract_q_columns(fc_table)

    data = fc_table[["target", *l2_cols.values(), *q_cols.values()]].copy()
    data = data.sort_values("target").reset_index(drop=True)

    # Paginación simple incremental (sin controles): 90 → 180 → ...
    page_key = f"simple_page::{key}"
    show_n = int(st.session_state.get(page_key, 90))
    view = data.head(show_n)

    cid = f"simple-{_slug(key)}"

    # Build pills HTML
    pills_html: List[str] = []
    for _, row in view.iterrows():
        gene = str(row["target"]) if pd.notna(row["target"]) else ""
        icons_html: List[str] = []
        hover_rows: List[str] = []
        for m in ("promedio", "advanced", "gen_ref"):
            if m not in l2_cols:
                continue
            val = row[l2_cols[m]]
            try:
                fv = float(val)
            except Exception:
                fv = 0.0
            # color by direction; intensity by |log2FC|
            ratio = max(0.0, min(1.0, abs(fv) / float(max_cap)))
            alpha = 0.25 + 0.7 * ratio
            if _sgn(fv, eps) > 0:
                color = f"rgba(46,125,50,{alpha:.3f})"      # green
            elif _sgn(fv, eps) < 0:
                color = f"rgba(94,53,177,{alpha:.3f})"      # purple
            else:
                color = f"rgba(176,176,176,{max(0.25, 0.4*alpha):.3f})"  # neutral
            # icon shape element
            if _METHOD_META[m]["shape"] == "square":
                icon = f"<span class='sq' style='background:{color}'></span>"
            elif _METHOD_META[m]["shape"] == "circle":
                icon = f"<span class='ci' style='background:{color}'></span>"
            else:
                icon = f"<span class='di' style='background:{color}'></span>"
            # detailed row for hover card
            qtxt = ""
            qcol = q_cols.get(m)
            if qcol and qcol in row and pd.notna(row[qcol]):
                try:
                    qtxt = f" · q={float(row[qcol]):.3f}"
                except Exception:
                    qtxt = ""
            hover_rows.append(
                f"<div class='trow'><span class='shape'>{icon}</span><span class='mlabel'>{_METHOD_META[m]['label']}</span><span class='mval'>{fv:+.2f}</span><span class='mq'>{qtxt}</span></div>"
            )
            title_attr = _escape(f"{_METHOD_META[m]['label']}: {fv:+.2f}{qtxt}", quote=True)
            icons_html.append(f"<span class='micon' title=\"{title_attr}\">{icon}</span>")
        pills_html.append(
            f"""
            <div class='pill'>
              <div class='header'><span class='gene'>{_escape(gene)}</span></div>
              <div class='ind'>{''.join(icons_html)}</div>
              <div class='tinfo'>{''.join(hover_rows)}</div>
            </div>
            """
        )

    html = "".join(pills_html)
    from streamlit.components.v1 import html as st_html
    approx_rows = max(1, int(math.ceil(max(1, len(view)) / 6.0)))
    height = int(90 + approx_rows * 42)
    st_html(
        f"""
        <style>
        #{cid} .grid {{ display:grid; grid-template-columns: repeat(auto-fill, minmax(170px, 1fr)); gap:6px; }}
        #{cid} .pill {{ position:relative; background: linear-gradient(180deg, rgba(241,245,249,0.97), rgba(241,245,249,0.82));
                        border:1px solid rgba(15,23,42,0.10); border-radius:8px; padding:5px 6px; box-shadow:0 5px 12px rgba(15,23,42,0.06); }}
        #{cid} .header {{ display:flex; align-items:center; justify-content:flex-start; margin-bottom:3px; }}
        #{cid} .gene {{ display:inline-block; padding:2px 6px; border-radius:999px; font-weight:700; letter-spacing:.02em; font-size:10.5px; background: linear-gradient(180deg, rgba(203,213,225,0.45), rgba(203,213,225,0.28)); border:1px solid rgba(148,163,184,0.22); color:#0b1220; }}
        #{cid} .ind {{ display:flex; align-items:center; gap:6px; justify-content:flex-start; padding:1px 0 0 1px; min-height:12px; }}
        #{cid} .micon {{}}
        #{cid} .sq {{ width:8px; height:8px; display:inline-block; border-radius:2px; }}
        #{cid} .ci {{ width:8px; height:8px; display:inline-block; border-radius:999px; }}
        #{cid} .di {{ width:8px; height:8px; display:inline-block; transform: rotate(45deg); border-radius:2px; }}
        /* Hover card */
        #{cid} .tinfo {{
            position:absolute; left:6px; right:6px; bottom:6px;
            background: rgba(255,255,255,0.98); border:1px solid rgba(15,23,42,0.10);
            border-radius:8px; box-shadow: 0 10px 24px rgba(15,23,42,0.10);
            padding:6px 8px; opacity:0; transform: translateY(4px);
            pointer-events:none; transition: opacity .16s ease, transform .16s ease;
            z-index:5;
        }}
        #{cid} .pill:hover .tinfo {{ opacity:1; transform: translateY(0); }}
        #{cid} .trow {{ display:flex; align-items:center; gap:8px; font-size:11px; line-height:1.15; color:#0b1220; }}
        #{cid} .trow + .trow {{ margin-top:4px; }}
        #{cid} .trow .shape {{ display:inline-flex; width:10px; justify-content:center; }}
        #{cid} .trow .mlabel {{ flex: 0 0 72px; color:#475569; }}
        #{cid} .trow .mval {{ font-weight:700; }}
        #{cid} .trow .mq {{ color:#475569; }}
        #{cid} .legend {{ margin-top:6px; font-size:12px; color:#334155; display:flex; gap:10px; flex-wrap:wrap; align-items:center; }}
        #{cid} .lg-chip {{ display:inline-flex; align-items:center; gap:6px; }}
        #{cid} .lg-swatch {{ width:12px; height:6px; border-radius:4px; display:inline-block; }}
        </style>
        <div id='{cid}'>
          <div class='grid'>{html}</div>
          <div class='legend'>
            <span class='lg-chip'><span class='sq' style='background:rgba(15,23,42,0.65)'></span> Promedio</span>
            <span class='lg-chip'><span class='di' style='background:rgba(15,23,42,0.65)'></span> Avanzada</span>
            <span class='lg-chip'><span class='ci' style='background:rgba(15,23,42,0.65)'></span> Gen ref</span>
            <span class='lg-chip'><span class='lg-swatch' style='background:rgba(46,125,50,0.95)'></span> Sobreexpresión (verde)</span>
            <span class='lg-chip'><span class='lg-swatch' style='background:rgba(94,53,177,0.95)'></span> Subexpresión (morado)</span>
            <span class='lg-chip'><span class='lg-swatch' style='background:rgba(176,176,176,0.9)'></span> Neutro (gris)</span>
            <span class='lg-chip'>Intensidad ≈ |log2FC|</span>
          </div>
        </div>
        """,
        height=height,
    )

    # "Mostrar más" discreto (no es un filtro; aumenta elementos renderizados)
    if len(data) > show_n:
        if st.button("Mostrar más", key=f"more::{key}"):
            st.session_state[page_key] = show_n + 90
