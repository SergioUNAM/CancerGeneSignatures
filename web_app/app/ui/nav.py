from __future__ import annotations

import streamlit as st


def enable_nav_css() -> None:
    st.markdown(
        """
        <style>
        .cgs-nav { display:flex; gap:6px; border:1px solid #e5e5ea; padding:4px; border-radius:12px; background:#f2f2f7; width: fit-content; }
        .cgs-nav .seg { border:none; background:transparent; padding:6px 14px; border-radius:10px; cursor:pointer; color:#1c1c1e; }
        .cgs-nav .seg.active { background:#1f6feb; color:white; }
        .cgs-nav-container { display:flex; justify-content:center; margin: 8px 0 18px 0; }
        .cgs-next { display:flex; justify-content:flex-end; margin-top: 16px; }
        .cgs-next .btn { background:#1f6feb; color:white; border:none; padding:8px 16px; border-radius:10px; cursor:pointer; }
        .cgs-next .btn[disabled] { background:#a8b3cf; cursor:not-allowed; }
        </style>
        """,
        unsafe_allow_html=True,
    )


def render_nav_picker(active: str = "entrada") -> None:
    """Renderiza un selector de navegación (segmentado) para páginas multipágina.

    active: "entrada" | "normalizacion"
    """
    enable_nav_css()
    st.markdown('<div class="cgs-nav-container"><div class="cgs-nav">', unsafe_allow_html=True)
    col1, col2, _ = st.columns([1, 1, 4])
    with col1:
        if st.button("Entrada", key="nav_ent", help="Entrada y Clasificación", width="stretch"):
            if hasattr(st, "switch_page"):
                try:
                    st.switch_page("pages/1_Entrada_y_Clasificacion.py")
                except Exception:
                    pass
    with col2:
        if st.button("Normalización", key="nav_norm", help="Normalización y Fold Change", width="stretch"):
            if hasattr(st, "switch_page"):
                try:
                    st.switch_page("pages/2_Normalizacion_y_FoldChange.py")
                except Exception:
                    pass
    st.markdown('</div></div>', unsafe_allow_html=True)

    # Fallback si no existe switch_page
    if not hasattr(st, "switch_page"):
        st.caption("Navegación rápida")
        st.page_link("pages/1_Entrada_y_Clasificacion.py", label="Entrada y Clasificación")
        st.page_link("pages/2_Normalizacion_y_FoldChange.py", label="Normalización y Fold Change")


def render_next_button(current: str, *, enabled: bool = True) -> None:
    """Muestra un botón 'Siguiente' que avanza a la próxima vista del pipeline.

    current: "entrada" | "normalizacion"
    enabled: habilita o deshabilita el botón según precondiciones.
    """
    enable_nav_css()
    col = st.container()
    with col:
        st.markdown('<div class="cgs-next">', unsafe_allow_html=True)
        clicked = st.button(
            "Siguiente",
            key=f"next_{current}",
            disabled=not enabled,
            help="Avanza a la siguiente etapa del pipeline",
        )
        st.markdown('</div>', unsafe_allow_html=True)
        if clicked and enabled and hasattr(st, "switch_page"):
            try:
                target = (
                    "pages/2_Normalizacion_y_FoldChange.py" if current == "entrada" else "pages/2_Normalizacion_y_FoldChange.py"
                )
                st.switch_page(target)
            except Exception:
                pass


__all__ = ["enable_nav_css", "render_nav_picker", "render_next_button"]
