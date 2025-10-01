# Guía de Estilos (UI) — CancerGeneSignatures

Este documento define criterios visuales y de interacción para mantener consistencia y accesibilidad en la web app.

## Principios
- Consistencia tipográfica: usar la tipografía por defecto del sistema (Streamlit) y evitar mezclar familias sin necesidad.
- Jerarquía visual: títulos claros (st.title, st.subheader), secciones separadas por `st.divider()` y textos de apoyo concisos (`st.caption`).
- Accesibilidad: contrastes suficientes, botones identificables, estados deshabilitados explícitos, y sin dependencia exclusiva de color para transmitir estado.
- Navegación clara: siempre visible y consistente entre páginas.
- Prohibición de emojis: no usar emojis en textos, títulos o botones. Preferir rótulos descriptivos y, en su caso, iconografía accesible.

## Componentes
- Navegación (segmented control): definida en `app/ui/nav.py` como `render_nav_picker`.
  - Estilos CSS: bordes suaves, fondo neutro, estado activo diferenciado por color y fondo.
  - Botón “Siguiente”: `render_next_button`, se coloca al final de la vista y se deshabilita si faltan precondiciones.
- Formularios y selectores: usar controles nativos de Streamlit; sin estilos adicionales para evitar fricción y mantener accesibilidad.
- Tablas y gráficos: usar `st.dataframe` y Plotly con plantillas legibles; títulos informativos y leyendas claras.

## Colores
- Primario: `#1f6feb` (estado activo/acción principal).
- Neutros: contornos `#e5e5ea`, fondos `#f2f2f7`, textos principales `#1c1c1e`.
- Estados deshabilitados: `#a8b3cf`.

## Patrón de Navegación
- Landing (Home): `web_app/Home.py` con descripción general y enlaces a páginas.
- Entrada y Clasificación: `pages/1_Entrada_y_Clasificacion.py`.
- Normalización y Fold Change: `pages/2_Normalizacion_y_FoldChange.py`.
- La navegación superior (`render_nav_picker`) debe mostrarse en todas las páginas del flujo.
- El botón “Siguiente” sólo se activa cuando se cumplen precondiciones (p. ej., clasificación válida en Entrada).

## Trazabilidad
- Estilos implementados en: `web_app/app/ui/nav.py`.
- Páginas que consumen estos estilos: `web_app/Home.py`, `web_app/pages/1_Entrada_y_Clasificacion.py`, `web_app/pages/2_Normalizacion_y_FoldChange.py`.
- Cambios de estilo deben documentarse en este archivo con una nota de versión y el motivo.

## Buenas prácticas adicionales
- Evitar columnas excesivas; preferir 1–2 columnas para formularios.
- Mensajes accionables en errores/avisos (`st.error/warning/info`).
- Usar `st.download_button` para salidas clave (CSV) con nombres claros.
