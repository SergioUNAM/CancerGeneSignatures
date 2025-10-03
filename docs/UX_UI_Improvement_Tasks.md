# Propuestas de mejoras UX/UI para la web app

Este documento resume oportunidades de mejora de primer nivel en la experiencia de usuario (UX) y la interfaz (UI) de la aplicación Streamlit descrita en `web_app/streamlit_app.py`.

## 1. Guiado por etapas y reducción de carga cognitiva inicial
- **Problema**: La pantalla principal presenta un título y múltiples bloques de texto extensos (resumen, "Ruta rápida", política de imputación) antes de que el usuario ejecute alguna acción, lo que sobrecarga la atención y desplaza el contenido clave hacia abajo.【F:web_app/streamlit_app.py†L645-L704】
- **Propuesta**: Introducir un encabezado tipo "wizard" o un componente de progreso visible que organice las etapas (carga, clasificación, normalización, exportación) y mover la ayuda detallada a tooltips o acordeones colapsados por defecto. Esto facilita que perfiles científicos se concentren en el flujo sin perder el contexto.

## 2. Sidebar con estados y acciones dependientes
- **Problema**: La barra lateral muestra todos los selectores y botones simultáneamente, incluso cuando no se ha cargado archivo, y el botón "Procesar archivo" está disponible sin validaciones visuales previas.【F:web_app/streamlit_app.py†L666-L734】
- **Propuesta**: Convertir la barra lateral en un panel de pasos dependientes: deshabilitar controles hasta que el archivo se cargue correctamente, mostrar estados (p. ej., iconos de éxito/error por sección) y proveer confirmaciones visuales tras procesar el archivo.

## 3. Clasificación de controles/muestras más visual
- **Problema**: El módulo de clasificación combina múltiples `multiselect`, advertencias y métricas distribuidas verticalmente, requiriendo mucho scroll y relectura para validar la separación.【F:web_app/streamlit_app.py†L814-L875】【F:web_app/streamlit_app.py†L909-L975】
- **Propuesta**: Reorganizar el bloque en tarjetas o tabla editable con filtros interactivos, incorporar badges de conteo junto a cada prefijo y habilitar una vista previa sticky (p. ej., panel lateral derecho) que muestre diferencias en tiempo real sin depender de expanders.

## 4. Parámetros avanzados con jerarquía clara
- **Problema**: El expander "Parámetros" incluye explicaciones muy largas en formato Markdown plano, lo que dificulta identificar rápidamente campos críticos como FDR, candidatos o bootstrap.【F:web_app/streamlit_app.py†L909-L968】
- **Propuesta**: Agrupar parámetros en sub-secciones visuales (p. ej., tarjetas "Sensibilidad" vs. "Reproducibilidad"), resumir las descripciones con tooltips contextuales y ofrecer presets (Rápido, Equilibrado, Exhaustivo) para usuarios que no requieran ajustar cada control.

## 5. Descargas y resultados más accesibles
- **Problema**: Los resultados (heatmaps, datasets, genes DE) generan múltiples botones de descarga repartidos en columnas y expanders, lo que puede provocar que el usuario pierda archivos clave o repita descargas.【F:web_app/streamlit_app.py†L1065-L1167】【F:web_app/streamlit_app.py†L1170-L1275】
- **Propuesta**: Incorporar un "panel de exportación" fijo al final o en la barra lateral con enlaces consolidando las descargas más relevantes, acompañado de etiquetas de método y metadatos (fecha, parámetros usados).

## 6. Sección de anotación Ensembl con resumen compacto
- **Problema**: La consulta Ensembl requiere recordar ejecutar manualmente un botón y la tabla resultante desplaza el resumen numérico al lateral, dificultando evaluar rápidamente la cobertura de IDs obtenidos.【F:web_app/streamlit_app.py†L521-L603】
- **Propuesta**: Mostrar un estado principal (pendiente/completado) con tarjeta de KPIs, activar automáticamente la consulta cuando se detecta una nueva lista de genes y ofrecer búsqueda/filtrado directo sobre la tabla para agilizar exploración.

Cada tarea debe abordarse manteniendo el rigor científico del flujo, pero priorizando la claridad visual y la reducción de fricción durante la carga, configuración y exportación de resultados.

## Seguimiento de implementación

### Etapa 1 — Guiado por etapas y reducción de carga cognitiva inicial
- ✅ Se añadió un encabezado de progreso tipo wizard y las ayudas introductorias se trasladaron a expanders colapsados por defecto para aliviar la primera carga visual en la vista principal (`web_app/streamlit_app.py`).

### Etapa 2 — Sidebar con estados y acciones dependientes
- ✅ La barra lateral ahora renderiza el estado de cada etapa con indicadores compactos, bloquea selectores dependientes hasta que la carga de datos se completa y confirma visualmente el dataset activo tras cada procesamiento (`web_app/streamlit_app.py`).
- ✅ Los controles de imputación y etiquetado del estudio se trasladaron al cuerpo principal para configurarse con contexto y explicación inmediata, manteniendo la barra lateral enfocada en carga y estado global (`web_app/streamlit_app.py`).
- ✅ La carga del archivo también se movió al flujo principal, dejando la barra lateral exclusivamente como tablero de seguimiento del pipeline (`web_app/streamlit_app.py`).

### Etapa 3 — Clasificación de controles/muestras más visual
- ✅ La sección de clasificación ahora muestra un catálogo compacto con tabla resumen y selectores en tarjetas paralelas, mostrando badges de conteo por prefijo y reorganizando las asignaciones manuales para reducir el scroll (`web_app/app/ui/sections/classification.py`).
- ✅ Se incorporó un panel de vista previa sticky con métricas vivas y chips de ejemplos por grupo, sustituyendo los expanders por tabs directas para revisar los dataframes aplicados (`web_app/app/ui/sections/classification.py`).
