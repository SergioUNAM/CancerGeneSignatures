# Seguimiento de Progreso — Web App CancerGeneSignatures

Fecha: 2025-09-10
Ámbito: `web_app/` (Streamlit) + `src/core/` (lógica qPCR/FC)

**Resumen**
- Objetivo: Analizar datos qPCR desde Excel, calcular ΔΔCt y Fold Change (promedios vs gen de referencia), y visualizar resultados.
- Estado: MVP funcional con carga Excel, clasificación por prefijos, cálculo FC, y visualizaciones Plotly.

**Arquitectura Breve**
- UI: `web_app/streamlit_app.py`
- Config: `web_app/config/menu.json`
- IO/Parsing: `src/core/io.py`
- qPCR utils: `src/core/qpcr.py`
- Limpieza: `src/core/cleaning.py`
- Fold Change: `src/core/fold_change.py`
- Tablas/Gráficas: `src/core/tables.py`

**Cómo Ejecutar**
- Crear entorno: `pip install -r web_app/requirements.txt`
- Iniciar app: `streamlit run web_app/streamlit_app.py`
- Menú alternativo (opcional): exportar `CGS_MENU_PATH` a un JSON compatible.

**Hitos y Estado**
- [x] Carga Excel y selección de hoja (`web_app/streamlit_app.py:118`)
- [x] Vista previa y metadatos extraídos (`web_app/streamlit_app.py:151`)
- [x] Ancho→largo y filtro básicos de controles de máquina (`web_app/streamlit_app.py:166`, `src/core/cleaning.py:6`)
- [x] Clasificación por prefijos y estado persistente (`web_app/streamlit_app.py:202`)
- [x] Imputación NaN→máximo Ct global (`web_app/streamlit_app.py:254`)
- [x] Cálculo FC (promedios y gen de referencia) (`src/core/fold_change.py:22`)
- [x] Visualizaciones comparativas y tabla (`web_app/streamlit_app.py:283`, `src/core/tables.py:8`)
- [x] Menú configurable con cache y fallback (`web_app/streamlit_app.py:88`, `web_app/config/menu.json`)
- [x] Botones de descarga para CSV generados
- [x] Mejora clasificación (case-insensitive)
- [ ] Validaciones y mensajes más específicos
- [ ] Documentación de plantillas Excel soportadas
- [ ] Soporte opcional `.xls` (o guía para convertir a `.xlsx`)

**Cambios Recientes (Changelog)**
- 2025-09-10
  - Añadidos botones de descarga (`st.download_button`) para los CSV generados: controles/muestras limpios, consolidado FC, y expresión categorizada.
  - Clasificación de controles/muestras ahora case-insensitive usando `classify_tests` de `src/core/qpcr.py`.
  - Documentación inicial de seguimiento en `web_app/PROGRESO.md`.
- 2025-08-31
  - MVP de Streamlit con:
    - Carga Excel y selección de hoja (`list_excel_sheets`, `parse_qpcr_wide`).
    - Clasificación por prefijos manuales y estado por archivo.
    - Cálculo Fold Change por promedios y por gen de referencia (menor desviación estándar).
    - Gráfica comparativa ΔΔCt y FC, y tabla Plotly.
    - Menú externo `menu.json` con tipos de cáncer/contextos/método preferido.

**Decisiones Técnicas**
- Gen de referencia: seleccionado por mínima media de desviaciones estándar entre grupos (`src/core/fold_change.py:44`).
- Imputación: valores no determinados/NaN se convierten a máximo Ct global antes de FC (`web_app/streamlit_app.py:254-263`).
- Filtro controles de máquina: en UI se pasan PPC/RTC; función soporta lista por defecto más amplia (`src/core/cleaning.py:18`).
- Parsing Excel robusto (cabeceras y fila de nombres de prueba) (`src/core/io.py:76`).

**Novedades de Uso**
- Descargas: en la sección de resultados aparecen botones para descargar CSVs:
  - `controles_limpios.csv`, `muestras_limpias.csv`
  - `fold_change_consolidado.csv`, `expresion_categorizada.csv`
- Clasificación case-insensitive: los prefijos de controles y muestras se comparan sin distinción de mayúsculas/minúsculas.

**Riesgos y Pendientes**
- Posible omisión del primer nombre de prueba en resumen (`web_app/streamlit_app.py:190`).
- `idxmin` puede fallar si todas las std son NaN; faltan checks previos.
- `.xls` sin `xlrd`; README sugiere convertir a `.xlsx` (aceptable, documentar mejor).

**Próximos Pasos (Prioridad sugerida)**
- Alta: Validaciones previas a FC y mensajes claros si grupos vacíos o sin solapamiento de targets.
- Media: Integrar sugerencias de prefijos/sufijos con `suggest_name_affixes`.
- Media: Unificar filtros de controles de máquina usando la lista por defecto + extras (PPC/RTC).
- Baja: Documentar formato esperado de plantillas y políticas para “Undetermined/ND”.
- Baja: Soporte opcional `.xls` (añadir `xlrd`) o guía de conversión.

**Convenciones para Actualizar Este Documento**
- Añadir entradas al Changelog con fecha ISO y bullets concisos.
- Mantener “Hitos y Estado” actualizado (marcar con [x]/[ ] según progreso).
- Referenciar archivos y líneas cuando aplique, formato `ruta:línea`.

Plantilla para una nueva entrada:
- YYYY-MM-DD
  - Cambio 1 breve.
  - Cambio 2 breve.
  - Impacto/nota si aplica.

