# Seguimiento de Progreso — Web App CancerGeneSignatures

Fecha: 2025-09-10
Ámbito: `web_app/` (Streamlit) + `app/core/` (lógica qPCR/FC)

**Resumen**
- Objetivo: Analizar datos qPCR desde Excel, calcular ΔΔCt y Fold Change (promedios vs gen de referencia), y visualizar resultados.
- Estado: MVP funcional con carga Excel, clasificación por prefijos, cálculo FC, y visualizaciones Plotly.

**Arquitectura Breve**
- UI: `web_app/streamlit_app.py`
- Config: `web_app/config/menu.json`
- IO/Parsing: `app/core/io.py`
- qPCR utils: `app/core/qpcr.py`
- Limpieza: `app/core/cleaning.py`
- Fold Change: `app/core/fold_change.py`
- Tablas/Gráficas: `app/core/tables.py`

**Cómo Ejecutar**
- Crear entorno: `pip install -r web_app/requirements.txt`
- Iniciar app: `streamlit run web_app/streamlit_app.py`
- Menú alternativo (opcional): exportar `CGS_MENU_PATH` a un JSON compatible.

**Hitos y Estado**
- [x] Carga Excel y selección de hoja (`web_app/streamlit_app.py:118`)
- [x] Vista previa y metadatos extraídos (`web_app/streamlit_app.py:151`)
- [x] Ancho→largo y filtro básicos de controles de máquina (`web_app/streamlit_app.py:166`, `app/core/cleaning.py:6`)
- [x] Clasificación por prefijos y estado persistente (`web_app/streamlit_app.py:202`)
- [x] Imputación NaN→máximo Ct global (`web_app/streamlit_app.py:254`)
- [x] Cálculo FC (promedios y gen de referencia) (`app/core/fold_change.py:22`)
- [x] Visualizaciones comparativas y tabla (`web_app/streamlit_app.py:283`, `app/core/tables.py:8`)
- [x] Menú configurable con cache y fallback (`web_app/streamlit_app.py:88`, `web_app/config/menu.json`)
- [x] Botones de descarga para CSV generados
- [x] Mejora clasificación (case-insensitive)
- [x] Enriquecimiento STRING (UI + API, filtros GO/KEGG, descarga)
- [ ] Validaciones y mensajes más específicos
- [ ] Documentación de plantillas Excel soportadas
- [ ] Soporte opcional `.xls` (o guía para convertir a `.xlsx`)

**Cambios Recientes (Changelog)**
- 2025-09-12
  - Se introdujo `web_app/app_state.py` para centralizar preferencias de sesión (política ND, contexto, bandera de genes estables).
  - La barra lateral de Streamlit ahora recuerda selecciones previas usando el contenedor tipado en lugar de claves sueltas.
- 2025-09-10
  - Merge a `master` de `feature/webapp-ensembl-integration` (Ensembl interactivo, UX de clasificación mejorada, parser A4/B4 + heurística).
- 2025-09-10
- Parser Excel por coordenadas: la app intenta primero encabezado fijo en A4/B4 (fila 4) y cae a detección automática si falla.
  - Compatibilidad: si el entorno aún carga una versión previa del parser sin nuevos parámetros, la app detecta `TypeError` y reintenta con la firma antigua.
  - Integración Ensembl en la UI: anotación de `target` con `ensembl_id` y `description`, tabla y descarga `ensembl_anotado.csv`. Se añadió `requests` a `web_app/requirements.txt`.
  - Mejora de anotación Ensembl: fallback para descripción usando `lookup/id/{id}`, `lookup/symbol`, y `mygene.info` si es necesario; cache por símbolo.
  - UX de clasificación mejorada: pestañas "Por prefijos" (con sugerencias automáticas) y "Selección manual" (multiselect de pruebas) con persistencia en sesión.
- 2025-09-10
  - Merge a `master` de la rama `feature/webapp-downloads-caseinsensitive` (descargas CSV, clasificación case-insensitive, fix import `src`).
- 2025-09-10
  - Añadidos botones de descarga (`st.download_button`) para los CSV generados: controles/muestras limpios, consolidado FC, y expresión categorizada.
  - Clasificación de controles/muestras ahora case-insensitive usando `classify_tests` de `app/core/qpcr.py`.
  - Documentación inicial de seguimiento en `web_app/PROGRESO.md`.
  - Corrección de importaciones de `src`: se añadieron `src/__init__.py` y `app/core/__init__.py`, y se inyectó la raíz del proyecto en `sys.path` en `web_app/streamlit_app.py` para permitir `from app.core...` al ejecutar desde `web_app/`.
 - 2025-09-10
 - Enriquecimiento STRING integrado en la UI: selección por niveles de expresión, fuentes (GO/KEGG/Reactome), filtros por FDR/tamaño/top-N, tabla, barras (-log10 FDR) y descarga CSV.
 - Nuevo módulo `app/core/string_enrichment.py` con helpers `run_string_enrichment` y `filter_enrichment` (sin dependencias adicionales, usa `requests`).
- 2025-09-10
  - Bibliografía (PubMed) integrada en la web app:
    - Nuevo módulo `app/core/bibliography.py` con `search_pubmed_by_genes`, `classify_bibliography`, `aggregate_counts_by_level_and_cancer`.
    - UI con inputs para `NCBI Email` (obligatorio) y `NCBI API Key` (opcional) directamente en la app (demo-friendly, sin necesidad de exportar variables).
    - Progreso visible por gen (barra + estado) y logs de avance/errores.
    - Descargas CSV de resultados y de bibliografía clasificada; gráficos por tipo de cáncer y nivel.
    - Dependencia agregada: `biopython` (para Entrez/Medline).
  - Logging y trazabilidad:
    - Logs informativos en puntos largos: anotación Ensembl, enriquecimiento STRING y búsqueda PubMed (inicio/fin, conteos). Se controla por `CGS_LOGLEVEL` (INFO por defecto).
 - 2025-09-11
  - Firmas genéticas (MVP) en la web app:
    - Nuevo módulo `app/core/signatures.py` para construir firmas desde bibliografía clasificada y enriquecimiento Hallmarks.
    - Sección "Firmas genéticas" en la app: selección de GMTs, generación de DataFrame y descarga CSV.
    - Visualización básica con sunburst por gen→hallmark y color por -log10(p), separado por nivel de expresión.
    - Requiere `gseapy` (añadido a requirements).

**Ramas**
- `master`: actualizado con descargas CSV, clasificación case-insensitive y fix de importaciones.
- `feature/webapp-ensembl-integration`: fusionada en `master`.
- `feature/webapp-string-gsea`: nueva rama para integrar enriquecimiento STRING/GSEA y visualización GO.
- 2025-08-31
  - MVP de Streamlit con:
    - Carga Excel y selección de hoja (`list_excel_sheets`, `parse_qpcr_wide`).
    - Clasificación por prefijos manuales y estado por archivo.
    - Cálculo Fold Change por promedios y por gen de referencia (menor desviación estándar).
    - Gráfica comparativa ΔΔCt y FC, y tabla Plotly.
    - Menú externo `menu.json` con tipos de cáncer/contextos/método preferido.

**Decisiones Técnicas**
- Gen de referencia: seleccionado por mínima media de desviaciones estándar entre grupos (`app/core/fold_change.py:44`).
- Imputación: valores no determinados/NaN se convierten a máximo Ct global antes de FC (`web_app/streamlit_app.py:254-263`).
- Filtro controles de máquina: en UI se pasan PPC/RTC; función soporta lista por defecto más amplia (`app/core/cleaning.py:18`).
- Parsing Excel robusto (cabeceras y fila de nombres de prueba) (`app/core/io.py:76`).
- Encabezados por coordenadas: `parse_qpcr_wide(..., header_mode="coords", header_row_idx=3, well_col_idx=0, target_col_idx=1)` con fallback a `header_mode="auto"`.
- Hot-reload en UI para `app.core.io`: se recarga el módulo al iniciar la app para asegurar la firma más reciente del parser.
 - Heurística de fila de nombres de prueba: se prefieren filas con mayor número de códigos con guion (p.ej., `4GB-001`) y dígitos, para escoger la fila 2 cuando aplica (sobre etiquetas como `C1`, `DGC1`).

**Novedades de Uso**
- Descargas: en la sección de resultados aparecen botones para descargar CSVs:
  - `controles_limpios.csv`, `muestras_limpias.csv`
  - `fold_change_consolidado.csv`, `expresion_categorizada.csv`
- Clasificación case-insensitive: los prefijos de controles y muestras se comparan sin distinción de mayúsculas/minúsculas.
- Ensembl: nueva sección anota genes con `ensembl_id` y `description`; incluye descarga `ensembl_anotado.csv` (requiere conexión a internet).
  - Nota: si falta descripción en la primera consulta, se intenta completarla con métodos alternativos.
  - Explorar interactivo: pestañas Resumen (métricas + gráfico), Explorar (filtros por gen/descr./nivel + descarga filtrada) y Enlaces (accesos directos a Ensembl).

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

**Cobertura del notebook `notebooks/uimeo_data_analysis.py`**
- Carga Excel y extracción básica: [x]
  - Notebook: `cargar_datos_excel`, `extraer_informacion_cg`
  - App: `app/core/io.parse_qpcr_wide` + vista previa (`web_app/streamlit_app.py:151`), extracción de tests/genes/pozos (`web_app/streamlit_app.py:176`)
- Clasificación por prefijos (controles/muestras): [x]
  - Notebook: `procesar_ct` (case-insensitive)
  - App: `app/core/qpcr.classify_tests` (`web_app/streamlit_app.py:226`)
- Filtrado de controles de máquina (PPC/RTC): [x]
  - Notebook: `filtrar_controles_maquina` (PPC/RTC)
  - App: `app/core/cleaning.drop_machine_controls` (`web_app/streamlit_app.py:169`)
- Imputación de Ct (NaN → máximo Ct global): [x]
  - Notebook: `procesar_ct_column` (valor máximo)
  - App: bloque de imputación (`web_app/streamlit_app.py:254`)
- Cálculo ΔΔCt y Fold Change (promedios y gen ref): [x]
  - Notebook: `calcular_fold_change`
  - App: `app/core/fold_change.compute_fold_change` (`web_app/streamlit_app.py:268`)
- Categorización de expresión (sub/estable/sobre): [x]
  - Notebook: `categorizar_expresion`
  - App: `pd.cut` sobre FC escogido (`web_app/streamlit_app.py:300`)
- Visualización básica (tabla/plots comparación FC): [x]
  - Notebook: tablas y gráficos de resumen
  - App: `app/core/tables.fc_comparison_table` + barras/series Plotly (`web_app/streamlit_app.py:279`, `:283`)
- Descarga de resultados: [x] (CSV en app; Excel en notebook)
  - Notebook: `export_dfs_to_excel`
  - App: `st.download_button` para CSV (`web_app/streamlit_app.py:329`)
- Enriquecimiento Ensembl (IDs/descr.): [x]
  - Notebook: `add_ensembl_info_batch` (requests + fallback)
  - App: Integrado con `app/core/ensembl.add_ensembl_info_batch` (tabla + descarga)
- Enriquecimiento STRING y/o GSEA (gseapy): [ ]
  - Notebook: llamadas a STRING/gseapy, preparación y gráficos
  - App: no integrado aún (requiere endpoints/red y UI)
- Filtrado y visualización GO: [ ]
  - Notebook: `filtrar_enriquecimiento_terminos_GO` y gráficos
  - App: no integrado
- Redes/interacción (networkx/cytoscape): [ ]
  - Notebook: construcción de redes y exportaciones
  - App: no integrado

Siguiente foco propuesto (por orden):
- Visualización GO detallada (filtros por BP/MF/CC con tooltips, ranking por NES/OR) y enlaces.
- Explorar integración GSEA (gseapy) usando conjuntos locales en `gen-sets_GSEA_MSigDB/`.

**Notas — Ensembl (UI)**
- Implementado tras la categorización en `web_app/streamlit_app.py`.
- Consideraciones: requiere red; se captura error y se notifica al usuario.
