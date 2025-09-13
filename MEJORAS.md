# Plan de Mejora por Fases — CancerGeneSignatures

Este documento organiza un backlog de mejoras para la app y los módulos núcleo, priorizado y agrupado por fases incrementales. Cada fase incluye objetivo, tareas, criterios de aceptación, referencias y notas.

Estado inicial evaluado sobre `web_app/streamlit_app.py` y `src/core/*` a fecha 2025-09-11.

---

## Fase 0 — Hotfixes y Bugs (alta prioridad)

Objetivo: Corregir fallos evidentes que afectan la veracidad de resultados o UX inmediata.

Tareas
- [x] Resumen de pruebas omite el primer nombre: usar todos los `sample_names`.
  - Referencia: `web_app/streamlit_app.py:279`
- [x] Robustecer elección del gen de referencia cuando no hay datos suficientes.
  - Referencia: `src/core/fold_change.py:61`, `src/core/fold_change.py:64`, `src/core/fold_change.py:65`
  - Criterios: si toda la estabilidad es NaN o no hay intersección de genes entre grupos, mostrar mensaje claro en la UI y abortar el cálculo en lugar de fallar.
- [x] Evitar SettingWithCopy en imputación de Ct; operar sobre copias seguras o usar `.loc`.
  - Referencia: `web_app/streamlit_app.py:384`, `web_app/streamlit_app.py:385`
- [x] Unificar filtros de controles de máquina: aplicar lista por defecto + extras (PPC/RTC) en vez de reemplazarla.
  - Referencias: `web_app/streamlit_app.py:258`, `src/core/cleaning.py:25`

Criterios de aceptación
- No hay errores al cargar y procesar archivos con grupos sin solapamiento o con muchas NaN.
- El resumen lista todas las pruebas detectadas.
- La imputación no genera warnings de pandas.
- Los controles de máquina estándar y PPC/RTC se filtran correctamente por defecto.

Notas
- Mantener mensajes de error/aviso en español, breves y accionables.

---

## Fase 1 — Validaciones y Robustez de Flujo

Objetivo: Prevenir rutas inválidas con validaciones tempranas y configurables.

Tareas
- [x] Pre-checks antes de FC: grupos vacíos, intersección de `target`, ratio de NaN por grupo, n mínimo por gen/grupo.
  - Mostrar panel de calidad con métricas y avisos.
- [x] Exponer en la UI la política de “Undetermined/ND” (nan | ctmax | value) y valor por defecto.
  - Conectar con `parse_qpcr_wide(..., undetermined_policy=...)` en lugar de imputar después.
  - Referencia: `src/core/io.py:79`, `src/core/io.py:181`
- [x] Ampliar tokens de “Undetermined”: incluir variantes como `na`, `n.a.`, `n/d`, `und.` (normalizar minúsculas y espacios).
- [x] Detección y resolución de nombres de muestra duplicados (renombrar con sufijos incrementales).

Criterios de aceptación
- Si una validación falla, la UI muestra causa y sugerencia (sin traceback).
- El usuario puede cambiar la política “Undetermined” desde la barra lateral y ver su efecto.

Notas
- Registrar en logs por qué se tomó cada decisión (e.g., fila de nombres elegida, política de ND aplicada).

---

## Fase 2 — UX de Clasificación y Navegación

Objetivo: Hacer la clasificación de controles/muestras más flexible y separar pasos pesados en páginas.

Tareas
- [x] Clasificación por prefijos, sufijos y regex (apoyarse en `classify_by_prefixes`/`classify_by_suffixes`).
  - Mostrar vista previa de cuántas y cuáles pruebas caen en cada bucket antes de aplicar.
  - Referencia: `src/core/qpcr.py:84`, `src/core/qpcr.py:92`
- [x] Detectar y advertir colisiones (mismas pruebas en ambos grupos) con opción de auto-resolver.
  - Referencia: `web_app/streamlit_app.py:349`
- [x] Mover secciones pesadas a `web_app/pages/` (STRING, PubMed, Firmas) para un flujo por pasos.
- [x] Gráficos FC: toggle escala lineal/log, resaltar el gen de referencia, tooltips claros.

Criterios de aceptación
- La pestaña de clasificación permite elegir por prefijo/sufijo/regex y previsualiza resultados.
- La app arranca en “Análisis qPCR” y tiene páginas separadas para “Enriquecimiento”, “Bibliografía” y “Firmas”.

Notas
- Mantener estado entre páginas usando `st.session_state` con claves namespaced.

---

## Fase 3 — Rendimiento, Caché y APIs

Objetivo: Reducir latencia y llamadas redundantes a servicios externos; parametrizar endpoints.

Tareas
- [ ] CACHE: resultados de Ensembl/STRING/PubMed con claves deterministas y TTL (p. ej., `requests_cache` o `@st.cache_data(ttl=...)`).
  - Referencias: `web_app/streamlit_app.py:456`, `web_app/streamlit_app.py:463`, `web_app/streamlit_app.py:470`, `web_app/streamlit_app.py:748`
- [ ] Ensembl: limitar `max_workers` por UI y aplicar backoff ante rate limiting; métrica de acierto (IDs/descr. encontrados).
  - Referencia: `src/core/ensembl.py:1`
- [ ] STRING: parametrizar `species` (lista común: 9606, 10090, 10116) y `sources`; base URL via env `CGS_STRING_URL`.
  - Referencias: `web_app/streamlit_app.py:578`, `web_app/streamlit_app.py:582`, `web_app/streamlit_app.py:593`, `src/core/string_enrichment.py:11`
- [ ] PubMed: usar búsqueda en `[TIAB]`, paginar con `retstart`, y respetar guías NCBI (throttling con/sin API key).
  - Referencias: `src/core/bibliography.py:114`
- [ ] Credenciales: migrar a `st.secrets` y eliminar escritura a `os.environ` en runtime.
  - Referencias: `web_app/streamlit_app.py:748`, `web_app/streamlit_app.py:750`

Criterios de aceptación
- Re-ejecutar una consulta reciente usa caché; tiempos de respuesta se reducen notablemente.
- No se observan errores por rate limiting; hay mensajes de progreso y métricas de resultados.

Notas
- Añadir variables de entorno documentadas (`CGS_STRING_URL`, `CGS_LOGLEVEL`).

---

## Fase 4 — Parsing qPCR y Formatos

Objetivo: Soportar variaciones de plantillas y documentar claramente los requisitos.

Tareas
- [ ] Auto-detector de cabeceras más robusto: aceptar “Ct/CT/Ct Mean”, “Well Position” y normalizar unicode.
  - Referencia: `src/core/io.py:104`
- [ ] Documentar y exponer override manual de fila/columnas de cabecera (UI avanzada).
- [ ] Soporte opcional `.xls` o guía clara para convertir a `.xlsx`.
- [ ] Ampliar heurística de nombres de prueba y registrar decisión (por qué se eligió esa fila).

Criterios de aceptación
- Casos comunes de plantillas alternativas procesan sin errores o con instrucciones claras para el usuario.

Notas
- Mantener compatibilidad con la firma de `parse_qpcr_wide` y fallback a `header_mode="auto"`.

---

## Fase 5 — Firmas Genéticas y Enriquecimiento Hallmarks

Objetivo: Hacer configurables los parámetros de construcción de firmas y mejorar visualización.

Tareas
- [ ] UI: sliders para ponderar `recency` y `freq_ratio`, mínimo de artículos por gen, filtro por contexto.
  - Referencia: `src/core/signatures.py:87`
- [ ] Hallmarks (gseapy): validación de rutas GMT, subida de ficheros y mensajes claros cuando faltan.
- [ ] Visualizaciones: sunburst/heatmap de scores y -log10(p); enlaces a PMIDs desde tooltips.

Criterios de aceptación
- El usuario puede ajustar parámetros y observar cambios en la tabla/gráficos de firmas.
- Descargas CSV/XLSX contienen columnas de p-valor por hallmark y mapping gene→PMIDs.

Notas
- Mantener `gseapy` como dependencia opcional; mensajes si no está instalado.

---

## Fase 6 — Calidad, Tests y Estilo

Objetivo: Mejorar confiabilidad con pruebas automáticas y estilo consistente.

Tareas
- [ ] Unit tests para: `parse_qpcr_wide` (cabeceras, undetermined), `compute_fold_change` (casos borde), `string_enrichment` y `bibliography` (mocks de red).
- [ ] Lint/format con `ruff` y `black`; pre-commit.
- [ ] Empaquetado editable (`pip install -e .`) para `src/` y evitar `sys.path` en la app.

Criterios de aceptación
- Test suite pasa en CI y local; sin flake/ruff críticos.

Notas
- No ampliar alcance de tests a notebooks por ahora.

---

## Fase 7 — Documentación y Ejemplos

Objetivo: Alinear documentación con la app real y facilitar onboarding.

Tareas
- [ ] Actualizar `README.md` (root) para reflejar módulos reales y flujo web.
- [ ] Ampliar `web_app/README.md` con troubleshooting (cabeceras, dependencias), límites de tiempo y buenas prácticas de PubMed/Ensembl.
- [ ] Añadir un Excel de ejemplo anonimizado en `raw_data/` + capturas del flujo.

Criterios de aceptación
- Cualquier persona puede ejecutar la app y reproducir el flujo con el dataset de ejemplo.

---

## Fase 8 — Despliegue

Objetivo: Facilitar ejecución reproducible y segura.

Tareas
- [ ] Dockerfile con `web_app/requirements.txt` y `streamlit run web_app/streamlit_app.py`.
- [ ] Variables seguras vía `st.secrets` y/o `.env` (no commits).
- [ ] Config de Streamlit (recursos, ancho, límites de subida) y healthcheck.

Criterios de aceptación
- La app levanta en Docker localmente; imágenes reproducibles.

---

## Apéndice — Lista de mejoras rápidas (checklist)

- [x] Corregir omisión de primer test en resumen (`web_app/streamlit_app.py:279`).
- [x] Guardas en FC para `NaN` y solapamiento (`src/core/fold_change.py:61`, `:64`, `:65`).
- [x] Unificar filtro de controles de máquina (usar defaults + PPC/RTC) (`web_app/streamlit_app.py:258`, `src/core/cleaning.py:25`).
- [x] Imputación segura de Ct (sin SettingWithCopy) (`web_app/streamlit_app.py:384`, `:385`).
- [x] UI para política “Undetermined” y tokens extendidos (`src/core/io.py:79`, `:181`).
- [ ] Clasificar por sufijos/regex con previsualización (`src/core/qpcr.py:84`, `:92`).
- [ ] Pages de Streamlit para modularizar flujo.
- [ ] Caché a disco/TTL para Ensembl/STRING/PubMed (`web_app/streamlit_app.py:456`, `:463`, `:470`, `:748`).
- [ ] STRING configurable (species/sources/base URL) (`web_app/streamlit_app.py:578`, `:582`, `:593`, `src/core/string_enrichment.py:11`).
- [ ] PubMed: TIAB + paginación + throttling (`src/core/bibliography.py:114`).
- [ ] Secrets en `st.secrets` en vez de `os.environ` (`web_app/streamlit_app.py:748`, `:750`).
- [ ] Tests, linting, empaquetado editable y README actualizados.

---

## Estimación y Orden de Ataque sugerido

1) Fase 0 (0.5–1 día) → 2) Fase 1 (0.5–1.5 días) → 3) Fase 2 (1–2 días) → 4) Fase 3 (1–2 días) → 5) Fase 4 (1 día) → 6) Fase 5 (1–2 días) → 7) Fase 6 (1–2 días) → 8) Fase 7 (0.5–1 día) → 9) Fase 8 (0.5–1 día)

Las fases 0–3 cubren la mayor parte del valor para un MVP robusto.

---

## Convenciones

- Mantener referencias de archivo con `ruta:línea` para trazabilidad.
- Mensajes en español, concisos y con acción sugerida.
- No introducir dependencias nuevas sin justificar en la fase correspondiente.
