# Filosofía de Ingeniería de Software — CancerGeneSignatures

Este documento establece las prácticas de ingeniería que guiarán la evolución de la web app. Su objetivo es asegurar un producto mantenible, verificable y fácil de extender a medida que reincorporemos módulos avanzados.

## 1. Principios de diseño
- **Modularidad explícita**: `app.core` provee lógica pura, `app.services` orquesta dominios y la capa Streamlit se limita a presentación. Cualquier nueva característica debe respetar esta separación.
- **Dependencias opcionales desacopladas**: integraciones externas (STRING, PubMed, NLP…) deben cargarse de forma perezosa (`try/except` o wrappers) y declararse como extras en `requirements` cuando se reactiven.
- **Single source of truth para configuraciones**: parámetros viven en `config/menu.json` o variables de entorno (`CGS_*`). Evitar valores mágicos en el código UI.
- **Persistencia mínima de estado**: `st.session_state` solo almacena artefactos imprescindibles (datos cargados, selección de usuario). Limpiar entradas auxiliares tras cada flujo.

## 2. Calidad y validación
- **Pruebas unitarias** para funciones puras (e.g., `parse_qpcr_wide`, imputación, selección de referencias). Las pruebas deben residir en `tests/` y usar datos sintéticos.
- **Pruebas de integración ligeras**: ejercen servicios orquestadores con fixtures reducidas (fold change, normalización). Evitar tests que dependan de red.
- **Lint + formato**: `ruff` (lint) y `black` (formato) en modo estricto antes de cada PR. Configuración compartida en `pyproject.toml`.
- **Revisión por pares**: todo cambio significativo entra via PR, acompañado de descripción, capturas o enlaces a runs de Streamlit cuando aplique.

## 3. Observabilidad y errores
- **Logging estructurado**: usar `logging.getLogger("cgs.<modulo>")`, mensajes concisos en español y niveles apropiados (`info`, `warning`, `error`). NUNCA silenciar excepciones sin registrar contexto.
- **Mensajes para usuarios**: en UI preferir `st.info/warning/error` con instrucciones accionables; evitar mostrar stack traces.
- **Feature flags**: al reintroducir módulos, rodear comportamiento inestable con toggles (`CGS_ENABLE_*` o configuración en menú) para facilitar rollback.

## 4. Datos y entrada
- **Validaciones tempranas**: toda carga de Excel, CSV o JSON debe pasar por funciones validadoras antes de su uso aguas abajo. Registrar razones de aceptación/rescate.
- **Inmutabilidad por defecto**: retornar copias (`.copy()`) al exponer DataFrames para impedir efectos colaterales. Mutaciones explícitas deben documentarse.
- **Compatibilidad hacia atrás**: si se modifican estructuras de datos, ofrecer adaptadores o scripts de migración.

## 5. Seguridad y secretos
- **Credenciales fuera del código**: utilizar `st.secrets` o variables de entorno. Prohibido subir claves a repositorio.
- **Manejo de errores de red**: envolver requests en timeouts y reintentos limitados. Nunca propagar respuestas crudas sin sanitizar.

## 6. Entrega continua
- **CI/CD**: pipeline mínimo con pasos `lint → tests → build (streamlit smoke)` antes de desplegar.
- **Versionado semántico**: etiquetar releases (`vX.Y.Z`) y mantener changelog (`docs/registro_cambios.md`).
- **Backups de configuración**: cambios en `config/menu.json` deben acompañarse de validación automática (`app.config.models.MenuConfig`).

## 7. Documentación viva
- **README funcional** (web_app) siempre acorde al alcance real.
- **Guías específicas** por módulo cuando el comportamiento no sea obvio (e.g., explicación de heurísticas de nombre de muestra).
- **Notebook ≠ especificación**: decisiones tomadas en notebooks deben migrar a documentación o código productivo antes de incorporarse.

## 8. Roadmap y deuda técnica
- Mantener una lista priorizada (e.g., `MEJORAS.md`) con tareas pendientes y responsables.
- Registrar decisiones clave (ADR ligero) en `docs/` cuando impacten arquitectura o dependencias.
- Revisar deuda en cada iteración; si no se aborda una tarea crítica, documentar el impacto.

Este README actúa como contrato de ingeniería: cualquier excepción debe justificarse y planificarse para su regularización en el backlog.

---

## Anexo A — Diseño de la Normalización Avanzada (qPCR)

Objetivo: por exploración exhaustiva (fuerza bruta) elegir el conjunto de genes de referencia que mejor actúa como regla común (estable intra-grupo y sin sesgo entre grupos), normalizar Ct a expresión relativa (log2), y evaluar separación Control vs Muestra con evidencia estadística robusta. El set elegido se reutiliza para fold change posterior.

Implementación actual (producida en `app.core.reference_normalization` y orquestada desde `app.services.normalization`):

- Candidatos y búsqueda K (1–3)
  - Se calculan candidatos por baja desviación estándar de Ct por gen.
  - Se prueban combinaciones de tamaño K y se puntúan con `stability_score_for_refs` (promedio ponderado: variabilidad intra-grupo + diferencia inter-grupo sobre el Ct de referencia medio por test).
- Normalización con referencias elegidas
  - Para cada test se promedia el Ct de las referencias → `ref_ct`; se computa `ΔCt = Ct - ref_ct` y `log2_rel_expr = -ΔCt`.
- Evaluación de evidencia y robustez
  - Estadística diferencial por gen (t-test, p y q Benjamini–Hochberg).
  - Frecuencia bootstrap de significancia (re-muestreo de tests).
  - Tasa empírica de falsos positivos por permutación de etiquetas Control/Muestra.
- Salidas
  - `df_norm` (tabla normalizada), `df_heatmap` (genes×tests), `reference_result` (refs elegidas + ranking), `differential` (tabla con p, q, efectos, bootstrap, FPR).
  - En la UI se ofrecen heatmaps con dendrogramas y un comparativo con métodos básicos (promedio global, gen de referencia único) para todo el panel o subconjuntos.

Estado de adecuación al objetivo: CUMPLE. La versión actual explora combinaciones, elige el set más estable y aporta métricas de robustez. Se integra con el cálculo de fold change y visualizaciones comparativas.

Mejoras recomendadas (no bloqueantes):
- Selección con criterio compuesto: además de estabilidad de referencias, incluir una métrica de separación post‑normalización (p.ej., mediana de |log2FC| o varianza explicada entre grupos) con validación cruzada para evitar sobreajuste.
- Penalizar referencias con alta fracción imputada o cobertura incompleta por test.
- Permitir opción de mediana (en vez de media) al promediar Ct de referencias o al promedio global por test.
- Umbrales de elegibilidad: replicados mínimos por grupo, |Δmedia(Ct)| entre grupos < 0.5–1.0 Ct para considerar un gen como referencia candidata.
- Afinar supuestos del método de promedios: calcular la línea base sobre genes comunes entre grupos y permitir seleccionar media/mediana.

Trazabilidad de módulos:
- Selección/normalización/estadística: `web_app/app/core/reference_normalization.py`
- Orquestación y matrices comparativas: `web_app/app/services/normalization.py`
- Heatmaps con dendrogramas: `web_app/app/services/heatmap_visuals.py`
- Comparativos ΔΔCT/FC (métodos básicos): `web_app/app/core/fold_change.py`, `web_app/app/services/fold_change_visuals.py`
