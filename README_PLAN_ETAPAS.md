# Roadmap de Mejora por Etapas — Mejores Prácticas de Ingeniería

Guía para la adopción progresiva de mejores prácticas en **CancerGeneSignatures**. Cada etapa define objetivos, acciones clave, entregables y criterios de salida. El plan se centra en robustecer calidad, mantenibilidad, observabilidad y operaciones sin introducir nuevas funciones de negocio.

## Resumen de estado
| Etapa | Objetivo principal | Estado | Comentario breve |
|-------|-------------------|--------|------------------|
| 0 | Diagnóstico y bases compartidas | ✅ Cerrada | Inventario, estándares y prioridades consensuadas. |
| 1 | Modularización y configuración | 🔄 En curso | UI multipágina, plantillas y documentación final. |
| 1.5 | Preparativos de calidad/pruebas | 🔄 En curso | Carpeta `tests/` + workbook sintético, pendiente ampliar escenarios. |
| 2 | Calidad de código y seguridad | ⏳ Pendiente | Pre-commit, auditoría de dependencias y políticas. |
| 3 | Estrategia de pruebas | ⏳ Pendiente | Pruebas unitarias/contrato y smoke UI. |
| 4 | Observabilidad y manejo de errores | ⏳ Pendiente | Trazabilidad, métricas y alertas. |
| 5 | Automatización CI/CD | ⏳ Pendiente | Pipelines de build/test/deploy consistentes. |
| 6 | Gobernanza y onboarding | ⏳ Pendiente | Documentación final, procesos de incorporación. |

---

## Etapa 0 — Diagnóstico y Bases Compartidas *(Cerrada)*
**Objetivo**: Alinear al equipo con el estado actual del código y acordar estándares mínimos.

**Logros clave**
- Inventario de módulos críticos (`src/core`, `web_app`) y dependencias externas.
- Documentación de decisiones arquitectónicas y deudas técnicas visibles.
- Definición de estilo y tipado (PEP 8 + hints con `mypy`/`pyright`).
- Tablero Kanban con backlog priorizado de mejoras.

**Entregables**
- Informe de diagnóstico y prioridades.
- Documento de estándares de código y convenciones de logs/errores.

**Criterios de salida alcanzados**
- Equipo alineado en alcance y limitaciones.
- Glosario común listo para etapas posteriores.

---

## Etapa 1 — Modularización y Gestión de Configuración *(En progreso)*
**Objetivo**: Reducir acoplamiento y centralizar configuraciones para facilitar pruebas y despliegues.

**Avances**
- `web_app/streamlit_app.py` reducido a bootstrapping; lógica en `app/services` y `app/ui`.
- Configuración tipada (`AppConfig`, `ServicesConfig`) con lectura resiliente.
- Servicios especializados (qPCR, fold change, STRING, bibliografía, heurística, NLP).
- Secciones UI reutilizables en `app/ui/sections.py`.

**Acciones restantes prioritarias**
- Crear `app/ui/pages/*` para desacoplar qPCR, enriquecimiento, firmas y reportes.
- Añadir adaptadores ligeros para nuevas fuentes (cacheo Ensembl/STRING).
- Publicar plantillas de configuración por entorno (`dev`, `staging`, `prod`).

**Entregables pendientes**
- README principal actualizado (ver `README.md`).
- Plantillas de configuración + documentación de variables ambiente.

**Criterios de salida actualizados**
- Cada flujo de la web app invocado desde secciones o páginas dedicadas.
- Configuraciones/credenciales resueltas exclusivamente vía `AppConfig`/`ServicesConfig`.

---

## Etapa 1.5 — Preparativos para Calidad y Pruebas *(En curso)*
Puente entre la modularización (Etapa 1) y la ejecución de pruebas (Etapa 2/3).

**Progreso reciente**
- Creada `tests/` con `pytest.ini`, fixtures compartidos y workbook qPCR sintético (ver `tests/conftest.py`).
- Pruebas de humo para ingesta qPCR (`parse_qpcr_wide`, `melt_wide_to_long`) y `compute_fold_change` activas en CI local (`pytest`).

**Acciones clave pendientes**
- Documentar escenarios prioritarios (qPCR, heurística, NLP, enriquecimiento) y criterios de aceptación.
- Incorporar datasets sintéticos adicionales para bibliografía y enriquecimiento.
- Evaluar `tox`/`nox` o flujos alternativos para ejecutar matrices de pruebas y cobertura.

**Entregables**
- Carpeta `tests/` inicial con fixtures y ejemplos *(completado)*.
- Guía corta de ejecución de pruebas y expectativas de cobertura *(documentada en `README.md`)*.

**Criterios de salida**
- Fixtures listos para uso en etapas siguientes.
- Ejecución de `pytest` local sin configuración manual adicional.

---

## Etapa 2 — Calidad de Código y Seguridad *(Pendiente)*
**Objetivo**: Establecer controles automáticos que garanticen calidad mínima y manejo correcto de datos.

**Acciones**
- Añadir `pre-commit` con `ruff`, `black`, `isort` y verificación de tipos (`mypy`/`pyright`).
- Ejecutar auditoría de dependencias (`pip-audit` o `safety`), fijar versiones críticas.
- Implantar revisiones de pares con checklist de seguridad (manejo de secretos, validación de entradas).
- Documentar políticas de gestión de datos (retención, anonimización, exportaciones).

**Entregables**
- `.pre-commit-config.yaml` + guía rápida de uso.
- Registro de hallazgos de seguridad y plan de remediación.

**Criterios de salida**
- Commits bloqueados si no pasan linters, formato o type-checks.
- Dependencias críticas monitoreadas por CVEs y actualizaciones planificadas.

---

## Etapa 3 — Estrategia de Pruebas y Datos de Referencia *(Pendiente)*
**Objetivo**: Asegurar el comportamiento de los componentes principales con suites reproducibles.

**Acciones**
- Consolidar `tests/` con `pytest` y fixtures qPCR/bibliografía/heurística.
- Cubrir funciones núcleo (`parse_qpcr_wide`, servicios `fold_change`, `string_enrichment`, `heuristics`, `nlp`).
- Implementar pruebas de contrato para servicios externos (PubMed, STRING, Ensembl) con mocks.
- Añadir smoke tests para la UI (Streamlit-testing/Playwright headless).

**Entregables**
- Reporte de cobertura inicial y metas de mejora (>70% en módulos críticos).
- Guía para generar/actualizar datasets de prueba.

**Criterios de salida**
- `pytest` estable y integrado al flujo de desarrollo.
- Defectos críticos reproducibles mediante pruebas automatizadas.

---

## Etapa 4 — Observabilidad y Manejo de Errores *(Pendiente)*
**Objetivo**: Garantizar trazabilidad y diagnósticos rápidos en ejecución.

**Acciones sugeridas**
- Centralizar logging estructurado (JSON/`structlog`) con niveles configurables.
- Añadir tracking de métricas clave (tiempos de API, cache hit ratio, errores por módulo).
- Implementar capa de errores controlados en servicios y UI (mensajes accionables + IDs de incidente).
- Documentar procedimientos de soporte y tableros básicos (Grafana/Streamlit metrics).

**Entregables**
- Config de logging unificada y guía de monitoreo.
- Checklist de respuesta ante incidentes y plantillas de reporte.

**Criterios de salida**
- Incidentes reproducibles con logs/metricas suficientes.
- Alertas tempranas cuando fallan integraciones externas o cacheos.

---

## Etapa 5 — Automatización CI/CD y Despliegues Confiables *(Pendiente)*
**Objetivo**: Automatizar build, pruebas y despliegue con controles consistentes.

**Acciones sugeridas**
- Definir pipelines para lint/test/build (GitHub Actions, GitLab CI u otro).
- Empaquetar app (`pip install -e .`) y generar artefactos de despliegue reproducibles (Docker).
- Definir estrategias de despliegue (staging → prod) con aprobaciones manuales cuando aplique.
- Integrar escaneo de dependencias y publicación de cambios en release notes.

**Entregables**
- Workflow CI/CD versionado con pasos claros.
- Dockerfile y documentación de ejecución en contenedor.

**Criterios de salida**
- Cada merge desencadena pipeline automática y controlada.
- Despliegues reproducibles con rollback definido.

---

## Etapa 6 — Gobernanza, Documentación y Onboarding *(Pendiente)*
**Objetivo**: Consolidar documentación, procesos y onboarding para nuevos integrantes.

**Acciones sugeridas**
- Actualizar `README.md`, `web_app/README.md` y guías de troubleshooting.
- Crear manual de onboarding con checklist de setup local, estándares y flujos de trabajo.
- Estandarizar plantillas de PR/issue y políticas de versionado.
- Añadir ejemplos anonimizados (`raw_data/`) y capturas de flujo completo.

**Entregables**
- Paquete de documentación actualizado (docs + assets).
- Checklist de onboarding y de handoff para releases.

**Criterios de salida**
- Nuevos integrantes completan onboarding en <1 semana.
- Documentación cubre ejecución, pruebas, despliegue y soporte.

---

## Recomendaciones transversales
- Revisar retrospectivamente el cierre de cada etapa para ajustar alcance de la siguiente.
- Etiquetar tareas del backlog con la etapa correspondiente para visibilidad.
- Mantener métricas de progreso (porcentaje completado, bloqueos, deuda técnica restante).
- Usar feature flags o toggles para introducir cambios estructurales sin interrumpir usuarios finales.
- Consultar `MEJORAS.md` para granularidad adicional y referencias por archivo/línea.

---

## Próximos pasos sugeridos (inicio Q4)
1. Cerrar Etapa 1 formalmente: mover secciones a `app/ui/pages`, publicar plantillas de configuración y actualizar documentación asociada.
2. Ejecutar Etapa 1.5: crear `tests/` con fixtures sintéticos y pipelines básicos (`pytest`, cobertura).
3. Preparar infraestructura para Etapa 2: establecer `pre-commit`, definir responsabilidades de revisión y plan de auditoría de dependencias.
