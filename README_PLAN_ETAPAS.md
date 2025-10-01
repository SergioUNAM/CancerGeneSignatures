# Roadmap de Mejora por Etapas — Mejores Prácticas de Ingeniería

Documento orientado a guiar la adopción progresiva de mejores prácticas de ingeniería de software en **CancerGeneSignatures**. Cada etapa define objetivos, acciones clave, entregables y criterios de finalización. El enfoque está en robustecer calidad, mantenibilidad, observabilidad y operaciones sin introducir nuevas funciones de negocio.

---

## Etapa 0 — Diagnóstico y Bases Compartidas
- **Objetivo**: Alinear al equipo con el estado actual del código y acordar estándares mínimos.
- **Acciones**
  - Inventariar módulos críticos (`src/core`, `web_app`) y dependencias externas.
  - Documentar decisiones arquitectónicas actuales y deudas técnicas visibles.
  - Definir un estándar de estilo y tipado (PEP 8 + hints tipo `mypy`/`pyright`).
  - Configurar un tablero Kanban con backlog priorizado de mejoras.
- **Entregables**
  - Informe corto de diagnóstico con riesgos y prioridades.
  - Documento de estándares de código y convenciones de logs/errores.
- **Criterios de salida**
  - Todo el equipo comprende el alcance y las limitaciones actuales.
  - Existe un glosario común y criterios listos para etapas posteriores.

## Etapa 1 — Modularización y Gestión de Configuración (EN PROGRESO)
- **Objetivo**: Reducir acoplamiento y centralizar configuraciones para facilitar pruebas y despliegues.
- **Avances clave**
  - `web_app/streamlit_app.py` reducido a bootstrapping; la lógica reside en `app/services` y `app/ui`.
  - Configuración tipada (`AppConfig`, `ServicesConfig`) y lectura tolerante a entornos sin `secrets`.
  - Servicios especializados: `qpcr`, `fold_change`, `string_enrichment`, `bibliography`, `heuristics`, `nlp`.
  - Secciones UI reutilizables (`app/ui/sections.py`).
- **Acciones restantes**
  - Crear `app/ui/pages` para separar análisis qPCR, enriquecimiento, firmas y reportes.
  - Añadir adaptadores ligeros si se requiere integrar nuevas fuentes (ej. caching de Ensembl).
- **Entregables pendientes**
  - Actualización del README principal con la nueva estructura modular.
  - Plantillas de configuración por entorno (`dev`, `staging`, `prod`).
- **Criterios de salida (actualizados)**
  - Cada flujo de la web app se invoca via `app/ui/sections` o páginas dedicadas (sin lógica extensa en Streamlit).
  - Configuraciones y credenciales se resuelven exclusivamente via `AppConfig` / `ServicesConfig`.

## Etapa 1.5 — Preparativos para Calidad y Pruebas *(nueva sub-etapa propuesta)*
Puente entre la modularización y la implementación de pruebas (Etapa 2/3).
- Definir estructura `tests/` y datos sintéticos para servicios extraídos.
- Documentar escenarios de prueba prioritarios (qPCR, heurística, NLP).
- Configurar `pytest` básico con cobertura e integración inicial en CI (job “test” opcional).

## Etapa 2 — Calidad de Código y Seguridad
- **Objetivo**: Establecer controles automáticos que garanticen calidad mínima y manejo correcto de datos.
- **Acciones**
  - Añadir `pre-commit` con `ruff`, `black`, `isort` y verificación de tipos.
  - Ejecutar auditoría de dependencias (`pip-audit`/`safety`) y definir política de actualizaciones.
  - Implantar revisiones de pares obligatorias con checklist de seguridad (manejo de secretos, validación de entradas).
  - Documentar políticas de gestión de datos (retención, anonimización, exportaciones).
- **Entregables**
  - Archivo `.pre-commit-config.yaml` y guía rápida de uso.
  - Registro de hallazgos de seguridad y plan de remediación.
- **Criterios de salida**
  - Commits rechazados si no pasan linters, formatos o type-checks.
  - Las dependencias críticas tienen versiones fijadas y se monitorean CVEs.

## Etapa 3 — Estrategia de Pruebas y Datos de Referencia
- **Objetivo**: Asegurar el comportamiento de los componentes principales con suites de pruebas reproducibles.
- **Acciones**
  - Configurar estructura `tests/` con `pytest` y fixtures qPCR/bibliografía/heurística.
  - Cubrir funciones núcleo (`parse_qpcr_wide`, servicios `fold_change`, `string_enrichment`, `heuristics`, `nlp`).
  - Implementar pruebas de contrato para servicios externos simulando respuestas (PubMed, STRING, Ensembl).
  - Añadir pruebas de smoke para la UI (`streamlit-testing`, Playwright headless).
- **Entregables**
  - Reporte de cobertura inicial y metas de mejora (+70% en módulos críticos).
  - Guía para generar/actualizar datasets de prueba.
- **Criterios de salida**
  - `pytest` se ejecuta sin flakiness y forma parte del flujo de desarrollo.
  - Defectos críticos reproducibles mediante pruebas automatizadas.

## Etapa 4 — Observabilidad y Manejo de Errores
- ... *(sin cambios respecto al plan original)*

## Etapa 5 — Automatización CI/CD y Despliegues Confiables
- ...

## Etapa 6 — Gobernanza, Documentación y Onboarding
- ...

---

## Recomendaciones Transversales
- Programar retrospectivas al cierre de cada etapa para ajustar alcance siguiente.
- Etiquetar tareas del backlog con la etapa correspondiente para visibilidad.
- Mantener métricas de progreso (porcentaje completado, bloqueos, deuda técnica pendiente).
- Usar feature flags para introducir cambios estructurales sin interrumpir usuarios finales.

---

## Siguiente Etapa Propuesta (Etapa 1 a finalizar)
1. **UI multipágina**: crear `app/ui/pages` (`analysis.py`, `enrichment.py`, `bibliography.py`, `insights.py`) y habilitar `streamlit.navigation` o patrón de sidebar para seleccionar vistas.
2. **Documentación**: actualizar `README.md` con la nueva arquitectura, flujo de ejecución y dependencias.
3. **Pruebas preparatorias**: generar fixtures mínimos para servicios qPCR/NLP y definir objetivos de cobertura de cara a la Etapa 3.
