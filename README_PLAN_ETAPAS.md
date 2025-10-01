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

## Etapa 1 — Modularización y Gestión de Configuración
- **Objetivo**: Reducir acoplamiento y centralizar configuraciones para facilitar pruebas y despliegues.
- **Acciones**
  - Extraer lógica de `web_app/streamlit_app.py` en paquetes `ui`, `services`, `domain` y `adapters`.
  - Introducir modelos tipados para configuración (p. ej. `pydantic`/`dataclasses`).
  - Unificar carga de settings desde `config/`, variables de entorno y `st.secrets`.
  - Establecer rutas de importación relativas y evitar `sys.path` dinámicos.
- **Entregables**
  - Nueva estructura de carpetas documentada en el README principal.
  - Plantillas de configuración por entorno (`dev`, `staging`, `prod`).
- **Criterios de salida**
  - Ningún módulo excede 400 líneas y la lógica de negocio es invocable sin UI.
  - Los parámetros críticos se resuelven vía objetos de configuración validados.

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
  - Configurar estructura `tests/` con `pytest` y fixtures de datos sintéticos qPCR/bibliografía.
  - Cubrir funciones núcleo (`parse_qpcr_wide`, `compute_fold_change`, enriquecimientos) con pruebas unitarias y de integración.
  - Implementar pruebas de contrato para servicios externos simulando respuestas (PubMed, STRING, Ensembl).
  - Añadir pruebas de smoke para la UI con `streamlit-testing` o Playwright headless.
- **Entregables**
  - Reporte de cobertura inicial y metas de mejora (+70% en módulos críticos).
  - Guía para generar/actualizar datasets de prueba.
- **Criterios de salida**
  - `pytest` se ejecuta sin flakiness y forma parte del flujo de desarrollo.
  - Los defectos críticos son reproducibles mediante alguna prueba automatizada.

## Etapa 4 — Observabilidad y Manejo de Errores
- **Objetivo**: Mejorar la trazabilidad y resiliencia para diagnósticos rápidos en producción.
- **Acciones**
  - Estandarizar logging estructurado (JSON) con correlación por sesión/usuario.
  - Centralizar manejo de excepciones y mensajes de error amigables en la UI.
  - Instrumentar métricas clave (tiempo de carga, llamadas externas, fallos) y exponerlas vía Prometheus o logs.
  - Configurar alertas básicas (p. ej. en Slack/Email) para errores críticos o degradaciones.
- **Entregables**
  - Configuración de logging por entorno y guía de interpretación.
  - Tablero de observabilidad o reporte periódico.
- **Criterios de salida**
  - Los incidentes se pueden rastrear con timestamp + ID de sesión.
  - Alertas se disparan antes de que el usuario reporte degradaciones severas.

## Etapa 5 — Automatización CI/CD y Despliegues Confiables
- **Objetivo**: Garantizar pipelines reproducibles y despliegues controlados.
- **Acciones**
  - Definir `pyproject.toml` y lockfile (`uv`, `poetry` o `pip-tools`) para dependencias reproducibles.
  - Crear Dockerfile multi-stage y scripts de entrada (`make`, `task`, `nox`) para comandos comunes.
  - Configurar CI (GitHub Actions/GitLab) para ejecutar lint, tests, build de imagen y escaneo de seguridad.
  - Establecer estrategias de despliegue (staging→prod) con checklist de rollback y backups.
- **Entregables**
  - Pipeline CI/CD documentado con diagramas de flujo.
  - Imagen base publicada en registro interno o público.
- **Criterios de salida**
  - Cada merge en `main` genera artefactos/versiones etiquetadas automáticamente.
  - Los despliegues se ejecutan con una orden y tienen rollback probado.

## Etapa 6 — Gobernanza, Documentación y Onboarding
- **Objetivo**: Facilitar la incorporación de nuevos colaboradores y asegurar continuidad operativa.
- **Acciones**
  - Actualizar `README.md`, `CONTRIBUTING.md` y `docs/` con guías de arquitectura, decisiones ADR y troubleshooting.
  - Formalizar política de versiones (semver), branching model y definición de Done.
  - Implementar reviews post-mortem y métricas de calidad (lead time, MTTR) trimestrales.
  - Mantener changelog/Release notes automatizadas.
- **Entregables**
  - Documentación navegable (mkdocs o similar) con índice completo.
  - Playbook de soporte y plan de capacitación.
- **Criterios de salida**
  - Onboarding < 5 días hábiles para nuevos devs.
  - Cada release documenta cambios, riesgos y pasos de validación.

---

## Recomendaciones Transversales
- Programar retrospectivas al cierre de cada etapa para ajustar alcance siguiente.
- Etiquetar tareas en el backlog con la etapa correspondiente para visibilidad.
- Mantener métricas de progreso (porcentaje completado, bloqueos, deuda técnica pendiente).
- Usar feature flags para introducir cambios estructurales sin interrumpir usuarios finales.

