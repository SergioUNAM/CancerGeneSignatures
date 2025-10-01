# Plan de mejoras para la arquitectura de CancerGeneSignatures Web App

Este plan consolida acciones para robustecer la estructura actual (`app/core`, `app/services`, UI Streamlit) y facilitar la incorporación de nuevas funcionalidades.

## Objetivos generales
- Mantener una separación clara entre lógica de dominio (`app.core`) y orquestación/UI (`app.services`, `app.ui`).
- Simplificar la incorporación de nuevos módulos o funcionalidades reutilizando patrones existentes.
- Reducir deuda técnica previa (imports, estado de sesión, documentación y pruebas).

## Acciones priorizadas

### 1. Modularizar el flujo de UI (Alta prioridad)
- [ ] Extraer secciones principales de `web_app/streamlit_app.py` a helpers en `app.ui.sections` (e.g. `render_fold_change_section`, `render_advanced_normalization_section`).
- [ ] Documentar los parámetros y estructuras de retorno de cada helper.
- Resultado esperado: archivo principal más corto y fácil de mantener; posibilidad de reutilizar secciones en páginas futuras.

### 2. Formalizar interfaces de servicios (Alta)
- [ ] Agregar docstrings detalladas en `app/services/*.py` describiendo inputs/outputs y excepciones.
- [ ] Añadir anotaciones de tipos en puntos críticos y validar con `mypy` (al menos en modo `--strict-optional`).
- [ ] Revisar que toda manipulación de `st.session_state` se haga vía `AppSessionState` o funciones auxiliares centralizadas.

### 3. Pruebas unitarias para módulos core (Media)
- [ ] Crear suite de tests adicional: `tests/core/test_reference_normalization.py`, `tests/core/test_qpcr.py`, etc., apuntando a funciones de `app.core`.
- [ ] Configurar `pytest` para correr en CI/local con scripts sencillos (`make test` o `poetry run pytest`).
- [ ] Documentar cómo correr las pruebas en README o PROGRESO.

### 4. Documentación y guías de contribución (Media)
- [ ] Actualizar `web_app/PROGRESO.md` y README con la nueva estructura (`app.core`, `app.services`).
- [ ] Redactar una guía corta en `docs/` explicando patrones de nombres (core vs services) y flujo de datos.

### 5. Organización de constantes y textos (Baja)
- [ ] Definir un módulo `app.ui.texts` (o similar) con textos/etiquetas recurrentes para la UI.
- [ ] Usar enums o constantes para llaves de `st.session_state` y nombres de archivos descargables.

### 6. Seguimiento y validación continua
- [ ] Establecer un checklist de “nueva funcionalidad” (tests, docstrings, actualización de PROGRESO).
- [ ] Revisar periódicamente dependencias (`requirements.txt`) y eliminar paquetes no usados.

## Métricas sugeridas
- % de módulos en `app/core` cubiertos por tests unitarios.
- Tiempo promedio en agregar nuevos métodos de normalización (disminuir a medida que se modulariza).
- Número de warnings de linters/typing (objetivo: cero en cada merge).

## Nota final
Este documento se debe revisar después de cada entrega relevante. Actualiza los checkboxes y añade comentarios con fecha para mantener un historial de decisiones.
