# Etapa 1 — Plan de Modularización y Gestión de Configuración

Este documento desglosa el trabajo de la Etapa 1 con foco en desacoplar la web app, introducir un modelo de configuración tipado y sentar las bases para futuras automatizaciones. No implica aún modificar código; define el plan de refactor.

---

## 1. Objetivos Claros
- Separar responsabilidades de UI, orquestación y dominio actualmente concentradas en `web_app/streamlit_app.py`.
- Definir un sistema de configuración centralizado, validado y versionado para todos los entornos.
- Eliminar dependencias implícitas (por ejemplo, `sys.path` manual) y preparar empaquetado limpio.

---

## 2. Diagnóstico Resumido
- `web_app/streamlit_app.py` (~2k líneas) mezcla:
  - Renderizado Streamlit.
  - Acceso a datos, caches, orquestación de procesos y llamadas a integraciones.
  - Definiciones de utilidades y fallbacks.
- Los parámetros de la app se dispersan entre `streamlit_app.py`, `config/menu.json`, variables de entorno y valores "duros".
- `AppSessionState` maneja pocas piezas del estado; no existe un objeto de configuración global.
- No hay paquete instalable para la web app (el código depende de rutas relativas).

---

## 3. Arquitectura Objetivo (alto nivel)

```
web_app/
├── app/                      # Paquete instalable (renombrar a opcional "cgs_web")
│   ├── __init__.py
│   ├── config/               # Definiciones de configuración y cargadores
│   │   ├── models.py         # dataclasses/pydantic para settings
│   │   └── loader.py         # merge menú JSON + env vars + defaults
│   ├── state/
│   │   └── session.py        # Manejo de session state tipado
│   ├── services/
│   │   ├── qpcr.py           # Orquestación de carga/pre-procesamiento
│   │   ├── signatures.py     # Generación de firmas
│   │   ├── enrichment.py     # STRING / Hallmarks
│   │   └── bibliography.py   # PubMed + interpretaciones
│   ├── ui/
│   │   ├── layout.py         # Layout general y navegación
│   │   ├── components.py     # Componentes reutilizables (tabs, tablas, gráficos)
│   │   └── pages/
│   │       ├── qpcr.py
│   │       ├── enrichment.py
│   │       ├── bibliography.py
│   │       └── signatures.py
│   ├── adapters/
│   │   ├── io_adapter.py     # Conecta servicios con `src.core.io`
│   │   ├── string_adapter.py # Configuración de cliente STRING
│   │   └── pubmed_adapter.py # Solicitudes y cache PubMed
│   └── app.py                # Punto de entrada Streamlit (mínimo)
├── config/
│   ├── menu.json             # Sigue siendo fuente de datos, validada
│   └── menu.schema.json      # (nuevo) especificación para validación
└── streamlit_app.py          # Wrapper que invoca `app.app.main()`
```

- `src/core/*` permanece como librería de dominio; sólo se crean adaptadores si hace falta transformar inputs/outputs.
- Las páginas pueden registrarse mediante un router simple (`st.Page`, `multipage`) manteniendo estado centralizado.

---

## 4. Sistema de Configuración

### 4.1 Modelos propuestos (`app/config/models.py`)
- `AppConfig`
  - `project_name`, `log_level`, `cache_ttl`
  - `menu: MenuConfig`
  - `services: ServicesConfig`
- `MenuConfig`
  - `version`, `cancer_types`, `contexts`, `normalization_methods`
  - Validaciones: listas no vacías, claves presentes, labels únicos.
- `ServicesConfig`
  - `ensembl: EnsemblConfig`
  - `string: StringConfig`
  - `pubmed: PubMedConfig`
  - `google_nlp: GoogleNLPConfig`
- Fijar defaults en un módulo `defaults.py` para facilitar actualización.

### 4.2 Cargador (`app/config/loader.py`)
1. Leer defaults.
2. Combinar con variables de entorno (`os.getenv`) y `st.secrets` si disponible.
3. Validar `config/menu.json` contra `MenuConfig` (p. ej., usando `pydantic` o `jsonschema`).
4. Generar instancia inmutable (`AppConfig`) cacheada en `st.session_state`.

### 4.3 Uso
- Servicios y páginas reciben `AppConfig` (o sub-secciones) por inyección explícita.
- `AppSessionState` almacena selecciones del usuario ligadas a `AppConfig` (p. ej., `context_label` debe pertenecer a `MenuConfig.contexts`).

---

## 5. Plan de Migración Incremental

1. **Preparación**
   - Crear paquete `web_app/app/` con estructura vacía.
   - Mover `AppSessionState` a `app/state/session.py`; adaptar importaciones sin cambiar comportamiento.
   2. **Config**
      - Implementar `config/models.py` y `loader.py` con tests unitarios.
      - Reemplazar uso directo de `MENU = load_menu()` por `config = load_app_config()`.
      - [x] Unificar credenciales vía `AppConfig.services` (PubMed, Google NLP) eliminando dependencias de `os.getenv` en la UI.
   3. **Servicios y Adaptadores**
      - Extraer funciones puras de `streamlit_app.py` hacia `services/` empezando por flujo qPCR (carga + clasificación + fold change).
        - [x] `app/services/qpcr.py` contiene `build_long_table`, `summarize_extraction` y clasificadores reutilizables (prefijos, sufijos, regex, selección, resolución de colisiones).
        - [x] `app/services/fold_change.py` encapsula imputación de Ct, métricas de calidad y generación de tablas de expresión.
        - [x] `app/services/string_enrichment.py` concentra la ejecución y filtrado de enriquecimiento STRING.
        - [x] `app/services/bibliography.py` expone helpers para PubMed y anexo de niveles de expresión.
        - [x] `app/services/heuristics.py` centraliza la interpretación heurística y datos derivados para visualizaciones.
        - [x] `app/services/nlp.py` prepara corpus y orquesta llamadas a Google NLP para la sección de insights.
   - Crear adaptadores finos para integraciones externas (PubMed, STRING, Google).
4. **UI y Páginas**
   - Dividir `streamlit_app.py` en páginas bajo `ui/pages/` manteniendo primero el flujo qPCR.
   - Simplificar `streamlit_app.py` a Bootstrapping (`from app.app import main; main()`).
5. **Refactor final**
   - Eliminar `sys.path` manual al convertir `web_app/app` en paquete instalable (`pip install -e web_app`).
   - Actualizar documentación (`README`, `MEJORAS`) con nueva estructura.

> Ejecutar pruebas `pytest` y una smoke manual en cada paso para garantizar regresión mínima.

---

## 6. Consideraciones Técnicas
- **Tipado**: preferible usar `pydantic` v2 para validaciones rápidas; si se evita dependencia adicional, usar `dataclasses` + validaciones manuales.
- **Compatibilidad**: mantener nombres de claves en `st.session_state` actual para no romper sesiones guardadas; introducir migrações si se renombra.
- **Empaquetado**: opcional crear `pyproject.toml` dentro de `web_app/` en etapas posteriores (anticipar import absoluto `from app ...`).
- **Testing**: cada servicio extraído debe acompañarse de pruebas (definidas en Etapa 3, pero conviene preparar archivos de ejemplo y fixtures).

---

## 7. Próximos Pasos
1. Revisar este plan con el equipo (arquitectura + datos) y ajustar nombres/alcances.
2. Priorizar qué servicios extraer primero (sugerencia: `qpcr` → `enrichment` → `bibliography` → `signatures`).
3. Crear issues/tareas por subpaso con criterios de aceptación claros.
4. Programar la implementación iterativa (feature flags si es necesario) dentro del sprint.

---

## 8. Riesgos y Mitigaciones
- **Rupturas por refactor masivo** → Migrar por etapas con wrappers temporales y pruebas de humo frecuentes.
- **Desalineación de configuraciones** → Mantener `menu.schema.json` y validación automática en CI.
- **Sobrecarga cognitiva** → Documentar diagramas simples (flujo de datos, dependencias) y mantener README actualizado.
