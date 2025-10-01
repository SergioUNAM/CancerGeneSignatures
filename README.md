# CancerGeneSignatures
Generador de firmas genéticas asociadas a tipos de cáncer mediante análisis de expresión génica, integración bibliográfica y construcción de redes de interacción. La meta es ofrecer un flujo reproducible que ayude a equipos oncológicos a priorizar genes candidatos, explorar literatura relevante y producir reportes consistentes.

## Panorama rápido
- Pipeline completo: importación de qPCR, cálculos de Fold Change, enriquecimiento STRING/Hallmarks y síntesis bibliográfica.
- Servicios especializados desacoplados de Streamlit (`web_app/app/services/*`) para favorecer pruebas y reutilización.
- UI modular multipágina en construcción, apoyada en secciones reutilizables (`web_app/app/ui/sections.py`).
- Configuración centralizada (`AppConfig`, `ServicesConfig`) con soporte para distintos entornos y credenciales externas.

## Flujo funcional
1. **Ingesta qPCR** (`app/services/qpcr.py`): limpieza, clasificación de muestras/controles, imputación de Ct y panel de calidad.
2. **Fold Change y métricas** (`app/services/fold_change.py`): reglas de `Undetermined`, control de NaN por grupo y cálculo de métricas.
3. **Anotación y enriquecimiento** (`app/services/string_enrichment.py`, `app/services/signatures.py`): integración Ensembl/STRING, hallmarks, filtros y consolidación.
4. **Bibliografía y heurística** (`app/services/bibliography.py`, `app/services/heuristics.py`): consultas PubMed, etiquetado heurístico y visualizaciones.
5. **Insights NLP** (`app/services/nlp.py`): generación de corpus, etiquetado por tema y agregación con resultados de expresión.
6. **UI multipágina (WIP)** (`app/ui/sections.py`, `app/ui/pages/*`): vistas independientes para qPCR, enriquecimiento, bibliografía e insights.

## Arquitectura del repositorio
```
CancerGeneSignatures/
├── web_app/                 # Streamlit app (bootstrap, configuración, servicios, UI, estado)
│   ├── streamlit_app.py
│   └── app/
│       ├── config/          # Modelos de configuración y cargadores tipados
│       ├── services/        # Lógica de negocio desacoplada
│       ├── ui/              # Componentes y páginas (en expansión)
│       └── state/           # Manejo de `st.session_state`
├── src/                     # Núcleo reusable (core analytics, integraciones)
│   ├── core/
│   └── integrations/
├── docs/                    # ADRs, especificaciones y planes de modularización
├── notebooks/               # Experimentación y validación exploratoria
├── raw_data/                # Datasets fuente (mantener anonimizado)
├── resultados/              # Salidas generadas durante el análisis
├── MEJORAS.md               # Backlog priorizado por fases
├── README_PLAN_ETAPAS.md    # Roadmap detallado
└── README.md                # Este documento
```

## Requisitos previos
- Python 3.10+ (recomendado 3.11 para aprovechar mejoras recientes).
- `pip` o `uv`; opcionalmente `conda` para gestionar dependencias científicas.
- Credenciales opcionales para APIs externas: PubMed (email + api key), STRING, Google Natural Language.

## Instalación rápida
```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
pip install -r requirements.txt          # dependencias núcleo
pip install -r web_app/requirements.txt  # dependencias específicas de la app Streamlit
```

## Ejecución de la app
```bash
streamlit run web_app/streamlit_app.py
```
- Se puede definir `CGS_MENU_PATH` para apuntar a un menú alternativo (`web_app/config/menu.json` por defecto).
- `CGS_LOGLEVEL`, `CGS_PUBMED_EMAIL`, `CGS_PUBMED_API_KEY` y `CGS_GOOGLE_NLP_API_KEY` permiten ajustar logs y credenciales.
- Para despliegues, sobrescribir variables mediante `st.secrets` o ficheros `.env` gestionados fuera del control de versiones.

## Trabajo con datos
- Colocar ficheros qPCR en `raw_data/` o cargarlos desde la UI; mantener identificadores anonimizados.
- Las llamadas a PubMed/STRING se cachean mediante `@st.cache_data` (TTL configurable en `app/services/*`).
- Referencias a datasets de ejemplo y fixtures se documentarán en la sub-etapa 1.5 (`README_PLAN_ETAPAS.md`).

## Pruebas rápidas
- Ejecutar `pytest` en la raíz para lanzar las pruebas básicas (`tests/`).
- Los fixtures sintéticos viven en `tests/conftest.py` e incluyen un workbook qPCR en memoria para los servicios de ingesta.
- `pytest.ini` fija rutas y filtros mínimos; se puede ampliar con marcadores según crezca la batería.

## Documentación relacionada
- `docs/etapa1_modularizacion_config.md`: decisiones sobre configuración y modularización actual.
- `MEJORAS.md`: backlog ordenado por fases con estado granular.
- `README_PLAN_ETAPAS.md`: roadmap detallado y criterios de salida.

## Roadmap activo
- **Etapa 1 (en curso)**: separación completa en `app/ui/pages`, plantillas de configuración y documentación final de la arquitectura.
- **Etapa 1.5 (en curso)**: ampliación de fixtures sintéticos, estructura `tests/` y automatización inicial.
- **Etapas 2-3**: calidad, pruebas y cobertura con CI. Detalles completos en el README de etapas.

## Contribuir
1. Crear rama descriptiva (`feature/`, `fix/`, `docs/`).
2. Mantener la lógica en `app/services` o `src/core` sin dependencias explícitas de Streamlit.
3. Añadir pruebas (cuando existan fixtures) y documentar decisiones clave en `docs/`.
4. Enviar PR incluyendo referencias a tareas del roadmap o checklist en `MEJORAS.md`.

## Licencia
Proyecto bajo Licencia MIT (`LICENSE`).
