# CancerGeneSignatures
Generador de firmas genéticas asociadas a tipos de cáncer mediante análisis de expresión génica, integración bibliográfica y construcción de redes de interacción. Facilita la identificación de genes clave y su aplicación en investigación oncológica, diagnóstico o pronóstico.

## 🚀 Flujo General

1. **Carga y análisis qPCR** (`web_app/app/services/qpcr.py`): clasificación de controles/muestras, imputación de Ct y cálculo de Fold Change.
2. **Anotación y enriquecimiento** (`app/services/string_enrichment.py`, `app/services/signatures.py`): integración con Ensembl, STRING y Hallmarks/MSigDB.
3. **Bibliografía y heurística** (`app/services/bibliography.py`, `app/services/heuristics.py`): búsqueda PubMed, síntesis heurística por gen y visualizaciones.
4. **Insights NLP** (`app/services/nlp.py`): preparación de corpus y análisis con Google Natural Language API.
5. **UI modular** (`app/ui/sections.py`, próximamente `app/ui/pages/*`): cada sección se renderiza a partir de los servicios anteriores.

## 📁 Estructura Principal
```
web_app/
├── streamlit_app.py        # Bootstrap, navegación y composición de secciones
├── app/
│   ├── config/             # Modelos y loader de configuración (AppConfig, ServicesConfig)
│   ├── services/           # Lógica de negocio (qpcr, fold_change, string_enrichment, bibliography, heuristics, nlp, visuals)
│   ├── ui/                 # Secciones reutilizables y (próximamente) páginas multipage
│   └── state/              # Estado persistente de sesión Streamlit
└── config/menu.json        # Menú por defecto (tipos de cáncer, contextos, normalización)
```

## 🛠️ Instalación Rápida
```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
pip install -r web_app/requirements.txt
streamlit run web_app/streamlit_app.py
```

### Variables útiles
- `CGS_MENU_PATH`: ruta alternativa al menú JSON.
- `CGS_LOGLEVEL`: nivel de logging (`INFO`, `DEBUG`, etc.).
- `CGS_PUBMED_EMAIL`, `CGS_PUBMED_API_KEY`: credenciales por defecto para PubMed.
- `CGS_GOOGLE_NLP_API_KEY`: clave de Google NLP.

## 📊 Servicios Destacados
- `app/services/qpcr.py`: orquestación del flujo qPCR (clasificación, imputación, resumen).
- `app/services/fold_change.py`: políticas de “Undetermined”, métricas de calidad y cálculo de FC.
- `app/services/string_enrichment.py`: llamadas a STRING y filtrado final.
- `app/services/bibliography.py`: consulta a PubMed y fusión con niveles de expresión.
- `app/services/heuristics.py`: síntesis heurística de la bibliografía clasificada.
- `app/services/nlp.py`: preparación de corpus y agregación de resultados de Google NLP.

## 🧭 Roadmap (resumen)
- **Etapa actual (1)**: modularización + UI distribuida. Pendiente mover cada sección a `app/ui/pages`, documentar estructura y publicar plantillas de configuración.
- **Próxima sub-etapa (1.5)**: preparar fixtures/datasets sintéticos y base `tests/` para abordar Calidad (Etapa 2) y Pruebas (Etapa 3).
- **Etapas posteriores**: ver `README_PLAN_ETAPAS.md` para el detalle de CI/CD, observabilidad y gobernanza.

## 🤝 Contribuciones
Las contribuciones son bienvenidas. Recomendaciones actuales:
- Mantener servicios libres de dependencias Streamlit.
- Crear pruebas para cada módulo nuevo.
- Documentar decisiones clave (ADR) en `docs/`.

## 📄 Licencia
Proyecto bajo Licencia MIT (`LICENSE`).
