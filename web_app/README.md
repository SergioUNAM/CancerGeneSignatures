Ejecución local de la Web App (MVP)

Requisitos
- Python 3.9+
- Entorno virtual recomendado

Instalación
```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\\Scripts\\activate
pip install -r web_app/requirements.txt
```

Ejecutar
```bash
streamlit run web_app/streamlit_app.py
```

Flujo
- Sube Excel (.xlsx/.xls), selecciona hoja (si aplica).
- (Opcional) Introduce prefijos de controles/muestras o usa los detectados de la plantilla.
- Revisa clasificación, Fold Change (promedios y gen de referencia) y gráficas.
- Anotación Ensembl (IDs/descripciones) y descarga de resultados.
- Enriquecimiento STRING por nivel de expresión (GO/KEGG/Reactome), filtros y descargas.
- Bibliografía (PubMed): ingresa tu NCBI Email (obligatorio) y API Key (opcional) directamente en la sección. La app mostrará progreso por gen, tablas y gráficos, y permitirá descargar CSV.
- Firmas genéticas: genera firmas por tipo de cáncer/nivel desde la bibliografía clasificada y enriquece Hallmarks (MSigDB). Incluye visualización Sunburst y descarga CSV. Requiere `gseapy` y archivos GMT locales (ruta por defecto en `gen-sets_GSEA_MSigDB/`).

Notas
- La app lee `web_app/config/menu.json` para los parámetros (contexto, tipo de cáncer, método preferido).
- Para Excel se usa openpyxl; si hay formatos antiguos, conviértelos a .xlsx.
- Logs: se puede ajustar el nivel con la variable de entorno `CGS_LOGLEVEL` (INFO por defecto).
