Ejecución local de la Web App (MVP)

Requisitos
- Python 3.9+
- Entorno virtual recomendado

Instalación
```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
pip install -r web_app/requirements.txt
```

Ejecutar
```bash
streamlit run web_app/streamlit_app.py
```

Flujo principal
- Sube un Excel (.xlsx/.xls) y selecciona la hoja qPCR si aplica.
- Clasifica controles y muestras por prefijo, sufijo, regex o selección manual; revisa la previsualización y limpia colisiones si es necesario.
- Aplica la política "Undetermined/ND" (nan, ctmax o valor fijo) y revisa la calidad de datos antes de continuar.
- Calcula fold change (promedios o gen de referencia) y genera la tabla de expresión categorizada.
- Ejecuta la normalización avanzada: ajusta parámetros (α, candidatos, K referencias, bootstrap, permutaciones) y obtén df_norm, ranking de referencias y estadísticas diferenciales.
- Anota genes vía Ensembl (IDs y descripciones) y descarga los resultados.

Notas
- La app usa `web_app/config/menu.json` para poblar tipos de cáncer, contextos y otras listas.
- Si trabajas con Excel heredados, conviértelos a `.xlsx` para una importación estable (se usa `openpyxl`).
- `CGS_LOGLEVEL` permite ajustar el nivel de logging (INFO por defecto).
- Dependencias opcionales como STRING, PubMed, Google NLP o firmas se reintroducirán en módulos separados cuando se reactiven.
- Para mantener estado entre ejecuciones se persiste la sesión en `st.session_state`; usa "Limpiar clasificación" si necesitas empezar de cero.
