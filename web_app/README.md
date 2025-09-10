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

Notas
- La app lee `web_app/config/menu.json` para los parámetros (contexto, tipo de cáncer, método preferido).
- Para Excel se usa openpyxl; si hay formatos antiguos, conviértelos a .xlsx.
