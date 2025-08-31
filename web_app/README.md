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
- Sube Excel/CSV, selecciona hoja (si aplica).
- Define columnas ID (opcional) y normalización.
- Revisa resultados y descarga ZIP.

Notas
- La app guarda resultados en memoria para descarga; no escribe en disco por defecto.
- Para Excel se usa openpyxl; si hay formatos antiguos, conviértelos a .xlsx.

