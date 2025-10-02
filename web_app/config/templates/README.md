# Plantillas de configuración por entorno

Este directorio contiene ejemplos de configuración para distintos entornos de despliegue.
Las plantillas ayudan a homogeneizar valores sensibles y evitar que `streamlit_app.py`
se configure directamente en código.

## Contenido

- `.env.dev.example`: variables pensadas para desarrollo local.
- `.env.staging.example`: configuración recomendada para un entorno de pruebas
  compartido.
- `.env.prod.example`: valores sugeridos para producción.
- `menu.dev.json`: menú de navegación acotado para experimentación.
- `menu.shared.json`: menú base para `staging` y `prod`.

## Uso

1. Copia la plantilla adecuada y elimina la terminación `.example`.
2. Ajusta los valores según el entorno, sin comprometer secretos en el repositorio.
3. Exporta las variables (o súbelas a `st.secrets`) antes de ejecutar la app.
4. Opcionalmente establece `CGS_MENU_PATH` apuntando al menú deseado.

```bash
cp web_app/config/templates/.env.dev.example .env
source .env
streamlit run web_app/streamlit_app.py
```

Para `staging`/`prod`, sube el contenido del `.env` al gestor de secretos de tu
plataforma y coloca el archivo `menu.shared.json` en un bucket accesible.
