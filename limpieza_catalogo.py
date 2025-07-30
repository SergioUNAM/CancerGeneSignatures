import json

# Rutas de entrada y salida
archivo_entrada = 'catalogo.json'
archivo_salida = 'catalogo_florida.json'

# Cargar el catálogo
with open(archivo_entrada, 'r', encoding='utf-8') as f:
    catalogo = json.load(f)

# Filtrar artículos que NO contienen "no usar" en la descripción (ignorando mayúsculas)
catalogo_filtrado = [
    item for item in catalogo
    if "no usar" not in item.get("descripcion", "").lower()
]

# Guardar el nuevo archivo filtrado
with open(archivo_salida, 'w', encoding='utf-8') as f:
    json.dump(catalogo_filtrado, f, ensure_ascii=False, indent=2)

print(f"Filtrado completado. Artículos restantes: {len(catalogo_filtrado)}")