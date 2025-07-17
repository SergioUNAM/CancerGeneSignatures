import pandas as pd
import json

# Leer los archivos Excel
areas_df = pd.read_excel('areas_filtradas.xlsx', dtype={'ean': str})
catalog_df = pd.read_excel('catalogo_florida.xlsx', dtype={'ean': str})

# Renombrar columna para adaptarse al JSON de salida
catalog_df = catalog_df.rename(columns={'product_name': 'descripcion'})

# Seleccionar solo las columnas necesarias
catalog_df = catalog_df[['ean', 'descripcion', 'stock']]
areas_df = areas_df[['ean', 'area']]

# Unir los DataFrames por 'ean'
merged_df = pd.merge(catalog_df, areas_df, on='ean', how='left')

# Convertir a lista de diccionarios (JSON)
productos = merged_df.to_dict(orient='records')

# Exportar a JSON
with open('catalogo.florida.json', 'w', encoding='utf-8') as f:
    json.dump(productos, f, ensure_ascii=False, indent=2)

print("JSON generado correctamente con", len(productos), "productos.")