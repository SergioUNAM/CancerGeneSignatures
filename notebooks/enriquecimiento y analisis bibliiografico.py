# %% [markdown]
# # Funciones de enriquecimiento genetico: API de STRING y gseapy
# 
# ## Enriquecimiento STRING
# 
# ### Agregar Información de Ensembl en Lotes
# 
# Este código agrega información de **Ensembl** (IDs y descripciones) a un DataFrame que contiene nombres de genes. La consulta se realiza en lotes para evitar sobrecargar el servidor de Ensembl.  
# 
# **Funcionalidad**:  
# 1. **Consulta en Lotes**:  
#    - Los genes se dividen en lotes de tamaño `batch_size` (por defecto, 50 genes por lote).  
#    - Se realiza una solicitud HTTP POST a la API de Ensembl para obtener los IDs y descripciones de los genes en cada lote.  
# 
# 2. **Manejo de Respuestas**:  
#    - Si la consulta es exitosa, se extraen los IDs y descripciones de los genes.  
#    - Si falla, se asigna "Not found" a los genes del lote.  
# 
# 3. **Resultados**:  
#    - Se agregan dos columnas al DataFrame: `ensembl_id` (ID de Ensembl) y `description` (descripción del gen).  
# 
# **Uso**:  
# - Aplica la función `add_ensembl_info_batch` a DataFrames que contengan genes subexpresados, estables y sobreexpresados para enriquecerlos con información de Ensembl.  

# %%
import requests
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
from IPython.display import display, Markdown

def get_ensembl_id(gene_symbol):
    """Busca en Ensembl con tolerancia a alias y nombres alternativos"""
    base_url = "https://rest.ensembl.org"
    endpoints = [
        f"/lookup/symbol/homo_sapiens/{gene_symbol}?expand=1",
        f"/xrefs/symbol/homo_sapiens/{gene_symbol}?"
    ]

    for endpoint in endpoints:
        try:
            response = requests.get(base_url + endpoint, headers={"Content-Type": "application/json"})
            if response.status_code == 200:
                data = response.json()
                if isinstance(data, list):  # Para el endpoint xrefs
                    for item in data:
                        if item['type'] == 'gene':
                            return item['id'], item.get('description', 'No description')
                else:  # Para el endpoint lookup
                    return data['id'], data.get('description', 'No description')
        except Exception as e:
            continue

    # Si falla Ensembl, intentar con MyGene.Info
    mygene_response = requests.get(f"https://mygene.info/v3/query?q={gene_symbol}&species=human")
    if mygene_response.status_code == 200:
        results = mygene_response.json().get('hits', [])
        for hit in results:
            if 'ensembl' in hit:
                ensembl_id = hit['ensembl'].get('gene', 'Not found')
                description = hit.get('summary', 'No description')
                return ensembl_id, description

    return 'Not found', 'No description'

def process_gene(gene):
    """Función wrapper para procesamiento paralelo"""
    ensembl_id, description = get_ensembl_id(gene)

    # Búsqueda aproximada si no se encuentra
    if ensembl_id == 'Not found':
        try:
            fuzzy_response = requests.get(
                f"https://rest.ensembl.org/lookup/search/homo_sapiens?q={gene}",
                headers={"Content-Type": "application/json"},
                timeout=10
            )
            if fuzzy_response.status_code == 200:
                fuzzy_data = fuzzy_response.json()
                if fuzzy_data.get('hits'):
                    best_match = fuzzy_data['hits'][0]
                    return best_match['id'], best_match.get('description', 'No description')
        except:
            pass

    return ensembl_id, description

def add_ensembl_info_batch(df, max_workers=5):
    """Versión con procesamiento paralelo"""
    targets = df['target'].tolist()

    # Procesamiento paralelo
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        results = list(executor.map(process_gene, targets))

    ids, descriptions = zip(*results)

    df['ensembl_id'] = ids
    df['description'] = descriptions

    # Estadísticas detalladas
    not_found = df[df['ensembl_id'] == 'Not found']
    if not not_found.empty:
        alert_message(f"Genes no encontrados ({len(not_found)}): {', '.join(not_found['target'].tolist()[:5])}" +
                     ("..." if len(not_found) > 5 else ""))
    else:
        success_message("¡Todos los genes fueron encontrados exitosamente!")

    result_message(f"Resumen:\n- Total genes: {len(df)}\n- Encontrados: {len(df) - len(not_found)}\n- No encontrados: {len(not_found)}")

    return df


# Función de alerta modificada
def display_alert_if_genes_not_found(df, dataset_name):
    not_found_genes = df[df['ensembl_id'] == 'Not found']
    if not not_found_genes.empty:
        alert_message(f"Genes no encontrados en {dataset_name}: {', '.join(not_found_genes['target'].unique()[:3])}..." +
                     f" ({len(not_found_genes)} totales)")
    else:
        success_message(f"Todos los genes en {dataset_name} fueron encontrados en Ensembl.")

# Ejemplo de uso
display(Markdown("### Agregar Información de Ensembl a Genes Subexpresados"))
df_sobreexpresados = add_ensembl_info_batch(df_sobreexpresados, max_workers=8)
display_alert_if_genes_not_found(df_sobreexpresados, "genes subexpresados")

display(Markdown("### Agregar Información de Ensembl a Genes Sobreexpresados"))
df_subexpresados = add_ensembl_info_batch(df_subexpresados, max_workers=8)
display_alert_if_genes_not_found(df_subexpresados, "genes sobreexpresados")

df_sobreexpresados['nivel_expresion'] = 'Sobreexpresados'
df_subexpresados['nivel_expresion'] = 'Subexpresados'

df_emsamble_info = pd.concat([df_sobreexpresados, df_subexpresados], ignore_index=True)

# Exportar el DataFrame con información de Ensembl

export_dfs_to_excel(
    [df_emsamble_info],
    ["informacion_ensembl_genes"],
    "IDs y descripciones por gen: Archivo que contiene la información de Ensembl para los genes categorizados.")

# %% [markdown]
# ### Visualización de la distribución de niveles de expresión
# 
# > Aplica para las lista de genes
# 
# Antes de iniciar el enriquecimiento, veamos la distribución de los niveles de expresión

# %%
import plotly.graph_objects as go
from plotly.colors import sample_colorscale

# Orden deseado de las categorías
orden_niveles = ['estable', 'subexpresado',  'sobreexpresado']

# Calcular frecuencias y reindexar
expression_counts = df_combined_expresion['nivel_expresion'].value_counts()
expression_counts = expression_counts.reindex(orden_niveles, fill_value=0)

# Generar colores usando escala divergente nativa de Plotly
colors = sample_colorscale('RdYlBu_r', [0.0, 0.4, 0.8])  # Azul -> Amarillo -> Rojo

# Crear gráfico interactivo
fig = go.Figure()

fig.add_trace(
    go.Bar(
        y=expression_counts.index,
        x=expression_counts.values,
        orientation='h',
        marker=dict(
            color=colors,
            line=dict(color='#333', width=1)
        ),
        text=expression_counts.values,
        texttemplate='<b>%{text}</b>',
        textposition='outside',
        hovertemplate=(
            "<b>%{y}</b><br>"
            "Frecuencia: %{x}<br>"
            "<extra></extra>"
        )
    )
)

# Personalización avanzada
fig.update_layout(
    title=dict(
        text='<b>Distribución de Niveles de Expresión</b><br><sub>Análisis de expresión génica</sub>',
        x=0.5,
        font=dict(size=24, family='Times New Roman', color='#2C3E50'),
    ),
    xaxis=dict(
        title='Frecuencia',
        showgrid=True,
        gridcolor='rgba(150,150,150,0.2)',
        linecolor='#333',
        mirror=True,
        zeroline=False
    ),
    yaxis=dict(
        title='Nivel de Expresión',
        type='category',
        linecolor='#333',
        mirror=True,
        tickfont=dict(size=14)
    ),
    plot_bgcolor='white',
    height=500,
    width=1000,
    margin=dict(t=120, b=80, l=120),
    hoverlabel=dict(
        bgcolor='white',
        font_size=14,
        bordercolor='#333'
    )
)

fig.show()

guardar_grafico(fig, 'distribucion_niveles_expresion.png')

# %% [markdown]
# ### Visualización de fold change por nivel de expresión
# 
# > datos: listas de genes
# 
# Este código genera una visualización en subplots para cada nivel de expresión génica (**Subexpresado**, **Estable**, **Sobreexpresado**), mostrando el **Fold Change** de los genes en cada categoría.
# 
# 
# **Uso*:
# - Esta visualización es útil para comparar el `fold_change` de los genes en cada nivel de expresión y identificar patrones o valores atípicos.  

# %%
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plotly.colors import sample_colorscale, qualitative


# Filtrar solo los niveles que existen en los datos
unique_levels = [level for level in orden_niveles if level in df_combined_expresion['nivel_expresion'].unique()]
num_levels = len(unique_levels)

# Crear figura con subplots
fig = make_subplots(
    rows=num_levels,
    cols=1,
    subplot_titles=[f'Nivel de expresión: {level}' for level in unique_levels],
    vertical_spacing=0.15,
    specs=[[{"secondary_y": False}] for _ in range(num_levels)]
)

# Configurar estilo común
axis_style = dict(
    showgrid=True,
    gridcolor='rgba(200,200,200,0.3)',
    linecolor='#444',
    mirror=True,
    zeroline=False
)

for i, level in enumerate(unique_levels, 1):
    # Filtrar y ordenar datos
    data_subset = df_combined_expresion[df_combined_expresion['nivel_expresion'] == level]
    data_subset = data_subset.sort_values('fold_change', ascending=False)

    if data_subset.empty:
        continue

    num_bars = len(data_subset)

    # Generar colores con escala secuencial de Plotly
    if num_bars > 1:
        colors = sample_colorscale('Agsunset', [n/(num_bars-1) for n in range(num_bars)])
    else:
        colors = [qualitative.Plotly[0]]  # Color base de Plotly

    # Crear gráfico de barras
    fig.add_trace(
        go.Bar(
            x=data_subset['target'],
            y=data_subset['fold_change'],
            marker=dict(
                color=colors,
                line=dict(color='rgba(0,0,0,0.5)', width=1)
            ),
            name=level,
            hovertemplate=(
                "<b>%{x}</b><br>"
                "Fold Change: %{y:.2f}<br>"
                "Nivel: %{fullData.name}"
                "<extra></extra>"
            )
        ),
        row=i,
        col=1
    )

    # Configurar escalas
    min_val = data_subset['fold_change'].min()
    max_val = data_subset['fold_change'].max()

    yaxis_config = dict(
        type='log' if (max_val / min_val > 100 and min_val > 0) else 'linear',
        title='Fold Change (log)' if (max_val / min_val > 100) else 'Fold Change',
        tickformat=".1f",
        **axis_style
    )

    fig.update_yaxes(yaxis_config, row=i, col=1)
    fig.update_xaxes(
        title='Genes',
        tickangle=90,
        categoryorder='total descending',
        **axis_style,
        row=i,
        col=1
    )

# Configuración global del layout
fig.update_layout(
    height=400 * num_levels + 200,
    width=1200,
    title=dict(
        text='<b>Análisis de Expresión Génica</b>',
        x=0.5,
        font=dict(size=24, family='Arial')
    ),
    plot_bgcolor='white',
    hoverlabel=dict(
        bgcolor='white',
        font_size=12,
        font_family="Arial"
    ),
    margin=dict(t=150, b=100),
    showlegend=False,
    bargap=0.4  # Espaciado entre barras
)

# Añadir anotación de crédito
fig.add_annotation(
    text="Visualización generada con Plotly",
    xref="paper",
    yref="paper",
    x=0.5,
    y=-0.15,
    showarrow=False,
    font=dict(color='#666', size=10)
)

fig.show()

guardar_grafico(fig, 'analisis_expresion_genica.png')

# %% [markdown]
# ### Obtención de los datos de enriquecimiento funcional con API STRING
# 
# 
# > datos: listas de genes
# 
# Este código realiza un **análisis de enriquecimiento funcional** utilizando la API de STRING para identificar procesos biológicos, rutas y funciones moleculares asociadas a los genes en cada nivel de expresión (**Subexpresado**, **Estable**, **Sobreexpresado**).
# 
# **Uso**:
# - Este análisis es útil para identificar procesos biológicos clave asociados a los genes en cada nivel de expresión, lo que puede ayudar a interpretar los resultados del análisis diferencial.  

# %%
from IPython.display import display, Markdown, HTML
import requests
import pandas as pd

# Definir la URL base y el formato de salida para la API de STRING
url_api_base = "https://version-12-0.string-db.org/api"
formato_salida = "json"
metodo = "enrichment"

# Construir la URL de solicitud
url_solicitud = "/".join([url_api_base, formato_salida, metodo])

# Definir una función para realizar el enriquecimiento funcional
def realizar_enriquecimiento_STRING(lista_genes, descripcion):
    """
    Realiza un análisis de enriquecimiento funcional utilizando la API de STRING.

    Parámetros:
    lista_genes (list): Lista de genes a analizar.
    descripcion (str): Descripción del nivel de expresión (subexpresados, estables, sobreexpresados).

    Retorna:
    pd.DataFrame: DataFrame con los resultados del enriquecimiento funcional.
    """
    if not lista_genes:
        alert_message(f"La lista de genes está vacía. No se puede realizar el enriquecimiento para {descripcion}")
        return pd.DataFrame()

    # Convertir lista en string separado por '%0d' como lo requiere la API
    genes_principales = "%0d".join(lista_genes)
    parametros = {
        "identifiers": genes_principales,
        "species": 9606,  # Código taxonómico para Homo sapiens
        "caller_identity": "UIMEO"
    }

    try:
        respuesta = requests.post(url_solicitud, data=parametros)
        respuesta.raise_for_status()
        datos = respuesta.json()
        success_message(f"Enriquecimiento funcional completado para {descripcion}")
        return pd.DataFrame(datos)
    except requests.exceptions.HTTPError as errh:
        error_message(f"Error HTTP en {descripcion}:{errh}")
    except requests.exceptions.ConnectionError as errc:
        error_message(f"Error de conexión en {descripcion}: {errc}")
    except requests.exceptions.Timeout as errt:
        error_message(f"Tiempo de espera agotado en {descripcion}: {errt}")
    except requests.exceptions.RequestException as err:
        error_message(f"Error en la solicitud para {descripcion}: {err}")
    except ValueError:
        error_message(f"Error al decodificar la respuesta JSON para {descripcion}")
    return pd.DataFrame()

# Realizar el enriquecimiento para cada nivel de expresión
display(Markdown("### Análisis de Enriquecimiento Funcional con STRING"))

display(Markdown("#### Enriquecimiento Funcional para Genes Subexpresados"))
enriquecimiento_sobreexpresados_STRING = realizar_enriquecimiento_STRING(lista_sobreexpresados, "genes sobreexpresados")
enriquecimiento_sobreexpresados_STRING['nivel_expresion'] = "Sobreexpresados"

display(Markdown("#### Enriquecimiento Funcional para Genes Estables"))
enriquecimiento_estables_STRING = realizar_enriquecimiento_STRING(lista_estables, "genes estables")
enriquecimiento_estables_STRING['nivel_expresion'] = "Estables"

display(Markdown("#### Enriquecimiento Funcional para Genes Subexpresados"))
enriquecimiento_subexpresados_STRING = realizar_enriquecimiento_STRING(lista_subexpresados, "genes subexpresados")
enriquecimiento_subexpresados_STRING['nivel_expresion'] = "Subexpresados"

# Concatenar los tres DataFrames en uno solo
df_enriquecimiento_completo_STRING = pd.concat([
    enriquecimiento_sobreexpresados_STRING,
    enriquecimiento_estables_STRING,
    enriquecimiento_subexpresados_STRING
])

# Exportar el DataFrame concatenado a un archivo Excel

display(Markdown("### Exportación de Resultados de Enriquecimiento Funcional a Excel"))

export_dfs_to_excel(
    [df_enriquecimiento_completo_STRING, enriquecimiento_sobreexpresados_STRING, enriquecimiento_estables_STRING, enriquecimiento_subexpresados_STRING],
    ["enriquecimiento_completo_STRING", "enriquecimiento_sobreexpresados_STRING", "enriquecimiento_estables_STRING", "enriquecimiento_subexpresados_STRING"],
    "Enriquecimiento funcional completo API-STRING: Archivo que contiene los resultados del enriquecimiento funcional con STRING, cada subgrupo de expresión en una hoja separada."
)



# %% [markdown]
# ### Construcción de un diccionario  de significados de categorías de enriquecimiento funcional
# 
# Este diccionario proporciona una descripción detallada de las categorías utilizadas en el análisis de enriquecimiento funcional. Cada categoría está asociada con un significado específico que ayuda a interpretar los resultados del análisis.

# %%
from IPython.display import display, Markdown, HTML

# Diccionario de significados de las categorías ampliado y detallado
categorias_significados = {
    "HPO": (
        "Human Phenotype Ontology: Un sistema de clasificación que describe los fenotipos observados en humanos. "
    ),
    "Process": (
        "Procesos biológicos generales asociados a los genes: Actividades celulares, fisiológicas o moleculares llevadas a cabo por uno o más genes."
    ),
    "PMID": (
        "Publicaciones científicas referenciadas mediante PubMed ID: Proporciona enlaces a estudios y artículos revisados por pares disponibles en PubMed."
    ),
    "RCTM": (
        "Reactome: Base de datos centrada en vías metabólicas y de señalización celular."
    ),
    "COMPARTMENTS": (
        "Localizaciones subcelulares y estructuras biológicas: Describe las ubicaciones dentro de la célula, como el núcleo, citoplasma, membrana plasmática y organelos."
    ),
    "WikiPathways": (
        "Rutas metabólicas y biológicas de la base WikiPathways: Una base de datos colaborativa que documenta y organiza rutas biológicas. "
    ),
    "KEGG": (
        "Rutas y procesos biológicos anotados en KEGG (Kyoto Encyclopedia of Genes and Genomes)."
    ),
    "Component": (
        "Componentes celulares asociados a los genes: Identifica las estructuras específicas dentro de la célula donde los genes o sus productos desempeñan un papel funcional."
    ),
    "TISSUES": (
        "Asociaciones específicas con tejidos biológicos: Describe la relación entre genes y su expresión preferencial o específica en tejidos como cerebro, hígado, pulmón o tejidos tumorales."
    ),
    "Keyword": (
        "Palabras clave relacionadas con las funciones genéticas: Etiquetas utilizadas para clasificar y agrupar genes según sus roles funcionales."
    ),
    "DISEASES": (
        "Relación con enfermedades humanas conocidas: Vincula genes y términos con patologías específicas como cáncer, enfermedades cardiovasculares o trastornos genéticos."
    ),
    "Function": (
        "Funciones moleculares de los genes: Define las actividades bioquímicas llevadas a cabo por los productos génicos, como proteínas y ARN."
    ),
    "NetworkNeighborAL": (
        "Genes conectados en redes de interacción cercanas: Representa la proximidad funcional o física de genes en redes biológicas."
    ),
    "SMART": (
        "Dominios específicos de proteínas anotados en SMART (Simple Modular Architecture Research Tool): Identifica regiones funcionales dentro de las proteínas."
    ),
    "InterPro": (
        "Clasificación y anotación de familias de proteínas: Integra múltiples bases de datos para identificar relaciones evolutivas y dominios funcionales en proteínas."
    ),
    "Pfam": (
        "Familias de proteínas definidas en la base Pfam: Ofrece información detallada sobre agrupaciones de proteínas relacionadas por secuencia y estructura."
    )
}

# Mostrar el diccionario de significados
display(Markdown("### Diccionario de Significados de Categorías"))

result_message("Diccionario de categorias cargado con exito")

# %% [markdown]
# ### Procesamiento y visualización de los resultados de enriquecimiento Funcional
# 
# > Uso exclusivo enriquecimiento STRING
# 
# Este código define una función `graficar_data_enriquecimiento` que toma un DataFrame de resultados de enriquecimiento funcional y realiza lo siguiente:
# 
# **Uso**:
# - La función se aplica a los DataFrames de enriquecimiento funcional para genes **sobreexpresados**, **estables** y **subexpresados**.  
# - Proporciona una visión clara y detallada de las categorías más relevantes en cada nivel de expresión, facilitando la interpretación de los resultados.  
# 

# %%
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go


# Preparar datos
def graficar_data_enriquecimiento(df):
    df_grouped = df.groupby(['category', 'nivel_expresion']).size().reset_index(name='count')
    df_categorias = pd.DataFrame.from_dict(categorias_significados, orient='index', columns=['descripcion']).reset_index().rename(columns={'index':'category'})
    df_plot = df_grouped.merge(df_categorias, on='category')

    # Ordenar por conteo total
    total_counts = df_plot.groupby('category')['count'].sum().sort_values(ascending=False)
    df_plot['category'] = pd.Categorical(df_plot['category'], categories=total_counts.index, ordered=True)
    df_plot = df_plot.sort_values('category')

    color_map = {
        "Sobreexpresados": "#FBB4AE",
        "Estables": "#D2E8B0",
        "Subexpresados": "#B3CDE3",
    }

    # Crear gráfico base
    fig = px.bar(
        df_plot,
        x='category',
        y='count',
        color='nivel_expresion',
        color_discrete_map=color_map,
        title='<b>Distribución de Categorías de Enriquecimiento por Nivel de Expresión</b>',
        labels={'count': 'Número de Términos', 'category': 'Categoría'},
        hover_data=['descripcion'],
        height=700,
    )

    # Personalización avanzada
    fig.update_layout(
        plot_bgcolor='rgba(245,245,245,1)',
        paper_bgcolor='white',
        font=dict(family='Arial', size=12, color='#2c3e50'),
        xaxis=dict(
            title='',
            tickangle=45,
            showgrid=False,
            tickfont=dict(size=14)
        ),
        yaxis=dict(
            title='Número de Términos',
            showgrid=True,
            gridcolor='#e0e0e0',
            title_font=dict(size=16)
        ),
        hoverlabel=dict(
            bgcolor='white',
            font_size=12,
            font_family='Arial'
        ),
        legend=dict(
            title='Nivel de Expresión',
            orientation='v',
            yanchor='bottom',
            y=0.7,
            xanchor='center',
            font=dict(size=14),
            # Mover a un costado del grafico
            x=0.9,
        ),
        margin=dict(l=50, r=50, b=200, t=100),  # Margen superior aumentado
        title=dict(
            x=0.5,
            y=0.95,
            font=dict(size=24, color='#2c3e50', family='Arial Black')
        ),
    )

    # Mejorar tooltips
    fig.update_traces(
        hovertemplate="<b>%{x}</b><br>"
                    "Nivel: %{customdata[0]}<br>"
                    "Términos: %{y}<br>"
                    "<extra></extra>"
    )

    # Añadir anotación explicativa
    fig.add_annotation(
        x=0.5,
        y=-0.4,
        xref='paper',
        yref='paper',
        text="* Los datos de significancia (FDR < 0.05) fueron ajustados mediante el método de Benjamini-Hochberg",
        showarrow=False,
        font=dict(size=10, color='gray')
    )

    # Configurar botones de filtrado
    botones = [
        dict(
            label="Mostrar Todos",
            method="update",
            args=[{"visible": [True]*len(fig.data)}, {"title": "<b>Todos los niveles de expresión</b>"}]
        ),
        dict(
            label="Solo Sobreexpresados",
            method="update",
            args=[{"visible": [d.name == 'Sobreexpresados' for d in fig.data]},
                {"title": "<b>Sobreexpresados</b>"}]
        ),
        dict(
            label="Solo Subexpresados",
            method="update",
            args=[{"visible": [d.name == 'Subexpresados' for d in fig.data]},
                {"title": "<b>Subexpresados</b>"}]
        )
    ]

    # Añadir interactividad mejorada
    fig.update_layout(
        updatemenus=[
            dict(
                type="buttons",
                direction="right",
                active=0,
                buttons=botones,
                pad={"r": 20, "t": 10},
                x=0.5,
                xanchor="center",
                y=1,  # Posición ajustada
                yanchor="top",
                font=dict(size=14)
            )
        ]
    )

    # Asegurar visibilidad inicial correcta
    for trace in fig.data:
        trace.visible = True

    fig.show()

    guardar_grafico(fig, 'distribucion_categorias_enriquecimiento.png')

# Graficar los datos de enriquecimiento
graficar_data_enriquecimiento(df_enriquecimiento_completo_STRING)

# %% [markdown]
# ### Filtrado y ordenado de datos de enriquecimiento GO
# 
# > Uso exclusivo enriquecimiento STRING
# 
# Este código define una función `filtrar_enriquecimiento_go` que procesa un DataFrame de resultados de enriquecimiento funcional para incluir solo términos relacionados con **Gene Ontology (GO)**. Luego, ordena los datos y elimina columnas innecesarias.
# 
# **Uso**:
# - Esta función es útil para preparar los datos de enriquecimiento GO para análisis posteriores, como la visualización o la integración con otros resultados.

# %%
def filtrar_enriquecimiento_terminos_GO(df):
    """
    Filtra un DataFrame para incluir solo las filas donde la columna 'term' comienza con 'GO:',
    ordena el DataFrame por el valor de 'fdr' en orden ascendente y elimina las columnas
    'ncbiTaxonId' y 'preferredNames'.

    Parámetros:
    dataframe (pd.DataFrame): DataFrame que contiene los datos de enriquecimiento con al menos las columnas 'term', 'fdr', 'ncbiTaxonId' y 'preferredNames'.
    category (str): Categoría de genes (Sobreexpresados, Estables, Subexpresados).

    Retorna:
    pd.DataFrame: DataFrame filtrado, ordenado y con las columnas especificadas eliminadas.

    Excepciones:
    ValueError: Si el DataFrame no contiene las columnas necesarias.
    """
    # Validar que las columnas requeridas están presentes
    columnas_requeridas = ['term', 'fdr', 'ncbiTaxonId', 'preferredNames']
    columnas_faltantes = [col for col in columnas_requeridas if col not in df.columns]
    if columnas_faltantes:
        raise ValueError(f"El DataFrame no contiene las columnas requeridas: {', '.join(columnas_faltantes)}")

    # Filtrar filas donde 'term' comienza con "GO:"
    dataframe_filtrado = df[df['term'].str.startswith("GO:")]

    # Ordenar el DataFrame por 'fdr' en orden ascendente
    dataframe_ordenado = dataframe_filtrado.sort_values(by='fdr', ascending=True)

    # Eliminar las columnas 'ncbiTaxonId' y 'preferredNames'
    dataframe_ordenado = dataframe_ordenado.drop(['ncbiTaxonId', 'preferredNames'], axis=1)
    # Reiniciar los índices y descartar los antiguos
    dataframe_ordenado.reset_index(drop=True, inplace=True)

    return dataframe_ordenado



# %% [markdown]
# #### Aplicación de la función: `filtrar_enriquecimiento_terminos_GO` sobre los datos obtenidos de STRING

# %%
# # Aplicar la función de filtrado a los diferentes conjuntos de datos
# display(Markdown("### Filtrado y Ordenamiento de Términos GO de los datos obtenidos vía STRING"))
# try:
#     enriquecimiento_GO_STRING_sobreexpresados = filtrar_enriquecimiento_terminos_GO(enriquecimiento_sobreexpresados_STRING)
#     result_message("Filtrado del enriquecimiento GO sobreexpresados")
#     enriquecimiento_GO_STRING_estables = filtrar_enriquecimiento_terminos_GO(enriquecimiento_estables_STRING)
#     result_message("Filtrado del enriquecimiento GO estables")
#     enriquecimiento_GO_STRING_subexpresados = filtrar_enriquecimiento_terminos_GO(enriquecimiento_subexpresados_STRING)
#     result_message("Filtrado del enriquecimiento GO subexpresados")
# except ValueError as e:
#     error_message(f"Error al procesar los datos: {e}")

# enriquecimiento_GO_STRING_completo = pd.concat([
#     enriquecimiento_GO_STRING_sobreexpresados,
#     enriquecimiento_GO_STRING_estables,
#     enriquecimiento_GO_STRING_subexpresados
# ], ignore_index=True)


# export_dfs_to_excel(
#     [enriquecimiento_GO_STRING_completo, enriquecimiento_GO_STRING_sobreexpresados, enriquecimiento_GO_STRING_estables, enriquecimiento_GO_STRING_subexpresados],
#     ["enriquecimiento_GO_completo", "enriquecimiento_GO_sobreexpresados", "enriquecimiento_GO_estables", "enriquecimiento_GO_subexpresados"],
#     "Filtrado términos GO API STRING: Archivo que contiene los resultados de enriquecimiento funcional con STRING para términos GO, cada subgrupo de expresión en una hoja separada.")



# %% [markdown]
# ### Visualización de la distribución de  las categorias de términos GO por nivel de expresión
# 
# > Uso enriquecimiento STRING y gseapy
# 
# Este código define una función `graficar_distribucion_terminos_GO` que genera un gráfico de barras agrupadas para visualizar el número de términos **Gene Ontology (GO)** enriquecidos en cada categoría de genes (**Sobreexpresados**, **Estables**, **Subexpresados**) y por tipo de GO.
# 
# **Uso**:
# - Esta función es útil para comparar visualmente la cantidad de términos GO enriquecidos en cada categoría de genes y tipo de GO, lo que ayuda a identificar patrones o diferencias significativas.

# %%
# import plotly.express as px
# import pandas as pd

# def graficar_distribucion_terminos_GO(df):
#     """
#     Función para graficar el número de términos GO enriquecidos por categoría de genes y tipo de GO usando Plotly.

#     Parámetros:
#     - df: DataFrame con los datos de enriquecimiento GO.
#     """

#     # Contar términos GO por categoría y tipo
#     df_counts = (
#         df
#         .groupby(['nivel_expresion', 'category'])
#         .size()
#         .reset_index(name='Número de términos GO enriquecidos')
#     )

#     # Obtener categorías únicas y ordenar
#     go_categories = df['category'].unique()
#     categoria_orden = ['Sobreexpresados', 'Estables', 'Subexpresados']

#     # Crear MultiIndex para todas las combinaciones posibles
#     multi_index = pd.MultiIndex.from_product(
#         [categoria_orden, go_categories],
#         names=['nivel_expresion', 'category']
#     )

#     # Reindexar para incluir combinaciones faltantes
#     df_counts = (
#         df_counts
#         .set_index(['nivel_expresion', 'category'])
#         .reindex(multi_index, fill_value=0)
#         .reset_index()
#     )

#     # Ordenar categorías
#     df_counts['nivel_expresion'] = pd.Categorical(
#         df_counts['nivel_expresion'],
#         categories=categoria_orden,
#         ordered=True
#     )
#     df_counts.sort_values(['nivel_expresion', 'category'], inplace=True)

#     # Paleta de colores tipo magma
#     magma_colores = ['#0d0887', '#7e03a8', '#cc4778'][:len(go_categories)]

#     # Crear gráfico interactivo
#     fig = px.bar(
#         df_counts,
#         x='nivel_expresion',
#         y='Número de términos GO enriquecidos',
#         color='category',
#         barmode='group',
#         color_discrete_sequence=magma_colores,
#         labels={'category': 'Tipo de GO', 'nivel_expresion': 'Categoría de genes'},
#         title='Términos GO Enriquecidos por Categoría de Genes y Tipo de GO',
#     )

#     # Personalizar diseño
#     fig.update_layout(
#         title_x=0.5,
#         title_font_size=20,
#         xaxis_title_font_size=14,
#         yaxis_title_font_size=14,
#         legend_title_text='Tipo de GO',
#         legend_title_font_size=12,
#         template='plotly_white',
#         hovermode='x unified',
#         bargap=0.15,
#         bargroupgap=0.1
#     )

#     # Ajustar etiquetas
#     fig.update_traces(
#         textfont_size=12,
#         textposition='outside',
#         cliponaxis=False
#     )

#     guardar_grafico(fig, 'distribucion_terminos_GO.png')

#     fig.show()


# %% [markdown]
# #### Aplicación de la función `graficar_distribucion_terminos_GO` sobre los datos obtenidos de STRING
# 
# 

# %%
# # Graficar el enriquecimiento GO
# graficar_distribucion_terminos_GO(enriquecimiento_GO_STRING_completo)

# %% [markdown]
# ### Visualización de términos GO más significativos
# 
# > Uso enriquecimiento STRING y gseapy
# 
# Este código define una función `graficar_top_terminos_GO` que genera gráficos de barras horizontales para visualizar los **términos GO más significativos** (basados en el valor p) en cada categoría de genes (**Sobreexpresados**, **Estables**, **Subexpresados**).
# 

# %%
# import plotly.express as px
# import plotly.graph_objects as go
# from plotly.subplots import make_subplots
# import numpy as np

# def graficar_top_terminos_GO(df, expression_group, n_terms=25):
#     """
#     Grafica los términos GO más significativos usando Plotly con interactividad y estilo moderno.

#     Parámetros:
#     df (pd.DataFrame): DataFrame con los términos GO y sus valores p.
#     expression_group (str): Nombre del grupo de expresión.
#     n_terms (int): Número máximo de términos a mostrar por categoría.
#     """
#     # Filtrar por grupo de expresión
#     group_df = df[df["nivel_expresion"] == expression_group]

#     if group_df.empty:
#         print(f"No hay datos para {expression_group}.")
#         return

#     # Calcular -log10(p-value)
#     group_df = group_df.assign(**{"-log10(p_value)": -np.log10(group_df["p_value"])})

#     # Obtener categorías únicas
#     categories = group_df["category"].unique()

#     # Crear subgráficos
#     fig = make_subplots(
#         rows=len(categories),
#         cols=1,
#         subplot_titles=[f"Categoría: {cat}" for cat in categories],
#         vertical_spacing=0.15
#     )

#     # Configuración de estilo
#     title_font = dict(family="Arial", size=24, color='#2e2e2e')
#     axis_font = dict(family="Arial", size=12, color='#5e5e5e')
#     hover_font = dict(family="Arial", size=12)

#     for i, category in enumerate(categories, 1):
#         # Filtrar y ordenar datos
#         category_df = (group_df[group_df["category"] == category]
#                        .nlargest(n_terms, "-log10(p_value)"))

#         if category_df.empty:
#             continue

#         # Acortar descripciones
#         max_length = 40
#         category_df = category_df.assign(
#             short_desc=category_df["description"].apply(
#                 lambda x: (x[:max_length] + "...") if len(x) > max_length else x
#             )
#         )

#         # Generar paleta de colores
#         n_barras = len(category_df)
#         colores = px.colors.sample_colorscale(
#             px.colors.sequential.PuBu,
#             np.linspace(1, 0, n_barras)  # Invertir escala para rojo = más significativo
#         )

#         # Crear gráfico de barras
#         fig.add_trace(
#             go.Bar(
#                 x=category_df["-log10(p_value)"],
#                 y=category_df["short_desc"],
#                 orientation='h',
#                 marker_color=colores,
#                 hovertext=category_df["description"],  # Texto completo en hover
#                 hoverinfo="text+x",
#                 showlegend=False,
#                 dy=0.5  # Espaciado entre barras
#             ),
#             row=i, col=1,
#         )

#         # Configurar ejes para cada subplot
#         fig.update_yaxes(
#             title_text="Término GO",
#             autorange="reversed",
#             row=i, col=1,
#             tickfont=axis_font,
#             titlefont=axis_font,
#         )

#         fig.update_xaxes(
#             title_text="-log10(p-value)",
#             row=i, col=1,
#             tickfont=axis_font,
#             titlefont=axis_font
#         )

#     # Configuración general del layout
#     fig.update_layout(
#         title={
#             'text': f"Términos GO más significativos - {expression_group}",
#             'y':0.95,
#             'x':0.5,
#             'xanchor': 'center',
#             'yanchor': 'top',
#             'font': title_font
#         },
#         height=600 * len(categories),
#         width=1200,
#         margin=dict(t=150, b=40, l=80, r=40),
#         hoverlabel=dict(
#             bgcolor="white",
#             font=hover_font
#         ),
#         plot_bgcolor='white'
#     )

#     # Ajustar títulos de subplots
#     fig.update_annotations(
#         font=dict(family="Arial", size=14, color='#4a4a4a'),
#         yshift=10
#     )

#     guardar_grafico(fig, f'top_terminos_GO_{expression_group}.png')

#     fig.show()

# %% [markdown]
# #### Aplicación de la función `graficar_top_terminos_GO` sobre los datos obtenidos de STRING

# %%
# display(Markdown("## Términos GO más Significativos por Categoría de Genes Sobreexpresados"))
# graficar_top_terminos_GO(enriquecimiento_GO_STRING_completo, "Sobreexpresados")

# display(Markdown("## Términos GO más Significativos por Categoría de Genes Estables"))
# graficar_top_terminos_GO(enriquecimiento_GO_STRING_completo, "Estables")

# display(Markdown("## Términos GO más Significativos por Categoría de Genes Subexpresados"))
# graficar_top_terminos_GO(enriquecimiento_GO_STRING_completo, "Subexpresados")

# %% [markdown]
# ### Filtrado y ordenado de datos de enriquecimiento PMID (bibliografía pubmed)
# 
# > Uso datos enriquecimiento STRING
# 
# Este segmento de código filtra un **DataFrame** de enriquecimiento funcional para incluir solo los términos referenciados por **PubMed ID (PMID)**, ordena los resultados según el valor de **FDR** y elimina columnas irrelevantes para el análisis.
# 
# Los datos resultantes nos permite obtener un landscape de la tématica sobre la cual se enriquecieron nuestras listas de genes, ya que estos datos solo se obtienen de la API  de STRING solo se aplicará a estos dataframes para el enriquecimiento posterior realizado mediante GSEAPY se obtendran otros datos

# %% [markdown]
# 

# %%
#def filtrar_enriquecimiento_pmid(df):

    # Validar que las columnas requeridas están presentes
    #columnas_requeridas = ['term', 'fdr', 'ncbiTaxonId', 'preferredNames']
    #columnas_faltantes = [col for col in columnas_requeridas if col not in df.columns]
    #if columnas_faltantes:
        #raise ValueError(f"El dataframe no contiene las columnas requeridas: {', '.join(columnas_faltantes)}")

    # Filtrar filas donde 'term' comienza con "PMID:"
    #df_filtrado = df[df['term'].str.startswith("PMID:")]

    # Ordenar el df por 'fdr' en orden ascendente
    #df_ordenado = df_filtrado.sort_values(by='fdr', ascending=True)

    # Eliminar las columnas 'ncbiTaxonId' y 'preferredNames'
    #df_ordenado = df_ordenado.drop(['ncbiTaxonId', 'preferredNames'], axis=1)
    # Reiniciar los índices y descartar los antiguos
    #df_ordenado.reset_index(drop=True, inplace=True)

    # Renombramos a la columna description a 'article_title'
    #df_ordenado.rename(columns={'description': 'article_title'}, inplace=True)

    #return df_ordenado
# Aplicar la función de filtrado a los diferentes conjuntos de datos
#try:
    #enriquecimiento_pmid_sobreexpresados = filtrar_enriquecimiento_pmid(enriquecimiento_sobreexpresados_STRING)
    #enriquecimiento_pmid_sobreexpresados['nivel_expresion'] = "Sobreexpresados"
    #enriquecimiento_pmid_subexpresados = filtrar_enriquecimiento_pmid(enriquecimiento_subexpresados_STRING)
#except ValueError as e:
    #print(f"Error al procesar los datos: {e}")

#enriquecimiento_pmid_completo = pd.concat([enriquecimiento_pmid_sobreexpresados, enriquecimiento_pmid_subexpresados], ignore_index=True)

#success_message("Filtrado de términos PMID completado con éxito")
#alert_message("Se genero el dataframe enriquecimiento_pmid_completo")


# %% [markdown]
# ### Análisis y exportación de estudios relacionados con cáncer y TEM o micro RNAs
# 
# Este bloque de código realiza el análisis de un conjunto de datos bibliográficos, permitiendo elegir la relación a investigar: **cáncer y transición epitelio-mesénquima (TEM)** o **cáncer y micro RNAs** (no de forma simultánea). A partir de la selección, se identifican los estudios relevantes y se exportan los resultados a un archivo Excel.
# 
# > **Nota**: Este código es útil para explorar tendencias en la literatura científica y facilita la organización y visualización de datos complejos en estudios genómicos.

# %% [markdown]
# #### Diccionario de tipos de cáncer y criterio adicional de búsqueda
# 
# A continuación, se define un diccionario que será utilizado en la clasificación automática basada en coincidencias dentro de los títulos o resúmenes de la bibliografía obtenida. Además, se incluye un criterio de búsqueda adicional que permite elegir entre **cáncer y TEM** o **cáncer y micro RNAs** para refinar los resultados.

# %%

cancer_type_keywords = {
    'Breast Cancer': [
        'breast', 'mammary', 'triple negative breast', 'er-positive',
        'pr-positive', 'ductal carcinoma', 'lobular carcinoma', 'estrogen receptor', 'progesterone receptor'
    ],
    'Melanoma': [
        'melanoma', 'cutaneous melanoma', 'skin melanoma', 'uveal melanoma', 'ocular melanoma',
        'acral melanoma', 'braf', 'braf mutation'
    ],
    'Colon Cancer': [
        'colon', 'colorectal', 'rectal', 'intestinal', 'intestine', 'colon adenocarcinoma',
        'rectal adenocarcinoma', 'lynch syndrome', 'hereditary colon cancer'
    ],
    'Hepatocellular Carcinoma': [
        'hepatocellular', 'hcc', 'liver cancer', 'hepatic carcinoma', 'hepatoma', 'cirrhosis',
        'hepatitis b', 'hbv', 'hepatitis c', 'hcv', 'fibrolamellar carcinoma', 'afp', 'alpha-fetoprotein'
    ],
    'Prostate Cancer': [
        'prostate', 'prostatic', 'prostate carcinoma', 'androgen receptor', 'ar-positive', 'psa',
        'prostate-specific antigen', 'gleason score', 'crpc', 'castration-resistant'
    ],
    'Lung Cancer': [
        'lung', 'pulmonary', 'bronchial', 'nsclc', 'non-small cell lung cancer', 'sclc',
        'small cell lung cancer', 'egfr mutation', 'alk rearrangement', 'lung adenocarcinma', 'lung squamous cell carcinoma', 'lung squamous'
    ],
    'Pancreatic Cancer': [
        'pancreatic', 'pancreas', 'pdac', 'ductal adenocarcinoma', 'pancreatic carcinoma',
        'neuroendocrine tumor', 'pancreatic neuroendocrine', 'pancreatic islet cell tumor', 'whipple procedure'
    ],
    'Leukemia': [
        'leukemia', 'leukaemia', 'aml', 'acute myeloid leukemia', 'acute lymphoblastic leukemia',
        'cll', 'chronic lymphocytic leukemia', 'cml', 'chronic myeloid leukemia', 'philadelphia chromosome',
        'bcr-abl', 'mpn', 'myeloproliferative neoplasm'
    ],
    'Lymphoma': [
        'lymphoma', 'hodgkin', 'non-hodgkin', 'nhl', 'b-cell lymphoma', 't-cell lymphoma',
        'dlbcl', 'diffuse large b-cell lymphoma', 'follicular lymphoma', 'mantle cell lymphoma',
        'burkitt lymphoma', 'cutaneous lymphoma'
    ],
    'Ovarian Cancer': [
        'ovarian', 'ovary', 'ovarian carcinoma', 'high-grade serous carcinoma', 'hgsc', 'brca1',
        'brca2', 'ca125', 'clear cell carcinoma', 'granulosa cell tumor'
    ],
    'Cervical Cancer': [
        'cervical', 'cervix', 'cervical carcinoma', 'cervical squamous cell carcinoma', 'cervical adenocarcinoma',
        'hpv', 'human papillomavirus', 'pap smear', 'cervical dysplasia', 'cervical squamous'
    ],
    'Renal Cancer': [
        'renal', 'kidney', 'rcc', 'renal carcinoma', 'clear cell carcinoma', 'papillary renal carcinoma',
        'chromophobe renal carcinoma', 'von hippel-lindau', 'vhl'
    ],
    'Bladder Cancer': [
        'bladder', 'urothelial', 'bladder carcinoma', 'transitional cell carcinoma', 'non-muscle invasive',
        'muscle invasive', 'nmibc', 'mibc'
    ],
    'Glioblastoma': [
        'glioblastoma', 'gbm', 'brain tumor', 'glioma', 'astrocytoma', 'oligodendroglioma', 'idh',
        'idh mutation', 'mgmt methylation', 'temozolomide'
    ],
    'Thyroid Cancer': [
        'thyroid', 'thyroid carcinoma', 'papillary thyroid', 'follicular thyroid', 'medullary thyroid',
        'anaplastic thyroid', 'braf', 'braf mutation', 'ret', 'ret mutation'
    ],
    'Gastric Cancer': [
        'stomach', 'gastric', 'gastric carcinoma', 'diffuse gastric cancer', 'intestinal metaplasia',
        'epstein-barr virus', 'helicobacter pylori', 'h. pylori', 'pylori'
    ],
    'Esophageal Cancer': [
        'esophageal', 'esophagus', 'esophageal carcinoma', 'barrett esophagus', 'esophageal adenocarcinoma',
        'esophageal squamous cell carcinoma', 'squamous Cell Carcinoma of the Esophagus.'
    ],
    'Skin Cancer': [
        'skin cancer', 'basal cell carcinoma', 'bcc', 'squamous cell carcinoma of the skin', 'scc',
        'cutaneous carcinoma', 'actinic keratosis', 'uv radiation', 'squamous cell carcinoma of skin', 'squamous-cell carcinoma of the skin'
    ],
    'Testicular Cancer': [
        'testicular', 'testis', 'germ cell tumor', 'seminoma', 'non-seminoma', 'beta-hcg',
        'alpha-fetoprotein', 'afp'
    ],
    'Endometrial Cancer': [
        'endometrial', 'uterine', 'endometrial carcinoma', 'serous carcinoma', 'lynch syndrome',
        'womb cancer'
    ],
    'Mesothelioma': [
        'mesothelioma', 'pleural mesothelioma', 'peritoneal mesothelioma', 'asbestos', 'asbestos exposure',
        'asbestosis'
    ],
    'Head and Neck Cancer': [
        'head and neck', 'hnc', 'oral squamous cell carcinoma', 'oscc', 'pharyngeal carcinoma',
        'laryngeal carcinoma', 'hpv-positive', 'hpv-negative', 'nasopharyngeal carcinoma' , 'laryngeal squamous cell carcinoma', 'oral squamous', 'naso-oropharyngeal carcinoma', 'tongue squamous'
    ],
    'Bone Cancer': [
    'bone cancer', 'osteosarcoma', 'chondrosarcoma', 'ewing sarcoma', 'primary bone tumor',
    'metastatic bone cancer', 'paget disease', 'bone lesion', 'osteogenic sarcoma'
]
}


emt_keywords = [
    'epithelial-to-mesenchymal transition',
    'epithelial mesenchymal transition',
    'emt',
    'mesenchymal-epithelial transition',
    'transición epitelio-mesenquima',
    'transición epitelio mesénquima',
    'epitelio-mesenquimal',
    'epitelio mesenquima', 'mesenchymal stemstromal', 'cvmscs',
    'epithelial-mesenchymal', 'epithelial to mesenchymal', 'mesenchymal stromal cells',
    'mesenchymalstem', 'mesenchymal transition', 'Epithelial-mesenchymal', 'mesenchymal'
]

microrna_keywords = [
    'microrna', 'mirna', 'mirnas', 'mir', 'mir-', 'miRNA', 'miRNAs',
    'microRNA', 'microRNAs', 'micro-RNA', 'micro-rnas',
    'small non-coding rna', 'sncrna', 'sncrnas',
    'non-coding rna', 'noncoding rna', 'ncRNA', 'ncRNAs',
    'hsa-mir', 'hsa-miR', 'microRNA expression',
    'exosomal micrornas', 'exosomal mirnas', 'mRNA'
]

keywords_cancer = [
    'cancer', 'tumor', 'tumours', 'neoplasm', 'neoplasia', 'carcinoma',
    'malignancy', 'malignant', 'oncology', 'oncogene', 'oncogenic',
    'tumorigenesis', 'tumorigenic', 'tumor suppressor', 'tumour suppressor',
    'metastasis', 'metastatic', 'metastasize', 'metastasizing',
    'carcinogenesis', 'carcinogenic'
]

general_cancer_keywords = ['cancer', 'tumor', 'metastasis']


# %% [markdown]
# ### Clasificación automática de tipos de cáncer en descripciones bibliográficas
# 
# > Uso enriquecimiento STRING y gseapy
# 
# Clasificación de tipos de cáncer según términos de búsqueda en titulos o resúmenes de artículos de PUBMED
# 

# %%
import pandas as pd
import re

def clasificar_tipos_cancer_bibliografia(df, cancer_types=cancer_type_keywords, selected_context="Cáncer y TEM"):
    """
    Clasifica artículos médicos según si se menciona cáncer, añadiendo las siguientes columnas:
      - cancer_relation: True si en el título se detecta mención (ya sea específica o genérica) de cáncer.
      - cancer_type: Tipo(s) de cáncer detectado(s) en el título.
      - emt_relation: True si en el abstract (o título, en ausencia de abstract) se detecta relación con TEM.
      - mRNAs_relation: True si en el abstract (o título) se detecta relación con micro RNAs.

    Además, se filtran los artículos para conservar solo aquellos relacionados con cáncer, y se retorna
    un DataFrame que conserva todas las columnas originales más las nuevas:
      cancer_relation, cancer_type, emt_relation, mRNAs_relation.

    Parámetros:
      - df: DataFrame con datos bibliográficos.
      - cancer_types: Diccionario de tipos de cáncer y sus palabras clave.
      - selected_context: ("Cáncer y TEM" o "Cáncer y micro RNAs"). (No se utiliza para condicionar el cálculo,
                          ya que se generan ambas relaciones).

    Retorna:
      - DataFrame filtrado y con las columnas originales y las nuevas columnas indicadas.
    """
    # Asegurar que el título esté en 'article_title'
    if 'article_title' not in df.columns:
        if 'Title' in df.columns:
            df = df.rename(columns={'Title': 'article_title'})
        else:
            raise ValueError("El DataFrame debe contener la columna 'article_title' o 'Title'.")
    title_column = 'article_title'

    # Usar el abstract si está disponible
    abstract_column = 'Abstract' if 'Abstract' in df.columns else None

    # Función para detectar mención y tipo de cáncer a partir del título
    def classify_cancer_type(text: str, type_dict: dict) -> tuple:
        text_lower = text.lower()
        # Buscar tipos específicos de cáncer según las palabras clave
        detected = [cancer for cancer, keywords in type_dict.items()
                    if any(re.search(rf'\b{keyword}\b', text_lower) for keyword in keywords)]
        if detected:
            return True, ", ".join(detected)
        # Si no se detecta un tipo específico, buscar términos genéricos
        generic_terms = ["cancer", "tumor", "neoplasm", "metastasis", "carcinoma"]
        if any(re.search(rf'\b{term}\b', text_lower) for term in generic_terms):
            return True, "General cancer mention"
        return False, None

    # Función para detectar relación en el texto (abstract o título)
    def detect_relation(text: str, keywords_list: list) -> bool:
        text_lower = text.lower()
        return any(keyword.lower() in text_lower for keyword in keywords_list)

    # Aplicar la clasificación al título para obtener 'cancer_relation' y 'cancer_type'
    df[['cancer_relation', 'cancer_type']] = df[title_column].apply(
        lambda x: pd.Series(classify_cancer_type(x, cancer_types))
    )

    # Seleccionar la fuente de texto: se prefiere el abstract, sino se usa el título
    source_text = df[abstract_column] if abstract_column else df[title_column]

    # Calcular ambas relaciones: para TEM y para micro RNAs
    df['emt_relation'] = source_text.apply(lambda x: detect_relation(x, emt_keywords))
    df['mRNAs_relation'] = source_text.apply(lambda x: detect_relation(x, microrna_keywords))

    # Filtrar solo los artículos que están relacionados con cáncer
    df = df[df['cancer_relation'] == True]

    return df

# %%
import plotly.express as px
import pandas as pd

def graficar_distribucion_cancer_bibliografia_barras(df):

    df = df.assign(
        cancer_type=df['cancer_type'].str.split(', ')
    ).explode('cancer_type')

    dataset_cancer_counts = df.groupby(['nivel_expresion', 'cancer_type']).size().reset_index(name='study_count')

    # Ordenar los tipos de cáncer por frecuencia total
    cancer_order = dataset_cancer_counts.groupby('cancer_type')['study_count'].sum().sort_values(ascending=False).index
    dataset_cancer_counts['cancer_type'] = pd.Categorical(dataset_cancer_counts['cancer_type'],
                                                        categories=cancer_order,
                                                        ordered=True)

    color_discrete_map = {
        'Sobreexpresados': '#FBB4AE',
        'Subexpresados': '#B3CDE3'}

    # Crear gráfico interactivo con Plotly
    fig = px.bar(
        dataset_cancer_counts.sort_values('cancer_type'),
        x='cancer_type',
        y='study_count',
        color='nivel_expresion',
        # barmode p
        barmode='group',
        color_discrete_map=color_discrete_map,
        title=f'<b>Número de estudios por tipo de cáncer y conjuntos de expresión relacionados con {contexto_biologico}</b>',
        labels={
            'study_count': '<b>Número de estudios</b>',
            'cancer_type': '<b>Tipo de cáncer</b>',
            'dataset': '<b>Conjunto de datos</b>'
        },
        hover_data={'study_count': ':.0f'},
        height=600,
        width=1200
    )

    # Personalizar layout
    fig.update_layout(
        template='plotly_white',
        title_font_size=20,
        title_x=0.5,
        xaxis_tickangle=-45,
        hoverlabel=dict(
            bgcolor="white",
            font_size=12,
            font_family="Arial"
        ),
        legend=dict(
            title_text='<b>Conjunto de datos</b>',
            orientation="v",
            yanchor="bottom",
            y=0.8,
            xanchor="right",
            x=0.9
        ),
        margin=dict(l=50, r=50, b=150, t=100),
        xaxis=dict(
            tickfont=dict(size=12),
            title_standoff=25
        ),
        yaxis=dict(
            range=[0, dataset_cancer_counts['study_count'].max() * 1.15],
            gridcolor='lightgrey'
        )
    )

    # Mejorar tooltips
    fig.update_traces(
        hovertemplate="<br>".join([
            "<b>Tipo de cáncer:</b> %{x}",
            "<b>Conjunto de datos:</b> %{fullData.name}",
            "<b>Número de estudios:</b> %{y}"
        ])
    )

    guardar_grafico(fig, 'distribucion_cancer_bibliografia.png')
    fig.show()


# %%
import pandas as pd
import plotly.express as px

def graficar_distribucion_cancer_bibliografia(df: pd.DataFrame, selected_context) -> None:
    """
    Genera y muestra un gráfico interactivo que visualiza la distribución de artículos
    de bibliografía categorizados por tipo de cáncer y nivel de expresión,
    únicamente para los artículos relacionados con el contexto seleccionado.

    El DataFrame de entrada debe contener, al menos, las siguientes columnas:
        - 'cancer_type': Cadena de texto con tipos de cáncer, separados por comas y un espacio (", ").
        - Una columna indicadora de relación con el contexto, que será:
              * 'emt_relation' para "Cáncer y TEM"
              * 'mRNAs_relation' para "Cáncer y micro RNAs"
        - 'nivel_expresion': Nivel de expresión asociado a cada artículo.

    Parámetros:
        df (pd.DataFrame): DataFrame con la información de la bibliografía.
        selected_context (str): Contexto a utilizar ("Cáncer y TEM" o "Cáncer y micro RNAs").

    Ejemplo de uso:
        graficar_distribucion_cancer_bibliografia(mi_dataframe, selected_context="Cáncer y micro RNAs")
    """
    # Determinar la columna a utilizar y el nombre descriptivo del contexto
    if selected_context == "Cáncer y TEM":
        relation_column = "emt_relation"
        context_name = "Transición Epitelio Mesénquima (EMT)"
    elif selected_context == "Cáncer y micro RNAs":
        relation_column = "mRNAs_relation"
        context_name = "micro RNAs"
    else:
        raise ValueError(f"Contexto no reconocido: '{selected_context}'. Use 'Cáncer y TEM' o 'Cáncer y micro RNAs'.")

    # Verificar que el DataFrame contenga las columnas necesarias
    columnas_requeridas = {'cancer_type', relation_column, 'nivel_expresion'}
    columnas_faltantes = columnas_requeridas - set(df.columns)
    if columnas_faltantes:
        raise ValueError(f"El DataFrame no contiene las siguientes columnas requeridas para el contexto '{selected_context}': {columnas_faltantes}")

    # Filtrar solo los artículos relacionados con el contexto seleccionado
    df_filtrado = df[df[relation_column] == True]

    # Transformar el DataFrame: separar la columna 'cancer_type' en múltiples filas y agrupar para contar
    df_transformado = (
        df_filtrado.copy()
        .assign(cancer_type=lambda d: d['cancer_type'].str.split(', '))
        .explode('cancer_type')
        .groupby(['cancer_type', 'nivel_expresion'], as_index=False)
        .size()
        .rename(columns={'size': 'count'})
    )

    # Definir mapa de colores para cada nivel de expresión
    color_discrete_map = {
        "Sobreexpresados": "#FBB4AE",
        "Subexpresados": "#B3CDE3",
    }

    # Creación del gráfico de dispersión usando fstrings para incluir el contexto en el título
    fig = px.scatter(
        df_transformado,
        x='count',
        y='cancer_type',
        color='nivel_expresion',
        color_discrete_map=color_discrete_map,
        size='count',
        size_max=25,
        title=f"<b>Artículos relacionados con {context_name} y Nivel de Expresión</b>",
        labels={'count': 'Número de Artículos', 'cancer_type': 'Tipo de Cáncer'},
        height=600,
        width=1000,
        template='plotly_white'
    )

    # Ajuste del layout para mejorar la visualización
    fig.update_layout(
        yaxis={'categoryorder': 'total ascending'},
        hovermode='closest'
    )

    # Guardar el gráfico como imagen
    guardar_grafico(fig, f'distribucion_cancer_bibliografia_{contexto_biologico}.png')

    # Mostrar el gráfico interactivo
    fig.show()

# %% [markdown]
# # Búsqueda automatizada de bibliografía en PubMed para genes individuales

# %%
from Bio import Entrez, Medline
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor
import time
from urllib.error import HTTPError
import random
import os
import pandas as pd
import re

import ssl
ssl._create_default_https_context = ssl._create_unverified_context

# Configuración segura de la API (usar variables de entorno)
Entrez.email = os.getenv('NCBI_EMAIL', 'sergio.castelarf@gmail.com')
Entrez.api_key = os.getenv('NCBI_API_KEY', 'cbcc86f288ac31140c7e0978fd743cfb0808')

# Se asume que las siguientes variables globales están definidas:
# emt_keywords: lista de palabras clave para TEM (ej: ['EMT', 'transición epitelio mesénquima'])
# microrna_keywords: lista de palabras clave para micro RNAs (ej: ['microRNA', 'miRNA'])
# general_cancer_keywords: lista de palabras clave para cáncer (ej: ['cancer', 'tumor', 'carcinoma'])

def fetch_pubmed_data(gene_info, selected_context="Cáncer y TEM"):
    gene, ensembl_id = gene_info

    if selected_context == "Cáncer y TEM":
        context_keywords = emt_keywords
    elif selected_context == "Cáncer y micro RNAs":
        context_keywords = microrna_keywords
    else:
        raise ValueError(f"Contexto no reconocido: {selected_context}")

    combined_keywords = context_keywords + general_cancer_keywords
    query = f'("{gene}"[Title] OR "{ensembl_id}"[Title]) AND ("{" OR ".join(combined_keywords)}")'
    # Se elimina la impresión de la query para no saturar la salida

    article_data = []

    for attempt in range(10):
        try:
            handle = Entrez.esearch(db="pubmed", term=query, retmax=200)
            record = Entrez.read(handle)
            handle.close()
            article_ids = record["IdList"]

            if not article_ids:
                return []

            handle = Entrez.efetch(db="pubmed", id=article_ids, rettype="medline", retmode="text")
            records = Medline.parse(handle)

            for record in records:
                title = record.get('TI', 'No title').lower()
                abstract = record.get('AB', 'No abstract').lower()

                gene_in_title = re.search(
                    rf'\b{re.escape(gene.lower())}\b|\b{re.escape(ensembl_id.lower())}\b', title
                )
                if not gene_in_title:
                    continue

                year = None
                date_str = record.get('DP', '')
                if date_str:
                    year_match = re.search(r'\b(\d{4})\b', date_str)
                    if year_match:
                        year = int(year_match.group(1))

                text_for_search = f"{title} {abstract}"

                if selected_context == "Cáncer y TEM":
                    found_context = any(
                        re.search(rf'\b{re.escape(kw.lower())}\b', text_for_search) for kw in emt_keywords
                    )
                elif selected_context == "Cáncer y micro RNAs":
                    found_context = any(
                        re.search(rf'\b{re.escape(kw.lower())}\b', text_for_search) for kw in microrna_keywords
                    )

                found_cancer = any(
                    re.search(rf'\b{re.escape(kw.lower())}\b', text_for_search) for kw in general_cancer_keywords
                )

                if found_context or found_cancer:
                    if found_context and found_cancer:
                        category = "Ambos"
                    elif found_context:
                        category = "TEM" if selected_context == "Cáncer y TEM" else "micro RNAs"
                    elif found_cancer:
                        category = "Cáncer"
                else:
                    category = "Sin coincidencia"

                if category != "Sin coincidencia":
                    article_data.append({
                        "Gene": gene,
                        "Ensembl_ID": ensembl_id,
                        "Title": record.get('TI', 'No title'),
                        "Year": year,
                        "Abstract": record.get('AB', 'No abstract'),
                        "Link": f"https://pubmed.ncbi.nlm.nih.gov/{record.get('PMID', '')}",
                        "Keyword_Category": category
                    })

            handle.close()
            break

        # CAPTURA ESPECÍFICA DEL ERROR RuntimeError
        except RuntimeError as e:
            if "Database is not supported" in str(e):
                wait = random.uniform(1, 5)
                time.sleep(wait)
                print(f"Error temporal en PubMed, reintentando en {wait:.1f} segundos...")
            else:
                print(f"Error inesperado: {str(e)}")
                break
                
        except HTTPError as e:
            if e.code == 429:
                wait = random.uniform(1, 5)
                time.sleep(wait)
            else:
                print(f"Error HTTP ({e.code}): {str(e)}")
                break
                
        except Exception as e:
            print(f"Error inesperado: {type(e).__name__}: {str(e)}")
            break

    return article_data

def procesar_genes_y_guardar(df, selected_context="Cáncer y TEM"):
    genes = df['target'].tolist()
    ensembl_ids = df['ensembl_id'].tolist()

    with ThreadPoolExecutor(max_workers=2) as executor:
        results = list(tqdm(
            executor.map(lambda gene_info: fetch_pubmed_data(gene_info, selected_context), zip(genes, ensembl_ids)),
            total=len(genes),
            desc=f"Buscando artículos por gen para {selected_context}"
        ))

    articles = [art for sublist in results for art in sublist]
    if articles:
        df_resultado = pd.DataFrame(articles)
        df_resultado = df_resultado.drop_duplicates(subset=['Link'], keep='first')
        df_final = pd.merge(
            df_resultado,
            df[['target', 'ensembl_id', 'nivel_expresion']],
            left_on=['Gene', 'Ensembl_ID'],
            right_on=['target', 'ensembl_id']
        )
        df_final = df_final.drop(columns=['target', 'ensembl_id'])
    else:
        df_final = pd.DataFrame()

    print(f"Búsqueda completada para {selected_context}. Artículos encontrados: {len(df_final)}")
    return df_final

# Uso del código
df_resultado_final = procesar_genes_y_guardar(df_emsamble_info, contexto_biologico)
success_message(f"Búsqueda de artículos en PubMed completada con éxito para el contexto: {contexto_biologico}")

# %% [markdown]
# ## Clasificación de la bibliografía

# %% [markdown]
# Aplicamos la función para clasificar la bibliografía encontrada y la asigna al tipo o tipos de cancer mencionados en cada articulo

# %%
try:
    # Intentar generar el DataFrame normalmente
    bibliografia_genes_clasificada = clasificar_tipos_cancer_bibliografia(df_resultado_final)
except Exception as e:
    print(f"Error al generar el DataFrame, cargando respaldo desde Excel: {e}")
    bibliografia_genes_clasificada = pd.read_excel("../bibliografia_clasif.xlsx")

# Exportar el DataFrame (sea generado o cargado)
export_dfs_to_excel(
    [bibliografia_genes_clasificada],
    ["bibliografia_genes_clasificada"],
    "Bibliografía de genes clasificada: Archivo que contiene los resultados de la búsqueda en PubMed por genes, con información adicional sobre la clasificación de tipos de cáncer y relación con EMT."
)

# %% [markdown]
# ## Grafico de la distribución de los tipos de cáncer con respecto a la cantidad de bibligrafia separados en grupos de expresión

# %%
graficar_distribucion_cancer_bibliografia_barras(bibliografia_genes_clasificada)

# %%
graficar_distribucion_cancer_bibliografia(bibliografia_genes_clasificada, contexto_biologico)