# %% [markdown]
# # Creación de firmas
# 
# ## Etapa 1: Firmas analíticas

# %% [markdown]
# ### Generamos las firmas de la etapa 1
# 

# %%
import pandas as pd
import numpy as np
import re
import datetime
import gseapy as gp



def hallmarks_enrichment(lista_de_expresion, background_genes):
    """
    Función externa de enriquecimiento
    """
    hallmark_results = gp.enrich(
        gene_list=lista_de_expresion,
        gene_sets="./gen-sets_GSEA_MSigDB/gsea_hallmarks_formatted.gmt",
        background=background_genes,
        outdir=None,
        verbose=False
    )

    # Crear un diccionario plano con términos como claves
    hallmark_dict = {}
    if not hallmark_results.results.empty:
        for _, row in hallmark_results.results.iterrows():
            term = row["Term"].replace(" ", "_")  # Ej: HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
            hallmark_dict[f"{term}_pvalue"] = row["Adjusted P-value"]

    return hallmark_dict

def crear_firma_genetica(df):
    """
    Genera una firma genética con características listas para Machine Learning,
    incluyendo las puntuaciones de Hallmarks y una nueva métrica robusta de
    significancia bibliográfica que integra:
      - Frecuencia relativa del gen en el cáncer versus global.
      - Número de artículos (suavizbado con logaritmo).
      - Ponderación por recencia de los artículos.
    """
    import re
    import datetime

    # Copia del DataFrame para no modificar el original
    df = df.copy()

    # 1. Extracción del PubMed ID a partir de la URL
    def extraer_pubmed_id(link):
        match = re.search(r'/(\d+)$', str(link))
        return match.group(1) if match else None
    df['pubmed_id'] = df['Link'].apply(extraer_pubmed_id)

    # 2. Procesamiento del campo 'cancer_type'
    df['cancer_type'] = df['cancer_type'].astype(str).str.split(',').apply(lambda x: [i.strip() for i in x])
    df = df.explode('cancer_type')

    # 3. Carga del conjunto de fondo (background) desde el archivo GMT
    background_genes = set()
    with open("../gen-sets_GSEA_MSigDB/C5- ontology gene sets.gmt", "r") as file:
        for line in file:
            genes = line.strip().split("\t")[2:]
            background_genes.update(genes)
    background_genes = sorted(background_genes)

    # 4. Cálculos globales:
    # - Conteo total de artículos por tipo de cáncer y nivel de expresión.
    cancer_total_counts = df.groupby(['cancer_type', 'nivel_expresion']).size().to_dict()

    # - Conteo de artículos relacionados con el contexto biológico,
    #   separados por niveles de expresión.
    if contexto_biologico == "Cáncer y TEM":
        biological_context = df[df['emt_relation'] == True]
        cancer_biological_context_counts = biological_context.groupby(['cancer_type', 'nivel_expresion']).size().to_dict()
    elif contexto_biologico == "Cáncer y micro RNAs":
        biological_context = df[df['mRNAs_relation'] == True]
        cancer_biological_context_counts = biological_context.groupby(['cancer_type', 'nivel_expresion']).size().to_dict()

    # Total de artículos en el DataFrame
    total_articulos = len(df)

    # Conteo global de artículos por gen en todo el DataFrame
    global_gene_counts = df['Gene'].value_counts().to_dict()

    # Parámetro de decaimiento para ponderar la recencia y año actual
    decay = 0.05
    current_year = datetime.datetime.now().year


    # 5. Función interna para el enriquecimiento de Hallmarks con manejo de excepciones
    # 5. Función interna para el enriquecimiento de Hallmarks con manejo de excepciones
    def hallmarks_enrichment_local(gene_list):
        """Devuelve un diccionario plano con características por Hallmark."""
        if not gene_list:
            return {}
        try:
            enr = gp.enrich(
                gene_list=gene_list,
                gene_sets="../gen-sets_GSEA_MSigDB/gsea_hallmarks_formatted.gmt",
                background=background_genes,
                outdir=None,
                verbose=False
            )

            hallmark_features = {}
            for _, row in enr.results.iterrows():
                term = row['Term'].replace(' ', '_').replace('HALLMARK_', '')
                prefix = f"hallmark_{term}"
                genes_overlap = row['Genes'].split(';')
                # Crear una entrada única para los genes de cada hallmark
                hallmark_features[f"{prefix}_genes"] = genes_overlap  # Clave dinámica
                hallmark_features[f"{prefix}_pvalue"] = row['Adjusted P-value']
            return hallmark_features
        except Exception as e:
            print(f"Error en enriquecimiento: {str(e)}")
            return {}

    # 6. Procesamiento por grupos (agrupando por 'cancer_type' y 'nivel_expresion')
    result_rows = []
    grouped = df.groupby(['cancer_type', 'nivel_expresion'])

    for (cancer, nivel), group in grouped:
        # Genes únicos en el grupo
        genes_unicos = group['Gene'].unique().tolist()

        # Inicializar listas para almacenar resultados
        genes_list = []
        conteo_list = []
        score_list = []
        pubmed_dict = {}

        for gene in genes_unicos:
            gene_df = group[group['Gene'] == gene]
            n_gene = gene_df.shape[0]

            # Registrar gene y su conteo en el grupo
            genes_list.append(gene)
            conteo_list.append(n_gene)

            # Total de artículos para este tipo de cáncer
            total_cancer = cancer_total_counts.get(cancer, 1)
            # Conteo global para el gen en todo el DataFrame
            n_gene_global = global_gene_counts.get(gene, 1)

            # Frecuencia relativa:
            # (n_gene / total_cancer) en este cáncer frente a (n_gene_global / total_articulos) en general
            freq_ratio = (n_gene / total_cancer) / (n_gene_global / total_articulos)

            # Factor de recencia: promedio del decaimiento exponencial basado en la diferencia de años
            recency_values = np.exp(-decay * (current_year - gene_df['Year'].astype(float)))
            recency_factor = recency_values.mean()

            # Métrica robusta combinada
            gene_score = freq_ratio * np.log(1 + n_gene) * recency_factor
            score_list.append(gene_score)

            # Recopilar PubMed IDs asociados a este gen
            pubmed_dict[gene] = gene_df['pubmed_id'].dropna().unique().tolist()

        # Enriquecimiento de Hallmarks para la lista de genes del grupo
        hallmark_features = hallmarks_enrichment_local(genes_list)

        # Construir la fila de resultados para el grupo
        result_row = {
            'cancer_type': cancer,
            'nivel_expresion': nivel,
            'genes': genes_list,
            'conteo_articulos_por_gene': conteo_list,
            'estadistica_significancia_bibliografica': score_list,
            'conteo_global_contexto_biologico': cancer_biological_context_counts.get(cancer, 0),
            'pubmed_ids': pubmed_dict,
            **hallmark_features  # Se integran las características obtenidas del enriquecimiento
        }

        result_rows.append(result_row)

    # 7. Creación del DataFrame final
    result_df = pd.DataFrame(result_rows)

    # Rellenar NaN para las columnas de hallmarks (importante para ML)
    hallmark_cols = [col for col in result_df.columns if col.startswith('hallmark_')]
    result_df[hallmark_cols] = result_df[hallmark_cols].fillna(0)

    return result_df

# %%
df_firmas_geneticas = crear_firma_genetica(bibliografia_genes_clasificada)

export_dfs_to_excel(
    [df_firmas_geneticas],
    ["firmas_por_tipo"],
    "Firmas genéticas por tipo de cáncer: Archivo que contiene la firma genética por tipo de cáncer y nivel de expresión, con información sobre los genes asociados, estadísticas y PubMed IDs."
)


# %% [markdown]
# ### Gráfico de distribución de los hallmarks de las firmas por tipo de cáncer y grupo de expresión

# %%
import pandas as pd
import plotly.express as px
import numpy as np
from typing import List, Dict



def generar_sunburst_para_tipo_cancer(
    df_firmas_geneticas: pd.DataFrame,
    hallmarks: List[str] = [
        "EVADING_GROWTH_SUPPRESSORS",
        "EVADING_IMMUNE_DESTRUCTION",
        "GENOME_INSTABILITY",
        "REPLICATIVE_IMMORTALITY",
        "REPROGRAMMING_ENERGY_METABOLISM",
        "RESISTING_CELL_DEATH",
        "SUSTAINED_ANGIOGENESIS",
        "SUSTAINING_PROLIFERATIVE_SIGNALING",
        "TISSUE_INVASION_AND_METASTASIS",
        "TUMOR-PROMOTING_INFLAMMATION"
    ],
    fig_width: int = 1200,
    fig_height: int = 900,
    color_scale: str = 'RdBu_r'
) -> Dict[str, px.sunburst]:

    # Configuraciones
    NON_SIG_COLOR = '#4a4a4a'
    PVAL_THRESHOLD = 0.05
    LOG_THRESHOLD = -np.log10(PVAL_THRESHOLD)

    # Diccionario de alias
    hallmark_aliases = {
        "EVADING_GROWTH_SUPPRESSORS": "Evade Growth Suppression",
        "EVADING_IMMUNE_DESTRUCTION": "Evade Immune Destruction",
        "GENOME_INSTABILITY": "Genome Instability",
        "REPLICATIVE_IMMORTALITY": "Replicative Immortality",
        "REPROGRAMMING_ENERGY_METABOLISM": "Metabolic Reprogramming",
        "RESISTING_CELL_DEATH": "Resist Cell Death",
        "SUSTAINED_ANGIOGENESIS": "Sustain Angiogenesis",
        "SUSTAINING_PROLIFERATIVE_SIGNALING": "Sustain Proliferation",
        "TISSUE_INVASION_AND_METASTASIS": "Invasion & Metastasis",
        "TUMOR-PROMOTING_INFLAMMATION": "Tumor Inflammation"
    }

    # Procesamiento de datos
    registros = []
    for _, row in df_firmas_geneticas.iterrows():
        cancer = row.get("cancer_type")
        nivel = row.get("nivel_expresion")
        genes_firma = row.get("genes", [])
        counts_firma = row.get("conteo_articulos_por_gene", [])

        # Solo procesar el tipo de cáncer especificado
        if cancer != tipo_cancer:
            continue

        if not (isinstance(genes_firma, list) and isinstance(counts_firma, list) and len(genes_firma) == len(counts_firma)):
            continue

        gene_to_count = dict(zip(genes_firma, counts_firma))

        for gene in genes_firma:
            articles_gene = gene_to_count.get(gene, 0)

            for hallmark in hallmarks:
                col_genes = f"hallmark_{hallmark}_genes"
                col_pvalue = f"hallmark_{hallmark}_pvalue"

                if col_genes in row and col_pvalue in row:
                    genes_hallmark = row[col_genes] if isinstance(row[col_genes], list) else []
                    pvalue = row[col_pvalue]

                    if gene in genes_hallmark:
                        registros.append({
                            "cancer_type": cancer,
                            "nivel_expresion": nivel,
                            "gene": gene,
                            "hallmark": hallmark_aliases.get(hallmark, hallmark),
                            "articles": articles_gene,
                            "pvalue": pvalue
                        })

    if not registros:
        raise ValueError(f"No se encontraron registros válidos para {tipo_cancer}.")

    df_plano = pd.DataFrame(registros)
    df_plano["log_pvalue"] = -np.log10(df_plano["pvalue"].replace(0, 1e-300))

    # Funciones auxiliares
    def create_custom_colorscale(max_log_p: float) -> list:
        """Crea la escala de color personalizada"""
        if max_log_p <= LOG_THRESHOLD:
            return [[0, NON_SIG_COLOR], [1, NON_SIG_COLOR]]

        threshold_pos = LOG_THRESHOLD / max_log_p
        RdBu_r_colors = px.colors.diverging.RdBu_r

        return [
            [0, NON_SIG_COLOR],
            [threshold_pos, NON_SIG_COLOR],
            [threshold_pos, RdBu_r_colors[0]],
            *[[threshold_pos + (i/(len(RdBu_r_colors)-1)*(1-threshold_pos)), color]
              for i, color in enumerate(RdBu_r_colors)],
            [1, RdBu_r_colors[-1]]
        ]

    def add_custom_annotations(fig):
        """Añade las anotaciones al gráfico"""
        fig.add_annotation(
            x=0.98, y=0.05,
            xref="paper", yref="paper",
            text=(
                "<b>📚 Tamaño: Evidencia bibliográfica<br>"
                "🎨 Color: Significancia estadística<br>"
                "⚪ p ≥ 0.05 (No significativo)</b>"
            ),
            showarrow=False,
            font=dict(size=11, color="#404040"),
            align="right",
            bordercolor="#a0a0a0",
            borderwidth=1,
            bgcolor="#ffffff"
        )

    # Generación de gráficos
    figuras_por_nivel = {}
    for nivel in df_plano["nivel_expresion"].unique():
        df_nivel = df_plano[df_plano["nivel_expresion"] == nivel].copy()

        max_log_p = df_nivel["log_pvalue"].max()
        color_scale_custom = create_custom_colorscale(max_log_p)

        fig = px.sunburst(
            df_nivel,
            path=['gene', 'hallmark'],  # Eliminamos cancer_type ya que solo hay uno
            values='articles',
            color='log_pvalue',
            color_continuous_scale=color_scale_custom,
            range_color=[0, max_log_p],
            title=f"Sunburst - {tipo_cancer} - Nivel: {nivel}",
            hover_data={'pvalue': True, 'log_pvalue': True}
        )

        add_custom_annotations(fig)

        fig.update_layout(
            width=fig_width,
            height=fig_height,
            margin=dict(t=60, b=20, l=20, r=20),
            coloraxis_colorbar=dict(
                title="<b>-log10(p)</b>",
                thickness=15,
                yanchor="middle",
                y=0.5
            )
        )

        figuras_por_nivel[nivel] = fig

    return figuras_por_nivel

# Ejecución
figs_por_nivel = generar_sunburst_para_tipo_cancer(df_firmas_geneticas)
for nivel, fig in figs_por_nivel.items():
    guardar_grafico(fig, f"sunburst_{tipo_cancer.lower().replace(' ', '_')}_nivel_{nivel}.png")
    fig.show()

# %%
import pandas as pd
import plotly.express as px
import numpy as np
from typing import List, Dict

def generar_sunburst_top5_canceres_por_nivel(
    df_firmas_geneticas: pd.DataFrame,
    hallmarks: List[str] = [
        "EVADING_GROWTH_SUPPRESSORS",
        "EVADING_IMMUNE_DESTRUCTION",
        "GENOME_INSTABILITY",
        "REPLICATIVE_IMMORTALITY",
        "REPROGRAMMING_ENERGY_METABOLISM",
        "RESISTING_CELL_DEATH",
        "SUSTAINED_ANGIOGENESIS",
        "SUSTAINING_PROLIFERATIVE_SIGNALING",
        "TISSUE_INVASION_AND_METASTASIS",
        "TUMOR-PROMOTING_INFLAMMATION"
    ],
    fig_width: int = 1200,
    fig_height: int = 900,
    color_scale: str = 'RdBu_r'
) -> Dict[str, px.sunburst]:

    NON_SIG_COLOR = '#4a4a4a'
    PVAL_THRESHOLD = 0.05
    LOG_THRESHOLD = -np.log10(PVAL_THRESHOLD)

    hallmark_aliases = {
        "EVADING_GROWTH_SUPPRESSORS": "Evade Growth Suppression",
        "EVADING_IMMUNE_DESTRUCTION": "Evade Immune Destruction",
        "GENOME_INSTABILITY": "Genome Instability",
        "REPLICATIVE_IMMORTALITY": "Replicative Immortality",
        "REPROGRAMMING_ENERGY_METABOLISM": "Metabolic Reprogramming",
        "RESISTING_CELL_DEATH": "Resist Cell Death",
        "SUSTAINED_ANGIOGENESIS": "Sustain Angiogenesis",
        "SUSTAINING_PROLIFERATIVE_SIGNALING": "Sustain Proliferation",
        "TISSUE_INVASION_AND_METASTASIS": "Invasion & Metastasis",
        "TUMOR-PROMOTING_INFLAMMATION": "Tumor Inflammation"
    }

    registros = []
    for _, row in df_firmas_geneticas.iterrows():
        cancer = row.get("cancer_type")
        nivel = row.get("nivel_expresion")
        genes_firma = row.get("genes", [])
        counts_firma = row.get("conteo_articulos_por_gene", [])

        if not (isinstance(genes_firma, list) and isinstance(counts_firma, list) and len(genes_firma) == len(counts_firma)):
            continue

        gene_to_count = dict(zip(genes_firma, counts_firma))

        for gene in genes_firma:
            articles_gene = gene_to_count.get(gene, 0)

            for hallmark in hallmarks:
                col_genes = f"hallmark_{hallmark}_genes"
                col_pvalue = f"hallmark_{hallmark}_pvalue"

                if col_genes in row and col_pvalue in row:
                    genes_hallmark = row[col_genes] if isinstance(row[col_genes], list) else []
                    pvalue = row[col_pvalue]

                    if gene in genes_hallmark:
                        registros.append({
                            "cancer_type": cancer,
                            "nivel_expresion": nivel,
                            "gene": gene,
                            "hallmark": hallmark_aliases.get(hallmark, hallmark),
                            "articles": articles_gene,
                            "pvalue": pvalue
                        })

    if not registros:
        raise ValueError("No se encontraron registros válidos en el DataFrame.")

    df_plano = pd.DataFrame(registros)
    df_plano["log_pvalue"] = -np.log10(df_plano["pvalue"].replace(0, 1e-300))

    # Seleccionar top 5 cánceres más frecuentes
    top5_canceres = (
        df_plano['cancer_type']
        .value_counts()
        .nlargest(6)
        .index.tolist()
    )
    df_plano = df_plano[df_plano['cancer_type'].isin(top5_canceres)]

    def create_custom_colorscale(max_log_p: float) -> list:
        if max_log_p <= LOG_THRESHOLD:
            return [[0, NON_SIG_COLOR], [1, NON_SIG_COLOR]]

        threshold_pos = LOG_THRESHOLD / max_log_p
        RdBu_r_colors = px.colors.diverging.RdBu_r

        return [
            [0, NON_SIG_COLOR],
            [threshold_pos, NON_SIG_COLOR],
            [threshold_pos, RdBu_r_colors[0]],
            *[[threshold_pos + (i/(len(RdBu_r_colors)-1)*(1-threshold_pos)), color]
              for i, color in enumerate(RdBu_r_colors)],
            [1, RdBu_r_colors[-1]]
        ]

    def add_custom_annotations(fig):
        fig.add_annotation(
            x=0.98, y=0.05,
            xref="paper", yref="paper",
            text=(
                "<b>📚 Tamaño: Evidencia bibliográfica<br>"
                "🎨 Color: Significancia estadística<br>"
                "⚪ p ≥ 0.05 (No significativo)</b>"
            ),
            showarrow=False,
            font=dict(size=11, color="#404040"),
            align="right",
            bordercolor="#a0a0a0",
            borderwidth=1,
            bgcolor="#ffffff"
        )

    figuras_por_nivel = {}
    for nivel in df_plano["nivel_expresion"].unique():
        df_nivel = df_plano[df_plano["nivel_expresion"] == nivel].copy()

        max_log_p = df_nivel["log_pvalue"].max()
        color_scale_custom = create_custom_colorscale(max_log_p)

        fig = px.sunburst(
            df_nivel,
            path=['cancer_type', 'gene', 'hallmark'],
            values='articles',
            color='log_pvalue',
            color_continuous_scale=color_scale_custom,
            range_color=[0, max_log_p],
            title=f"Sunburst – Nivel de expresión: {nivel} (Top 5 Cánceres)",
            hover_data={'pvalue': True, 'log_pvalue': True}
        )

        add_custom_annotations(fig)

        fig.update_layout(
            width=fig_width,
            height=fig_height,
            margin=dict(t=60, b=20, l=20, r=20),
            coloraxis_colorbar=dict(
                title="<b>-log10(p)</b>",
                thickness=15,
                yanchor="middle",
                y=0.5
            )
        )

        figuras_por_nivel[nivel] = fig

    return figuras_por_nivel

# Ejemplo de ejecución:
figs_por_nivel = generar_sunburst_top5_canceres_por_nivel(df_firmas_geneticas)

for nivel, fig in figs_por_nivel.items():
    fig.show()

# %%
# print(len(lista_sobre))

# %%
import pandas as pd
import pandas as pd
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px


# Diccionario de alias
hallmark_aliases = {
    "EVADING_GROWTH_SUPPRESSORS": "Evade Growth Suppression",
    "EVADING_IMMUNE_DESTRUCTION": "Evade Immune Destruction",
    "GENOME_INSTABILITY": "Genome Instability",
    "REPLICATIVE_IMMORTALITY": "Replicative Immortality",
    "REPROGRAMMING_ENERGY_METABOLISM": "Metabolic Reprogramming",
    "RESISTING_CELL_DEATH": "Resist Cell Death",
    "SUSTAINED_ANGIOGENESIS": "Sustain Angiogenesis",
    "SUSTAINING_PROLIFERATIVE_SIGNALING": "Sustain Proliferation",
    "TISSUE_INVASION_AND_METASTASIS": "Invasion & Metastasis",
    "TUMOR-PROMOTING_INFLAMMATION": "Tumor Inflammation"
}

# Función para generar el gráfico
def generar_grafico(cancer_seleccionado, nivel_expresion):
    # Filtrar datos
    filtro = (df_firmas_geneticas['cancer_type'] == cancer_seleccionado) & \
             (df_firmas_geneticas['nivel_expresion'] == nivel_expresion)

    if df_firmas_geneticas[filtro].empty:
        return go.Figure()

    datos = df_firmas_geneticas[filtro].iloc[0]

    # Procesar p-values
    pvalue_cols = [c for c in datos.index if c.startswith('hallmark_') and c.endswith('_pvalue')]
    hallmarks = ['_'.join(c.split('_')[1:-1]) for c in pvalue_cols]
    pvalues = datos[pvalue_cols].tolist()

    # Transformar valores
    transformed_values = [-np.log10(p) if p > 0 else 0 for p in pvalues]
    significant = [p < 0.05 for p in pvalues]

# Configurar colores con paleta Pastel1
    colores = []
    palette = px.colors.sequential.Aggrnyl
    color_idx = 0
    for sig in significant:
        if sig:
            colores.append(palette[color_idx % len(palette)])
            color_idx += 1
        else:
            colores.append('#e0e0e0')  # Gris pastel

    # Crear figura
    fig = go.Figure()

    # Añadir barras traslúcidas
    fig.add_trace(go.Barpolar(
        r=transformed_values,
        theta=[hallmark_aliases[h] for h in hallmarks],
        marker_color=[c.replace('rgb', 'rgba').replace(')', ',0.4)') for c in colores],
        marker_line_color='rgba(100,100,100,0.5)',
        marker_line_width=1,
        hoverinfo='text',
        hovertext=[f"<b>{hallmark_aliases[h]}</b><br>-log10(p): {v:.2f}<br>p-value: {p:.4e}"
                  for h, v, p in zip(hallmarks, transformed_values, pvalues)]
    ))

    # Añadir marcadores circulares encima
    fig.add_trace(go.Scatterpolar(
        r=transformed_values,
        theta=[hallmark_aliases[h] for h in hallmarks],
        mode='markers',
        marker=dict(
            color=colores,
            size=6,
            symbol='circle',
            line=dict(color='rgba(0,0,0,0.5)', width=1)
        ),
        hoverinfo='none'
    ))

    # Diseño mejorado
    fig.update_layout(
        title= f"Hallmarks en {cancer_seleccionado} - {nivel_expresion}",
        polar={
            'bgcolor': 'white',
            'radialaxis': {
                'title': {'text': '<b>-log10(p-value)</b>', 'font': {'size': 14, 'color':'#000'}},
                'gridcolor': '#dfe6e9',
                'linecolor': '#b2bec3',
                'tickfont': {'color': '#000'},
                'range': [0, max(transformed_values)*1.2]
            },
            'angularaxis': {
                'rotation': 90,
                'direction': 'clockwise',
                'gridcolor': '#dfe6e9',
                'linecolor': '#b2bec3',
                'tickfont': {'size': 12, 'color':'#000'}
            }
        },
        paper_bgcolor='white',
        plot_bgcolor='white',
        height=700,
        margin={'t': 100, 'b': 80, 'l': 80, 'r': 80},
        showlegend=False,
        annotations=[
            dict(
                x=0.98,
                y=0.98,
                xref='paper',
                yref='paper',
                text='<b style="color:#2c3e50">● Hallmarks significativos (p < 0.05)</b>',
                showarrow=False,
                bgcolor='rgba(255,255,255,0.9)',
                bordercolor='#dfe6e9',
                borderwidth=1,
                font={'size': 12}
            )
        ]
    )

    return fig

from IPython.display import display

# Lista de combinaciones que queremos graficar

niveles_expresion = ['Sobreexpresados', 'Subexpresados']

# Mostrar gráficos para cada combinación

for nivel in niveles_expresion:
    fig = generar_grafico(tipo_cancer, nivel)
    guardar_grafico(fig, f"hallmarks_{tipo_cancer}_{nivel}.png")
    display(fig)

# %%


# %%
def preparar_y_visualizar_fingerprint(df, tipo_cancer):
    """
    """
    import ast
    import pandas as pd
    import numpy as np
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    import plotly.figure_factory as ff
    from scipy.cluster.hierarchy import linkage, leaves_list

    # Función de conversión segura
    def convertir_a_lista_si_es_str(x):
        if isinstance(x, str) and x.startswith("[") and x.endswith("]"):
            try:
                return ast.literal_eval(x)
            except Exception:
                return x
        return x

    # Limpieza de datos
    df = df.copy()
    df["cancer_type"] = df["cancer_type"].astype(str).str.strip()
    df["nivel_expresion"] = df["nivel_expresion"].astype(str).str.strip()
    df["genes"] = df["genes"].apply(convertir_a_lista_si_es_str)
    
    cols_genes = [col for col in df.columns if col.startswith("hallmark_") and col.endswith("_genes")]
    for col in cols_genes:
        df[col] = df[col].apply(convertir_a_lista_si_es_str)

    pval_cols = [col for col in df.columns if col.startswith("hallmark_") and col.endswith("_pvalue")]

    for nivel_expresion in ["Sobreexpresados", "Subexpresados"]:
        firma = df[
            (df["cancer_type"] == tipo_cancer.strip()) &
            (df["nivel_expresion"] == nivel_expresion.strip())
        ]

        if firma.empty:
            print(f"❌ No se encontraron datos para: {tipo_cancer} - {nivel_expresion}")
            continue

        fila = firma.iloc[0]
        genes_firma = fila["genes"]

        registros = []
        for gcol in cols_genes:
            hallmark_base = gcol.replace("hallmark_", "").replace("_genes", "")
            pval_col = f"hallmark_{hallmark_base}_pvalue"
            hallmark_name = hallmark_base.replace("_", " ").title()

            genes_hallmark = fila[gcol]
            pval = fila[pval_col] if fila[pval_col] > 0 else 1.0
            logp = -np.log10(pval)

            for gen in genes_firma:
                if isinstance(genes_hallmark, list) and gen in genes_hallmark:
                    registros.append({
                        "gen": gen,
                        "hallmark": hallmark_name,
                        "-log10(p)": logp,
                        "p-value": 10**-logp
                    })

        df_plot = pd.DataFrame(registros)
        if df_plot.empty:
            print(f"⚠️ No hay genes coincidentes para: {tipo_cancer} - {nivel_expresion}")
            continue

        # Crear matriz y aplicar clustering
        matriz = df_plot.pivot(index="gen", columns="hallmark", values="-log10(p)")
        matriz_filled = matriz.fillna(0)
        
        # Clustering jerárquico
        if len(matriz) > 1:
            row_linkage = linkage(matriz_filled, method='ward')
            row_order = leaves_list(row_linkage)
            matriz = matriz.iloc[row_order[::-1]]  # Invertir para mejor visualización

        if len(matriz.columns) > 1:
            col_linkage = linkage(matriz_filled.T, method='ward')
            col_order = leaves_list(col_linkage)
            matriz = matriz.iloc[:, col_order]

        # Generar texto para hover
        hover_text = np.array([[
            f"<b>Gen:</b> {gen}<br><b>Hallmark:</b> {hm}<br>"
            f"<b>-log10(p):</b> {matriz.loc[gen, hm]:.2f}<br>"
            f"<b>p-value:</b> {10**-matriz.loc[gen, hm]:.2e}"
            if not np.isnan(matriz.loc[gen, hm]) else ""
            for hm in matriz.columns
        ] for gen in matriz.index])

        # Crear visualización interactiva
        fig = go.Figure(go.Heatmap(
            z=matriz.values,
            x=matriz.columns,
            y=matriz.index,
            colorscale="Viridis",
            colorbar={"title": "-log10(p)", "tickvals": np.linspace(0, matriz.max().max(), 5)},
            hoverinfo="text",
            text=hover_text,
            texttemplate="",
            zmin=0,
            zmid=matriz.max().max()/2,
            zmax=matriz.max().max()
        ))

        # Ajustes de diseño profesional
        fig.update_layout(
            title={
                'text': f"<b>Fingerprint Genómico</b><br>{tipo_cancer} ({nivel_expresion})",
                'y':0.95,
                'x':0.5,
                'xanchor': 'center',
                'font': {'size': 24}
            },
            xaxis={
                'title': 'Hallmarks',
                'tickangle': 45,
                'tickfont': {'size': 12},
                'automargin': True
            },
            yaxis={
                'title': 'Genes',
                'tickfont': {'size': 12},
                'automargin': True
            },
            margin={'l': 150, 'r': 50, 't': 150, 'b': 150},
            width=900,
            height=700,
            hoverlabel={
                'bgcolor': 'white',
                'font_size': 14,
                'font_family': 'Arial'
            },
            plot_bgcolor='rgba(245,245,245,0.9)'
        )

        # Añadir interactividad avanzada
        fig.update_xaxes(fixedrange=False)
        fig.update_yaxes(fixedrange=False)
        
        # Mostrar y guardar
        fig.show()
        guardar_grafico(fig, f"fingerprint_mejorado_{tipo_cancer.replace(' ', '_')}_{nivel_expresion}.png")

# %%
# Ejecuta para cualquier tipo de cáncer
preparar_y_visualizar_fingerprint(df_firmas_geneticas, "Breast Cancer")

# %%


# %%


# %% [markdown]
# ### Niveles de expresión relativa
# 
# 

# %%
import pandas as pd
import plotly.express as px
import plotly.figure_factory as ff
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# ------------------------------
# 1. Filtrar firmas genéticas y extraer las listas
# ------------------------------
# Se asume que 'df_firmas_geneticas' ya está definido y contiene:
# 'cancer_type', 'nivel_expresion' y 'genes'

# Definir el tipo de cáncer a analizar (por ejemplo, 'Breast Cancer')
df_cancer = df_firmas_geneticas[df_firmas_geneticas['cancer_type'] == tipo_cancer]

# Función para asegurar que 'genes' se procese como lista
def ensure_list(val):
    if isinstance(val, list):
        return [item.strip() for item in val]
    elif isinstance(val, str):
        return [item.strip() for item in val.split(',')]
    else:
        return list(val)

# Extraer genes según el nivel de expresión
cadena_sobre = df_cancer.loc[df_cancer['nivel_expresion'] == 'Sobreexpresados', 'genes'].iloc[0]
cadena_sub   = df_cancer.loc[df_cancer['nivel_expresion'] == 'Subexpresados', 'genes'].iloc[0]

lista_sobre = ensure_list(cadena_sobre)
lista_sub   = ensure_list(cadena_sub)

print(f"Genes sobreexpresados: {lista_sobre}")
print(f"Genes subexpresados: {lista_sub}")

# ------------------------------
# 2. Concatenar los DataFrames de expresión (controles y muestras)
# ------------------------------
# Se asume que 'controles_df' y 'muestras_df' ya están definidos y contienen:
# 'test', 'tipo', 'target' y 'ct'
df_total = pd.concat([controles_df, muestras_df])

# ------------------------------
# 3. Normalización de datos según el método elegido
# ------------------------------

if metodo_elegido == "gen de referencia":
    # Se extrae el Ct del gen de referencia para cada test
    ref_ct = df_total[df_total['target'] == gen_referencia][['test', 'ct']].rename(
        columns={'ct': f'{gen_referencia}_ct'}
    )
    df_norm = df_total.merge(ref_ct, on='test')
    df_norm['delta_ct'] = df_norm['ct'] - df_norm[f'{gen_referencia}_ct']
    df_norm['log2_rel_expr'] = -df_norm['delta_ct']  # log2(2^(-ΔCt)) = -ΔCt

elif metodo_elegido == "promedios":
    # Se calcula el promedio de Ct para cada test y se normaliza respecto a ese promedio
    df_norm = df_total.copy()
    df_norm['ct_promedio_test'] = df_norm.groupby('test')['ct'].transform('mean')
    df_norm['delta_ct_promedio'] = df_norm['ct'] - df_norm['ct_promedio_test']
    df_norm['log2_rel_expr'] = -df_norm['delta_ct_promedio']

else:
    raise ValueError("El método de normalización debe ser 'gen de referencia' o 'promedios'.")

# Crear la matriz para el heatmap (tabla pivote)
# Si se normaliza por gen de referencia, se elimina dicho gen de la matriz
if metodo_elegido == "gen de referencia":
    heatmap_data = df_norm.pivot_table(
        index='target',
        columns='test',
        values='log2_rel_expr',
        aggfunc='mean'
    ).drop(gen_referencia, errors='ignore')
else:
    heatmap_data = df_norm.pivot_table(
        index='target',
        columns='test',
        values='log2_rel_expr',
        aggfunc='mean'
    )

# ------------------------------
# 4. Filtrar la matriz para cada nivel de expresión
# ------------------------------
heatmap_data_sobre = heatmap_data[heatmap_data.index.isin(lista_sobre)]
heatmap_data_sub   = heatmap_data[heatmap_data.index.isin(lista_sub)]

# ------------------------------
# 5. Crear el clustergrama (heatmap con dendrogramas)
# ------------------------------
from scipy.cluster.hierarchy import linkage
import plotly.figure_factory as ff
from plotly.subplots import make_subplots
import plotly.graph_objects as go

def create_clustered_heatmap(heatmap_data, title):
    # Clustering y dendrogramas para columnas (tests)
    dendro_col = ff.create_dendrogram(
        heatmap_data.T,
        orientation='bottom',
        labels=heatmap_data.columns.tolist(),
        linkagefun=lambda x: linkage(x, method='average', metric='euclidean')
    )
    col_order = dendro_col['layout']['xaxis']['ticktext']

    # Clustering y dendrogramas para filas (genes)
    dendro_row = ff.create_dendrogram(
        heatmap_data,
        orientation='right',
        labels=heatmap_data.index.tolist(),
        linkagefun=lambda x: linkage(x, method='average', metric='euclidean')
    )
    row_order = dendro_row['layout']['yaxis']['ticktext']

    # Reordenar la matriz según el clustering obtenido
    clustered_data = heatmap_data.loc[row_order, col_order]

    # Definir la escala de colores personalizada:
    # - Valores bajos (subexpresados): azul
    # - 0 (expresión similar al método de normalización): amarillo pastel
    # - Valores altos (sobreexpresados): rojo
    custom_colorscale = [
        [0.0, '#2c7bb6'],   # Azul
        [0.5, '#ffffb2'],   # Amarillo pastel
        [1.0, '#d7191c']    # Rojo
    ]

    # Crear subplots para combinar dendrogramas y heatmap
    fig = make_subplots(
        rows=2, cols=2,
        shared_xaxes=True,
        shared_yaxes=True,
        vertical_spacing=0.02,
        horizontal_spacing=0.02,
        column_widths=[0.8, 0.2],
        row_heights=[0.10, 1],
        specs=[[{"type": "scatter", "colspan": 2}, None],
               [{"type": "heatmap"}, {"type": "scatter"}]]
    )

    # Añadir dendrograma de columnas
    for trace in dendro_col['data']:
        fig.add_trace(trace, row=1, col=1)

    # Añadir dendrograma de filas
    for trace in dendro_row['data']:
        fig.add_trace(trace, row=2, col=2)

    # Añadir heatmap utilizando la escala de colores personalizada y centrándose en 0
    heatmap = go.Heatmap(
        z=clustered_data.values,
        x=col_order,
        y=row_order,
        colorscale=custom_colorscale,
        colorbar=dict(
            title=f"log2(Exp. Rel. { 'vs ' + gen_referencia if metodo_elegido=='gen de referencia' else 'por promedios' })"
        ),
        showscale=True,
        zmid=0
    )
    fig.add_trace(heatmap, row=2, col=1)

    # Ajustes de layout
    fig.update_layout(
        title_text=title,
        width=1000,
        height=800,
        showlegend=False
    )

    # Configuración de ejes
    fig.update_xaxes(showticklabels=False, row=1, col=1)
    fig.update_xaxes(showticklabels=True, row=2, col=1, title="Tests")
    fig.update_yaxes(showticklabels=True, row=2, col=1, title="Genes")
    fig.update_yaxes(showticklabels=False, row=2, col=2)

    return fig

# Definir el título según el método de normalización
if metodo_elegido == "gen de referencia":
    norm_title = f"Normalizado vs {gen_referencia}"
elif metodo_elegido == "promedios":
    norm_title = "Normalizado por promedios"
else:
    norm_title = "Normalización desconocida"

# ------------------------------
# 6. Clustergrama interactivo para genes sobreexpresados
# ------------------------------
fig_sobre = create_clustered_heatmap(
    heatmap_data_sobre,
    title=f"Expresión de Genes Sobreexpresados en {tipo_cancer} con Dendrograma ({norm_title})"
)
guardar_grafico(fig_sobre, f"clustergrama_sobreexpresados_{tipo_cancer}.png")

fig_sobre.show()

# ------------------------------
# 7. Clustergrama interactivo para genes subexpresados
# ------------------------------
fig_sub = create_clustered_heatmap(
    heatmap_data_sub,
    title=f"Expresión de Genes Subexpresados en {tipo_cancer} con Dendrograma ({norm_title})"
)
guardar_grafico(fig_sub, f"clustergrama_subexpresados_{tipo_cancer}.png")
fig_sub.show()
