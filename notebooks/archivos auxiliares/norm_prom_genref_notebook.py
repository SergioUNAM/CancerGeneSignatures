# ============================================================
# CANCER GENE SIGNATURES — PIPELINE ΔΔCT
# Archivo: norm_prom_genref_notebook.py
# Propósito: Cargar raw Excel, clasificar CT en controles/muestras,
#            imputar CT, filtrar controles de máquina, calcular ΔΔCT
#            por dos métodos (promedio global y gen de referencia),
#            exportar resultados y visualizar.
# ============================================================

# ============================== SECCIÓN 0. IMPORTS ===============================
import pandas as pd
import numpy as np
from typing import Optional, Dict, Tuple
import plotly.graph_objects as go
from IPython.display import display, Markdown

# ==================== SECCIÓN 1. CARGA Y VALIDACIÓN DE ARCHIVO ====================
# Entradas: ruta (str)
# Salidas: DataFrame crudo
# - Lee Excel crudo (header=None)
# - Maneja errores comunes y notifica
def cargar_datos_excel(ruta):
    display(Markdown("### Resultados de la extracción"))
    try:
        df = pd.read_excel(ruta, engine='openpyxl', header=None)
        print("Archivo cargado correctamente")
        return df
    except FileNotFoundError:
        print(f"Archivo '{ruta}' no encontrado")
    except Exception as e:
        print(f"Error al cargar: {e}")
    return None

# ==================== SECCIÓN 2. METADATOS: TESTS, POZOS Y GENES ==================
# Entradas: df (DataFrame crudo)
# Salidas: test_names (list[str]), pozos (list[str]), genes (list[str])
# - Extrae nombres de pruebas (fila 5), pozos (col 0) y genes (col 1)
# - Valida forma mínima de la tabla
def extraer_informacion_cg(df: pd.DataFrame) -> tuple:
    rows, cols = df.shape
    if rows < 9 or cols < 3:
        print("El archivo no tiene el formato esperado.")
        return [], [], []
    test_names = df.iloc[5, 2:].dropna().astype(str).tolist()
    pozos = df.iloc[8:, 0].dropna().astype(str).tolist()
    genes = df.iloc[8:, 1].dropna().astype(str).tolist()
    return test_names, pozos, genes

# ====================== SECCIÓN 3. CLASIFICACIÓN DE CT POR TIPO ===================
# Entradas: df (crudo), genes (list[str]), nombre_controles (str), nombre_muestras (str)
# Salidas: controles_df, muestras_df
# - Detecta columnas 'CT' y separa en controles vs muestras
# - Empata longitud con lista de genes
def procesar_ct(df: pd.DataFrame, genes: list, nombre_controles: str, nombre_muestras: str) -> tuple:
    if len(df) <= 8:
        print("El DataFrame no tiene suficientes filas.")
        return pd.DataFrame(), pd.DataFrame()

    fila_ct = df.iloc[7]
    if not any(str(val).strip().upper() == "CT" for val in fila_ct.values):
        print("No se encontró ninguna columna con 'CT' en la fila 7.")
        return pd.DataFrame(), pd.DataFrame()

    # Normalizar prefijos
    nombre_controles = nombre_controles.strip().upper()
    nombre_muestras = nombre_muestras.strip().upper()

    controles_list = []
    muestras_list = []

    for idx in np.where(fila_ct == "CT")[0]:
        raw_test_name = str(df.iloc[5, idx])
        nombre_test = raw_test_name.strip().upper()
        print(f"[DEBUG] idx={idx}, nombre_test='{nombre_test}'")

        ct_values = df.iloc[8:, idx].tolist()
        if len(ct_values) < len(genes):
            ct_values.extend([np.nan] * (len(genes) - len(ct_values)))
        elif len(ct_values) > len(genes):
            ct_values = ct_values[:len(genes)]

        df_temp = pd.DataFrame({
            'test': [nombre_test] * len(genes),
            'target': genes,
            'ct': ct_values
        })

        if nombre_test.startswith(nombre_controles):
            controles_list.append(df_temp)
            print(f"Control detectado: {nombre_test}")
        elif nombre_test.startswith(nombre_muestras):
            muestras_list.append(df_temp)
            print(f"Muestra detectada: {nombre_test}")
        else:
            print(f"'{nombre_test}' no coincide con control ni muestra. Se omite.")

    controles_df = pd.concat(controles_list, ignore_index=True) if controles_list else pd.DataFrame()
    muestras_df = pd.concat(muestras_list, ignore_index=True) if muestras_list else pd.DataFrame()
    return controles_df, muestras_df

# ===================== SECCIÓN 4. ORQUESTACIÓN DEL PIPELINE BÁSICO =================
# Entradas: ruta_archivo (str), nombre_controles (str), nombre_muestras (str)
# Salidas: dict con test_names, pozos, genes, controles_df, muestras_df
# - Carga -> extracción de info -> partición CT -> resúmenes
def flujo_procesamiento(ruta_archivo: str, nombre_controles: str, nombre_muestras: str) -> dict:
    df = cargar_datos_excel(ruta_archivo)
    if df is None:
        raise SystemExit(f"❌ Error crítico: No se pudo cargar el archivo '{ruta_archivo}'")

    test_names, pozos, genes = extraer_informacion_cg(df)

    for titulo, datos, filtro in [
        ("### Nombres de las pruebas realizadas", test_names, lambda x: x[1:]),
        ("### Genes objetivo analizados", genes, lambda x: [g for g in x if pd.notna(g) and g.strip()])
    ]:
        display(Markdown(titulo))
        datos_filtrados = filtro(datos)
        print(f"Total {len(datos_filtrados)}: \n{', '.join(datos_filtrados)}")

    controles_df, muestras_df = procesar_ct(df, genes, nombre_controles, nombre_muestras)

    display(Markdown("### Resumen de clasificación"))
    for tipo, df_tipo in [('Controles', controles_df), ('Muestras', muestras_df)]:
        if not df_tipo.empty:
            unique_tests = df_tipo['test'].unique()
            print(f"{tipo}: {len(unique_tests)} pruebas → {', '.join(unique_tests)}")
        else:
            print(f"No se detectaron {tipo.lower()} en el archivo.")

    print("Proceso de CT completado exitosamente")
    return {
        "test_names": test_names,
        "pozos": pozos,
        "genes": genes,
        "controles_df": controles_df,
        "muestras_df": muestras_df
    }

# ======================= SECCIÓN 5. PARÁMETROS Y EJECUCIÓN RÁPIDA ==================
# Entradas: ruta_archivo, nombre_controles, nombre_muestras
# Salidas: controles_df, muestras_df
# - Configuración de archivo y prefijos
# - Ejecución de flujo y obtención de dataframes base
ruta_archivo = "../raw_data/mir-29_PRUEBA_PLANTILLA.xlsx"
nombre_controles = "7BSC"
nombre_muestras = "3CM"
resultados = flujo_procesamiento(ruta_archivo, nombre_controles, nombre_muestras)
controles_df = resultados['controles_df']
muestras_df = resultados['muestras_df']

# =================== SECCIÓN 6. LIMPIEZA E IMPUTACIÓN DE VALORES CT =================
# Entradas: controles_df, muestras_df
# Salidas: controles_df, muestras_df imputados + flags
# - Normaliza tokens de 'Undetermined' y NaN
# - Imputa por límite técnico (recomendado) o máximo observado
def procesar_ct_column(
    df1: pd.DataFrame,
    df2: pd.DataFrame,
    columna_ct: str = 'ct',
    *,
    max_ct: Optional[float] = None,          # p.ej. 40 o 45
    by: Optional[list] = None,               # p.ej. ['plate_id','target']
    undetermined_tokens: tuple = ('undetermined','undet','nd','na'),
    inplace: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Imputa 'Undetermined' y faltantes en Ct con un valor límite (preferible al máximo observado).
    - Si max_ct se da: usa ese límite (recomendado).
    - Si no: usa el máximo observado *estratificado por 'by'*. Evita mezclar controles/muestras.
    Añade una columna booleana '<columna_ct>_imputed'.
    """
    def _clean(df: pd.DataFrame) -> pd.DataFrame:
        d = df if inplace else df.copy()
        # normaliza strings y convierte a numérico
        d[columna_ct] = (
            d[columna_ct]
            .astype(str)
            .str.strip()
            .str.lower()
            .replace({t: pd.NA for t in undetermined_tokens})
        )
        d[columna_ct] = pd.to_numeric(d[columna_ct], errors='coerce')
        return d

    a = _clean(df1)
    b = _clean(df2)

    # construir máscara de NaN antes de imputar
    a_flag = a[columna_ct].isna()
    b_flag = b[columna_ct].isna()

    if max_ct is not None:
        # imputación fija por límite técnico
        a.loc[a_flag, columna_ct] = float(max_ct)
        b.loc[b_flag, columna_ct] = float(max_ct)
        used_info = f"Imputación por límite técnico max_ct={max_ct}"
    else:
        # imputación por máximo observado, pero estratificada por 'by' si existe
        if by:
            # concatenar para calcular máximos por estrato sin mezclar controles/muestras
            comb = pd.concat([
                a[by + [columna_ct]],
                b[by + [columna_ct]]
            ], axis=0)
            grp_max = comb.groupby(by, dropna=False)[columna_ct].max(min_count=1)

            def fill_grouped(df, flag):
                if not by:
                    return df
                # join del máximo por estrato
                m = df.merge(grp_max.rename('ct_max_grp').reset_index(), on=by, how='left')
                df.loc[flag, columna_ct] = m.loc[flag, 'ct_max_grp'].values
                return df

            a = fill_grouped(a, a_flag)
            b = fill_grouped(b, b_flag)
            used_info = f"Imputación por máximo observado estratificado por {by}"
        else:
            # global, pero evitando todo-NaN
            valor_maximo = pd.concat([a[columna_ct], b[columna_ct]], ignore_index=True).max()
            if pd.isna(valor_maximo):
                raise ValueError(
                    "No hay valores numéricos válidos de Ct; provee max_ct o revisa los datos."
                )
            a.loc[a_flag, columna_ct] = valor_maximo
            b.loc[b_flag, columna_ct] = valor_maximo
            used_info = f"Imputación por máximo observado global={valor_maximo}"

    # flags
    a[f'{columna_ct}_imputed'] = a_flag
    b[f'{columna_ct}_imputed'] = b_flag

    # validaciones suaves
    for d, name in [(a,'df1'), (b,'df2')]:
        if not pd.api.types.is_float_dtype(d[columna_ct]):
            d[columna_ct] = d[columna_ct].astype(float)
        # rangos plausibles 0-50
        out_of_range = d.loc[(d[columna_ct] < 0) | (d[columna_ct] > 50), columna_ct]
        if len(out_of_range):
            raise ValueError(f"{name}: Ct fuera de rango plausible: {out_of_range.unique()}")

    # podrías registrar used_info en logs en vez de display
    a.attrs['imputation_info'] = used_info
    b.attrs['imputation_info'] = used_info

    return a, b

# Aplicar imputación básica (ajusta max_ct según tu plataforma)
controles_df, muestras_df = procesar_ct_column(controles_df, muestras_df, max_ct=40)
print(controles_df.attrs['imputation_info'])
print("Controles imputados:", controles_df['ct_imputed'].sum())
print("Muestras imputadas:", muestras_df['ct_imputed'].sum())

# ============= SECCIÓN 7. FILTRADO DE CONTROLES DE MÁQUINA (PPC/RTC) ==============
# Entradas: controles_df, muestras_df
# Salidas: controles_df, muestras_df sin PPC/RTC
# - Elimina targets que comienzan con prefijos de control de plataforma
def filtrar_controles_maquina(df, columna_target='target', controles_maquina=None):
    """
    Filtra controles de máquina cuyo target inicia con los prefijos dados (por defecto PPC/RTC).
    """
    if controles_maquina is None:
        controles_maquina = ['PPC', 'RTC']

    print(f"Controles de máquina a eliminar: {', '.join(controles_maquina)}")
    df_filtrado = df[~df[columna_target].str.startswith(tuple(controles_maquina))]
    print("Filas eliminadas con controles de máquina")
    return df_filtrado

display(Markdown("### Filtrado de controles de máquina para los datos en controles"))
controles_df = filtrar_controles_maquina(controles_df)
display(Markdown("### Filtrado de controles de máquina para los datos en muestras"))
muestras_df = filtrar_controles_maquina(muestras_df)

# Añadir columna para identificar el tipo de cada fila
controles_df = controles_df.assign(tipo='Control')
muestras_df = muestras_df.assign(tipo='Muestra')

# ================= SECCIÓN 8. EXPORTACIÓN DE TABLAS LIMPIAS PRE-NORMALIZACIÓN =================
# [EXPORT] Entradas: controles_df, muestras_df → generar archivos "datos_control_limpio" y "datos_muestras_limpio" en ruta_resultados

# ====================== SECCIÓN 9. DOCUMENTACIÓN DE METODOLOGÍAS ====================
# - Derivaciones matemáticas usadas en el cálculo
"""
## Procesamiento de datos

### Análisis de expresión génica mediante CT y Fold Change utilizando las metodologías promedio CT de un gen de referencia y promedio global

#### Método con gen de referencia

1. Promedio de Ct para cada gen:

$$\text{Promedio } Ct_{\text{gen de referencia}} \hspace{5pt}= \frac{Ct_1 + Ct_2 + Ct_3 + ... +Ct_n}{n}$$

2. ΔCt para cada condición:

$$\Delta Ct = \text{Promedio } Ct_{\text{gen de interés}} \hspace{5pt}- \text{Promedio } Ct_{\text{gen de referencia}}$$

3. ΔΔCt:

$$\Delta\Delta Ct = \Delta Ct_{\text{muestra}} - \Delta Ct_{\text{control}}
$$

4. Fold Change:

$$\text{Fold Change} = 2^{-\Delta\Delta Ct}$$

#### Métdodo con promedio global

1. Promedio de Ct para cada gen:

$$
\text{Promedio } Ct_{\text{gen}} = \frac{Ct_1 + Ct_2 + Ct_3 + ... + Ct_n}{n}
$$

2. Promedio global de Ct para cada condición

$$
\text{Promedio global } Ct_{\text{condición}} \hspace{5pt}= \frac{\sum \text{Promedio } Ct_{\text{gen}}}{N}
$$

3. ΔCt para cada gen:

$$
\Delta Ct_{\text{gen}} = \text{Promedio } Ct_{\text{gen}} - \text{Promedio global } Ct_{\text{condición}}
$$

4. ΔΔCt para cada gen:

$$
\Delta\Delta Ct_{\text{gen}} = \Delta Ct_{\text{gen, tratamiento}} \hspace{5pt}- \Delta Ct_{\text{gen, control}}
$$

5. Fold Change:

$$
\text{Fold Change} = 2^{-\Delta\Delta Ct_{\text{gen}}}
$$
"""

# ================== SECCIÓN 10. CÁLCULO DE ΔΔCT Y FOLD CHANGE ======================
# Entradas: controles_df, muestras_df
# Salidas: df_consolidado, df_promedios, df_gen_ref, gen_referencia
# - Método promedio global y método con gen de referencia estable
# - Selección de gen más estable por menor desviación promedio
def calcular_fold_change(
    controles: pd.DataFrame,
    muestras: pd.DataFrame
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, str]:
    """
    Calcula el fold change usando los dataframes de controles y muestras.

    1. Agrupamos los datos por 'target' y calculamos la media y la desviación estándar (CT).
    2. Calculamos los delta CT comparando cada gen con el promedio general.
    3. Seleccionamos el gen más estable: el que tiene la menor variabilidad (promedio de desviaciones).
    4. Utilizamos este gen de referencia para normalizar y obtener el fold change.
    """

    # Función auxiliar para procesar cada grupo (controles o muestras)
    def procesar_grupo(df: pd.DataFrame, grupo: str) -> pd.DataFrame:
        stats = df.groupby('target', as_index=False).agg(
            ct_promedio=('ct', 'mean'),
            ct_std=('ct', 'std')
        ).rename(columns={
            'ct_promedio': f'ct_promedio_{grupo}',
            'ct_std': f'ct_std_{grupo}'
        })
        return stats

    # Procesar controles y muestras
    stats_controles = procesar_grupo(controles, 'controles')
    stats_muestras = procesar_grupo(muestras, 'muestras')

    # Promedios generales
    stats_controles['promedio_general_controles'] = stats_controles['ct_promedio_controles'].mean()
    stats_muestras['promedio_general_muestras'] = stats_muestras['ct_promedio_muestras'].mean()

    # ΔCt respecto al promedio del grupo
    stats_controles['delta_ct_promedio_controles'] = (
        stats_controles['ct_promedio_controles'] - stats_controles['promedio_general_controles']
    )
    stats_muestras['delta_ct_promedio_muestras'] = (
        stats_muestras['ct_promedio_muestras'] - stats_muestras['promedio_general_muestras']
    )

    # Merge inicial
    df_consolidado = pd.merge(
        stats_controles, stats_muestras,
        on='target', how='outer', suffixes=('_controles', '_muestras')
    )

    # Gen más estable
    df_consolidado['stability'] = (df_consolidado['ct_std_controles'] + df_consolidado['ct_std_muestras']) / 2
    gen_mas_estable = df_consolidado.loc[df_consolidado['stability'].idxmin(), 'target']
    print(f"Gen de referencia seleccionado: {gen_mas_estable}")

    # ΔCt respecto al gen de referencia
    ref_controles = stats_controles.query("target == @gen_mas_estable")['ct_promedio_controles'].values[0]
    ref_muestras = stats_muestras.query("target == @gen_mas_estable")['ct_promedio_muestras'].values[0]
    stats_controles['delta_ct_gen_ref_controles'] = stats_controles['ct_promedio_controles'] - ref_controles
    stats_muestras['delta_ct_gen_ref_muestras'] = stats_muestras['ct_promedio_muestras'] - ref_muestras

    # Merge final para FC por ambos métodos
    df_consolidado = pd.merge(
        stats_controles, stats_muestras,
        on='target', how='outer', suffixes=('_controles', '_muestras')
    )

    # ΔΔCt y FC — método promedio
    df_consolidado['delta_delta_ct_promedio'] = (
        df_consolidado['delta_ct_promedio_muestras'] - df_consolidado['delta_ct_promedio_controles']
    )
    df_consolidado['fold_change_promedio'] = 2 ** (-df_consolidado['delta_delta_ct_promedio'])

    # ΔΔCt y FC — método gen de referencia
    df_consolidado['delta_delta_ct_gen_ref'] = (
        df_consolidado['delta_ct_gen_ref_muestras'] - df_consolidado['delta_ct_gen_ref_controles']
    )
    df_consolidado['fold_change_gen_ref'] = 2 ** (-df_consolidado['delta_delta_ct_gen_ref'])

    # DataFrames especializados
    cols_promedio = ['target'] + [c for c in df_consolidado.columns if 'promedio' in c]
    cols_gen_ref = ['target'] + [c for c in df_consolidado.columns if 'gen_ref' in c]
    df_promedios = df_consolidado[cols_promedio].copy().rename(columns={'fold_change_promedio': 'fold_change'})
    df_gen_ref = df_consolidado[cols_gen_ref].copy().rename(columns={'fold_change_gen_ref': 'fold_change'})

    print("Cálculo de fold change finalizado 🚀")
    return df_consolidado, df_promedios, df_gen_ref, gen_mas_estable

# Ejemplo de uso: cálculo principal
df_consolidado, df_promedios, df_gen_ref, gen_referencia = calcular_fold_change(controles_df, muestras_df)

# ================= SECCIÓN 11. EXPORTACIÓN DE RESULTADOS DE NORMALIZACIÓN ================
# [EXPORT] Entradas: df_consolidado, df_promedios, df_gen_ref → exportar tablas de resultados de normalización

# ===================== SECCIÓN 12. TABLA INTERACTIVA COMPARATIVA ====================
# Entradas: df_consolidado
# Salidas: fig (plotly)
# - Tabla Plotly: fold_change por método con codificación por color
SETTINGS = {
    "color_thresholds": [1, 2],
    "colors": ['#9fccff', '#FFFBCA', '#FFCCE1'],
    "header_style": {
        "fill_color": '#2C3E50',
        "font_color": 'white',
        "align": ['center', 'center', 'center'],
        "font_size": 14,
        "height": 40
    },
    "cell_style": {
        "font_size": 12,
        "align": ['left', 'center', 'center'],
        "height": 30
    },
    "decimals": 2,
    "char_width": 7,
    "min_column_width": 60,
    "max_column_width": 200,
    "row_height": 35,
    "header_height": 40,
    "base_height": 300
}

def prepare_data(df):
    return (df[['target', 'fold_change_promedio', 'fold_change_gen_ref']]
            .sort_values(by='target', ascending=True)
            .assign(
                fold_change_promedio=lambda x: x['fold_change_promedio'].round(SETTINGS['decimals']),
                fold_change_gen_ref=lambda x: x['fold_change_gen_ref'].round(SETTINGS['decimals'])
            ))

def calculate_colors(series):
    return np.select(
        [series > SETTINGS["color_thresholds"][1],
         series > SETTINGS["color_thresholds"][0]],
        SETTINGS["colors"][::-1][:2],
        SETTINGS["colors"][0]
    )

def calculate_column_widths(df):
    col_widths = []
    for col in df.columns:
        max_header_length = len(str(col))
        max_content_length = df[col].astype(str).apply(len).max()
        max_length = max(max_header_length, max_content_length)
        calculated_width = max_length * SETTINGS["char_width"]
        final_width = np.clip(calculated_width, SETTINGS["min_column_width"], SETTINGS["max_column_width"])
        col_widths.append(final_width)
    return col_widths

def create_interactive_table(df):
    df_display = df.copy()
    for col in ['fold_change_promedio', 'fold_change_gen_ref']:
        df_display[col] = df_display[col].astype(str)

    metadata = {
        "columns": df_display.columns.tolist(),
        "promedio_colors": calculate_colors(df['fold_change_promedio']),
        "gen_ref_colors": calculate_colors(df['fold_change_gen_ref']),
        "values": [df_display[col] for col in df_display.columns]
    }

    fig = go.Figure(data=[go.Table(
        columnwidth=calculate_column_widths(df_display),
        header=dict(
            values=[f"<b>{col}</b>" for col in metadata["columns"]],
            **SETTINGS["header_style"],
            font=dict(
                size=SETTINGS["header_style"]["font_size"],
                color=SETTINGS["header_style"]["font_color"]
            )
        ),
        cells=dict(
            values=metadata["values"],
            fill_color=[
                ['#EEEEEE']*len(df),
                metadata["promedio_colors"],
                metadata["gen_ref_colors"]
            ],
            **SETTINGS["cell_style"],
            font=dict(size=SETTINGS["cell_style"]["font_size"])
        )
    )])

    fig.update_layout(
        title={
            'text': f'<b>Tabla comparativa promedios vs gen de referencia ({gen_referencia})</b>',
            'y': 0.97,
            'x': 0.5,
            'xanchor': 'center',
            'font': dict(size=20, color='#2C3E50')
        },
        margin=dict(t=100, b=20, l=20, r=20),
        paper_bgcolor='white',
        autosize=False,
        height=300 + len(df)*4
    )

    return fig

df_fold_change = df_consolidado[['target', 'fold_change_promedio', 'fold_change_gen_ref']]
df_processed = prepare_data(df_fold_change)
fig = create_interactive_table(df_processed)
fig.show(config={'responsive': False})
# [EXPORT] Figura: fig → guardar como 'analisis_comparativo.png'

# ===================== SECCIÓN 13. GRÁFICO COMPARATIVO (ΔΔCT vs FC) ====================
# Entradas: df_consolidado
# Salidas: fig (plotly)
# - Barras: ΔΔCT por método (eje izquierdo)
# - Líneas: FC por método (eje derecho, log)
fig = go.Figure()

# 1. Barras para Delta-Delta CT (Eje Y Izquierdo)
fig.add_trace(go.Bar(
    x=df_consolidado['target'],
    y=df_consolidado['delta_delta_ct_promedio'],
    name='ΔΔCT (Promedios)',
    marker_color='#1f77b4',
    opacity=0.8,
    yaxis='y'
))
fig.add_trace(go.Bar(
    x=df_consolidado['target'],
    y=df_consolidado['delta_delta_ct_gen_ref'],
    name='ΔΔCT (Gen Ref)',
    marker_color='#ff7f0e',
    opacity=0.8,
    yaxis='y'
))

# 2. Líneas + Marcadores para Fold Change (Eje Y Derecho - Escala Log)
fig.add_trace(go.Scatter(
    x=df_consolidado['target'],
    y=df_consolidado['fold_change_promedio'],
    name='Fold Change (Promedios)',
    mode='markers+lines',
    marker=dict(color='#2ca02c', size=10, symbol='diamond'),
    line=dict(color='#2ca02c', width=2, dash='dot'),
    yaxis='y2'
))
fig.add_trace(go.Scatter(
    x=df_consolidado['target'],
    y=df_consolidado['fold_change_gen_ref'],
    name='Fold Change (Gen Ref)',
    mode='markers+lines',
    marker=dict(color='#d62728', size=10, symbol='diamond'),
    line=dict(color='#d62728', width=2, dash='dot'),
    yaxis='y2'
))

# 3. Ajustes de Layout
fig.update_layout(
    title=dict(
        text='Análisis Comparativo de Métodos de Cálculo',
        font=dict(size=24, family='Arial Black'),
        x=0.5
    ),
    template='plotly_white',
    barmode='group',
    yaxis=dict(
        title=dict(text='ΔΔCT Value', font=dict(color='#1f77b4', size=14)),
        tickfont=dict(color='#1f77b4'),
        showgrid=True,
        gridcolor='lightgray'
    ),
    yaxis2=dict(
        title=dict(text='Fold Change (Log Scale)', font=dict(color='#d62728', size=14)),
        tickfont=dict(color='#d62728'),
        overlaying='y',
        side='right',
        type='log',
        showgrid=False
    ),
    hovermode='x unified',
    legend=dict(
        orientation="h",
        yanchor="bottom",
        y=1.02,
        xanchor="right",
        x=1,
        bgcolor='rgba(255,255,255,0.9)'
    ),
    height=650,
    width=1200,
    margin=dict(b=100, t=100),
    plot_bgcolor='rgba(245,245,245,0.9)'
)

# 4. Anotaciones Clave
fig.add_annotation(
    text="Fold Change = 1",
    xref="paper", yref="y2",
    x=0.05, y=1,
    showarrow=False,
    font=dict(color="#d62728", size=12)
)

fig.show()
# [EXPORT] Figura: fig → guardar como 'analisis_comparativo.png'

# =================== SECCIÓN 16. CATEGORIZACIÓN DE NIVEL DE EXPRESIÓN ====================
# Entradas: df_promedios o df_gen_ref (elige uno para categorizar)
# Salidas: df_combined_expresion y listas por categoría
# - Etiquetas: subexpresado (<1), estable [1,2), sobreexpresado (>=2)
def categorizar_expresion(df, umbral_bajo=1, umbral_alto=2):
    # Validación de columnas
    required_cols = ['target', 'fold_change']
    if not all(col in df.columns for col in required_cols):
        raise ValueError("Columnas requeridas no encontradas")

    # Categorización vectorizada
    condiciones = [
        (df['fold_change'] < umbral_bajo),
        (df['fold_change'] >= umbral_bajo) & (df['fold_change'] < umbral_alto),
        (df['fold_change'] >= umbral_alto)
    ]

    categorias = ['subexpresado', 'estable', 'sobreexpresado']

    df_categorizado = df[required_cols].copy()
    df_categorizado['nivel_expresion'] = np.select(condiciones, categorias, default='sin_categorizar')
    df_categorizado.sort_values('fold_change', ascending=False)

    return df_categorizado

# Elegir método para categorizar: aquí por gen de referencia
df_combined_expresion = categorizar_expresion(df_gen_ref)

# Extraer listas por categoría
df_sobreexpresados = df_combined_expresion[df_combined_expresion['nivel_expresion'] == 'sobreexpresado'][['target']]
lista_sobreexpresados = df_sobreexpresados['target'].tolist()
df_estables = df_combined_expresion[df_combined_expresion['nivel_expresion'] == 'estable'][['target']]
lista_estables = df_estables['target'].tolist()
df_subexpresados = df_combined_expresion[df_combined_expresion['nivel_expresion'] == 'subexpresado'][['target']]
lista_subexpresados = df_subexpresados['target'].tolist()

# ===================== SECCIÓN 17. EXPORTACIÓN DE CATEGORÍAS POR GEN =====================
# [EXPORT] Entradas: df_combined_expresion, df_sobreexpresados, df_estables, df_subexpresados → exportar categorías por gen

# =============================== FIN DEL PIPELINE =================================
# Notas:
# - Este script evita funciones auxiliares decorativas; usa print() para trazas mínimas.
# - Marcadores [EXPORT] indican qué objetos deben exportarse externamente.