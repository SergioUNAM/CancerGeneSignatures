import pandas as pd
import numpy as np


def cargar_datos_excel(ruta: str) -> pd.DataFrame:
    """Carga un archivo de Excel devolviendo un DataFrame."""
    return pd.read_excel(ruta, engine="openpyxl")


def extraer_informacion_cg(df: pd.DataFrame) -> tuple[list, list, list]:
    """Obtiene nombres de pruebas, pozos y genes de un DataFrame."""
    test_names = df.iloc[0].dropna().astype(str).tolist()
    pozos = df.iloc[2:, 0].dropna().astype(str).tolist()
    genes = df.iloc[3:, 1].dropna().astype(str).tolist()
    return test_names, pozos, genes


def procesar_ct(df: pd.DataFrame, genes: list, nombre_controles: str, nombre_muestras: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Separa los datos en controles y muestras a partir de las columnas CT."""
    fila_ct = df.iloc[2]
    indices_ct = np.where(fila_ct == "CT")[0]
    controles = []
    muestras = []
    for idx in indices_ct:
        nombre_test = str(df.iloc[0, idx])
        bloque = pd.DataFrame({
            "test": [nombre_test] * len(genes),
            "target": genes,
            "ct": df.iloc[3:, idx].tolist(),
        })
        if nombre_test.startswith(nombre_controles):
            controles.append(bloque)
        elif nombre_test.startswith(nombre_muestras):
            muestras.append(bloque)
    return pd.concat(controles, ignore_index=True), pd.concat(muestras, ignore_index=True)


def filtrar_controles_maquina(df: pd.DataFrame, columna_target: str = "target", controles_maquina: list | None = None) -> pd.DataFrame:
    """Elimina filas cuyo target coincide con controles de máquina."""
    if controles_maquina is None:
        controles_maquina = ["PPC", "RTC", "miRTC"]
    return df[~df[columna_target].str.startswith(tuple(controles_maquina))]


def calcular_fold_change(controles: pd.DataFrame, muestras: pd.DataFrame) -> pd.DataFrame:
    """Calcula Fold Change usando el método del gen de referencia más estable."""

    def resumen(df: pd.DataFrame, grupo: str) -> pd.DataFrame:
        stats = df.groupby("target", as_index=False).agg(ct_promedio=("ct", "mean"), ct_std=("ct", "std"))
        stats["promedio_%s" % grupo] = stats["ct_promedio"].mean()
        stats["delta_ct_%s" % grupo] = stats["ct_promedio"] - stats["promedio_%s" % grupo]
        return stats

    stats_c = resumen(controles, "controles")
    stats_m = resumen(muestras, "muestras")

    df_consolidado = pd.merge(stats_c, stats_m, on="target", suffixes=("_controles", "_muestras"))
    df_consolidado["stability"] = (df_consolidado["ct_std_controles"] + df_consolidado["ct_std_muestras"]) / 2
    gen_ref = df_consolidado.loc[df_consolidado["stability"].idxmin(), "target"]

    ref_c = stats_c.query("target == @gen_ref")["ct_promedio"].values[0]
    ref_m = stats_m.query("target == @gen_ref")["ct_promedio"].values[0]

    stats_c["delta_ct_ref_controles"] = stats_c["ct_promedio"] - ref_c
    stats_m["delta_ct_ref_muestras"] = stats_m["ct_promedio"] - ref_m

    df_final = pd.merge(stats_c, stats_m, on="target", suffixes=("_controles", "_muestras"))
    df_final["delta_delta_ct"] = df_final["delta_ct_ref_muestras"] - df_final["delta_ct_ref_controles"]
    df_final["fold_change"] = 2 ** (-df_final["delta_delta_ct"])
    return df_final
