import argparse
from pathlib import Path
from src.processing import (
    cargar_datos_excel,
    extraer_informacion_cg,
    procesar_ct,
    filtrar_controles_maquina,
    calcular_fold_change,
)
from src.visualization import tabla_fold_change


def run_pipeline(ruta_excel: str, nombre_controles: str, nombre_muestras: str, salida: str | None = None):
    df = cargar_datos_excel(ruta_excel)
    tests, pozos, genes = extraer_informacion_cg(df)
    controles_df, muestras_df = procesar_ct(df, genes, nombre_controles, nombre_muestras)
    controles_df = filtrar_controles_maquina(controles_df)
    muestras_df = filtrar_controles_maquina(muestras_df)

    resultados = calcular_fold_change(controles_df, muestras_df)
    fig = tabla_fold_change(resultados)

    if salida:
        fig.write_html(Path(salida))
    else:
        fig.show()


def main():
    parser = argparse.ArgumentParser(description="Procesamiento de datos qPCR")
    parser.add_argument("archivo", help="Ruta al archivo Excel")
    parser.add_argument("--controles", default="Control", help="Prefijo para las columnas de controles")
    parser.add_argument("--muestras", default="Sample", help="Prefijo para las columnas de muestras")
    parser.add_argument("--salida", help="Ruta al archivo HTML de salida")
    args = parser.parse_args()
    run_pipeline(args.archivo, args.controles, args.muestras, args.salida)


if __name__ == "__main__":
    main()
