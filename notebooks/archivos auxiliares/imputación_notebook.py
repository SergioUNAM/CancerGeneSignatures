import pandas as pd
from typing import Optional, Dict

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

controles_df, muestras_df = procesar_ct_column(controles_df, muestras_df)

# recuperar el texto guardado en atributos
print(controles_df.attrs['imputation_info'])
print("Controles imputados:", controles_df['ct_imputed'].sum())
print("Muestras imputadas:", muestras_df['ct_imputed'].sum())