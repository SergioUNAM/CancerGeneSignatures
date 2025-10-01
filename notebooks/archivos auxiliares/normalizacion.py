# =========================================================================================
# ETAPA DE NORMALIZADORES (SIN IMPUTACIÓN): ESTABILIDAD + NORMALIZACIÓN + FDR + ROBUSTEZ
# =========================================================================================
# REQUISITOS DE ENTRADA (ya resueltos por tus pasos previos):
#   - df_total con columnas obligatorias:
#       'test'    -> identificador de muestra (único por muestra biológica/técnica dentro del análisis)
#       'target'  -> gen (string)
#       'ct'      -> Ct numérico (ya imputado/limpio en tus pasos previos)
#       'tipo'    -> 'Control' o cualquier otra etiqueta considerada 'Muestra'
#
# OBJETIVO EN ESTA ETAPA:
#   1) Seleccionar referencias por ESTABILIDAD (no por “qué separa más”).
#      - Minimizar varianza intra-grupo y diferencia inter-grupo del Ct de la referencia.
#   2) Normalizar respecto a las referencias elegidas:
#      - ref_ct (por muestra) = promedio de Ct de las referencias
#      - delta_ct = ct - ref_ct
#      - log2_rel_expr = -delta_ct
#   3) Evaluar diferencia por gen (Welch t-test) con corrección BH-FDR y tamaño de efecto (Cohen's d).
#   4) Robustez:
#      - Bootstrap por muestra: frecuencia con la que un gen vuelve a salir significativo.
#      - Permutación de etiquetas: FPR empírica del pipeline.
#   5) Entregar:
#      - refs elegidas + score de estabilidad
#      - df_norm (con columnas agregadas)
#      - df_stats (genes con p, q, t_abs, d, medias y SD por grupo, y bootstrap_freq)
#      - df_heat (matriz para tu clustergrama)
#
# NOTA: No toco tus funciones de lectura ni tu create_clustered_heatmap(...).
# =========================================================================================

import numpy as np
import pandas as pd
from itertools import combinations
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

# -----------------------------------------------------------------------------------------
# 0) CAPA DEL “DF DE ENTRADA”
# -----------------------------------------------------------------------------------------
# Usa el df_total que ya existe en tu pipeline.
df_base = df_total.copy()

# Homogeneiza tipos básicos por si quedaron inconsistencias
df_base['test'] = df_base['test'].astype(str)
df_base['target'] = df_base['target'].astype(str)
df_base['tipo'] = df_base['tipo'].astype(str)

# Validación mínima
cols_req = {'test','target','ct','tipo'}
faltan = cols_req - set(df_base.columns)
if faltan:
    raise ValueError(f"Faltan columnas obligatorias en df_total: {faltan}")

# -----------------------------------------------------------------------------------------
# 1) FUNCIÓN DE ESTABILIDAD PARA UN CONJUNTO DE REFERENCIAS
# -----------------------------------------------------------------------------------------
def stability_score_for_refs(df, refs):
    """
    Calcula un score de estabilidad para un conjunto de genes de referencia 'refs'.
    Menor score = mejor (más estable).

    Lógica:
      - ref_ct por muestra (test) = promedio de Ct de las referencias.
      - s_intra: desviación estándar intra-grupo del ref_ct (promedio sobre tipos).
      - inter_diff: diferencia de medias de ref_ct entre grupos (Control vs Muestra).
      - score = 0.7*s_intra + 0.3*inter_diff (pesos ajustables).

    Recomendación:
      - Mantener refs pequeñas (K=2 o 3). K grande diluye biología y amplifica ruido técnico.
    """
    if not refs:
        return np.inf

    tmp = df.copy()
    ref_ct = (tmp[tmp['target'].isin(refs)]
              .groupby('test')['ct']
              .mean()
              .rename('ref_ct'))
    tmp = tmp.join(ref_ct, on='test')

    # Solo filas de las referencias para medir estabilidad de la referencia en sí
    ref_rows = tmp[tmp['target'].isin(refs)].copy()

    # SD intra-grupo de ref_ct (promedio de SDs por 'tipo')
    s_intra = (ref_rows.groupby(['tipo','test'])['ref_ct'].mean()
                        .groupby(level=0).std(ddof=1))
    s_intra_mean = s_intra.mean() if s_intra.notna().any() else np.inf

    # Diferencia inter-grupo de la media de ref_ct
    m_by_grp = ref_rows.groupby('tipo')['ref_ct'].mean()
    inter_diff = m_by_grp.max() - m_by_grp.min() if len(m_by_grp) >= 2 else np.inf

    score = 0.7 * s_intra_mean + 0.3 * inter_diff
    return float(score)

# -----------------------------------------------------------------------------------------
# 2) CANDIDATOS Y BÚSQUEDA DE REFERENCIAS (BRUTE FORCE HONESTA)
# -----------------------------------------------------------------------------------------
# Candidatos por menor SD de Ct (no uses t para esto).
N_CANDIDATOS = 20
K_REFS = 2     # Puedes probar 3 si tienes replicación suficiente

sd_por_gen = (df_base.groupby('target')['ct'].std()
              .replace([np.inf, -np.inf], np.nan)
              .dropna()
              .sort_values())

candidatos = sd_por_gen.index.tolist()[:N_CANDIDATOS]
if len(candidatos) < K_REFS:
    raise ValueError(f"No hay suficientes candidatos ({len(candidatos)}) para K={K_REFS}.")

mejor_refs, mejor_score = None, np.inf
for refs in combinations(candidatos, K_REFS):
    sc = stability_score_for_refs(df_base, refs)
    if sc < mejor_score:
        mejor_refs, mejor_score = refs, sc

# -----------------------------------------------------------------------------------------
# 3) NORMALIZACIÓN CON LAS REFERENCIAS ELEGIDAS
# -----------------------------------------------------------------------------------------
def normaliza_con_refs(df, refs):
    """
    Añade columnas:
      - ref_ct: promedio Ct de las referencias por muestra
      - delta_ct = ct - ref_ct
      - log2_rel_expr = -delta_ct
    """
    out = df.copy()
    ref_ct = (out[out['target'].isin(refs)]
              .groupby('test')['ct']
              .mean()
              .rename('ref_ct'))
    out = out.join(ref_ct, on='test')
    out['delta_ct'] = out['ct'] - out['ref_ct']
    out['log2_rel_expr'] = -out['delta_ct']
    return out

df_norm = normaliza_con_refs(df_base, mejor_refs)

# -----------------------------------------------------------------------------------------
# 4) EVALUACIÓN DIFERENCIAL POR GEN (Welch + BH-FDR + Cohen's d)
# -----------------------------------------------------------------------------------------
def evaluar_dif_por_gen(temp):
    """
    Para cada gen, compara Control vs Muestra usando log2_rel_expr:
      - Welch t-test (no asume varianzas iguales)
      - p-value + BH-FDR (q)
      - Cohen's d (tamaño de efecto aproximado con SD combinada)
      - medias y SD por grupo (para interpretar magnitud)
    Devuelve DataFrame ordenado por q ascendente.
    """
    rows = []
    for g, sub in temp.groupby('target', sort=False):
        c = sub.loc[sub['tipo'] == 'Control', 'log2_rel_expr']
        m = sub.loc[sub['tipo'] != 'Control', 'log2_rel_expr']
        if len(c) > 1 and len(m) > 1:
            t, p = ttest_ind(c, m, equal_var=False)
            sd_p = np.sqrt((c.var(ddof=1) + m.var(ddof=1)) / 2.0)
            d = (m.mean() - c.mean()) / sd_p if sd_p > 0 else np.nan
            rows.append((g,
                         float(np.abs(t)), float(p), float(d),
                         float(c.mean()), float(m.mean()),
                         float(c.std(ddof=1) if len(c)>1 else np.nan),
                         float(m.std(ddof=1) if len(m)>1 else np.nan)))
    res = pd.DataFrame(rows, columns=[
        'gene','t_abs','p','cohen_d','mean_ctrl','mean_case','sd_ctrl','sd_case'
    ]).sort_values('p')
    if not res.empty:
        res['q'] = multipletests(res['p'], method='fdr_bh')[1]
        res = res.sort_values(['q','t_abs'])
    return res

df_stats = evaluar_dif_por_gen(df_norm)

# Criterio de significancia
ALPHA_Q = 0.05
genes_signif = df_stats.loc[df_stats['q'] < ALPHA_Q, 'gene'].tolist()

# -----------------------------------------------------------------------------------------
# 5) ROBUSTEZ: BOOTSTRAP POR MUESTRA Y PERMUTACIÓN PARA FPR
# -----------------------------------------------------------------------------------------
rng = np.random.default_rng(123)

def bootstrap_significancia(temp, refs, B=300, alpha=0.05):
    """
    Bootstrap a nivel 'test' con reemplazo.
    En cada iteración:
      - toma subconjunto re-muestreado
      - re-normaliza con las MISMAS refs
      - recalcula diferencias por gen
    Devuelve proporción de iteraciones con q<alpha por gen.
    """
    tests = temp['test'].unique()
    genes = temp['target'].unique()
    hits = dict.fromkeys(genes, 0)

    for _ in range(B):
        boots = rng.choice(tests, size=len(tests), replace=True)
        sub = temp[temp['test'].isin(boots)].copy()
        sub = normaliza_con_refs(sub, refs)
        r = evaluar_dif_por_gen(sub)
        if not r.empty:
            sig = set(r.loc[r['q'] < alpha, 'gene'])
            for g in sig:
                hits[g] += 1

    return pd.Series(hits, name='bootstrap_freq').div(B).sort_values(ascending=False)

def permutation_fpr(temp, refs, R=200, alpha=0.05):
    """
    Permuta etiquetas 'tipo' por muestra para romper asociación.
    Calcula la tasa de falsos positivos promedio del pipeline (empírico).
    """
    tests = temp['test'].unique()
    fprs = []

    for _ in range(R):
        # etiquetas permutadas por test
        tipo_por_test = temp.groupby('test')['tipo'].first()
        tipo_perm = pd.Series(rng.permutation(tipo_por_test.values), index=tipo_por_test.index)
        sub = temp.copy()
        sub['tipo_perm'] = sub['test'].map(tipo_perm)
        # Normalización con MISMAS referencias
        sub = normaliza_con_refs(sub, refs)

        # Evalúa usando 'tipo_perm'
        pvals = []
        for g, gsub in sub.groupby('target', sort=False):
            c = gsub.loc[gsub['tipo_perm'] == 'Control', 'log2_rel_expr']
            m = gsub.loc[gsub['tipo_perm'] != 'Control', 'log2_rel_expr']
            if len(c) > 1 and len(m) > 1:
                _, p = ttest_ind(c, m, equal_var=False)
                pvals.append(p)
        if pvals:
            q_perm = multipletests(pvals, method='fdr_bh')[1]
            fprs.append((np.array(q_perm) < alpha).mean())

    return float(np.nanmean(fprs)) if fprs else np.nan

# Ejecuta robustez (opcional; descomenta si necesitas rapidez)
bootstrap_freq = bootstrap_significancia(df_base, mejor_refs, B=300, alpha=ALPHA_Q)
fpr_empirica = permutation_fpr(df_base, mejor_refs, R=200, alpha=ALPHA_Q)

# Integra bootstrap al ranking
if not df_stats.empty:
    df_stats = df_stats.merge(bootstrap_freq.rename('bootstrap_freq'),
                              left_on='gene', right_index=True, how='left')
    df_stats['bootstrap_freq'] = df_stats['bootstrap_freq'].fillna(0.0)

# -----------------------------------------------------------------------------------------
# 6) MATRIZ PARA CLUSTERGRAMA Y RESUMEN
# -----------------------------------------------------------------------------------------
genes_para_heatmap = genes_signif if genes_signif else df_stats.head(20)['gene'].tolist()
df_heat = df_norm[df_norm['target'].isin(genes_para_heatmap)].pivot_table(
    index='target', columns='test', values='log2_rel_expr', aggfunc='mean'
)

resumen_refs = {
    "refs_elegidas": list(mejor_refs),
    "score_estabilidad": round(mejor_score, 4),
    "n_genes_signif_q<%.2f" % ALPHA_Q: int((df_stats['q'] < ALPHA_Q).sum()),
    "fpr_empirica_aprox": round(fpr_empirica, 4)
}

print("===== REFERENCIAS ELEGIDAS (ESTABILIDAD) =====")
print(resumen_refs)
print("\n===== TOP 20 GENES (por q) =====")
display(df_stats.head(20))

# Visualización si ya tienes tu función
# fig = create_clustered_heatmap(df_heat, title=f"Clustergrama: q<{ALPHA_Q} | refs: {', '.join(mejor_refs)}")
# fig.show()