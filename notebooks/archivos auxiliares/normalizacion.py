# ================================================================
# 2. Cálculo y visualización usando el gen más estable como referencia
# ================================================================
def calcular_fold_change(controles, muestras):
    def proc(df, grupo):
        stats = df.groupby('target', as_index=False).agg(
            ct_promedio=('ct', 'mean'), ct_std=('ct', 'std')
        ).rename(columns={'ct_promedio': f'ct_promedio_{grupo}', 'ct_std': f'ct_std_{grupo}'})
        return stats
    stats_cont = proc(controles, 'controles')
    stats_mues = proc(muestras, 'muestras')
    stats_cont['promedio_general_controles'] = stats_cont['ct_promedio_controles'].mean()
    stats_mues['promedio_general_muestras'] = stats_mues['ct_promedio_muestras'].mean()
    stats_cont['delta_ct_promedio_controles'] = stats_cont['ct_promedio_controles'] - stats_cont['promedio_general_controles']
    stats_mues['delta_ct_promedio_muestras'] = stats_mues['ct_promedio_muestras'] - stats_mues['promedio_general_muestras']
    df_consolidado = pd.merge(stats_cont, stats_mues, on='target', how='outer')
    df_consolidado['stability'] = (df_consolidado['ct_std_controles'] + df_consolidado['ct_std_muestras']) / 2
    gen_mas_estable = df_consolidado.loc[df_consolidado['stability'].idxmin(), 'target']
    ref_cont = stats_cont.query("target == @gen_mas_estable")['ct_promedio_controles'].values[0]
    ref_mues = stats_mues.query("target == @gen_mas_estable")['ct_promedio_muestras'].values[0]
    stats_cont['delta_ct_gen_ref_controles'] = stats_cont['ct_promedio_controles'] - ref_cont
    stats_mues['delta_ct_gen_ref_muestras'] = stats_mues['ct_promedio_muestras'] - ref_mues
    df_consolidado = pd.merge(stats_cont, stats_mues, on='target', how='outer')
    df_consolidado['delta_delta_ct_promedio'] = (
        df_consolidado['delta_ct_promedio_muestras'] - df_consolidado['delta_ct_promedio_controles']
    )
    df_consolidado['fold_change_promedio'] = 2 ** (-df_consolidado['delta_delta_ct_promedio'])
    df_consolidado['delta_delta_ct_gen_ref'] = (
        df_consolidado['delta_ct_gen_ref_muestras'] - df_consolidado['delta_ct_gen_ref_controles']
    )
    df_consolidado['fold_change_gen_ref'] = 2 ** (-df_consolidado['delta_delta_ct_gen_ref'])
    return df_consolidado, gen_mas_estable

df_consolidado, gen_referencia = calcular_fold_change(controles_df, muestras_df)

df_total = pd.concat([controles_df, muestras_df])
ref_ct = df_total[df_total['target'] == gen_referencia][['test', 'ct']].rename(
    columns={'ct': f'{gen_referencia}_ct'}
)
df_norm = df_total.merge(ref_ct, on='test')
df_norm['delta_ct'] = df_norm['ct'] - df_norm[f'{gen_referencia}_ct']
df_norm['log2_rel_expr'] = -df_norm['delta_ct']

df_separadores = df_consolidado.copy()
df_separadores['abs_delta'] = df_separadores['delta_delta_ct_promedio'].abs()
df_separadores = df_separadores.sort_values('abs_delta', ascending=False)
umbral = 2
genes_separadores = df_separadores[df_separadores['abs_delta'] > umbral]['target'].tolist()

# ====== Función para crear el clustergrama ======
def create_clustered_heatmap(heatmap_data, title):
    dendro_col = ff.create_dendrogram(
        heatmap_data.T,
        orientation='bottom',
        labels=heatmap_data.columns.tolist(),
        linkagefun=lambda x: linkage(x, method='average', metric='euclidean')
    )
    col_order = dendro_col['layout']['xaxis']['ticktext']
    dendro_row = ff.create_dendrogram(
        heatmap_data,
        orientation='right',
        labels=heatmap_data.index.tolist(),
        linkagefun=lambda x: linkage(x, method='average', metric='euclidean')
    )
    row_order = dendro_row['layout']['yaxis']['ticktext']
    clustered_data = heatmap_data.loc[row_order, col_order]
    custom_colorscale = [
        [0.0, '#2c7bb6'], [0.5, '#ffffb2'], [1.0, '#d7191c']
    ]
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
    for trace in dendro_col['data']:
        fig.add_trace(trace, row=1, col=1)
    for trace in dendro_row['data']:
        fig.add_trace(trace, row=2, col=2)
    heatmap = go.Heatmap(
        z=clustered_data.values,
        x=col_order,
        y=row_order,
        colorscale=custom_colorscale,
        colorbar=dict(title=f"log2(Exp. Rel. vs {gen_referencia})"),
        showscale=True,
        zmid=0
    )
    fig.add_trace(heatmap, row=2, col=1)
    fig.update_layout(title_text=title, width=1000, height=max(800, len(row_order) * 40), showlegend=False)
    fig.update_xaxes(showticklabels=False, row=1, col=1)
    fig.update_xaxes(showticklabels=True, row=2, col=1, title="Tests")
    fig.update_yaxes(showticklabels=True, row=2, col=1, title="Genes")
    fig.update_yaxes(showticklabels=False, row=2, col=2)
    return fig

# ==== 1. Clustergrama con el gen más estable ====
if genes_separadores:
    heatmap_data = df_norm[df_norm['target'].isin(genes_separadores)].pivot_table(
        index='target', columns='test', values='log2_rel_expr', aggfunc='mean'
    )
    fig1 = create_clustered_heatmap(
        heatmap_data,
        title=f"Clustergrama: Genes separadores (normalización vs {gen_referencia})"
    )
    fig1.show()
    heatmap_data1 = heatmap_data.copy()
else:
    print("No hay suficientes genes separadores para el clustergrama (gen más estable).")

# ================================================================
# 3. Algoritmo de fuerza bruta para selección de candidatos normalizadores
# ================================================================
print("\nIniciando búsqueda de sets óptimos de genes normalizadores (fuerza bruta)...\n")
df = df_total.copy()

def separation_score_for_refs(refs):
    temp = df.copy()
    ref_cts = temp[temp['target'].isin(refs)].groupby('test')['ct'].mean().reset_index()
    ref_cts = ref_cts.rename(columns={'ct': 'ref_ct'})
    temp = temp.merge(ref_cts, on='test')
    temp['delta_ct'] = temp['ct'] - temp['ref_ct']
    temp['log2_rel_expr'] = -temp['delta_ct']
    vals_ctrl = temp[temp['tipo'] == 'Control']['log2_rel_expr']
    vals_case = temp[temp['tipo'] != 'Control']['log2_rel_expr']
    if len(vals_ctrl) < 2 or len(vals_case) < 2:
        return (refs, np.nan)
    t_stat, _ = ttest_ind(vals_ctrl, vals_case, equal_var=False)
    return (refs, abs(t_stat))

scores_indiv = []
for gene in df['target'].unique():
    temp = df.copy()
    ref_ct = temp[temp['target'] == gene][['test', 'ct']].rename(columns={'ct': f'{gene}_ct'})
    temp = temp.merge(ref_ct, on='test')
    temp['delta_ct'] = temp['ct'] - temp[f'{gene}_ct']
    temp['log2_rel_expr'] = -temp['delta_ct']
    vals_ctrl = temp[temp['tipo'] == 'Control']['log2_rel_expr']
    vals_case = temp[temp['tipo'] != 'Control']['log2_rel_expr']
    if len(vals_ctrl) > 1 and len(vals_case) > 1:
        t_stat, _ = ttest_ind(vals_ctrl, vals_case, equal_var=False)
        scores_indiv.append((gene, abs(t_stat)))
scores_indiv = sorted(scores_indiv, key=lambda x: -x[1])
top_N = 12
top_genes = [x[0] for x in scores_indiv[:top_N]]

K = 2  # Cambia K para buscar sets de 2, 3 genes normalizadores
comb_list = list(combinations(top_genes, K))

results = []
with ThreadPoolExecutor() as executor:
    for res in tqdm(executor.map(separation_score_for_refs, comb_list), total=len(comb_list), desc=f"Probando sets de {K} genes"):
        results.append(res)

df_results = pd.DataFrame(results, columns=['genes', 'score']).sort_values('score', ascending=False)
display(df_results.head(10))

# Normalización global usando el/los mejores normalizadores seleccionados
best_genes = df_results.iloc[0]['genes']
temp = df.copy()
ref_cts = temp[temp['target'].isin(best_genes)].groupby('test')['ct'].mean().reset_index()
ref_cts = ref_cts.rename(columns={'ct': 'ref_ct'})
temp = temp.merge(ref_cts, on='test')
temp['delta_ct'] = temp['ct'] - temp['ref_ct']
temp['log2_rel_expr'] = -temp['delta_ct']

# Cálculo de score para todos los genes bajo esta normalización
genes_all = temp['target'].unique()
scores_all = []
for gene in genes_all:
    sub = temp[temp['target'] == gene]
    vals_ctrl = sub[sub['tipo'] == 'Control']['log2_rel_expr']
    vals_case = sub[sub['tipo'] != 'Control']['log2_rel_expr']
    if len(vals_ctrl) > 1 and len(vals_case) > 1:
        t_stat, _ = ttest_ind(vals_ctrl, vals_case, equal_var=False)
        scores_all.append((gene, abs(t_stat)))
df_scores_all = pd.DataFrame(scores_all, columns=['gene', 'score']).sort_values('score', ascending=False)

umbral = 2
genes_separadores = df_scores_all[df_scores_all['score'] > umbral]['gene'].tolist()

df_sep = temp[temp['target'].isin(genes_separadores)]
heatmap_data = df_sep.pivot_table(
    index='target',
    columns='test',
    values='log2_rel_expr',
    aggfunc='mean'
)

# ==== 2. Clustergrama usando el/los mejores genes normalizadores de la búsqueda de fuerza bruta ====
if len(genes_separadores) > 0:

    fig2 = create_clustered_heatmap(
        heatmap_data,
        title=f"Clustergrama: Genes separadores (score > {umbral}) usando normalización por {', '.join(best_genes)}"
    )
    fig2.show()
    heatmap_data2 = heatmap_data.copy()
    display(df_scores_all[df_scores_all['gene'].isin(genes_separadores)])
else:
    print("No hay suficientes genes separadores para el clustergrama (normalización óptima).")