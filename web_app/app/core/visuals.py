from __future__ import annotations

from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

try:
    from plotly.subplots import make_subplots  # type: ignore
    import plotly.figure_factory as ff  # type: ignore
    from scipy.cluster.hierarchy import linkage, leaves_list  # type: ignore
except Exception:
    # Optional dependencies; functions that use them will raise if missing
    make_subplots = None  # type: ignore
    ff = None  # type: ignore
    linkage = None  # type: ignore
    leaves_list = None  # type: ignore


# Aliases para mostrar nombres amigables de Hallmarks
HALLMARK_ALIASES: Dict[str, str] = {
    "EVADING_GROWTH_SUPPRESSORS": "Evade Growth Suppression",
    "EVADING_IMMUNE_DESTRUCTION": "Evade Immune Destruction",
    "GENOME_INSTABILITY": "Genome Instability",
    "REPLICATIVE_IMMORTALITY": "Replicative Immortality",
    "REPROGRAMMING_ENERGY_METABOLISM": "Metabolic Reprogramming",
    "RESISTING_CELL_DEATH": "Resist Cell Death",
    "SUSTAINED_ANGIOGENESIS": "Sustain Angiogenesis",
    "SUSTAINING_PROLIFERATIVE_SIGNALING": "Sustain Proliferation",
    "TISSUE_INVASION_AND_METASTASIS": "Invasion & Metastasis",
    "TUMOR-PROMOTING_INFLAMMATION": "Tumor Inflammation",
}


def _ensure_list(val) -> List[str]:
    if isinstance(val, list):
        return [str(x).strip() for x in val]
    if isinstance(val, str):
        if val.startswith("[") and val.endswith("]"):
            # podría ser una representación de lista
            import ast
            try:
                lst = ast.literal_eval(val)
                return [str(x).strip() for x in lst] if isinstance(lst, list) else [str(val).strip()]
            except Exception:
                return [s.strip() for s in val.split(",") if s.strip()]
        return [s.strip() for s in val.split(",") if s.strip()]
    try:
        return [str(x).strip() for x in list(val)]
    except Exception:
        return [str(val)] if pd.notna(val) else []


def hallmarks_polar_chart(
    df_signatures: pd.DataFrame,
    cancer_type: str,
    nivel_expresion: str,
    pvalue_threshold: float = 0.05,
) -> go.Figure:
    """Gráfico polar de Hallmarks para un tipo de cáncer y nivel específico.

    Requiere que df_signatures tenga columnas hallmark_*_pvalue y hallmark_*_genes,
    más 'cancer_type', 'nivel_expresion', 'genes' y 'conteo_articulos_por_gene'.
    """
    t = df_signatures.copy()
    filt = (t.get("cancer_type") == cancer_type) & (t.get("nivel_expresion") == nivel_expresion)
    if filt.sum() == 0:
        return go.Figure()
    row = t.loc[filt].iloc[0]

    # Columnas hallmark
    pval_cols = [c for c in row.index if str(c).startswith("hallmark_") and str(c).endswith("_pvalue")]
    hallmarks = [c[len("hallmark_"):-len("_pvalue")] for c in pval_cols]
    pvals = [row[c] for c in pval_cols]
    # Transformaciones y significancia
    logp = [(-np.log10(float(p)) if (isinstance(p, (int, float)) and float(p) > 0) else 0.0) for p in pvals]
    significant = [float(p) < pvalue_threshold if isinstance(p, (int, float)) else False for p in pvals]

    # Colores: significativo → paleta; no significativo → gris
    palette = px.colors.sequential.Aggrnyl
    colors: List[str] = []
    idx = 0
    for sig in significant:
        if sig:
            colors.append(palette[idx % len(palette)])
            idx += 1
        else:
            colors.append("#e0e0e0")

    theta_labels = [HALLMARK_ALIASES.get(h, h).replace("HALLMARK_", "").replace("_", " ") for h in hallmarks]

    fig = go.Figure()
    # Barras traslúcidas
    rgba = [c.replace("rgb", "rgba").replace(")", ",0.4)") if c.startswith("rgb") else c for c in colors]
    fig.add_trace(go.Barpolar(
        r=logp,
        theta=theta_labels,
        marker_color=rgba,
        marker_line_color='rgba(100,100,100,0.5)',
        marker_line_width=1,
        hoverinfo='text',
        hovertext=[f"<b>{lbl}</b><br>-log10(p): {v:.2f}<br>p-value: {p:.4e}" for lbl, v, p in zip(theta_labels, logp, pvals)],
        opacity=1.0,
    ))
    # Marcadores
    fig.add_trace(go.Scatterpolar(
        r=logp,
        theta=theta_labels,
        mode='markers',
        marker=dict(color=colors, size=7, symbol='circle', line=dict(color='rgba(0,0,0,0.5)', width=1)),
        hoverinfo='skip',
    ))
    fig.update_layout(
        title=f"Hallmarks en {cancer_type} — {nivel_expresion}",
        polar=dict(
            bgcolor='white',
            radialaxis=dict(title=dict(text='-log10(p-value)'), gridcolor='#dfe6e9', linecolor='#b2bec3', range=[0, max(logp) * 1.2 if logp else 1.0]),
            angularaxis=dict(rotation=90, direction='clockwise', gridcolor='#dfe6e9', linecolor='#b2bec3'),
        ),
        template='plotly_white',
        showlegend=False,
        height=650,
        margin=dict(t=80, b=60, l=60, r=60),
    )
    return fig


def fingerprint_heatmap(
    df_signatures: pd.DataFrame,
    cancer_type: str,
    nivel_expresion: str,
    colorscale: str | List = "Viridis",
) -> go.Figure:
    """Heatmap genes × hallmarks con z = -log10(p) para un nivel.

    Construye una matriz donde cada gen asociado a un hallmark toma el valor -log10(p) del término.
    """
    if df_signatures is None or df_signatures.empty:
        return go.Figure()
    base = df_signatures[(df_signatures.get('cancer_type') == cancer_type) & (df_signatures.get('nivel_expresion') == nivel_expresion)]
    if base.empty:
        return go.Figure()
    row = base.iloc[0]
    genes_firma = _ensure_list(row.get('genes'))

    registros: List[Dict[str, object]] = []
    for col in row.index:
        s = str(col)
        if s.startswith('hallmark_') and s.endswith('_genes'):
            base_term = s[len('hallmark_'):-len('_genes')]
            pval = row.get(f"hallmark_{base_term}_pvalue", 1.0)
            pval = float(pval) if isinstance(pval, (int, float, str)) else 1.0
            try:
                pval = float(pval)
            except Exception:
                pval = 1.0
            logp = -np.log10(pval) if pval > 0 else 0.0
            genes_hm = _ensure_list(row.get(col))
            for g in genes_firma:
                if g in genes_hm:
                    registros.append({"gen": g, "hallmark": base_term.replace('HALLMARK_', '').replace('_', ' ').title(), "-log10(p)": logp})
    if not registros:
        return go.Figure()
    df_plot = pd.DataFrame(registros)
    matriz = df_plot.pivot(index="gen", columns="hallmark", values="-log10(p)")
    matriz = matriz.fillna(0)

    # Ordenamiento jerárquico si hay múltiples filas/columnas
    try:
        if linkage is not None and leaves_list is not None and ff is not None:
            if matriz.shape[0] > 1:
                row_linkage = linkage(matriz.values, method='ward')
                row_order = leaves_list(row_linkage)
                matriz = matriz.iloc[row_order[::-1]]
            if matriz.shape[1] > 1:
                col_linkage = linkage(matriz.values.T, method='ward')
                col_order = leaves_list(col_linkage)
                matriz = matriz.iloc[:, col_order]
    except Exception:
        pass

    hover_text = np.array([[
        (f"<b>Gen:</b> {r}<br><b>Hallmark:</b> {c}<br><b>-log10(p):</b> {matriz.loc[r, c]:.2f}")
        for c in matriz.columns
    ] for r in matriz.index])

    fig = go.Figure(go.Heatmap(
        z=matriz.values,
        x=matriz.columns,
        y=matriz.index,
        colorscale=colorscale,
        colorbar=dict(title="-log10(p)"),
        hoverinfo="text",
        text=hover_text,
        zmin=0,
        zmax=float(np.nanmax(matriz.values)) if matriz.size else 1.0,
    ))
    fig.update_layout(
        title=dict(text=f"Fingerprint Genómico — {cancer_type} ({nivel_expresion})", x=0.5),
        template='plotly_white',
        height=700,
        margin=dict(l=120, r=40, t=80, b=100),
        xaxis=dict(title='Hallmarks', tickangle=45),
        yaxis=dict(title='Genes')
    )
    return fig


def _create_clustered_heatmap(heatmap_data: pd.DataFrame, title: str) -> go.Figure:
    """Replica del clustergrama con dendrogramas para una matriz genes×tests.

    heatmap_data: índice=genes, columnas=tests, valores=log2_rel_expr (centrado en 0).
    """
    if ff is None or linkage is None:
        raise RuntimeError("Se requieren scipy y plotly.figure_factory para el clustergrama.")

    dendro_col = ff.create_dendrogram(
        heatmap_data.T,
        orientation='bottom',
        labels=heatmap_data.columns.tolist(),
        linkagefun=lambda x: linkage(x, method='average', metric='euclidean'),
    )
    col_order = dendro_col['layout']['xaxis']['ticktext']

    dendro_row = ff.create_dendrogram(
        heatmap_data,
        orientation='right',
        labels=heatmap_data.index.tolist(),
        linkagefun=lambda x: linkage(x, method='average', metric='euclidean'),
    )
    row_order = dendro_row['layout']['yaxis']['ticktext']

    clustered = heatmap_data.loc[row_order, col_order]
    custom_colorscale = [
        [0.0, '#2c7bb6'],  # azul
        [0.5, '#ffffb2'],  # amarillo pastel
        [1.0, '#d7191c'],  # rojo
    ]

    fig = make_subplots(
        rows=2, cols=2,
        shared_xaxes=True, shared_yaxes=True,
        vertical_spacing=0.02, horizontal_spacing=0.02,
        column_widths=[0.8, 0.2], row_heights=[0.10, 1],
        specs=[[{"type": "scatter", "colspan": 2}, None],
               [{"type": "heatmap"}, {"type": "scatter"}]],
    )
    for tr in dendro_col['data']:
        fig.add_trace(tr, row=1, col=1)
    for tr in dendro_row['data']:
        fig.add_trace(tr, row=2, col=2)
    fig.add_trace(go.Heatmap(
        z=clustered.values,
        x=col_order,
        y=row_order,
        colorscale=custom_colorscale,
        showscale=True,
        zmid=0,
        colorbar=dict(title="log2(Exp. Rel.)"),
    ), row=2, col=1)
    fig.update_layout(title_text=title, width=1000, height=800, showlegend=False)
    fig.update_xaxes(showticklabels=False, row=1, col=1)
    fig.update_xaxes(showticklabels=True, row=2, col=1, title="Tests")
    fig.update_yaxes(showticklabels=True, row=2, col=1, title="Genes")
    fig.update_yaxes(showticklabels=False, row=2, col=2)
    return fig


def clustered_expression_by_level(
    controles_df: pd.DataFrame,
    muestras_df: pd.DataFrame,
    lista_sobre: List[str],
    lista_sub: List[str],
    metodo_elegido: str = "gen de referencia",
    gen_referencia: Optional[str] = None,
) -> Dict[str, go.Figure]:
    """Genera clustergramas para 'Sobreexpresados' y 'Subexpresados'.

    - metodo_elegido: "gen de referencia" o "promedios".
    - gen_referencia: requerido si se usa "gen de referencia" y debe existir en los datos.
    """
    if controles_df is None or muestras_df is None or (controles_df.empty and muestras_df.empty):
        return {}

    df_total = pd.concat([controles_df, muestras_df], ignore_index=True)
    df_total = df_total.copy()
    df_total['ct'] = pd.to_numeric(df_total['ct'], errors='coerce')

    metodo = metodo_elegido.strip().lower()
    if metodo == "gen de referencia":
        if not gen_referencia:
            raise ValueError("Se requiere 'gen_referencia' para normalizar por gen de referencia.")
        ref_ct = df_total[df_total['target'] == gen_referencia][['test', 'ct']].rename(columns={'ct': f'{gen_referencia}_ct'})
        if ref_ct.empty:
            # fallback a promedios si no existe el gen de referencia
            metodo = "promedios"
        else:
            df_norm = df_total.merge(ref_ct, on='test', how='left')
            df_norm['delta_ct'] = df_norm['ct'] - df_norm[f'{gen_referencia}_ct']
            df_norm['log2_rel_expr'] = -df_norm['delta_ct']
    if metodo == "promedios":
        df_norm = df_total.copy()
        df_norm['ct_promedio_test'] = df_norm.groupby('test')['ct'].transform('mean')
        df_norm['delta_ct_promedio'] = df_norm['ct'] - df_norm['ct_promedio_test']
        df_norm['log2_rel_expr'] = -df_norm['delta_ct_promedio']

    heatmap_base = df_norm.pivot_table(index='target', columns='test', values='log2_rel_expr', aggfunc='mean')

    figs: Dict[str, go.Figure] = {}
    if lista_sobre:
        sub = heatmap_base.loc[[g for g in heatmap_base.index if g in set(lista_sobre)]]
        if not sub.empty:
            figs['Sobreexpresados'] = _create_clustered_heatmap(sub, title="Expresión relativa — Sobreexpresados")
    if lista_sub:
        sub2 = heatmap_base.loc[[g for g in heatmap_base.index if g in set(lista_sub)]]
        if not sub2.empty:
            figs['Subexpresados'] = _create_clustered_heatmap(sub2, title="Expresión relativa — Subexpresados")
    return figs

