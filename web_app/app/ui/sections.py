from __future__ import annotations

import io
from typing import Optional

import pandas as pd
import streamlit as st

from app.services.heuristics import (
    compute_heuristic_summary,
    build_heatmap_data,
    build_sankey_data,
    merge_with_expression,
    build_function_long,
    compute_function_counts,
    map_functions_to_hallmarks,
)
from app.services.visuals import (
    build_heuristic_heatmap,
    build_heuristic_sankey,
)
from app.services.nlp import (
    prepare_corpus,
    extract_texts,
    analyze_texts,
)

try:
    from src.core.bibliography import filter_bibliography_by_cancer  # type: ignore
except Exception:  # pragma: no cover - dependencia opcional
    filter_bibliography_by_cancer = None  # type: ignore


@st.cache_data(show_spinner=False, ttl=3600)
def _prepare_corpus_cached(
    df_csv: str,
    cancer_label: str,
    gene: Optional[str],
    apply_filter: bool,
) -> dict:
    try:
        df = pd.read_csv(io.StringIO(df_csv))  # type: ignore[name-defined]
    except Exception:
        return {}
    try:
        info = prepare_corpus(
            df,
            cancer_label=cancer_label,
            gene=gene,
            apply_filter=apply_filter,
        )
    except ValueError as exc:
        return {"error": str(exc)}
    return {
        "filtered_csv": info.filtered.to_csv(index=False),
        "gene_csv": info.gene_filtered.to_csv(index=False),
        "title_col": info.title_column,
        "abstract_col": info.abstract_column or "",
        "total_all": info.total_all,
        "total_gene": info.total_gene,
        "total_filtered": info.total_filtered,
    }


def render_heuristics_section(app_state, cancer_type: str) -> None:
    st.markdown("---")
    st.header("Heurística bibliográfica y relación con niveles/firmas")

    focused = st.session_state.get("bibliografia_clasificada")
    if not (isinstance(focused, pd.DataFrame) and not focused.empty):
        st.info("Primero ejecuta la búsqueda en PubMed y la clasificación para generar la heurística.")
        return

    summary = compute_heuristic_summary(focused)
    if summary.empty:
        st.info("No se pudieron etiquetar relaciones con las heurísticas actuales.")
        return

    ordered_summary = summary.sort_values('heuristic_summary').reset_index(drop=True)
    st.dataframe(ordered_summary)
    st.download_button(
        label="Descargar resumen heurístico (CSV)",
        data=ordered_summary.to_csv(index=False),
        file_name="resumen_heuristico_relaciones.csv",
        mime="text/csv",
        use_container_width=True,
    )

    heatmap_fig = build_heuristic_heatmap(build_heatmap_data(summary))
    if heatmap_fig is not None:
        st.plotly_chart(heatmap_fig, use_container_width=True)

    sankey_fig = build_heuristic_sankey(build_sankey_data(summary))
    if sankey_fig is not None:
        st.plotly_chart(sankey_fig, use_container_width=True)

    try:
        import numpy as np
        import plotly.express as px

        df_levels = st.session_state.get('df_expr')
        merged = merge_with_expression(summary, df_levels, bool(app_state.exclude_stable))
        long_df = build_function_long(merged)
        if long_df.empty:
            st.info("Sin datos combinados de heurística y niveles de expresión.")
        else:
            st.subheader("Niveles × funciones")
            colm1, colm2, colm3 = st.columns([1, 1, 1])
            with colm1:
                agg_mode = st.selectbox("Métrica", ["conteo", "score"], index=1, key="agg_mode_heu")
            with colm2:
                top_funcs = st.number_input("Top funciones", min_value=3, max_value=20, value=8, step=1, key="top_funcs_heu")
            with colm3:
                norm_scores = st.checkbox("Normalizar por columna", value=True, key="norm_scores_heu")

            if agg_mode == 'conteo':
                pivot = long_df.pivot_table(index='nivel_expresion', columns='funcion', values='flag', aggfunc='sum', fill_value=0)
            else:
                pivot = long_df.pivot_table(index='nivel_expresion', columns='funcion', values='score', aggfunc='sum', fill_value=0.0)
            totals = pivot.sum(axis=0).sort_values(ascending=False)
            top_cols = totals.head(int(top_funcs)).index.tolist()
            pivot_top = pivot[top_cols]
            if norm_scores:
                pivot_top = (pivot_top - pivot_top.min()) / (pivot_top.max() - pivot_top.min() + 1e-9)
            fig_lv = px.imshow(
                pivot_top.values,
                labels=dict(x='Función', y='Nivel expresión', color='Valor'),
                x=top_cols,
                y=pivot_top.index.tolist(),
                aspect='auto',
                title=f"Niveles × Funciones ({agg_mode})",
            )
            st.plotly_chart(fig_lv, use_container_width=True)

            view_mode_lvl = st.radio("Vista complementaria", ["Tabla detallada", "Treemap"], horizontal=True, index=0)
            if view_mode_lvl == "Tabla detallada":
                st.dataframe(long_df)
            else:
                trea = long_df.merge(merged[['Gene', 'nivel_expresion', 'fold_change']].drop_duplicates('Gene'), on='Gene', how='left')
                trea['log2fc_abs'] = np.log2(trea['fold_change'].clip(lower=1e-12)).abs()
                trea['weight'] = trea['score'] * (1.0 + trea['log2fc_abs'])
                fig_tree = px.treemap(
                    trea,
                    path=['nivel_expresion', 'funcion', 'Gene'],
                    values='weight',
                    color='score',
                    color_continuous_scale='Viridis',
                    title='Treemap: Nivel → Función → Gen (peso = score×|log2FC|)'
                )
                st.plotly_chart(fig_tree, use_container_width=True)

            df_sigs_v = st.session_state.get('df_signatures')
            if isinstance(df_sigs_v, pd.DataFrame) and not df_sigs_v.empty:
                st.subheader("Funciones ↔ Hallmarks (firmas)")
                try:
                    map_hm = map_functions_to_hallmarks(summary, df_sigs_v, cancer_type)
                    if not map_hm.empty:
                        mf = pd.merge(long_df, map_hm, on='Gene', how='inner')
                        agg_hm = mf.groupby(['funcion', 'hallmark'])['flag'].sum().reset_index(name='n')
                        import plotly.graph_objects as go
                        funcs = agg_hm['funcion'].unique().tolist()
                        terms = agg_hm['hallmark'].unique().tolist()
                        nodes = funcs + terms
                        idx = {n: i for i, n in enumerate(nodes)}
                        src = [idx[f] for f in agg_hm['funcion']]
                        tgt = [idx[t] for t in agg_hm['hallmark']]
                        val = agg_hm['n'].astype(float).tolist()
                        fig_fh = go.Figure(go.Sankey(node=dict(label=nodes, pad=12, thickness=12), link=dict(source=src, target=tgt, value=val)))
                        fig_fh.update_layout(margin=dict(l=10, r=10, t=10, b=10), title='Funciones → Hallmarks (conteos)')
                        st.plotly_chart(fig_fh, use_container_width=True)
                except Exception:
                    pass
    except Exception as exc:
        st.warning(f"No se pudo ejecutar la heurística: {exc}")


def render_google_nlp_section(app_state, cancer_type: str, google_cfg) -> None:
    st.markdown("---")
    st.header("Insights de la literatura (Google NLP)")

    try:
        sec_gkey = st.secrets["GOOGLE_NLP_API_KEY"]  # type: ignore[attr-defined]
    except Exception:
        sec_gkey = ""
    default_google_key = sec_gkey or (google_cfg.api_key or "")

    ins_col1, ins_col2, ins_col3 = st.columns([2, 1, 1])
    with ins_col1:
        g_api_key = st.text_input(
            "Google NLP API Key",
            value=default_google_key,
            type="password",
            help="Se recomienda usar st.secrets['GOOGLE_NLP_API_KEY'] o la variable de entorno GOOGLE_NLP_API_KEY.",
        )
    with ins_col2:
        n_docs = st.number_input("Máx. artículos a analizar", min_value=5, max_value=200, value=30, step=5)
    with ins_col3:
        lang_sel = st.selectbox("Idioma (opcional)", options=["auto", "es", "en"], index=0, help="Déjalo en 'auto' para detección automática.")

    apply_cancer_filter = st.checkbox(
        "Filtrar por el tipo de cáncer seleccionado",
        value=True,
        help="Desactívalo si tu CSV incluye artículos de otros cánceres y quieres analizarlos igualmente.",
    )
    run_insights = st.button("Generar insights (Google NLP)")

    uploaded_bib_ins = st.file_uploader(
        "Subir bibliografía (CSV) — columnas esperadas: Title, Abstract, Gene (opcional)",
        type=["csv"],
        key="upl_bib_nlp_csv",
        help="Si subes un archivo aquí, se usará para los insights en lugar de la bibliografía de la sesión.",
    )

    gene_selected_input: Optional[str] = None
    gene_source_df = None

    def _read_uploaded_csv(uploaded_file):
        if uploaded_file is None:
            return None
        for attempt in range(3):
            try:
                uploaded_file.seek(0)
            except Exception:
                pass
            try:
                if attempt == 0:
                    return pd.read_csv(uploaded_file)
                if attempt == 1:
                    return pd.read_csv(uploaded_file, sep=None, engine='python')
                return pd.read_csv(uploaded_file, sep=';', encoding='utf-8-sig')
            except Exception:
                continue
        return None

    try:
        if uploaded_bib_ins is not None:
            gene_source_df = _read_uploaded_csv(uploaded_bib_ins)
        else:
            gene_source_df = st.session_state.get("bibliografia_clasificada")
    except Exception:
        gene_source_df = None

    if isinstance(gene_source_df, pd.DataFrame) and not gene_source_df.empty:
        gene_col = None
        for c in ["Gene", "gene", "Symbol", "symbol", "target"]:
            if c in gene_source_df.columns:
                gene_col = c
                break
        if gene_col:
            uniq_genes = (
                gene_source_df[gene_col].dropna().astype(str).str.strip().unique().tolist()
            )
            uniq_genes = sorted([g for g in uniq_genes if g])
            if uniq_genes:
                gene_selected_input = st.selectbox("Gen a analizar (desde archivo/datos)", options=uniq_genes, index=0)
        if not gene_selected_input:
            gene_selected_input = st.text_input("Gen a analizar (manual si no hay columna)", value="", placeholder="Ej.: TP53")
    else:
        gene_selected_input = st.text_input("Gen a analizar", value="", placeholder="Ej.: TP53")

    if not run_insights:
        return

    if uploaded_bib_ins is not None:
        if isinstance(gene_source_df, pd.DataFrame) and not gene_source_df.empty:
            bib_df = gene_source_df.copy()
        else:
            try:
                bib_df = _read_uploaded_csv(uploaded_bib_ins)
            except Exception as exc:
                st.error(f"No se pudo leer el CSV subido: {exc}")
                bib_df = None
        if bib_df is None:
            st.error("No se pudo leer el CSV subido (prueba con separador ';' o codificación UTF-8).")
            return
    else:
        bib_df = st.session_state.get("bibliografia_clasificada")

    if not (isinstance(bib_df, pd.DataFrame) and not bib_df.empty):
        st.info("Primero genera o carga la bibliografía clasificada en la sección anterior.")
        return

    if not (g_api_key or google_cfg.api_key):
        st.warning("Falta la API Key de Google NLP. Ingresa una o configúrala en el entorno/secrets.")
        return

    with st.spinner("Analizando textos con Google NLP…"):
        cache = _prepare_corpus_cached(
            bib_df.to_csv(index=False),
            cancer_label=cancer_type,
            gene=gene_selected_input or None,
            apply_filter=bool(apply_cancer_filter),
        )

    if not cache:
        st.info("No hay textos válidos en el CSV (requiere columnas Title/Abstract con contenido).")
        return

    if cache.get("error"):
        st.error(cache["error"])
        return

    filtered_df = pd.read_csv(io.StringIO(cache["filtered_csv"])) if cache.get("filtered_csv") else pd.DataFrame()
    total_docs_all = int(cache.get("total_all", 0))
    total_docs_gene = int(cache.get("total_gene", 0))
    total_docs_final = int(cache.get("total_filtered", 0))

    st.caption(
        f"Artículos: {total_docs_all} en total → {total_docs_gene} para el gen seleccionado"
        + (f" → {total_docs_final} tras filtrar por '{cancer_type}'" if apply_cancer_filter else "")
    )

    texts = extract_texts(
        filtered_df,
        title_column=cache.get("title_col", "Title"),
        abstract_column=cache.get("abstract_col") or None,
        limit=int(n_docs),
    )

    if not texts:
        if apply_cancer_filter:
            st.info("No hay textos para analizar con el filtro de cáncer activo. Desactívalo para analizar todo tu CSV.")
        else:
            st.info("No hay textos válidos en el CSV (requiere columnas Title/Abstract con contenido).")
        return

    try:
        glang = None if lang_sel == "auto" else lang_sel
        res = analyze_texts(
            texts,
            api_key=g_api_key or google_cfg.api_key,
            language=glang,
        ) or {}
    except Exception as exc:
        st.error(f"Error llamando Google NLP: {exc}")
        res = {}

    if not res:
        return

    try:
        gname = (gene_selected_input or '').strip()
        sent = res.get("sentiment", {}) or {}
        cats = res.get("categories", []) or []
        ents = res.get("entities", []) or []
        pos_list = [c['category'] for c in cats if c.get('confidence_sum', 0) > 0.2]
        top_cats = [c['category'] for c in cats][:5]
        top_ent_names = [e['name'] for e in ents[:5]] if ents else []
        avg_s = float(sent.get('avg_score', 0))
        avg_m = float(sent.get('avg_magnitude', 0))
        n_docs_used = int(sent.get('n', 0))
        partes = []
        partes.append(
            f"Se revisaron {n_docs_used} artículos" + (f" sobre {gname}" if gname else "") + (f" en {cancer_type}" if apply_cancer_filter else "") + "."
        )
        if top_cats:
            partes.append(f"Los temas dominantes fueron: {', '.join(top_cats)}.")
        if top_ent_names:
            partes.append(f"Las entidades más relevantes incluyeron {', '.join(top_ent_names)}.")
        if pos_list:
            partes.append(f"Se observaron asociaciones positivas destacadas con {', '.join(pos_list)}.")
        neg_list = [e['name'] for e in ents if e.get('sentiment_score', 0) < -0.1]
        if neg_list:
            partes.append(f"También emergieron señales negativas en torno a {', '.join(neg_list)}.")
        tono = 'neutral'
        if avg_s > 0.1:
            tono = 'ligeramente positivo'
        elif avg_s < -0.1:
            tono = 'ligeramente negativo'
        partes.append(f"El tono global fue {tono} (score {avg_s:.2f}, magnitud {avg_m:.2f}).")
        st.info(" ".join(partes))
    except Exception:
        pass

    try:
        df_sub = bib_df.copy()
        gcand = None
        for c in ["Gene", "gene", "Symbol", "symbol", "target"]:
            if c in df_sub.columns:
                gcand = c
                break
        if gcand and gene_selected_input:
            df_sub = df_sub[df_sub[gcand].astype(str).str.strip().str.lower() == gene_selected_input.strip().lower()]
        if apply_cancer_filter and filter_bibliography_by_cancer is not None:
            try:
                df_sub = filter_bibliography_by_cancer(df_sub, cancer_type)
            except Exception:
                pass
        summ = compute_heuristic_summary(df_sub)
        if not summ.empty:
            row = summ.iloc[0]
            effects = []
            if float(row.get('upregulated_score', 0)) > float(row.get('downregulated_score', 0)):
                effects.append('tendencia a sobreexpresión')
            elif float(row.get('downregulated_score', 0)) > float(row.get('upregulated_score', 0)):
                effects.append('tendencia a subexpresión')
            prog = 'desfavorable' if float(row.get('prognosis_bad_score', 0)) > float(row.get('prognosis_good_score', 0)) else ('favorable' if float(row.get('prognosis_good_score', 0)) > 0 else 'incierta')
            flags = [k for k in ['proliferation','invasion','migration','metastasis','drug_resistance','drug_sensitivity','apoptosis','emt_related'] if row.get(k, 0) > 0]
            texto = (
                f"Heurística complementaria: para {gene_selected_input or 'el gen seleccionado'} se observa {', '.join(effects) if effects else 'un efecto de expresión incierto'}, "
                f"y una señal pronóstica {prog}. Funciones implicadas: {', '.join(flags) if flags else 'no concluyentes'}."
            )
            st.caption(texto)
    except Exception:
        pass

    sent = res.get("sentiment", {}) or {}
    s1, s2, s3 = st.columns(3)
    with s1:
        st.metric("Sentimiento promedio", f"{float(sent.get('avg_score', 0.0)):.2f}")
    with s2:
        st.metric("Magnitud promedio", f"{float(sent.get('avg_magnitude', 0.0)):.2f}")
    with s3:
        st.metric("Artículos analizados", int(sent.get('n', 0)))

    ents = pd.DataFrame(res.get("entities") or [])
    if not ents.empty:
        st.subheader("Entidades más relevantes")
        show = ents.rename(columns={"salience_sum": "saliencia_acum", "mentions": "menciones"})
        st.dataframe(show.head(50))

    entsent = pd.DataFrame(res.get("entity_sentiment") or [])
    if not entsent.empty:
        st.subheader("Entidades con sentimiento (promedios)")
        show2 = entsent.rename(columns={
            "avg_sentiment": "sent_prom",
            "avg_magnitude": "magn_prom",
            "mentions": "menciones",
        })
        st.dataframe(show2.head(50))

    cats = pd.DataFrame(res.get("categories") or [])
    if not cats.empty:
        st.subheader("Categorías (classifyText)")
        try:
            import plotly.express as px

            figc = px.bar(cats.head(20), x="confidence_sum", y="category", orientation="h", title="Top categorías")
            st.plotly_chart(figc, use_container_width=True)
        except Exception:
            pass
        st.dataframe(cats)
