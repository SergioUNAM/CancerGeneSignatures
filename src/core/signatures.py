from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional
import os
import re

import pandas as pd
import numpy as np


@dataclass
class HallmarkConfig:
    hallmark_gmt: str = "gen-sets_GSEA_MSigDB/gsea_hallmarks_formatted.gmt"
    background_gmt: Optional[str] = "gen-sets_GSEA_MSigDB/C5- ontology gene sets.gmt"
    pval_col: str = "Adjusted P-value"


def _load_background_genes(background_gmt: Optional[str]) -> List[str]:
    if not background_gmt:
        return []
    path = os.fspath(background_gmt)
    if not os.path.exists(path):
        return []
    genes: set = set()
    with open(path, "r") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                genes.update(parts[2:])
    return sorted(genes)


def _extract_pmid(link: str) -> Optional[str]:
    m = re.search(r"/(\d+)$", str(link) if link is not None else "")
    return m.group(1) if m else None


def _hallmarks_enrichment(genes: List[str], cfg: HallmarkConfig, background: Optional[List[str]] = None) -> Dict[str, object]:
    if not genes:
        return {}
    try:
        import gseapy as gp  # local import to avoid hard dependency in non-UI flows
        enr = gp.enrich(
            gene_list=genes,
            gene_sets=os.fspath(cfg.hallmark_gmt),
            background=background or None,
            outdir=None,
            verbose=False,
        )
        out: Dict[str, object] = {}
        res = getattr(enr, "results", None)
        if res is None or res.empty:
            return {}
        for _, row in res.iterrows():
            term = str(row.get("Term", "")).replace(" ", "_")
            pval = row.get(cfg.pval_col)
            genes_overlap = str(row.get("Genes", "")).split(";") if isinstance(row.get("Genes"), str) else []
            # keys:
            key_base = f"hallmark_{term}"
            out[f"{key_base}_pvalue"] = pval
            out[f"{key_base}_genes"] = genes_overlap
        return out
    except Exception:
        return {}


def create_signatures(
    df_bib: pd.DataFrame,
    contexto_biologico: str,
    hallmark_cfg: Optional[HallmarkConfig] = None,
) -> pd.DataFrame:
    """
    Construye firmas por (cancer_type, nivel_expresion) desde bibliografía clasificada.
    Requiere columnas: Gene, Link, cancer_type, nivel_expresion, Year, emt_relation, mRNAs_relation.
    """
    if df_bib is None or df_bib.empty:
        return pd.DataFrame()
    hallmark_cfg = hallmark_cfg or HallmarkConfig()

    t = df_bib.copy()
    # pmid
    t["pubmed_id"] = t["Link"].apply(_extract_pmid)
    # expandir cancer_type por coma/espacio
    t["cancer_type"] = t["cancer_type"].astype(str).str.split(", ")
    t = t.explode("cancer_type")
    t["cancer_type"] = t["cancer_type"].astype(str).str.strip()

    # Contadores globales
    total_articulos = len(t)
    global_gene_counts: Dict[str, int] = t["Gene"].value_counts().to_dict()
    # Conteo por cancer y nivel (para denominar en ratio)
    cancer_level_counts = t.groupby(["cancer_type", "nivel_expresion"]).size().to_dict()
    # Filtrado por contexto
    if contexto_biologico == "Cáncer y TEM":
        filt_ctx = t[t.get("emt_relation", False) == True]
    elif contexto_biologico == "Cáncer y micro RNAs":
        filt_ctx = t[t.get("mRNAs_relation", False) == True]
    else:
        filt_ctx = t
    ctx_counts = filt_ctx.groupby(["cancer_type", "nivel_expresion"]).size().to_dict()

    # recency
    decay = 0.05
    current_year = pd.Timestamp.now().year

    # background genes
    background = _load_background_genes(hallmark_cfg.background_gmt)

    rows: List[Dict[str, object]] = []
    for (cancer, nivel), grp in t.groupby(["cancer_type", "nivel_expresion"], dropna=False):
        genes_unicos = grp["Gene"].dropna().astype(str).unique().tolist()
        conteo_list: List[int] = []
        score_list: List[float] = []
        pubmed_map: Dict[str, List[str]] = {}

        for gene in genes_unicos:
            sub = grp[grp["Gene"] == gene]
            n_gene = int(sub.shape[0])
            conteo_list.append(n_gene)
            total_cancer_level = int(cancer_level_counts.get((cancer, nivel), 1))
            n_gene_global = int(global_gene_counts.get(gene, 1))
            # frecuencia relativa contextualizada
            freq_ratio = (n_gene / max(1, total_cancer_level)) / (n_gene_global / max(1, total_articulos))
            # recency factor
            years = pd.to_numeric(sub.get("Year"), errors="coerce").dropna()
            if not years.empty:
                recency = np.exp(-decay * (current_year - years.astype(float))).mean()
            else:
                recency = 1.0
            score = float(freq_ratio) * float(np.log(1 + n_gene)) * float(recency)
            score_list.append(score)
            # pmids
            pmids = sub["pubmed_id"].dropna().astype(str).unique().tolist()
            pubmed_map[gene] = pmids

        # hallmarks
        hm = _hallmarks_enrichment(genes_unicos, hallmark_cfg, background)

        rows.append({
            "cancer_type": cancer,
            "nivel_expresion": nivel,
            "genes": genes_unicos,
            "conteo_articulos_por_gene": conteo_list,
            "estadistica_significancia_bibliografica": score_list,
            "conteo_global_contexto_biologico": int(ctx_counts.get((cancer, nivel), 0)),
            "pubmed_ids": pubmed_map,
            **hm,
        })

    result = pd.DataFrame(rows)
    if not result.empty:
        # fill NaN in hallmark pvalues to 0 (for ML friendliness)
        hm_cols = [c for c in result.columns if c.startswith("hallmark_")]
        if hm_cols:
            result[hm_cols] = result[hm_cols].fillna(0)
    return result

