from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List, Dict, Tuple, Optional, Callable, Any
import re
import difflib
from functools import lru_cache
import os
import time
import re

import pandas as pd

try:
    # Lazy import; Biopython is required for PubMed access
    from Bio import Entrez, Medline  # type: ignore
except Exception as _e:
    Entrez = None  # type: ignore
    Medline = None  # type: ignore


# --- Keyword dictionaries (concise, extensible) ---

GENERAL_CANCER_TERMS = [
    "cancer", "tumor", "tumour", "neoplasm", "metastasis", "carcinoma"
]

from typing import Dict, List

DEFAULT_CANCER_TYPES: Dict[str, List[str]] = {
    "Breast Cancer": [
        "breast", "mammary", "mamma", "triple negative breast",
        "er-positive", "pr-positive", "ductal carcinoma",
        "lobular carcinoma", "estrogen receptor", "progesterone receptor"
    ],
    "Melanoma": [
        "melanoma", "cutaneous melanoma", "skin melanoma",
        "uveal", "uveal melanoma", "ocular melanoma",
        "acral melanoma", "braf", "braf mutation"
    ],
    "Colon Cancer": [
        "colon", "colorectal", "rectal", "crc", "intestinal",
        "intestine", "colon adenocarcinoma", "rectal adenocarcinoma",
        "lynch syndrome", "hereditary colon cancer"
    ],
    "Lung Cancer": [
        "lung", "pulmonary", "bronchial", "nsclc", "sclc",
        "non-small cell lung cancer", "small cell lung cancer",
        "egfr mutation", "alk rearrangement",
        "lung adenocarcinoma", "lung squamous cell carcinoma", "lung squamous"
    ],
    "Prostate Cancer": [
        "prostate", "prostatic", "psa", "androgen",
        "androgen receptor", "ar-positive",
        "prostate-specific antigen", "gleason score",
        "crpc", "castration-resistant", "prostate carcinoma"
    ],
    "Pancreatic Cancer": [
        "pancreatic", "pancreas", "pdac", "ductal adenocarcinoma",
        "pancreatic carcinoma", "neuroendocrine tumor",
        "pancreatic neuroendocrine", "islet cell tumor", "whipple procedure"
    ],
    "Leukemia": [
        "leukemia", "leukaemia", "aml", "cml", "cll",
        "acute myeloid leukemia", "acute lymphoblastic leukemia",
        "chronic lymphocytic leukemia", "chronic myeloid leukemia",
        "philadelphia chromosome", "bcr-abl",
        "mpn", "myeloproliferative neoplasm"
    ],
    "Lymphoma": [
        "lymphoma", "hodgkin", "non-hodgkin", "b-cell", "t-cell",
        "nhl", "dlbcl", "diffuse large b-cell lymphoma",
        "follicular lymphoma", "mantle cell lymphoma",
        "burkitt lymphoma", "cutaneous lymphoma"
    ],
    "Ovarian Cancer": [
        "ovarian", "ovary", "ovarian carcinoma",
        "high-grade serous carcinoma", "hgsc", "brca1", "brca2",
        "ca125", "clear cell carcinoma", "granulosa cell tumor"
    ],
    "Cervical Cancer": [
        "cervical", "cervix", "hpv", "human papillomavirus",
        "cervical carcinoma", "cervical squamous cell carcinoma",
        "cervical adenocarcinoma", "pap smear",
        "cervical dysplasia", "cervical squamous"
    ],
    "Renal Cancer": [
        "renal", "kidney", "rcc", "renal carcinoma",
        "clear cell carcinoma", "papillary renal carcinoma",
        "chromophobe renal carcinoma", "von hippel-lindau", "vhl"
    ],
    "Bladder Cancer": [
        "bladder", "urothelial", "bladder carcinoma",
        "transitional cell carcinoma", "non-muscle invasive",
        "muscle invasive", "nmibc", "mibc"
    ],
    "Glioblastoma": [
        "glioblastoma", "gbm", "glioma", "astrocytoma",
        "oligodendroglioma", "brain tumor",
        "idh", "idh mutation", "mgmt methylation", "temozolomide"
    ],
    "Thyroid Cancer": [
        "thyroid", "thyroid carcinoma", "papillary thyroid",
        "follicular thyroid", "medullary thyroid",
        "anaplastic thyroid", "braf", "braf mutation",
        "ret", "ret mutation"
    ],
    "Gastric Cancer": [
        "gastric", "stomach", "gastric carcinoma",
        "diffuse gastric cancer", "intestinal metaplasia",
        "epstein-barr virus", "helicobacter pylori",
        "h. pylori", "pylori"
    ],
    "Esophageal Cancer": [
        "esophageal", "esophagus", "esophageal carcinoma",
        "barrett esophagus", "esophageal adenocarcinoma",
        "esophageal squamous cell carcinoma"
    ],
    "Skin Cancer": [
        "skin cancer", "basal cell carcinoma", "bcc",
        "squamous cell carcinoma of the skin", "scc",
        "cutaneous carcinoma", "actinic keratosis",
        "uv radiation", "squamous-cell carcinoma of the skin"
    ],
    "Testicular Cancer": [
        "testicular", "testis", "germ cell tumor",
        "seminoma", "non-seminoma", "beta-hcg",
        "alpha-fetoprotein", "afp"
    ],
    "Endometrial Cancer": [
        "endometrial", "uterine", "endometrial carcinoma",
        "serous carcinoma", "lynch syndrome", "womb cancer"
    ],
    "Mesothelioma": [
        "mesothelioma", "pleural mesothelioma",
        "peritoneal mesothelioma", "asbestos",
        "asbestos exposure", "asbestosis"
    ],
    "Head and Neck Cancer": [
        "head and neck", "laryngeal", "pharyngeal",
        "oral squamous", "oral squamous cell carcinoma", "oscc",
        "pharyngeal carcinoma", "laryngeal carcinoma",
        "hpv-positive", "hpv-negative", "nasopharyngeal carcinoma",
        "laryngeal squamous cell carcinoma",
        "naso-oropharyngeal carcinoma", "tongue squamous"
    ],
    "Bone Cancer": [
        "bone cancer", "osteosarcoma", "chondrosarcoma",
        "ewing", "ewing sarcoma", "primary bone tumor",
        "metastatic bone cancer", "paget disease",
        "bone lesion", "osteogenic sarcoma"
    ],
    "Hepatocellular Carcinoma": [
        "hepatocellular", "hcc", "liver cancer", "hepatic carcinoma",
        "hepatoma", "cirrhosis", "hepatitis b", "hbv",
        "hepatitis c", "hcv", "fibrolamellar carcinoma",
        "afp", "alpha-fetoprotein"
    ]
}

EMT_KEYWORDS = [
    "epithelial-to-mesenchymal transition", "epithelial mesenchymal transition", "emt",
    "mesenchymal-epithelial transition", "epithelial-mesenchymal",
    "transición epitelio mesénquima", "epitelio-mesenquimal",
]

MICRORNA_KEYWORDS = [
    "microrna", "mirna", "mirnas", "miRNA", "miRNAs", "microRNA", "microRNAs",
    "non-coding rna", "noncoding rna", "ncRNA", "ncRNAs",
]


@dataclass
class PubMedConfig:
    email: str
    api_key: Optional[str] = None
    max_per_gene: int = 100
    sleep_between: float = 0.34  # ~3 req/s without API key


def _require_entrez_config(email: Optional[str] = None, api_key: Optional[str] = None) -> PubMedConfig:
    """Obtiene configuración de Entrez desde argumentos o variables de entorno.
    Da prioridad a parámetros explícitos para evitar depender de `os.environ` en la app.
    """
    em = email or os.getenv("NCBI_EMAIL")
    if not em:
        raise RuntimeError("NCBI_EMAIL no configurado (parámetro o entorno).")
    key = api_key or os.getenv("NCBI_API_KEY")
    return PubMedConfig(email=em, api_key=key)


def _setup_entrez(cfg: PubMedConfig) -> None:
    if Entrez is None:
        raise RuntimeError("Biopython no está instalado. Agrega 'biopython' a requirements.")
    Entrez.email = cfg.email
    if cfg.api_key:
        Entrez.api_key = cfg.api_key


@lru_cache(maxsize=1024)
def _compile_word_regex(words_t: tuple) -> re.Pattern:
    parts = [rf"\b{re.escape(w).replace(r'\-', '[- ]?')}\b" for w in words_t if w]
    return re.compile("|".join(parts), flags=re.IGNORECASE) if parts else re.compile(r"$.^")


def _compile_kw_regex(words: List[str]) -> re.Pattern:
    return _compile_word_regex(tuple(w for w in words if w))


def _detect_any(text: str, rx: re.Pattern) -> bool:
    return bool(rx.search((text or "").lower()))


def search_pubmed_by_genes(
    df_genes: pd.DataFrame,
    symbol_col: str = "target",
    ensembl_col: str = "ensembl_id",
    selected_context: str = "Cáncer y TEM",
    max_per_gene: int = 100,
    progress: Optional[Callable[[int, int, str], None]] = None,
    logger: Optional[Any] = None,
    email: Optional[str] = None,
    api_key: Optional[str] = None,
) -> pd.DataFrame:
    """
    Busca artículos en PubMed por gen (símbolo) y filtra por contexto/cáncer de forma ligera.
    Devuelve columnas: Gene, Ensembl_ID, Title, Year, Abstract, Link
    """
    cfg = _require_entrez_config(email=email, api_key=api_key)
    _setup_entrez(cfg)
    max_n = int(max_per_gene or cfg.max_per_gene)

    if selected_context == "Cáncer y TEM":
        ctx_kw = EMT_KEYWORDS
    else:
        ctx_kw = MICRORNA_KEYWORDS

    rx_ctx = _compile_kw_regex(ctx_kw)
    rx_cancer = _compile_kw_regex(GENERAL_CANCER_TERMS)

    rows: List[Dict[str, object]] = []
    total = len(df_genes)
    for idx, (_, row) in enumerate(df_genes.iterrows(), start=1):
        gene = str(row.get(symbol_col, "")).strip()
        ensembl_id = str(row.get(ensembl_col, "")).strip()
        if not gene:
            continue
        if progress:
            try:
                progress(idx, total, gene)
            except Exception:
                pass
        if logger:
            try:
                logger.info(f"PubMed: [{idx}/{total}] {gene}")
            except Exception:
                pass
        # ESearch query: gene AND (cancer terms OR context terms)
        # Keep the query simple and broad; refine on client side
        query = f"{gene} AND ({' OR '.join(GENERAL_CANCER_TERMS + ctx_kw)})"
        try:
            ids: List[str] = []
            retmax = min(200, max_n)
            retstart = 0
            while len(ids) < max_n:
                h = Entrez.esearch(db="pubmed", term=query, retmax=retmax, retstart=retstart)
                record = Entrez.read(h)
                h.close()
                batch = record.get("IdList", [])
                if not batch:
                    break
                ids.extend(batch)
                retstart += len(batch)
                if len(batch) < retmax:
                    break
            if not ids:
                time.sleep(cfg.sleep_between)
                continue
            # Limitar a max_n ids
            ids = ids[:max_n]
            h2 = Entrez.efetch(db="pubmed", id=ids, rettype="medline", retmode="text")
            parsed = list(Medline.parse(h2))
            h2.close()
            for rec in parsed:
                title = rec.get("TI", "")
                abstract = rec.get("AB", "")
                year = None
                dp = rec.get("DP", "")
                m = re.search(r"(19|20)\d{2}", dp)
                if m:
                    try:
                        year = int(m.group(0))
                    except Exception:
                        year = None
                text = f"{title} {abstract}".lower()
                has_ctx = _detect_any(text, rx_ctx)
                has_cancer = _detect_any(text, rx_cancer)
                if not (has_ctx or has_cancer):
                    continue
                pmid = rec.get("PMID", "")
                link = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}" if pmid else ""
                rows.append({
                    "Gene": gene,
                    "Ensembl_ID": ensembl_id,
                    "Title": title,
                    "Year": year,
                    "Abstract": abstract,
                    "Link": link,
                })
        except Exception as e:
            if logger:
                try:
                    logger.warning(f"PubMed error for gene {gene}: {type(e).__name__}: {e}")
                except Exception:
                    pass
            # be gentle with NCBI
            time.sleep(1.0)
        time.sleep(cfg.sleep_between)

    out = pd.DataFrame(rows)
    if not out.empty:
        out = out.drop_duplicates(subset=["Link"], keep="first")
    return out


def classify_bibliography(
    df: pd.DataFrame,
    cancer_types: Optional[Dict[str, List[str]]] = None,
) -> pd.DataFrame:
    """
    Clasifica artículos por 'cancer_relation', 'cancer_type', 'emt_relation', 'mRNAs_relation'.
    Acepta columnas 'Title' y opcional 'Abstract'.
    """
    if df is None or df.empty:
        return pd.DataFrame()
    cancer_types = cancer_types or DEFAULT_CANCER_TYPES

    def classify_cancer_type(title: str) -> Tuple[bool, Optional[str]]:
        t = (title or "").lower()
        # tipos específicos
        detected = [label for label, kws in cancer_types.items()
                    if any(re.search(rf"\b{re.escape(k.lower())}\b", t) for k in kws)]
        if detected:
            return True, ", ".join(sorted(set(detected)))
        # genérico
        if any(re.search(rf"\b{re.escape(term)}\b", t) for term in GENERAL_CANCER_TERMS):
            return True, "General cancer mention"
        return False, None

    def detect_any(text: str, words: List[str]) -> bool:
        return any(re.search(rf"\b{re.escape(w.lower())}\b", (text or "").lower()) for w in words)

    out = df.copy()
    title_col = "Title" if "Title" in out.columns else ("article_title" if "article_title" in out.columns else None)
    if not title_col:
        raise ValueError("Se requiere columna 'Title' o 'article_title' en el DataFrame de bibliografía.")
    abs_col = "Abstract" if "Abstract" in out.columns else None
    out[["cancer_relation", "cancer_type"]] = out[title_col].apply(lambda x: pd.Series(classify_cancer_type(x)))
    text_src = out[abs_col] if abs_col else out[title_col]
    out["emt_relation"] = text_src.apply(lambda t: detect_any(t, EMT_KEYWORDS))
    out["mRNAs_relation"] = text_src.apply(lambda t: detect_any(t, MICRORNA_KEYWORDS))
    # filtrar solo relacionados a cáncer
    out = out[out["cancer_relation"] == True]
    return out.reset_index(drop=True)


def aggregate_counts_by_level_and_cancer(df: pd.DataFrame) -> pd.DataFrame:
    if df is None or df.empty:
        return pd.DataFrame()
    if "nivel_expresion" not in df.columns or "cancer_type" not in df.columns:
        return pd.DataFrame()
    t = df.copy()
    t["cancer_type"] = t["cancer_type"].fillna("")
    # expandir múltiples tipos en filas
    t = t.assign(cancer_type=t["cancer_type"].astype(str).str.split(", ")).explode("cancer_type")
    g = t.groupby(["nivel_expresion", "cancer_type"], as_index=False).size().rename(columns={"size": "count"})
    return g


# ---------------- Heurística de interpretación por cáncer -----------------

def _norm(s: str) -> str:
    return (s or "").strip().lower()


def _best_fuzzy_key(query: str, keys: List[str], cutoff: float = 0.82) -> Optional[str]:
    matches = difflib.get_close_matches(query, keys, n=1, cutoff=cutoff)
    return matches[0] if matches else None


def _match_cancer_key(label: str) -> Optional[str]:
    """Mejor coincidencia: exacta, por sinónimos/inclusión, luego fuzzy por keys."""
    if not label:
        return None
    lab = _norm(label)
    keys = list(DEFAULT_CANCER_TYPES.keys())
    # exacta por key
    for k in keys:
        if lab == k.lower():
            return k
    # sinónimos / inclusión
    for k in keys:
        syns = [k] + DEFAULT_CANCER_TYPES.get(k, [])
        if any(lab == _norm(s) for s in syns):
            return k
        if any(lab in _norm(s) or _norm(s) in lab for s in syns):
            return k
    # fuzzy por keys
    fk = _best_fuzzy_key(lab, [k.lower() for k in keys], cutoff=0.82)
    if fk:
        return next((k for k in keys if k.lower() == fk), None)
    return None


HEUR_UP = [
    "overexpress", "over-expression", "overexpression", "upregulated", "up-regulated",
    "high expression", "elevated", "increased", "up regulation", "gain of", "amplified",
    "induces expression", "activates", "enhances", "promotes expression",
]
HEUR_DOWN = [
    "downregulated", "down-regulated", "low expression", "decreased", "reduced",
    "silenced", "loss of", "suppressed", "knockdown reduces", "knockout reduces",
    "inhibits expression", "attenuates",
]
HEUR_PROGNOSIS_BAD = [
    "poor prognosis", "worse survival", "shorter survival", "adverse outcome", "high risk",
    "lower overall survival", "lower os", "lower dfs", "higher hazard", "hazard ratio >",
]
HEUR_PROGNOSIS_GOOD = [
    "better survival", "good prognosis", "favorable prognosis", "longer survival",
    "higher overall survival", "higher os", "beneficial prognosis",
]
HEUR_FUNCTIONS = {
    "proliferation": ["proliferation", "cell growth", "cell cycle progression", "s phase entry"],
    "apoptosis": ["apoptosis", "apoptotic", "anoikis", "caspase activation"],
    "invasion": ["invasion", "invasive", "matrigel invasion"],
    "migration": ["migration", "migratory", "wound healing assay"],
    "metastasis": ["metastasis", "metastatic", "metastatic potential"],
    "drug_resistance": ["drug resistance", "chemoresistance", "resistant", "multidrug resistance"],
    "drug_sensitivity": ["chemosensitivity", "drug sensitivity", "sensitizes", "resensitizes"],
    "emt_related": ["emt", "epithelial-mesenchymal transition", "vimentin", "e-cadherin loss"],
}
NEGATIONS = [
    "not", "no", "lack of", "did not", "does not", "fails to", "without", "little to no",
    "non-", "absence of", "insufficient", "no association", "not correlated",
]
HEDGES = [
    "may", "might", "could", "suggests", "appears to", "potentially", "trend towards", "likely",
]


def _contains_any(text: str, words: List[str]) -> bool:
    rx = _compile_kw_regex(words)
    return bool(rx.search(text or ""))


def filter_bibliography_by_cancer(df: pd.DataFrame, cancer_label: str) -> pd.DataFrame:
    """Filtra artículos cuyo título/abstract mencionan el cáncer indicado o cuyo
    campo 'cancer_type' lo contiene (cuando está disponible)."""
    if df is None or df.empty:
        return pd.DataFrame()
    k = _match_cancer_key(cancer_label)
    terms = DEFAULT_CANCER_TYPES.get(k, GENERAL_CANCER_TERMS)
    t = df.copy()
    title_col = "Title" if "Title" in t.columns else ("article_title" if "article_title" in t.columns else None)
    abs_col = "Abstract" if "Abstract" in t.columns else None
    title_s = (t[title_col] if title_col else pd.Series("", index=t.index)).astype(str)
    abs_s = (t[abs_col] if abs_col else pd.Series("", index=t.index)).astype(str)
    text = (title_s + " " + abs_s)
    mask_terms = text.apply(lambda s: bool(_compile_kw_regex(terms).search(s)))
    mask_ct = (
        t.get("cancer_type", pd.Series("", index=t.index)).astype(str)
         .str.contains(str(k or cancer_label), case=False, na=False)
    )
    out = t[mask_terms | mask_ct]
    return out.reset_index(drop=True)


def _span_contains_negation(txt: str, start: int, end: int) -> bool:
    left = max(0, start - 200)
    right = min(len(txt), end + 200)
    return bool(_compile_kw_regex(NEGATIONS).search(txt[left:right]))


def _span_contains_hedge(txt: str, start: int, end: int) -> bool:
    left = max(0, start - 200)
    right = min(len(txt), end + 200)
    return bool(_compile_kw_regex(HEDGES).search(txt[left:right]))


def _score_hits(txt: str, gene: str, lexicon: List[str], title_txt: str = "") -> float:
    if not txt:
        return 0.0
    rx = _compile_kw_regex(lexicon)
    score = 0.0
    for m in rx.finditer(txt):
        s = 1.0
        if gene:
            try:
                gpos = [gm.start() for gm in re.finditer(re.escape(gene), txt, flags=re.IGNORECASE)]
            except re.error:
                gpos = []
            if gpos and any(abs(g - m.start()) <= 80 for g in gpos):
                s *= 1.3
        if _span_contains_negation(txt, m.start(), m.end()):
            continue
        if _span_contains_hedge(txt, m.start(), m.end()):
            s *= 0.7
        score += s
    if title_txt and rx.search(title_txt):
        score *= 1.5
    return score


def interpret_gene_relations(
    df: pd.DataFrame,
    gene_col: str = "Gene",
    title_col: str = "Title",
    abs_col: str = "Abstract",
) -> pd.DataFrame:
    """Calcula scores heurísticos por relación (+ flags por umbral) y efecto neto de expresión."""
    if df is None or df.empty:
        return pd.DataFrame()
    t = df.copy()
    if title_col not in t.columns and "article_title" in t.columns:
        title_col = "article_title"
    if abs_col not in t.columns:
        abs_col = None

    def _txt(row):
        title = str(row.get(title_col, "") or "")
        abstract = str(row.get(abs_col, "") or "") if abs_col else ""
        return title, (title + " " + abstract)

    def row_scores(row):
        gene = str(row.get(gene_col, "") or "").strip()
        title, full = _txt(row)
        out: Dict[str, float | bool] = {}
        out["upregulated_score"] = _score_hits(full, gene, HEUR_UP, title)
        out["downregulated_score"] = _score_hits(full, gene, HEUR_DOWN, title)
        out["prognosis_bad_score"] = _score_hits(full, gene, HEUR_PROGNOSIS_BAD, title)
        out["prognosis_good_score"] = _score_hits(full, gene, HEUR_PROGNOSIS_GOOD, title)
        for k, words in HEUR_FUNCTIONS.items():
            out[f"{k}_score"] = _score_hits(full, gene, words, title)
        TH = 1.2
        out["upregulated"] = out["upregulated_score"] >= TH
        out["downregulated"] = out["downregulated_score"] >= TH
        out["prognosis_bad"] = out["prognosis_bad_score"] >= TH
        out["prognosis_good"] = out["prognosis_good_score"] >= TH
        for k in HEUR_FUNCTIONS.keys():
            out[k] = out[f"{k}_score"] >= TH
        return out

    scores_df = t.apply(row_scores, axis=1, result_type="expand")
    out = pd.concat([t.reset_index(drop=True), scores_df], axis=1)

    def net_direction(r) -> str:
        up = float(r.get("upregulated_score", 0.0))
        dn = float(r.get("downregulated_score", 0.0))
        if max(up, dn) < 1.2:
            return "uncertain"
        if up > dn * 1.25:
            return "up"
        if dn > up * 1.25:
            return "down"
        return "mixed"

    out["expression_effect"] = out.apply(net_direction, axis=1)
    return out


def summarize_relations_by_gene(df: pd.DataFrame, gene_col: str = "Gene") -> pd.DataFrame:
    if df is None or df.empty or gene_col not in df.columns:
        return pd.DataFrame()
    rel_flags = [
        "upregulated", "downregulated", "prognosis_bad", "prognosis_good",
        "proliferation", "apoptosis", "invasion", "migration", "metastasis",
        "drug_resistance", "drug_sensitivity", "emt_related",
    ]
    rel_scores = [f"{c}_score" for c in rel_flags]
    present_flags = [c for c in rel_flags if c in df.columns]
    present_scores = [c for c in rel_scores if c in df.columns]
    agg_dict: Dict[str, str] = {c: "sum" for c in present_flags}
    agg_dict.update({c: "sum" for c in present_scores})
    g = df.groupby(gene_col).agg(agg_dict).reset_index()

    def _mk_summary(r) -> str:
        parts = []
        for c in present_flags:
            try:
                val = int(r.get(c, 0))
                if val > 0:
                    parts.append(f"{c}({val})")
            except Exception:
                pass
        return ", ".join(parts)

    g["heuristic_summary"] = g.apply(_mk_summary, axis=1)

    def _net_gene(r) -> str:
        up = float(r.get("upregulated_score", 0.0))
        dn = float(r.get("downregulated_score", 0.0))
        if max(up, dn) < 2.5:
            return "uncertain"
        if up > dn * 1.2:
            return "up"
        if dn > up * 1.2:
            return "down"
        return "mixed"

    g["net_expression_effect"] = g.apply(_net_gene, axis=1)
    return g
