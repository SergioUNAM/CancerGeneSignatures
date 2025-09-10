from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List, Dict, Tuple, Optional, Callable, Any
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


def _require_entrez_config() -> PubMedConfig:
    email = os.getenv("NCBI_EMAIL")
    if not email:
        raise RuntimeError("NCBI_EMAIL no configurado en el entorno.")
    api_key = os.getenv("NCBI_API_KEY")
    return PubMedConfig(email=email, api_key=api_key)


def _setup_entrez(cfg: PubMedConfig) -> None:
    if Entrez is None:
        raise RuntimeError("Biopython no está instalado. Agrega 'biopython' a requirements.")
    Entrez.email = cfg.email
    if cfg.api_key:
        Entrez.api_key = cfg.api_key


def _compile_kw_regex(words: List[str]) -> re.Pattern:
    words_escaped = [re.escape(w.lower()) for w in words if w]
    return re.compile(r"\b(" + "|".join(words_escaped) + r")\b") if words_escaped else re.compile(r"$.^")


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
) -> pd.DataFrame:
    """
    Busca artículos en PubMed por gen (símbolo) y filtra por contexto/cáncer de forma ligera.
    Devuelve columnas: Gene, Ensembl_ID, Title, Year, Abstract, Link
    """
    cfg = _require_entrez_config()
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
            h = Entrez.esearch(db="pubmed", term=query, retmax=max_n)
            record = Entrez.read(h)
            h.close()
            ids = record.get("IdList", [])
            if not ids:
                time.sleep(cfg.sleep_between)
                continue
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
