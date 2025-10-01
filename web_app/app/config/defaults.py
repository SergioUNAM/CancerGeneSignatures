from __future__ import annotations

from typing import Any, Dict

# Valores por defecto que se usaban históricamente en la app.
DEFAULT_PROJECT_NAME = "CancerGeneSignatures - Web App"
DEFAULT_LOG_LEVEL = "INFO"

DEFAULT_MENU_DATA: Dict[str, Any] = {
    "version": 1,
    "menu": {
        "cancer_types": [
            "Breast Cancer",
            "Melanoma",
            "Colon Cancer",
        ],
        "contexts": [
            {"key": "TEM", "label": "Cáncer y TEM"},
            {"key": "micro_rnas", "label": "Cáncer y micro RNAs"},
        ],
        "normalization_methods": [
            {"key": "reference_gene", "label": "gen de referencia"},
            {"key": "means", "label": "promedios"},
        ],
    },
}

DEFAULT_SERVICES_SETTINGS: Dict[str, Any] = {
    "ensembl": {
        "max_workers": 3,
    },
    "string": {
        "base_url": None,
        "species": 9606,
    },
    "pubmed": {
        "default_email": None,
        "api_key": None,
        "max_per_gene": 20,
    },
    "google_nlp": {
        "api_key": None,
    },
}
