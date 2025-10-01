from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Mapping, Optional

import pandas as pd

from src.core.string_enrichment import enrich_by_levels, filter_enrichment


@dataclass(frozen=True)
class EnrichmentResult:
    """Estructura con los resultados filtrados del enriquecimiento STRING."""

    combined: pd.DataFrame
    by_level: Dict[str, pd.DataFrame]

    def is_empty(self) -> bool:
        return (self.combined is None or self.combined.empty) and all(
            df.empty for df in self.by_level.values()
        )


def perform_string_enrichment(
    expression_df: pd.DataFrame,
    *,
    levels: Iterable[str],
    species: int,
    sources: Optional[Iterable[str]] = None,
    max_fdr: float = 0.05,
    min_term_genes: int = 3,
    top_n: int = 25,
) -> EnrichmentResult:
    """Ejecuta y filtra el enriquecimiento STRING por nivel de expresión."""

    levels_list = list(levels)
    sources_list = list(sources) if sources is not None else None

    res = enrich_by_levels(
        expression_df,
        symbol_col="target",
        level_col="nivel_expresion",
        levels=levels_list,
        species=int(species),
        sources=sources_list,
    )
    combined = res.get("combined") if isinstance(res, Mapping) else None
    by_level_raw = res.get("by_level") if isinstance(res, Mapping) else {}

    filtered_combined = filter_enrichment(
        combined,
        include_categories=sources_list,
        max_fdr=float(max_fdr),
        min_term_genes=int(min_term_genes),
        top_n=int(top_n),
    ) if combined is not None else pd.DataFrame()

    filtered_by_level: Dict[str, pd.DataFrame] = {}
    if isinstance(by_level_raw, Mapping):
        for lvl, df_lvl in by_level_raw.items():
            filtered_by_level[str(lvl)] = filter_enrichment(
                df_lvl,
                include_categories=sources_list,
                max_fdr=float(max_fdr),
                min_term_genes=int(min_term_genes),
                top_n=int(top_n),
            )

    return EnrichmentResult(
        combined=filtered_combined,
        by_level=filtered_by_level,
    )


__all__ = [
    "EnrichmentResult",
    "perform_string_enrichment",
]
