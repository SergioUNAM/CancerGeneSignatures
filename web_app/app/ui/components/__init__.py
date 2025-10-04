"""Componentes UI generales reutilizables."""

from .exports import ExportRegistry, render_export_panel  # noqa: F401
from . import highlights as _highlights
from .insight_board import Badge, Card, render_insight_board  # noqa: F401
from .insights_overview import (  # noqa: F401
    BoardSpec,
    InsightsData,
    build_extraction_specs,
    build_ddct_specs,
    build_ddct_comparative_spec,
    build_fc_specs,
    build_fc_comparative_spec,
    build_concordance_spec_from_summary,
    build_default_insight_specs,
    render_insights_overview,
)
from .classification_preview import render_classification_preview  # noqa: F401
from .gene_pills import render_gene_pills_fc, render_gene_methods_view  # noqa: F401
from .cards import InfoCard, render_info_cards  # noqa: F401
from .progress import (  # noqa: F401
    build_step_sequence,
    render_pipeline_progress,
    render_sidebar_progress,
)

Highlight = _highlights.Highlight  # noqa: F401
render_highlight_pills = _highlights.render_highlight_pills  # noqa: F401

__all__ = (
    "build_step_sequence",
    "render_pipeline_progress",
    "render_sidebar_progress",
    "ExportRegistry",
    "render_export_panel",
    "Highlight",
    "render_highlight_pills",
    "Badge",
    "Card",
    "render_insight_board",
    "render_classification_preview",
    "render_gene_pills_fc",
    "render_gene_methods_view",
    "BoardSpec",
    "InsightsData",
    "build_extraction_specs",
    "build_ddct_specs",
    "build_ddct_comparative_spec",
    "build_fc_specs",
    "build_fc_comparative_spec",
    "build_concordance_spec_from_summary",
    "build_default_insight_specs",
    "render_insights_overview",
    "InfoCard",
    "render_info_cards",
)
