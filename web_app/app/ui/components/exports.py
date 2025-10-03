"""Helpers para consolidar descargas en un panel de exportación."""

from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass
from datetime import datetime
from typing import Any, Iterable, Mapping, MutableMapping, Optional

import pandas as pd
import streamlit as st


@dataclass
class ExportArtifact:
    """Representa un artefacto descargable generado durante la sesión."""

    key: str
    section: str
    label: str
    file_name: str
    payload: bytes | str
    method_label: Optional[str] = None
    description: Optional[str] = None
    parameters: Optional[Mapping[str, Any]] = None
    mime: str = "text/csv"

    def as_payload(self) -> bytes | str:
        """Retorna la carga normalizada para `st.download_button`."""

        return self.payload


class ExportRegistry:
    """Acumula los artefactos generados para mostrarlos en un panel unificado."""

    def __init__(self) -> None:
        self._artifacts: MutableMapping[str, ExportArtifact] = {}
        self._order: list[str] = []

    def register(self, artifact: ExportArtifact) -> None:
        """Registra o actualiza un artefacto en la colección."""

        if artifact.key not in self._artifacts:
            self._order.append(artifact.key)
        self._artifacts[artifact.key] = artifact

    def register_dataframe(
        self,
        *,
        key: str,
        section: str,
        label: str,
        file_name: str,
        dataframe: pd.DataFrame,
        method_label: Optional[str] = None,
        description: Optional[str] = None,
        parameters: Optional[Mapping[str, Any]] = None,
        include_index: bool = False,
    ) -> None:
        """Convierte un DataFrame a CSV y lo agrega al registro."""

        csv_data = dataframe.to_csv(index=include_index)
        artifact = ExportArtifact(
            key=key,
            section=section,
            label=label,
            file_name=file_name,
            payload=csv_data,
            method_label=method_label,
            description=description,
            parameters=parameters,
            mime="text/csv",
        )
        self.register(artifact)

    def register_text(
        self,
        *,
        key: str,
        section: str,
        label: str,
        file_name: str,
        text: str,
        method_label: Optional[str] = None,
        description: Optional[str] = None,
        parameters: Optional[Mapping[str, Any]] = None,
        mime: str = "text/plain",
    ) -> None:
        """Registra contenido plano en el panel de exportación."""

        artifact = ExportArtifact(
            key=key,
            section=section,
            label=label,
            file_name=file_name,
            payload=text,
            method_label=method_label,
            description=description,
            parameters=parameters,
            mime=mime,
        )
        self.register(artifact)

    def __iter__(self) -> Iterable[ExportArtifact]:
        for key in self._order:
            artifact = self._artifacts.get(key)
            if artifact is not None:
                yield artifact

    def group_by_section(self) -> Mapping[str, list[ExportArtifact]]:
        grouped: "OrderedDict[str, list[ExportArtifact]]" = OrderedDict()
        for artifact in self:
            grouped.setdefault(artifact.section, []).append(artifact)
        return grouped

    def is_empty(self) -> bool:
        return not self._order


def render_export_panel(
    registry: ExportRegistry,
    *,
    study_context: Optional[Mapping[str, Any]] = None,
    title: str = "Panel de exportación consolidado",
) -> None:
    """Renderiza un panel con enlaces de descarga consolidados."""

    if registry.is_empty():
        return

    st.divider()
    st.subheader(title)
    st.caption(
        "Descarga desde un solo lugar los artefactos generados en esta sesión. "
        "Cada entrada incluye el método o presets activos para favorecer la trazabilidad."
    )

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")
    st.caption(f"Generado el {timestamp} (hora local del servidor)")

    if study_context:
        context_lines = [
            f"• **{k}**: {v}" for k, v in study_context.items() if v not in (None, "")
        ]
        if context_lines:
            st.markdown("**Contexto del análisis**")
            for line in context_lines:
                st.markdown(line)

    grouped = registry.group_by_section()
    for section, artifacts in grouped.items():
        with st.expander(f"{section} · {len(artifacts)} descargas", expanded=True):
            for artifact in artifacts:
                header = artifact.label
                if artifact.method_label:
                    header = f"{header} — {artifact.method_label}"
                st.markdown(f"**{header}**")
                if artifact.description:
                    st.caption(artifact.description)
                if artifact.parameters:
                    params = ", ".join(
                        f"{k}: {v}" for k, v in artifact.parameters.items() if v not in (None, "")
                    )
                    if params:
                        st.caption(f"Parámetros clave → {params}")
                st.download_button(
                    f"Descargar {artifact.label}",
                    artifact.as_payload(),
                    file_name=artifact.file_name,
                    mime=artifact.mime,
                    key=f"export::{artifact.key}",
                    use_container_width=True,
                )

    st.caption(
        "Tip: conserva este panel como evidencia de configuración adjuntándolo al cuaderno de análisis."
    )
