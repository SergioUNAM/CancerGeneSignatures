from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import streamlit as st

from src.core.io import LoadResult


@dataclass
class UndeterminedSettings:
    """Preferencias sobre el tratamiento de valores 'Undetermined/ND'."""

    policy: str = "nan"
    value: float = 40.0


@dataclass
class AppSessionState:
    """Contenedor tipado para valores persistidos en ``st.session_state``."""

    undetermined: UndeterminedSettings = field(default_factory=UndeterminedSettings)
    context_label: str = ""
    exclude_stable: bool = False
    df_loaded: Optional[LoadResult] = None

    @classmethod
    def load(cls) -> "AppSessionState":
        """Recupera el estado almacenado o crea uno nuevo."""
        existing = st.session_state.get("_app_state")
        if isinstance(existing, cls):
            return existing
        state = cls()
        st.session_state["_app_state"] = state
        return state

    def persist(self) -> None:
        """Guarda el estado actual y sincroniza claves legadas."""
        st.session_state["_app_state"] = self
        st.session_state["und_policy"] = self.undetermined.policy
        st.session_state["und_value"] = self.undetermined.value
        st.session_state["exclude_stable"] = self.exclude_stable
        st.session_state["df_loaded"] = self.df_loaded
