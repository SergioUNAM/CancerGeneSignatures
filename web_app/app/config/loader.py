from __future__ import annotations

import json
import logging
import os
from pathlib import Path
from typing import Any, Callable, Mapping, Optional

from .defaults import (
    DEFAULT_LOG_LEVEL,
    DEFAULT_MENU_DATA,
    DEFAULT_PROJECT_NAME,
    DEFAULT_SERVICES_SETTINGS,
)
from .models import (
    AppConfig,
    ConfigValidationError,
    GoogleNLPConfig,
    MenuConfig,
    PubMedConfig,
    ServicesConfig,
    StringConfig,
    EnsemblConfig,
)

LOGGER = logging.getLogger(__name__)


class ConfigError(Exception):
    """Error genérico al cargar la configuración de la aplicación."""


def load_app_config(
    menu_path: Optional[Path] = None,
    *,
    env: Optional[Mapping[str, str]] = None,
    secrets: Optional[Mapping[str, Any]] = None,
    warn: Optional[Callable[[str], None]] = None,
) -> AppConfig:
    """Carga la configuración de la app consolidando defaults, archivo y entorno."""

    env_map = env or os.environ
    warn_fn = warn or (lambda msg: LOGGER.warning(msg))

    project_name = env_map.get("CGS_PROJECT_NAME", DEFAULT_PROJECT_NAME)
    log_level = env_map.get("CGS_LOGLEVEL", DEFAULT_LOG_LEVEL).upper()

    menu_config = _load_menu_config(menu_path=menu_path, env=env_map, secrets=secrets, warn=warn_fn)
    services_config = _load_services_config(env=env_map, secrets=secrets, warn=warn_fn)

    app_config = AppConfig(
        project_name=project_name,
        log_level=log_level,
        menu=menu_config,
        services=services_config,
    )
    return app_config.sanitized()


def _load_menu_config(
    *,
    menu_path: Optional[Path],
    env: Mapping[str, str],
    secrets: Optional[Mapping[str, Any]],
    warn: Callable[[str], None],
) -> MenuConfig:
    menu_data: Optional[Mapping[str, Any]] = None

    # Prioridad: argumento explícito > env > secrets > archivo por defecto > defaults en código
    if menu_path and menu_path.exists():
        menu_data = _read_json(menu_path, warn)
    else:
        env_path = env.get("CGS_MENU_PATH")
        if env_path:
            path = Path(env_path)
            if path.exists():
                menu_data = _read_json(path, warn)
            else:
                warn(f"No se encontró el archivo definido en CGS_MENU_PATH: {path}")

    if menu_data is None and secrets is not None:
        menu_secret = secrets.get("menu")
        if isinstance(menu_secret, Mapping):
            menu_data = menu_secret

    if menu_data is None:
        default_menu_path = Path(__file__).resolve().parents[2] / "config" / "menu.json"
        if default_menu_path.exists():
            menu_data = _read_json(default_menu_path, warn)
        else:
            warn("No se encontró config/menu.json; usando menú por defecto embebido.")
            menu_data = DEFAULT_MENU_DATA

    try:
        return MenuConfig.from_dict(menu_data)
    except ConfigValidationError as exc:
        raise ConfigError(f"Menú inválido: {exc}") from exc


def _load_services_config(
    *,
    env: Mapping[str, str],
    secrets: Optional[Mapping[str, Any]],
    warn: Callable[[str], None],
) -> ServicesConfig:
    defaults = DEFAULT_SERVICES_SETTINGS

    ensembl_max_workers = _read_int(env.get("CGS_ENSEMBL_MAX_WORKERS"), defaults["ensembl"]["max_workers"], warn)

    string_base_url = env.get("CGS_STRING_BASE_URL") or defaults["string"]["base_url"]
    string_species = _read_int(env.get("CGS_STRING_SPECIES"), defaults["string"]["species"], warn)

    pubmed_email = (
        env.get("CGS_PUBMED_EMAIL")
        or env.get("NCBI_EMAIL")
        or defaults["pubmed"].get("default_email")
    )
    pubmed_api_key = _get_secret(
        "pubmed_api_key",
        env_key="CGS_PUBMED_API_KEY",
        env=env,
        secrets=secrets,
        alt_env_keys=["NCBI_API_KEY"],
        alt_secret_keys=["NCBI_API_KEY"],
    )
    pubmed_max_per_gene = _read_int(env.get("CGS_PUBMED_MAX_PER_GENE"), defaults["pubmed"]["max_per_gene"], warn)

    google_api_key = _get_secret(
        "google_nlp_api_key",
        env_key="CGS_GOOGLE_NLP_API_KEY",
        env=env,
        secrets=secrets,
        alt_env_keys=["GOOGLE_NLP_API_KEY"],
        alt_secret_keys=["GOOGLE_NLP_API_KEY"],
    )

    ensembl_cfg = EnsemblConfig(max_workers=ensembl_max_workers)
    string_cfg = StringConfig(base_url=string_base_url, species=string_species)
    pubmed_cfg = PubMedConfig(default_email=pubmed_email, api_key=pubmed_api_key, max_per_gene=pubmed_max_per_gene)
    google_cfg = GoogleNLPConfig(api_key=google_api_key)

    return ServicesConfig(
        ensembl=ensembl_cfg,
        string=string_cfg,
        pubmed=pubmed_cfg,
        google_nlp=google_cfg,
    ).sanitized()


def _read_json(path: Path, warn: Callable[[str], None]) -> Mapping[str, Any]:
    try:
        with path.open("r", encoding="utf-8") as handler:
            return json.load(handler)
    except FileNotFoundError:
        warn(f"Archivo de configuración no encontrado: {path}")
    except json.JSONDecodeError as exc:
        raise ConfigError(f"No se pudo parsear {path}: {exc}") from exc
    return DEFAULT_MENU_DATA


def _read_int(value: Optional[str], default: int, warn: Callable[[str], None]) -> int:
    if value is None:
        return int(default)
    try:
        return int(value)
    except (TypeError, ValueError):
        warn(f"Valor entero inválido '{value}', usando {default}.")
        return int(default)


def _get_secret(
    key: str,
    *,
    env_key: str,
    env: Mapping[str, str],
    secrets: Optional[Mapping[str, Any]],
) -> Optional[str]:
    if env_key in env and env[env_key]:
        return env[env_key]
    alt_env_keys = []
    alt_secret_keys = []
    if key == "pubmed_api_key":
        alt_env_keys.append("NCBI_API_KEY")
        alt_secret_keys.append("NCBI_API_KEY")
    if key == "google_nlp_api_key":
        alt_env_keys.append("GOOGLE_NLP_API_KEY")
        alt_secret_keys.append("GOOGLE_NLP_API_KEY")
    for alt in alt_env_keys:
        if alt in env and env[alt]:
            return env[alt]
    if secrets is not None and key in secrets and secrets[key]:
        value = secrets[key]
        return str(value)
    if secrets is not None:
        for alt in alt_secret_keys:
            if alt in secrets and secrets[alt]:
                return str(secrets[alt])
    return None


__all__ = [
    "load_app_config",
    "ConfigError",
]
