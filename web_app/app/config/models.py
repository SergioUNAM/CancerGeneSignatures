from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Mapping, Optional, Sequence


class ConfigValidationError(Exception):
    """Se produce cuando un archivo de configuración presenta valores inválidos."""


@dataclass(frozen=True)
class MenuOption:
    key: str
    label: str

    def __post_init__(self) -> None:
        if not self.key:
            raise ConfigValidationError("Los elementos de menú requieren un 'key'.")
        if not self.label:
            raise ConfigValidationError("Los elementos de menú requieren un 'label'.")


NormalizationMethod = MenuOption


@dataclass(frozen=True)
class MenuConfig:
    version: int
    cancer_types: List[str] = field(default_factory=list)
    contexts: List[MenuOption] = field(default_factory=list)
    normalization_methods: List[NormalizationMethod] = field(default_factory=list)
    created_at: Optional[str] = None

    @classmethod
    def from_dict(cls, data: Mapping[str, object]) -> "MenuConfig":
        menu_section = data.get("menu") if isinstance(data.get("menu"), Mapping) else {}
        cancer_types_raw = list(menu_section.get("cancer_types", [])) if isinstance(menu_section, Mapping) else []
        contexts_raw = list(menu_section.get("contexts", [])) if isinstance(menu_section, Mapping) else []
        norm_methods_raw = list(menu_section.get("normalization_methods", [])) if isinstance(menu_section, Mapping) else []

        contexts = [MenuOption(**ctx) for ctx in contexts_raw if isinstance(ctx, Mapping)]
        normalization_methods = [NormalizationMethod(**nm) for nm in norm_methods_raw if isinstance(nm, Mapping)]

        config = cls(
            version=int(data.get("version", 1)),
            cancer_types=[str(item) for item in cancer_types_raw if item],
            contexts=contexts,
            normalization_methods=normalization_methods,
            created_at=str(data.get("created_at")) if data.get("created_at") else None,
        )
        config.validate()
        return config

    def validate(self) -> None:
        if not self.cancer_types:
            raise ConfigValidationError("La configuración de menú requiere al menos un tipo de cáncer.")
        if not self.contexts:
            raise ConfigValidationError("La configuración de menú requiere al menos un contexto disponible.")
        if not self.normalization_methods:
            raise ConfigValidationError("La configuración de menú requiere métodos de normalización.")
        if len(set(self.cancer_types)) != len(self.cancer_types):
            raise ConfigValidationError("Los tipos de cáncer deben ser únicos.")
        _ensure_unique_labels(self.contexts, "contextos")
        _ensure_unique_labels(self.normalization_methods, "métodos de normalización")


@dataclass(frozen=True)
class EnsemblConfig:
    max_workers: int = 3

    def sanitized(self) -> "EnsemblConfig":
        workers = max(1, int(self.max_workers))
        return EnsemblConfig(max_workers=workers)


@dataclass(frozen=True)
class StringConfig:
    base_url: Optional[str] = None
    species: int = 9606

    def sanitized(self) -> "StringConfig":
        species = int(self.species) if self.species else 9606
        return StringConfig(base_url=self.base_url or None, species=species)


@dataclass(frozen=True)
class PubMedConfig:
    default_email: Optional[str] = None
    api_key: Optional[str] = None
    max_per_gene: int = 20

    def sanitized(self) -> "PubMedConfig":
        max_per_gene = max(1, int(self.max_per_gene))
        return PubMedConfig(
            default_email=self.default_email or None,
            api_key=self.api_key or None,
            max_per_gene=max_per_gene,
        )


@dataclass(frozen=True)
class GoogleNLPConfig:
    api_key: Optional[str] = None

    def sanitized(self) -> "GoogleNLPConfig":
        return GoogleNLPConfig(api_key=self.api_key or None)


@dataclass(frozen=True)
class ServicesConfig:
    ensembl: EnsemblConfig = field(default_factory=EnsemblConfig)
    string: StringConfig = field(default_factory=StringConfig)
    pubmed: PubMedConfig = field(default_factory=PubMedConfig)
    google_nlp: GoogleNLPConfig = field(default_factory=GoogleNLPConfig)

    def sanitized(self) -> "ServicesConfig":
        return ServicesConfig(
            ensembl=self.ensembl.sanitized(),
            string=self.string.sanitized(),
            pubmed=self.pubmed.sanitized(),
            google_nlp=self.google_nlp.sanitized(),
        )


@dataclass(frozen=True)
class AppConfig:
    project_name: str
    log_level: str
    menu: MenuConfig
    services: ServicesConfig = field(default_factory=ServicesConfig)

    def sanitized(self) -> "AppConfig":
        return AppConfig(
            project_name=self.project_name,
            log_level=self.log_level.upper(),
            menu=self.menu,
            services=self.services.sanitized(),
        )


def _ensure_unique_labels(options: Sequence[MenuOption], scope: str) -> None:
    labels = [opt.label for opt in options]
    if len(labels) != len(set(labels)):
        raise ConfigValidationError(f"Los {scope} deben tener etiquetas únicas.")
