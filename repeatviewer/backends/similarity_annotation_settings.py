"""User-configurable parameters for similarity-based annotation."""

from __future__ import annotations

from dataclasses import dataclass

from PySide6.QtCore import QSettings


DEFAULT_SIMILARITY_EVALUE = 1e-4
SETTINGS_VERSION = 8


@dataclass
class SimilarityAnnotationSettings:
    """Search thresholds for similarity annotation."""

    blastn_evalue: float = DEFAULT_SIMILARITY_EVALUE
    blastp_evalue: float = DEFAULT_SIMILARITY_EVALUE
    hmm_evalue: float = DEFAULT_SIMILARITY_EVALUE
    cdd_evalue: float = DEFAULT_SIMILARITY_EVALUE
    max_hits_per_database: int = 10
    show_database_name_on_track: bool = True
    show_process_window: bool = True
    close_process_window_when_finished: bool = True


def format_similarity_evalue_suffix(evalue: float | str | None) -> str:
    if evalue is None:
        return ""
    text = str(evalue).strip()
    if not text:
        return ""
    try:
        value = float(text)
    except (TypeError, ValueError):
        return ""
    if value < 0:
        return ""
    if value >= 0.01:
        return f" (E={value:g})"
    return f" (E={value:.2g})"


def similarity_feature_name(
    database_name: str,
    target_name: str,
    *,
    show_database_name: bool,
    evalue: float | str | None = None,
) -> str:
    target = target_name.strip() or "hit"
    if show_database_name and database_name.strip():
        label = f"{database_name}: {target}"
    else:
        label = target
    suffix = format_similarity_evalue_suffix(evalue)
    return f"{label}{suffix}" if suffix else label


def similarity_feature_track_label(
    name: str,
    feature_type: str,
    *,
    show_database_name: bool,
) -> str:
    if not name:
        return feature_type
    if show_database_name or feature_type != "domain":
        return name
    _prefix, sep, rest = name.partition(": ")
    if sep and rest:
        return rest
    return name


def default_similarity_annotation_settings() -> SimilarityAnnotationSettings:
    return SimilarityAnnotationSettings()


def _get_float(s: QSettings, key: str, default: float) -> float:
    return float(s.value(key, default))


def _get_int(s: QSettings, key: str, default: int) -> int:
    return int(s.value(key, default))


def _get_bool(s: QSettings, key: str, default: bool) -> bool:
    value = s.value(key, default)
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.lower() in ("1", "true", "yes")
    return bool(value)


def load_similarity_annotation_settings() -> SimilarityAnnotationSettings:
    s = QSettings()
    prefix = "similarity_annotation/"
    version = _get_int(s, f"{prefix}settings_version", 0)
    max_hits = _get_int(s, f"{prefix}max_hits_per_database", 10)
    show_db_name = _get_bool(s, f"{prefix}show_database_name_on_track", True)
    show_process = _get_bool(s, f"{prefix}show_process_window", True)
    close_process = _get_bool(s, f"{prefix}close_process_window_when_finished", True)
    # v6 and earlier stored a single blast_evalue for both BLASTN and BLASTP.
    legacy_blast = _get_float(s, f"{prefix}blast_evalue", DEFAULT_SIMILARITY_EVALUE)
    blastn = _get_float(s, f"{prefix}blastn_evalue", legacy_blast)
    blastp = _get_float(s, f"{prefix}blastp_evalue", legacy_blast)
    hmm = _get_float(s, f"{prefix}hmm_evalue", DEFAULT_SIMILARITY_EVALUE)
    cdd = _get_float(s, f"{prefix}cdd_evalue", DEFAULT_SIMILARITY_EVALUE)
    settings = SimilarityAnnotationSettings(
        blastn_evalue=blastn,
        blastp_evalue=blastp,
        hmm_evalue=hmm,
        cdd_evalue=cdd,
        max_hits_per_database=max_hits,
        show_database_name_on_track=show_db_name,
        show_process_window=show_process,
        close_process_window_when_finished=close_process,
    )
    if version < SETTINGS_VERSION:
        save_similarity_annotation_settings(settings)
    return settings


def save_similarity_annotation_settings(settings: SimilarityAnnotationSettings) -> None:
    s = QSettings()
    prefix = "similarity_annotation/"
    s.setValue(f"{prefix}settings_version", SETTINGS_VERSION)
    s.setValue(f"{prefix}blastn_evalue", settings.blastn_evalue)
    s.setValue(f"{prefix}blastp_evalue", settings.blastp_evalue)
    s.setValue(f"{prefix}hmm_evalue", settings.hmm_evalue)
    s.setValue(f"{prefix}cdd_evalue", settings.cdd_evalue)
    s.setValue(f"{prefix}max_hits_per_database", settings.max_hits_per_database)
    s.setValue(f"{prefix}show_database_name_on_track", settings.show_database_name_on_track)
    s.setValue(f"{prefix}show_process_window", bool(settings.show_process_window))
    s.setValue(
        f"{prefix}close_process_window_when_finished",
        bool(settings.close_process_window_when_finished),
    )
