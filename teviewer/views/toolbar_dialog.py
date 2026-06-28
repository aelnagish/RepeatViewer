"""Toolbar appearance settings."""

from __future__ import annotations

TOOLBAR_SETTINGS_GROUP = "toolbar"

DEFAULT_TOOLBAR_ORDER: list[str] = [
    "aliview",
    "_sep1",
    "align",
    "blast",
    "extend",
    "teaid",
    "msa_clean",
    "consensus",
    "flexidot",
    "te_genome_map",
    "te_benchmark",
]

TOOLBAR_LABELS: dict[str, str] = {
    "aliview": "AliView",
    "_sep1": "— separator —",
    "align": "Align (MAFFT)",
    "blast": "Blast",
    "te_genome_map": "TE visualization",
    "te_benchmark": "Benchmarking",
    "extend": "Extend",
    "teaid": "TEAid",
    "msa_clean": "MSA clean",
    "consensus": "Consensus",
    "flexidot": "Flexidot",
}


def load_toolbar_font_size() -> int:
    from PySide6.QtCore import QSettings

    s = QSettings()
    s.beginGroup(TOOLBAR_SETTINGS_GROUP)
    size = int(s.value("font_size", 11))
    s.endGroup()
    return max(8, min(16, size))


def save_toolbar_font_size(size: int) -> None:
    from PySide6.QtCore import QSettings

    s = QSettings()
    s.beginGroup(TOOLBAR_SETTINGS_GROUP)
    s.setValue("font_size", max(8, min(16, int(size))))
    s.endGroup()
