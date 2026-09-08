from __future__ import annotations

from dataclasses import dataclass

MINIMAP2_PRESETS = (
    "map-ont",
    "map-hifi",
    "map-pb",
    "map-hifi-pb",
    "asm5",
    "asm10",
    "asm20",
    "splice",
    "sr",
)


@dataclass
class Minimap2Options:
    preset: str = "map-ont"
    threads: int = 4
    output_format: str = "paf"
    reference_fasta: str = ""
