"""Example plugin: write a small stats file next to each selected FASTA."""

from __future__ import annotations

from pathlib import Path

from Bio import SeqIO

from teviewer.plugins.api import TEviewerPlugin
from teviewer.project import TEviewerProject


class SequenceStatsPlugin(TEviewerPlugin):
    id = "sequence_stats"
    label = "Seq stats"
    tooltip = "Write sequence count and total length beside each FASTA"
    fasta_only = True

    def validate(self, files: list[Path], project: TEviewerProject) -> str | None:
        if not files:
            return "Select at least one FASTA file."
        return None

    def run(self, file_path: Path, project: TEviewerProject, **job_options) -> list[str]:
        path = file_path.expanduser().resolve()
        records = list(SeqIO.parse(str(path), "fasta"))
        if not records:
            raise RuntimeError("No sequences in file.")
        total_len = sum(len(r.seq) for r in records)
        out = path.parent / f"{path.name}_stats.txt"
        out.write_text(
            f"file: {path.name}\n"
            f"sequences: {len(records)}\n"
            f"total_length: {total_len}\n",
            encoding="utf-8",
        )
        return [str(out)]
