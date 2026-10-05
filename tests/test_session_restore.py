"""Session restore settings round-trip."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QSettings

from repeatviewer.database_settings import SETTINGS_GROUP
from repeatviewer.session_restore import (
    load_annotation_views,
    load_folder_file_memory,
    lookup_annotation_view,
    save_annotation_view,
    save_folder_file_memory,
)


def test_folder_file_memory_round_trip(tmp_path: Path) -> None:
    folder = tmp_path / "dir"
    folder.mkdir()
    file_path = folder / "seq.fa"
    file_path.write_text(">a\nACGT\n")

    s = QSettings()
    s.beginGroup(SETTINGS_GROUP)
    s.remove("folder_last_files_json")
    s.remove("folder_last_scroll_json")
    s.endGroup()

    key = str(folder.resolve())
    save_folder_file_memory(
        {key: str(file_path.resolve())},
        {key: 42},
    )
    files, scrolls = load_folder_file_memory()
    assert files[key].endswith("seq.fa")
    assert scrolls[key] == 42


def test_annotation_view_round_trip(tmp_path: Path) -> None:
    fasta = tmp_path / "x.fa"
    fasta.write_text(">q\nACGT\n")

    s = QSettings()
    s.beginGroup(SETTINGS_GROUP)
    s.remove("annotation_views_by_path_json")
    s.endGroup()

    state = {"start": 100, "end": 500, "ppb": 2.5, "scroll_px": 80.0, "sequence_id": "q"}
    save_annotation_view(fasta, state, sequence_id="q")
    hit = lookup_annotation_view(fasta, sequence_id="q")
    assert hit is not None
    assert hit["start"] == 100
    assert hit["end"] == 500
    views = load_annotation_views()
    assert len(views) == 1
