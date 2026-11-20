"""PathBrowseRow shows the end of long absolute paths when unfocused."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import QApplication

from repeatviewer.widgets.path_browse import PathBrowseRow


def test_path_browse_shows_end_of_long_path(qapp: QApplication, tmp_path: Path) -> None:
    deep = tmp_path
    for part in ("Documents", "RepeatViewer", "Epigenomics", "project_data"):
        deep = deep / part
        deep.mkdir()
    target = deep / "reads_nanopore_pass.fastq.gz"
    target.write_text("x")

    row = PathBrowseRow(minimum_width=180)
    row.resize(320, 32)
    row.show()
    qapp.processEvents()
    row.set_path(target)
    qapp.processEvents()

    edit = row.line_edit()
    assert not edit.hasFocus()
    assert edit.cursorPosition() == len(edit.text())
    assert edit.text().endswith(target.name)
    row.deleteLater()


def test_elide_middle_still_shows_basename(qapp: QApplication, tmp_path: Path) -> None:
    target = tmp_path / "very_long_alignment_file_name.fasta"
    target.write_text(">a\nACGT\n")
    row = PathBrowseRow(elide_middle=True, minimum_width=120)
    row.resize(200, 32)
    row.show()
    qapp.processEvents()
    row.set_path(target)
    qapp.processEvents()
    assert target.name in row.line_edit().text() or "…" in row.line_edit().text()
    assert row.path() == target.resolve()
    row.deleteLater()
