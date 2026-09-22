"""File color persistence and Files-table Color column."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QSettings
from PySide6.QtWidgets import QApplication

from repeatviewer.widgets.file_colors import (
    get_file_color,
    replace_file_color_path,
    set_file_color,
)
from repeatviewer.widgets.file_list_table import (
    COL_COLOR,
    _COL_NAME,
    FileListTableWidget,
)


def test_file_color_persists_and_colors_name(tmp_path: Path, monkeypatch) -> None:
    app = QApplication.instance() or QApplication([])
    settings = QSettings()
    settings.clear()
    monkeypatch.setenv("QT_QPA_PLATFORM", "offscreen")

    fasta = tmp_path / "demo.fa"
    fasta.write_text(">a\nACGT\n")
    set_file_color(fasta, "#ea580c")
    assert get_file_color(fasta) == "#ea580c"

    table = FileListTableWidget()
    table.add_file_row(fasta)
    name = table.item(0, _COL_NAME)
    color_cell = table.item(0, COL_COLOR)
    assert name is not None
    assert color_cell is not None
    assert color_cell.text() == "Orange"
    assert name.foreground().color().name() == "#ea580c"
    assert color_cell.foreground().color().name() == "#ea580c"

    settings.clear()


def test_file_color_column_sort(tmp_path: Path, monkeypatch) -> None:
    app = QApplication.instance() or QApplication([])
    settings = QSettings()
    settings.clear()

    blue = tmp_path / "blue.fa"
    orange = tmp_path / "orange.fa"
    plain = tmp_path / "plain.fa"
    for path in (blue, orange, plain):
        path.write_text(">a\nACGT\n")
    set_file_color(blue, "#2563eb")
    set_file_color(orange, "#ea580c")

    table = FileListTableWidget()
    table.add_file_row(orange)
    table.add_file_row(plain)
    table.add_file_row(blue)
    assert table.item(1, COL_COLOR).text() == "—"
    table._sort_column = COL_COLOR
    table._sort_ascending = True
    table.apply_default_sort()

    names = [table.item(row, _COL_NAME).text() for row in range(table.rowCount())]
    # Blue (#2563eb) before Orange (#ea580c); empty last.
    assert names == ["blue.fa", "orange.fa", "plain.fa"]

    settings.clear()


def test_replace_file_color_path(tmp_path: Path) -> None:
    QApplication.instance() or QApplication([])
    settings = QSettings()
    settings.clear()

    old = tmp_path / "old.fa"
    new = tmp_path / "new.fa"
    old.write_text(">a\nACGT\n")
    set_file_color(old, "#16a34a")
    replace_file_color_path(old, new)
    assert get_file_color(old) is None
    assert get_file_color(new) == "#16a34a"

    settings.clear()
