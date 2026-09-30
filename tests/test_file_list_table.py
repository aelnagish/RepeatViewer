"""File list table widget behavior."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QApplication

from repeatviewer.widgets.file_list_table import (
    _COL_DESCRIPTION,
    _COL_GC,
    _COL_NAME,
    FileListTableWidget,
)


def test_clear_rows_preserves_headers_and_refills(tmp_path: Path) -> None:
    app = QApplication.instance() or QApplication([])
    table = FileListTableWidget()
    first = tmp_path / "alpha.fa"
    second = tmp_path / "beta.fa"
    first.write_text(">a\nACGT\n")
    second.write_text(">b\nACGT\n")

    table.add_file_row(first)
    table.add_file_row(second)
    assert table.rowCount() == 2
    assert table.horizontalHeaderItem(_COL_NAME).text() == "Name"

    table.clear_rows()
    assert table.rowCount() == 0
    assert table.horizontalHeaderItem(_COL_NAME).text() == "Name"
    assert table.horizontalHeaderItem(_COL_GC).text() == "% GC"

    row_item = table.add_file_row(first)
    assert table.rowCount() == 1
    assert str(row_item.data(Qt.ItemDataRole.UserRole)) == str(first.resolve())
    row_item.setFont(row_item.font())


def test_clear_delegates_to_clear_rows(tmp_path: Path) -> None:
    app = QApplication.instance() or QApplication([])
    table = FileListTableWidget()
    first = tmp_path / "alpha.fa"
    first.write_text(">a\nACGT\n")
    table.add_file_row(first)
    table.clear()
    assert table.rowCount() == 0
    assert table.horizontalHeaderItem(_COL_NAME).text() == "Name"
    assert table.horizontalHeaderItem(_COL_GC).text() == "% GC"


def test_repair_chrome_restores_wiped_headers(tmp_path: Path) -> None:
    app = QApplication.instance() or QApplication([])
    table = FileListTableWidget()
    fasta = tmp_path / "x.fa"
    fasta.write_text(">a\nACGT\n")
    table.add_file_row(fasta)
    # Simulate the Qt bug: QTableWidget.clear() replaces labels with section numbers.
    super(FileListTableWidget, table).clear()
    assert table.horizontalHeaderItem(_COL_NAME) is None

    table.repair_chrome()
    assert table.horizontalHeaderItem(_COL_NAME).text() == "Name"
    assert table.horizontalHeaderItem(_COL_DESCRIPTION).text() == "Description"


def test_checkbox_selection_preferred_over_row_selection(tmp_path: Path) -> None:
    app = QApplication.instance() or QApplication([])
    table = FileListTableWidget()
    one = tmp_path / "one.fa"
    two = tmp_path / "two.fa"
    one.write_text(">a\nACGT\n")
    two.write_text(">b\nACGT\n")
    row_one = table.add_file_row(one)
    table.add_file_row(two)
    table.selectRow(1)
    row_one.setChecked(True)

    selected = table.selectedItems()
    assert len(selected) == 1
    assert str(selected[0].data(Qt.ItemDataRole.UserRole)) == str(one.resolve())


def test_header_sort_reorders_rows(tmp_path: Path) -> None:
    app = QApplication.instance() or QApplication([])
    table = FileListTableWidget()
    zebra = tmp_path / "zebra.fa"
    alpha = tmp_path / "alpha.fa"
    zebra.write_text(">z\nACGT\n")
    alpha.write_text(">a\nACGT\n")
    table.add_file_row(zebra)
    table.add_file_row(alpha)
    table.apply_default_sort()

    names = [
        str(table.item(row, _COL_NAME).text())
        for row in range(table.rowCount())
    ]
    assert names == ["alpha.fa", "zebra.fa"]

    header = table.horizontalHeader()
    header.sectionClicked.emit(_COL_NAME)
    names_desc = [
        str(table.item(row, _COL_NAME).text())
        for row in range(table.rowCount())
    ]
    assert names_desc == ["zebra.fa", "alpha.fa"]


def test_sorting_with_a_missing_sort_key_does_not_recurse(qapp):
    """A cell without a sort key (e.g. a just-pasted file) used to recurse until
    the app segfaulted: super().__lt__ dispatches back to the Python override."""
    import sys

    from PySide6.QtWidgets import QTableWidget

    from repeatviewer.widgets.file_list_table import _SortableTableItem

    table = QTableWidget(3, 1)
    table.setItem(0, 0, _SortableTableItem("b.fa", sort_key="b"))
    table.setItem(1, 0, _SortableTableItem("new.fa", sort_key=None))
    table.setItem(2, 0, _SortableTableItem("a.fa", sort_key="a"))
    old_limit = sys.getrecursionlimit()
    sys.setrecursionlimit(200)
    try:
        table.sortItems(0)
    finally:
        sys.setrecursionlimit(old_limit)
    assert sorted(table.item(r, 0).text() for r in range(3)) == ["a.fa", "b.fa", "new.fa"]
