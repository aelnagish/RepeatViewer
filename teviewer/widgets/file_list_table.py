"""Multi-column file table for the current folder."""

from __future__ import annotations

import platform
from pathlib import Path

from PySide6.QtCore import QMimeData, QSize, Qt, Signal
from PySide6.QtGui import QBrush, QDrag, QFont
from PySide6.QtWidgets import (
    QAbstractItemView,
    QHeaderView,
    QMenu,
    QTableView,
    QTableWidget,
    QTableWidgetItem,
    QTreeView,
    QWidget,
)

from teviewer.backends.blast_genome_plot import find_blast_table_for_query, is_te_hit_table_file
from teviewer.backends.separate_fasta import count_fasta_sequences
from teviewer.browser_file_metadata import BrowserFileMetadata, metadata_for_file
from teviewer.widgets.browser_icons import browser_file_icon

_FASTA_SUFFIXES = {".fa", ".fasta", ".fas"}
_ROLE_PATH = Qt.ItemDataRole.UserRole
_ROLE_SORT = Qt.ItemDataRole.UserRole + 1

_COL_CHECK = 0
_COL_NAME = 1
_COL_LENGTH = 2
_COL_SEQUENCES = 3
_COL_GC = 4
_COL_TYPE = 5
_COL_DESCRIPTION = 6
_COL_MODIFIED = 7
_COL_SIZE = 8

_HEADERS = (
    "",
    "Name",
    "Length",
    "# Sequences",
    "% GC",
    "Molecule Type",
    "Description",
    "Modified",
    "Size",
)
_CHECKBOX_COLUMN_WIDTH = 30
_DEFAULT_COLUMN_WIDTHS = {
    _COL_NAME: 200,
    _COL_LENGTH: 80,
    _COL_SEQUENCES: 96,
    _COL_GC: 64,
    _COL_TYPE: 110,
    _COL_DESCRIPTION: 220,
    _COL_MODIFIED: 150,
    _COL_SIZE: 72,
}
_LEFT_ALIGN = Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter


class _SortableTableItem(QTableWidgetItem):
    """Table cell that sorts by an explicit key stored in ``_ROLE_SORT``."""

    def __init__(self, text: str, *, sort_key: object) -> None:
        super().__init__(text)
        self.setData(_ROLE_SORT, sort_key)

    def __lt__(self, other: QTableWidgetItem) -> bool:  # type: ignore[override]
        if other is None:
            return super().__lt__(other)
        left = self.data(_ROLE_SORT)
        right = other.data(_ROLE_SORT)
        if left is None or right is None:
            return super().__lt__(other)
        try:
            return bool(left < right)
        except TypeError:
            return str(left).lower() < str(right).lower()


class FileListRow:
    """Row handle compatible with legacy QListWidgetItem styling helpers."""

    def __init__(self, table: "FileListTableWidget", row: int) -> None:
        self._table = table
        self._row = row

    def data(self, role: Qt.ItemDataRole) -> object | None:
        if role == _ROLE_PATH:
            return self._table._path_for_row(self._row)
        return None

    def setBackground(self, brush: QBrush) -> None:
        self._table._style_row(self._row, background=brush)

    def setForeground(self, brush: QBrush) -> None:
        self._table._style_row(self._row, foreground=brush)

    def setFont(self, font: QFont) -> None:
        self._table._style_row(self._row, font=font)

    def font(self) -> QFont:
        name_item = self._table.item(self._row, _COL_NAME)
        if name_item is not None:
            return name_item.font()
        return self._table.font()

    def setToolTip(self, text: str) -> None:
        name_item = self._table.item(self._row, _COL_NAME)
        if name_item is not None:
            name_item.setToolTip(text)

    def setSelected(self, selected: bool) -> None:
        self._table._set_row_selected(self._row, selected)

    def setChecked(self, checked: bool) -> None:
        self._table._set_row_checked(self._row, checked)

    def isChecked(self) -> bool:
        return self._table._is_row_checked(self._row)


class FileListTableWidget(QTableWidget):
    """Multi-column file list for the selected folder."""

    currentItemChanged = Signal(object, object)
    itemClicked = Signal(object)
    itemDoubleClicked = Signal(object)
    open_in_preview_requested = Signal(list)
    open_in_text_editor_requested = Signal(list)
    delete_files_requested = Signal(list)
    duplicate_files_requested = Signal(list)
    cut_files_requested = Signal(list)
    copy_files_requested = Signal(list)
    paste_files_requested = Signal()
    rename_file_requested = Signal(Path)
    separate_fasta_requested = Signal(Path)
    separate_fasta_subfolders_requested = Signal(Path)
    blast_genome_map_requested = Signal(Path)
    blast_divergence_requested = Signal(Path)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("FileList")
        self._sort_column = _COL_NAME
        self._sort_ascending = True
        self._ensure_headers()
        self.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.setShowGrid(False)
        self.setAlternatingRowColors(True)
        self.setWordWrap(False)
        self.setSortingEnabled(True)
        self.setDragEnabled(True)
        self.setDragDropMode(QAbstractItemView.DragDropMode.DragOnly)
        self.verticalHeader().setVisible(False)

        header = self.horizontalHeader()
        header.setStretchLastSection(False)
        header.setSectionsMovable(True)
        header.setSectionsClickable(True)
        header.setHighlightSections(True)
        header.setMinimumSectionSize(48)
        header.setDefaultAlignment(_LEFT_ALIGN)
        header.setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
        header.setSectionResizeMode(_COL_CHECK, QHeaderView.ResizeMode.Fixed)
        header.resizeSection(_COL_CHECK, _CHECKBOX_COLUMN_WIDTH)
        for column, width in _DEFAULT_COLUMN_WIDTHS.items():
            header.resizeSection(column, width)
        header.sectionClicked.connect(self._on_header_clicked)

        self.setIconSize(QSize(18, 18))
        self.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.customContextMenuRequested.connect(self._show_context_menu)
        self.cellClicked.connect(self._on_cell_clicked)
        self.cellDoubleClicked.connect(self._on_cell_double_clicked)
        self.currentCellChanged.connect(self._on_current_cell_changed)

    def _ensure_headers(self, *, force: bool = False) -> None:
        expected = len(_HEADERS)
        if self.columnCount() != expected:
            self.setColumnCount(expected)
        if force:
            self.setHorizontalHeaderLabels(list(_HEADERS))
            return
        for column, label in enumerate(_HEADERS):
            item = self.horizontalHeaderItem(column)
            if item is None or item.text() != label:
                self.setHorizontalHeaderLabels(list(_HEADERS))
                break

    def _reapply_cell_alignment(self) -> None:
        for row in range(self.rowCount()):
            for column in range(self.columnCount()):
                cell = super().item(row, column)
                if cell is None:
                    continue
                if column == _COL_CHECK:
                    cell.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                else:
                    cell.setTextAlignment(_LEFT_ALIGN)

    def _legitimate_header_views(self) -> set[QHeaderView]:
        """Header views still owned by table/tree widgets in this window."""
        headers: set[QHeaderView] = {self.horizontalHeader(), self.verticalHeader()}
        window = self.window()
        if window is None:
            return headers
        for table in window.findChildren(QTableView):
            headers.add(table.horizontalHeader())
            headers.add(table.verticalHeader())
        for tree in window.findChildren(QTreeView):
            headers.add(tree.header())
        return headers

    def repair_chrome(self) -> None:
        """Restore header labels and remove duplicate header views after style changes."""
        self._ensure_headers(force=True)
        header = self.horizontalHeader()
        vertical = self.verticalHeader()
        self.setHorizontalHeader(header)
        header.setVisible(True)
        vertical.setVisible(False)
        window = self.window()
        if window is not None:
            legit = self._legitimate_header_views()
            for orphan in window.findChildren(QHeaderView):
                if orphan in legit:
                    continue
                orphan.hide()
                orphan.deleteLater()
        self._reapply_cell_alignment()
        header.updateGeometries()
        self.updateGeometry()
        self.viewport().update()

    def clear(self) -> None:  # noqa: A003 (Qt API name)
        """Never call QTableWidget.clear() — it wipes headers and leaves ghost rows."""
        self.clear_rows()

    def clearContents(self) -> None:
        """Clear cell data without destroying column headers."""
        self.clear_rows()

    def clear_rows(self) -> None:
        """Remove all file rows without clearing column headers (QTableWidget.clear() wipes them)."""
        self.setRowCount(0)
        self._ensure_headers()

    def count(self) -> int:
        return self.rowCount()

    def item(self, row: int, column: int | None = None) -> FileListRow | QTableWidgetItem | None:
        """Row handle (legacy API) when ``column`` is omitted; else the table cell."""
        if column is None:
            if row < 0 or row >= self.rowCount():
                return None
            return FileListRow(self, row)
        return super().item(row, column)

    def currentItem(self) -> FileListRow | None:
        row = self.currentRow()
        return self.item(row) if row >= 0 else None

    def checkedItems(self) -> list[FileListRow]:
        rows: list[FileListRow] = []
        for row in range(self.rowCount()):
            if self._is_row_checked(row):
                rows.append(FileListRow(self, row))
        return rows

    def selectedItems(self) -> list[FileListRow]:
        checked = self.checkedItems()
        if checked:
            return checked
        rows = sorted({index.row() for index in self.selectionModel().selectedRows()})
        return [FileListRow(self, row) for row in rows]

    def setCurrentItem(self, row_item: FileListRow | None) -> None:
        if row_item is None:
            self.clearSelection()
            return
        self.selectRow(row_item._row)
        self.setCurrentCell(row_item._row, _COL_NAME)

    def scrollToItem(
        self,
        row_item: FileListRow,
        hint: QAbstractItemView.ScrollHint = QAbstractItemView.ScrollHint.EnsureVisible,
    ) -> None:
        cell = super().item(row_item._row, _COL_NAME)
        if cell is not None:
            super().scrollToItem(cell, hint)

    def add_file_row(self, path: Path) -> FileListRow:
        path = path.expanduser().resolve()
        was_sorting = self.isSortingEnabled()
        self.setSortingEnabled(False)
        row = self.rowCount()
        self.insertRow(row)
        meta = metadata_for_file(path)
        path_text = str(path)

        check_item = QTableWidgetItem()
        check_item.setFlags(
            Qt.ItemFlag.ItemIsUserCheckable
            | Qt.ItemFlag.ItemIsEnabled
            | Qt.ItemFlag.ItemIsSelectable
        )
        check_item.setCheckState(Qt.CheckState.Unchecked)
        check_item.setData(_ROLE_PATH, path_text)
        check_item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
        self.setItem(row, _COL_CHECK, check_item)

        name_item = _SortableTableItem(path.name, sort_key=meta.sort_name)
        name_item.setData(_ROLE_PATH, path_text)
        name_item.setIcon(browser_file_icon(self, path))
        name_item.setTextAlignment(_LEFT_ALIGN)
        self.setItem(row, _COL_NAME, name_item)

        for column, text, sort_key in (
            (_COL_LENGTH, meta.length_text, meta.sort_length),
            (_COL_SEQUENCES, meta.sequence_count_text, meta.sort_sequence_count),
            (_COL_GC, meta.gc_text, meta.sort_gc),
            (_COL_TYPE, meta.molecule_type, meta.molecule_type.lower()),
            (_COL_DESCRIPTION, meta.description_text, meta.sort_description),
            (_COL_MODIFIED, meta.modified_text, meta.sort_modified),
            (_COL_SIZE, meta.size_text, meta.sort_size),
        ):
            cell = _SortableTableItem(text, sort_key=sort_key)
            cell.setData(_ROLE_PATH, path_text)
            cell.setTextAlignment(_LEFT_ALIGN)
            self.setItem(row, column, cell)

        self.setSortingEnabled(was_sorting)
        return FileListRow(self, row)

    def apply_default_sort(self) -> None:
        order = Qt.SortOrder.AscendingOrder if self._sort_ascending else Qt.SortOrder.DescendingOrder
        self.sortItems(self._sort_column, order)

    def add_placeholder_row(self, text: str) -> None:
        was_sorting = self.isSortingEnabled()
        self.setSortingEnabled(False)
        row = self.rowCount()
        self.insertRow(row)
        cell = QTableWidgetItem(text)
        cell.setFlags(Qt.ItemFlag.NoItemFlags)
        self.setItem(row, _COL_NAME, cell)
        self.setSpan(row, _COL_CHECK, 1, self.columnCount())
        self.setSortingEnabled(was_sorting)

    def _set_row_selected(self, row: int, selected: bool) -> None:
        if row < 0 or row >= self.rowCount():
            return
        index = self.model().index(row, 0)
        flags = (
            QAbstractItemView.SelectionFlag.Select
            | QAbstractItemView.SelectionFlag.Rows
        )
        if selected:
            self.selectionModel().select(index, flags)
        else:
            self.selectionModel().select(index, QAbstractItemView.SelectionFlag.Deselect)

    def _set_row_checked(self, row: int, checked: bool) -> None:
        cell = self.item(row, _COL_CHECK)
        if cell is None:
            return
        cell.setCheckState(
            Qt.CheckState.Checked if checked else Qt.CheckState.Unchecked
        )

    def _is_row_checked(self, row: int) -> bool:
        if row < 0 or row >= self.rowCount():
            return False
        cell = self.item(row, _COL_CHECK)
        return cell is not None and cell.checkState() == Qt.CheckState.Checked

    def _path_for_row(self, row: int) -> str | None:
        if row < 0 or row >= self.rowCount():
            return None
        item = self.item(row, _COL_NAME)
        if item is None:
            item = self.item(row, _COL_CHECK)
        if item is None:
            return None
        raw = item.data(_ROLE_PATH)
        return str(raw) if raw else None

    def _style_row(
        self,
        row: int,
        *,
        background: QBrush | None = None,
        foreground: QBrush | None = None,
        font: QFont | None = None,
    ) -> None:
        for column in range(self.columnCount()):
            cell = self.item(row, column)
            if cell is None:
                continue
            if background is not None:
                cell.setBackground(background)
            if foreground is not None:
                cell.setForeground(foreground)
            if font is not None:
                cell.setFont(font)

    def _selected_paths(self) -> list[Path]:
        paths: list[Path] = []
        for row_item in self.selectedItems():
            raw = row_item.data(_ROLE_PATH)
            if not raw:
                continue
            path = Path(str(raw))
            if path.is_file():
                paths.append(path)
        return paths

    def _on_header_clicked(self, column: int) -> None:
        if column == _COL_CHECK:
            return
        if self._sort_column == column:
            self._sort_ascending = not self._sort_ascending
        else:
            self._sort_column = column
            self._sort_ascending = True
        self.apply_default_sort()

    def _on_cell_clicked(self, row: int, column: int) -> None:
        if column == _COL_CHECK:
            return
        row_item = self.item(row)
        if row_item is not None:
            self.itemClicked.emit(row_item)

    def _on_cell_double_clicked(self, row: int, column: int) -> None:
        if column == _COL_CHECK:
            return
        row_item = self.item(row)
        if row_item is not None:
            self.itemDoubleClicked.emit(row_item)

    def _on_current_cell_changed(
        self,
        current_row: int,
        _current_col: int,
        previous_row: int,
        _previous_col: int,
    ) -> None:
        current = self.item(current_row) if current_row >= 0 else None
        previous = self.item(previous_row) if previous_row >= 0 else None
        self.currentItemChanged.emit(current, previous)

    def _show_context_menu(self, position) -> None:
        index = self.indexAt(position)
        if index.isValid() and index.row() not in {r._row for r in self.selectedItems()}:
            self.selectRow(index.row())
        paths = self._selected_paths()
        window = self.window()
        can_paste = (
            hasattr(window, "_clipboard_can_paste")
            and window._clipboard_can_paste()
            and hasattr(window, "_current_folder")
            and window._current_folder is not None
        )
        if not paths and not can_paste:
            return
        menu = QMenu(self)
        if paths:
            if len(paths) == 1:
                path = paths[0]
                menu.addAction("Open in preview", lambda: self.open_in_preview_requested.emit(paths))
                menu.addAction(
                    "Open in text editor",
                    lambda: self.open_in_text_editor_requested.emit(paths),
                )
                menu.addAction("Rename", lambda: self.rename_file_requested.emit(path))
                if path.suffix.lower() in _FASTA_SUFFIXES:
                    try:
                        seq_n = count_fasta_sequences(path)
                    except Exception:  # noqa: BLE001
                        seq_n = 0
                    if seq_n > 1:
                        menu.addAction(
                            f"Separate sequence into {seq_n} files",
                            lambda: self.separate_fasta_requested.emit(path),
                        )
                        menu.addAction(
                            f"Separate sequence into {seq_n} folders",
                            lambda: self.separate_fasta_subfolders_requested.emit(path),
                        )
                    blast_txt = find_blast_table_for_query(path)
                    if blast_txt is not None:
                        menu.addAction(
                            "TE density plot",
                            lambda bt=blast_txt: self.blast_genome_map_requested.emit(bt),
                        )
                if is_te_hit_table_file(path):
                    menu.addAction(
                        "TE density plot",
                        lambda bt=path: self.blast_genome_map_requested.emit(bt),
                    )
                    menu.addAction(
                        "TE divergence plot",
                        lambda bt=path: self.blast_divergence_requested.emit(bt),
                    )
            else:
                menu.addAction(
                    f"Open {len(paths)} files in tabs",
                    lambda: self.open_in_preview_requested.emit(paths),
                )
                menu.addAction(
                    f"Open {len(paths)} files in text editor",
                    lambda: self.open_in_text_editor_requested.emit(paths),
                )
            menu.addSeparator()
            menu.addAction("Cut", lambda: self.cut_files_requested.emit(paths))
            menu.addAction("Copy", lambda: self.copy_files_requested.emit(paths))
            if can_paste:
                menu.addAction("Paste", lambda: self.paste_files_requested.emit())
            menu.addSeparator()
            if len(paths) == 1:
                menu.addAction(
                    "Duplicate",
                    lambda: self.duplicate_files_requested.emit(paths),
                )
            else:
                menu.addAction(
                    f"Duplicate {len(paths)} files",
                    lambda: self.duplicate_files_requested.emit(paths),
                )
            menu.addSeparator()
            menu.addAction("Open in AliView", self._context_open_aliview)
            if len(paths) == 1:
                label = (
                    "Reveal in Finder"
                    if platform.system() == "Darwin"
                    else "Show in file manager"
                )
                menu.addAction(label, lambda: self._reveal_paths(paths))
            menu.addSeparator()
            if len(paths) == 1:
                menu.addAction("Delete", lambda: self.delete_files_requested.emit(paths))
            else:
                menu.addAction(
                    f"Delete {len(paths)} files",
                    lambda: self.delete_files_requested.emit(paths),
                )
        elif can_paste:
            menu.addAction("Paste", lambda: self.paste_files_requested.emit())
        menu.exec(self.mapToGlobal(position))

    def _context_open_aliview(self) -> None:
        window = self.window()
        if hasattr(window, "_open_selected_in_aliview"):
            window._open_selected_in_aliview()

    def _reveal_paths(self, paths: list[Path]) -> None:
        window = self.window()
        if hasattr(window, "_reveal_in_file_manager"):
            window._reveal_in_file_manager(paths)

    def mimeData(self, items: list[FileListRow]) -> QMimeData:  # noqa: N802 (Qt API)
        mime = QMimeData()
        urls = []
        for row_item in items:
            raw = row_item.data(_ROLE_PATH)
            if raw:
                urls.append(Path(str(raw)).as_uri())
        if urls:
            from PySide6.QtCore import QUrl

            mime.setUrls([QUrl(u) for u in urls])
        return mime

    def startDrag(self, supportedActions):  # noqa: N802 (Qt API)
        items = self.selectedItems()
        if not items:
            return
        drag = QDrag(self)
        drag.setMimeData(self.mimeData(items))
        drag.exec(supportedActions)
