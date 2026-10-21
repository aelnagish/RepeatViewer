"""Multi-column file table for the current folder."""

from __future__ import annotations

import platform
from pathlib import Path

from PySide6.QtCore import QItemSelectionModel, QMimeData, QSize, Qt, QTimer, Signal
from PySide6.QtGui import QBrush, QColor, QDrag, QFont, QFontMetrics, QMouseEvent, QPainter, QPixmap
from PySide6.QtWidgets import (
    QAbstractItemView,
    QApplication,
    QColorDialog,
    QDialog,
    QHeaderView,
    QMenu,
    QTableView,
    QTableWidget,
    QTableWidgetItem,
    QTreeView,
    QWidget,
)

from repeatviewer.backends.blast_genome_plot import (
    find_blast_table_for_query,
    is_blast_table_file,
    is_te_hit_table_file,
)
from repeatviewer.backends.cyto_images import image_count_for_sequence
from repeatviewer.backends.separate_fasta import count_fasta_sequences
from repeatviewer.browser_file_metadata import BrowserFileMetadata
from repeatviewer.ui.theme import item_view_checkbox_stylesheet
from repeatviewer.widgets.browser_icons import browser_file_icon

_FASTA_SUFFIXES = {".fa", ".fasta", ".fas", ".faa"}
_ROLE_PATH = Qt.ItemDataRole.UserRole
_ROLE_SORT = Qt.ItemDataRole.UserRole + 1
_ROLE_FAMILY_BAND = Qt.ItemDataRole.UserRole + 2  # True = tinted name-group band

_COL_CHECK = 0
_COL_NAME = 1
_COL_LENGTH = 2
_COL_SEQUENCES = 3
_COL_GC = 4
_COL_TYPE = 5
_COL_DESCRIPTION = 6
_COL_MODIFIED = 7
_COL_SIZE = 8
_COL_COLOR = 9
_COL_IMAGES = 10

# Public aliases for column visibility settings.
COL_CHECK = _COL_CHECK
COL_NAME = _COL_NAME
COL_LENGTH = _COL_LENGTH
COL_SEQUENCES = _COL_SEQUENCES
COL_GC = _COL_GC
COL_TYPE = _COL_TYPE
COL_DESCRIPTION = _COL_DESCRIPTION
COL_MODIFIED = _COL_MODIFIED
COL_SIZE = _COL_SIZE
COL_COLOR = _COL_COLOR
COL_IMAGES = _COL_IMAGES

_HEADERS = (
    "",
    "Name",
    "Length",
    "# Sequences",
    "% GC",
    "Type",
    "Description",
    "Modified",
    "Size",
    "Color",
    "\U0001F52C",  # microscope: sequences with cytogenetic images
)
_CHECKBOX_COLUMN_WIDTH = 22
_DEFAULT_COLUMN_WIDTHS = {
    _COL_NAME: 200,
    _COL_LENGTH: 80,
    _COL_SEQUENCES: 96,
    _COL_GC: 64,
    _COL_TYPE: 110,
    _COL_DESCRIPTION: 220,
    _COL_MODIFIED: 150,
    _COL_SIZE: 72,
    _COL_COLOR: 80,
    _COL_IMAGES: 34,
}
_LEFT_ALIGN = Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter


class _SortableTableItem(QTableWidgetItem):
    """Table cell that sorts by an explicit key stored in ``_ROLE_SORT``."""

    def __init__(self, text: str, *, sort_key: object) -> None:
        super().__init__(text)
        self.setData(_ROLE_SORT, sort_key)

    def __lt__(self, other: QTableWidgetItem) -> bool:  # type: ignore[override]
        # Never call super().__lt__ here: in PySide6 Qt's operator< is virtual
        # and dispatches straight back to this method, so it recursed until the
        # app crashed (segfault) whenever a cell had no sort key yet, e.g. a
        # file that was just pasted. Compare the visible text instead.
        if other is None:
            return False
        left = self.data(_ROLE_SORT)
        right = other.data(_ROLE_SORT)
        if left is None or right is None:
            return self.text().lower() < other.text().lower()
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
    # Emitted when a Files-list drag starts (copy/move to Folders); cancel Preview.
    dragStarted = Signal()
    open_in_preview_requested = Signal(list)
    open_in_preview_tab_requested = Signal(list)
    open_in_text_editor_requested = Signal(list)
    delete_files_requested = Signal(list)
    duplicate_files_requested = Signal(list)
    cut_files_requested = Signal(list)
    copy_files_requested = Signal(list)
    paste_files_requested = Signal()
    rename_file_requested = Signal(Path)
    batch_rename_requested = Signal(list)  # Edit > Batch Rename… for these files
    separate_fasta_requested = Signal(Path)
    separate_fasta_subfolders_requested = Signal(Path)
    extract_sequences_requested = Signal(Path)      # pick records out of one file
    extract_hit_subjects_requested = Signal(Path)  # hit table -> subject sequences
    group_sequences_requested = Signal(list)        # combine several files
    new_folder_with_files_requested = Signal(list)  # copy selection into a new folder
    blast_genome_map_requested = Signal(Path)
    blast_divergence_requested = Signal(Path)
    add_blast_database_requested = Signal(Path)
    add_annotation_requested = Signal(Path)
    export_annotation_requested = Signal(Path)
    cyto_images_requested = Signal(Path)             # open the image(s) of a sequence
    attach_cyto_image_requested = Signal(Path, list)  # FASTA, image paths ([] = ask)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("FileList")
        self._sort_column = _COL_NAME
        self._sort_ascending = True
        self._suppress_selection_preview = False
        # After a drag, Qt often fires cellClicked on mouse release — skip Preview.
        self._skip_click_preview_once = False
        # Name sort: tint source + tool-output families. Other sorts: even/odd rows.
        self._family_banding_enabled = False
        self._row_zebra_enabled = False
        # While True, metadata chunk updates must not re-enable sorting (avoids
        # scroll jumps when Type/Length values change mid-fill on large folders).
        self._metadata_sort_frozen = False
        self._metadata_saved_scroll: int | None = None
        # When set, skip per-row QSettings reads during bulk listing.
        self._bulk_file_colors: dict[str, str] | None = None
        self._ensure_headers()
        self.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.setShowGrid(False)
        self.setAlternatingRowColors(False)
        self.setWordWrap(False)
        # Long names: keep start + end (extension), elide the middle — like Finder.
        self.setTextElideMode(Qt.TextElideMode.ElideMiddle)
        # Selected rows use white text; keep checkbox ticks black on the light box.
        self.setStyleSheet(item_view_checkbox_stylesheet())
        self.setSortingEnabled(True)
        self.setDragEnabled(True)
        self.setDragDropMode(QAbstractItemView.DragDropMode.DragOnly)
        # Image files dropped from Finder / Explorer onto a FASTA row are attached
        # as cytogenetic images (see dragEnterEvent / dropEvent).
        self.setAcceptDrops(True)
        self.viewport().setAcceptDrops(True)
        self._cyto_drop_row: int | None = None
        self.verticalHeader().setVisible(False)

        header = self.horizontalHeader()
        header.setStretchLastSection(False)
        header.setSectionsMovable(True)
        header.setSectionsClickable(True)
        header.setHighlightSections(True)
        header.setMinimumSectionSize(22)
        header.setFixedHeight(18)
        header.setDefaultAlignment(_LEFT_ALIGN)
        header.setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
        header.setSectionResizeMode(_COL_CHECK, QHeaderView.ResizeMode.Fixed)
        # Allow a narrow checkbox column; other columns keep readable widths below.
        header.setMinimumSectionSize(22)
        header.resizeSection(_COL_CHECK, _CHECKBOX_COLUMN_WIDTH)
        for column, width in _DEFAULT_COLUMN_WIDTHS.items():
            header.resizeSection(column, width)
        header.sectionClicked.connect(self._on_header_clicked)
        header.sectionMoved.connect(self._schedule_save_column_layout)
        header.sectionResized.connect(self._schedule_save_column_layout)

        self._column_layout_save_timer = QTimer(self)
        self._column_layout_save_timer.setSingleShot(True)
        self._column_layout_save_timer.setInterval(400)
        self._column_layout_save_timer.timeout.connect(self._save_column_layout)
        self._restoring_column_layout = False

        self.setIconSize(QSize(18, 18))
        self.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.customContextMenuRequested.connect(self._show_context_menu)
        self.cellClicked.connect(self._on_cell_clicked)
        self.cellDoubleClicked.connect(self._on_cell_double_clicked)
        self.currentCellChanged.connect(self._on_current_cell_changed)

        from repeatviewer.widgets.file_list_columns import load_visible_file_columns

        self.apply_visible_columns(load_visible_file_columns())
        self._restore_column_layout()
        self._place_images_column_once()

    def _schedule_save_column_layout(self, *_args) -> None:
        if self._restoring_column_layout:
            return
        self._column_layout_save_timer.start()

    def _save_column_layout(self) -> None:
        from repeatviewer.widgets.file_list_columns import save_column_order, save_column_widths

        header = self.horizontalHeader()
        order: list[int] = []
        widths: dict[int, int] = {}
        for visual in range(header.count()):
            logical = header.logicalIndex(visual)
            order.append(logical)
            widths[logical] = header.sectionSize(logical)
        save_column_order(order)
        save_column_widths(widths)

    def _restore_column_layout(self) -> None:
        from repeatviewer.widgets.file_list_columns import (
            load_column_order,
            load_column_widths,
        )

        header = self.horizontalHeader()
        order = load_column_order()
        widths = load_column_widths()
        if not order and not widths:
            return

        self._restoring_column_layout = True
        try:
            header.blockSignals(True)
            if order:
                for target_visual, logical in enumerate(order):
                    if logical < 0 or logical >= header.count():
                        continue
                    current_visual = header.visualIndex(logical)
                    if current_visual >= 0 and current_visual != target_visual:
                        header.moveSection(current_visual, target_visual)
            for logical, width in widths.items():
                if 0 <= logical < header.count() and logical != _COL_CHECK:
                    header.resizeSection(logical, width)
            header.resizeSection(_COL_CHECK, _CHECKBOX_COLUMN_WIDTH)
        finally:
            header.blockSignals(False)
            self._restoring_column_layout = False

    def _place_images_column_once(self) -> None:
        """First run with the images column: put it right after Name, where it is seen.

        One-shot (QSettings flag), so a user who later drags it elsewhere keeps
        their order.
        """
        from PySide6.QtCore import QSettings

        from repeatviewer.widgets.file_list_columns import FILES_SETTINGS_GROUP

        s = QSettings()
        s.beginGroup(FILES_SETTINGS_GROUP)
        done = str(s.value("images_column_placed_v2", False)).lower() in {"1", "true"}
        s.endGroup()
        if done:
            return
        header = self.horizontalHeader()
        target = header.visualIndex(_COL_NAME) + 1
        current = header.visualIndex(_COL_IMAGES)
        if current >= 0 and current != target:
            header.moveSection(current, target)
        s = QSettings()
        s.beginGroup(FILES_SETTINGS_GROUP)
        s.setValue("images_column_placed_v2", True)
        s.endGroup()

    def apply_visible_columns(self, visible: set[int]) -> None:
        """Show/hide optional columns (checkbox + name always stay visible)."""
        from repeatviewer.widgets.file_list_columns import (
            ALWAYS_VISIBLE_COLUMNS,
            TOGGLEABLE_COLUMNS,
        )

        allowed = ALWAYS_VISIBLE_COLUMNS | set(TOGGLEABLE_COLUMNS.keys())
        show = ALWAYS_VISIBLE_COLUMNS | {column for column in visible if column in allowed}
        for column in range(self.columnCount()):
            self.setColumnHidden(column, column not in show)

    def _ensure_headers(self, *, force: bool = False) -> None:
        expected = len(_HEADERS)
        if self.columnCount() != expected:
            self.setColumnCount(expected)
        if force:
            self.setHorizontalHeaderLabels(list(_HEADERS))
            images_header = self.horizontalHeaderItem(_COL_IMAGES)
            if images_header is not None:
                images_header.setToolTip("Cytogenetic images attached (click a tick to open)")
            return
        for column, label in enumerate(_HEADERS):
            item = self.horizontalHeaderItem(column)
            if item is None or item.text() != label:
                self.setHorizontalHeaderLabels(list(_HEADERS))
                break
        images_header = self.horizontalHeaderItem(_COL_IMAGES)
        if images_header is not None:
            images_header.setToolTip("Cytogenetic images attached (click a tick to open)")
            images_header.setTextAlignment(Qt.AlignmentFlag.AlignCenter)

    def _reapply_cell_alignment(self) -> None:
        for row in range(self.rowCount()):
            for column in range(self.columnCount()):
                cell = super().item(row, column)
                if cell is None:
                    continue
                if column in (_COL_CHECK, _COL_IMAGES):
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
            self.setCurrentCell(-1, -1)
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

    def begin_bulk_file_insert(self, folder: Path | None = None) -> None:
        """Cache file colors (and cyto counts) once for a folder rebuild."""
        from repeatviewer.widgets.file_colors import load_file_colors

        self._bulk_file_colors = load_file_colors()
        self._bulk_cyto_counts: dict[str, int] | None = {}
        if folder is not None:
            try:
                from repeatviewer.backends.cyto_images import image_counts_by_sequence

                self._bulk_cyto_counts = image_counts_by_sequence(folder)
            except Exception:  # noqa: BLE001 - listing must not fail on cyto
                self._bulk_cyto_counts = {}

    def end_bulk_file_insert(self) -> None:
        self._bulk_file_colors = None
        self._bulk_cyto_counts = None

    def add_file_row(
        self,
        path: Path,
        *,
        metadata: BrowserFileMetadata | None = None,
        fill_metadata_columns: bool | None = None,
    ) -> FileListRow:
        path = path.expanduser()
        # During bulk folder rebuild paths are already absolute from scandir —
        # skip resolve() (N syscalls) so ROLE_PATH matches cached stats keys.
        if self._bulk_file_colors is None or not path.is_absolute():
            try:
                path = path.resolve()
            except OSError:
                pass
        # Caller may already have sorting disabled for a bulk insert; avoid
        # toggling it per row (expensive when adding hundreds of files).
        was_sorting = self.isSortingEnabled()
        if was_sorting:
            self.setSortingEnabled(False)
        row = self.rowCount()
        self.insertRow(row)
        if fill_metadata_columns is None:
            from repeatviewer.widgets.file_list_columns import metadata_columns_visible

            fill_metadata_columns = metadata_columns_visible()
        if metadata is None:
            if fill_metadata_columns:
                from repeatviewer.browser_file_metadata import quick_metadata_for_file

                meta = quick_metadata_for_file(path, resolve_path=False)
            else:
                from repeatviewer.browser_file_metadata import minimal_name_metadata

                meta = minimal_name_metadata(path)
        else:
            meta = metadata
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
        name_item.setToolTip(path.name)
        self.setItem(row, _COL_NAME, name_item)

        self._apply_file_color_cells(row, path)
        # Small folders: fill metadata now. Large folders: name only here;
        # deferred fill creates Modified/Type/Size (and later Length/GC).
        if fill_metadata_columns:
            self._apply_metadata_cells(row, path_text, meta)
        self._apply_cyto_cell(row, path)
        if was_sorting:
            self.setSortingEnabled(True)
        return FileListRow(self, row)

    def _apply_file_color_cells(
        self,
        row: int,
        path: Path,
        *,
        resolved: bool = False,
    ) -> None:
        """Apply stored file color to Name, Color, and all metadata cells.

        Also clears job highlight backgrounds/bold so previous green/blue tints
        do not linger after a new tool run.

        Pass ``resolved=True`` when ``path`` is already absolute/resolved
        (e.g. ROLE_PATH) to avoid a redundant ``resolve()`` syscall.
        """
        from repeatviewer.widgets.file_colors import (
            color_label,
            get_file_color,
        )

        if not resolved:
            path = path.expanduser().resolve()
        path_text = str(path)
        if self._bulk_file_colors is not None:
            color = self._bulk_file_colors.get(path_text)
        else:
            color = get_file_color(path)
        brush = QBrush(QColor(color)) if color else QBrush()
        # Keep zebra / family band tint under the colored text.
        band_bg = self._family_band_brush(row)

        label = color_label(color) if color else "—"
        # Empty colors sort last when ascending.
        sort_key = color.lower() if color else "\uffff"
        cell = _SortableTableItem(label, sort_key=sort_key)
        cell.setData(_ROLE_PATH, path_text)
        cell.setTextAlignment(_LEFT_ALIGN)
        cell.setForeground(brush)
        cell.setBackground(band_bg)
        self.setItem(row, _COL_COLOR, cell)

        # Name + every other populated cell (Length, GC, …) update immediately.
        for column in range(self.columnCount()):
            existing = self.item(row, column)
            if existing is None:
                continue
            existing.setBackground(band_bg)
            existing.setForeground(brush)
            font = existing.font()
            if font.bold():
                font.setBold(False)
                existing.setFont(font)
        name_item = self.item(row, _COL_NAME)
        if name_item is not None:
            name_item.setToolTip("")

    def set_family_banding_enabled(self, enabled: bool) -> None:
        """Enable/disable name-family zebra banding (online-BLAST style)."""
        enabled = bool(enabled)
        if enabled == self._family_banding_enabled:
            if enabled:
                self.apply_family_group_banding()
            return
        self._family_banding_enabled = enabled
        if enabled:
            self._row_zebra_enabled = False
            self.apply_family_group_banding()
            return
        self._clear_family_band_flags()
        # Only clear band tints — do not reload every row's file color (slow on
        # large folders and raced with deferred metadata fill).
        self._clear_band_backgrounds()

    def family_banding_enabled(self) -> bool:
        return self._family_banding_enabled

    def row_zebra_enabled(self) -> bool:
        return self._row_zebra_enabled

    def set_list_banding(self, *, family: bool, row_zebra: bool) -> None:
        """Choose name-family banding, simple row zebra, or neither."""
        family = bool(family)
        row_zebra = bool(row_zebra) and not family
        changed = (
            family != self._family_banding_enabled
            or row_zebra != self._row_zebra_enabled
        )
        self._family_banding_enabled = family
        self._row_zebra_enabled = row_zebra
        if not family:
            self._clear_family_band_flags()
        if not family and not row_zebra:
            if changed:
                self._clear_band_backgrounds()
            return
        self.apply_list_banding()

    def apply_list_banding(self) -> None:
        if self._family_banding_enabled:
            self.apply_family_group_banding()
        elif self._row_zebra_enabled:
            self.apply_row_zebra_banding()

    def _clear_family_band_flags(self) -> None:
        for row in range(self.rowCount()):
            name_item = self.item(row, _COL_NAME)
            if name_item is not None:
                name_item.setData(_ROLE_FAMILY_BAND, None)

    def _clear_band_backgrounds(self) -> None:
        clear = QBrush()
        for row in range(self.rowCount()):
            for column in range(self.columnCount()):
                cell = self.item(row, column)
                if cell is not None:
                    cell.setBackground(clear)

    def _band_tint_hex(self, *, family: bool) -> str:
        from PySide6.QtGui import QPalette
        from repeatviewer.browser_file_groups import (
            FAMILY_BAND_DARK,
            FAMILY_BAND_LIGHT,
            ROW_ZEBRA_DARK,
            ROW_ZEBRA_LIGHT,
        )

        window = self.palette().color(QPalette.ColorRole.Window)
        dark = window.lightness() < 128
        if family:
            return FAMILY_BAND_DARK if dark else FAMILY_BAND_LIGHT
        return ROW_ZEBRA_DARK if dark else ROW_ZEBRA_LIGHT

    def _family_band_brush(self, row: int) -> QBrush:
        if self._family_banding_enabled:
            name_item = self.item(row, _COL_NAME)
            if name_item is None or not name_item.data(_ROLE_FAMILY_BAND):
                return QBrush()
            return QBrush(QColor(self._band_tint_hex(family=True)))
        if self._row_zebra_enabled and row % 2 == 1:
            return QBrush(QColor(self._band_tint_hex(family=False)))
        return QBrush()

    def apply_family_group_banding(self) -> None:
        """Tint alternating source+output families when name banding is on."""
        if not self._family_banding_enabled:
            return
        from repeatviewer.browser_file_groups import family_band_flags

        names: list[str] = []
        rows: list[int] = []
        for row in range(self.rowCount()):
            name_item = self.item(row, _COL_NAME)
            if name_item is None:
                continue
            text = name_item.text()
            if not text or text.startswith("("):
                continue
            names.append(text)
            rows.append(row)
        flags = family_band_flags(names)
        for row, tinted in zip(rows, flags):
            name_item = self.item(row, _COL_NAME)
            if name_item is not None:
                name_item.setData(_ROLE_FAMILY_BAND, bool(tinted))
            brush = self._family_band_brush(row)
            for column in range(self.columnCount()):
                cell = self.item(row, column)
                if cell is None:
                    continue
                # Background only — keep user file-color text if set.
                cell.setBackground(brush)

    def apply_row_zebra_banding(self) -> None:
        """Tint every other row when the list is not sorted by Name."""
        if not self._row_zebra_enabled:
            return
        for row in range(self.rowCount()):
            name_item = self.item(row, _COL_NAME)
            if name_item is None:
                continue
            text = name_item.text()
            if not text or text.startswith("("):
                continue
            brush = self._family_band_brush(row)
            for column in range(self.columnCount()):
                cell = self.item(row, column)
                if cell is not None:
                    # Background only — keep user file-color text if set.
                    cell.setBackground(brush)

    def reapply_file_colors(self) -> None:
        """Refresh Name/Color styling from stored colors for every row."""
        for row in range(self.rowCount()):
            raw = self._path_for_row(row)
            if not raw:
                continue
            # ROLE_PATH is already resolved — skip per-row resolve().
            self._apply_file_color_cells(row, Path(raw), resolved=True)

    def reapply_file_color_at_row(self, row: int) -> None:
        """Clear job tint for one row without rebuilding the path index."""
        raw = self._path_for_row(row)
        if raw:
            self._apply_file_color_cells(row, Path(raw), resolved=True)

    def reapply_file_color_for_path(self, path: Path) -> None:
        key = str(path.expanduser().resolve())
        row = self._path_row_index().get(key)
        if row is not None:
            self._apply_file_color_cells(row, Path(key), resolved=True)

    def reapply_file_colors_for_keys(self, keys: set[str]) -> None:
        """Clear job tint for specific path keys (one index build)."""
        if not keys:
            return
        row_index = self._path_row_index()
        for key in keys:
            row = row_index.get(key)
            if row is not None:
                self._apply_file_color_cells(row, Path(key), resolved=True)

    # -- cytogenetic images ------------------------------------------------

    def _apply_cyto_cell(self, row: int, path: Path) -> None:
        count = 0
        if path.suffix.lower() in _FASTA_SUFFIXES:
            try:
                count = image_count_for_sequence(path)
            except Exception:  # noqa: BLE001 - never break the listing over images
                count = 0
        cell = _SortableTableItem("\u2713" if count else "", sort_key=count)
        cell.setData(_ROLE_PATH, str(path))
        cell.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
        if count:
            cell.setToolTip(
                f"{count} cytogenetic image{'s' if count != 1 else ''} \u2014 click to open"
            )
            bold = cell.font()
            bold.setBold(True)
            cell.setFont(bold)
        self.setItem(row, _COL_IMAGES, cell)

    def cyto_image_count(self, row: int) -> int:
        cell = super().item(row, _COL_IMAGES)
        if cell is None:
            return 0
        value = cell.data(_ROLE_SORT)
        return int(value) if isinstance(value, int) else 0

    def refresh_cyto_images(self) -> None:
        """Re-read the tick column (after an image was attached or removed)."""
        was_sorting = self.isSortingEnabled()
        self.setSortingEnabled(False)
        try:
            for row in range(self.rowCount()):
                raw = self._path_for_row(row)
                if raw:
                    self._apply_cyto_cell(row, Path(raw))
        finally:
            self.setSortingEnabled(was_sorting)

    def _fasta_row_at(self, pos) -> int | None:
        index = self.indexAt(pos)
        if not index.isValid():
            return None
        raw = self._path_for_row(index.row())
        if not raw or Path(raw).suffix.lower() not in _FASTA_SUFFIXES:
            return None
        return index.row()

    def dragEnterEvent(self, event) -> None:  # noqa: N802 (Qt API)
        from repeatviewer.widgets.cyto_images_ui import image_paths_from_mime

        if event.source() is None and image_paths_from_mime(event.mimeData()):
            event.acceptProposedAction()
            return
        event.ignore()

    def dragMoveEvent(self, event) -> None:  # noqa: N802 (Qt API)
        from repeatviewer.widgets.cyto_images_ui import image_paths_from_mime

        if event.source() is not None or not image_paths_from_mime(event.mimeData()):
            event.ignore()
            return
        row = self._fasta_row_at(event.position().toPoint())
        if row is None:
            self._cyto_drop_row = None
            event.ignore()
            return
        if row != self._cyto_drop_row:
            self._cyto_drop_row = row
            self.selectRow(row)
        event.acceptProposedAction()

    def dragLeaveEvent(self, event) -> None:  # noqa: N802 (Qt API)
        self._cyto_drop_row = None
        super().dragLeaveEvent(event)

    def dropEvent(self, event) -> None:  # noqa: N802 (Qt API)
        from repeatviewer.widgets.cyto_images_ui import image_paths_from_mime

        self._cyto_drop_row = None
        paths = image_paths_from_mime(event.mimeData()) if event.source() is None else []
        row = self._fasta_row_at(event.position().toPoint())
        if not paths or row is None:
            event.ignore()
            return
        event.acceptProposedAction()
        fasta = Path(str(self._path_for_row(row)))
        # Open the dialog after the drop returns (Finder waits for the drop otherwise).
        QTimer.singleShot(0, lambda: self.attach_cyto_image_requested.emit(fasta, list(paths)))

    def _row_decoration(self, row: int) -> tuple[QBrush | None, QBrush | None, QFont | None]:
        """Background/foreground/font already applied to this row (if any)."""
        for column in (_COL_NAME, _COL_CHECK, *_DEFAULT_COLUMN_WIDTHS.keys()):
            cell = super().item(row, column)
            if cell is None:
                continue
            bg = cell.background()
            fg = cell.foreground()
            font = cell.font()
            if (
                (bg.style() != Qt.BrushStyle.NoBrush and bg.color().alpha() > 0)
                or (fg.style() != Qt.BrushStyle.NoBrush and fg.color().alpha() > 0)
                or font.bold()
            ):
                return bg, fg, font
        return None, None, None

    def _apply_metadata_cells(
        self,
        row: int,
        path_text: str,
        meta: BrowserFileMetadata,
    ) -> None:
        bg, fg, font = self._row_decoration(row)
        for column, text, sort_key in (
            (_COL_LENGTH, meta.length_text, meta.sort_length),
            (_COL_SEQUENCES, meta.sequence_count_text, meta.sort_sequence_count),
            (_COL_GC, meta.gc_text, meta.sort_gc),
            (_COL_TYPE, meta.molecule_type, meta.molecule_type.lower()),
            (_COL_DESCRIPTION, meta.description_text, meta.sort_description),
            (_COL_MODIFIED, meta.modified_text, meta.sort_modified),
            (_COL_SIZE, meta.size_text, meta.sort_size),
        ):
            cell = self.item(row, column)
            if cell is None:
                cell = _SortableTableItem(text, sort_key=sort_key)
                cell.setData(_ROLE_PATH, path_text)
                cell.setTextAlignment(_LEFT_ALIGN)
                self.setItem(row, column, cell)
            else:
                cell.setText(text)
                cell.setData(_ROLE_SORT, sort_key)
                cell.setData(_ROLE_PATH, path_text)
            if bg is not None:
                cell.setBackground(bg)
            if fg is not None:
                cell.setForeground(fg)
            if font is not None:
                cell.setFont(font)

    def update_rows_metadata(
        self,
        items: list[tuple[Path, BrowserFileMetadata]],
    ) -> None:
        """Apply several metadata updates without re-sorting between rows."""
        if not items:
            return
        row_index = self._path_row_index()
        frozen = self._metadata_sort_frozen
        was_sorting = self.isSortingEnabled()
        if not frozen:
            self.setSortingEnabled(False)
        self.setUpdatesEnabled(False)
        try:
            for path, meta in items:
                # listed_paths / ROLE_PATH are already absolute — avoid resolve().
                key = str(path)
                row = row_index.get(key)
                if row is None:
                    try:
                        key = str(path.expanduser().resolve())
                    except OSError:
                        continue
                    row = row_index.get(key)
                if row is None:
                    continue
                self._apply_metadata_cells(row, key, meta)
        finally:
            self.setUpdatesEnabled(True)
            # Keep sorting off while a large-folder fill is in progress.
            if not frozen:
                self.setSortingEnabled(was_sorting)

    def update_row_metadata(self, path: Path, meta: BrowserFileMetadata) -> bool:
        """Refresh metadata columns for an existing row. Returns True if found."""
        key = str(path.expanduser().resolve())
        row = self._path_row_index().get(key)
        if row is None:
            return False
        self._apply_metadata_cells(row, key, meta)
        return True

    def rename_listed_path(
        self,
        old_path: Path,
        new_path: Path,
        *,
        metadata: BrowserFileMetadata | None = None,
        fill_metadata_columns: bool | None = None,
    ) -> bool:
        """Update an existing row after a same-folder rename (no remove/add).

        Returns True when a matching row was found and updated.
        """
        old_key = str(old_path.expanduser().resolve())
        new_resolved = new_path.expanduser().resolve()
        new_key = str(new_resolved)
        row_index = self._path_row_index()
        row = row_index.get(old_key)
        if row is None:
            return False

        was_sorting = self.isSortingEnabled()
        if was_sorting:
            self.setSortingEnabled(False)
        try:
            for column in range(self.columnCount()):
                cell = self.item(row, column)
                if cell is not None:
                    cell.setData(_ROLE_PATH, new_key)

            name_item = self.item(row, _COL_NAME)
            if name_item is not None:
                name_item.setText(new_resolved.name)
                name_item.setData(_ROLE_SORT, new_resolved.name.lower())
                name_item.setIcon(browser_file_icon(self, new_resolved))
                name_item.setToolTip(new_resolved.name)

            self._apply_file_color_cells(row, new_resolved, resolved=True)
            self._apply_cyto_cell(row, new_resolved)

            if fill_metadata_columns is None:
                from repeatviewer.widgets.file_list_columns import metadata_columns_visible

                fill_metadata_columns = metadata_columns_visible()
            if fill_metadata_columns:
                if metadata is None:
                    from repeatviewer.browser_file_metadata import quick_metadata_for_file

                    metadata = quick_metadata_for_file(new_resolved)
                self._apply_metadata_cells(row, new_key, metadata)
        finally:
            if was_sorting:
                self.setSortingEnabled(True)
        return True

    def begin_metadata_fill(self) -> None:
        """Freeze sorting for the duration of deferred metadata updates."""
        if self._metadata_sort_frozen:
            return
        self._metadata_sort_frozen = True
        self._metadata_saved_scroll = int(self.verticalScrollBar().value())
        self.setSortingEnabled(False)

    def end_metadata_fill(self, *, resort: bool = True) -> None:
        """Unfreeze sorting; optionally sort once and restore scroll."""
        if not self._metadata_sort_frozen:
            return
        self._metadata_sort_frozen = False
        scroll = self._metadata_saved_scroll
        self._metadata_saved_scroll = None
        if resort and self.rowCount() > 0:
            self.setSortingEnabled(True)
            self.apply_default_sort()
        else:
            self.setSortingEnabled(True)
        if scroll is not None:
            vsb = self.verticalScrollBar()
            vsb.setValue(min(int(scroll), vsb.maximum()))

    def update_rows_metadata(
        self,
        items: list[tuple[Path, BrowserFileMetadata]],
    ) -> None:
        """Apply several metadata updates without re-sorting between rows."""
        if not items:
            return
        row_index = self._path_row_index()
        frozen = self._metadata_sort_frozen
        was_sorting = self.isSortingEnabled()
        if not frozen:
            self.setSortingEnabled(False)
        self.setUpdatesEnabled(False)
        try:
            for path, meta in items:
                key = str(path.expanduser().resolve())
                row = row_index.get(key)
                if row is None:
                    continue
                self._apply_metadata_cells(row, key, meta)
        finally:
            self.setUpdatesEnabled(True)
            # Keep sorting off while a large-folder fill is in progress.
            if not frozen:
                self.setSortingEnabled(was_sorting)

    def _path_row_index(self) -> dict[str, int]:
        index: dict[str, int] = {}
        for row in range(self.rowCount()):
            key = self._path_for_row(row)
            if key:
                index[key] = row
        return index

    def apply_default_sort(self) -> None:
        order = Qt.SortOrder.AscendingOrder if self._sort_ascending else Qt.SortOrder.DescendingOrder
        self.sortItems(self._sort_column, order)
        self.apply_list_banding()

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
            QItemSelectionModel.SelectionFlag.Select
            | QItemSelectionModel.SelectionFlag.Rows
        )
        if selected:
            self.selectionModel().select(index, flags)
        else:
            self.selectionModel().select(index, QItemSelectionModel.SelectionFlag.Deselect)

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

    def listed_paths(self) -> list[Path]:
        paths: list[Path] = []
        for row in range(self.rowCount()):
            raw = self._path_for_row(row)
            if raw:
                paths.append(Path(raw))
        return paths

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

    def suppress_selection_preview(self) -> bool:
        """True while a right-click context menu is active (must not swap Preview)."""
        return self._suppress_selection_preview

    def set_suppress_selection_preview(self, suppressed: bool) -> None:
        self._suppress_selection_preview = bool(suppressed)

    def clear_suppress_selection_preview(self) -> None:
        self._suppress_selection_preview = False

    def _clear_suppress_selection_preview(self) -> None:
        self.clear_suppress_selection_preview()

    def mousePressEvent(self, event: QMouseEvent) -> None:  # noqa: N802
        # Qt changes the current row on right-click *before* customContextMenuRequested.
        # Block Preview reload until the context menu finishes (see main_window handler).
        if event.button() == Qt.MouseButton.RightButton:
            self._suppress_selection_preview = True
        elif event.button() == Qt.MouseButton.LeftButton:
            # Never leave click-to-preview stuck after a modal dialog from a
            # context-menu action (e.g. Add annotation…).
            self._suppress_selection_preview = False
        super().mousePressEvent(event)

    def _on_header_clicked(self, column: int) -> None:
        if column == _COL_CHECK:
            return
        if self._sort_column == column:
            self._sort_ascending = not self._sort_ascending
        else:
            self._sort_column = column
            self._sort_ascending = True
        # Family banding only when Name order groups related files; otherwise
        # use simple even/odd row zebra.
        want_family = column == _COL_NAME
        self.set_list_banding(family=want_family, row_zebra=not want_family)
        self.apply_default_sort()
        window = self.window()
        if hasattr(window, "_refresh_file_list_highlights"):
            window._refresh_file_list_highlights()
        # Keep deferred metadata moving after a large-folder sort.
        if hasattr(window, "_ensure_file_list_metadata_fill"):
            window._ensure_file_list_metadata_fill()

    def _on_cell_clicked(self, row: int, column: int) -> None:
        if column == _COL_CHECK:
            return
        if self._skip_click_preview_once:
            # Mouse-up after drag must not open the dragged file in Preview.
            self._skip_click_preview_once = False
            self._suppress_selection_preview = False
            return
        if column == _COL_IMAGES and self.cyto_image_count(row) > 0:
            raw = self._path_for_row(row)
            if raw:
                self.cyto_images_requested.emit(Path(raw))
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
        self._suppress_selection_preview = True
        try:
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
                    menu.addAction(
                        "Open in Preview",
                        lambda: self.open_in_preview_requested.emit(paths),
                    )
                    menu.addAction(
                        "Open in Additional Tab",
                        lambda: self.open_in_preview_tab_requested.emit(paths),
                    )
                    menu.addAction(
                        "Open in Text Editor",
                        lambda: self.open_in_text_editor_requested.emit(paths),
                    )
                    menu.addAction("Rename", lambda: self.rename_file_requested.emit(path))
                    if path.suffix.lower() in _FASTA_SUFFIXES:
                        menu.addAction(
                            "Add as BLAST Database…",
                            lambda p=path: self.add_blast_database_requested.emit(p),
                        )
                        menu.addAction(
                            "Add Annotation…",
                            lambda p=path: self.add_annotation_requested.emit(p),
                        )
                        from repeatviewer.views.export_annotation_dialog import (
                            fasta_has_exportable_annotations,
                        )

                        export_act = menu.addAction(
                            "Export Annotation…",
                            lambda p=path: self.export_annotation_requested.emit(p),
                        )
                        export_act.setEnabled(fasta_has_exportable_annotations(path))
                        menu.addAction(
                            "Attach Cytogenetic Image…",
                            lambda p=path: self.attach_cyto_image_requested.emit(p, []),
                        )
                        try:
                            n_images = image_count_for_sequence(path)
                        except Exception:  # noqa: BLE001
                            n_images = 0
                        if n_images:
                            menu.addAction(
                                f"Show Cytogenetic Images ({n_images})",
                                lambda p=path: self.cyto_images_requested.emit(p),
                            )
                        try:
                            seq_n = count_fasta_sequences(path)
                        except Exception:  # noqa: BLE001
                            seq_n = 0
                        if seq_n > 1:
                            menu.addAction(
                                f"Separate Sequence into {seq_n} Files",
                                lambda: self.separate_fasta_requested.emit(path),
                            )
                            menu.addAction(
                                f"Separate Sequence into {seq_n} Folders",
                                lambda: self.separate_fasta_subfolders_requested.emit(path),
                            )
                        extract_act = menu.addAction(
                            "Extract Sequences from File",
                            lambda p=path: self.extract_sequences_requested.emit(p),
                        )
                        extract_act.setEnabled(seq_n > 1)
                        blast_txt = find_blast_table_for_query(path)
                        if blast_txt is not None:
                            menu.addAction(
                                "TE Density Plot",
                                lambda bt=blast_txt: self.blast_genome_map_requested.emit(bt),
                            )
                    if is_blast_table_file(path):
                        menu.addAction(
                            "Extract Subject Sequences",
                            lambda bt=path: self.extract_hit_subjects_requested.emit(bt),
                        )
                    if is_te_hit_table_file(path):
                        menu.addAction(
                            "TE Density Plot",
                            lambda bt=path: self.blast_genome_map_requested.emit(bt),
                        )
                        menu.addAction(
                            "TE Divergence Plot",
                            lambda bt=path: self.blast_divergence_requested.emit(bt),
                        )
                else:
                    menu.addAction(
                        "Open in Preview",
                        lambda: self.open_in_preview_requested.emit([paths[0]]),
                    )
                    menu.addAction(
                        f"Open {len(paths)} Files in Preview",
                        lambda: self.open_in_preview_tab_requested.emit(paths),
                    )
                    menu.addAction(
                        f"Open {len(paths)} Files in Text Editor",
                        lambda: self.open_in_text_editor_requested.emit(paths),
                    )
                    menu.addAction(
                        "Batch Rename…",
                        lambda ps=list(paths): self.batch_rename_requested.emit(ps),
                    )
                menu.addSeparator()
                if len(paths) == 1:
                    new_folder_label = "New Folder with Selected Items"
                else:
                    new_folder_label = f"New Folder with {len(paths)} Selected Items"
                menu.addAction(
                    new_folder_label,
                    lambda ps=list(paths): self.new_folder_with_files_requested.emit(ps),
                )
                menu.addSeparator()
                menu.addAction("Cut", lambda: self.cut_files_requested.emit(paths))
                menu.addAction("Copy", lambda: self.copy_files_requested.emit(paths))
                if len(paths) == 1:
                    menu.addAction(
                        "Copy Path",
                        lambda ps=list(paths): self._copy_file_paths_to_clipboard(ps),
                    )
                else:
                    menu.addAction(
                        f"Copy {len(paths)} Paths",
                        lambda ps=list(paths): self._copy_file_paths_to_clipboard(ps),
                    )
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
                        f"Duplicate {len(paths)} Files",
                        lambda: self.duplicate_files_requested.emit(paths),
                    )
                menu.addSeparator()
                menu.addAction("Open in AliView", self._context_open_aliview)
                if len(paths) == 1:
                    label = (
                        "Reveal in Finder"
                        if platform.system() == "Darwin"
                        else "Show in File Manager"
                    )
                    menu.addAction(label, lambda: self._reveal_paths(paths))
                menu.addSeparator()
                self._add_file_color_menu(menu, paths)
                menu.addSeparator()
                if len(paths) == 1:
                    menu.addAction("Delete", lambda: self.delete_files_requested.emit(paths))
                else:
                    menu.addAction(
                        f"Delete {len(paths)} Files",
                        lambda: self.delete_files_requested.emit(paths),
                    )
            elif can_paste:
                menu.addAction("Paste", lambda: self.paste_files_requested.emit())

            # Grouping needs several FASTA files, so it belongs on a multi-selection.
            fasta_paths = [p for p in paths if p.suffix.lower() in _FASTA_SUFFIXES]
            if len(fasta_paths) > 1:
                menu.addSeparator()
                menu.addAction(
                    f"Combine {len(fasta_paths)} Sequences to One File",
                    lambda ps=list(fasta_paths): self.group_sequences_requested.emit(ps),
                )
            menu.exec(self.mapToGlobal(position))
        finally:
            # Keep suppress through this event-loop turn (deferred preview opens).
            QTimer.singleShot(0, self._clear_suppress_selection_preview)

    def _add_file_color_menu(self, menu: QMenu, paths: list[Path]) -> None:
        from repeatviewer.widgets.file_colors import FILE_COLOR_PRESETS

        color_menu = menu.addMenu("File Color")
        color_menu.addAction(
            "Choose Color…",
            lambda ps=list(paths): self._choose_files_color(ps),
        )
        for label, hex_color in FILE_COLOR_PRESETS:
            color_menu.addAction(
                label,
                lambda checked=False, c=hex_color, ps=list(paths): self._set_files_color(
                    ps, c
                ),
            )
        color_menu.addSeparator()
        color_menu.addAction(
            "Clear Color",
            lambda ps=list(paths): self._clear_files_color(ps),
        )

    def _set_files_color(self, paths: list[Path], color: str) -> None:
        from repeatviewer.widgets.file_colors import set_file_color

        row_index = self._path_row_index()
        was_sorting = self.isSortingEnabled()
        self.setSortingEnabled(False)
        try:
            for path in paths:
                resolved = path.expanduser().resolve()
                set_file_color(resolved, color)
                row = row_index.get(str(resolved))
                if row is not None:
                    self._apply_file_color_cells(row, resolved)
        finally:
            self.setSortingEnabled(was_sorting)
        if self._sort_column == _COL_COLOR:
            self.apply_default_sort()

    def _choose_files_color(self, paths: list[Path]) -> None:
        from repeatviewer.widgets.file_colors import (
            FILE_COLOR_PRESETS,
            get_file_color,
        )

        if not paths:
            return
        stored = get_file_color(paths[0])
        initial = QColor(stored) if stored else QColor(FILE_COLOR_PRESETS[0][1])
        if not initial.isValid():
            initial = QColor(FILE_COLOR_PRESETS[0][1])

        dialog = QColorDialog(initial, self.window())
        dialog.setWindowTitle("File Color")
        dialog.setStyleSheet("")
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        picked = dialog.selectedColor()
        if not picked.isValid():
            return
        self._set_files_color(paths, picked.name())

    def _clear_files_color(self, paths: list[Path]) -> None:
        from repeatviewer.widgets.file_colors import set_file_color

        row_index = self._path_row_index()
        was_sorting = self.isSortingEnabled()
        self.setSortingEnabled(False)
        try:
            for path in paths:
                resolved = path.expanduser().resolve()
                set_file_color(resolved, None)
                row = row_index.get(str(resolved))
                if row is not None:
                    self._apply_file_color_cells(row, resolved)
        finally:
            self.setSortingEnabled(was_sorting)
        if self._sort_column == _COL_COLOR:
            self.apply_default_sort()

    def _context_open_aliview(self) -> None:
        window = self.window()
        if hasattr(window, "_open_selected_in_aliview"):
            window._open_selected_in_aliview()

    def _copy_file_paths_to_clipboard(self, paths: list[Path]) -> None:
        """Copy absolute file path(s) to the system clipboard (one per line)."""
        lines: list[str] = []
        for path in paths:
            try:
                lines.append(str(Path(path).expanduser().resolve()))
            except OSError:
                lines.append(str(Path(path).expanduser()))
        if not lines:
            return
        QApplication.clipboard().setText("\n".join(lines))
        window = self.window()
        status = getattr(window, "_set_status_message", None)
        if callable(status):
            if len(lines) == 1:
                status("File path copied", 3000)
            else:
                status(f"{len(lines)} file paths copied", 3000)

    def _reveal_paths(self, paths: list[Path]) -> None:
        window = self.window()
        if hasattr(window, "_reveal_in_file_manager"):
            window._reveal_in_file_manager(paths)

    def mimeData(self, items: list[FileListRow]) -> QMimeData:  # noqa: N802 (Qt API)
        mime = QMimeData()
        urls = []
        names: list[str] = []
        for row_item in items:
            raw = row_item.data(_ROLE_PATH)
            if raw:
                path = Path(str(raw))
                urls.append(path.as_uri())
                names.append(path.name)
        if urls:
            from PySide6.QtCore import QUrl

            mime.setUrls([QUrl(u) for u in urls])
        # Prefer filenames for any text-based drag preview (not file:// URIs).
        if names:
            mime.setText("\n".join(names))
        return mime

    def _drag_filename_pixmap(self, items: list[FileListRow]) -> QPixmap:
        """Compact drag ghost: short elided name, trails right of the cursor.

        Long FASTA names must not cover Folders drop targets under the pointer.
        Built at the widget device-pixel ratio so Retina/HiDPI stays sharp.
        """
        names: list[str] = []
        for row_item in items:
            raw = row_item.data(_ROLE_PATH)
            if raw:
                names.append(Path(str(raw)).name)
        if not names:
            return QPixmap()
        if len(names) == 1:
            text = names[0]
        else:
            text = f"{len(names)} files"
        font = self.font()
        font.setPointSize(max(9, font.pointSize() - 1))
        fm = QFontMetrics(font)
        # Keep the badge narrow so Folders rows stay readable under the cursor.
        text = fm.elidedText(text, Qt.TextElideMode.ElideMiddle, 140)
        pad_x, pad_y = 8, 4
        width = max(36, fm.horizontalAdvance(text) + pad_x * 2)
        height = max(20, fm.height() + pad_y * 2)
        dpr = max(1.0, float(self.devicePixelRatioF()))
        pixmap = QPixmap(int(round(width * dpr)), int(round(height * dpr)))
        pixmap.setDevicePixelRatio(dpr)
        pixmap.fill(QColor(0, 0, 0, 0))
        painter = QPainter(pixmap)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        painter.setRenderHint(QPainter.RenderHint.TextAntialiasing, True)
        # Solid white so the short label stays readable over Folders/Files.
        painter.setBrush(QColor(255, 255, 255))
        painter.setPen(QColor(160, 160, 160))
        painter.drawRoundedRect(0, 0, width - 1, height - 1, 5, 5)
        painter.setPen(self.palette().color(self.foregroundRole()))
        painter.setFont(font)
        # Integer baseline avoids sub-pixel blur on the text.
        painter.drawText(pad_x, pad_y + fm.ascent(), text)
        painter.end()
        return pixmap

    def startDrag(self, supportedActions):  # noqa: N802 (Qt API)
        items = self.selectedItems()
        if not items:
            return
        # Selecting a row to drag must not open / switch Preview.
        self._suppress_selection_preview = True
        self._skip_click_preview_once = True
        self.dragStarted.emit()
        drag = QDrag(self)
        drag.setMimeData(self.mimeData(items))
        pixmap = self._drag_filename_pixmap(items)
        if not pixmap.isNull():
            drag.setPixmap(pixmap)
            # Anchor at the left so the label sits to the right of the cursor
            # and does not cover the Folders row being hovered.
            drag.setHotSpot(pixmap.rect().topLeft())
        try:
            drag.exec(supportedActions)
        finally:
            # Delay past the mouse-release cellClicked that Qt fires after a drag.
            QTimer.singleShot(150, self._finish_drag_preview_guard)

    def _finish_drag_preview_guard(self) -> None:
        """Clear suppress if no post-drag click arrived (e.g. Escape cancel)."""
        if self._skip_click_preview_once:
            self._skip_click_preview_once = False
        self._suppress_selection_preview = False
