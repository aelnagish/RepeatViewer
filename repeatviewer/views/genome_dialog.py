"""Choose or change the project genome FASTA (File → Define genome)."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QEvent, Qt, Signal
from PySide6.QtGui import QDragEnterEvent, QDropEvent
from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QFrame,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QRadioButton,
    QScrollArea,
    QStackedWidget,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.pan_genome import (
    GENOME_MODE_PAN,
    GENOME_MODE_SINGLE,
    MAX_PAN_GENOME_NAME_LEN,
    PanGenomeEntry,
    normalize_pan_genome_name_strict,
)
from repeatviewer.project import RepeatViewerProject
from repeatviewer.views.settings_dialog import (
    MAX_GENOME_HISTORY,
    _genome_history_label,
    add_genome_to_history,
    load_genome_history,
    save_workspace,
)

_GENOME_SUFFIXES = {".fa", ".fasta", ".fna"}
_FASTA_FILTER = "FASTA (*.fa *.fasta *.fna);;All files (*)"
MIN_PAN_GENOMES = 2
_NAME_FIELD_WIDTH = 120
_SINGLE_STACK_HEIGHT = 44
_PAN_ROW_BLOCK_HEIGHT = 32
_PAN_SCROLL_BASE_HEIGHT = 68
_PAN_SCROLL_MAX_HEIGHT = 360
_PAN_LIST_EXTRA_HEIGHT = 40


def is_genome_fasta(path: Path) -> bool:
    try:
        resolved = path.expanduser().resolve()
    except OSError:
        return False
    return resolved.is_file() and resolved.suffix.lower() in _GENOME_SUFFIXES


def _paths_from_mime(mime) -> list[Path]:
    if not mime.hasUrls():
        return []
    paths: list[Path] = []
    for url in mime.urls():
        if not url.isLocalFile():
            continue
        path = Path(url.toLocalFile())
        if is_genome_fasta(path):
            paths.append(path.expanduser().resolve())
    return paths


def _default_genome_name(path: Path) -> str:
    return path.stem[:MAX_PAN_GENOME_NAME_LEN]


class GenomeFastaDropLineEdit(QLineEdit):
    """Path field that accepts FASTA files dragged from the Files panel or Finder."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setAcceptDrops(True)
        self.setPlaceholderText("Drop a genome FASTA here, browse, or pick from Files")
        self.setClearButtonEnabled(True)
        self._on_path_dropped = None

    def set_on_path_dropped(self, callback) -> None:
        self._on_path_dropped = callback

    def dragEnterEvent(self, event: QDragEnterEvent) -> None:  # noqa: N802
        if _paths_from_mime(event.mimeData()):
            event.acceptProposedAction()
            return
        event.ignore()

    def dragMoveEvent(self, event) -> None:  # noqa: N802
        if _paths_from_mime(event.mimeData()):
            event.acceptProposedAction()
            return
        event.ignore()

    def dropEvent(self, event: QDropEvent) -> None:  # noqa: N802
        paths = _paths_from_mime(event.mimeData())
        if not paths:
            event.ignore()
            return
        self.set_path(paths[0])
        event.acceptProposedAction()

    def set_path(self, path: Path | str) -> None:
        self.setText(str(path))
        if self._on_path_dropped is not None:
            self._on_path_dropped(Path(str(path)))


class _PanGenomeRow(QWidget):
    """One pan-genome line: Name, Genome FASTA, Browse, Remove."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(8)

        self.name = QLineEdit()
        self.name.setMaxLength(MAX_PAN_GENOME_NAME_LEN)
        self.name.setPlaceholderText("Name (max 10)")
        self.name.setFixedWidth(_NAME_FIELD_WIDTH)

        self.path = GenomeFastaDropLineEdit()
        self.path.set_on_path_dropped(self._maybe_fill_name_from_path)

        browse = QPushButton("Browse")
        browse.clicked.connect(self._browse)

        self._remove_btn = QPushButton("Remove")
        self._remove_btn.setToolTip("Remove this genome from the list")
        self._remove_btn.clicked.connect(self._request_remove)

        layout.addWidget(QLabel("Name:"))
        layout.addWidget(self.name)
        layout.addWidget(QLabel("Genome FASTA:"))
        layout.addWidget(self.path, stretch=1)
        layout.addWidget(browse)
        layout.addWidget(self._remove_btn)

        self._on_remove = None
        self._can_remove = True

    def set_can_remove(self, allowed: bool) -> None:
        self._can_remove = allowed
        self._remove_btn.setEnabled(allowed)

    def set_on_remove(self, callback) -> None:
        self._on_remove = callback

    def _maybe_fill_name_from_path(self, path: Path) -> None:
        if not self.name.text().strip():
            self.name.setText(_default_genome_name(path))

    def _request_remove(self) -> None:
        if self._can_remove and self._on_remove is not None:
            self._on_remove(self)

    def _browse(self) -> None:
        start = self.path.text().strip() or str(Path.home())
        chosen, _ = QFileDialog.getOpenFileName(
            self,
            "Select genome FASTA",
            start,
            _FASTA_FILTER,
        )
        if chosen:
            self.path.set_path(chosen)

    def focus_path_field(self) -> None:
        self.path.setFocus(Qt.FocusReason.OtherFocusReason)

    def has_path(self) -> bool:
        return bool(self.path.text().strip())

    def to_entry(self) -> PanGenomeEntry | None:
        name = self.name.text().strip()
        path_text = self.path.text().strip()
        if not name and not path_text:
            return None
        if not name or not path_text:
            raise ValueError("Each pan-genome entry needs both a name and a FASTA path.")
        return PanGenomeEntry(name=name, path=Path(path_text))


class PanGenomeListWidget(QWidget):
    size_changed = Signal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._rows: list[_PanGenomeRow] = []
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(8)

        self._scroll = QScrollArea()
        self._scroll.setWidgetResizable(True)
        self._scroll.setFrameShape(QFrame.Shape.NoFrame)
        self._container = QWidget()
        self._rows_layout = QVBoxLayout(self._container)
        self._rows_layout.setContentsMargins(0, 0, 0, 0)
        self._rows_layout.setSpacing(2)
        self._scroll.setWidget(self._container)
        outer.addWidget(self._scroll)

        self._add_btn = QPushButton("Add another genome")
        self._add_btn.clicked.connect(self._on_add_genome)
        add_row_layout = QHBoxLayout()
        add_row_layout.setContentsMargins(0, 0, 0, 0)
        add_row_layout.addWidget(self._add_btn)
        add_row_layout.addStretch(1)
        outer.addLayout(add_row_layout)

        self._focused_row: _PanGenomeRow | None = None
        self._path_field_focused = False

    def _on_add_genome(self) -> None:
        self.add_row(focus=True, expand=True)

    def _update_scroll_height(self) -> None:
        n = len(self._rows)
        height = _PAN_SCROLL_BASE_HEIGHT + max(0, n - 2) * _PAN_ROW_BLOCK_HEIGHT
        height = min(height, _PAN_SCROLL_MAX_HEIGHT)
        self._scroll.setMinimumHeight(height)
        self._scroll.setMaximumHeight(height)
        self.setMinimumHeight(height + _PAN_LIST_EXTRA_HEIGHT)

    def add_row(
        self,
        entry: PanGenomeEntry | None = None,
        *,
        focus: bool = False,
        expand: bool = False,
    ) -> _PanGenomeRow:
        row = _PanGenomeRow(self._container)
        row.set_on_remove(self._remove_row)
        row.path.installEventFilter(self)
        row.name.installEventFilter(self)
        if entry is not None:
            row.name.setText(entry.name)
            row.path.setText(str(entry.path))
        self._rows.append(row)
        self._rows_layout.addWidget(row)
        self._sync_remove_buttons()
        self._update_scroll_height()
        if focus:
            row.focus_path_field()
            self._focused_row = row
        if expand:
            self.size_changed.emit()
        return row

    def eventFilter(self, watched, event) -> bool:  # noqa: N802
        if event.type() in {QEvent.Type.FocusIn, QEvent.Type.MouseButtonPress}:
            for row in self._rows:
                if watched is row.path:
                    self._focused_row = row
                    self._path_field_focused = True
                    break
                if watched is row.name:
                    self._focused_row = row
                    self._path_field_focused = False
                    break
        return super().eventFilter(watched, event)

    def _sync_remove_buttons(self) -> None:
        can_remove = len(self._rows) > MIN_PAN_GENOMES
        for row in self._rows:
            row.set_can_remove(can_remove)

    def _remove_row(self, row: _PanGenomeRow) -> None:
        if len(self._rows) <= MIN_PAN_GENOMES:
            return
        if row not in self._rows:
            return
        self._rows.remove(row)
        self._rows_layout.removeWidget(row)
        row.deleteLater()
        self._sync_remove_buttons()
        self._update_scroll_height()
        self.size_changed.emit()
        if self._focused_row is row:
            self._focused_row = self._rows[-1] if self._rows else None

    def snapshot_rows(self) -> list[tuple[str, str]]:
        """Capture all row text, including partial or not-yet-validated entries."""
        return [(row.name.text(), row.path.text()) for row in self._rows]

    def restore_snapshot(self, snapshot: list[tuple[str, str]]) -> None:
        for row in list(self._rows):
            self._rows.remove(row)
            self._rows_layout.removeWidget(row)
            row.deleteLater()
        items = list(snapshot)
        while len(items) < MIN_PAN_GENOMES:
            items.append(('', ''))
        for name, path in items:
            row = self.add_row()
            row.name.setText(name)
            row.path.setText(path)
        self._focused_row = None
        self._path_field_focused = False
        self._update_scroll_height()

    def set_entries(self, entries: list[PanGenomeEntry]) -> None:
        for row in list(self._rows):
            self._rows.remove(row)
            self._rows_layout.removeWidget(row)
            row.deleteLater()
        items = list(entries)
        while len(items) < MIN_PAN_GENOMES:
            items.append(None)
        for entry in items:
            self.add_row(entry if isinstance(entry, PanGenomeEntry) else None)
        self._focused_row = None
        self._path_field_focused = False
        self._update_scroll_height()

    def apply_path(self, path: Path) -> bool:
        if not is_genome_fasta(path):
            return False
        row = self._focused_row
        if row is None or not self._path_field_focused:
            return False
        row.path.set_path(path)
        return True

    def apply_paths(self, paths: list[Path]) -> None:
        for path in paths:
            self.apply_path(path)

    def activate_panel_path_pick(self) -> None:
        """Focus the active pan-genome path field for Files-panel selection."""
        self._path_field_focused = True
        row = self._focused_row
        if row is None and self._rows:
            row = self._rows[0]
            self._focused_row = row
        if row is not None:
            row.path.setFocus(Qt.FocusReason.OtherFocusReason)

    def entries(self) -> list[PanGenomeEntry]:
        out: list[PanGenomeEntry] = []
        seen: set[str] = set()
        for row in self._rows:
            entry = row.to_entry()
            if entry is None:
                continue
            if entry.name in seen:
                raise ValueError(f"Duplicate genome name: {entry.name}")
            if not entry.path.is_file():
                raise ValueError(f"Genome file not found for '{entry.name}': {entry.path}")
            seen.add(entry.name)
            out.append(entry)
        return out


class GenomeDialog(QDialog):
    """Select single or pan-genome mode for Blast, Extend, TEAid, and Flexidot."""

    def __init__(
        self,
        project: RepeatViewerProject,
        parent: QWidget | None = None,
        *,
        selected_paths: list[Path] | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Define genome")
        self.setMinimumWidth(760)
        self._project = project
        self.setAcceptDrops(True)
        self.setWindowFlag(Qt.WindowType.WindowStaysOnTopHint, True)

        self._mode_single = QRadioButton("Single genome mode")
        self._mode_pan = QRadioButton("Pan-genome mode (Beta)")
        self._mode_single.setToolTip(
            "One reference genome for Blast, MSA Extend, TEAid, Flexidot, and Visualization."
        )
        self._mode_pan.setToolTip(
            "Several genomes with short names (max 10 characters). "
            "Blast, TEAid, and Flexidot run against each genome; Extend uses the genome you pick."
        )

        self._genome = QComboBox()
        self._genome.setEditable(True)
        self._genome.setInsertPolicy(QComboBox.InsertPolicy.NoInsert)
        self._genome.setMinimumWidth(360)
        line = self._genome.lineEdit()
        if line is not None:
            line.setPlaceholderText("Drop a FASTA, browse, or click a file in Files")
            line.setClearButtonEnabled(True)
            line.setAcceptDrops(True)
            line.installEventFilter(self)
        self._populate_genome_history(project.genome_file)
        self._genome.activated.connect(self._on_genome_history_selected)
        self._genome.installEventFilter(self)

        self._pan_list = PanGenomeListWidget()
        self._pan_list.size_changed.connect(self._adjust_dialog_size)
        self._single_genome_path_focused = False

        self._pan_snapshot: list[tuple[str, str]] = [
            (entry.name, str(entry.path)) for entry in project.pan_genomes
        ]
        self._single_genome_snapshot = (
            self._path_str(project.genome_file) if project.genome_file is not None else ''
        )

        if project.genome_mode == GENOME_MODE_PAN:
            self._mode_pan.setChecked(True)
        else:
            self._mode_single.setChecked(True)

        self._pan_list.restore_snapshot(self._pan_snapshot)
        if self._single_genome_snapshot:
            self._genome.lineEdit().setText(self._single_genome_snapshot)

        mode_row = QHBoxLayout()
        mode_row.setContentsMargins(0, 0, 0, 0)
        mode_row.setSpacing(16)
        mode_row.addWidget(self._mode_single)
        mode_row.addWidget(self._mode_pan)
        mode_row.addStretch(1)
        self._mode_bar = QWidget()
        self._mode_bar.setLayout(mode_row)

        single_page = QWidget()
        single_line = QHBoxLayout(single_page)
        single_line.setContentsMargins(0, 0, 0, 0)
        single_line.setSpacing(8)
        single_line.addWidget(QLabel("Genome FASTA:"))
        single_line.addWidget(self._genome, stretch=1)
        browse_btn = QPushButton("Browse")
        browse_btn.clicked.connect(self._browse_single_genome)
        single_line.addWidget(browse_btn)

        pan_page = QWidget()
        pan_layout = QVBoxLayout(pan_page)
        pan_layout.setContentsMargins(0, 0, 0, 0)
        pan_layout.addWidget(self._pan_list)

        self._content_stack = QStackedWidget()
        self._content_stack.addWidget(single_page)
        self._content_stack.addWidget(pan_page)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self._accept)
        buttons.rejected.connect(self.reject)

        self._hint = QLabel()
        self._hint.setWordWrap(True)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 12, 12, 12)
        layout.setSpacing(10)
        layout.addWidget(
            self._mode_bar,
            alignment=Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop,
        )
        layout.addWidget(self._content_stack)
        layout.addWidget(self._hint)
        layout.addWidget(buttons)

        self._mode_single.toggled.connect(self._sync_mode_widgets)
        self._sync_mode_widgets()

        if selected_paths and not self._project.has_genome_configured():
            self._apply_initial_paths(selected_paths)

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        self.raise_()
        self.activateWindow()
        self._adjust_dialog_size()
        self._activate_panel_path_pick()

    def _activate_panel_path_pick(self) -> None:
        """Focus the active path field so the next Files-panel click can fill it."""
        if self._mode_pan.isChecked():
            self._pan_list.activate_panel_path_pick()
        else:
            self._single_genome_path_focused = True
            line = self._genome.lineEdit()
            if line is not None:
                line.setFocus(Qt.FocusReason.OtherFocusReason)

    def _adjust_dialog_size(self) -> None:
        if self._mode_pan.isChecked():
            pan_h = self._pan_list.sizeHint().height()
            self._content_stack.setMinimumHeight(pan_h)
            self._content_stack.setMaximumHeight(_PAN_SCROLL_MAX_HEIGHT + _PAN_LIST_EXTRA_HEIGHT)
        else:
            self._content_stack.setMinimumHeight(_SINGLE_STACK_HEIGHT)
            self._content_stack.setMaximumHeight(_SINGLE_STACK_HEIGHT)
        self.adjustSize()
        hint = self.sizeHint()
        self.resize(max(self.minimumWidth(), self.width()), hint.height())

    def apply_file_from_panel(self, path: Path) -> bool:
        """Fill the focused genome path field from a Files panel click."""
        if not is_genome_fasta(path):
            return False
        resolved = path.expanduser().resolve()
        if self._mode_pan.isChecked():
            return self._pan_list.apply_path(resolved)
        if not self._single_genome_path_focused:
            return False
        line = self._genome.lineEdit()
        if line is None:
            return False
        line.setText(str(resolved))
        return True

    def _apply_initial_paths(self, paths: list[Path]) -> None:
        genomes = [p.expanduser().resolve() for p in paths if is_genome_fasta(p)]
        if not genomes:
            return
        if self._mode_pan.isChecked() or len(genomes) > 1:
            self._mode_pan.setChecked(True)
            entries = [
                PanGenomeEntry(name=_default_genome_name(path), path=path)
                for path in genomes
            ]
            self._pan_list.set_entries(entries)
            self._pan_snapshot = self._pan_list.snapshot_rows()
            self._adjust_dialog_size()
        else:
            self._mode_single.setChecked(True)
            self._genome.lineEdit().setText(str(genomes[0]))

    def _sync_mode_widgets(self) -> None:
        if self._mode_pan.isChecked():
            self._content_stack.setCurrentIndex(1)
            self._pan_list.restore_snapshot(self._pan_snapshot)
        else:
            self._pan_snapshot = self._pan_list.snapshot_rows()
            self._single_genome_snapshot = self._genome.lineEdit().text()
            self._content_stack.setCurrentIndex(0)
            if self._single_genome_snapshot:
                self._genome.lineEdit().setText(self._single_genome_snapshot)
        self._update_hint()
        self._adjust_dialog_size()
        if self.isVisible():
            self._activate_panel_path_pick()

    def _update_hint(self) -> None:
        if self._mode_pan.isChecked():
            self._hint.setText(
                "Pan-genome mode: at least two genomes, each with a unique name (≤10 characters). "
                "Click in a Genome FASTA field, then click a file in the Files panel, or drag "
                "FASTA files onto a path field or browse. "
                "Blast, TEAid, and Flexidot run against every genome; Extend uses the genome you pick. "
                "Visualization is available only in single-genome mode."
            )
        else:
            self._hint.setText(
                "Single-genome mode: click in the path field, then click a file in the Files panel, "
                "or drag a FASTA onto the path field or browse. Chromosome sizes and the BLAST "
                "index are prepared after you choose a genome."
            )

    def dragEnterEvent(self, event: QDragEnterEvent) -> None:  # noqa: N802
        if _paths_from_mime(event.mimeData()):
            event.acceptProposedAction()
            return
        event.ignore()

    def dragMoveEvent(self, event) -> None:  # noqa: N802
        if _paths_from_mime(event.mimeData()):
            event.acceptProposedAction()
            return
        event.ignore()

    def dropEvent(self, event: QDropEvent) -> None:  # noqa: N802
        paths = _paths_from_mime(event.mimeData())
        if not paths:
            event.ignore()
            return
        if self._mode_pan.isChecked():
            self._pan_list.apply_paths(paths)
        else:
            self._genome.lineEdit().setText(str(paths[0]))
        event.acceptProposedAction()

    def eventFilter(self, watched, event) -> bool:  # noqa: N802
        line = self._genome.lineEdit()
        if watched is self._genome or (line is not None and watched is line):
            if event.type() in {QEvent.Type.FocusIn, QEvent.Type.MouseButtonPress}:
                self._single_genome_path_focused = True
            if line is not None and watched is line:
                if event.type() == QEvent.Type.DragEnter:
                    if _paths_from_mime(event.mimeData()):
                        event.acceptProposedAction()
                        return True
                    event.ignore()
                    return True
                if event.type() == QEvent.Type.DragMove:
                    if _paths_from_mime(event.mimeData()):
                        event.acceptProposedAction()
                        return True
                    event.ignore()
                    return True
                if event.type() == QEvent.Type.Drop:
                    paths = _paths_from_mime(event.mimeData())
                    if paths:
                        line.setText(str(paths[0]))
                        event.acceptProposedAction()
                        return True
                    event.ignore()
                    return True
        return super().eventFilter(watched, event)

    def _path_str(self, path: Path | None) -> str:
        return str(path) if path else ""

    def _populate_genome_history(self, current: Path | None) -> None:
        history = load_genome_history()
        if current is not None:
            try:
                resolved = current.expanduser().resolve()
            except OSError:
                resolved = None
            if resolved is not None and resolved.is_file() and resolved not in history:
                history.insert(0, resolved)
        self._genome.blockSignals(True)
        self._genome.clear()
        for path in history[:MAX_GENOME_HISTORY]:
            self._genome.addItem(_genome_history_label(path), str(path))
        self._genome.blockSignals(False)
        if current is not None:
            self._genome.lineEdit().setText(self._path_str(current))
        elif self._genome.count() > 0:
            self._genome.setCurrentIndex(0)
            self._on_genome_history_selected(0)

    def _on_genome_history_selected(self, index: int) -> None:
        if index < 0:
            return
        path = self._genome.itemData(index)
        if path:
            self._genome.lineEdit().setText(str(path))

    def _current_genome_path(self) -> Path | None:
        text = self._genome.lineEdit().text().strip()
        if not text:
            return None
        return Path(text).expanduser()

    def _browse_single_genome(self) -> None:
        start = self._genome.lineEdit().text().strip() or str(Path.home())
        chosen, _ = QFileDialog.getOpenFileName(
            self,
            "Select genome FASTA",
            start,
            _FASTA_FILTER,
        )
        if chosen:
            self._genome.lineEdit().setText(chosen)

    def _accept(self) -> None:
        try:
            if self._mode_pan.isChecked():
                entries = self._pan_list.entries()
                if len(entries) < MIN_PAN_GENOMES:
                    raise ValueError(f"Pan-genome mode requires at least {MIN_PAN_GENOMES} genomes.")
                for entry in entries:
                    # Strict validation only for the definition dialog.
                    normalize_pan_genome_name_strict(entry.name)
            else:
                genome = self._current_genome_path()
                if genome is None or not genome.is_file():
                    raise ValueError("Choose a valid genome FASTA file.")
        except ValueError as exc:
            QMessageBox.warning(self, "Define genome", str(exc))
            return
        self.accept()

    def apply_to_project(self) -> RepeatViewerProject:
        if self._mode_pan.isChecked():
            entries = self._pan_list.entries()
            self._project.genome_mode = GENOME_MODE_PAN
            self._project.pan_genomes = entries
            self._project.genome_file = entries[0].path if entries else None
            for entry in entries:
                add_genome_to_history(entry.path)
        else:
            genome = self._current_genome_path()
            self._project.genome_mode = GENOME_MODE_SINGLE
            self._project.genome_file = genome
            self._pan_snapshot = self._pan_list.snapshot_rows()
            complete_pan: list[PanGenomeEntry] = []
            for name, path_text in self._pan_snapshot:
                name = name.strip()
                path_text = path_text.strip()
                if not name or not path_text:
                    continue
                try:
                    complete_pan.append(PanGenomeEntry(name=name, path=Path(path_text)))
                except ValueError:
                    continue
            if complete_pan:
                self._project.pan_genomes = complete_pan
            if genome is not None and genome.is_file():
                add_genome_to_history(genome)
        save_workspace(self._project)
        return self._project
