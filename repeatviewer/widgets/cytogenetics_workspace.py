"""Cytogenetics workspace page (opens from Tools > Cytogenetics).

Document tabs along the top (Library, Image, Signals, Idiogram, Figure) and the options
in collapsible sections on the right, laid out like the Mapping contig view. The tabs
are filled in step by step; each tab is a ``_TabPage`` that later steps replace.
"""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QSettings, Qt, Signal
from PySide6.QtWidgets import (
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QScrollArea,
    QSizePolicy,
    QSplitter,
    QTabWidget,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends.cytogenetics import workspace as ws
from repeatviewer.widgets.cytogenetics_library import LibraryTab

_SETTINGS = "cytogenetics_workspace"
SIDE_W = 270

HELP_TEXT = (
    "<b>Cytogenetics</b><br>"
    "Images of this folder (FISH, GISH and other stainings) with their slides, "
    "channels, marked chromosomes and signals.<br><br>"
    "<b>Library</b>: all images, grouped by slide, with filters.<br>"
    "<b>Image</b>: one cell with its channels, levels and scale bar.<br>"
    "<b>Signals</b>: number chromosomes and mark signals.<br>"
    "<b>Idiogram</b>: signals on the karyotype, compared with the assembly.<br>"
    "<b>Figure</b>: panels with letters, labels and a legend.<br><br>"
    "Use <b>Options</b> to show or hide the side panel."
)


def _settings() -> QSettings:
    return QSettings()


def _read(key: str, default, kind=str):
    try:
        value = _settings().value(f"{_SETTINGS}/{key}", default)
        if kind is bool:
            return value if isinstance(value, bool) else str(value).lower() in ("1", "true", "yes")
        return kind(value)
    except (TypeError, ValueError):
        return default


def _write(key: str, value) -> None:
    try:
        _settings().setValue(f"{_SETTINGS}/{key}", value)
    except Exception:  # noqa: BLE001 - a remembered setting is never worth an error
        pass


class _TabPage(QWidget):
    """One workspace tab. Later steps replace the placeholder body with the real view."""

    def __init__(self, name: str, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.name = name
        lay = QVBoxLayout(self)
        lay.setContentsMargins(16, 16, 16, 16)
        self._hint = QLabel(self)
        self._hint.setWordWrap(True)
        self._hint.setAlignment(Qt.AlignmentFlag.AlignTop | Qt.AlignmentFlag.AlignLeft)
        self._hint.setText(f"The {name} tab is not available yet.")
        lay.addWidget(self._hint, 1)

    def refresh(self, folder: Path | None) -> None:
        """Reload from the folder's store (called when the images change)."""


class CytogeneticsWorkspace(QWidget):
    """The workspace of one data folder (see ``backends.cytogenetics.workspace``)."""

    openRequested = Signal(object)  # Path of a file to open in the preview
    fileWritten = Signal(object)    # Path of a file written from the workspace
    tabChanged = Signal(str)

    def __init__(self, path: Path | None = None, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        from repeatviewer.widgets.cyto_images_ui import cyto_notifier
        from repeatviewer.widgets.help_button import make_help_button

        self.setAutoFillBackground(True)
        self.setBackgroundRole(self.palette().ColorRole.Window)
        self._path: Path | None = None
        self._folder: Path | None = None
        self._sequence_filter = ""

        self._title = QLabel(self)
        self._title.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        font = self._title.font()
        font.setBold(True)
        self._title.setFont(font)

        self._side_btn = QToolButton(self)
        self._side_btn.setText("Options")
        self._side_btn.setCheckable(True)
        self._side_btn.setToolTip("Show or hide the options panel")
        self._side_btn.toggled.connect(self._on_side_toggled)
        self._help = make_help_button(self, "Cytogenetics", HELP_TEXT)

        header = QHBoxLayout()
        header.setContentsMargins(10, 6, 10, 2)
        header.addWidget(self._title, 1)
        header.addWidget(self._side_btn)
        header.addWidget(self._help)

        self._tabs = QTabWidget(self)
        self._tabs.setDocumentMode(True)
        self._pages: dict[str, QWidget] = {}
        for name in ws.TAB_NAMES:
            page = LibraryTab(self) if name == "Library" else _TabPage(name, self)
            self._pages[name] = page
            self._tabs.addTab(page, name)
        self._tabs.currentChanged.connect(self._on_tab_changed)

        self._data_form = QFormLayout()
        self._side = self._build_side()
        self._side.setMinimumWidth(SIDE_W)

        self._split = QSplitter(Qt.Orientation.Horizontal, self)
        self._split.setChildrenCollapsible(False)
        self._split.addWidget(self._tabs)
        self._split.addWidget(self._side)
        self._split.setStretchFactor(0, 1)
        self._split.setStretchFactor(1, 0)
        self._split.setSizes([900, SIDE_W])

        root = QVBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)
        root.addLayout(header)
        root.addWidget(self._split, 1)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

        self._side_btn.setChecked(_read("side_panel", True, bool))
        self._on_side_toggled(self._side_btn.isChecked())
        self._tabs.setCurrentIndex(self._tabs.indexOf(self._pages[ws.normalize_tab(_read("tab", ""))]))
        cyto_notifier().changed.connect(self._on_store_changed)
        if path is not None:
            self.load(path)

    # ---------------------------------------------------------------- side panel

    def _build_side(self) -> QWidget:
        from repeatviewer.widgets.sequence_annotation_rail import _CollapsibleSection

        body = QWidget(self)
        col = QVBoxLayout(body)
        col.setContentsMargins(10, 8, 10, 8)
        col.setSpacing(12)
        self._sections: dict[str, QWidget] = {}

        def section(key: str, title: str, default: bool = True):
            s = _CollapsibleSection(title, body, expanded=_read(f"section/{key}", default, bool))
            s.expandedChanged.connect(lambda on, k=key: _write(f"section/{k}", bool(on)))
            self._sections[key] = s
            col.addWidget(s)
            return s

        s = section("data", "Data")
        self._data_form.setContentsMargins(0, 0, 0, 0)
        self._data_form.setHorizontalSpacing(8)
        self._data_form.setVerticalSpacing(6)
        s.add_layout(self._data_form)
        col.addStretch(1)

        scroll = QScrollArea(self)
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QScrollArea.Shape.NoFrame)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setWidget(body)
        return scroll

    def _fill_data(self) -> None:
        while self._data_form.rowCount():
            self._data_form.removeRow(0)
        if self._folder is None:
            return
        for label, value in ws.summary_lines(ws.folder_summary(self._folder)):
            self._data_form.addRow(QLabel(label, self), QLabel(value, self))

    def _on_side_toggled(self, on: bool) -> None:
        self._side.setVisible(bool(on))
        _write("side_panel", bool(on))

    # ---------------------------------------------------------------- loading

    def load(self, path: Path) -> None:
        """Open the workspace of the folder that ``path`` (a ``.rvcyto`` marker) belongs to."""
        self._path = Path(path)
        self._folder = ws.workspace_folder(self._path)
        self._title.setText(f"Cytogenetics: {self._folder.name}")
        self._title.setToolTip(str(self._folder))
        self.refresh()

    def refresh(self) -> None:
        self._fill_data()
        for page in self._pages.values():
            page.refresh(self._folder)

    def _on_store_changed(self, folder: str) -> None:
        if self._folder is not None and Path(folder) == self._folder:
            self.refresh()

    # ---------------------------------------------------------------- tabs

    def folder(self) -> Path | None:
        return self._folder

    def tab_name(self) -> str:
        page = self._tabs.currentWidget()
        return getattr(page, "name", ws.TAB_NAMES[0])

    def show_tab(self, name: str) -> None:
        self._tabs.setCurrentWidget(self._pages[ws.normalize_tab(name)])

    def sequence_filter(self) -> str:
        """FASTA name the Library is filtered to ("" = all images)."""
        return self._sequence_filter

    def set_sequence_filter(self, fasta_name: str) -> None:
        self._sequence_filter = str(fasta_name or "")
        self._pages["Library"].set_sequence(self._sequence_filter)
        self.refresh()

    def library(self) -> LibraryTab:
        return self._pages["Library"]

    def _on_tab_changed(self, index: int) -> None:
        name = self.tab_name()
        _write("tab", name)
        self.tabChanged.emit(name)
