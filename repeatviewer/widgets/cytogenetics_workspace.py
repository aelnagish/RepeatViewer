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
    QPushButton,
    QScrollArea,
    QSizePolicy,
    QSplitter,
    QTabWidget,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends.cytogenetics import workspace as ws
from repeatviewer.widgets.cytogenetics_figure import FigureTab
from repeatviewer.widgets.cytogenetics_idiogram import IdiogramTab
from repeatviewer.widgets.cytogenetics_image import ImageTab
from repeatviewer.widgets.cytogenetics_library import LibraryTab
from repeatviewer.widgets.cytogenetics_signals import SignalsTab

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

        from repeatviewer.widgets.cytogenetics_style import workspace_css

        self.setObjectName("cytoWorkspace")
        self.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)
        self.setStyleSheet(workspace_css())
        self._path: Path | None = None
        self._folder: Path | None = None
        self._sequence_filter = ""
        self._current_image = ""

        self._title = QLabel(self)
        self._title.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        font = self._title.font()
        font.setBold(True)
        font.setPointSizeF(font.pointSizeF() + 1)
        self._title.setFont(font)

        self._side_btn = QToolButton(self)
        self._side_btn.setText("Options")
        self._side_btn.setCheckable(True)
        self._side_btn.setToolTip("Show or hide the options panel")
        self._side_btn.toggled.connect(self._on_side_toggled)
        self._help = make_help_button(self, "Cytogenetics", HELP_TEXT)

        self._crumb = QLabel(self)
        self._crumb.setObjectName("cytoCrumb")
        self._crumb.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        self._import_btn = QPushButton("Import slide session", self)
        self._import_btn.setObjectName("cytoPrimary")
        self._import_btn.clicked.connect(lambda: self._pages["Library"].start_import())

        chrome = QWidget(self)
        chrome.setObjectName("cytoChrome")
        chrome.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)
        header = QHBoxLayout(chrome)
        header.setContentsMargins(16, 10, 16, 10)
        header.setSpacing(12)
        header.addWidget(self._title)
        header.addWidget(self._crumb, 1)
        header.addWidget(self._import_btn)
        header.addWidget(self._help)
        self._side_btn.hide()  # the mockup has no options panel; Data counts stay reachable from the code

        self._tabs = QTabWidget(self)
        self._tabs.setDocumentMode(True)
        self._pages: dict[str, QWidget] = {}
        for name in ws.TAB_NAMES:
            page = {"Library": LibraryTab, "Image": ImageTab, "Signals": SignalsTab, "Idiogram": IdiogramTab, "Figure": FigureTab}.get(name, None)
            page = page(self) if page else _TabPage(name, self)
            self._pages[name] = page
            self._tabs.addTab(page, name)
        self._tabs.currentChanged.connect(self._on_tab_changed)
        self._pages["Library"].imageOpened.connect(self.open_image)
        self._pages["Library"].signalsRequested.connect(self.open_signals)
        self._pages["Library"].imageSelected.connect(self._on_image_shown)
        self._pages["Image"].imageChanged.connect(lambda i: self._on_image_shown(i, "Image"))
        self._pages["Signals"].imageChanged.connect(lambda i: self._on_image_shown(i, "Signals"))
        self._pages["Signals"].idiogramRequested.connect(lambda: self.show_tab("Idiogram"))
        self._pages["Signals"].karyotypesChanged.connect(self.refresh)
        self._pages["Idiogram"].karyotypesChanged.connect(self.refresh)
        self._pages["Idiogram"].figureRequested.connect(self.add_idiogram_to_figure)
        self._pages["Library"].figureRequested.connect(self.add_to_figure)
        self._pages["Figure"].crumbChanged.connect(self._on_figure_crumb)

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
        root.addWidget(chrome)
        root.addWidget(self._split, 1)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

        self._side_btn.setChecked(False)
        self._on_side_toggled(self._side_btn.isChecked())
        self._tabs.setCurrentIndex(self._tabs.indexOf(self._pages[ws.normalize_tab(_read("tab", ""))]))
        self._import_btn.setVisible(self.tab_name() == "Library")
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
        _write("side_panel_v2", bool(on))

    # ---------------------------------------------------------------- loading

    def load(self, path: Path) -> None:
        """Open the workspace of the folder that ``path`` (a ``.rvcyto`` marker) belongs to."""
        self._path = Path(path)
        self._folder = ws.workspace_folder(self._path)
        self._title.setText("Cytogenetics")
        self._title.setToolTip(str(self._folder))
        self._crumb.setText(ws.breadcrumb(self._folder, self._current_image))
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

    def image_tab(self) -> ImageTab:
        return self._pages["Image"]

    def open_image(self, image_id: str) -> None:
        """Show one image in the Image tab."""
        self.show_tab("Image")
        self._pages["Image"].show_image(image_id)

    def open_signals(self, image_id: str) -> None:
        """Show one image in the Signals tab."""
        self.show_tab("Signals")
        self._pages["Signals"].show_image(image_id)

    def add_to_figure(self, image_id: str) -> None:
        """Add the channels and the merge of one cell to the current figure."""
        if not image_id:
            return
        self.show_tab("Figure")
        self._pages["Figure"].add_cell(image_id)

    def add_idiogram_to_figure(self) -> None:
        self.show_tab("Figure")
        self._pages["Figure"].add_idiogram()

    def _on_figure_crumb(self, text: str) -> None:
        if self.tab_name() == "Figure" and text:
            self._crumb.setText(text)

    def _sync_crumb(self) -> None:
        if self._folder is None:
            return
        name = self.tab_name()
        if name == "Figure":
            self._crumb.setText(self._pages["Figure"].crumb())
        elif name == "Idiogram":
            self._crumb.setText(f"{self._folder.name} › all slides")
        else:
            self._crumb.setText(ws.breadcrumb(self._folder, self._current_image))

    def _on_image_shown(self, image_id: str, source: str = "Library") -> None:
        """Keep the header and the other image tabs on the same cell."""
        self._current_image = image_id
        self._pages["Figure"].current_image = image_id
        if self._folder is not None:
            self._sync_crumb()
        for name in ("Image", "Signals"):
            if name != source:
                self._pages[name].follow(image_id)

    def library(self) -> LibraryTab:
        return self._pages["Library"]

    def _on_tab_changed(self, index: int) -> None:
        name = self.tab_name()
        _write("tab", name)
        self._import_btn.setVisible(name == "Library")
        self._sync_crumb()
        self.tabChanged.emit(name)
