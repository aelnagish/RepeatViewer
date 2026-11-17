"""The Primers and Probes workspace: one tab for the primer work on a sequence.

Tools > Primers and Probes used to show its panel inside the Annotation view's
side rail, where a click on Display, Annotation, Find or Others replaced it and
nothing brought it back. The workspace is a tab of its own: a row of sections on
top (the menu item that was clicked is the one selected), the sequence it works
on, the sequence map, and the settings of the section on the right. Nothing in
it can be switched away by the rail, because the rail's tab strip is not shown.

The sequence view and the panel are the ones the Annotation view uses
(``SequenceMapPanel`` and its ``PrimersPanel``), so every primer function works
as before.

Developer: Abdullah Yusuf
"""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Qt, QTimer, Signal
from PySide6.QtWidgets import (
    QApplication,
    QComboBox,
    QHBoxLayout,
    QLabel,
    QMessageBox,
    QPushButton,
    QSplitter,
    QStackedWidget,
    QTabBar,
    QVBoxLayout,
    QWidget,
)

TAB_NAME = "Primers and Probes"

SECTION_DESIGN = "design"
SECTION_PCR = "pcr"
SECTION_DIGEST = "digest"
SECTIONS: tuple[tuple[str, str], ...] = (
    (SECTION_DESIGN, "Design"),
    (SECTION_PCR, "PCR"),
    (SECTION_DIGEST, "Digest"),
)

# Pages of the PrimersPanel behind each section. Design holds two for now: the
# design settings and the primers on the sequence.
PANEL_DESIGN, PANEL_PRIMER, PANEL_PCR, PANEL_ENZYMES = range(4)
_PANEL_PAGES: dict[str, tuple[int, ...]] = {
    SECTION_DESIGN: (PANEL_DESIGN, PANEL_PRIMER),
    SECTION_PCR: (PANEL_PCR,),
    SECTION_DIGEST: (PANEL_ENZYMES,),
}

FASTA_FILTER = (
    "FASTA (*.fa *.fasta *.fna *.fas *.fsa *.ffn *.seq *.aln *.mfa *.gz);;All files (*)"
)


SEQUENCE_HELP = (
    "<b>Sequence</b>: the file the workspace works on. Click a FASTA's name in Files to "
    "use it, or press Browse. Ticking a file's box does not change it.<br>"
    "<b>Record</b>: for a file with several sequences, the one to work on.<br>"
    "To open a file in a tab of its own instead, right-click it in Files and choose "
    "Open in Additional Tab."
)
SETTINGS_HEADER_HELP = (
    "The design settings are shown as one line while a primer is selected, to make room "
    "for its details.<br>"
    "<b>Show</b>: opens all settings again.<br>"
    "<b>Restore Defaults</b>: puts every setting back to its default.<br>"
    "<b>Design</b>: runs the design again with these settings. Primers already on the "
    "sequence are kept."
)
LIST_HELP = (
    "One list for what Design found and for the primers already on this sequence.<br>"
    "<b>candidate</b>: found by the last Design run. Click it to see it on the sequence, "
    "tick the ones to keep and press <b>Keep Ticked</b>. Nothing is written before that.<br>"
    "<b>on the sequence</b>: a primer of this sequence. Click it to see its details on "
    "the right and its binding site on the map; the tick marks the primers that PCR uses. "
    "Right-click it to edit, copy, save to the Primer Database or remove it.<br>"
    "<b>Add Primer</b>, <b>From Database</b>: put a primer on the sequence by hand, or "
    "from your saved primers.<br>"
    "<b>Open Report</b>: writes the usual primer report for the last Design run."
)


def _path_row_of(widget):
    """The PathBrowseRow a focused widget belongs to, or None."""
    from repeatviewer.widgets.path_browse import PathBrowseRow

    while widget is not None:
        if isinstance(widget, PathBrowseRow):
            return widget
        widget = widget.parentWidget()
    return None


def section_of_panel_page(index: int) -> str:
    """The section a page of the PrimersPanel belongs to."""
    for key, pages in _PANEL_PAGES.items():
        if int(index) in pages:
            return key
    return SECTION_DESIGN


def section_key(name: str | None) -> str:
    """A section key from what a caller passes ("enzymes" and "primer" are old names)."""
    text = str(name or "").strip().lower()
    if text in ("enzymes", "digest"):
        return SECTION_DIGEST
    if text == "pcr":
        return SECTION_PCR
    return SECTION_DESIGN


def same_file(a: Path | str | None, b: Path | str | None) -> bool:
    if not a or not b:
        return False
    try:
        return Path(a).expanduser().resolve() == Path(b).expanduser().resolve()
    except OSError:
        return str(a) == str(b)


class PrimersWorkspacePage(QWidget):
    """A tool tab: sections, the sequence, its map and the section's settings."""

    statusMessage = Signal(str)
    annotationsSaved = Signal(object)    # path of the annotation file just written
    sequenceExtracted = Signal(object)   # path of a FASTA written from the map
    sequenceChanged = Signal(object)     # Path of the sequence now worked on

    def __init__(self, fasta_path: Path, parent: QWidget | None = None, *,
                 section: str = SECTION_DESIGN) -> None:
        super().__init__(parent)
        from repeatviewer.widgets.path_browse import PathBrowseRow
        from repeatviewer.widgets.sequence_map_viewer import SequenceMapPanel

        # A tool page, not the preview of a file: the tabs must never treat it
        # as the page of the FASTA it works on.
        self.path = None
        self._syncing = False
        self._sequence_path = Path(fasta_path).expanduser().resolve()
        # Inside a tab a page does not paint its own background, and the white
        # Preview page would show through.
        self.setAutoFillBackground(True)

        root = QVBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)

        self.sections = QTabBar(self)
        self.sections.setObjectName("PrimersWorkspaceSections")
        self.sections.setExpanding(False)
        self.sections.setDrawBase(False)
        for key, label in SECTIONS:
            index = self.sections.addTab(label)
            self.sections.setTabData(index, key)
        root.addWidget(self.sections)

        row = QHBoxLayout()
        row.setContentsMargins(12, 6, 12, 6)
        row.setSpacing(8)
        row.addWidget(QLabel("Sequence", self))
        self.sequence_row = PathBrowseRow(
            self,
            initial=self._sequence_path,
            browse_title="Choose the Sequence to Work On",
            file_filter=FASTA_FILTER,
            start_dir=str(self._sequence_path.parent),
            elide_middle=True,
            minimum_width=260,
        )
        self.sequence_row.setMaximumWidth(460)
        self.sequence_row.set_default_tooltip(
            "The file the workspace works on. Click a FASTA in Files to use it, "
            "or press Browse. Double-click a file in Files to open it in its own tab."
        )
        row.addWidget(self.sequence_row, 1)
        self.record_label = QLabel("Record", self)
        row.addWidget(self.record_label)
        self.record_combo = QComboBox(self)
        self.record_combo.setMinimumWidth(180)
        self.record_combo.setToolTip("For a file with several sequences: the one to work on.")
        row.addWidget(self.record_combo)
        row.addStretch(2)
        from repeatviewer.widgets.help_button import make_help_button

        row.addWidget(make_help_button(self, "Sequence", SEQUENCE_HELP))
        root.addLayout(row)

        self.panel = SequenceMapPanel(self._sequence_path, self, autoload=False)
        root.addWidget(self.panel, 1)
        self.panel.statusMessage.connect(self._on_panel_status)
        self.panel.annotationsSaved.connect(self.annotationsSaved.emit)
        self.panel.sequenceExtracted.connect(self.sequenceExtracted.emit)
        self.panel.activeSequenceChanged.connect(self._on_panel_record_changed)

        rail = self.panel._rail
        # The workspace shows one page of the rail and never another, so the
        # strip of rail buttons has no job here.
        rail._tab_strip.hide()
        rail.set_panel_collapsed(False, persist=False)

        self.primers = self.panel.primers_panel()
        self._armed_row = None
        self._build_list_area()
        self._build_settings_header()
        app = QApplication.instance()
        if app is not None:
            app.focusChanged.connect(self._on_focus_changed)
        self.sections.currentChanged.connect(self._on_section_changed)
        self.primers.tabs.currentChanged.connect(self._on_panel_page_changed)
        self.record_combo.activated.connect(self._on_record_chosen)
        self.sequence_row.path_changed.connect(self._on_sequence_field_changed)

        self._load(self._sequence_path)
        self.set_section(section, page=PANEL_DESIGN)

    def _on_panel_status(self, text: str) -> None:
        # Changing the sequence ends a Transfer session the view never had here;
        # its notice would only confuse.
        if str(text).startswith("Transfer stopped"):
            return
        self.statusMessage.emit(str(text))

    # ------------------------------------------------------------------ naming

    def tab_name(self) -> str:
        return TAB_NAME

    def sequence_path(self) -> Path:
        return self._sequence_path

    def primers_panel(self):
        return self.primers

    # ---------------------------------------------------------------- sections

    def section(self) -> str:
        return str(self.sections.tabData(self.sections.currentIndex()) or SECTION_DESIGN)

    def set_section(self, section: str, *, page: int | None = None) -> None:
        """Show ``section``.

        Design has two pages (settings, primers): it keeps the one it was on,
        unless ``page`` names one.
        """
        key = section_key(section)
        if page is not None and int(page) in _PANEL_PAGES[key]:
            self._syncing = True
            try:
                self.primers.tabs.setCurrentIndex(int(page))
            finally:
                self._syncing = False
        for index in range(self.sections.count()):
            if self.sections.tabData(index) == key:
                if self.sections.currentIndex() != index:
                    self.sections.setCurrentIndex(index)   # runs _on_section_changed
                else:
                    self._apply_section(key)
                return

    def _on_section_changed(self, _index: int) -> None:
        if self._syncing:
            return
        self._apply_section(self.section())

    def _apply_section(self, key: str) -> None:
        pages = _PANEL_PAGES[key]
        tabs = self.primers.tabs
        before = tabs.currentIndex()
        target = before if before in pages else pages[0]
        if target == PANEL_PRIMER and not self.primers._selected:
            target = PANEL_DESIGN           # the primer page needs a primer
        self._syncing = True
        try:
            tabs.setCurrentIndex(target)
            # The sections above and the list under the map do the panel's
            # own tab bar's job.
            tabs.setVisible(False)
        finally:
            self._syncing = False
        if target == before:
            # The page stayed the same, so its signal did not fire: let it draw
            # its overlays (region and candidates, products, enzyme sites).
            self.primers._on_workbench_tab(target)
        self.bottom.setVisible(key in (SECTION_DESIGN, SECTION_PCR))
        self.bottom_stack.setCurrentIndex(1 if key == SECTION_PCR else 0)
        self._sync_settings_summary()

    def _on_panel_page_changed(self, index: int) -> None:
        """The panel moved by itself (kept pairs open the primer page, say)."""
        if self._syncing:
            return
        key = section_of_panel_page(index)
        if key != self.section():
            self.set_section(key)
        elif int(index) == PANEL_PRIMER and not self.primers._selected:
            # nothing to show there: stay on the settings
            QTimer.singleShot(0, lambda: self._show_design_page(PANEL_DESIGN))
        self._sync_settings_summary()

    # ------------------------------------------------- the list under the map

    def _build_list_area(self) -> None:
        """Candidates and primers on the sequence under the map; PCR products too."""
        from repeatviewer.widgets.help_button import make_help_button
        from repeatviewer.widgets.primers_workspace_list import WorkspaceList

        primers = self.primers
        design = primers.design_pane
        # The two tables stay the source of the list, but are never shown.
        self._kept = QWidget(self)
        self._kept.hide()
        candidates_section = design.table.parentWidget()
        self.list = WorkspaceList(
            design.table, primers.table, self,
            primer_columns=(primers.COL_TICK, primers.COL_NAME, primers.COL_TM,
                            primers.COL_SITE))
        self.list.context_menu_for_primer = primers._on_context_menu
        self.list.setMinimumHeight(96)

        self.bottom = QWidget(self)
        self.bottom.setAutoFillBackground(True)
        outer = QVBoxLayout(self.bottom)
        outer.setContentsMargins(12, 6, 12, 8)
        outer.setSpacing(0)
        self.bottom_stack = QStackedWidget(self.bottom)
        outer.addWidget(self.bottom_stack)

        page = QWidget(self.bottom_stack)
        lay = QVBoxLayout(page)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(6)
        lay.addWidget(self.list, 1)
        bar = QHBoxLayout()
        bar.setContentsMargins(0, 0, 0, 0)
        bar.setSpacing(8)
        primers.add_button.setText("Add Primer")
        primers.database_button.setText("From Database")
        self.list_note = design.cand_note
        for button in (primers.add_button, primers.database_button):
            button.setSizePolicy(QPushButton().sizePolicy())
            bar.addWidget(button)
        bar.addWidget(self.list_note, 1)
        bar.addWidget(make_help_button(page, "Primers", LIST_HELP))
        bar.addWidget(design.report_button)
        bar.addWidget(design.keep_button)
        lay.addLayout(bar)
        self.bottom_stack.addWidget(page)

        products = QWidget(self.bottom_stack)
        pl = QVBoxLayout(products)
        pl.setContentsMargins(0, 0, 0, 0)
        primers.pcr_pane.table.setMinimumHeight(96)
        pl.addWidget(primers.pcr_pane.table)
        self.bottom_stack.addWidget(products)

        design.table.setParent(self._kept)
        primers.table.setParent(self._kept)
        candidates_section.hide()
        # the two buttons of the old primer page that the workspace replaces
        primers.design_button.hide()
        primers.pcr_button.hide()

        # under the map, left of the settings; the user can drag the border
        layout = self.panel.layout()
        index = layout.indexOf(self.panel._map)
        self.splitter = QSplitter(Qt.Orientation.Vertical, self.panel)
        self.splitter.setChildrenCollapsible(False)
        layout.insertWidget(max(0, index), self.splitter, 1)
        self.splitter.addWidget(self.panel._map)
        self.splitter.addWidget(self.bottom)
        self.splitter.setStretchFactor(0, 1)
        self.splitter.setStretchFactor(1, 0)
        self.splitter.setSizes([600, 190])

        self.list.picked.connect(self._on_list_picked)
        self.list.cellClicked.connect(self._on_list_cell_clicked)

    def _on_list_picked(self, kind: str) -> None:
        from repeatviewer.widgets.primers_workspace_list import KIND_PRIMER

        if self.section() != SECTION_DESIGN:
            return
        self._show_design_page(PANEL_PRIMER if kind == KIND_PRIMER else PANEL_DESIGN)

    def _on_list_cell_clicked(self, row: int, column: int) -> None:
        """A click on the row that is already selected opens its page again."""
        from repeatviewer.widgets.primers_workspace_list import COL_TICK

        kind = self.list.kind_of_row(row)
        if kind and column != COL_TICK:
            self._on_list_picked(kind)

    def _show_design_page(self, page: int) -> None:
        """Design settings, or the details of the selected primer."""
        if self.section() != SECTION_DESIGN:
            return
        if int(page) == PANEL_PRIMER and not self.primers._selected:
            page = PANEL_DESIGN
        self._syncing = True
        try:
            self.primers.tabs.setCurrentIndex(int(page))
        finally:
            self._syncing = False
        self._sync_settings_summary()

    # ------------------------------------------- settings folded to one line

    def _build_settings_header(self) -> None:
        """On top of the primer details: the design settings as one line."""
        from repeatviewer.widgets.help_button import make_help_button

        design = self.primers.design_pane
        page = self.primers.pages.widget(PANEL_PRIMER)
        header = QWidget(page)
        lay = QVBoxLayout(header)
        lay.setContentsMargins(0, 0, 0, 6)
        lay.setSpacing(6)
        top = QHBoxLayout()
        top.setContentsMargins(0, 0, 0, 0)
        top.addWidget(QLabel("<b>Design Settings</b>", header), 1)
        self.show_settings_button = QPushButton("Show", header)
        top.addWidget(self.show_settings_button)
        top.addWidget(make_help_button(header, "Design Settings", SETTINGS_HEADER_HELP))
        lay.addLayout(top)
        self.settings_summary = QLabel(header)
        self.settings_summary.setWordWrap(True)
        self.settings_summary.setStyleSheet("color: palette(placeholder-text);")
        lay.addWidget(self.settings_summary)
        row = QHBoxLayout()
        row.setContentsMargins(0, 0, 0, 0)
        row.addStretch(1)
        self.header_defaults_button = QPushButton("Restore Defaults", header)
        self.header_design_button = QPushButton("Design", header)
        row.addWidget(self.header_defaults_button)
        row.addWidget(self.header_design_button)
        lay.addLayout(row)
        page.layout().insertWidget(0, header)
        self.settings_header = header

        self.show_settings_button.clicked.connect(
            lambda: self._show_design_page(PANEL_DESIGN))
        self.header_defaults_button.clicked.connect(design.defaults_button.click)
        self.header_defaults_button.clicked.connect(self._sync_settings_summary)

        def design_again() -> None:
            self._show_design_page(PANEL_DESIGN)
            design.design_button.click()

        self.header_design_button.clicked.connect(design_again)

    def settings_summary_text(self) -> str:
        design = self.primers.design_pane
        parts = [design.mode.currentText()]
        try:
            region = design.region()
        except Exception:  # noqa: BLE001 - a half-typed range is no region
            region = None
        if region is not None:
            parts.append(f"region {region.start:,} to {region.end:,}")
        if not design.is_probe_mode():
            parts.append(f"Tm {design.tm_min.value():g} to {design.tm_max.value():g} °C")
        return ", ".join(parts)

    def _sync_settings_summary(self, *_args) -> None:
        self.settings_summary.setText(self.settings_summary_text())

    # ------------------------------------------ path fields filled from Files

    def _on_focus_changed(self, _old, new) -> None:
        """Remember the path field the cursor was in when the user goes to Files."""
        if new is None:
            return
        try:
            if not self.isAncestorOf(new):
                return                      # focus left the workspace: keep the field
        except RuntimeError:
            return
        row = _path_row_of(new)
        self._armed_row = row if row is not None and row is not self.sequence_row else None

    def fill_armed_path(self, path: Path | str) -> bool:
        """Put a file clicked in Files into the path field the cursor was in."""
        row = self._armed_row
        if row is None:
            return False
        try:
            if not row.isVisible():
                return False
            row.set_path(Path(path))
        except RuntimeError:
            self._armed_row = None
            return False
        self._armed_row = None
        return True

    # ---------------------------------------------------------------- sequence

    def _load(self, path: Path) -> None:
        self._sequence_path = Path(path).expanduser().resolve()
        self.panel.load_path(self._sequence_path)
        rail = self.panel._rail
        rail._tab_strip.hide()
        rail.set_panel_collapsed(False, persist=False)
        self.panel.show_primers_tab()
        self._sync_records()

    def _sync_records(self) -> None:
        source = self.panel._sequence_combo
        self.record_combo.blockSignals(True)
        try:
            self.record_combo.clear()
            for index in range(source.count()):
                self.record_combo.addItem(str(source.itemData(index) or source.itemText(index)),
                                          source.itemData(index))
            self.record_combo.setCurrentIndex(max(0, source.currentIndex()))
        finally:
            self.record_combo.blockSignals(False)
        several = self.record_combo.count() > 1
        self.record_label.setVisible(several)
        self.record_combo.setVisible(several)

    def _on_record_chosen(self, index: int) -> None:
        seq_id = self.record_combo.itemData(index)
        if seq_id:
            self.panel.ensure_active_sequence(str(seq_id))

    def _on_panel_record_changed(self, seq_id: str) -> None:
        index = self.record_combo.findData(seq_id)
        if index >= 0 and index != self.record_combo.currentIndex():
            self.record_combo.blockSignals(True)
            self.record_combo.setCurrentIndex(index)
            self.record_combo.blockSignals(False)

    def _show_sequence_in_field(self) -> None:
        self.sequence_row.set_path(self._sequence_path, emit=False)
        self.sequence_row.set_start_dir(str(self._sequence_path.parent))

    def _on_sequence_field_changed(self) -> None:
        """Browse, or a path typed into the field."""
        chosen = self.sequence_row.path()
        if chosen is None or same_file(chosen, self._sequence_path):
            self._show_sequence_in_field()      # an emptied field shows the file again
            return
        if not self.set_sequence(chosen):
            self._show_sequence_in_field()

    def set_sequence(self, path: Path | str) -> bool:
        """Work on another FASTA. False when it is not one, or the user stayed."""
        target = Path(path).expanduser()
        if not target.is_file():
            self.statusMessage.emit(f"{target.name} is not a file.")
            return False
        if same_file(target, self._sequence_path):
            self._show_sequence_in_field()
            return True
        if not self.confirm_leave():
            return False
        self._load(target)
        self._show_sequence_in_field()
        # The page in front draws its overlays for the new sequence.
        self._apply_section(self.section())
        self.sequenceChanged.emit(self._sequence_path)
        return True

    def reload_annotations(self) -> None:
        """The sequence's annotations were saved somewhere else: show them."""
        if self.has_unsaved_changes():
            return
        record = self.record_combo.currentData()
        self.panel.load_path(self._sequence_path, fit_after_load=False)
        if record:
            self.panel.ensure_active_sequence(str(record))
        self.panel.show_primers_tab()
        self._sync_records()
        self._apply_section(self.section())

    # ------------------------------------------------------------ leaving safely

    def has_unsaved_changes(self) -> bool:
        try:
            return bool(self.panel.has_unsaved_annotation_changes())
        except RuntimeError:
            return False

    def confirm_leave(self) -> bool:
        """Unsaved annotation changes: Save / Don't Save / Cancel. True to go on."""
        if not self.has_unsaved_changes():
            return True
        box = QMessageBox(self.window())
        box.setIcon(QMessageBox.Icon.Warning)
        box.setWindowTitle("Unsaved Changes")
        box.setText(f'Save changes to "{self._sequence_path.name}"?')
        box.setStandardButtons(
            QMessageBox.StandardButton.Save
            | QMessageBox.StandardButton.Discard
            | QMessageBox.StandardButton.Cancel
        )
        box.setDefaultButton(QMessageBox.StandardButton.Save)
        box.setEscapeButton(QMessageBox.StandardButton.Cancel)
        choice = box.exec()
        if choice == QMessageBox.StandardButton.Save:
            return bool(self.panel.save_annotations(silent=True))
        return choice == QMessageBox.StandardButton.Discard
