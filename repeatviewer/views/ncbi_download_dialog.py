"""Download from NCBI window (Tools > Download from NCBI).

Tabs: **Search** | **Options**. Type accessions (one or many) or a search
term, press Search, scroll through the results (more load as you scroll),
tick what you want and press Run to download into the folder shown in
"Save to". A click on a folder or file in the panels fills "Save to".

Searching runs on a Python thread; results come back through a Qt signal, so
the window never freezes while NCBI answers.
"""

from __future__ import annotations

import threading
from dataclasses import replace
from pathlib import Path

from PySide6.QtCore import QObject, Qt, QTimer, QUrl, Signal
from PySide6.QtGui import QDesktopServices, QGuiApplication
from PySide6.QtWidgets import (
    QAbstractItemView,
    QButtonGroup,
    QCheckBox,
    QDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QMenu,
    QMessageBox,
    QPushButton,
    QRadioButton,
    QSizePolicy,
    QStackedWidget,
    QTabBar,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends import ncbi_download as nd
from repeatviewer.job_options import (
    DIALOG_TAB_OPTIONS,
    NCBI_DIALOG_TAB_SEARCH,
    NCBI_DIALOG_TABS,
    NCBI_RIGHT_CLICK_CLIPBOARD,
    NCBI_RIGHT_CLICK_DIALOG,
    NcbiDownloadOptions,
    opening_dialog_tab,
)
from repeatviewer.views.tool_dialogs import (
    TOOL_PARAMS_DIALOG_WIDTH_WIDE,
    dialog_should_run,
    make_cancel_save_run_bar,
    make_process_output_box,
)
from repeatviewer.widgets.combo_boxes import NoWheelComboBox as QComboBox
from repeatviewer.widgets.help_button import field_with_help, make_help_button
from repeatviewer.widgets.path_browse import PathBrowseRow
from repeatviewer.widgets.spin_boxes import NoWheelSpinBox as QSpinBox

_DIALOG_HEIGHT = 640
_DIALOG_WIDTH = max(TOOL_PARAMS_DIALOG_WIDTH_WIDE, 860)
# Load the next page when the scroll bar is this close to the bottom (rows).
_SCROLL_MARGIN_ROWS = 8

SEARCH_HELP = (
    "<b>Accessions</b>: type or paste one or more, separated by spaces, commas or new "
    "lines, e.g. <b>NC_045512.2 MG815809</b>. They are found directly and ticked.<br>"
    "<b>Genome assemblies</b>: <b>GCA_</b> / <b>GCF_</b> accessions switch the database "
    "to Genome Assembly.<br>"
    "<b>Search term</b>: any words, e.g. <b>CroSat satellite</b>, or NCBI search syntax "
    "such as <b>ITS2[Title] AND Crocus[Organism]</b>.<br>"
    "Press <b>Enter</b> or <b>Search</b>. More results load when you scroll to the end of "
    "the list."
)
DATABASE_HELP = (
    "<b>Nucleotide</b>: GenBank, RefSeq and WGS DNA / RNA records.<br>"
    "<b>Protein</b>: protein records (saved as <b>.faa</b>).<br>"
    "<b>Genome Assembly</b>: whole genomes from NCBI Datasets (genome FASTA, optionally "
    "with its GFF3 annotation). Genomes can be several GB; the download runs in the "
    "background with a Stop button."
)
FILTERS_HELP = (
    "Filters are added to the search term (ignored when you search by accession).<br>"
    "<b>Organism</b>: scientific name or taxon, e.g. <b>Crocus sativus</b> or "
    "<b>Iridaceae</b>.<br>"
    "<b>Length</b>: sequence length in bp (aa for protein); 0 means no limit.<br>"
    "<b>Molecule</b>: genomic DNA, mRNA or rRNA records only.<br>"
    "<b>RefSeq only</b>: curated NCBI Reference Sequences.<br>"
    "<b>Latest only</b>: skip replaced versions of an assembly."
)
FORMAT_HELP = (
    "<b>FASTA</b>: sequences only.<br>"
    "<b>FASTA with NCBI annotations</b>: FASTA plus the record's NCBI features (genes, "
    "CDS, repeats, ...) added to the Annotation track.<br>"
    "<b>GenBank / GenPept</b>: the full NCBI flat file, as text.<br>"
    "<b>Genome FASTA + GFF3</b>: the GFF3 is saved next to the genome with the same "
    "name, so tools that use a companion annotation pick it up."
)
SAVE_AS_HELP = (
    "<b>One file</b>: all ticked records in one FASTA / GenBank file, named in "
    "<b>File name</b>.<br>"
    "<b>One file per record</b>: each record in its own file named after its accession, "
    "e.g. <b>NC_045512.2.fa</b>.<br>"
    "Genome assemblies are always one file per assembly."
)
SAVE_TO_HELP = (
    "Folder the files are written to. It starts as the folder open in the Files panel; "
    "click a folder or file in the panels to change it."
)
ACCOUNT_HELP = (
    "NCBI asks programs to send an <b>email</b> so they can contact you if something "
    "goes wrong; it is optional.<br>"
    "An <b>API key</b> raises the limit from 3 to 10 requests per second, which makes "
    "large downloads faster. Create one for free under <b>Account settings</b> on "
    "ncbi.nlm.nih.gov."
)
RIGHT_CLICK_HELP = (
    "<b>Download the accessions on the clipboard</b>: copy accessions from a paper or a "
    "table, then right-click the toolbar button. They are downloaded with the saved "
    "format into the folder open in the Files panel.<br>"
    "<b>Open this window</b>: right-click does the same as a click."
)

_COLUMNS = {
    nd.DB_NUCLEOTIDE: ("", "Accession", "Description", "Organism", "Length (bp)", "Type", "Updated"),
    nd.DB_PROTEIN: ("", "Accession", "Description", "Organism", "Length (aa)", "Type", "Updated"),
    nd.DB_ASSEMBLY: ("", "Accession", "Assembly", "Organism", "Size (bp)", "Level", "Released"),
}


def format_length(n: int) -> str:
    return f"{int(n):,}" if n else ""


def results_summary(count: int, shown: int, *, accessions: bool = False) -> str:
    if count == 0:
        return "No records found."
    noun = "record" if count == 1 else "records"
    if accessions or shown >= count:
        return f"{count:,} {noun}."
    return f"{count:,} {noun}, showing {shown:,}. Scroll down to load more."


def _wide_with_help(field: QWidget, title: str, help_text: str) -> QWidget:
    """Like ``field_with_help`` but the field takes the row's free width."""
    wrap = QWidget()
    row = QHBoxLayout(wrap)
    row.setContentsMargins(0, 0, 0, 0)
    row.setSpacing(6)
    row.addWidget(field, 1)
    row.addWidget(make_help_button(wrap, title, help_text), 0, Qt.AlignmentFlag.AlignVCenter)
    return wrap


class _SearchBridge(QObject):
    """Carries search results from the Python thread back to the GUI thread."""

    page_ready = Signal(int, object)
    page_failed = Signal(int, str)


class NcbiDownloadDialog(QDialog):
    """Search NCBI, tick records, Run downloads them."""

    _TAB_SEARCH = 0
    _TAB_OPTIONS = 1

    def __init__(
        self,
        *,
        initial: NcbiDownloadOptions,
        folder: Path | None,
        query_text: str = "",
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Download from NCBI")
        self.setObjectName("NcbiDownloadDialog")
        self._opts = initial
        self._database = initial.database
        self._formats = {
            nd.DB_NUCLEOTIDE: initial.nucleotide_format,
            nd.DB_PROTEIN: initial.protein_format,
            nd.DB_ASSEMBLY: initial.assembly_format,
        }
        self._records: list[nd.NcbiRecord] = []
        self._ticked: set[str] = set()
        self._total = 0
        self._query: nd.NcbiQuery | None = None
        self._accession_mode = False
        self._loading = False
        self._request_id = 0
        self._bridge = _SearchBridge(self)
        self._bridge.page_ready.connect(self._on_page_ready)
        self._bridge.page_failed.connect(self._on_page_failed)

        # ---------------------------------------------------------- Search tab
        self._db_combo = QComboBox(self)
        for key in nd.DATABASES:
            self._db_combo.addItem(nd.DATABASE_LABELS[key], key)
        self._set_combo(self._db_combo, self._database)
        self._db_combo.currentIndexChanged.connect(self._on_database_changed)

        self._search_edit = QLineEdit(self)
        self._search_edit.setClearButtonEnabled(True)
        self._search_edit.setPlaceholderText("Accession(s) or search term, e.g. NC_045512.2 or CroSat")
        self._search_edit.returnPressed.connect(self.start_search)
        self._search_edit.textChanged.connect(self._on_text_changed)
        self._search_btn = QPushButton("Search", self)
        self._search_btn.setAutoDefault(False)
        self._search_btn.clicked.connect(self.start_search)
        search_row = QWidget(self)
        srl = QHBoxLayout(search_row)
        srl.setContentsMargins(0, 0, 0, 0)
        srl.setSpacing(6)
        srl.addWidget(self._search_edit, 1)
        srl.addWidget(self._search_btn)

        self._organism = QLineEdit(self)
        self._organism.setClearButtonEnabled(True)
        self._organism.setPlaceholderText("e.g. Crocus sativus")
        self._organism.returnPressed.connect(self.start_search)
        self._min_len = QSpinBox(self)
        self._max_len = QSpinBox(self)
        for box in (self._min_len, self._max_len):
            box.setRange(0, 2_000_000_000)
            box.setSingleStep(100)
            box.setSpecialValueText("Any")
            box.setFixedWidth(110)
        self._molecule = QComboBox(self)
        for key, label in nd.MOLECULE_LABELS:
            self._molecule.addItem(label, key)
        self._molecule.setFixedWidth(130)
        self._refseq_only = QCheckBox("RefSeq only", self)
        self._latest_only = QCheckBox("Latest only", self)
        self._latest_only.setChecked(True)

        filters = QWidget(self)
        fl = QHBoxLayout(filters)
        fl.setContentsMargins(0, 0, 0, 0)
        fl.setSpacing(6)
        fl.addWidget(self._organism, 2)
        self._len_label = QLabel("Length:", self)
        fl.addWidget(self._len_label)
        fl.addWidget(self._min_len)
        self._len_dash = QLabel("to", self)
        fl.addWidget(self._len_dash)
        fl.addWidget(self._max_len)
        self._mol_label = QLabel("Molecule:", self)
        fl.addWidget(self._mol_label)
        fl.addWidget(self._molecule)
        fl.addWidget(self._refseq_only)
        fl.addWidget(self._latest_only)

        top = QFormLayout()
        top.setLabelAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        top.setHorizontalSpacing(10)
        top.setVerticalSpacing(8)
        db_row = QWidget(self)
        dbl = QHBoxLayout(db_row)
        dbl.setContentsMargins(0, 0, 0, 0)
        dbl.addWidget(field_with_help(self._db_combo, "Database", DATABASE_HELP))
        dbl.addStretch(1)
        top.addRow("Database:", db_row)
        top.addRow("Search:", _wide_with_help(search_row, "Search", SEARCH_HELP))
        top.addRow("Organism:", _wide_with_help(filters, "Filters", FILTERS_HELP))

        self._summary = QLabel("Type an accession or a search term, then press Search.", self)
        self._summary.setWordWrap(True)
        self._summary.setStyleSheet("color: palette(placeholder-text);")
        self._tick_all = QCheckBox("Tick all shown", self)
        self._tick_all.setTristate(False)
        self._tick_all.clicked.connect(self._on_tick_all)
        self._ticked_label = QLabel("", self)
        summary_row = QHBoxLayout()
        summary_row.setContentsMargins(0, 0, 0, 0)
        summary_row.addWidget(self._tick_all)
        summary_row.addSpacing(12)
        summary_row.addWidget(self._summary, 1)
        summary_row.addWidget(self._ticked_label)

        self._table = QTableWidget(0, 7, self)
        self._table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self._table.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self._table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self._table.setAlternatingRowColors(True)
        self._table.verticalHeader().setVisible(False)
        self._table.setWordWrap(False)
        self._table.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self._table.customContextMenuRequested.connect(self._show_table_menu)
        self._table.itemChanged.connect(self._on_item_changed)
        self._table.cellDoubleClicked.connect(self._toggle_row)
        self._table.verticalScrollBar().valueChanged.connect(self._maybe_load_more)
        self._table.setMinimumHeight(220)
        self._apply_columns()

        self._format = QComboBox(self)
        self._format.currentIndexChanged.connect(self._on_format_changed)
        self._one_file = QRadioButton("One file", self)
        self._per_record = QRadioButton("One file per record", self)
        save_group = QButtonGroup(self)
        save_group.addButton(self._one_file)
        save_group.addButton(self._per_record)
        (self._one_file if initial.one_file else self._per_record).setChecked(True)
        self._one_file.toggled.connect(self._sync_file_name)
        save_as = QWidget(self)
        sal = QHBoxLayout(save_as)
        sal.setContentsMargins(0, 0, 0, 0)
        sal.addWidget(self._one_file)
        sal.addWidget(self._per_record)
        sal.addStretch(1)
        self._file_name = QLineEdit(self)
        self._file_name.setPlaceholderText("Named after the search")
        self._file_name_edited = False
        self._file_name.textEdited.connect(lambda _t: setattr(self, "_file_name_edited", True))
        self._save_to = PathBrowseRow(
            self,
            initial=folder,
            placeholder="Click a folder in the Folders panel",
            browse_title="Save NCBI downloads to",
            mode="directory",
        )

        out = QFormLayout()
        out.setLabelAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        out.setHorizontalSpacing(10)
        out.setVerticalSpacing(8)
        fmt_row = QWidget(self)
        frl = QHBoxLayout(fmt_row)
        frl.setContentsMargins(0, 0, 0, 0)
        frl.addWidget(field_with_help(self._format, "Format", FORMAT_HELP))
        frl.addStretch(1)
        out.addRow("Format:", fmt_row)
        self._save_as_label = QLabel("Save as:", self)
        self._save_as_field = field_with_help(save_as, "Save as", SAVE_AS_HELP, stretch=True)
        out.addRow(self._save_as_label, self._save_as_field)
        self._file_name_label = QLabel("File name:", self)
        out.addRow(self._file_name_label, self._file_name)
        out.addRow("Save to:", _wide_with_help(self._save_to, "Save to", SAVE_TO_HELP))

        search_page = QWidget(self)
        spl = QVBoxLayout(search_page)
        spl.setContentsMargins(4, 12, 4, 6)
        spl.setSpacing(8)
        spl.addLayout(top)
        spl.addLayout(summary_row)
        spl.addWidget(self._table, 1)
        spl.addLayout(out)

        # --------------------------------------------------------- Options tab
        self._right_click_buttons: dict[str, QRadioButton] = {}
        rc_group = QButtonGroup(self)
        rc_box = QGroupBox("Right-click on Download from NCBI", self)
        rcl = QVBoxLayout(rc_box)
        for label, value in (
            ("Download the accessions on the clipboard with saved settings", NCBI_RIGHT_CLICK_CLIPBOARD),
            ("Open this window", NCBI_RIGHT_CLICK_DIALOG),
        ):
            radio = QRadioButton(label, self)
            radio.setChecked(initial.right_click_action == value)
            rc_group.addButton(radio)
            self._right_click_buttons[value] = radio
        rcl.addWidget(field_with_help(
            self._right_click_buttons[NCBI_RIGHT_CLICK_CLIPBOARD], "Right-click", RIGHT_CLICK_HELP
        ))
        rcl.addWidget(self._right_click_buttons[NCBI_RIGHT_CLICK_DIALOG])

        self._email = QLineEdit(initial.email, self)
        self._email.setPlaceholderText("optional")
        self._api_key = QLineEdit(initial.api_key, self)
        self._api_key.setPlaceholderText("optional")
        self._api_key.setEchoMode(QLineEdit.EchoMode.PasswordEchoOnEdit)
        self._page_size = QSpinBox(self)
        self._page_size.setRange(20, nd.MAX_PAGE_SIZE)
        self._page_size.setSingleStep(20)
        self._page_size.setValue(int(initial.page_size))
        self._page_size.setFixedWidth(110)
        account = QGroupBox("NCBI account", self)
        al = QFormLayout(account)
        al.setLabelAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        al.addRow("Email:", _wide_with_help(self._email, "NCBI account", ACCOUNT_HELP))
        al.addRow("API key:", self._api_key)
        al.addRow("Results per page:", self._page_size)

        process_box, self._show_process_window, self._close_process_window = (
            make_process_output_box(
                self,
                tool_name="the download",
                checked=initial.show_process_window,
                close_when_finished=initial.close_process_window_when_finished,
            )
        )
        options_page = QWidget(self)
        ol = QVBoxLayout(options_page)
        ol.setContentsMargins(16, 14, 16, 14)
        ol.setSpacing(14)
        ol.setAlignment(Qt.AlignmentFlag.AlignTop)
        ol.addWidget(rc_box)
        ol.addWidget(account)
        ol.addWidget(process_box)

        # -------------------------------------------------------------- frame
        self._tab_bar = QTabBar(self)
        self._tab_bar.setExpanding(False)
        self._tab_bar.setDocumentMode(True)
        self._tab_bar.addTab("Search")
        self._tab_bar.addTab("Options")
        self._stack = QStackedWidget(self)
        self._stack.addWidget(search_page)
        self._stack.addWidget(options_page)
        self._tab_bar.currentChanged.connect(self._stack.setCurrentIndex)
        open_tab = opening_dialog_tab(initial.last_dialog_tab, allowed=NCBI_DIALOG_TABS,
                                      first_tab=NCBI_DIALOG_TAB_SEARCH)
        self._tab_bar.setCurrentIndex(
            self._TAB_OPTIONS if open_tab == DIALOG_TAB_OPTIONS else self._TAB_SEARCH
        )
        self._stack.setCurrentIndex(self._tab_bar.currentIndex())

        restore = QPushButton("Restore Defaults", self)
        restore.setAutoDefault(False)
        restore.clicked.connect(self._restore_defaults)
        self._button_box = make_cancel_save_run_bar(
            self, left_widgets=[restore], on_accept=self._on_accept
        )
        # Enter searches; it must never start a download by accident.
        for btn in (getattr(self, "_save_button", None), getattr(self, "_run_button", None)):
            if btn is not None:
                btn.setDefault(False)
                btn.setAutoDefault(False)
        if self._run_button is not None:
            self._run_button.setToolTip("Download the ticked records")

        layout = QVBoxLayout(self)
        layout.setContentsMargins(22, 18, 22, 16)
        layout.setSpacing(8)
        layout.addWidget(self._tab_bar, 0)
        layout.addWidget(self._stack, 1)
        layout.addWidget(self._button_box, 0)

        self._refresh_formats()
        self._sync_filters()
        self._sync_ticked()
        self.resize(_DIALOG_WIDTH, _DIALOG_HEIGHT)
        if query_text:
            self._search_edit.setText(query_text)
            QTimer.singleShot(0, self.start_search)
        else:
            QTimer.singleShot(0, self._search_edit.setFocus)

    # ------------------------------------------------------------- helpers

    @staticmethod
    def _set_combo(combo: QComboBox, value) -> None:
        idx = combo.findData(value)
        if idx >= 0:
            combo.setCurrentIndex(idx)

    def database(self) -> str:
        return str(self._db_combo.currentData() or nd.DB_NUCLEOTIDE)

    def current_dialog_tab(self) -> str:
        return DIALOG_TAB_OPTIONS if self._tab_bar.currentIndex() == self._TAB_OPTIONS \
            else NCBI_DIALOG_TAB_SEARCH

    def _apply_columns(self) -> None:
        cols = _COLUMNS[self.database()]
        self._table.setColumnCount(len(cols))
        self._table.setHorizontalHeaderLabels(list(cols))
        hh = self._table.horizontalHeader()
        hh.setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
        hh.setSectionResizeMode(0, QHeaderView.ResizeMode.Fixed)
        hh.setSectionResizeMode(2, QHeaderView.ResizeMode.Stretch)
        for col, width in enumerate((28, 120, 0, 140, 95, 105, 85)):
            if width:
                self._table.setColumnWidth(col, width)

    def _sync_filters(self) -> None:
        db = self.database()
        seq_db = db in (nd.DB_NUCLEOTIDE, nd.DB_PROTEIN)
        for w in (self._len_label, self._min_len, self._len_dash, self._max_len):
            w.setVisible(seq_db)
        for w in (self._mol_label, self._molecule):
            w.setVisible(db == nd.DB_NUCLEOTIDE)
        self._latest_only.setVisible(db == nd.DB_ASSEMBLY)
        is_assembly = db == nd.DB_ASSEMBLY
        self._save_as_label.setVisible(not is_assembly)
        self._save_as_field.setVisible(not is_assembly)
        self._file_name_label.setVisible(not is_assembly)
        self._file_name.setVisible(not is_assembly)
        self._sync_file_name()

    def _refresh_formats(self) -> None:
        db = self.database()
        self._format.blockSignals(True)
        self._format.clear()
        for key, label in nd.FORMATS_BY_DB[db]:
            self._format.addItem(label, key)
        self._set_combo(self._format, self._formats.get(db))
        self._format.blockSignals(False)

    def _on_format_changed(self, _index: int) -> None:
        self._formats[self.database()] = str(self._format.currentData())
        self._sync_file_name()

    def _on_database_changed(self, _index: int) -> None:
        self._database = self.database()
        self._refresh_formats()
        self._sync_filters()
        self._clear_results("Press Search to look in the new database.")
        self._apply_columns()

    def _on_text_changed(self, text: str) -> None:
        guess = nd.guess_database(text)
        if guess and guess != self.database():
            self._set_combo(self._db_combo, guess)

    # ------------------------------------------------------------- search

    def current_query(self) -> nd.NcbiQuery:
        return nd.NcbiQuery(
            database=self.database(),
            text=self._search_edit.text().strip(),
            organism=self._organism.text().strip(),
            min_length=self._min_len.value(),
            max_length=self._max_len.value(),
            molecule=str(self._molecule.currentData() or ""),
            refseq_only=self._refseq_only.isChecked(),
            latest_only=self._latest_only.isChecked(),
        )

    def _client(self) -> nd.NcbiClient:
        return nd.NcbiClient(email=self._email.text().strip(), api_key=self._api_key.text().strip())

    def _clear_results(self, message: str) -> None:
        self._request_id += 1
        self._records = []
        self._ticked = set()
        self._total = 0
        self._loading = False
        self._query = None
        self._table.setRowCount(0)
        self._summary.setText(message)
        self._sync_ticked()

    def start_search(self) -> None:
        query = self.current_query()
        try:
            nd.build_term(query)
        except nd.NcbiError as exc:
            self._summary.setText(str(exc))
            return
        self._clear_results("Searching NCBI…")
        self._query = query
        self._accession_mode = bool(nd.looks_like_accessions(query.text))
        if not self._file_name_edited:
            self._file_name.clear()
        self._fetch_page(0)

    def _page_size_value(self) -> int:
        return int(self._page_size.value())

    def _fetch_page(self, retstart: int) -> None:
        if self._query is None:
            return
        self._loading = True
        self._search_btn.setEnabled(False)
        request_id = self._request_id
        query = self._query
        client = self._client()
        size = max(self._page_size_value(), len(nd.looks_like_accessions(query.text)))
        bridge = self._bridge

        def _work() -> None:
            try:
                page = nd.search(client, query, retstart=retstart, retmax=size)
            except Exception as exc:  # noqa: BLE001 - shown in the window
                message = str(exc) if isinstance(exc, nd.NcbiError) else f"{type(exc).__name__}: {exc}"
                try:
                    bridge.page_failed.emit(request_id, message)
                except RuntimeError:
                    pass  # window closed meanwhile
                return
            try:
                bridge.page_ready.emit(request_id, page)
            except RuntimeError:
                pass

        threading.Thread(target=_work, name="ncbi-search", daemon=True).start()

    def _on_page_ready(self, request_id: int, page) -> None:
        if request_id != self._request_id:
            return
        self._loading = False
        self._search_btn.setEnabled(True)
        self._total = int(page.count)
        self.add_records(page.records, tick=self._accession_mode)
        self._summary.setText(results_summary(self._total, len(self._records),
                                              accessions=self._accession_mode))
        if self._accession_mode:
            wanted = nd.looks_like_accessions(self._query.text if self._query else "")
            missing = missing_accessions(wanted, self._records)
            if missing:
                self._summary.setText(
                    self._summary.text() + " Not found: " + ", ".join(missing[:8])
                    + ("…" if len(missing) > 8 else "")
                )
        QTimer.singleShot(0, self._fill_short_view)

    def _on_page_failed(self, request_id: int, message: str) -> None:
        if request_id != self._request_id:
            return
        self._loading = False
        self._search_btn.setEnabled(True)
        self._summary.setText(message)

    def has_more(self) -> bool:
        return not self._accession_mode and len(self._records) < self._total

    def _fill_short_view(self) -> None:
        """Keep loading while the list is too short to scroll (nothing to scroll to)."""
        if self._table.verticalScrollBar().maximum() == 0:
            self._maybe_load_more()

    def _maybe_load_more(self, *_args) -> None:
        if self._loading or not self.has_more() or self._query is None:
            return
        bar = self._table.verticalScrollBar()
        row_h = max(1, self._table.verticalHeader().defaultSectionSize())
        if bar.maximum() == 0 or bar.value() >= bar.maximum() - _SCROLL_MARGIN_ROWS * row_h:
            self._summary.setText(f"Loading more… ({len(self._records):,} of {self._total:,})")
            self._fetch_page(len(self._records))

    def add_records(self, records: list[nd.NcbiRecord], *, tick: bool = False) -> None:
        """Append result rows (also used by tests without network)."""
        known = {r.uid for r in self._records}
        self._table.blockSignals(True)
        for rec in records:
            if rec.uid in known:
                continue
            known.add(rec.uid)
            self._records.append(rec)
            row = self._table.rowCount()
            self._table.insertRow(row)
            box = QTableWidgetItem()
            box.setFlags(Qt.ItemFlag.ItemIsUserCheckable | Qt.ItemFlag.ItemIsEnabled
                         | Qt.ItemFlag.ItemIsSelectable)
            if tick:
                self._ticked.add(rec.uid)
            box.setCheckState(Qt.CheckState.Checked if rec.uid in self._ticked
                              else Qt.CheckState.Unchecked)
            self._table.setItem(row, 0, box)
            values = (rec.accession, rec.title, rec.organism, format_length(rec.length),
                      rec.kind, rec.updated)
            for col, value in enumerate(values, start=1):
                item = QTableWidgetItem(value)
                if col == 4:
                    item.setTextAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
                if col == 2 and rec.title:
                    item.setToolTip(rec.title)
                self._table.setItem(row, col, item)
        self._table.blockSignals(False)
        self._sync_ticked()

    # -------------------------------------------------------------- ticks

    def _on_item_changed(self, item: QTableWidgetItem) -> None:
        if item.column() != 0:
            return
        row = item.row()
        if not 0 <= row < len(self._records):
            return
        uid = self._records[row].uid
        if item.checkState() == Qt.CheckState.Checked:
            self._ticked.add(uid)
        else:
            self._ticked.discard(uid)
        self._sync_ticked()

    def _toggle_row(self, row: int, _col: int = 0) -> None:
        item = self._table.item(row, 0)
        if item is None:
            return
        item.setCheckState(Qt.CheckState.Unchecked if item.checkState() == Qt.CheckState.Checked
                           else Qt.CheckState.Checked)

    def _on_tick_all(self, checked: bool) -> None:
        self._table.blockSignals(True)
        state = Qt.CheckState.Checked if checked else Qt.CheckState.Unchecked
        for row, rec in enumerate(self._records):
            item = self._table.item(row, 0)
            if item is not None:
                item.setCheckState(state)
            if checked:
                self._ticked.add(rec.uid)
            else:
                self._ticked.discard(rec.uid)
        self._table.blockSignals(False)
        self._sync_ticked()

    def _sync_ticked(self) -> None:
        n = len(self._ticked)
        self._ticked_label.setText(f"{n:,} ticked" if n else "")
        self._tick_all.blockSignals(True)
        self._tick_all.setChecked(bool(self._records) and n == len(self._records))
        self._tick_all.blockSignals(False)
        self._sync_file_name()

    def chosen_records(self) -> list[nd.NcbiRecord]:
        """Ticked records; with nothing ticked, the selected rows."""
        ticked = [r for r in self._records if r.uid in self._ticked]
        if ticked:
            return ticked
        rows = sorted({i.row() for i in self._table.selectedIndexes()})
        return [self._records[r] for r in rows if 0 <= r < len(self._records)]

    def _sync_file_name(self, *_args) -> None:
        records = self.chosen_records() if hasattr(self, "_records") else []
        combined = (self.database() != nd.DB_ASSEMBLY and self._one_file.isChecked()
                    and len(records) > 1)
        self._file_name.setEnabled(combined)
        self._file_name_label.setEnabled(combined)
        if combined and not self._file_name_edited:
            self._file_name.setText(nd.default_file_name(
                self.database(), self.output_format(), records, self._search_edit.text()))
        elif not combined and not self._file_name_edited:
            self._file_name.clear()

    # -------------------------------------------------------- table menu

    def _show_table_menu(self, pos) -> None:
        records = self.chosen_records()
        row = self._table.rowAt(pos.y())
        if not records and 0 <= row < len(self._records):
            records = [self._records[row]]
        if not records:
            return
        menu = QMenu(self)
        db = self.database()
        menu.addAction("Open on NCBI Website",
                       lambda: [QDesktopServices.openUrl(QUrl(r.web_url(db))) for r in records[:10]])
        menu.addAction("Copy Accession" if len(records) == 1 else "Copy Accessions",
                       lambda: QGuiApplication.clipboard().setText(
                           "\n".join(r.accession for r in records)))
        menu.exec(self._table.viewport().mapToGlobal(pos))

    # ------------------------------------------------------- panel click

    def apply_path_from_panel(self, path: Path) -> bool:
        path = Path(path)
        folder = path if path.is_dir() else path.parent
        if not folder.is_dir():
            return False
        self._save_to.set_path(folder)
        return True

    # ------------------------------------------------------------ results

    def output_format(self) -> str:
        return str(self._format.currentData() or nd.FORMAT_FASTA)

    def one_file(self) -> bool:
        return self.database() != nd.DB_ASSEMBLY and self._one_file.isChecked()

    def file_name(self) -> str:
        text = self._file_name.text().strip()
        if not text:
            return ""
        suffix = nd.output_suffix(self.database(), self.output_format())
        return text if Path(text).suffix else text + suffix

    def out_dir(self) -> Path | None:
        return self._save_to.path()

    def query_text(self) -> str:
        return self._search_edit.text().strip()

    def should_run(self) -> bool:
        return dialog_should_run(self)

    def options(self) -> NcbiDownloadOptions:
        rc = next((k for k, b in self._right_click_buttons.items() if b.isChecked()),
                  NCBI_RIGHT_CLICK_DIALOG)
        return replace(
            self._opts,
            database=self.database(),
            nucleotide_format=self._formats[nd.DB_NUCLEOTIDE],
            protein_format=self._formats[nd.DB_PROTEIN],
            assembly_format=self._formats[nd.DB_ASSEMBLY],
            one_file=self._one_file.isChecked(),
            page_size=int(self._page_size.value()),
            email=self._email.text().strip(),
            api_key=self._api_key.text().strip(),
            right_click_action=rc,
            show_process_window=self._show_process_window.isChecked(),
            close_process_window_when_finished=self._close_process_window.isChecked(),
            last_dialog_tab=self.current_dialog_tab(),
        )

    def _on_accept(self) -> None:
        if getattr(self, "_close_action", "") == "run":
            if not self.chosen_records():
                QMessageBox.information(self, "Download from NCBI",
                                        "Tick the records to download first.")
                self._close_action = "cancel"
                return
            folder = self.out_dir()
            if folder is None or not folder.is_dir():
                QMessageBox.information(self, "Download from NCBI",
                                        "Choose the folder to save to: click a folder in the "
                                        "Folders panel or use Browse.")
                self._close_action = "cancel"
                return
        self._request_id += 1  # drop late search answers
        self.accept()

    def reject(self) -> None:  # noqa: D401 - Qt override
        self._request_id += 1
        super().reject()

    def _restore_defaults(self) -> None:
        d = NcbiDownloadOptions()
        self._formats = {
            nd.DB_NUCLEOTIDE: d.nucleotide_format,
            nd.DB_PROTEIN: d.protein_format,
            nd.DB_ASSEMBLY: d.assembly_format,
        }
        self._refresh_formats()
        (self._one_file if d.one_file else self._per_record).setChecked(True)
        self._organism.clear()
        self._min_len.setValue(0)
        self._max_len.setValue(0)
        self._molecule.setCurrentIndex(0)
        self._refseq_only.setChecked(False)
        self._latest_only.setChecked(True)
        self._page_size.setValue(d.page_size)
        self._right_click_buttons[d.right_click_action].setChecked(True)
        self._show_process_window.setChecked(d.show_process_window)
        self._close_process_window.setChecked(d.close_process_window_when_finished)
        self._file_name_edited = False
        self._sync_file_name()


def missing_accessions(wanted: list[str], records: list[nd.NcbiRecord]) -> list[str]:
    """Typed accessions NCBI did not return (version-insensitive)."""
    found = set()
    for r in records:
        acc = r.accession.upper()
        found.add(acc)
        found.add(acc.split(".")[0])
    return [w for w in wanted if w.upper() not in found and w.upper().split(".")[0] not in found]
