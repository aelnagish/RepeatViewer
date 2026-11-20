"""Sequence list for a large FASTA (Preview page shown instead of "too big to show").

Mockup: Design canvas "Large genomes mockup", screen 1 "Pick a sequence from a
large genome". The file is indexed once on a worker thread (backends/large_fasta),
then every sequence is listed with its length; the user opens one sequence, a
region of it, or writes it to its own FASTA. Nothing is loaded into memory.

Step 1 (16 Nov 2026): Open writes the chosen region next to the genome
(`<genome>_<name>_<start>-<end>.fa`, reused while newer than the genome) and opens
that file. Step 2 replaces this with a window view in genome positions.
"""

from __future__ import annotations

import threading
from pathlib import Path

from PySide6.QtCore import QObject, QRectF, Qt, Signal
from PySide6.QtGui import QColor, QFont, QPainter, QPalette
from PySide6.QtWidgets import (
    QAbstractItemView,
    QFileDialog,
    QFrame,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QMessageBox,
    QProgressBar,
    QPushButton,
    QSizePolicy,
    QStyle,
    QStyledItemDelegate,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends import large_fasta as lf

_COL_NAME, _COL_LENGTH, _COL_BAR, _COL_ANNOT = range(4)
_LENGTH_ROLE = Qt.ItemDataRole.UserRole + 1
_FRACTION_ROLE = Qt.ItemDataRole.UserRole + 2

_HELP = (
    "<b>Sequences</b>: this file is too large to load at once, so RepeatViewer "
    "indexed it (a <b>.fai</b> next to the file, as samtools does) and lists its "
    "sequences. Only the part you open is read from disk.<br>"
    "<b>Double-click</b> a sequence: up to 100 Mb opens whole; a longer one opens "
    "its first 5 Mb.<br>"
    "<b>Region</b>: type <b>scaffold_5:120,000,000-125,000,000</b>, a name alone "
    "(whole sequence) or <b>1000-2000</b> (selected sequence). At most 100 Mb.<br>"
    "<b>Open Region</b> / <b>Open Whole Sequence</b>: write the part next to the "
    "genome as <b>&lt;genome&gt;_&lt;name&gt;_&lt;start&gt;-&lt;end&gt;.fa</b> and "
    "open it. The header keeps the genome positions.<br>"
    "<b>Extract to File</b>: the same, to a file name you choose.<br>"
    "<b>.gz files</b>: bgzip files open as fast as plain ones. A plain gzip file "
    "is read from its start for every region; <b>Convert to bgzip</b> makes "
    "<b>&lt;name&gt;_bgzip.fa.gz</b> once."
)


class _Relay(QObject):
    """Carries worker-thread results to the GUI thread (no module-level QObject)."""

    progress = Signal(int, int)
    indexed = Signal(object)
    counted = Signal(object)
    written = Signal(object, bool)  # path, open it after
    converted = Signal(object)
    failed = Signal(str)


class _BarDelegate(QStyledItemDelegate):
    """Length bar relative to the longest sequence."""

    def paint(self, painter: QPainter, option, index) -> None:
        super().paint(painter, option, index)
        fraction = index.data(_FRACTION_ROLE)
        if fraction is None:
            return
        rect = option.rect.adjusted(6, 0, -6, 0)
        h = 10.0
        top = rect.center().y() - h / 2 + 1
        width = max(4.0, rect.width() * float(fraction))
        colour = option.palette.color(QPalette.ColorRole.Highlight)
        if not option.state & QStyle.StateFlag.State_Selected:
            colour = QColor(colour)
            colour.setAlpha(110)
        painter.save()
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        painter.setPen(Qt.PenStyle.NoPen)
        painter.setBrush(colour)
        painter.drawRoundedRect(QRectF(rect.left(), top, width, h), h / 2, h / 2)
        painter.restore()


class _LengthItem(QTableWidgetItem):
    def __lt__(self, other: QTableWidgetItem) -> bool:  # numeric sort
        a = self.data(_LENGTH_ROLE)
        b = other.data(_LENGTH_ROLE)
        if a is None or b is None:
            return super().__lt__(other)
        return int(a) < int(b)


class SequenceListPage(QWidget):
    """Preview page listing the sequences of a large FASTA."""

    openRequested = Signal(object)  # Path of a written region to open in Preview
    fileWritten = Signal(object)  # Path of a new file (Files panel turns it green)
    statusMessage = Signal(str)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._path: Path | None = None
        self._index: lf.FastaIndex | None = None
        self._cancel = threading.Event()
        self._busy = False
        self._relay = _Relay(self)
        self._relay.progress.connect(self._on_progress)
        self._relay.indexed.connect(self._on_indexed)
        self._relay.counted.connect(self._on_counted)
        self._relay.written.connect(self._on_written)
        self._relay.converted.connect(self._on_converted)
        self._relay.failed.connect(self._on_failed)
        self._build()

    # ------------------------------------------------------------------ UI

    def _build(self) -> None:
        from repeatviewer.widgets.help_button import make_help_button

        root = QVBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)

        body = QWidget(self)
        lay = QVBoxLayout(body)
        lay.setContentsMargins(20, 16, 20, 10)
        lay.setSpacing(10)

        top = QHBoxLayout()
        top.setSpacing(12)
        self._title = QLabel("Sequences", body)
        font = QFont(self._title.font())
        font.setPointSizeF(font.pointSizeF() * 1.3)
        font.setBold(True)
        self._title.setFont(font)
        self._summary = self._muted(QLabel("", body))
        self._filter = QLineEdit(body)
        self._filter.setPlaceholderText("Filter names")
        self._filter.setClearButtonEnabled(True)
        self._filter.setFixedWidth(200)
        self._filter.textChanged.connect(self._apply_filter)
        top.addWidget(self._title, 0, Qt.AlignmentFlag.AlignBaseline)
        top.addWidget(self._summary, 0, Qt.AlignmentFlag.AlignBaseline)
        top.addStretch(1)
        top.addWidget(self._filter)
        top.addWidget(make_help_button(body, "Sequences", _HELP))
        lay.addLayout(top)

        # Indexing progress / gzip advice, shown only when needed.
        self._notice = QFrame(body)
        nl = QHBoxLayout(self._notice)
        nl.setContentsMargins(0, 0, 0, 0)
        nl.setSpacing(10)
        self._notice_text = QLabel("", self._notice)
        self._notice_text.setWordWrap(True)
        self._notice_text.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred)
        self._bar = QProgressBar(self._notice)
        self._bar.setFixedWidth(220)
        self._bar.setTextVisible(False)
        self._notice_button = QPushButton("Stop", self._notice)
        self._notice_button.clicked.connect(self._on_notice_button)
        nl.addWidget(self._notice_text, 1)
        nl.addWidget(self._bar)
        nl.addWidget(self._notice_button)
        self._notice.hide()
        lay.addWidget(self._notice)

        self._table = QTableWidget(0, 4, body)
        self._table.setHorizontalHeaderLabels(["Name", "Length", "Size", "Annotations"])
        self._table.verticalHeader().hide()
        self._table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self._table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self._table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self._table.setAlternatingRowColors(False)
        self._table.setShowGrid(False)
        self._table.setItemDelegateForColumn(_COL_BAR, _BarDelegate(self._table))
        header = self._table.horizontalHeader()
        header.setDefaultAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        header.setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
        header.setStretchLastSection(False)
        header.setSortIndicator(-1, Qt.SortOrder.AscendingOrder)  # file order until a click
        self._table.setSortingEnabled(True)
        for col, width in ((_COL_NAME, 180), (_COL_LENGTH, 110), (_COL_BAR, 420), (_COL_ANNOT, 110)):
            self._table.setColumnWidth(col, width)
        self._table.itemSelectionChanged.connect(self._on_selection)
        self._table.itemDoubleClicked.connect(self._on_double_click)
        lay.addWidget(self._table, 1)

        self._hint = self._muted(QLabel(
            "Double-click a sequence to open it (longer than 100 Mb: its first 5 Mb).", body))
        lay.addWidget(self._hint)
        root.addWidget(body, 1)

        line = QFrame(self)
        line.setFrameShape(QFrame.Shape.HLine)
        line.setFrameShadow(QFrame.Shadow.Sunken)
        root.addWidget(line)

        bottom = QWidget(self)
        bottom.setAutoFillBackground(True)
        bl = QHBoxLayout(bottom)
        bl.setContentsMargins(20, 10, 20, 10)
        bl.setSpacing(10)
        region_label = QLabel("Region", bottom)
        bold = QFont(region_label.font())
        bold.setBold(True)
        region_label.setFont(bold)
        self._region = QLineEdit(bottom)
        self._region.setFixedWidth(300)
        self._region.setPlaceholderText("scaffold_5:120,000,000-125,000,000")
        self._region.textChanged.connect(self._sync_region_note)
        self._region.returnPressed.connect(self._open_region)
        region_label.setBuddy(self._region)
        self._region_note = self._muted(QLabel("", bottom))
        self._extract_button = QPushButton("Extract to File…", bottom)
        self._whole_button = QPushButton("Open Whole Sequence", bottom)
        self._open_button = QPushButton("Open Region", bottom)
        self._open_button.setDefault(True)
        self._extract_button.clicked.connect(self._extract)
        self._whole_button.clicked.connect(self._open_whole)
        self._open_button.clicked.connect(self._open_region)
        bl.addWidget(region_label)
        bl.addWidget(self._region)
        bl.addWidget(self._region_note)
        bl.addStretch(1)
        bl.addWidget(self._extract_button)
        bl.addWidget(self._whole_button)
        bl.addWidget(self._open_button)
        root.addWidget(bottom)
        self._set_actions_enabled(False)

    def _muted(self, label: QLabel) -> QLabel:
        pal = label.palette()
        pal.setColor(QPalette.ColorRole.WindowText,
                     self.palette().color(QPalette.ColorRole.PlaceholderText))
        label.setPalette(pal)
        return label

    # ---------------------------------------------------------------- loading

    def path(self) -> Path | None:
        return self._path

    def load(self, path: Path) -> None:
        path = Path(path).expanduser().resolve()
        if self._path == path and self._index is not None:
            return
        self._cancel.set()
        self._cancel = threading.Event()
        self._path = path
        self._index = None
        self._table.setRowCount(0)
        self._region.clear()
        self._title.setText("Sequences")
        self._summary.setText(f"{path.name} · {self._disk_size(path)} on disk")
        index = lf.load_index(path)
        if index is not None:
            self._on_indexed(index)
            return
        self._show_notice("Indexing the file once (it is read from start to end)…",
                          button="Stop", progress=True)
        cancel = self._cancel
        relay = self._relay

        def work() -> None:
            try:
                result = lf.build_index(path, progress=lambda d, t: relay.progress.emit(d, t),
                                        cancel=cancel.is_set)
            except lf.IndexCancelled:
                return
            except Exception as exc:  # shown to the user
                if not cancel.is_set():
                    relay.failed.emit(f"Could not index {path.name}: {exc}")
                return
            if not cancel.is_set():
                relay.indexed.emit(result)

        threading.Thread(target=work, name="fasta-index", daemon=True).start()

    def stop(self) -> None:
        """Stop background work (page closed or another file shown)."""
        self._cancel.set()

    @staticmethod
    def _disk_size(path: Path) -> str:
        try:
            size = path.stat().st_size
        except OSError:
            return "?"
        for unit in ("B", "KB", "MB", "GB", "TB"):
            if size < 1024 or unit == "TB":
                return f"{size:.0f} {unit}" if unit == "B" else f"{size:.1f} {unit}"
            size /= 1024
        return ""

    def _on_progress(self, done: int, total: int) -> None:
        if total > 0:
            self._bar.setRange(0, 1000)
            self._bar.setValue(int(1000 * done / total))
            self._notice_text.setText(
                f"Indexing the file once: {self._gb(done)} of {self._gb(total)} read…")

    @staticmethod
    def _gb(n: int) -> str:
        return f"{n / 1024 ** 3:.1f} GB" if n >= 1024 ** 3 else f"{n / 1024 ** 2:.0f} MB"

    def _on_indexed(self, index: lf.FastaIndex) -> None:
        if self._path is None or Path(index.path) != self._path:
            return
        self._index = index
        self._hide_notice()
        n = len(index.entries)
        self._title.setText(f"{n:,} sequence" + ("" if n == 1 else "s"))
        kind = {lf.KIND_BGZF: " (bgzip)", lf.KIND_GZIP: " (gzip)"}.get(index.kind, "")
        self._summary.setText(
            f"{lf.format_bp(index.total_length)} · {self._disk_size(index.path)} on disk{kind}"
            " · index ready")
        self._fill_table(index)
        if index.kind == lf.KIND_GZIP:
            self._show_notice(
                "This file is gzip, not bgzip: every region is read from the start of the "
                "file, which is slow for a genome. Convert it once to bgzip to open regions "
                "at once.", button="Convert to bgzip", progress=False)
        self._count_annotations(index)

    def _fill_table(self, index: lf.FastaIndex) -> None:
        longest = max(1, index.longest)
        self._table.setSortingEnabled(False)
        self._table.setRowCount(len(index.entries))
        for row, entry in enumerate(index.entries):
            name = QTableWidgetItem(entry.name)
            if entry.description:
                name.setToolTip(f"{entry.name} {entry.description}")
            length = _LengthItem(lf.format_bp(entry.length))
            length.setData(_LENGTH_ROLE, entry.length)
            length.setToolTip(f"{entry.length:,} bp")
            bar = _LengthItem("")
            bar.setData(_LENGTH_ROLE, entry.length)
            bar.setData(_FRACTION_ROLE, entry.length / longest)
            annot = QTableWidgetItem("")
            for col, item in ((_COL_NAME, name), (_COL_LENGTH, length), (_COL_BAR, bar),
                              (_COL_ANNOT, annot)):
                item.setTextAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
                self._table.setItem(row, col, item)
        self._table.setSortingEnabled(True)
        self._apply_filter(self._filter.text())
        if index.entries:
            self._table.selectRow(0)
        self._set_actions_enabled(True)

    def _count_annotations(self, index: lf.FastaIndex) -> None:
        relay = self._relay
        names = index.names()
        path = Path(index.path)

        def work() -> None:
            try:
                from repeatviewer.backends.sequence_features import load_features_by_id_from_json

                feats = load_features_by_id_from_json(path, names)
            except Exception:
                feats = {}
            relay.counted.emit((path, {k: len(v) for k, v in feats.items()}))

        threading.Thread(target=work, name="fasta-annotations", daemon=True).start()

    def _on_counted(self, payload) -> None:
        path, counts = payload
        if self._index is None or Path(self._index.path) != path:
            return
        for row in range(self._table.rowCount()):
            name = self._table.item(row, _COL_NAME).text()
            n = counts.get(name, 0)
            self._table.item(row, _COL_ANNOT).setText(f"{n:,}" if n else "None")

    # -------------------------------------------------------------- notices

    def _show_notice(self, text: str, *, button: str, progress: bool) -> None:
        self._notice_text.setText(text)
        self._notice_button.setText(button)
        self._bar.setVisible(progress)
        self._bar.setRange(0, 0)
        self._notice.show()

    def _hide_notice(self) -> None:
        self._notice.hide()

    def _on_notice_button(self) -> None:
        label = self._notice_button.text()
        if label == "Stop":
            self._cancel.set()
            self._busy = False
            self._hide_notice()
            if self._index is None:
                self._summary.setText("Indexing stopped. Open the file again to index it.")
            else:
                self._set_actions_enabled(True)
            return
        if label == "Convert to bgzip" and self._index is not None:
            self._convert()

    def _on_failed(self, message: str) -> None:
        self._busy = False
        self._hide_notice()
        self._set_actions_enabled(self._index is not None)
        QMessageBox.warning(self, "Sequences", message)

    # -------------------------------------------------------------- actions

    def _set_actions_enabled(self, on: bool) -> None:
        for w in (self._region, self._extract_button, self._whole_button, self._open_button):
            w.setEnabled(on and not self._busy)

    def _apply_filter(self, text: str) -> None:
        needle = (text or "").strip().lower()
        for row in range(self._table.rowCount()):
            item = self._table.item(row, _COL_NAME)
            self._table.setRowHidden(row, bool(needle) and needle not in item.text().lower())

    def _selected_name(self) -> str | None:
        rows = self._table.selectionModel().selectedRows() if self._table.selectionModel() else []
        if not rows:
            return None
        item = self._table.item(rows[0].row(), _COL_NAME)
        return item.text() if item else None

    def _on_selection(self) -> None:
        name = self._selected_name()
        if name is None or self._index is None:
            return
        region = lf.opening_region(self._index, name)
        self._region.setText(f"{name}:{region.start:,}-{region.end:,}")

    def _sync_region_note(self, _text: str = "") -> None:
        region = self._parse(quiet=True)
        self._region_note.setText(lf.format_bp(region.length) if region else "")

    def _parse(self, *, quiet: bool = False) -> lf.Region | None:
        if self._index is None:
            return None
        try:
            return lf.parse_region(self._region.text(), self._index, self._selected_name())
        except lf.LargeFastaError as exc:
            if not quiet:
                QMessageBox.information(self, "Region", str(exc))
            return None

    def _on_double_click(self, item: QTableWidgetItem) -> None:
        if self._index is None:
            return
        name = self._table.item(item.row(), _COL_NAME).text()
        self._write(lf.opening_region(self._index, name), open_after=True)

    def _open_region(self) -> None:
        region = self._parse()
        if region is not None:
            self._write(region, open_after=True)

    def _open_whole(self) -> None:
        name = self._selected_name()
        if name is None or self._index is None:
            return
        entry = self._index.entry(name)
        self._write(lf.Region(name, 1, max(1, entry.length)), open_after=True)

    def _extract(self) -> None:
        region = self._parse()
        if region is None or self._index is None:
            return
        default = lf.numbered_path(lf.region_output_path(self._index, region))
        target, _ = QFileDialog.getSaveFileName(
            self, "Extract to File", str(default), "FASTA (*.fa *.fasta *.fna)")
        if target:
            self._write(region, open_after=False, target=Path(target))

    def _write(self, region: lf.Region, *, open_after: bool, target: Path | None = None) -> None:
        index = self._index
        if index is None or self._busy:
            return
        if open_after and region.length > lf.WHOLE_SEQUENCE_LIMIT:
            QMessageBox.information(
                self, "Region",
                f"{region.label()} is {lf.format_bp(region.length)}. The Annotation view shows "
                "at most 100 Mb: choose a shorter region, or use Extract to File.")
            return
        out = target or lf.region_output_path(index, region)
        if target is None and self._fresh(out, Path(index.path)):
            self._on_written(out, open_after)
            return
        self._busy = True
        self._set_actions_enabled(False)
        self._show_notice(f"Writing {region.label()}…", button="Stop", progress=True)
        cancel = self._cancel
        relay = self._relay

        def work() -> None:
            try:
                lf.write_region(index, region, out,
                                progress=lambda d, t: relay.progress.emit(d, t),
                                cancel=cancel.is_set)
            except lf.IndexCancelled:
                return
            except Exception as exc:
                relay.failed.emit(f"Could not write {out.name}: {exc}")
                return
            relay.written.emit(out, open_after)

        threading.Thread(target=work, name="fasta-region", daemon=True).start()

    @staticmethod
    def _fresh(out: Path, source: Path) -> bool:
        try:
            return out.is_file() and out.stat().st_mtime >= source.stat().st_mtime
        except OSError:
            return False

    def _on_written(self, out, open_after: bool) -> None:
        self._busy = False
        self._hide_notice()
        if self._index is not None and self._index.kind == lf.KIND_GZIP:
            self._show_notice(
                "This file is gzip, not bgzip: every region is read from the start of the "
                "file, which is slow for a genome. Convert it once to bgzip to open regions "
                "at once.", button="Convert to bgzip", progress=False)
        self._set_actions_enabled(True)
        out = Path(out)
        self.fileWritten.emit(out)
        self.statusMessage.emit(f"Wrote {out.name}")
        if open_after:
            self.openRequested.emit(out)

    def _convert(self) -> None:
        index = self._index
        if index is None or self._busy:
            return
        self._busy = True
        self._set_actions_enabled(False)
        self._show_notice("Converting to bgzip (the file is read once)…", button="Stop",
                          progress=True)
        cancel = self._cancel
        relay = self._relay
        source = Path(index.path)

        def work() -> None:
            try:
                out = lf.convert_to_bgzip(source, progress=lambda d, t: relay.progress.emit(d, t),
                                          cancel=cancel.is_set)
            except lf.IndexCancelled:
                return
            except Exception as exc:
                relay.failed.emit(f"Could not convert {source.name}: {exc}")
                return
            relay.converted.emit(out)

        threading.Thread(target=work, name="fasta-bgzip", daemon=True).start()

    def _on_converted(self, out) -> None:
        self._busy = False
        self._hide_notice()
        self._set_actions_enabled(True)
        out = Path(out)
        self.fileWritten.emit(out)
        self.statusMessage.emit(f"Wrote {out.name}: open it to work with regions at once.")
        self.openRequested.emit(out)
