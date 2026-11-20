"""Sanger Check: .ab1 reads of a clone against the sequence it should carry.

Choose the expected sequence (the open file, or any FASTA), add the reads, and
the window says whether the insert is confirmed. Every difference is listed
with the quality behind it, and the chromatogram is drawn around the base so
the call can be judged by eye.

Developer: Abdullah Yusuf
"""

from __future__ import annotations

import html
from pathlib import Path

from PySide6.QtCore import QRectF, Qt
from PySide6.QtGui import QColor, QFont, QPainter, QPainterPath, QPen
from PySide6.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QDialog,
    QFileDialog,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends.primers import sanger as sg
from repeatviewer.widgets.help_button import make_help_button

BASE_COLORS = {"G": "#1a1a1a", "A": "#1a7f4b", "T": "#d0342c", "C": "#1f57c3"}
OK_COLOR = "#1a7f4b"
CHECK_COLOR = "#b26b00"
WINDOW_BASES = 21          # bases shown around the chosen one

SANGER_HELP = (
    "<b>Expected</b>: the sequence the clone should carry (the open file, or Browse). "
    "Pick the record when the file holds several.<br>"
    "<b>Reads</b>: the .ab1 files from sequencing (plain .seq or FASTA work too, but "
    "without a chromatogram or quality).<br>"
    "Each read is trimmed by quality (Mott's rule, as phred uses), then aligned to the "
    "expected sequence in both directions. A difference is <b>doubtful</b> below Q30: the "
    "machine, not the clone. The verdict is green only when the reads cover every base and "
    "nothing is left over.<br>"
    "Click a difference to see the chromatogram around it. <b>Save Report</b> writes an "
    "HTML page beside the reads."
)
TRACE_HELP = (
    "The four dyes around the chosen base: G black, A green, T red, C blue. The called "
    "base is on top, the expected one underneath when they differ. A clean single peak "
    "with a high quality score is a real difference in the clone; a double peak, a shoulder "
    "or a low score is a poor call."
)


class TraceView(QWidget):
    """The chromatogram around one base of one read."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._check: sg.ReadCheck | None = None
        self._centre = -1
        self._expected = ""
        self._position = 0
        self.setMinimumHeight(170)
        self.setToolTip("The dye signals around the chosen base.")

    def show_at(self, check: sg.ReadCheck | None, index: int, expected_base: str = "",
                position: int = 0) -> None:
        self._check, self._centre = check, int(index)
        self._expected, self._position = expected_base, int(position)
        self.update()

    def _message(self, painter: QPainter, text: str) -> None:
        painter.setPen(QColor("#6b7280"))
        painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, text)

    def paintEvent(self, event) -> None:  # noqa: N802
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        painter.fillRect(self.rect(), QColor("#ffffff"))
        trace = self._check.trace if self._check is not None else None
        if trace is None or self._centre < 0:
            self._message(painter, "Click a difference to see the chromatogram.")
            painter.end()
            return
        if not trace.has_trace:
            self._message(painter, "This file holds no chromatogram, only bases.")
            painter.end()
            return
        half = WINDOW_BASES // 2
        first = max(0, self._centre - half)
        last = min(len(trace.locations) - 1, self._centre + half)
        if last <= first:
            self._message(painter, "No trace around this base.")
            painter.end()
            return
        x0, x1 = trace.locations[first], trace.locations[last]
        span = max(1, x1 - x0)
        w, h = self.width(), self.height()
        top, bottom = 26.0, h - 26.0
        peak = max((max(v[x0:x1 + 1], default=1) for v in trace.channels.values()), default=1) or 1
        for base, values in trace.channels.items():
            path = QPainterPath()
            for i, x in enumerate(range(x0, min(x1 + 1, len(values)))):
                px = i / span * (w - 16) + 8
                py = bottom - (values[x] / peak) * (bottom - top)
                path.moveTo(px, py) if i == 0 else path.lineTo(px, py)
            painter.setPen(QPen(QColor(BASE_COLORS.get(base, "#888888")), 1.4))
            painter.drawPath(path)
        font = QFont(painter.font())
        font.setPointSizeF(max(8.0, font.pointSizeF()))
        painter.setFont(font)
        for i in range(first, last + 1):
            px = (trace.locations[i] - x0) / span * (w - 16) + 8
            base = trace.bases[i] if i < len(trace.bases) else "?"
            if i == self._centre:
                painter.fillRect(QRectF(px - 7, 2, 14, h - 4), QColor(255, 196, 0, 46))
            painter.setPen(QColor(BASE_COLORS.get(base, "#555555")))
            painter.drawText(QRectF(px - 7, 2, 14, 18), Qt.AlignmentFlag.AlignCenter, base)
            if i == self._centre and self._expected and self._expected != base:
                painter.setPen(QColor(CHECK_COLOR))
                painter.drawText(QRectF(px - 9, h - 22, 18, 18),
                                 Qt.AlignmentFlag.AlignCenter, self._expected.lower())
        painter.setPen(QColor("#6b7280"))
        note = f"read {self._check.name}, base {self._centre + 1:,}"
        if self._position:
            note += f" · expected {self._position:,}"
        painter.drawText(QRectF(8, h - 20, w - 16, 18),
                         Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter, note)
        painter.end()


class SangerDialog(QDialog):
    """``expected`` is the sequence the clone should carry; reads are added by the user."""

    def __init__(self, expected: str = "", *, expected_name: str = "",
                 expected_path: Path | None = None, folder: Path | None = None,
                 parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Sanger Check")
        self.setObjectName("SangerDialog")
        self._expected = str(expected or "").upper()
        self._expected_name = expected_name or "expected sequence"
        self._folder = Path(folder) if folder else (
            expected_path.parent if expected_path else Path.cwd())
        self._checks: list[sg.ReadCheck] = []
        self._rows: list[tuple[sg.ReadCheck, sg.Difference | None]] = []
        self._paths: list[Path] = []

        lay = QVBoxLayout(self)
        lay.setSpacing(8)
        head = QHBoxLayout()
        self.expected_label = QLabel(self)
        self.expected_label.setWordWrap(True)
        head.addWidget(self.expected_label, 1)
        self.expected_button = QPushButton("Expected Sequence", self)
        self.add_button = QPushButton("Add Reads", self)
        head.addWidget(self.expected_button)
        head.addWidget(self.add_button)
        head.addWidget(make_help_button(self, "Sanger Check", SANGER_HELP))
        lay.addLayout(head)

        self.record_combo = QComboBox(self)
        self.record_combo.hide()
        lay.addWidget(self.record_combo)

        self.verdict = QLabel(self)
        self.verdict.setWordWrap(True)
        self.verdict.setTextFormat(Qt.TextFormat.RichText)
        self.verdict.setObjectName("sangerVerdict")
        lay.addWidget(self.verdict)

        self.table = QTableWidget(0, 4, self)
        self.table.setHorizontalHeaderLabels(["Read", "Result", "Difference", "Quality"])
        self.table.verticalHeader().hide()
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setShowGrid(False)
        self.table.setAlternatingRowColors(True)
        hh = self.table.horizontalHeader()
        hh.setDefaultAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        hh.setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
        for c, w in ((0, 150), (2, 190), (3, 80)):
            self.table.setColumnWidth(c, w)
        self.table.setMinimumHeight(190)
        lay.addWidget(self.table, 1)

        trace_head = QHBoxLayout()
        trace_head.addWidget(QLabel("<b>Chromatogram</b>", self), 1)
        trace_head.addWidget(make_help_button(self, "Chromatogram", TRACE_HELP))
        lay.addLayout(trace_head)
        self.trace = TraceView(self)
        lay.addWidget(self.trace)

        buttons = QHBoxLayout()
        self.report_button = QPushButton("Save Report", self)
        close = QPushButton("Close", self)
        buttons.addStretch(1)
        buttons.addWidget(self.report_button)
        buttons.addWidget(close)
        lay.addLayout(buttons)

        self.add_button.clicked.connect(self._on_add_reads)
        self.expected_button.clicked.connect(self._on_pick_expected)
        self.record_combo.currentIndexChanged.connect(self._on_record)
        self.table.itemSelectionChanged.connect(self._on_row)
        self.report_button.clicked.connect(self._on_report)
        close.clicked.connect(self.accept)
        self.resize(820, 640)
        self._refresh()

    # ------------------------------------------------------------------ inputs

    def set_expected(self, sequence: str, name: str = "") -> None:
        self._expected = str(sequence or "").upper()
        if name:
            self._expected_name = name
        self._rerun()

    def _on_pick_expected(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Expected Sequence", str(self._folder),
            "Sequences (*.fa *.fasta *.fas *.fna *.txt *.seq);;All files (*)")
        if not path:
            return
        records = _read_fasta(Path(path))
        if not records:
            self.verdict.setText("That file holds no sequence.")
            return
        self._folder = Path(path).parent
        self.record_combo.blockSignals(True)
        self.record_combo.clear()
        for name, seq in records:
            self.record_combo.addItem(f"{name} ({len(seq):,} bp)", (name, seq))
        self.record_combo.setVisible(len(records) > 1)
        self.record_combo.setCurrentIndex(0)
        self.record_combo.blockSignals(False)
        name, seq = records[0]
        self.set_expected(seq, name)

    def _on_record(self, index: int) -> None:
        data = self.record_combo.itemData(index)
        if data:
            self.set_expected(data[1], data[0])

    def _on_add_reads(self) -> None:
        names = "Sanger reads (*.ab1 *.abi *.seq *.fa *.fasta *.txt);;All files (*)"
        paths, _ = QFileDialog.getOpenFileNames(self, "Sanger Reads", str(self._folder), names)
        if paths:
            self._folder = Path(paths[0]).parent
            self.add_reads([Path(p) for p in paths])

    def add_reads(self, paths) -> None:
        for p in paths:
            if Path(p) not in self._paths:
                self._paths.append(Path(p))
        self._rerun()

    def _rerun(self) -> None:
        self._checks = sg.check_reads(self._expected, self._paths) if self._paths else []
        self._refresh()

    # ------------------------------------------------------------------ view

    def _refresh(self) -> None:
        n = len(self._expected)
        self.expected_label.setText(
            f"<b>{html.escape(self._expected_name)}</b>, {n:,} bp · "
            f"{len(self._paths)} read{'s' if len(self._paths) != 1 else ''}"
            if n else "No expected sequence yet: press Expected Sequence.")
        self._fill_table()
        self._show_verdict()
        self.report_button.setEnabled(bool(self._checks) and bool(self._expected))

    def _show_verdict(self) -> None:
        if not self._expected or not self._checks:
            self.verdict.setText("")
            self.verdict.setStyleSheet("")
            return
        tone, line = sg.verdict(self._expected, self._checks)
        bg, border, fg = (("#E7F6EC", "#9BD3AE", "#14532D") if tone == "green"
                          else ("#FFF4E5", "#E8A33D", "#5A3700"))
        self.verdict.setStyleSheet(
            f"#sangerVerdict {{ background: {bg}; border: 1px solid {border};"
            f" border-radius: 6px; padding: 6px; color: {fg}; }}")
        self.verdict.setText(html.escape(line))

    def _fill_table(self) -> None:
        self._rows = []
        for check in self._checks:
            if check.differences:
                self._rows += [(check, d) for d in check.differences]
            else:
                self._rows.append((check, None))
        self.table.blockSignals(True)
        self.table.setRowCount(len(self._rows))
        for r, (check, diff) in enumerate(self._rows):
            if diff is None:
                cells = [check.name, check.summary(), "none", ""]
            else:
                cells = [check.name, check.summary(), diff.text().split(" (")[0],
                         "no quality" if diff.quality < 0 else f"Q{diff.quality}"]
            for c, text in enumerate(cells):
                item = QTableWidgetItem(text)
                if diff is not None and c >= 2:
                    item.setForeground(QColor("#9ca3af" if diff.doubtful else CHECK_COLOR))
                elif diff is None and c == 2:
                    item.setForeground(QColor(OK_COLOR))
                self.table.setItem(r, c, item)
        self.table.blockSignals(False)
        if self._rows:
            self.table.selectRow(0)
        else:
            self.trace.show_at(None, -1)

    def _on_row(self) -> None:
        rows = self.table.selectionModel().selectedRows() if self.table.selectionModel() else []
        if not rows or rows[0].row() >= len(self._rows):
            return
        check, diff = self._rows[rows[0].row()]
        if diff is None:
            middle = (check.trimmed[0] + check.trimmed[1]) // 2
            self.trace.show_at(check, middle)
        else:
            self.trace.show_at(check, diff.trace_index, diff.expected, diff.position)

    def _on_report(self) -> None:
        if not self._checks:
            return
        safe = "".join(c if c.isalnum() or c in "-_." else "_" for c in self._expected_name)[:40]
        out = self._folder / f"{safe or 'clone'}_sanger.html"
        n = 2
        while out.exists():
            out = self._folder / f"{safe or 'clone'}_sanger_{n}.html"
            n += 1
        try:
            sg.write_report(out, self._expected_name, self._expected, self._checks)
        except OSError as exc:
            self.verdict.setText(f"Could not write the report: {exc}")
            return
        self.report_path = out
        main = self.parent().window() if self.parent() else None
        show = getattr(main, "_show_new_output", None)
        if show is not None:
            show(out)


def _read_fasta(path: Path) -> list[tuple[str, str]]:
    """(name, sequence) of every record; a plain sequence file counts as one."""
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return []
    out: list[tuple[str, list[str]]] = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            out.append((line[1:].split()[0] or path.stem, []))
        elif out:
            out[-1][1].append(line)
        else:
            out.append((path.stem, [line]))
    return [(name, "".join(parts).upper()) for name, parts in out if any(parts)]
