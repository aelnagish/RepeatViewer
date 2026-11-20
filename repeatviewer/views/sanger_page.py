"""Sanger Check as a tab in the middle, with its settings on the right.

The sequence a clone should carry, the .ab1 reads, what differs, and the
chromatogram around each difference. No window: this is a page beside the
sequence tabs, and the panel on its right holds only its own settings.

Developer: Abdullah Yusuf
"""

from __future__ import annotations

import html
import json
from pathlib import Path

from PySide6.QtCore import QRectF, QSettings, Qt
from PySide6.QtGui import QColor, QPainter
from PySide6.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QComboBox,
    QFileDialog,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QPushButton,
    QSizePolicy,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends.primers import sanger as sg
from repeatviewer.views.sanger_dialog import SANGER_HELP, TRACE_HELP, TraceView, _read_fasta
from repeatviewer.widgets.help_button import make_help_button
from repeatviewer.widgets.primer_design_pane import (
    FIELD_W,
    LABEL_W,
    _row,
    _rule,
    _section,
    _spin,
)

_SETTINGS_KEY = "primers/sanger_page"
PANEL_W = 330
READ_COLORS = ("#1F57C3", "#C2410C", "#1a7f4b", "#7A3E9D")

COVERAGE_HELP = (
    "Which part of the expected sequence each read covers, and where the differences are. "
    "Amber marks what no read reached: the clone is only confirmed where a read went."
)


class CoverageBar(QWidget):
    """The expected sequence as a bar, one bar per read, differences marked."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._length = 0
        self._checks: list[sg.ReadCheck] = []
        self.setMinimumHeight(96)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.setToolTip("What each read covers, and where it differs.")

    def set_data(self, length: int, checks: list[sg.ReadCheck]) -> None:
        self._length = max(0, int(length))
        self._checks = [c for c in checks if c.end]
        self.setMinimumHeight(60 + 22 * max(1, len(self._checks)))
        self.update()

    def paintEvent(self, event) -> None:  # noqa: N802
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        n = self._length
        if n <= 0:
            painter.end()
            return
        left, right = 8.0, max(80.0, self.width() - 8.0)
        span = right - left
        font = painter.font()
        font.setPointSizeF(max(7.5, font.pointSizeF() * 0.82))
        painter.setFont(font)

        def x_of(pos: int) -> float:
            return left + span * (min(max(1, pos), n) - 1) / max(1, n - 1)

        painter.setPen(QColor("#5A5D63"))
        for pos in (1, n // 4, n // 2, 3 * n // 4, n):
            x = x_of(pos)
            align = (Qt.AlignmentFlag.AlignLeft if pos == 1
                     else Qt.AlignmentFlag.AlignRight if pos == n
                     else Qt.AlignmentFlag.AlignHCenter)
            painter.drawText(QRectF(x - 40, 0, 80, 14), align | Qt.AlignmentFlag.AlignVCenter,
                             f"{pos:,}")
        painter.fillRect(QRectF(left, 18, span, 10), QColor("#D9D9DE"))
        for a, b in sg.gaps_in_coverage(n, self._checks):
            painter.fillRect(QRectF(x_of(a), 18, max(2.0, x_of(b) - x_of(a)), 10),
                             QColor("#E8A33D"))
        gaps = sg.gaps_in_coverage(n, self._checks)
        if gaps:
            painter.setPen(QColor("#8A4B00"))
            first = ", ".join(f"{a:,} to {b:,}" for a, b in gaps[:2])
            painter.drawText(QRectF(left, 30, span, 14),
                             Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter,
                             f"not read: {first}")
        y = 48.0
        for i, check in enumerate(self._checks):
            colour = QColor(READ_COLORS[i % len(READ_COLORS)])
            x1, x2 = x_of(check.start), x_of(check.end)
            painter.fillRect(QRectF(x1, y, max(2.0, x2 - x1), 10), colour)
            painter.setPen(colour)
            arrow = "→" if check.strand == "+" else "←"
            text = f"{check.name} {arrow} {check.start:,} to {check.end:,}"
            width = painter.fontMetrics().horizontalAdvance(text) + 10
            if right - x2 >= width:                     # after the bar
                painter.drawText(QRectF(x2 + 6, y - 1, width, 12),
                                 Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter,
                                 text)
            elif x1 - left >= width:                    # before it
                painter.drawText(QRectF(x1 - width - 6, y - 1, width, 12),
                                 Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter,
                                 text)
            else:                                       # on it
                painter.setPen(QColor("#FFFFFF"))
                painter.drawText(QRectF(x1 + 6, y - 1, max(10.0, x2 - x1 - 12), 12),
                                 Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter,
                                 text)
            for d in check.differences:
                dx = x_of(d.position)
                painter.fillRect(QRectF(dx - 1, 16, 3, y - 16 + 12),
                                 QColor("#9CA3AF") if d.doubtful else QColor("#B26B00"))
            y += 22
        painter.end()


class SangerPage(QWidget):
    """One Sanger check: the reads of a clone against the sequence it should carry."""

    def __init__(self, expected: str = "", *, expected_name: str = "",
                 folder: Path | None = None, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.path = None                 # not a file page; the tab bar only needs a name
        self._expected = str(expected or "").upper()
        self._expected_name = expected_name or "expected sequence"
        self._folder = Path(folder) if folder else Path.cwd()
        self._checks: list[sg.ReadCheck] = []
        self._rows: list[tuple[sg.ReadCheck, sg.Difference | None]] = []
        self._paths: list[Path] = []
        saved = self._load()

        from repeatviewer.widgets.primer_design_pane import (
            _card, _results_surface, _scroll_panel)

        outer = QHBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)

        middle = QWidget(self)
        lay = QVBoxLayout(middle)
        lay.setContentsMargins(14, 12, 14, 12)
        lay.setSpacing(10)
        self.verdict = QLabel(middle)
        self.verdict.setWordWrap(True)
        self.verdict.setObjectName("sangerVerdict")
        lay.addWidget(self.verdict)

        self.table = QTableWidget(0, 6, middle)
        self.table.setHorizontalHeaderLabels(
            ["Read", "Direction", "Covers", "Identity", "Difference", "Quality"])
        self.table.verticalHeader().hide()
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setShowGrid(False)
        self.table.setAlternatingRowColors(True)
        head = self.table.horizontalHeader()
        head.setDefaultAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        head.setSectionResizeMode(4, QHeaderView.ResizeMode.Stretch)
        for c, w in ((0, 170), (1, 90), (2, 120), (3, 80), (5, 90)):
            self.table.setColumnWidth(c, w)
        self.table.setMinimumHeight(90)
        lay.addWidget(_card(self.table, middle))

        cov_head = QHBoxLayout()
        cov_head.addWidget(QLabel("<b>Coverage of the expected sequence</b>", middle), 1)
        cov_head.addWidget(make_help_button(middle, "Coverage", COVERAGE_HELP))
        lay.addLayout(cov_head)
        self.coverage = CoverageBar(middle)
        lay.addWidget(_card(self.coverage, middle))

        trace_head = QHBoxLayout()
        self.trace_title = QLabel("<b>Chromatogram</b>", middle)
        trace_head.addWidget(self.trace_title, 1)
        trace_head.addWidget(make_help_button(middle, "Chromatogram", TRACE_HELP))
        lay.addLayout(trace_head)
        self.trace = TraceView(middle)
        lay.addWidget(_card(self.trace, middle), 1)
        _results_surface(middle)
        outer.addWidget(middle, 1)

        panel = QWidget(self)
        panel.setFixedWidth(PANEL_W)
        panel.setObjectName("sangerPanel")
        pl = QVBoxLayout(panel)
        pl.setContentsMargins(10, 12, 10, 12)
        pl.setSpacing(10)
        title = QHBoxLayout()
        title.addWidget(QLabel("<b>Sanger Check</b>", panel), 1)
        title.addWidget(make_help_button(panel, "Sanger Check", SANGER_HELP))
        pl.addLayout(title)
        pl.addWidget(_rule(panel))

        sec, sl, _ = _section("Expected", panel, "Expected", SANGER_HELP)
        self.expected_button = QPushButton("Choose", sec)
        self.expected_button.setFixedWidth(FIELD_W)
        sl.addLayout(_row("Sequence", self.expected_button, parent=sec))
        self.record_combo = QComboBox(sec)
        self.record_combo.setFixedWidth(FIELD_W)
        sl.addLayout(_row("Record", self.record_combo, parent=sec))
        pl.addWidget(sec)
        pl.addWidget(_rule(panel))

        sec, sl, _ = _section("Reads", panel, "Reads", SANGER_HELP)
        self.add_button = QPushButton("Add Reads", sec)
        self.add_button.setFixedWidth(FIELD_W)
        sl.addLayout(_row("Files", self.add_button, parent=sec))
        self.reads_note = QLabel("None yet.", sec)
        self.reads_note.setWordWrap(True)
        self.reads_note.setIndent(LABEL_W + 2)
        self.reads_note.setStyleSheet("color: palette(placeholder-text);")
        sl.addWidget(self.reads_note)
        pl.addWidget(sec)
        pl.addWidget(_rule(panel))

        sec, sl, _ = _section("Quality", panel, "Quality", SANGER_HELP)
        self.cutoff = _spin(0.001, 0.5, float(saved.get("cutoff", sg.TRIM_CUTOFF)), decimals=3,
                            parent=sec)
        sl.addLayout(_row("Trim cutoff", self.cutoff, parent=sec))
        self.good_q = _spin(0, 60, int(saved.get("good_quality", sg.GOOD_QUALITY)), parent=sec)
        sl.addLayout(_row("Real from Q", self.good_q, parent=sec))
        pl.addWidget(sec)
        pl.addWidget(_rule(panel))

        sec, sl, _ = _section("Options", panel, "Options", SANGER_HELP)
        self.report_check = QCheckBox("Write a report beside the reads", sec)
        self.report_check.setChecked(bool(saved.get("report", False)))
        sl.addLayout(_row("", self.report_check, parent=sec))
        self.defaults_button = QPushButton("Restore Defaults", sec)
        self.save_button = QPushButton("Save", sec)
        self.check_button = QPushButton("Check", sec)
        self.check_button.setDefault(True)
        br = QHBoxLayout()
        br.addWidget(self.defaults_button)
        br.addStretch(1)
        br.addWidget(self.save_button)
        br.addWidget(self.check_button)
        sl.addLayout(br)
        self.saved_note = QLabel("", sec)
        self.saved_note.setStyleSheet("color: palette(placeholder-text);")
        sl.addWidget(self.saved_note)
        pl.addWidget(sec)
        pl.addStretch(1)
        outer.addWidget(_scroll_panel(panel, self))

        self.expected_button.clicked.connect(self._on_pick_expected)
        self.record_combo.currentIndexChanged.connect(self._on_record)
        self.add_button.clicked.connect(self._on_add_reads)
        self.table.itemSelectionChanged.connect(self._on_row)
        self.check_button.clicked.connect(self.recheck)
        self.save_button.clicked.connect(self._on_save)
        self.defaults_button.clicked.connect(self._on_defaults)
        self._refresh()

    # ------------------------------------------------------------------ settings

    def _load(self) -> dict:
        raw = QSettings().value(_SETTINGS_KEY, "")
        try:
            d = json.loads(raw) if raw else {}
        except (TypeError, ValueError):
            d = {}
        return d if isinstance(d, dict) else {}

    def _on_save(self) -> None:
        QSettings().setValue(_SETTINGS_KEY, json.dumps(
            {"cutoff": float(self.cutoff.value()), "good_quality": int(self.good_q.value()),
             "report": self.report_check.isChecked()}))
        self.saved_note.setText("Saved.")

    def _on_defaults(self) -> None:
        self.cutoff.setValue(sg.TRIM_CUTOFF)
        self.good_q.setValue(sg.GOOD_QUALITY)
        self.report_check.setChecked(False)
        self.recheck()

    # ------------------------------------------------------------------ inputs

    def tab_name(self) -> str:
        return f"Sanger Check · {self._expected_name}"

    def set_expected(self, sequence: str, name: str = "") -> None:
        self._expected = str(sequence or "").upper()
        if name:
            self._expected_name = name
        self.recheck()

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
        self.record_combo.setCurrentIndex(0)
        self.record_combo.blockSignals(False)
        self.set_expected(records[0][1], records[0][0])

    def _on_record(self, index: int) -> None:
        data = self.record_combo.itemData(index)
        if data:
            self.set_expected(data[1], data[0])

    def _on_add_reads(self) -> None:
        paths, _ = QFileDialog.getOpenFileNames(
            self, "Sanger Reads", str(self._folder),
            "Sanger reads (*.ab1 *.abi *.seq *.fa *.fasta *.txt);;All files (*)")
        if paths:
            self._folder = Path(paths[0]).parent
            self.add_reads([Path(p) for p in paths])

    def add_reads(self, paths) -> None:
        for p in paths:
            if Path(p) not in self._paths:
                self._paths.append(Path(p))
        self.recheck()

    # ------------------------------------------------------------------ run

    def recheck(self) -> None:
        self._checks = (sg.check_reads(self._expected, self._paths,
                                       cutoff=float(self.cutoff.value()))
                        if self._paths and self._expected else [])
        for check in self._checks:
            for d in check.differences:
                d.quality = d.quality          # kept as read; the threshold is applied below
        self._refresh()
        if self.report_check.isChecked() and self._checks:
            self._write_report()

    def _doubtful(self, d: sg.Difference) -> bool:
        return 0 <= d.quality < int(self.good_q.value())

    def _refresh(self) -> None:
        n = len(self._expected)
        self.reads_note.setText(
            f"{len(self._paths)} read{'s' if len(self._paths) != 1 else ''}: "
            + ", ".join(p.stem for p in self._paths[:3])
            + (" …" if len(self._paths) > 3 else "") if self._paths else "None yet.")
        self._fill_table()
        self.coverage.set_data(n, self._checks)
        self._show_verdict()

    def _show_verdict(self) -> None:
        if not self._expected:
            self.verdict.setText("Choose the sequence this clone should carry.")
            self.verdict.setStyleSheet("")
            return
        if not self._checks:
            self.verdict.setText(f"{self._expected_name}, {len(self._expected):,} bp. "
                                 "Add the .ab1 reads.")
            self.verdict.setStyleSheet("")
            return
        tone, line = sg.verdict(self._expected, self._checks)
        bg, border, fg = (("#E7F6EC", "#9BD3AE", "#14532D") if tone == "green"
                          else ("#FFF4E5", "#E8A33D", "#5A3700"))
        self.verdict.setStyleSheet(
            f"#sangerVerdict {{ background: {bg}; border: 1px solid {border};"
            f" border-radius: 6px; padding: 8px; color: {fg}; }}")
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
            way = "forward" if check.strand == "+" else "reverse"
            covers = f"{check.start:,} to {check.end:,}" if check.end else "no match"
            cells = [check.name, way, covers, f"{check.identity * 100:.1f} %",
                     "none" if diff is None else diff.text().split(" (")[0],
                     "" if diff is None else
                     ("no quality" if diff.quality < 0 else f"Q{diff.quality}")]
            for c, text in enumerate(cells):
                item = QTableWidgetItem(text)
                if c == 4:
                    item.setForeground(QColor("#1a7f4b") if diff is None else
                                       QColor("#9CA3AF") if self._doubtful(diff) else
                                       QColor("#B26B00"))
                self.table.setItem(r, c, item)
        self.table.blockSignals(False)
        rows = sum(self.table.rowHeight(r) for r in range(self.table.rowCount()))
        self.table.setMaximumHeight(self.table.horizontalHeader().height() + rows + 6)
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
            self.trace_title.setText(f"<b>Chromatogram</b> &nbsp;<span style='color:gray'>"
                                     f"{check.name}, middle of the read</span>")
        else:
            self.trace.show_at(check, diff.trace_index, diff.expected, diff.position)
            self.trace_title.setText(
                f"<b>Chromatogram</b> &nbsp;<span style='color:gray'>{check.name}, "
                f"base {diff.trace_index + 1:,} · expected {diff.position:,}</span>")

    def _write_report(self) -> None:
        safe = "".join(c if c.isalnum() or c in "-_." else "_" for c in self._expected_name)[:40]
        out = self._folder / f"{safe or 'clone'}_sanger.html"
        n = 2
        while out.exists():
            out = self._folder / f"{safe or 'clone'}_sanger_{n}.html"
            n += 1
        try:
            sg.write_report(out, self._expected_name, self._expected, self._checks)
        except OSError as exc:
            self.saved_note.setText(f"Could not write the report: {exc}")
            return
        self.report_path = out
        self.saved_note.setText(f"Report: {out.name}")
