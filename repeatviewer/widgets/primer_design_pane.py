"""Design tab of the primer workbench (Primers rail, beside the sequence).

Region from the map selection (or typed), Amplify / Inside / Avoid, the main
conditions, Design on the open sequence, candidates to tick and keep. Ladder,
FISH and multi-file design stay in the Design window ("Design Window…").

Developer: Abdullah Yusuf
"""

from __future__ import annotations

import json
import threading

from PySide6.QtCore import QObject, QSettings, Qt, Signal
from PySide6.QtWidgets import (
    QAbstractItemView,
    QButtonGroup,
    QComboBox,
    QDoubleSpinBox,
    QFrame,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QPushButton,
    QSizePolicy,
    QSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends.primers import design as dz
from repeatviewer.backends.primers import workbench as wb
from repeatviewer.widgets.help_button import make_help_button

_SETTINGS_KEY = "primers/workbench_design"
LABEL_W = 60
BOX_W = 60
JOIN_W = 26
GAP = 2

REGION_HELP = (
    "<b>Region</b>: drag on the map (or the bases) and the range appears here, or type it.<br>"
    "<b>Amplify</b>: primers sit outside the region and the product covers all of it.<br>"
    "<b>Inside</b>: both primers and the product stay inside the region.<br>"
    "<b>Avoid</b>: no primer may sit in the region (a SNP, a methylated site, a repeat).<br>"
    "When the region and the settings cannot work together, an amber note says why "
    "and offers the settings that fit. Nothing is blocked."
)
SETTINGS_HELP = (
    "<b>Mode</b>: Standard PCR or qPCR (short products). Ladder PCR and FISH probes, and "
    "design on many files, are in the Design window.<br>"
    "<b>Size</b>, <b>Tm</b>: min, optimum and max of each primer.<br>"
    "<b>Product</b>: shortest and longest product.<br>"
    "GC, tails and the salt conditions come from the Design window (Design Window…)."
)
CANDIDATES_HELP = (
    "Pairs Primer3 found around the region, best first. Click a pair to see its product on "
    "the sequence. Tick the pairs to keep and press <b>Keep Ticked</b>: they are added to "
    "this sequence's primers (Primer tab). Nothing is written before that.<br>"
    "<b>Open Report</b> writes the usual primer report (…_primers.html) for this run."
)


class _Relay(QObject):
    done = Signal(object, object)   # Primer3Run | None, error text | None


def _section(title: str, parent: QWidget, help_title: str, help_text: str):
    box = QWidget(parent)
    lay = QVBoxLayout(box)
    lay.setContentsMargins(0, 0, 0, 0)
    lay.setSpacing(6)
    head = QHBoxLayout()
    head.setContentsMargins(0, 0, 0, 0)
    lab = QLabel(f"<b>{title}</b>", box)
    head.addWidget(lab, 1)
    head.addWidget(make_help_button(box, help_title, help_text))
    lay.addLayout(head)
    return box, lay, lab


def _row(label: str, *widgets, parent: QWidget, joins=()) -> QHBoxLayout:
    row = QHBoxLayout()
    row.setContentsMargins(0, 0, 0, 0)
    row.setSpacing(GAP)
    lab = QLabel(label, parent)
    lab.setFixedWidth(LABEL_W)
    lab.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
    row.addWidget(lab)
    for i, w in enumerate(widgets):
        if i and i - 1 < len(joins) and joins[i - 1]:
            j = QLabel(joins[i - 1], parent)
            j.setEnabled(False)
            j.setFixedWidth(JOIN_W)
            j.setAlignment(Qt.AlignmentFlag.AlignCenter)
            row.addWidget(j)
        row.addWidget(w)
    row.addStretch(1)
    return row


def _spin(lo: float, hi: float, value: float, *, decimals: int = 0, suffix: str = "",
          parent: QWidget):
    box = QDoubleSpinBox(parent) if decimals else QSpinBox(parent)
    if decimals:
        box.setDecimals(decimals)
    box.setRange(lo, hi)
    box.setValue(value)
    if suffix:
        box.setSuffix(suffix)
    box.setFixedWidth(BOX_W)
    box.setKeyboardTracking(False)
    return box


class DesignPane(QWidget):
    """``keepRequested`` carries the Primer objects of the ticked pairs."""

    keepRequested = Signal(list)
    designWindowRequested = Signal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._host = None
        self._run: dz.Primer3Run | None = None
        self._run_region: wb.Region | None = None
        self._run_settings: dz.DesignSettings | None = None
        self._busy = False
        self._shown: tuple[int, int] | None = None   # range shown for a candidate
        self._relay = _Relay(self)
        self._relay.done.connect(self._on_done)

        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(10)
        s = self._load_settings()

        # ---- Region
        sec, sl, _ = _section("Region", self, "Region", REGION_HELP)
        self.region_from = _spin(1, 2_000_000_000, 1, parent=sec)
        self.region_to = _spin(1, 2_000_000_000, 1, parent=sec)
        self.region_from.setFixedWidth(BOX_W + 16)
        self.region_to.setFixedWidth(BOX_W + 16)
        sl.addLayout(_row("From", self.region_from, self.region_to, parent=sec, joins=("to",)))
        self.use_group = QButtonGroup(self)
        self.use_group.setExclusive(True)
        use_row = []
        for i, use in enumerate(wb.USES):
            b = QPushButton(wb.USE_LABELS[use], sec)
            b.setCheckable(True)
            b.setProperty("use", use)
            b.setMinimumWidth(0)
            b.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
            b.setFixedWidth((3 * BOX_W + 2 * JOIN_W + 4 * GAP - 2 * GAP) // 3)
            self.use_group.addButton(b, i)
            use_row.append(b)
        use_row[0].setChecked(True)
        sl.addLayout(_row("Use it to", *use_row, parent=sec))
        self.region_note = QLabel("Drag on the map, or type the range.", sec)
        self.region_note.setWordWrap(True)
        self.region_note.setStyleSheet("color: palette(placeholder-text);")
        self.region_note.setIndent(LABEL_W + GAP)
        sl.addWidget(self.region_note)
        self.advice = QFrame(sec)
        self.advice.setObjectName("workbenchAdvice")
        self.advice.setStyleSheet(
            "#workbenchAdvice { background: #FFF4E5; border: 1px solid #E8A33D;"
            " border-radius: 6px; } #workbenchAdvice QLabel { color: #5A3700; }")
        al = QVBoxLayout(self.advice)
        al.setContentsMargins(8, 6, 8, 6)
        al.setSpacing(6)
        self.advice_label = QLabel(self.advice)
        self.advice_label.setWordWrap(True)
        al.addWidget(self.advice_label)
        self.fit_button = QPushButton("Set Product Range to Fit", self.advice)
        fr = QHBoxLayout()
        fr.addStretch(1)
        fr.addWidget(self.fit_button)
        al.addLayout(fr)
        self.advice.hide()
        sl.addWidget(self.advice)
        lay.addWidget(sec)
        lay.addWidget(_rule(self))

        # ---- Settings
        sec, sl, _ = _section("Settings", self, "Settings", SETTINGS_HELP)
        self.mode = QComboBox(sec)
        for m in wb.WORKBENCH_MODES:
            self.mode.addItem(dz.MODE_LABELS[m], m)
        self.mode.setFixedWidth(3 * BOX_W + 2 * JOIN_W + 4 * GAP)
        i = self.mode.findData(s.mode)
        self.mode.setCurrentIndex(max(0, i))
        sl.addLayout(_row("Mode", self.mode, parent=sec))
        self.size_min = _spin(5, 60, s.size_min, parent=sec)
        self.size_opt = _spin(5, 60, s.size_opt, parent=sec)
        self.size_max = _spin(5, 60, s.size_max, parent=sec)
        sl.addLayout(_row("Size", self.size_min, self.size_opt, self.size_max, parent=sec,
                          joins=("opt", "max")))
        self.tm_min = _spin(20, 95, s.tm_min, decimals=1, parent=sec)
        self.tm_opt = _spin(20, 95, s.tm_opt, decimals=1, parent=sec)
        self.tm_max = _spin(20, 95, s.tm_max, decimals=1, parent=sec)
        sl.addLayout(_row("Tm (°C)", self.tm_min, self.tm_opt, self.tm_max, parent=sec,
                          joins=("opt", "max")))
        self.prod_min = _spin(10, 100_000, s.product_min, parent=sec)
        self.prod_max = _spin(10, 100_000, s.product_max, parent=sec)
        sl.addLayout(_row("Product", self.prod_min, self.prod_max, parent=sec, joins=("to",)))
        self.window_button = QPushButton("Design Window…", sec)
        self.window_button.setToolTip("GC, tails, salt, Ladder PCR, FISH probes, many files")
        self.defaults_button = QPushButton("Restore Defaults", sec)
        self.design_button = QPushButton("Design", sec)
        self.design_button.setDefault(True)
        br = QHBoxLayout()
        br.setContentsMargins(0, 2, 0, 0)
        br.addWidget(self.window_button)
        br.addStretch(1)
        br.addWidget(self.defaults_button)
        br.addWidget(self.design_button)
        sl.addLayout(br)
        lay.addWidget(sec)
        lay.addWidget(_rule(self))

        # ---- Candidates
        sec, sl, self.cand_title = _section("Candidates", self, "Candidates", CANDIDATES_HELP)
        self.table = QTableWidget(0, 4, sec)
        self.table.setHorizontalHeaderLabels(["", "Pair", "Product", "Tm F / R"])
        self.table.verticalHeader().hide()
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setShowGrid(False)
        self.table.setAlternatingRowColors(True)
        head = self.table.horizontalHeader()
        head.setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
        head.setDefaultAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        head.setStretchLastSection(True)
        for c, w in enumerate((24, 70, 70)):
            self.table.setColumnWidth(c, w)
        self.table.setMinimumHeight(150)
        sl.addWidget(self.table)
        self.cand_note = QLabel("None yet. Press Design.", sec)
        self.cand_note.setWordWrap(True)
        self.cand_note.setStyleSheet("color: palette(placeholder-text);")
        sl.addWidget(self.cand_note)
        self.report_button = QPushButton("Open Report", sec)
        self.keep_button = QPushButton("Keep Ticked", sec)
        kr = QHBoxLayout()
        kr.addStretch(1)
        kr.addWidget(self.report_button)
        kr.addWidget(self.keep_button)
        sl.addLayout(kr)
        lay.addWidget(sec)
        lay.addStretch(1)

        for w in (self.region_from, self.region_to, self.size_min, self.size_max,
                  self.prod_min, self.prod_max):
            w.valueChanged.connect(self._update_advice)
        self.use_group.idClicked.connect(lambda _i: self._update_advice())
        self.fit_button.clicked.connect(self._on_fit)
        self.design_button.clicked.connect(self.run_design)
        self.defaults_button.clicked.connect(self._on_defaults)
        self.window_button.clicked.connect(self.designWindowRequested.emit)
        self.keep_button.clicked.connect(self._on_keep)
        self.report_button.clicked.connect(self._on_report)
        self.table.itemSelectionChanged.connect(self._on_pick)
        self._sync_buttons()

    # ------------------------------------------------------------- host / state

    def set_host(self, host) -> None:
        self._host = host
        self.reset_for_sequence()

    def _sequence(self) -> tuple[str, str]:
        if self._host is None:
            return "", ""
        try:
            return self._host.primers_host_sequence()
        except Exception:  # noqa: BLE001 - no sequence loaded
            return "", ""

    def reset_for_sequence(self) -> None:
        _sid, seq = self._sequence()
        n = max(1, len(seq))
        for box in (self.region_from, self.region_to):
            box.blockSignals(True)
            box.setRange(1, n)
            box.blockSignals(False)
        self._clear_candidates()
        self._update_advice()

    def set_region(self, start: int, end: int) -> None:
        """From the map selection (1-based, inclusive). Single clicks are ignored."""
        a, b = sorted((int(start), int(end)))
        if b - a < 1 or (a, b) == self._shown:
            return   # a click, or the product we just showed for a candidate
        for box, v in ((self.region_from, a), (self.region_to, b)):
            box.blockSignals(True)
            box.setValue(v)
            box.blockSignals(False)
        self._update_advice()

    def region(self) -> wb.Region | None:
        a, b = self.region_from.value(), self.region_to.value()
        if b - a < 1:
            return None
        btn = self.use_group.checkedButton()
        use = btn.property("use") if btn is not None else wb.USE_AMPLIFY
        return wb.Region(a, b, use)

    def settings(self) -> dz.DesignSettings:
        base = self._load_settings()
        d = base.to_dict()
        d.update(mode=str(self.mode.currentData()),
                 size_min=self.size_min.value(), size_opt=self.size_opt.value(),
                 size_max=self.size_max.value(), tm_min=self.tm_min.value(),
                 tm_opt=self.tm_opt.value(), tm_max=self.tm_max.value(),
                 product_min=self.prod_min.value(), product_max=self.prod_max.value())
        return dz.DesignSettings.from_dict(d)

    def _load_settings(self) -> dz.DesignSettings:
        raw = QSettings().value(_SETTINGS_KEY, "")
        try:
            data = json.loads(raw) if raw else None
        except (TypeError, ValueError):
            data = None
        if not data:
            return dz.DesignSettings(**dz.mode_preset(dz.MODE_PCR))
        s = dz.DesignSettings.from_dict(data)
        if s.mode not in wb.WORKBENCH_MODES:
            s.mode = dz.MODE_PCR
        return s

    def _save_settings(self) -> None:
        QSettings().setValue(_SETTINGS_KEY, json.dumps(self.settings().to_dict()))

    def _on_defaults(self) -> None:
        p = dz.mode_preset(str(self.mode.currentData()))
        for box, key in ((self.size_min, "size_min"), (self.size_opt, "size_opt"),
                         (self.size_max, "size_max"), (self.tm_min, "tm_min"),
                         (self.tm_opt, "tm_opt"), (self.tm_max, "tm_max"),
                         (self.prod_min, "product_min"), (self.prod_max, "product_max")):
            box.setValue(p[key])
        self._update_advice()

    # ------------------------------------------------------------------ advice

    def _update_advice(self) -> None:
        _sid, seq = self._sequence()
        region = self.region()
        if region is None:
            self.region_note.setText("Drag on the map, or type the range.")
            self.advice.hide()
            return
        self.region_note.setText(region.label())
        adv = wb.region_advice(region, self.settings(), len(seq))
        if adv is None:
            self.advice.hide()
            return
        self.advice_label.setText(adv.message)
        self.fit_button.setVisible(adv.can_fix)
        self.advice.show()

    def _on_fit(self) -> None:
        _sid, seq = self._sequence()
        adv = wb.region_advice(self.region(), self.settings(), len(seq))
        if adv is None or not adv.can_fix:
            return
        self.prod_min.setValue(adv.product_min)
        self.prod_max.setValue(adv.product_max)
        self._update_advice()

    # ------------------------------------------------------------------ design

    def run_design(self) -> None:
        if self._busy:
            return
        sid, seq = self._sequence()
        if not seq:
            self.cand_note.setText("Open a sequence first.")
            return
        s = self.settings()
        region = self.region()
        self._save_settings()
        self._busy = True
        self._sync_buttons()
        self.cand_note.setText("Designing…")
        relay = self._relay

        def work() -> None:
            try:
                run = wb.design_candidates(sid, seq, s, region)
            except Exception as exc:  # noqa: BLE001 - shown in the pane
                relay.done.emit(None, str(exc))
                return
            relay.done.emit(run, None)

        self._run_region, self._run_settings = region, s
        threading.Thread(target=work, daemon=True).start()

    def _on_done(self, run, error) -> None:
        self._busy = False
        if error:
            self._clear_candidates()
            self.cand_note.setText(f"Design failed: {error}")
            self._sync_buttons()
            return
        self._run = run
        self._fill(run)
        self._sync_buttons()

    def _clear_candidates(self) -> None:
        self._run = None
        self.table.setRowCount(0)
        self.cand_title.setText("<b>Candidates</b>")
        self.cand_note.setText("None yet. Press Design.")
        self._sync_buttons()

    def _fill(self, run: dz.Primer3Run) -> None:
        pairs = run.pairs
        self.table.blockSignals(True)
        self.table.setRowCount(len(pairs))
        for r, p in enumerate(pairs):
            tick = QTableWidgetItem()
            tick.setFlags(Qt.ItemFlag.ItemIsUserCheckable | Qt.ItemFlag.ItemIsEnabled
                          | Qt.ItemFlag.ItemIsSelectable)
            tick.setCheckState(Qt.CheckState.Checked if r == 0 else Qt.CheckState.Unchecked)
            self.table.setItem(r, 0, tick)
            tms = (f"{p.forward.tm:.1f} / {p.reverse.tm:.1f}"
                   if p.forward is not None and p.reverse is not None else "")
            for c, v in enumerate((f"Pair {r + 1}", f"{p.product_size:,} bp", tms), start=1):
                self.table.setItem(r, c, QTableWidgetItem(v))
        self.table.blockSignals(False)
        n = len(pairs)
        self.cand_title.setText(f"<b>Candidates</b> &nbsp;<span style='color:gray'>"
                                f"{n} pair{'s' if n != 1 else ''}</span>")
        if pairs:
            self.cand_note.setText("Click a pair to see it on the sequence. Tick the ones to keep.")
            self.table.selectRow(0)
        else:
            why = "; ".join(run.explain[:2]) if run.explain else ""
            self.cand_note.setText("No pair met the conditions." + (f" {why}" if why else ""))

    def _on_pick(self) -> None:
        rows = self.table.selectionModel().selectedRows()
        if not rows or self._run is None or self._host is None:
            return
        p = self._run.pairs[rows[0].row()]
        if p.forward is None or p.reverse is None:
            return
        self._shown = (int(p.forward.start), int(p.reverse.end))
        try:
            self._host.primers_host_show(p.forward.start, p.reverse.end)
        except Exception:  # noqa: BLE001 - showing is a convenience
            pass
        notes = "; ".join(p.notes[:2])
        self.cand_note.setText(
            f"Pair {rows[0].row() + 1}: {p.forward.start:,} – {p.reverse.end:,}, "
            f"{p.product_size:,} bp." + (f" {notes}." if notes else ""))

    def ticked_pairs(self) -> list[dz.DesignedPair]:
        if self._run is None:
            return []
        out = []
        for r, p in enumerate(self._run.pairs):
            item = self.table.item(r, 0)
            if item is not None and item.checkState() == Qt.CheckState.Checked:
                out.append(p)
        return out

    def _on_keep(self) -> None:
        primers = []
        for p in self.ticked_pairs():
            primers.extend(wb.candidate_primers(p))
        if primers:
            self.keepRequested.emit(primers)

    def _on_report(self) -> None:
        if self._run is None or self._host is None:
            return
        sid, seq = self._sequence()
        try:
            path = self._host.primers_host_path()
            report = wb.write_candidates_report(path, sid, seq, self._run_settings or self.settings(),
                                                self._run_region, self._run)
        except Exception as exc:  # noqa: BLE001
            self.cand_note.setText(f"Could not write the report: {exc}")
            return
        from PySide6.QtCore import QUrl
        from PySide6.QtGui import QDesktopServices

        QDesktopServices.openUrl(QUrl.fromLocalFile(str(report)))

    def _sync_buttons(self) -> None:
        has = self._run is not None and bool(self._run.pairs)
        self.design_button.setEnabled(not self._busy)
        self.keep_button.setEnabled(has and not self._busy)
        self.report_button.setEnabled(self._run is not None and not self._busy)


def _rule(parent: QWidget) -> QFrame:
    line = QFrame(parent)
    line.setFrameShape(QFrame.Shape.HLine)
    line.setFrameShadow(QFrame.Shadow.Plain)
    line.setStyleSheet("color: palette(midlight);")
    return line


__all__ = ["DesignPane"]
