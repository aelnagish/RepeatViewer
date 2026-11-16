"""Digest page of the Enzymes tab: the ticked enzymes on a tandem array, as a gel.

Monomer ladders from copies that lost the site, long fragments where the
user's data shows the site methylated, partial cutting (Completeness). One
lane per ticked enzyme, and a "no data" lane for enzymes blocked in the data,
so the two can be compared. Open in a Gel writes a gel document; Compare with
Gel Photo puts a real gel next to it.

Developer: Abdullah Yusuf
"""

from __future__ import annotations

import json
from pathlib import Path

from PySide6.QtCore import QSettings, Qt, QTimer, Signal
from PySide6.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QDoubleSpinBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends.primers import array_digest as ad
from repeatviewer.backends.primers import gel_doc as gd
from repeatviewer.widgets.primer_design_pane import BOX_W, _row, _rule, _section, _spin
from repeatviewer.widgets.primer_pcr_pane import FIELD_W, PCR_LABEL_W

_SETTINGS_KEY = "primers/workbench_digest"
SHOWN_BANDS = 10

DIGEST_HELP = (
    "<b>Lanes</b>: one per enzyme ticked on the map page. With methylation data, enzymes "
    "with blocked sites get a second lane without the data, to compare.<br>"
    "<b>Site loss</b>: <i>From the copies</i> uses the sites as they are in this sequence; "
    "copies that lost the site by mutation give the 2×, 3× … ladder. <i>A site in every "
    "copy</i> adds the missing sites at the array period.<br>"
    "<b>Completeness</b>: share of molecules cut at an open site. Each site is cut with "
    "completeness × (1 − share methylated in a blocking context) from your data.<br>"
    "The array period comes from monomer annotations (a name with three or more copies), "
    "otherwise it is estimated from the sequence."
)
GEL_HELP = (
    "The expected gel, brightness by DNA mass. Very many faint bands of nearby sizes run "
    "as a smear. <b>Open in a Gel</b> writes a gel document you can change; <b>Compare "
    "with Gel Photo</b> puts your real gel next to this one."
)
LANE_HELP = (
    "The brightest bands of the chosen lane, by size. <b>n×</b>: n monomers. "
    "<b>Brightness</b>: share of the lane's DNA in that band."
)


class DigestPage(QWidget):
    """``pane`` is the EnzymePane (ticked rows, methylation, the PrimersPanel)."""

    backRequested = Signal()

    def __init__(self, pane, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._pane = pane
        self._lanes: list[ad.DigestLane] = []
        self._info = ad.ArrayInfo()
        self._doc: gd.GelDoc | None = None
        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(150)
        self._timer.timeout.connect(self.recompute)
        s = self._load_settings()

        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(10)
        back = QPushButton("‹ Enzymes", self)
        back.setFlat(True)
        back.clicked.connect(self.backRequested)
        br = QHBoxLayout()
        br.addWidget(back)
        br.addStretch(1)
        lay.addLayout(br)

        sec, sl, self.digest_title = _section("Digest", self, "Digest", DIGEST_HELP)
        self.lanes_label = QLabel(sec)
        self.lanes_label.setWordWrap(True)
        self.lanes_label.setFixedWidth(FIELD_W)
        sl.addLayout(_row("Lanes", self.lanes_label, parent=sec, label_w=PCR_LABEL_W))
        self.array_label = QLabel(sec)
        self.array_label.setWordWrap(True)
        self.array_label.setFixedWidth(FIELD_W)
        sl.addLayout(_row("Array", self.array_label, parent=sec, label_w=PCR_LABEL_W))
        self.loss = QComboBox(sec)
        self.loss.setFixedWidth(FIELD_W)
        for k in (ad.LOSS_COPIES, ad.LOSS_EVERY):
            self.loss.addItem(ad.LOSS_LABELS[k], k)
        self.loss.setCurrentIndex(max(0, self.loss.findData(s.get("site_loss", ad.LOSS_COPIES))))
        sl.addLayout(_row("Site loss", self.loss, parent=sec, label_w=PCR_LABEL_W))
        self.meth_label = QLabel(sec)
        self.meth_label.setWordWrap(True)
        self.meth_label.setFixedWidth(FIELD_W)
        sl.addLayout(_row("Methylation", self.meth_label, parent=sec, label_w=PCR_LABEL_W))
        self.completeness = _spin(10, 100, int(s.get("completeness", 95)), suffix=" %",
                                  parent=sec)
        sl.addLayout(_row("Completeness", self.completeness, parent=sec, label_w=PCR_LABEL_W))
        lay.addWidget(sec)
        lay.addWidget(_rule(self))

        sec, sl, _ = _section("Gel", self, "Gel", GEL_HELP)
        self.agarose = QDoubleSpinBox(sec)
        self.agarose.setRange(0.3, 4.0)
        self.agarose.setSingleStep(0.1)
        self.agarose.setDecimals(1)
        self.agarose.setSuffix(" %")
        self.agarose.setValue(float(s.get("agarose", 1.2)))
        self.agarose.setFixedWidth(BOX_W + 10)
        self.agarose.setKeyboardTracking(False)
        self.ladder = QComboBox(sec)
        self.ladder.addItem("Automatic", "")
        for name in gd.LADDERS:
            self.ladder.addItem(name, name)
        self.ladder.setCurrentIndex(max(0, self.ladder.findData(s.get("ladder", ""))))
        self.ladder.setFixedWidth(FIELD_W - BOX_W - 10 - 4)
        sl.addLayout(_row("Agarose", self.agarose, self.ladder, parent=sec, label_w=PCR_LABEL_W))
        from repeatviewer.widgets.gel_view import GelView

        self.gel = GelView(sec, compact=True)
        self.gel.setMinimumHeight(260)
        sl.addWidget(self.gel, 1)
        lay.addWidget(sec, 1)
        lay.addWidget(_rule(self))

        sec, sl, self.lane_title = _section("Lane", self, "Lane", LANE_HELP)
        self.lane_combo = QComboBox(sec)
        self.lane_combo.setFixedWidth(FIELD_W)
        sl.addLayout(_row("Lane", self.lane_combo, parent=sec, label_w=PCR_LABEL_W))
        self.table = QTableWidget(0, 3, sec)
        self.table.setHorizontalHeaderLabels(["", "Size", "Brightness"])
        self.table.verticalHeader().hide()
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.NoSelection)
        self.table.setShowGrid(False)
        self.table.setAlternatingRowColors(True)
        head = self.table.horizontalHeader()
        head.setDefaultAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        head.setStretchLastSection(True)
        for c, w in enumerate((60, 90)):
            self.table.setColumnWidth(c, w)
        self.table.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.table.setMinimumHeight(150)
        sl.addWidget(self.table)
        self.why = QLabel(sec)
        self.why.setWordWrap(True)
        self.why.setEnabled(False)
        sl.addWidget(self.why)
        self.photo_button = QPushButton("Compare with Gel Photo", sec)
        self.open_button = QPushButton("Open in a Gel", sec)
        row = QHBoxLayout()
        row.addStretch(1)
        row.addWidget(self.photo_button)
        row.addWidget(self.open_button)
        sl.addLayout(row)
        lay.addWidget(sec)

        self.loss.currentIndexChanged.connect(self.schedule)
        self.completeness.valueChanged.connect(self.schedule)
        self.agarose.valueChanged.connect(self.schedule)
        self.ladder.currentIndexChanged.connect(self.schedule)
        self.lane_combo.currentIndexChanged.connect(self._fill_lane)
        self.photo_button.clicked.connect(self._on_photo)
        self.open_button.clicked.connect(self._on_open)

    # ------------------------------------------------------------------ settings

    def _load_settings(self) -> dict:
        raw = QSettings().value(_SETTINGS_KEY, "")
        try:
            d = json.loads(raw) if raw else {}
        except (TypeError, ValueError):
            d = {}
        return d if isinstance(d, dict) else {}

    def settings(self) -> ad.DigestSettings:
        m = self._pane.settings()
        return ad.DigestSettings(site_loss=str(self.loss.currentData()),
                                 completeness=float(self.completeness.value()) / 100.0,
                                 min_percent=m.min_percent, min_reads=m.min_reads,
                                 agarose=float(self.agarose.value()),
                                 ladder=str(self.ladder.currentData() or ""))

    def _save_settings(self, s: ad.DigestSettings) -> None:
        QSettings().setValue(_SETTINGS_KEY, json.dumps(
            {"site_loss": s.site_loss, "completeness": round(s.completeness * 100),
             "agarose": s.agarose, "ladder": s.ladder}))

    # ------------------------------------------------------------------ compute

    def schedule(self, *_a) -> None:
        self._timer.start()

    def recompute(self) -> None:
        panel = self._pane._panel
        seq = panel._template
        rows = self._pane.ticked_rows()
        meth = self._pane._meth
        has_data = meth is not None and bool(meth.calls)
        s = self.settings()
        self._save_settings(s)
        features = []
        host = panel._host
        if host is not None:
            try:
                features = list(host.primers_host_features())
            except Exception:  # noqa: BLE001
                features = []
        self._info = ad.array_info(seq, features) if seq else ad.ArrayInfo()
        self._lanes = ad.digest_lanes(seq, rows, self._info, s, has_data=has_data) if seq else []
        self.lanes_label.setText(" · ".join(ln.label for ln in self._lanes) or
                                 "Tick enzymes on the map page.")
        self.array_label.setText(self._info.label())
        loss_label = ad.LOSS_LABELS[ad.LOSS_COPIES]
        if self._info.copies:
            loss_label += f" ({self._info.copies})"
        self.loss.setItemText(0, loss_label)
        if has_data:
            self.meth_label.setText(f"{Path(meth.path).name} (≥ {s.min_percent:.0f} %, "
                                    f"≥ {s.min_reads} reads)")
        else:
            self.meth_label.setText("None (choose a file on the map page)")
        seq_id = (panel._seq_id or "sequence").split()[0]
        self.digest_title.setText(f"<b>Digest</b> &nbsp;<span style='color:gray'>{seq_id}, "
                                  f"{len(self._lanes)} lanes</span>")
        self._doc = ad.gel_doc(self._lanes, len(seq), s, title=f"{seq_id} digest") \
            if self._lanes else None
        self.gel.set_doc(self._doc)
        old = self.lane_combo.currentText()
        self.lane_combo.blockSignals(True)
        self.lane_combo.clear()
        for ln in self._lanes:
            self.lane_combo.addItem(ln.label)
        i = self.lane_combo.findText(old)
        self.lane_combo.setCurrentIndex(i if i >= 0 else 0)
        self.lane_combo.blockSignals(False)
        self._fill_lane()
        for b in (self.photo_button, self.open_button):
            b.setEnabled(self._doc is not None)

    def lanes(self) -> list[ad.DigestLane]:
        return list(self._lanes)

    def _fill_lane(self, *_a) -> None:
        i = self.lane_combo.currentIndex()
        if not (0 <= i < len(self._lanes)):
            self.table.setRowCount(0)
            self.why.setText("")
            self.lane_title.setText("<b>Lane</b>")
            return
        lane = self._lanes[i]
        bright = sorted(lane.bands, key=lambda b: -b.mass_share)[:SHOWN_BANDS]
        bright.sort(key=lambda b: b.size)
        self.lane_title.setText(f"<b>{lane.label} lane</b> &nbsp;<span style='color:gray'>"
                                f"{len(lane.bands)} bands</span>")
        self.table.setRowCount(len(bright))
        for r, b in enumerate(bright):
            label = b.label
            if label == "uncut":
                label = "array"
            vals = [label, f"{b.size:,} bp", f"{b.mass_share * 100:.0f} %"
                    + (" (uncut array)" if b.label == "uncut" else "")]
            for c, v in enumerate(vals):
                self.table.setItem(r, c, QTableWidgetItem(v))
        self.why.setText("Why: " + lane.why)

    # ------------------------------------------------------------------ actions

    def _gel_path(self) -> Path | None:
        path = self._pane._panel.fasta_path()
        if path is None:
            return None
        seq_id = (self._pane._panel._seq_id or "sequence").split()[0]
        safe = "".join(c if c.isalnum() or c in "-_." else "_" for c in seq_id)[:40]
        return gd.unique_gel_path(path.parent, f"{path.stem}_{safe}_digest")

    def _on_open(self) -> None:
        from repeatviewer.widgets.primers_panel import find_main_window

        if self._doc is None:
            return
        out = self._gel_path()
        if out is None:
            return
        try:
            self._doc.save(out)
        except OSError as exc:
            self.why.setText(f"Could not write the gel: {exc}")
            return
        mw = find_main_window()
        if mw is not None and hasattr(mw, "_show_new_output"):
            mw._show_new_output(out)

    def _on_photo(self) -> None:
        if self._doc is None:
            return
        from repeatviewer.views.gel_photo_dialog import GelPhotoDialog

        dlg = GelPhotoDialog(self._doc, gel_path=None, parent=self.window())
        dlg.setWindowFlag(Qt.WindowType.Tool, True)
        dlg.show()
