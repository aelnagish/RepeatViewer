"""Compare with Gel Photo: your real gel next to the in-silico one, band by band."""

from __future__ import annotations

import html
from pathlib import Path

import numpy as np
from PySide6.QtCore import QRectF, Qt
from PySide6.QtGui import QColor, QImage, QPainter, QPen
from PySide6.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QDialog,
    QFileDialog,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QMessageBox,
    QPushButton,
    QSizePolicy,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends.primers import gel_doc as gd
from repeatviewer.backends.primers import gel_photo as gp
from repeatviewer.widgets.help_button import make_help_button

PHOTO_HELP = (
    "<b>Open Photo</b>: a picture of your gel (PNG, JPG or TIFF; dark or light background).<br>"
    "Lanes and bands are found for you. Choose which photo lane holds the <b>ladder</b> (the "
    "gel's ladder, set in the gel window); its bands give the size scale.<br>"
    "For each photo lane choose the <b>gel lane</b> it should match. The table then names "
    "every band: <b>matched</b> (within 7 %), <b>off-target</b> (a predicted off-target "
    "product), <b>unexplained</b> (no product predicted: try the PCR with more mismatches), "
    "and expected bands <b>not seen</b>."
)

STATUS_COLORS = {"matched": "#1a7f4b", "off-target": "#b26b00", "unexplained": "#c2410c",
                 "not seen": "#8e8e93"}


def qimage_to_array(img: QImage) -> np.ndarray:
    g = img.convertToFormat(QImage.Format.Format_Grayscale8)
    w, h = g.width(), g.height()
    ptr = g.constBits()
    arr = np.frombuffer(ptr, dtype=np.uint8, count=g.bytesPerLine() * h).reshape(h, g.bytesPerLine())
    return arr[:, :w].astype(float)


class PhotoView(QWidget):
    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._img: QImage | None = None
        self._lanes: list[gp.Lane] = []
        self._ladder: int = -1
        self._ladder_ys: list[int] = []
        self._mark: tuple[int, int] | None = None   # lane, y
        self.setMinimumSize(360, 480)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

    def set_data(self, img: QImage | None, lanes, ladder: int, ladder_ys: list[int]) -> None:
        self._img, self._lanes, self._ladder, self._ladder_ys = img, list(lanes), ladder, list(ladder_ys)
        self.update()

    def mark(self, lane: int, y: int | None) -> None:
        self._mark = None if y is None else (lane, y)
        self.update()

    def _scale(self) -> tuple[float, float, float]:
        if self._img is None:
            return 1.0, 0.0, 0.0
        s = min(self.width() / self._img.width(), self.height() / self._img.height())
        ox = (self.width() - self._img.width() * s) / 2
        oy = (self.height() - self._img.height() * s) / 2
        return s, ox, oy

    def paintEvent(self, event) -> None:  # noqa: N802
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        if self._img is None:
            p.setPen(self.palette().color(self.foregroundRole()))
            p.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, "Open a photo of your gel")
            p.end()
            return
        s, ox, oy = self._scale()
        p.drawImage(QRectF(ox, oy, self._img.width() * s, self._img.height() * s), self._img)
        pen = QPen(QColor("#f5c542"), 1.5, Qt.PenStyle.DashLine)
        for i, ln in enumerate(self._lanes):
            r = QRectF(ox + (ln.x - ln.half) * s, oy + 2, 2 * ln.half * s, self._img.height() * s - 4)
            p.setPen(pen)
            p.setBrush(Qt.BrushStyle.NoBrush)
            p.drawRect(r)
            p.setPen(QColor("#f5c542"))
            p.drawText(QRectF(r.left(), oy + 4, r.width(), 16), Qt.AlignmentFlag.AlignHCenter,
                       "L" if i == self._ladder else str(i + 1))
        if 0 <= self._ladder < len(self._lanes):
            ln = self._lanes[self._ladder]
            p.setPen(Qt.PenStyle.NoPen)
            p.setBrush(QColor("#7fb0ff"))
            for y in self._ladder_ys:
                p.drawEllipse(QRectF(ox + (ln.x - ln.half) * s - 9, oy + y * s - 3, 6, 6))
        if self._mark is not None:
            lane, y = self._mark
            if 0 <= lane < len(self._lanes):
                ln = self._lanes[lane]
                p.setPen(QPen(QColor("#ff8a65"), 2))
                p.setBrush(Qt.BrushStyle.NoBrush)
                p.drawRoundedRect(QRectF(ox + (ln.x - ln.half) * s - 3, oy + y * s - 7,
                                         2 * ln.half * s + 6, 14), 3, 3)
        p.end()


class GelPhotoDialog(QDialog):
    def __init__(self, doc: gd.GelDoc, *, gel_path: Path | None = None,
                 parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Compare with Gel Photo")
        self.setMinimumSize(1180, 780)
        self._doc = doc
        self._gel_path = gel_path
        self._img: QImage | None = None
        self._arr: np.ndarray | None = None
        self._lanes: list[gp.Lane] = []
        self._bands: dict[int, list[gp.Band]] = {}
        self._fit: gp.LadderFit | None = None
        self._rows: list[tuple] = []

        root = QVBoxLayout(self)
        root.setContentsMargins(16, 12, 16, 12)
        root.setSpacing(10)
        top = QHBoxLayout()
        self.open_button = QPushButton("Open Photo…", self)
        top.addWidget(self.open_button)
        self.steps = QLabel(self)
        self.steps.setTextFormat(Qt.TextFormat.RichText)
        top.addWidget(self.steps, 1)
        top.addWidget(make_help_button(self, "Compare with Gel Photo", PHOTO_HELP))
        root.addLayout(top)

        body = QHBoxLayout()
        body.setSpacing(14)
        self.photo = PhotoView(self)
        body.addWidget(self.photo, 3)
        side = QVBoxLayout()
        side.setSpacing(8)
        row = QHBoxLayout()
        row.addWidget(QLabel("Ladder lane", self))
        self.ladder_combo = QComboBox(self)
        self.ladder_combo.setFixedWidth(120)
        row.addWidget(self.ladder_combo)
        row.addSpacing(8)
        row.addWidget(QLabel(f"Ladder: {html.escape(doc.ladder)}", self))
        row.addStretch(1)
        side.addLayout(row)
        side.addWidget(QLabel("<b>Photo lane → gel lane</b>", self))
        self.map_table = QTableWidget(0, 2, self)
        self.map_table.setHorizontalHeaderLabels(["Photo lane", "Gel lane"])
        self.map_table.verticalHeader().hide()
        self.map_table.setShowGrid(False)
        h = self.map_table.horizontalHeader()
        h.setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
        h.setDefaultAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        h.setStretchLastSection(True)
        self.map_table.setColumnWidth(0, 90)
        self.map_table.setFixedHeight(150)
        side.addWidget(self.map_table)
        self.fit_label = QLabel(self)
        self.fit_label.setWordWrap(True)
        side.addWidget(self.fit_label)
        side.addWidget(QLabel("<b>Bands</b>", self))
        self.table = QTableWidget(0, 4, self)
        self.table.setHorizontalHeaderLabels(["Lane", "Seen", "Matches", "What it is"])
        self.table.verticalHeader().hide()
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setShowGrid(False)
        h2 = self.table.horizontalHeader()
        h2.setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
        h2.setDefaultAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        h2.setStretchLastSection(True)
        for c, w in enumerate((60, 90, 90)):
            self.table.setColumnWidth(c, w)
        side.addWidget(self.table, 1)
        self.summary = QLabel(self)
        self.summary.setTextFormat(Qt.TextFormat.RichText)
        side.addWidget(self.summary)
        host = QWidget(self)
        host.setLayout(side)
        host.setFixedWidth(520)
        body.addWidget(host)
        root.addLayout(body, 1)
        bottom = QHBoxLayout()
        self.figure_button = QPushButton("Save Figure (Photo with Labels)…", self)
        self.table_button = QPushButton("Save Table (TSV)…", self)
        bottom.addWidget(self.figure_button)
        bottom.addWidget(self.table_button)
        bottom.addStretch(1)
        self.close_button = QPushButton("Close", self)
        bottom.addWidget(self.close_button)
        root.addLayout(bottom)

        self.open_button.clicked.connect(self._open)
        self.ladder_combo.currentIndexChanged.connect(lambda _i: self._analyse())
        self.table.itemSelectionChanged.connect(self._on_row)
        self.figure_button.clicked.connect(self._save_figure)
        self.table_button.clicked.connect(self._save_table)
        self.close_button.clicked.connect(self.close)
        self._sync_steps()

    # ---------------------------------------------------------------- steps

    def _sync_steps(self) -> None:
        done = '<span style="color:#1a7f4b">✓ {}</span>'
        todo = '<span style="color:#5f6368">{}</span>'
        cur = '<b style="color:#1f4e9a">{}</b>'
        has_img = self._img is not None
        has_fit = self._fit is not None
        parts = [
            (done if has_img else cur).format("1 Open the photo"),
            (done if has_img and self._lanes else (cur if has_img else todo)).format(
                f"2 Lanes found ({len(self._lanes)})" if self._lanes else "2 Lanes"),
            (done if has_fit else (cur if self._lanes else todo)).format(
                f"3 Ladder fitted ({len(self._fit.used)} bands)" if has_fit else "3 Ladder"),
            (cur if has_fit else todo).format("4 Check the bands"),
        ]
        self.steps.setText("&nbsp;&nbsp;›&nbsp;&nbsp;".join(parts))

    # ---------------------------------------------------------------- load

    def load_image(self, img: QImage) -> None:
        self._img = img
        self._arr = qimage_to_array(img)
        self._lanes = gp.find_lanes(self._arr)
        self._bands = {i: gp.find_bands(self._arr, ln) for i, ln in enumerate(self._lanes)}
        self.ladder_combo.blockSignals(True)
        self.ladder_combo.clear()
        for i in range(len(self._lanes)):
            self.ladder_combo.addItem(f"Lane {i + 1}", i)
        # ladder: the lane with the most bands
        if self._lanes:
            best = max(range(len(self._lanes)), key=lambda i: len(self._bands.get(i, [])))
            self.ladder_combo.setCurrentIndex(best)
        self.ladder_combo.blockSignals(False)
        self._fill_map()
        self._analyse()

    def _open(self) -> None:
        start = str(self._gel_path.parent) if self._gel_path else ""
        fn, _ = QFileDialog.getOpenFileName(self, "Open Gel Photo", start,
                                            "Images (*.png *.jpg *.jpeg *.tif *.tiff *.bmp)")
        if not fn:
            return
        img = QImage(fn)
        if img.isNull():
            QMessageBox.warning(self, "Gel Photo", "Could not read this image.")
            return
        self._photo_path = Path(fn)
        self.load_image(img)

    def _fill_map(self) -> None:
        sample_lanes = [(i, ln) for i, ln in enumerate(self._doc.lanes) if ln.kind != gd.LANE_LADDER]
        self.map_table.setRowCount(len(self._lanes))
        ladder = self.ladder_combo.currentData()
        k = 0
        for r in range(len(self._lanes)):
            self.map_table.setItem(r, 0, QTableWidgetItem(f"Lane {r + 1}"))
            combo = QComboBox(self.map_table)
            combo.addItem("none", -1)
            for i, ln in sample_lanes:
                combo.addItem(f"{i} · {ln.label}", i)
            if r != ladder and k < len(sample_lanes):
                combo.setCurrentIndex(k + 1)
                k += 1
            combo.currentIndexChanged.connect(lambda _i: self._analyse())
            self.map_table.setCellWidget(r, 1, combo)

    # ---------------------------------------------------------------- analyse

    def _analyse(self) -> None:
        if self._arr is None:
            self.photo.set_data(None, [], -1, [])
            self._sync_steps()
            return
        ladder = self.ladder_combo.currentData()
        ladder = int(ladder) if ladder is not None else -1
        sizes = [s for s, _m in gd.LADDERS.get(self._doc.ladder, [])]
        lad_ys = [b.y for b in self._bands.get(ladder, [])]
        self._fit = gp.fit_ladder(lad_ys, sizes) if ladder >= 0 else None
        self.photo.set_data(self._img, self._lanes, ladder, [y for y, _s in self._fit.used] if self._fit else lad_ys)
        self._rows = []
        if self._fit is None:
            self.fit_label.setText("The ladder lane needs at least 3 bands.")
        else:
            self.fit_label.setText(f"Ladder: {len(self._fit.used)} of {len(sizes)} bands used, "
                                   f"fit within {self._fit.max_error * 100:.0f} %.")
            for r in range(len(self._lanes)):
                if r == ladder:
                    continue
                combo = self.map_table.cellWidget(r, 1)
                gel_lane = int(combo.currentData()) if combo is not None else -1
                obs = [(self._fit.size_at(b.y), b.y) for b in self._bands.get(r, [])]
                expected = []
                if 0 <= gel_lane < len(self._doc.lanes):
                    for b in gd.lane_bands(self._doc.lanes[gel_lane], show_dimers=False):
                        expected.append((b.size, b.label, b.amount))
                calls = gp.compare([o for o, _y in obs], expected)
                by_size = {round(o): y for o, y in obs}
                for seen, exp, status, label in calls:
                    y = by_size.get(round(seen)) if seen is not None else None
                    self._rows.append((r, seen, exp, status, label, y))
        self._fill_table()
        self._sync_steps()

    def _fill_table(self) -> None:
        self.table.setRowCount(len(self._rows))
        counts: dict[str, int] = {}
        for i, (lane, seen, exp, status, label, _y) in enumerate(self._rows):
            counts[status] = counts.get(status, 0) + 1
            what = {"matched": label or "expected product",
                    "off-target": "predicted off-target product",
                    "unexplained": "no product predicted; try the PCR with more mismatches",
                    "not seen": f"expected {label}".strip() + ", not seen (faint or absent)"}[status]
            vals = [str(lane + 1), f"{seen:,.0f} bp" if seen else "–",
                    f"{exp:,} bp" if exp else "none", what]
            for c, v in enumerate(vals):
                item = QTableWidgetItem(v)
                if c in (2, 3):
                    item.setForeground(QColor(STATUS_COLORS.get(status, "#1d1d1f")))
                self.table.setItem(i, c, item)
        parts = [f'<span style="color:{STATUS_COLORS[k]}">● {n} {k}</span>'
                 for k, n in counts.items()]
        self.summary.setText(" &nbsp; ".join(parts))

    def _on_row(self) -> None:
        rows = sorted({i.row() for i in self.table.selectedIndexes()})
        if not rows:
            return
        lane, _seen, _exp, _status, _label, y = self._rows[rows[0]]
        self.photo.mark(lane, y)

    # ---------------------------------------------------------------- save

    def _save_figure(self) -> None:
        if self._img is None:
            return
        from repeatviewer.widgets.gel_view import save_view_image

        base = getattr(self, "_photo_path", None)
        start = str(base.with_name(base.stem + "_compared.png")) if base else "gel_compared.png"
        fn, _ = QFileDialog.getSaveFileName(self, "Save Figure", start, "PNG (*.png)")
        if fn:
            save_view_image(self.photo, Path(fn), "png")

    def _save_table(self) -> None:
        if not self._rows:
            return
        base = getattr(self, "_photo_path", None)
        start = str(base.with_name(base.stem + "_bands.tsv")) if base else "gel_bands.tsv"
        fn, _ = QFileDialog.getSaveFileName(self, "Save Table", start, "TSV (*.tsv)")
        if not fn:
            return
        with open(fn, "w", encoding="utf-8") as fh:
            fh.write("photo_lane\tseen_bp\texpected_bp\tstatus\tlabel\n")
            for lane, seen, exp, status, label, _y in self._rows:
                fh.write(f"{lane + 1}\t{'' if seen is None else round(seen)}\t{exp or ''}\t"
                         f"{status}\t{label}\n")
