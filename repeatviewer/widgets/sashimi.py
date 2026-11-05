"""Sashimi plot: read depth and splice junctions of RNA-seq alignment tracks.

One panel per alignment track over a region: depth as a filled profile, each
junction as an arc from intron start to end labelled with the reads that span
it (arcs alternate above and below the profile so the labels do not collide),
and the genes of the annotation tracks underneath."""

from __future__ import annotations

import numpy as np
from PySide6.QtCore import QPointF, QRectF, QSize, Qt
from PySide6.QtGui import QColor, QFont, QFontMetrics, QPainter, QPainterPath, QPen
from PySide6.QtWidgets import (
    QCheckBox,
    QDialog,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QScrollArea,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends import genome_tracks as gt

PANEL_H = 150
GENE_ROW_H = 18
LEFT_W = 120
RIGHT_PAD = 16
PALETTE = ("#c0504d", "#4f81bd", "#9bbb59", "#8064a2", "#f79646", "#4bacc6")


class SashimiView(QWidget):
    """The plot itself (also what Save image… renders)."""

    def __init__(self, ref_id: str, start: int, end: int, tracks: list, genes: list,
                 parent=None) -> None:
        super().__init__(parent)
        self.ref_id, self.start, self.end = ref_id, start, end
        self.tracks = list(tracks)
        self.genes = list(genes)
        self.min_junction = 1
        self.shared_scale = False
        self._data: dict[int, tuple] = {}
        self.setMinimumWidth(500)
        self._fit()

    # ---- data

    def _track_data(self, t, n: int):
        key = (id(t), n)
        got = self._data.get(key)
        if got is None:
            cov = t.coverage(self.ref_id, self.start, self.end, n)
            cov = np.zeros(n) if cov is None else np.nan_to_num(np.asarray(cov, float))
            j = t.junctions(self.ref_id, self.start, self.end)
            got = self._data[key] = (cov, j)
        return got

    def _gene_items(self):
        out = []
        for g in self.genes:
            try:
                out.extend(it for it, _r in g.fetch(self.ref_id, self.start, self.end)
                           if it.blocks or it.kind.lower() in ("gene", "mrna", "transcript"))
            except (gt.TrackError, OSError, ValueError, KeyError):
                continue
        # one row per transcript; plain genes only when there are no transcripts
        models = [it for it in out if it.blocks]
        return (models or out)[:12]

    def _fit(self) -> None:
        h = 24 + PANEL_H * max(1, len(self.tracks)) + 8
        h += GENE_ROW_H * len(self._gene_items()) + (16 if self.genes else 0)
        self.setMinimumHeight(h)
        self.resize(max(self.width(), 900), h)

    def sizeHint(self) -> QSize:  # noqa: N802
        return QSize(1000, self.minimumHeight())

    def set_min_junction(self, n: int) -> None:
        self.min_junction = max(1, int(n))
        self.update()

    def set_shared_scale(self, on: bool) -> None:
        self.shared_scale = bool(on)
        self.update()

    # ---- drawing

    def x_of(self, bp: float, w: float) -> float:
        span = max(1, self.end - self.start)
        return LEFT_W + (bp - self.start) * (w - LEFT_W - RIGHT_PAD) / span

    def paintEvent(self, _event) -> None:  # noqa: N802
        p = QPainter(self)
        self.render_to(p, self.width(), self.height())
        p.end()

    def render_to(self, p: QPainter, w: int, h: int) -> None:
        pal = self.palette()
        dark = pal.color(pal.ColorRole.Base).lightness() < 128
        p.fillRect(QRectF(0, 0, w, h), pal.color(pal.ColorRole.Base))
        fg = pal.color(pal.ColorRole.Text)
        mid = pal.color(pal.ColorRole.Mid)
        p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        font = QFont(self.font())
        small = QFont(font)
        small.setPointSize(max(7, font.pointSize() - 2))
        # title: the region
        p.setFont(font)
        p.setPen(fg)
        p.drawText(QRectF(LEFT_W, 2, w - LEFT_W, 18), Qt.AlignmentFlag.AlignVCenter,
                   f"{self.ref_id}:{self.start + 1:,}-{self.end:,}  "
                   f"({self.end - self.start:,} bp)")
        n = max(1, int(w - LEFT_W - RIGHT_PAD))
        data = [self._track_data(t, n) for t in self.tracks]
        peak_all = max([float(c.max()) for c, _j in data if len(c)] + [1.0])
        y = 24
        for k, (t, (cov, j)) in enumerate(zip(self.tracks, data)):
            color = QColor(PALETTE[k % len(PALETTE)])
            self._paint_panel(p, t, cov, j, y, w, color, fg, mid, small,
                              peak_all if self.shared_scale else None, dark)
            y += PANEL_H
        if self.genes:
            self._paint_genes(p, y + 8, w, fg, small)

    def _paint_panel(self, p, t, cov, j, top, w, color, fg, mid, small, peak, dark):
        fm = QFontMetrics(small)
        base = top + PANEL_H * 0.62         # depth above, lower arcs below
        cov_h = PANEL_H * 0.62 - 18
        peak = peak or max(1.0, float(cov.max()) if len(cov) else 1.0)
        p.setFont(small)
        p.setPen(fg)
        p.drawText(QRectF(6, top + 4, LEFT_W - 12, 14),
                   Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft,
                   fm.elidedText(t.name, Qt.TextElideMode.ElideMiddle, LEFT_W - 12))
        p.setPen(mid)
        p.drawText(QRectF(6, top + 18, LEFT_W - 12, 14),
                   Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft,
                   f"[0 – {peak:,.0f}]")
        # depth profile
        if len(cov):
            ys = base - np.clip(cov / peak, 0, 1) * cov_h
            path = QPainterPath(QPointF(LEFT_W, base))
            for i, yy in enumerate(ys.tolist()):
                path.lineTo(QPointF(LEFT_W + i, yy))
                path.lineTo(QPointF(LEFT_W + i + 1, yy))
            path.lineTo(QPointF(LEFT_W + len(ys), base))
            path.closeSubpath()
            fill = QColor(color)
            fill.setAlpha(170 if not dark else 200)
            p.setPen(Qt.PenStyle.NoPen)
            p.setBrush(fill)
            p.drawPath(path)
            p.setBrush(Qt.BrushStyle.NoBrush)
        p.setPen(QPen(mid, 1))
        p.drawLine(QPointF(LEFT_W, base), QPointF(w - RIGHT_PAD, base))
        # junctions, alternating above and below
        if not len(j):
            return
        depth = j[:, 2:].sum(axis=1)
        keep = depth >= self.min_junction
        j, depth = j[keep], depth[keep]
        if not len(j):
            return
        jpeak = max(1, int(depth.max()))
        order = np.argsort(j[:, 0], kind="stable")
        label_bg = QColor(self.palette().color(self.palette().ColorRole.Base))
        for rank, idx in enumerate(order.tolist()):
            s0, e0 = int(j[idx, 0]), int(j[idx, 1])
            d = int(depth[idx])
            x0, x1 = self.x_of(s0, w), self.x_of(e0, w)
            up = rank % 2 == 0
            # an arc leaves the profile at its height on each side (or the line)
            y0 = self._depth_y(cov, x0, base, cov_h, peak) if up else base
            y1 = self._depth_y(cov, x1, base, cov_h, peak) if up else base
            room = (min(y0, y1) - top - 30) if up else (top + PANEL_H - base - 6)
            hh = max(10.0, room * (0.45 + 0.55 * d / jpeak))
            apex = (min(y0, y1) - hh) if up else (base + hh)
            ctrl = (apex - min(y0, y1)) * 4 / 3 + min(y0, y1) if up else \
                base + hh * 4 / 3
            pen = QPen(color.darker(120), max(1.0, min(6.0, 1 + 5 * d / jpeak)))
            p.setPen(pen)
            path = QPainterPath(QPointF(x0, y0))
            path.cubicTo(QPointF(x0, ctrl), QPointF(x1, ctrl), QPointF(x1, y1))
            p.drawPath(path)
            text = f"{d:,}"
            tw = fm.horizontalAdvance(text) + 6
            r = QRectF((x0 + x1) / 2 - tw / 2, apex - 7, tw, 14)
            p.fillRect(r, label_bg)
            p.setPen(self.palette().color(self.palette().ColorRole.Text))
            p.drawText(r, Qt.AlignmentFlag.AlignCenter, text)

    def _depth_y(self, cov, x, base, cov_h, peak) -> float:
        i = int(x - LEFT_W)
        if not len(cov):
            return base
        i = min(max(i, 0), len(cov) - 1)
        lo, hi = max(0, i - 1), min(len(cov), i + 2)
        return base - min(1.0, float(cov[lo:hi].max()) / peak) * cov_h

    def _paint_genes(self, p, top, w, fg, small) -> None:
        fm = QFontMetrics(small)
        p.setFont(small)
        color = QColor("#2f4f8f")
        for r, it in enumerate(self._gene_items()):
            y = top + r * GENE_ROW_H
            mid_y = y + GENE_ROW_H / 2
            x0, x1 = self.x_of(it.start, w), self.x_of(it.end, w)
            p.setPen(fg)
            p.drawText(QRectF(6, y, LEFT_W - 12, GENE_ROW_H),
                       Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft,
                       fm.elidedText(it.name or it.kind, Qt.TextElideMode.ElideRight,
                                     LEFT_W - 12))
            p.setPen(QPen(color, 1))
            p.drawLine(QPointF(max(LEFT_W, x0), mid_y), QPointF(min(w - RIGHT_PAD, x1), mid_y))
            for s_, e_, _q in (it.blocks or [(it.start, it.end, 0)]):
                a, b = self.x_of(s_, w), self.x_of(e_, w)
                p.fillRect(QRectF(a, mid_y - 3, max(1.0, b - a), 6), color.lighter(130))
            for s_, e_ in it.thick:
                a, b = self.x_of(s_, w), self.x_of(e_, w)
                p.fillRect(QRectF(a, mid_y - 5, max(1.0, b - a), 10), color)


class SashimiDialog(QDialog):
    def __init__(self, parent, ref_id: str, start: int, end: int, tracks: list,
                 others: list, genes: list, *, min_junction: int = 1) -> None:
        super().__init__(parent)
        self.setWindowTitle(f"Sashimi plot – {ref_id}:{start + 1:,}-{end:,}")
        self.resize(1100, 700)
        self.view = SashimiView(ref_id, start, end, tracks, genes)
        self.view.set_min_junction(min_junction)
        lay = QVBoxLayout(self)
        bar = QHBoxLayout()
        bar.addWidget(QLabel("Junctions with at least"))
        self.min_box = QSpinBox()
        self.min_box.setRange(1, 100_000)
        self.min_box.setValue(max(1, min_junction))
        self.min_box.setSuffix(" reads")
        self.min_box.valueChanged.connect(self.view.set_min_junction)
        bar.addWidget(self.min_box)
        self.shared = QCheckBox("Same depth scale for all")
        self.shared.toggled.connect(self.view.set_shared_scale)
        bar.addWidget(self.shared)
        self._other_boxes = []
        for t in others:
            box = QCheckBox(t.name)
            box.setToolTip("Add this alignment track as a panel")
            box.toggled.connect(lambda on, t=t: self._toggle_track(t, on))
            bar.addWidget(box)
            self._other_boxes.append(box)
        bar.addStretch(1)
        save = QPushButton("Save image…")
        save.clicked.connect(self._save)
        bar.addWidget(save)
        lay.addLayout(bar)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(self.view)
        lay.addWidget(scroll, 1)

    def _toggle_track(self, t, on: bool) -> None:
        if on and t not in self.view.tracks:
            self.view.tracks.append(t)
        elif not on and t in self.view.tracks:
            self.view.tracks.remove(t)
        self.view._fit()
        self.view.update()

    def _save(self) -> None:
        from repeatviewer.widgets.image_export import save_widget_image

        save_widget_image(self, self.view, "sashimi",
                          lambda p, w, h: self.view.render_to(p, w, h))
