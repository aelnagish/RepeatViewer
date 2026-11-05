"""Drawing of genome-browser tracks (annotations, PAF / BLAST / SAM / BAM) in the
contig view. ``GenomeTracksMixin`` is mixed into ``ContigCanvas``: it only uses
the canvas geometry (``x_of``, ``bp_at``, ``ppb``, ``label_w``) and its model.

Data come from ``backends/genome_tracks.py``; nothing here reads files.
"""

from __future__ import annotations

import numpy as np
from PySide6.QtCore import QLineF, QPointF, QRectF, Qt
from PySide6.QtGui import QColor, QFont, QFontMetrics, QPainter, QPen, QPolygonF
from PySide6.QtWidgets import QWidget

from repeatviewer.backends import alignment_style as style
from repeatviewer.backends import genome_tracks as gt

ANN_ROW_H = 14
ANN_SQUISH_H = 6
ALN_COV_H = 26
ALN_ROW_H = 10
ALN_SQUISH_H = 4
ALN_COLLAPSE_H = 2          # collapsed: 2 px lines, no gap, a density picture
NOTE_H = 14
TRACK_GAP = 6
DEFAULT_ANN_ROWS = 8
DEFAULT_ALN_ROWS = 20
MAX_DRAWN = 30_000        # alignments in view above this: coverage only
MISMATCH_MIN_PPB = 0.5    # draw mismatches from 2 bp per pixel

GROUP_HEAD_H = 14         # band with a group's name, above its rows

MODE_EXPANDED = "expanded"
MODE_SQUISHED = "squished"
MODE_COLLAPSED = "collapsed"      # reads as thin touching lines (alignment tracks)

_SECONDARY_ALPHA = 90
_COV_FILL = QColor("#a9a9a9")
_INS = QColor("#8e44ad")
_CLIP = QColor("#d6452f")
_MATE_LINE = QColor("#9aa0a8")
_INTRON_LINE = QColor(110, 118, 130, 170)   # spliced gaps: one neutral colour (IGV-like)
_SEL_FILL = QColor(31, 95, 191, 22)
_SEL_EDGE = QColor("#1f5fbf")
DEFAULT_CLIP_MIN = 20     # mark soft clips of at least this many bases


def _pointed_bar(x0: float, x1: float, y: float, h: float, rev: bool,
                 tip: float | None = None) -> QPolygonF:
    """A bar from x0 to x1 whose strand end is an arrow point, kept inside x0..x1 so
    it never reaches into the next item on the row."""
    d = tip if tip is not None else min(max(3.0, h * 0.5), (x1 - x0) * 0.5)
    m = y + h / 2
    if rev:
        pts = [(x0, m), (x0 + d, y), (x1, y), (x1, y + h), (x0 + d, y + h)]
    else:
        pts = [(x0, y), (x1 - d, y), (x1, m), (x1 - d, y + h), (x0, y + h)]
    return QPolygonF([QPointF(a, b) for a, b in pts])


def _clear_of_chevrons(left: float, right: float, first: float,
                       step: float = 14) -> tuple[float, float]:
    """Widen a label patch so its edges fall between strand chevrons (drawn every
    ``step`` px from ``first``), never through one."""
    xc = first + step * round((left - first) / step)
    if xc - 3.5 < left < xc + 3.5:
        left = xc - 3.5
    xc = first + step * round((right - first) / step)
    if xc - 3.5 < right < xc + 3.5:
        right = xc + 3.5
    return left, right


def _edge_color(fill: QColor) -> QColor:
    """Outline of a read bar: the fill, darker and a little more saturated."""
    h, s_, v, _a = fill.getHsvF()
    return QColor.fromHsvF(max(0.0, h), min(1.0, s_ * 1.5 + 0.1), v * 0.72)


def _type_color(kind: str, custom: dict[str, str]) -> QColor:
    if kind in custom:
        return QColor(custom[kind])
    from repeatviewer.backends.gff_features import _feature_color

    return QColor(_feature_color(kind or "feature"))


class GenomeTracksMixin:
    """Extra tracks stacked under the reference, above the placed sequences."""

    def _init_tracks(self) -> None:
        self.tracks: list[gt.Track] = []
        self._extra_spans: list[tuple[int, int]] = []
        self._ruler_drag: tuple[float, float] | None = None   # (x0, x1) while selecting
        self._body_press: tuple[float, int, bool] | None = None  # (x, scroll, moved)
        self.region_sel: tuple[int, int] | None = None   # ruler selection [a, b) in bp
        self._sel_chip = None
        self._track_drag: int | None = None                   # track being moved
        self._track_drop: int | None = None                   # where it would land

    # ---------------------------------------------------------------- tracks

    def set_tracks(self, tracks: list) -> None:
        self.tracks = list(tracks)
        self._tracks_changed()

    def add_track(self, track, above_reads: bool = False) -> None:
        """Adds a track last; with ``above_reads`` (a file the user adds) an annotation
        track goes above the read tracks, as in IGV, so it never has to be dragged past
        a tall alignment track. Sessions load with their saved order."""
        at = len(self.tracks)
        if above_reads and track.kind == gt.KIND_ANNOTATION:
            at = next((i for i, t in enumerate(self.tracks)
                       if t.kind != gt.KIND_ANNOTATION), at)
        self.tracks.insert(at, track)
        self._tracks_changed()

    def remove_track(self, index: int) -> None:
        if 0 <= index < len(self.tracks):
            del self.tracks[index]
            self._tracks_changed()

    def set_track_mode(self, index: int, mode: str) -> None:
        if 0 <= index < len(self.tracks):
            self.tracks[index].settings["mode"] = mode
            self._tracks_changed()

    def move_track(self, src: int, dst: int) -> None:
        """Move track ``src`` so it ends up at position ``dst`` (0 = top)."""
        n = len(self.tracks)
        if not (0 <= src < n) or src == dst:
            return
        t = self.tracks.pop(src)
        self.tracks.insert(max(0, min(dst, n - 1)), t)
        self._tracks_changed()
        self.tracksReordered.emit()

    def _tracks_changed(self) -> None:
        self._tiles_epoch = getattr(self, "_tiles_epoch", 0) + 1
        if getattr(self, "_tile_cache", None):
            self._tile_cache.clear()
        self._update_scrollbars()
        self.viewport().update()
        self.viewChanged.emit()

    # ---------------------------------------------------------------- geometry

    @staticmethod
    def _mode(t) -> str:
        return t.settings.get("mode", MODE_EXPANDED)

    def _row_h(self, t) -> int:
        squished = self._mode(t) == MODE_SQUISHED
        if t.kind == gt.KIND_ANNOTATION:
            return ANN_SQUISH_H if squished else ANN_ROW_H
        if self._mode(t) == MODE_COLLAPSED:
            return ALN_COLLAPSE_H
        return ALN_SQUISH_H if squished else ALN_ROW_H

    def _max_rows(self, t) -> int:
        default = DEFAULT_ANN_ROWS if t.kind == gt.KIND_ANNOTATION else DEFAULT_ALN_ROWS
        n = int(t.settings.get("max_rows", default))
        mode = self._mode(t)
        if mode == MODE_COLLAPSED and t.kind != gt.KIND_ANNOTATION:
            return n * 4
        return n * 3 if mode == MODE_SQUISHED else n

    def _shown_rows(self, t) -> int:
        m = self.model
        if m is None:
            return 1
        if isinstance(t, gt.IndexedBamTrack):
            return self._max_rows(t)
        return max(1, min(self._max_rows(t), t.row_count(m.ref_id)))

    def _extra_track_h(self, t) -> int:
        if t.kind == gt.KIND_ANNOTATION:
            return self._shown_rows(t) * self._row_h(t) + TRACK_GAP
        return ALN_COV_H + 2 + self._layout(t)[1] + TRACK_GAP

    # ---------------------------------------------------------------- groups and sorting

    def _sync_grouping(self, t) -> None:
        """Apply the track's filter, "View as pairs" and "Group by" settings to its rows."""
        if t.kind != gt.KIND_ANNOTATION:
            t.set_filter(*style.read_filter(t.settings))
            t.set_pairs(bool(t.settings.get("pairs", False)))
        mode = t.settings.get("group_by", style.GROUP_NONE)
        tag_name = t.settings.get("group_tag", "")
        chrom = self.model.ref_id if self.model is not None else ""
        fn = style.group_function(mode, tag_name=tag_name, chrom=chrom,
                                  read_groups=getattr(t, "read_groups", {}))
        t.set_grouping((mode, tag_name, chrom), fn)

    def _layout(self, t):
        """(groups as [(label, top, rows shown)], height of the rows) of an alignment
        track; one unnamed group when the reads are not grouped."""
        rh = self._row_h(t)
        if t.kind == gt.KIND_ANNOTATION or self.model is None:
            return [], self._shown_rows(t) * rh
        self._sync_grouping(t)
        if not t.grouped:
            n = self._shown_rows(t)
            return [(None, 0, n)], n * rh
        counts = t.group_rows(self.model.ref_id)
        labels = sorted(counts, key=lambda g: (g.startswith(("no ", "unpaired")), g))
        if not labels:
            return [], rh
        cap = max(3, self._max_rows(t) // len(labels))
        fixed = isinstance(t, gt.IndexedBamTrack)    # rows grow while panning
        groups, y = [], 0
        for g in labels:
            shown = cap if fixed else max(1, min(cap, counts[g]))
            groups.append((g, y, shown))
            y += GROUP_HEAD_H + shown * rh
        return groups, y

    def _placed(self, t, a: int, b: int):
        """(item, y offset under the coverage or None when its row is not shown) of the
        alignments in [a, b], in groups and sorted as the track settings say."""
        m = self.model
        rh = self._row_h(t)
        groups, _h = self._layout(t)
        perm = self._sort_perm(t)
        if not t.grouped:
            shown = groups[0][2] if groups else 0
            out = []
            for it, row in t.fetch(m.ref_id, a, b):
                row = perm.get((None, row), row) if perm else row
                out.append((it, row * rh if row < shown else None))
            return out
        where = {g: (top + GROUP_HEAD_H, shown) for g, top, shown in groups}
        out = []
        for it, g, row in t.fetch_groups(m.ref_id, a, b):
            if perm:
                row = perm.get((g, row), row)
            top, shown = where.get(g, (0, 0))
            out.append((it, top + row * rh if row < shown else None))
        return out

    def _paint_group_heads(self, p: QPainter, t, rows_top: float, h: int, w: int) -> None:
        """Name band of each group (kept at the left edge of the data area), and a
        relayout when reading the file brought new groups."""
        groups, rows_h = self._layout(t)
        if ALN_COV_H + 2 + rows_h + TRACK_GAP != h:
            from PySide6.QtCore import QTimer

            def relayout():
                self._update_scrollbars()
                self.viewport().update()
            QTimer.singleShot(0, relayout)
        pal = self.palette()
        band = QColor(pal.color(pal.ColorRole.AlternateBase))
        text = pal.color(pal.ColorRole.Text)
        font = QFont(p.font())
        font.setPointSize(8)
        font.setBold(True)
        fm = QFontMetrics(font)
        p.setFont(font)
        counts = t.group_rows(self.model.ref_id)
        for g, top, shown in groups:
            y = rows_top + top
            p.fillRect(QRectF(self.label_w, y, w - self.label_w, GROUP_HEAD_H - 2), band)
            label = str(g)
            if counts.get(g, 0) > shown:
                # reads stack into rows so none overlap; only the first rows are drawn
                label += (f"  ·  showing {shown:,} of {counts[g]:,} rows of stacked reads"
                          "  ·  raise Rows in the side panel to see more")
            p.setPen(text)
            p.drawText(QRectF(self.label_w + 6, y, fm.horizontalAdvance(label) + 8,
                              GROUP_HEAD_H - 2),
                       Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft, label)

    def _sort_perm(self, t) -> dict:
        """(group, row) -> shown row: rows whose read covers the sort position come
        first, in the order of the sort key (the rest keep their order below)."""
        spec = t.settings.get("sort")
        m = self.model
        if not spec or m is None or spec.get("ref") != m.ref_id:
            return {}
        pos = int(spec.get("pos", 0))
        region = None
        if isinstance(t, gt.IndexedBamTrack):
            # only the loaded reads: asking the file for pos could move the region
            region = t._region
            if region is None or region.key != t.resolve(m.ref_id) or \
                    not region.lo * region.tile <= pos < region.hi * region.tile:
                return {}
        key = (m.ref_id, str(spec), region.serial if region is not None else 0,
               getattr(t, "_group_key", None))
        cache = getattr(t, "_sort_cache", None)
        if cache is not None and cache[0] == key:
            return cache[1]
        if t.grouped:
            hits = region.query_groups(pos, pos + 1) if region is not None else \
                t.fetch_groups(m.ref_id, pos, pos + 1)
            counts = t.group_rows(m.ref_id)
        else:
            hits = [(it, None, r) for it, r in (region.query(pos, pos + 1) if region is not None
                                                else t.fetch(m.ref_id, pos, pos + 1))]
            counts = {None: len(region.last_end) if region is not None
                      else t.row_count(m.ref_id)}
        hits = [h for h in hits if h[0].start <= pos < h[0].end]
        keyfn = style.sort_key_function(
            spec.get("by", "base"), pos, [h[0] for h in hits], reference=m.reference or "",
            tag_name=spec.get("tag", ""), read_groups=getattr(t, "read_groups", {}))
        perm: dict = {}
        by_group: dict = {}
        for it, g, r in hits:
            by_group.setdefault(g, []).append((keyfn(it), r))
        for g, rows in by_group.items():
            rows.sort(key=lambda kr: kr[0])
            first = [r for _k, r in rows]
            taken = set(first)
            rest = [r for r in range(max(counts.get(g, 0), max(taken) + 1)) if r not in taken]
            for new, old in enumerate(first + rest):
                if new != old:
                    perm[(g, old)] = new
        t._sort_cache = (key, perm)
        return perm

    def extra_tracks_h(self) -> int:
        return sum(self._extra_track_h(t) for t in self.tracks)

    def track_at_y(self, y: float) -> int:
        """Index of the extra track at header y (unwrapped coordinates), or -1."""
        for i, (top, h) in enumerate(self._extra_spans):
            if top <= y < top + h:
                return i
        return -1

    def _header_y(self, pos: QPointF) -> float | None:
        if self.wrap_active():
            loc = self._wrap_locate(pos)
            return None if loc is None else loc[1]
        if self.model is not None and self.model.browser:
            return pos.y() + self.verticalScrollBar().value()
        return pos.y()

    def in_header(self, pos: QPointF) -> bool:
        if self.wrap_active() or (self.model is not None and self.model.browser):
            return True
        return pos.y() < self.header_h()

    # ---------------------------------------------------------------- painting

    def _paint_extra_tracks(self, p: QPainter, bp0: int, bp1: int, w: int) -> None:
        m = self.model
        if m is None:
            return
        pal = self.palette()
        mid = pal.color(pal.ColorRole.Mid)
        for t, (top, h) in zip(self.tracks, self._extra_spans):
            p.setPen(QPen(pal.color(pal.ColorRole.Midlight), 1))
            p.drawLine(QPointF(self.label_w, top + h - 2), QPointF(w, top + h - 2))
            try:
                if t.kind == gt.KIND_ANNOTATION:
                    self._paint_body(p, t, top, h - TRACK_GAP, bp0, bp1)
                else:
                    span = self._reads_span(t, bp0, bp1, w)
                    if isinstance(span, tuple):
                        # coverage and reads come from the tiles, coverage scaled to
                        # the peak of the surrounding region so the tiles stay valid
                        self._tile_peak = self._coverage_peak(t, *span)
                        self._cov_max = getattr(self, "_cov_max", {})
                        self._cov_max[id(t)] = self._tile_peak
                        self._paint_body(p, t, top, h - TRACK_GAP, *span,
                                         extra_key=self._tile_peak)
                        if t.grouped:
                            self._paint_group_heads(p, t, top + ALN_COV_H + 2, h, w)
                    else:
                        self._paint_alignment_head(p, t, top, bp0, bp1, w, mid, note=span)
            except Exception as exc:  # noqa: BLE001 - a bad file must not break the view
                p.setPen(mid)
                p.drawText(QRectF(self.label_w + 6, top, w, NOTE_H),
                           Qt.AlignmentFlag.AlignVCenter, f"Cannot draw: {exc}")

    # Drawn rows are kept as pixmap tiles (TILE_PX wide, per zoom level): panning only
    # draws the tiles that come into view and copies the others.
    TILE_PX = 512
    TILE_CACHE_BYTES = 256_000_000

    def _reads_span(self, t, bp0: int, bp1: int, w: int):
        """(a, b) to draw reads for, or the reason (text) they are not drawn."""
        m = self.model
        a = max(0, self.bp_at(max(self.label_w, self.x_of(bp0))))
        b = max(a + 1, self.bp_at(min(w, self.x_of(bp1))))
        if not t.reads_visible(a, b):
            return f"Zoom in below {gt.BAM_READ_WINDOW // 1000:,} kb to see the reads"
        n_here = t.count_in(m.ref_id, a, b + 1)
        if n_here > MAX_DRAWN:
            return f"{n_here:,} alignments here: zoom in to see them"
        return a, b

    def _coverage_peak(self, t, a: int, b: int) -> float:
        """Depth scale for the coverage: the highest depth within a view on either
        side, rounded up to 1, 1.5, 2, 3, 5 or 7 × 10ⁿ so it rarely changes."""
        span = max(1, b - a)
        raw = t.coverage_peak(self.model.ref_id, max(0, a - span), b + span)
        if raw <= 0:
            return 1.0
        mag = 10 ** np.floor(np.log10(raw))
        for f in (1, 1.5, 2, 3, 5, 7, 10):
            if raw <= f * mag:
                return float(f * mag)
        return float(raw)

    def _paint_body(self, p: QPainter, t, top: int, height: int, bp0: int, bp1: int,
                    extra_key=None) -> None:
        m = self.model
        w = self.viewport().width()
        if self._wrap_origin is not None or height <= 0:
            self._cx0, self._cx1 = float(self.label_w), float(w)
            self._draw_body(p, t, top, bp0, bp1)
            return
        from collections import OrderedDict

        from PySide6.QtGui import QPixmap

        cache = getattr(self, "_tile_cache", None)
        if cache is None:
            cache = self._tile_cache = OrderedDict()
        tw = self.TILE_PX
        hval = self.horizontalScrollBar().value()
        dpr = self.viewport().devicePixelRatioF()
        k0 = max(0, int(hval // tw))
        k1 = int((hval + self._data_w()) // tw)
        colors = getattr(self, "type_colors", {}) or {}
        base_key = (id(t), m.ref_id, self.ppb, self._mode(t), self._shown_rows(t), height,
                    repr(self._layout(t)[0]) if t.kind != gt.KIND_ANNOTATION else None,
                    t.layout_version(m.ref_id, bp0, bp1 + 1), dpr,
                    hash(tuple(sorted(colors.items()))), getattr(self, "_tiles_epoch", 0),
                    extra_key)
        for k in range(k0, k1 + 1):
            x = self.label_w + k * tw - hval
            key = (*base_key, k)
            pm = cache.get(key)
            if pm is None:
                pm = QPixmap(int(tw * dpr), int(height * dpr))
                pm.setDevicePixelRatio(dpr)
                pm.fill(Qt.GlobalColor.transparent)
                tp = QPainter(pm)
                tp.setFont(p.font())
                tp.translate(-x, -top)
                self._cx0, self._cx1 = float(x), float(x + tw)
                a = max(0, int(k * tw / self.ppb) - 1)
                b = min(m.length, int((k + 1) * tw / self.ppb) + 2)
                # feature names are written once over the whole view, not per tile
                self._draw_body(tp, t, top, a, b, labels=False)
                tp.end()
                if t.tile_complete(m.ref_id, a, b):
                    cache[key] = pm
                    used = sum(v.width() * v.height() * 4 for v in cache.values())
                    while len(cache) > 8 and used > self.TILE_CACHE_BYTES:
                        _k, old = cache.popitem(last=False)
                        used -= old.width() * old.height() * 4
            else:
                cache.move_to_end(key)
            p.drawPixmap(QPointF(x, top), pm)
        self._cx0, self._cx1 = float(self.label_w), float(w)
        if t.kind == gt.KIND_ANNOTATION:
            pal = self.palette()
            self._paint_annotation_track(p, t, top, bp0, bp1, pal.color(pal.ColorRole.Text),
                                         shapes=False)

    def _draw_body(self, p: QPainter, t, top: int, a: int, b: int, labels: bool = True) -> None:
        pal = self.palette()
        if t.kind == gt.KIND_ANNOTATION:
            self._paint_annotation_track(p, t, top, a, b, pal.color(pal.ColorRole.Text),
                                         labels=labels)
        else:
            self._paint_coverage_bars(p, t, top, a, b, getattr(self, "_tile_peak", 1.0))
            self._paint_alignment_rows(p, t, top + ALN_COV_H + 2, a, b)

    def _paint_coverage_bars(self, p, t, top, a, b, peak) -> None:
        x0, x1 = self.x_of(a), self.x_of(b)
        n = max(1, int((x1 - x0) // 2))
        cov = t.coverage(self.model.ref_id, a, b, n)
        if cov is None or not len(cov):
            return
        from PySide6.QtGui import QPolygonF

        base = top + ALN_COV_H - 1
        step = (x1 - x0) / len(cov)
        ys = (base - np.minimum(cov, peak) * ((ALN_COV_H - 4) / peak)).tolist()
        xs = (x0 + np.arange(len(cov) + 1) * step).tolist()
        pts = [QPointF(x0, base)]
        for i, y in enumerate(ys):
            pts += (QPointF(xs[i], y), QPointF(xs[i + 1], y))
        pts.append(QPointF(xs[-1], base))
        p.setPen(Qt.PenStyle.NoPen)
        p.setBrush(_COV_FILL)
        p.drawPolygon(QPolygonF(pts))
        p.setBrush(Qt.BrushStyle.NoBrush)

    def _paint_annotation_track(self, p, t, top, bp0, bp1, fg, *, shapes: bool = True,
                                labels: bool = True) -> None:
        m = self.model
        rh = self._row_h(t)
        max_rows = self._shown_rows(t)
        squished = self._mode(t) == MODE_SQUISHED
        colors = getattr(self, "type_colors", {}) or {}
        font = QFont(p.font())
        font.setPointSize(7)
        p.setFont(font)
        fm = QFontMetrics(font)
        body_h = rh - 3 if not squished else rh - 1
        named: set = set()        # (row, name): one name per feature (and per gene's parts)
        for it, row in t.fetch(m.ref_id, bp0, bp1 + 1):
            row = min(row, max_rows - 1)
            y = top + row * rh
            x0, x1 = self.x_of(it.start), self.x_of(it.end)
            color = _type_color(it.kind, colors)
            if shapes:
                self._draw_feature(p, it, x0, x1, y, body_h, color, squished)
            if not labels or squished or not it.name or (row, it.name) in named:
                continue
            text_w = fm.horizontalAdvance(it.name)
            vis0 = max(x0, getattr(self, "_cx0", self.label_w))
            vis1 = min(x1, getattr(self, "_cx1", self.viewport().width()))
            if it.blocks:
                label = self._gene_label(it, color, y + 1, body_h, text_w, vis0, vis1)
                if label is None:
                    continue
                rect, pen, over_line, fill = label
                if over_line:
                    # clear the strand arrows under the name, keep the line and exons
                    p.save()
                    p.setClipRect(rect)
                    p.fillRect(rect, self.palette().color(self.palette().ColorRole.Base))
                    self._draw_feature(p, it, x0, x1, y, body_h, color, squished, ticks=False)
                    p.restore()
                    p.setFont(font)
                else:
                    # hide the chevrons behind the name, as on the other tracks
                    cx = rect.center().x()
                    left, right = _clear_of_chevrons(
                        cx - text_w / 2 - 4, cx + text_w / 2 + 4,
                        x0 + (8 if it.strand == "-" else 0) + 7)
                    patch = QRectF(left, rect.top(), right - left, rect.height()) & rect
                    p.fillRect(patch, fill)
                p.setPen(pen)
                p.drawText(rect, Qt.AlignmentFlag.AlignCenter, it.name)
            elif vis1 - vis0 > text_w + 10:
                cx = (vis0 + vis1) / 2
                left, right = _clear_of_chevrons(cx - text_w / 2 - 4, cx + text_w / 2 + 4,
                                                 x0 + 7)
                patch = QRectF(left, y + 1, right - left, body_h)
                p.fillRect(patch, color)
                p.setPen(QColor("#ffffff") if color.lightness() < 150 else QColor("#000000"))
                p.drawText(patch, Qt.AlignmentFlag.AlignCenter, it.name)
            else:
                continue
            named.add((row, it.name))

    def _draw_feature(self, p, it, x0: float, x1: float, y: float, body_h: float,
                      color: QColor, squished: bool, ticks: bool = True) -> None:
        if it.blocks:
            # gene model: line through introns, exons thin, CDS full height
            mid_y = y + 1 + body_h / 2
            p.setPen(QPen(color.darker(130), 1))
            p.drawLine(QPointF(x0, mid_y), QPointF(x1, mid_y))
            thin = max(2.0, body_h * 0.55)
            for s_, e_, _q in it.blocks:
                p.fillRect(QRectF(self.x_of(s_), mid_y - thin / 2,
                                  max(1.0, (e_ - s_) * self.ppb), thin), color.lighter(115))
            for s_, e_ in it.thick:
                p.fillRect(QRectF(self.x_of(s_), y + 1, max(1.0, (e_ - s_) * self.ppb),
                                  body_h), color)
            if ticks and not squished and x1 - x0 > 30:
                exons = [(self.x_of(s_), self.x_of(e_)) for s_, e_, _q in it.blocks]
                self._strand_ticks(p, it.strand, x0 + (8 if it.strand == "-" else 0),
                                   x1 - (8 if it.strand == "+" else 0), mid_y,
                                   color.darker(160), inside=exons,
                                   inside_color=QColor(255, 255, 255, 170))
            if it.strand in ("+", "-") and x1 - x0 > 6:
                # the 3' end is pointed, at the height of the part it ends on (a thin
                # UTR gets a thin point), inside the model's span
                rev = it.strand == "-"
                end = it.start if rev else it.end
                in_cds = any(s_ <= end <= e_ for s_, e_ in it.thick)
                h = body_h if in_cds else thin
                top = y + 1 if in_cds else mid_y - thin / 2
                fill = color if in_cds else color.lighter(115)
                d = min(max(2.5, h * 0.6), (x1 - x0) / 3)
                cut = QRectF(x0, top, d, h) if rev else QRectF(x1 - d, top, d, h)
                p.fillRect(cut, self.palette().color(self.palette().ColorRole.Base))
                p.setPen(Qt.PenStyle.NoPen)
                p.setBrush(fill)
                p.drawPolygon(_pointed_bar(cut.left(), cut.right(), top, h, rev, tip=d))
                p.setBrush(Qt.BrushStyle.NoBrush)
        elif it.strand in ("+", "-") and x1 - x0 > 4:
            p.setPen(Qt.PenStyle.NoPen)
            p.setBrush(color)
            p.drawPolygon(_pointed_bar(x0, x1, y + 1, body_h, it.strand == "-"))
            p.setBrush(Qt.BrushStyle.NoBrush)
            if ticks and not squished and x1 - x0 > 16:
                self._strand_ticks(p, it.strand, x0, x1 - 6, y + 1 + body_h / 2,
                                   QColor(255, 255, 255, 170))
        else:
            p.fillRect(QRectF(x0, y + 1, max(1.0, x1 - x0), body_h), color)
            if ticks and not squished and x1 - x0 > 16:
                self._strand_ticks(p, it.strand, x0, x1, y + 1 + body_h / 2,
                                   QColor(255, 255, 255, 170))

    def _gene_label(self, it, color: QColor, y: float, body_h: float, text_w: int,
                    vis0: float, vis1: float):
        """Where to write a gene model's name: (rect, pen, over the line?) or None.
        Inside its widest visible exon when it fits there, else over the intron line
        without a background, so the exon / intron structure stays visible."""
        best = None
        for parts, fill in ((it.thick, color), ([(s_, e_) for s_, e_, _q in it.blocks],
                                                color.lighter(115))):
            for s_, e_ in parts:
                a, b = max(self.x_of(s_), vis0), min(self.x_of(e_), vis1)
                if b - a > text_w + 6 and (best is None or b - a > best[1] - best[0]):
                    best = (a, b, fill)
            if best is not None:
                break               # CDS (full height) first, the thin UTR parts after
        if best is not None:
            a, b, fill = best
            pen = QColor("#ffffff") if fill.lightness() < 150 else QColor("#000000")
            if fill != color:
                # thin UTR: the name sits on a full-height patch of the feature colour,
                # so it reads the same (white on colour) as everywhere else
                fill = color
                pen = QColor("#ffffff") if fill.lightness() < 150 else QColor("#000000")
            return QRectF(a, y, b - a, body_h), pen, False, fill
        if vis1 - vis0 > text_w + 10:
            base = self.palette().color(self.palette().ColorRole.Base)
            pen = QColor("#000000") if base.lightness() > 128 else QColor("#ffffff")
            cx = (vis0 + vis1) / 2
            # no background: the intron line and thin exons show between the letters
            return QRectF(cx - text_w / 2 - 3, y, text_w + 6, body_h), pen, True, None
        return None

    @staticmethod
    def _strand_ticks(p, strand: str, x0: float, x1: float, y: float, color: QColor,
                      skip: QRectF | None = None, inside=None,
                      inside_color: QColor | None = None) -> None:
        """Chevrons along a feature; with ``inside`` (x ranges, e.g. exons) those
        within a range use ``inside_color`` so they read on the filled part."""
        if strand not in "+-" or strand == "":
            return
        step = 14
        d = 2.5 if strand == "+" else -2.5
        out_ln, in_ln = [], []
        x = x0 + step / 2
        while x < x1 - 3:
            if skip is not None and skip.left() - 3 <= x <= skip.right() + 3:
                x += step
                continue
            on = inside is not None and any(a <= x - abs(d) and x + abs(d) <= b
                                            for a, b in inside)
            if inside is not None and not on and any(a - abs(d) < x < b + abs(d)
                                                     for a, b in inside):
                x += step       # straddles an exon edge: skip it
                continue
            (in_ln if on else out_ln).extend((
                QLineF(x - d, y - 2.5, x, y), QLineF(x - d, y + 2.5, x, y)))
            x += step
        if out_ln:
            p.setPen(QPen(color, 1))
            p.drawLines(out_ln)
        if in_ln:
            p.setPen(QPen(inside_color, 1))
            p.drawLines(in_ln)

    def _paint_alignment_head(self, p, t, top, bp0, bp1, w, mid, note: str = "") -> None:
        """Wide views: coverage only (from the whole-genome cache when built), and why
        no reads are drawn."""
        m = self.model
        x_start = int(max(self.label_w, self.x_of(bp0)))
        x_end = int(min(w, self.x_of(bp1)))
        a = max(0, self.bp_at(x_start))
        b = max(a + 1, self.bp_at(x_end))
        cov = t.coverage(m.ref_id, a, b, max(1, (x_end - x_start) // 2))
        self._cov_max = getattr(self, "_cov_max", {})
        if cov is not None and len(cov):
            peak = float(cov.max()) or 1.0
            self._cov_max[id(t)] = peak
            from PySide6.QtGui import QPolygonF

            base = top + ALN_COV_H - 1
            step = (x_end - x_start) / len(cov)
            ys = (base - cov * ((ALN_COV_H - 4) / peak)).tolist()
            xs = (x_start + np.arange(len(cov) + 1) * step).tolist()
            pts = [QPointF(x_start, base)]
            for i, y in enumerate(ys):
                pts += (QPointF(xs[i], y), QPointF(xs[i + 1], y))
            pts.append(QPointF(xs[-1], base))
            p.setPen(Qt.PenStyle.NoPen)
            p.setBrush(_COV_FILL)
            p.drawPolygon(QPolygonF(pts))
            p.setBrush(Qt.BrushStyle.NoBrush)
        else:
            self._cov_max[id(t)] = 0.0
        if note:
            if cov is None and isinstance(t, gt.IndexedBamTrack):
                note += " · right-click the track name to compute whole-genome coverage"
            p.setPen(mid)
            p.drawText(QRectF(self.label_w + 6, top + ALN_COV_H + 2, w - self.label_w - 12,
                              NOTE_H),
                       Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft, note)

    def _paint_alignment_rows(self, p, t, rows_top, a, b) -> None:
        m = self.model
        rh = self._row_h(t)
        collapsed = self._mode(t) == MODE_COLLAPSED
        squished = self._mode(t) == MODE_SQUISHED or collapsed
        show_mism = self.ppb >= MISMATCH_MIN_PPB and not squished
        letters = self.ppb >= 7.0 and not squished
        ref = m.reference or ""
        font = QFont(p.font())
        font.setPointSize(7)
        # collapsed: rows touch, so the reads read as a density picture
        body = rh if collapsed else rh - (2 if not squished else 1)
        ppb = self.ppb
        x_base = self.x_of(0)
        fills: dict = {}
        rects: dict = {}      # style -> bars, drawn with one call per colour
        polys: dict = {}      # style -> read ends drawn pointed (strand)
        later = []            # reads whose mismatches are drawn on top
        clip_marks: list = []
        color_of = self._read_colorer(t)
        shade = bool(t.settings.get("shade_mapq", True))
        clip_min = int(t.settings.get("clip_min", DEFAULT_CLIP_MIN)) \
            if t.settings.get("flag_clips", True) else 0
        alleles = {} if show_mism else None     # position -> {base: reads}
        pairs = t.pairs and not squished and not collapsed
        chrom = t.resolve(m.ref_id)
        mate_lines: list = []
        intron_lines: list = []
        for it, yoff in self._placed(t, a, b + 1):
            if alleles is not None and ref:
                if it.cs:
                    gt.cs_details(it)
                mm = gt.mismatches_against(it, ref, a, b + 1, offsets=True)
                for pos, base, _k in mm:
                    d = alleles.setdefault(pos, {})
                    d[base] = d.get(base, 0) + 1
            else:
                mm = ()
            if yoff is None:
                continue
            y = rows_top + yoff
            rev = it.strand == "-"
            faint = it.secondary or (shade and it.mapq is not None and it.mapq == 0)
            key = (color_of(it), faint)
            if key not in fills:
                fill = QColor(key[0])
                edge = _edge_color(fill)
                if faint:
                    fill.setAlpha(_SECONDARY_ALPHA)
                    edge.setAlpha(_SECONDARY_ALPHA)
                fills[key] = (fill, QPen(edge, 1))
                rects[key] = []
                polys[key] = []
            parts = it.ref_blocks
            if pairs:
                span = gt._pair_span(it, chrom)
                if span is not None:
                    # the line to the mate, under the bars (each mate draws its side)
                    x0 = x_base + (it.end if span[1] >= it.start else span[1]) * ppb
                    x1 = x_base + (max(span[1], it.start)) * ppb
                    if x1 > x0:
                        mate_lines.append(QLineF(x0, y + body / 2, x1, y + body / 2))
            if len(parts) > 1:
                intron_lines.append(QLineF(x_base + parts[0][1] * ppb, y + body / 2,
                                 x_base + parts[-1][0] * ppb, y + body / 2))
            rc = rects[key]
            pointed = not squished and not collapsed and ppb * (it.end - it.start) > 10
            end_i = (0 if rev else len(parts) - 1) if pointed else -1
            for i, (s_, e_, _q) in enumerate(parts):
                bx0 = x_base + s_ * ppb
                bw = max(1.0, (e_ - s_) * ppb)
                if i == end_i and bw > 3:
                    # pointed end shows the strand, filled and inside the read's span
                    polys[key].append(_pointed_bar(bx0, bx0 + bw, y, body, rev))
                else:
                    rc.append(QRectF(bx0, y, bw, body))
            if clip_min:
                left, right = style.clips(it)
                if left >= clip_min:
                    x = x_base + it.start * ppb
                    clip_marks.append(QLineF(x, y - 0.5, x, y + body + 0.5))
                if right >= clip_min:
                    x = x_base + it.end * ppb
                    clip_marks.append(QLineF(x, y - 0.5, x, y + body + 0.5))
            if show_mism:
                later.append((it, y, mm))
        if mate_lines:
            p.setPen(QPen(_MATE_LINE, 1))
            p.drawLines(mate_lines)
        if intron_lines:
            # whatever the reads' colour, so a splice reads as one grey band
            p.setPen(QPen(_INTRON_LINE, 1))
            p.drawLines(intron_lines)
        p.setPen(Qt.PenStyle.NoPen)
        for key, (fill, _pen) in fills.items():
            if rects[key] or polys[key]:
                p.setBrush(fill)
                if rects[key]:
                    p.drawRects(rects[key])
                for poly in polys[key]:
                    p.drawPolygon(poly)
        p.setBrush(Qt.BrushStyle.NoBrush)
        if later:
            self._paint_read_details(p, t, later, a, b, body, font, letters)
        if clip_marks:
            p.setPen(QPen(_CLIP, 2 if ppb >= 0.5 else 1.5))
            p.drawLines(clip_marks)
        if alleles:
            self._paint_allele_bars(p, t, rows_top - 2 - ALN_COV_H, a, b, alleles)

    def _paint_read_details(self, p, t, later, a, b, body, font, letters) -> None:
        """Mismatches (faded by base quality) and insertions (with their length when
        there is room) over the read bars."""
        ppb = self.ppb
        x_base = self.x_of(0)
        by_quality = bool(t.settings.get("shade_quality", True))
        fm = QFontMetrics(font)
        p.setFont(font)
        boxes: dict = {}          # (base, opacity step) -> rects, one call per colour
        texts = []
        w = max(1.0, ppb)
        for it, y, mm in later:
            quals = style.qualities(it) if by_quality and mm else None
            nq = len(quals) if quals is not None else 0
            for pos, base, k in mm:
                step = 8
                if nq and 0 <= k < nq:
                    step = round(style.quality_alpha(quals[k]) * 8)
                rect = QRectF(x_base + pos * ppb, y, w, body)
                boxes.setdefault((base, step), []).append(rect)
                if letters:
                    texts.append((rect, base))
        p.setPen(Qt.PenStyle.NoPen)
        for (base, step), rs in boxes.items():
            c = QColor(self._colors.get(base, (QColor("#888888"), None))[0])
            if step < 8:
                c.setAlphaF(step / 8)
            p.setBrush(c)
            p.drawRects(rs)
        p.setBrush(Qt.BrushStyle.NoBrush)
        for rect, base in texts:
            p.setPen(self._colors.get(base, (None, QColor("#000000")))[1])
            p.drawText(rect, Qt.AlignmentFlag.AlignCenter, base)
        for it, y, _mm in later:
            for pos, n in style.insertions(it):
                if not a <= pos <= b:
                    continue
                x = x_base + pos * ppb
                label = f"{n:,}"
                tw = fm.horizontalAdvance(label) + 4
                if n > 1 and (ppb >= 2 or n >= 20) and body >= 7:
                    box = QRectF(x - tw / 2, y, tw, body)
                    p.fillRect(box, _INS)
                    p.setPen(QColor("#ffffff"))
                    p.drawText(box, Qt.AlignmentFlag.AlignCenter, label)
                else:
                    p.setPen(QPen(_INS, 2 if ppb > 2 else 1))
                    p.drawLine(QPointF(x, y), QPointF(x, y + body))

    def _paint_allele_bars(self, p, t, cov_top, a, b, alleles) -> None:
        """Coverage bars split by base where enough reads differ from the reference
        (IGV's allele-frequency threshold, 20 % by default)."""
        threshold = float(t.settings.get("allele_threshold", 0.2))
        lo, hi = max(a, 0), max(a + 1, b + 1)
        cov = t.coverage(self.model.ref_id, lo, hi, hi - lo)
        if cov is None or not len(cov):
            return
        peak = getattr(self, "_tile_peak", 1.0) or 1.0
        scale = (ALN_COV_H - 4) / peak
        base_y = cov_top + ALN_COV_H - 1
        ppb = self.ppb
        x_base = self.x_of(0)
        ref = self.model.reference or ""
        bars: dict = {}
        for pos, counts in alleles.items():
            k = pos - lo
            if not 0 <= k < len(cov):
                continue
            depth = float(cov[k])
            alt = sum(counts.values())
            if depth < 3 or alt / depth < threshold:
                continue
            x = x_base + pos * ppb
            w = max(1.0, ppb - (1 if ppb >= 4 else 0))
            y = base_y
            parts = sorted(counts.items(), key=lambda kv: -kv[1])
            ref_n = max(0.0, depth - alt)
            if pos < len(ref) and ref_n:
                parts.append((ref[pos].upper(), ref_n))
            for base, n in parts:
                h = min(depth, peak) * scale * (n / depth)
                bars.setdefault(base, []).append(QRectF(x, y - h, w, h))
                y -= h
        p.setPen(Qt.PenStyle.NoPen)
        for base, rs in bars.items():
            p.setBrush(self._colors.get(base, (QColor("#888888"), None))[0])
            p.drawRects(rs)
        p.setBrush(Qt.BrushStyle.NoBrush)

    def _paint_color_legend(self, p: QPainter, t, y: float, y_end: float,
                            font: QFont, mid: QColor) -> None:
        """'Colour by' mode and its key, in the name column under the track name."""
        mode = t.settings.get("color_by", style.COLOR_STRAND)
        tag_name = t.settings.get("color_tag", "")
        rows = [("", "Colour: " + style.mode_label(mode, tag_name))]
        rows += style.legend(mode, tag_name=tag_name)
        if mode in (style.COLOR_READ_GROUP, style.COLOR_SAMPLE, style.COLOR_LIBRARY):
            palette = style.group_colors(getattr(t, "read_groups", {}), mode)
            rows += [(c, v) for v, c in list(palette.items())[:8]]
        fm = QFontMetrics(font)
        p.setFont(font)
        line_h = 11
        for color, text in rows:
            if y + line_h > y_end:
                break
            x = 6.0
            if color:
                p.setPen(Qt.PenStyle.NoPen)
                p.setBrush(QColor(color))
                p.drawRoundedRect(QRectF(x, y + 2, 9, 7), 2, 2)
                p.setBrush(Qt.BrushStyle.NoBrush)
                x += 13
            p.setPen(mid)
            p.drawText(QRectF(x, y, self.label_w - x - 4, line_h),
                       Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft,
                       fm.elidedText(text, Qt.TextElideMode.ElideRight, int(self.label_w - x - 4)))
            y += line_h

    def _read_colorer(self, t):
        """Function item -> fill colour (hex) for the track's "Colour by" setting."""
        mode = t.settings.get("color_by", style.COLOR_STRAND)
        if mode == style.COLOR_STRAND:
            return lambda it: style.REV if it.strand == "-" else style.FWD
        chrom = self.model.ref_id if self.model is not None else ""
        kw = {"chrom": chrom, "tag_name": t.settings.get("color_tag", ""),
              "read_groups": getattr(t, "read_groups", {}),
              "palette": style.group_colors(getattr(t, "read_groups", {}), mode)}
        if mode == style.COLOR_INSERT:
            kw["insert_range"] = t.insert_range(chrom)
        return lambda it: style.read_color(it, mode, **kw)

    def _paint_extra_labels(self, p: QPainter, fg: QColor, mid: QColor) -> None:
        font = QFont(p.font())
        fm = QFontMetrics(font)
        small = QFont(font)
        small.setPointSize(7)
        for t, (top, h) in zip(self.tracks, self._extra_spans):
            p.setFont(font)
            p.setPen(fg)
            p.drawText(QRectF(6, top, self.label_w - 10, 16),
                       Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft,
                       fm.elidedText(t.name, Qt.TextElideMode.ElideMiddle, self.label_w - 12))
            note = ""
            if t.kind != gt.KIND_ANNOTATION:
                peak = getattr(self, "_cov_max", {}).get(id(t))
                if peak:
                    note = f"depth max {peak:,.0f}"
                if self.model is not None and \
                        not t.grouped and t.row_count(self.model.ref_id) > self._shown_rows(t):
                    note += " · more rows hidden" if note else "more rows hidden"
            elif self.model is not None:
                n = t.count(self.model.ref_id)
                note = f"{n:,} features" if n else "none on this sequence"
            if note and h > 28:
                p.setFont(small)
                p.setPen(mid)
                p.drawText(QRectF(6, top + 14, self.label_w - 10, 12),
                           Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft, note)
            if t.kind != gt.KIND_ANNOTATION:
                self._paint_color_legend(p, t, top + 28, top + h - 4, small, mid)
            p.setPen(QPen(mid, 1))
            p.drawLine(QPointF(0, top + h - 2), QPointF(self.label_w, top + h - 2))
        p.setFont(font)

    # ---------------------------------------------------------------- hover

    def track_items_at(self, pos: QPointF):
        """(track index, items under pos, depth text or None) in a track body, or None."""
        m = self.model
        y = self._header_y(pos)
        if m is None or y is None or not self.in_header(pos) or pos.x() < self.label_w:
            return None
        i = self.track_at_y(y)
        if i < 0:
            return None
        t = self.tracks[i]
        top, _h = self._extra_spans[i]
        bp = self.bp_at_pos(pos)
        slack = max(1, int(3 / max(self.ppb, 1e-9)))
        if t.kind == gt.KIND_ANNOTATION:
            row = int((y - top) // self._row_h(t))
            last = self._shown_rows(t) - 1
            hits = [it for it, r in t.fetch(m.ref_id, bp - slack, bp + slack + 1)
                    if min(r, last) == row]
        else:
            rows_top = top + ALN_COV_H + 2
            if y < rows_top:
                cov = t.coverage(m.ref_id, bp, bp + 1, 1)
                depth = float(cov[0]) if cov is not None and len(cov) else 0.0
                return i, [], f"Position {bp + 1:,}: depth {depth:,.0f}"
            if not t.reads_visible(bp - slack, bp + slack):
                return i, [], None
            dy = y - rows_top
            rh = self._row_h(t)
            hits = [it for it, yoff in self._placed(t, bp - slack, bp + slack + 1)
                    if yoff is not None and yoff <= dy < yoff + rh]
        # the item under the pointer first, then the near ones
        hits.sort(key=lambda it: 0 if it.start <= bp < it.end else 1)
        return i, hits, None

    def track_tooltip(self, pos: QPointF) -> str | None:
        y = self._header_y(pos)
        if self.model is None or y is None or not self.in_header(pos):
            return None
        i = self.track_at_y(y)
        if i < 0:
            return None
        t = self.tracks[i]
        if pos.x() < self.label_w:
            return f"{t.name}\n{t.path}\nRight-click for options"
        found = self.track_items_at(pos)
        if found is None:
            return None
        _i, hits, depth = found
        if depth:
            return f"{t.name}\n{depth}"
        if not hits:
            return None
        lines = []
        for it in hits[:4]:
            head = f"{it.name or it.kind}  {it.start + 1:,}–{it.end:,} ({it.strand})"
            extra = []
            if it.kind and it.kind != "alignment":
                extra.append(it.kind)
            if it.score is not None:
                extra.append(f"identity {it.score:.1%}")
            if it.mapq is not None:
                extra.append(f"MAPQ {it.mapq}")
            if it.secondary:
                extra.append("secondary")
            info = ", ".join(f"{k} {v}" for k, v in list(it.info.items())[:6]
                             if k not in ("ID", "Parent", "Name"))
            lines.append(head + ("\n  " + " · ".join(extra) if extra else "")
                         + ("\n  " + info if info else ""))
        return "\n".join(lines) + "\nClick for details"

    # ---------------------------------------------------------------- mouse: ruler, track order

    def _in_ruler(self, pos: QPointF) -> bool:
        if self.model is None or self.wrap_active() or pos.x() < self.label_w:
            return False
        y = self._header_y(pos)
        return y is not None and 0 <= y < self.RULER_ZONE

    RULER_ZONE = 24

    def tracks_mouse_press(self, event) -> bool:
        """Left button on the ruler starts a region selection; on a track name, a move."""
        if event.button() != Qt.MouseButton.LeftButton or self.model is None:
            return False
        pos = event.position()
        if self._in_ruler(pos):
            self._ruler_drag = (pos.x(), pos.x())
            self.setFocus(Qt.FocusReason.MouseFocusReason)   # Enter / Esc for the selection
            return True
        if self.tracks and pos.x() >= self.label_w and self.in_header(pos):
            y = self._header_y(pos)
            if y is not None and self.track_at_y(y) >= 0:
                # drag the track body to pan (as in IGV); a click without moving
                # opens the item under the pointer
                self._body_press = (pos.x(), self.horizontalScrollBar().value(), False)
                return True
        if self.tracks and pos.x() < self.label_w and self.in_header(pos):
            y = self._header_y(pos)
            k = self.track_at_y(y) if y is not None else -1
            if k >= 0 and len(self.tracks) > 1:
                self._track_drag = k
                self._track_drop = k        # the arrow stays: the drop line shows where
                return True
        return False

    def tracks_mouse_move(self, event) -> bool:
        pos = event.position()
        press = getattr(self, "_body_press", None)
        if press is not None:
            x0, h0, moved = press
            dx = pos.x() - x0
            if moved or abs(dx) > 3:
                if not moved:
                    self.viewport().setCursor(Qt.CursorShape.ClosedHandCursor)
                    self._body_press = (x0, h0, True)
                self.horizontalScrollBar().setValue(int(h0 - dx))
            return True
        if self._ruler_drag is not None:
            x = max(float(self.label_w), min(pos.x(), float(self.viewport().width())))
            self._ruler_drag = (self._ruler_drag[0], x)
            self.viewport().update()
            return True
        if self._track_drag is not None:
            self._drag_pos = QPointF(pos)
            self._update_track_drop()
            if not getattr(self, "_autoscrolling", False):
                self._drag_autoscroll()
            return True
        return False

    def _update_track_drop(self) -> None:
        y = self._header_y(self._drag_pos)
        if y is not None:
            self._track_drop = self._drop_index(y)
            self.viewport().update()

    def _drag_autoscroll(self) -> None:
        """While a track name is dragged near the top or bottom edge, the page scrolls
        (faster the closer to the edge), so a track can pass a tall read track."""
        self._autoscrolling = False
        if self._track_drag is None:
            return
        y = self._drag_pos.y()
        h = self.viewport().height()
        edge = 40
        step = 0
        if y < edge:
            step = -int(4 + (edge - y) * 0.8)
        elif y > h - edge:
            step = int(4 + (y - (h - edge)) * 0.8)
        if not step:
            return
        v = self.verticalScrollBar()
        before = v.value()
        v.setValue(before + step)
        if v.value() != before:
            self._update_track_drop()
            from PySide6.QtCore import QTimer

            self._autoscrolling = True
            QTimer.singleShot(30, self._drag_autoscroll)

    def tracks_hover_cursor(self, pos: QPointF):
        """Cursor hint: a text beam over the ruler (drag to select a region); None
        elsewhere (track names keep the arrow, and are still dragged to reorder)."""
        if self._in_ruler(pos):
            return Qt.CursorShape.IBeamCursor
        return None

    def tracks_mouse_release(self, event) -> bool:
        press = getattr(self, "_body_press", None)
        if press is not None:
            self._body_press = None
            if press[2]:
                self.viewport().unsetCursor()
            else:
                found = self.track_items_at(event.position())
                if found is not None and found[1]:
                    self.trackItemClicked.emit(found[0], found[1][0],
                                               event.globalPosition().toPoint())
            return True
        if self._ruler_drag is not None:
            x0, x1 = sorted(self._ruler_drag)
            self._ruler_drag = None
            if x1 - x0 >= 4:
                # the region stays selected; its Go button zooms there
                m = self.model
                a = max(0, self.bp_at(x0))
                b = min(m.length, max(a + 1, self.bp_at(x1)))
                self.set_region_selection(a, b)
            else:
                self.set_region_selection(None)
            self.viewport().update()
            return True
        if self._track_drag is not None:
            src, dst = self._track_drag, self._track_drop
            self._track_drag = self._track_drop = None
            self.viewport().unsetCursor()
            if dst is not None:
                # dst is an insertion slot 0..n; account for the removed track
                self.move_track(src, dst - 1 if dst > src else dst)
            self.viewport().update()
            return True
        return False

    # ---------------------------------------------------------------- region selection

    def set_region_selection(self, a: int | None, b: int | None = None) -> None:
        """Select [a, b) on the reference (None clears it). It is drawn as a band
        over every track, with a small bar to zoom there."""
        self.region_sel = None if a is None or b is None or b <= a else (int(a), int(b))
        self._place_sel_chip()
        self.viewport().update()
        if hasattr(self, "regionSelectionChanged"):
            self.regionSelectionChanged.emit()

    def go_to_selection(self) -> None:
        sel = self.region_sel
        if sel is not None:
            self.set_region_selection(None)
            self.show_range(*sel)

    def _sel_chip_widget(self):
        if self._sel_chip is None:
            self._sel_chip = _RegionChip(self)
        return self._sel_chip

    def _place_sel_chip(self) -> None:
        sel, m = self.region_sel, self.model
        if sel is None or m is None or self.wrap_active():
            if self._sel_chip is not None:
                self._sel_chip.hide()
            return
        chip = self._sel_chip_widget()
        chip.set_region(m.ref_id, *sel)
        vp = self.viewport()
        x0, x1 = self.x_of(sel[0]), self.x_of(sel[1])
        lo, hi = float(self.label_w) + 4, float(vp.width()) - chip.width() - 4
        x = min(max((x0 + x1) / 2 - chip.width() / 2, lo), max(lo, hi))
        y = self.RULER_ZONE + 4
        if m.browser:
            y -= self.verticalScrollBar().value()
        y = max(2, y)
        from PySide6.QtCore import QPoint

        target = QPoint(int(x), int(y))
        if chip.pos() != target:
            chip.move(target)
        if not chip.isVisible():
            chip.show()
            chip.raise_()

    def tracks_key_press(self, event) -> bool:
        if event.key() == Qt.Key.Key_Escape and self.region_sel is not None:
            self.set_region_selection(None)
            return True
        if event.key() in (Qt.Key.Key_Return, Qt.Key.Key_Enter) and self.region_sel:
            self.go_to_selection()
            return True
        return False

    def _paint_region_selection(self, p: QPainter, h: float) -> None:
        sel = self.region_sel
        if sel is None or self.model is None:
            return
        x0, x1 = self.x_of(sel[0]), self.x_of(sel[1])
        w = self.viewport().width()
        if x1 < self.label_w or x0 > w:
            return
        x0, x1 = max(x0, float(self.label_w)), min(x1, float(w))
        p.setClipping(False)
        p.fillRect(QRectF(x0, 0, max(1.0, x1 - x0), h), _SEL_FILL)
        p.fillRect(QRectF(x0, 0, max(1.0, x1 - x0), 3), _SEL_EDGE)
        p.setPen(QPen(_SEL_EDGE, 1))
        p.drawLine(QPointF(x0, 0), QPointF(x0, h))
        p.drawLine(QPointF(x1, 0), QPointF(x1, h))

    def _drop_index(self, y: float) -> int:
        """Insertion slot (0..n) for header y: before the track whose middle is below y."""
        for i, (top, h) in enumerate(self._extra_spans):
            if y < top + h / 2:
                return i
        return len(self._extra_spans)

    def _paint_drag_overlay(self, p: QPainter) -> None:
        """Selected ruler region, or where a dragged track would land."""
        h = self.viewport().height() + self.verticalScrollBar().value()
        self._paint_region_selection(p, h)
        if self._ruler_drag is not None:
            x0, x1 = sorted(self._ruler_drag)
            p.setClipping(False)
            p.fillRect(QRectF(x0, 0, max(1.0, x1 - x0), h), QColor(31, 95, 191, 40))
            p.setPen(QPen(QColor("#1f5fbf"), 1))
            p.drawLine(QPointF(x0, 0), QPointF(x0, h))
            p.drawLine(QPointF(x1, 0), QPointF(x1, h))
            if self.model is not None and x1 - x0 >= 4:
                a, b = max(0, self.bp_at(x0)), self.bp_at(x1)
                p.drawText(QPointF(x0 + 4, self.RULER_ZONE + 12), f"{a + 1:,}–{b:,} ({b - a:,} bp)")
        if self._track_drag is not None and self._track_drop is not None and self._extra_spans:
            spans = self._extra_spans
            k = self._track_drop
            y = spans[k][0] if k < len(spans) else spans[-1][0] + spans[-1][1]
            top, hh = spans[self._track_drag]
            p.setClipping(False)
            p.fillRect(QRectF(0, top, self.viewport().width(), hh), QColor(31, 95, 191, 25))
            p.fillRect(QRectF(0, y - 1.5, self.viewport().width(), 3), QColor("#1f5fbf"))


class _RegionChip(QWidget):
    """Floating bar over a ruler selection: where it is, Go (zoom there) and ×."""

    def __init__(self, canvas) -> None:
        from PySide6.QtWidgets import QHBoxLayout, QLabel, QToolButton

        super().__init__(canvas.viewport())
        self.canvas = canvas
        self.setObjectName("regionChip")
        self.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)
        self.setStyleSheet(
            "#regionChip { background: #1f5fbf; border-radius: 5px; }"
            "#regionChip QLabel { color: white; padding: 0 2px; }"
            "#regionChip QToolButton { color: white; border: none; padding: 1px 7px;"
            " border-radius: 4px; font-weight: 600; }"
            "#regionChip QToolButton:hover { background: rgba(255,255,255,0.22); }"
            "#regionChip QToolButton#goButton { background: white; color: #1f5fbf; }"
            "#regionChip QToolButton#goButton:hover { background: #e6eefb; }")
        lay = QHBoxLayout(self)
        lay.setContentsMargins(8, 3, 4, 3)
        lay.setSpacing(6)
        self.label = QLabel(self)
        lay.addWidget(self.label)
        self.go = QToolButton(self)
        self.go.setObjectName("goButton")
        self.go.setText("Go")
        self.go.setToolTip("Zoom to the selected region (Enter)")
        self.go.clicked.connect(canvas.go_to_selection)
        lay.addWidget(self.go)
        self.clear = QToolButton(self)
        self.clear.setText("×")
        self.clear.setToolTip("Clear the selection (Esc)")
        self.clear.clicked.connect(lambda: canvas.set_region_selection(None))
        lay.addWidget(self.clear)
        self.setCursor(Qt.CursorShape.ArrowCursor)

    def set_region(self, chrom: str, a: int, b: int) -> None:
        n = b - a
        size = f"{n / 1e6:.2f} Mb" if n >= 1e6 else f"{n / 1e3:.1f} kb" if n >= 1e3 \
            else f"{n} bp"
        text = f"{chrom}:{a + 1:,}–{b:,}  ·  {size}"
        if self.label.text() != text:
            self.label.setText(text)
            self.adjustSize()
