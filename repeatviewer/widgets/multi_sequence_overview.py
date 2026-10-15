"""Multi-sequence FASTA annotation overview (names | tracks, synced rows)."""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path

from PySide6.QtCore import (
    QEvent,
    QPoint,
    QPointF,
    QRect,
    QRectF,
    QSettings,
    QSize,
    Qt,
    QThread,
    QTimer,
    Signal,
    Slot,
)
from PySide6.QtGui import (
    QColor,
    QContextMenuEvent,
    QCursor,
    QFont,
    QFontMetrics,
    QGuiApplication,
    QImage,
    QKeyEvent,
    QPageSize,
    QPainter,
    QPalette,
    QPdfWriter,
    QPen,
    QPixmap,
    QPolygonF,
    QRegion,
    QWheelEvent,
)
from PySide6.QtWidgets import (
    QApplication,
    QCheckBox,
    QDialog,
    QFileDialog,
    QFrame,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMenu,
    QMessageBox,
    QScrollArea,
    QSizePolicy,
    QSlider,
    QSplitter,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends.annotation_colors import similarity_database_color
from repeatviewer.backends.auto_annotation import AutoAnnotationMode
from repeatviewer.backends.auto_annotation_settings import (
    load_auto_annotation_settings,
    save_auto_annotation_settings,
)
from repeatviewer.backends.database_prepare import database_search_readiness
from repeatviewer.backends.job_cancel import JobCancellation
from repeatviewer.backends.msa_coordinates import (
    DEFAULT_MSA_DIVERGENCE_THRESHOLD,
    MSA_COLOR_ALL_BACKGROUNDS,
    column_residue_frequencies,
    drop_features_overlapping_region,
    mask_rows_region,
    msa_residue_fill_background,
    normalize_msa_color_mode,
    normalize_msa_divergence_threshold,
    remap_features_after_column_delete,
    remove_alignment_columns,
    splice_rows_region,
)
from repeatviewer.backends.sequence_features import (
    SequenceFeature,
    feature_display_pack_key,
    is_aligned_fasta,
    load_aligned_fasta_rows,
    load_fasta_rows,
    load_all_features_by_id,
    filter_features_by_annotation_visibility,
    load_features_by_id_from_json,
    ordered_feature_pack_keys,
    remount_parent_ids,
    resolve_features_gff_path,
    resolve_sidecar_path,
    save_annotations_multi,
    write_aligned_fasta_rows,
    write_fasta_rows,
)
from repeatviewer.backends.similarity_annotation import (
    SimilarityDatabaseRun,
    load_use_system_blast_from_settings,
)
from repeatviewer.backends.similarity_annotation_settings import (
    load_similarity_annotation_settings,
    save_similarity_annotation_settings,
    similarity_feature_track_label,
)
from repeatviewer.backends.te_annotation_batch import (
    TeAnnotationOptions,
    run_te_annotation_on_fasta,
)
from repeatviewer.database_settings import load_similarity_annotation_databases
from repeatviewer.widgets.preview_controls import (
    ANNOTATION_TOOLBAR_RENDERED_HEIGHT,
    AnnotationZoomBar,
)
from repeatviewer.widgets.preview_zoom import (
    amino_acid_color_map,
    invalidate_residue_color_maps,
    nucleotide_base_color_map,
)
from repeatviewer.widgets.sequence_annotation_rail import (
    SequenceAnnotationRail,
    SimilarityDatabaseCheck,
)
from repeatviewer.workers.base import Worker

_AUTO_ANNOTATION_MODES: tuple[tuple[str, AutoAnnotationMode], ...] = (
    ("LTR", AutoAnnotationMode.LTR),
    ("TIR", AutoAnnotationMode.TIR),
    ("Tandem", AutoAnnotationMode.TANDEM),
    ("ORF", AutoAnnotationMode.ORF),
    ("TSD", AutoAnnotationMode.TSD),
)

_BAR_HEIGHT = 7
_FIND_HIT_FILL = QColor(245, 158, 11, 150)
_FIND_CURRENT_FILL = QColor(37, 99, 235, 220)
_FIND_HIT_BORDER = QColor(180, 83, 9, 230)
# Equal blank space above and below the sequence bar within its band.
_BAR_V_PAD = 4
_FEAT_HEIGHT = 12
_FEAT_LABEL_MIN_WIDTH = 28.0
_FEAT_LABEL_POINT_SIZE = 7
_FEAT_EDGE_HIT_PX = 6.0
_FEAT_GAP = 2
# Default gap between sequence track and first annotation (= Annotation distance).
_SEQ_FEATURE_GAP = 1
_MIN_ROW_HEIGHT = _BAR_V_PAD + _BAR_HEIGHT + _BAR_V_PAD
_NAME_COLUMN_DEFAULT = 200
_NAME_COLUMN_MAX = 220
# Narrowest splitter pane: enough for the row index ("12.") without the name.
_NAME_COLUMN_MIN = 28
_NAME_PAD = 8
# Point / 1-column feature bands stay visible when ppb is a fraction of a pixel.
_MIN_COLUMN_BAND_PX = 6.0


def _name_column_min_width(fm: QFontMetrics, n_sequences: int) -> int:
    """Width that still shows ``N.`` for the largest row index."""
    digits = max(1, len(str(max(1, int(n_sequences)))))
    sample = f"{'9' * digits}."
    return max(_NAME_COLUMN_MIN, int(fm.horizontalAdvance(sample) + 2 * _NAME_PAD))
_SELECTED_BG_LIGHT = QColor(11, 92, 171, 36)
_SELECTED_BG_DARK = QColor(80, 160, 255, 48)
_DEFAULT_FEAT = QColor("#c0392b")
# Only guards division/painting against a zero scale.  The actual zoom-out
# limit is calculated per data set so the longest sequence always fits.
_ABS_MIN_PPB = 1e-12
_MAX_PPB = 4.0
_MSA_MAX_PPB = 16.0
_MSA_LETTER_MIN_PPB = 7.5
_MSA_ROW_HEIGHT = 18
# Fixed letter size for MSA/multi Overview: scale double-click and zoom-slider
# 3/4 position target this DNA/protein extent; further zoom only widens ppb.
_MSA_LETTER_SIZE_PT = 10
_MSA_ZOOM_LETTER_SLIDER = 75  # 3/4 of 0–100 slider
# Letter-mode fill: full-band square at slider ¾; width grows to 3× by 4/4.
# Scale 1.0 → when Sequence distance is 0, row edges touch vertically.
_MSA_LETTER_BG_SCALE = 1.0
_MSA_LETTER_BG_MAX_WIDTH_SCALE = 3.0  # bg width at slider max / width at ¾
_DEFAULT_MSA_FONT_PT = _MSA_LETTER_SIZE_PT
_DEFAULT_TRACK_HEIGHT = 4  # Others → Track height (1–10); 10 == size-10 band
_DEFAULT_MSA_SEQUENCE_DISTANCE = 1
_DEFAULT_MULTI_SEQUENCE_DISTANCE = 2
_DEFAULT_MSA_COLUMN_DISTANCE = 0
_OVERVIEW_DISPLAY_SETTINGS_GROUP = "annotation_overview_display"
_MSA_DIFF_STANDARD_BY_PATH_GROUP = "msa_diff_standard_by_path"
_FEATURE_UNDO_LIMIT = 32


def _msa_letter_size_metrics(point_size: int = _MSA_LETTER_SIZE_PT) -> tuple[int, float]:
    """Return ``(residue_band_px, target_ppb)`` for a DNA/protein letter size.

    ``target_ppb`` is the horizontal scale where cells are square for that
    letter size (also the zoom-slider 3/4 / scale double-click anchor).
    """
    pt = max(2, min(36, int(point_size)))
    font = QFont()
    font.setFamily("Menlo")
    font.setStyleHint(QFont.StyleHint.Monospace)
    font.setPointSize(pt)
    fm = QFontMetrics(font)
    glyph_h = max(fm.height(), fm.ascent() + fm.descent())
    # Tight band so size-10 cells stay compact; fill uses ``_MSA_LETTER_BG_SCALE``.
    band = max(1, int(glyph_h + 2))
    ppb = float(max(_MSA_LETTER_MIN_PPB, band))
    return band, ppb


def _msa_letter_bg_rect(
    x: float,
    y: float,
    cell_w: float,
    cell_h: float,
    *,
    letter_side: float | None = None,
) -> QRectF:
    """Centered residue fill for letter zoom.

    At size-10 (slider ¾) the fill is a square of ``letter_side``. From ¾→max
    the fill spans the full column pitch so neighbouring columns stay flush
    (no white gutters); height stays at the size-10 side while width follows
    the cell as zoom grows toward ``_MSA_LETTER_BG_MAX_WIDTH_SCALE``×.
    """
    cell_w = max(0.5, float(cell_w))
    cell_h = max(0.5, float(cell_h))
    if letter_side is None:
        side = max(1.0, min(cell_w, cell_h) * _MSA_LETTER_BG_SCALE)
        return QRectF(
            float(x) + (cell_w - side) * 0.5,
            float(y) + (cell_h - side) * 0.5,
            side,
            side,
        )
    side = max(1.0, float(letter_side) * _MSA_LETTER_BG_SCALE)
    if cell_w <= side + 1e-9:
        # Approaching size-10: keep a square inside the cell.
        s = max(1.0, min(cell_w, cell_h) * _MSA_LETTER_BG_SCALE)
        return QRectF(
            float(x) + (cell_w - s) * 0.5,
            float(y) + (cell_h - s) * 0.5,
            s,
            s,
        )
    # Past ¾: fill the full column so max zoom has no inter-column gaps.
    bg_h = side
    return QRectF(
        float(x),
        float(y) + (cell_h - bg_h) * 0.5,
        cell_w,
        bg_h,
    )


def _clone_features_by_id(
    features_by_id: dict[str, list[SequenceFeature]],
) -> dict[str, list[SequenceFeature]]:
    return {
        rid: [SequenceFeature.from_dict(feat.to_dict()) for feat in feats]
        for rid, feats in features_by_id.items()
    }


@dataclass
class _OverviewUndoSnapshot:
    """Features plus optional MSA rows so Restore can undo column/row deletes."""

    features: dict[str, list[SequenceFeature]]
    sequences: list[str] | None = None
    lengths: dict[str, int] | None = None
    record_ids: list[str] | None = None


# MSA canvas / gap cells stay white (not palette Base grey).
_MSA_CANVAS_LIGHT = QColor(255, 255, 255)
_MSA_CANVAS_DARK = QColor(28, 28, 28)
# Tint only the selected region box (not whole rows / outside the band).
_REGION_BAND_LIGHT = QColor(11, 92, 171, 40)
_REGION_BAND_DARK = QColor(80, 160, 255, 55)
_MAX_ZOOM_OVER_FIT = 32.0
_QT_MAX_WIDGET_SIZE = 16_777_215
_BAR_RIGHT_PAD = 8
_ZOOM_STEP = 1.08
# Compact ruler: scale labels + short ticks (no selection position badges).
_SCALE_HEIGHT = 20
# Ignore tiny moves (click / accidental jitter); keep low for responsive feel.
_RUBBER_MIN_PX = 4
# Below this ppb, MSA overview paints one sample per screen pixel (like the
# single-sequence map / AliView fit view) instead of walking every column.
_MSA_PIXEL_SAMPLE_PPB = 1.0


def _msa_canvas_color(widget: QWidget) -> QColor:
    return _MSA_CANVAS_DARK if _is_dark(widget) else _MSA_CANVAS_LIGHT


def _guess_msa_is_protein(seqs: list[str]) -> bool:
    """Heuristic: protein when most alphabetic residues are not ACGTUN."""
    sample_parts: list[str] = []
    for seq in seqs[:30]:
        sample_parts.append(seq[:400])
        if sum(len(p) for p in sample_parts) >= 8_000:
            break
    letters = [c for c in "".join(sample_parts) if c.isalpha()]
    if not letters:
        return False
    nuc = sum(1 for c in letters if c in "ACGTUN")
    return (nuc / len(letters)) < 0.85


def _is_dark(widget: QWidget) -> bool:
    return widget.palette().color(QPalette.ColorRole.Window).lightness() < 128


def _wheel_delta_y(event: QWheelEvent) -> int:
    """Vertical pixels for a wheel/trackpad event (do not steal horizontal swipes)."""
    pixel = event.pixelDelta()
    if not pixel.isNull() and pixel.y() != 0:
        return int(pixel.y())
    angle = event.angleDelta()
    dy = int(angle.y())
    if dy == 0:
        return 0
    # Mouse wheel notches are typically ±120; map to ~one row, not 120px jumps.
    return int(round(dy * (_MIN_ROW_HEIGHT / 120.0)))


def _wheel_delta_x(event: QWheelEvent, *, scroll_step: int = 40) -> int:
    """Horizontal pixels for trackpad swipe (or Shift+wheel when not zooming)."""
    pixel = event.pixelDelta()
    if not pixel.isNull() and pixel.x() != 0:
        return int(pixel.x())
    angle = event.angleDelta()
    if angle.x() != 0:
        return int(round((angle.x() / 120.0) * max(1, scroll_step)))
    return 0


def _bar_color(widget: QWidget) -> QColor:
    # Soft slate for multi-sequence overview bars.
    return QColor("#b4bac2") if _is_dark(widget) else QColor("#d8dce1")


def _msa_plain_bar_color(widget: QWidget) -> QColor:
    """Light grey for MSA overview when scheme leaves a residue uncoloured."""
    return QColor("#b8b8b8") if _is_dark(widget) else QColor("#d0d0d0")


def _msa_gap_diff_color(widget: QWidget) -> QColor:
    """Fill for a gap that differs from the Color vs standard reference."""
    return QColor("#c0392b") if _is_dark(widget) else QColor("#e74c3c")


def _split_tracks_by_source_pref() -> bool:
    value = QSettings().value("annotation/split_tracks_by_source", True)
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.lower() in ("1", "true", "yes", "on")
    return bool(value)


def _pack_feature_tracks(
    features: list[SequenceFeature],
    *,
    split_by_source: bool | None = None,
) -> list[list[SequenceFeature]]:
    """Pack features into visual bands by (source, track) or track only.

    Within each band, overlapping features stack into additional rows.
    """
    if not features:
        return []
    if split_by_source is None:
        split_by_source = _split_tracks_by_source_pref()
    tracks: list[list[SequenceFeature]] = []
    keys = ordered_feature_pack_keys(features, split_by_source=split_by_source)
    for key in keys:
        group = [
            feat
            for feat in sorted(features, key=lambda f: (f.start, f.end, f.name))
            if feature_display_pack_key(feat, split_by_source=split_by_source) == key
        ]
        band_rows: list[list[SequenceFeature]] = []
        for feat in group:
            placed = False
            for row in band_rows:
                if all(feat.start > other.end or feat.end < other.start for other in row):
                    row.append(feat)
                    placed = True
                    break
            if not placed:
                band_rows.append([feat])
        tracks.extend(band_rows)
    return tracks


def _row_height(
    n_tracks: int,
    *,
    feat_height: int = _FEAT_HEIGHT,
    feat_gap: int = _FEAT_GAP,
) -> int:
    seq_band = _BAR_V_PAD + _BAR_HEIGHT + _BAR_V_PAD
    feat_block = 0
    if n_tracks > 0:
        # Same gap above the first annotation as between annotation rows.
        feat_block = (
            feat_gap
            + n_tracks * feat_height
            + (n_tracks - 1) * feat_gap
        )
    return max(_MIN_ROW_HEIGHT, seq_band + feat_block)


def _feature_display_label(
    feat: SequenceFeature,
    *,
    show_database_name: bool,
) -> str:
    base = str(feat.name or feat.feature_type or "").strip()
    return similarity_feature_track_label(
        base,
        str(feat.feature_type or ""),
        show_database_name=show_database_name,
    )


def _feature_with_bounds(
    feat: SequenceFeature, start: int, end: int
) -> SequenceFeature:
    a, b = int(start), int(end)
    if b < a:
        a, b = b, a
    return SequenceFeature(
        start=a,
        end=b,
        name=feat.name,
        feature_type=feat.feature_type,
        color=feat.color,
        strand=feat.strand,
        track=feat.track,
        parent=feat.parent,
        source=feat.source,
        pident=feat.pident,
        subject_id=feat.subject_id,
    )


def _draw_arrow_feature(
    painter: QPainter,
    rect: QRectF,
    color: QColor,
    *,
    forward: bool,
    label: str = "",
    label_point_size: int = _FEAT_LABEL_POINT_SIZE,
) -> None:
    if rect.width() < 1.5 or rect.height() < 2:
        return
    arrow_w = min(rect.height() * 0.85, max(3.0, rect.width() * 0.2))
    body = QRectF(rect)
    if forward:
        tip_x = rect.right()
        base_x = rect.right() - arrow_w
        body.setRight(base_x)
        arrow = QPolygonF(
            [
                QPointF(tip_x, rect.center().y()),
                QPointF(base_x, rect.top()),
                QPointF(base_x, rect.bottom()),
            ]
        )
    else:
        tip_x = rect.left()
        base_x = rect.left() + arrow_w
        body.setLeft(base_x)
        arrow = QPolygonF(
            [
                QPointF(tip_x, rect.center().y()),
                QPointF(base_x, rect.top()),
                QPointF(base_x, rect.bottom()),
            ]
        )
    painter.setPen(Qt.PenStyle.NoPen)
    painter.setBrush(color)
    if body.width() > 0.5:
        painter.drawRect(body)
    painter.drawPolygon(arrow)
    pt = max(0, min(24, int(label_point_size)))
    if label and pt > 0 and rect.width() > _FEAT_LABEL_MIN_WIDTH:
        painter.setPen(QColor("#ffffff"))
        label_font = QFont(painter.font())
        label_font.setPointSize(pt)
        painter.setFont(label_font)
        text_rect = rect.adjusted(3, 0, -3, 0)
        if not forward:
            text_rect = text_rect.adjusted(arrow_w * 0.4, 0, 0, 0)
        else:
            text_rect = text_rect.adjusted(0, 0, -arrow_w * 0.4, 0)
        painter.drawText(text_rect, Qt.AlignmentFlag.AlignVCenter, label)


def _tick_step_for_span(span_bp: float, bars_width_px: float, min_label_px: float) -> int:
    if bars_width_px <= 1 or span_bp <= 0:
        return max(1, int(span_bp) or 1)
    raw = max(1.0, float(span_bp) * float(min_label_px) / float(bars_width_px))
    exp = math.floor(math.log10(raw)) if raw > 0 else 0
    base = 10**exp
    for mult in (1, 2, 5, 10):
        step = int(mult * base)
        if step >= raw:
            return max(1, step)
    return max(1, int(10 * base))


def _elide_middle(fm: QFontMetrics, text: str, width: int) -> str:
    if width <= 0 or fm.horizontalAdvance(text) <= width:
        return text
    if width < fm.horizontalAdvance("…"):
        return ""
    ell = "…"
    left = max(1, len(text) // 2)
    right = left
    while left + right > 0:
        candidate = (
            text[:left] + ell + text[len(text) - right :] if right else text[:left] + ell
        )
        if fm.horizontalAdvance(candidate) <= width:
            return candidate
        if left >= right and left > 1:
            left -= 1
        elif right > 0:
            right -= 1
        else:
            left -= 1
    return ell


class _NamesColumn(QWidget):
    """Sequence names; row geometry matches the tracks canvas exactly."""

    rowClicked = Signal(int, object)  # row, keyboard modifiers
    rowsDragPreview = Signal(object)  # list[int] while dragging (live)
    rowsPixelPreview = Signal(float, float)  # content y0,y1 while dragging
    rowsDragSelected = Signal(object)  # list[int] on drag release (commit)
    selectAllRequested = Signal()
    deleteRequested = Signal()  # Delete / Backspace (macOS name focus)
    contextMenuRequested = Signal(QPoint, int)  # global pos, row under cursor

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._record_ids: list[str] = []
        self._row_tops: list[int] = []
        self._row_heights: list[int] = []
        self._selected_rows: set[int] = set()
        self._label_point_size = _DEFAULT_MSA_FONT_PT
        self._drag_origin: QPointF | None = None
        self._dragging = False
        self._press_row = -1
        self._press_modifiers: object = Qt.KeyboardModifier.NoModifier
        # Throttle name highlights during drag (~20 Hz); chrome follows every move.
        self._pending_drag_rows: set[int] | None = None
        self._drag_highlight_timer = QTimer(self)
        self._drag_highlight_timer.setSingleShot(True)
        self._drag_highlight_timer.setInterval(50)
        self._drag_highlight_timer.timeout.connect(self._flush_drag_highlight)
        # Elided name cache: invalidate on ids / width / font size.
        self._elided_labels: list[str] = []
        self._elided_labels_key: tuple | None = None
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.setAttribute(Qt.WidgetAttribute.WA_OpaquePaintEvent, True)
        self.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Minimum)

    def set_label_point_size(self, point_size: int) -> None:
        """Name labels use DNA/protein size + 1 (clamped).

        Cap matches zoomed-in MSA letter fonts (see ``_msa_letter_font_for_cell``).
        """
        self._label_point_size = max(2, min(96, int(point_size) + 1))
        self._elided_labels_key = None
        self.update()

    def label_point_size(self) -> int:
        return int(self._label_point_size)

    def set_rows(
        self,
        record_ids: list[str],
        tops: list[int],
        heights: list[int],
        selected: set[int],
        *,
        width: int,
        content_height: int,
    ) -> None:
        self._record_ids = list(record_ids)
        self._row_tops = list(tops)
        self._row_heights = list(heights)
        self._selected_rows = set(selected)
        self._elided_labels_key = None
        min_w = _name_column_min_width(QFontMetrics(self._label_font()), len(record_ids))
        w = max(min_w, int(width))
        h = max(_MIN_ROW_HEIGHT, int(content_height))
        self.setMinimumWidth(min_w)
        self.setMinimumHeight(h)
        self.resize(w, h)
        self.update()

    def set_selected_rows(self, rows: set[int]) -> None:
        new = set(rows)
        if new == self._selected_rows:
            return
        self._selected_rows = new
        self.update()

    def _label_font(self) -> QFont:
        """Name font = DNA/protein label size + 1."""
        font = QFont(self.font())
        font.setPointSize(max(2, int(self._label_point_size)))
        return font

    def _row_at_y(self, y: float) -> int:
        for index, (top, height) in enumerate(zip(self._row_tops, self._row_heights)):
            if top <= y < top + height:
                return index
        return -1

    def _rows_in_y_range(self, y0: float, y1: float) -> list[int]:
        tops = self._row_tops
        heights = self._row_heights
        n = len(tops)
        if n <= 0:
            return []
        lo_y, hi_y = (y0, y1) if y0 <= y1 else (y1, y0)
        lo, hi = 0, n - 1
        first = n
        while lo <= hi:
            mid = (lo + hi) // 2
            if tops[mid] + heights[mid] > lo_y:
                first = mid
                hi = mid - 1
            else:
                lo = mid + 1
        lo, hi = 0, n - 1
        last = -1
        while lo <= hi:
            mid = (lo + hi) // 2
            if tops[mid] < hi_y:
                last = mid
                lo = mid + 1
            else:
                hi = mid - 1
        if first > last:
            return []
        return list(range(first, last + 1))

    def _ensure_elided_labels(self, fm: QFontMetrics, avail: int) -> list[str]:
        key = (len(self._record_ids), int(avail), int(self._label_point_size), self.width())
        if self._elided_labels_key == key and len(self._elided_labels) == len(
            self._record_ids
        ):
            return self._elided_labels
        labels: list[str] = []
        for index, rid in enumerate(self._record_ids):
            prefix = f"{index + 1}."
            numbered = f"{prefix} {rid}"
            if fm.horizontalAdvance(numbered) <= avail:
                label = numbered
            elif fm.horizontalAdvance(prefix) <= avail:
                name_budget = max(0, avail - fm.horizontalAdvance(prefix))
                tail = fm.elidedText(
                    f" {rid}", Qt.TextElideMode.ElideRight, name_budget
                )
                label = prefix + (tail if tail and tail != "…" else "")
            else:
                label = fm.elidedText(prefix, Qt.TextElideMode.ElideRight, avail)
            labels.append(label)
        self._elided_labels = labels
        self._elided_labels_key = key
        return labels

    def mousePressEvent(self, event) -> None:  # noqa: N802
        if event.button() == Qt.MouseButton.LeftButton:
            self._drag_origin = QPointF(event.position())
            self._dragging = False
            self._press_row = self._row_at_y(event.position().y())
            self._press_modifiers = event.modifiers()
            # Do not commit on press — trackpad/three-finger drags would flash
            # the press row before the range selection. Click commits on release.
            self.grabMouse()
            event.accept()
            return
        super().mousePressEvent(event)

    def contextMenuEvent(self, event: QContextMenuEvent) -> None:  # noqa: N802
        row = self._row_at_y(event.pos().y())
        self.contextMenuRequested.emit(event.globalPos(), int(row))
        event.accept()

    def _flush_drag_highlight(self) -> None:
        pending = self._pending_drag_rows
        if pending is None:
            return
        if pending != self._selected_rows:
            self.set_selected_rows(pending)

    def mouseMoveEvent(self, event) -> None:  # noqa: N802
        if self._drag_origin is not None and event.buttons() & Qt.MouseButton.LeftButton:
            # Coalesce backlog: follow the latest cursor, not a stale event.
            current = self.mapFromGlobal(QCursor.pos())
            delta = QPointF(current) - self._drag_origin
            if (
                abs(delta.x()) >= _RUBBER_MIN_PX
                or abs(delta.y()) >= _RUBBER_MIN_PX
            ):
                self._dragging = True
                y0 = float(self._drag_origin.y())
                y1 = float(current.y())
                row_set = set(self._rows_in_y_range(y0, y1))
                # MSA chrome every move; name wash throttled like regional drag.
                self._pending_drag_rows = row_set
                if not self._drag_highlight_timer.isActive():
                    self._flush_drag_highlight()
                    self._drag_highlight_timer.start()
                self.rowsPixelPreview.emit(y0, y1)
            event.accept()
            return
        super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event) -> None:  # noqa: N802
        if event.button() == Qt.MouseButton.LeftButton:
            if QWidget.mouseGrabber() is self:
                self.releaseMouse()
            was_dragging = self._dragging
            origin = self._drag_origin
            press_row = self._press_row
            press_mods = self._press_modifiers
            self._drag_origin = None
            self._dragging = False
            self._press_row = -1
            self._drag_highlight_timer.stop()
            self._pending_drag_rows = None
            if was_dragging and origin is not None:
                current = self.mapFromGlobal(QCursor.pos())
                rubber = QRect(origin.toPoint(), current).normalized()
                rows = self._rows_in_y_range(rubber.top(), rubber.bottom())
                if rows:
                    self.set_selected_rows(set(rows))
                    self.rowsDragSelected.emit(rows)
            elif press_row >= 0:
                self.rowClicked.emit(press_row, press_mods)
            event.accept()
            return
        super().mouseReleaseEvent(event)

    def keyPressEvent(self, event: QKeyEvent) -> None:  # noqa: N802
        from PySide6.QtGui import QKeySequence

        if event.matches(QKeySequence.StandardKey.SelectAll) or (
            event.key() == Qt.Key.Key_A
            and event.modifiers()
            & (Qt.KeyboardModifier.ControlModifier | Qt.KeyboardModifier.MetaModifier)
        ):
            self.selectAllRequested.emit()
            event.accept()
            return
        # Plain Delete/Backspace (and macOS ⌘⌫ = StandardKey.Delete).
        if event.matches(QKeySequence.StandardKey.Delete) or event.key() in (
            Qt.Key.Key_Delete,
            Qt.Key.Key_Backspace,
        ):
            self.deleteRequested.emit()
            event.accept()
            return
        if event.matches(QKeySequence.StandardKey.Undo) or (
            event.key() == Qt.Key.Key_Z
            and event.modifiers()
            & (Qt.KeyboardModifier.ControlModifier | Qt.KeyboardModifier.MetaModifier)
            and not (event.modifiers() & Qt.KeyboardModifier.ShiftModifier)
        ):
            node = self.parent()
            while node is not None:
                restore = getattr(node, "_restore_features", None)
                if callable(restore):
                    restore(True)
                    event.accept()
                    return
                node = node.parent()
        # Forward Shift+arrows to the tracks canvas (feature edge expand).
        if event.modifiers() & Qt.KeyboardModifier.ShiftModifier and event.key() in (
            Qt.Key.Key_Left,
            Qt.Key.Key_Right,
        ):
            node = self.parent()
            while node is not None:
                bars = getattr(node, "_bars", None)
                apply = getattr(bars, "apply_shift_edge_expand", None) if bars else None
                if callable(apply):
                    direction = (
                        "left" if event.key() == Qt.Key.Key_Left else "right"
                    )
                    if apply(direction):
                        event.accept()
                        return
                    break
                node = node.parent()
        super().keyPressEvent(event)

    def paintEvent(self, event) -> None:  # noqa: N802
        painter = QPainter(self)
        clip = event.rect()
        base = self.palette().color(QPalette.ColorRole.Base)
        painter.fillRect(clip, base)
        sel_bg = _SELECTED_BG_DARK if _is_dark(self) else _SELECTED_BG_LIGHT
        text_color = self.palette().color(QPalette.ColorRole.Text)
        mid = self.palette().color(QPalette.ColorRole.Mid)

        font = self._label_font()
        painter.setFont(font)
        fm = QFontMetrics(font)
        avail = max(0, self.width() - 2 * _NAME_PAD)
        labels = self._ensure_elided_labels(fm, avail)
        text_h = fm.height()

        for index in self._rows_in_y_range(float(clip.top()), float(clip.bottom()) + 1.0):
            top = self._row_tops[index]
            height = self._row_heights[index]
            if index in self._selected_rows:
                painter.fillRect(0, top, self.width(), height, sel_bg)
            label = labels[index] if index < len(labels) else ""
            # Center on the row band (MSA residue / bar strip), not the gap below.
            text_y = int(top + max(0, (height - text_h) / 2))
            painter.setPen(text_color)
            painter.drawText(
                QRect(_NAME_PAD, text_y, avail, text_h),
                int(Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft),
                label,
            )

        painter.setPen(QPen(mid, 1))
        painter.drawLine(self.width() - 1, clip.top(), self.width() - 1, clip.bottom())
        painter.end()


class _SelectionChrome:
    """Selection outline/wash painted on the MSA canvas (not a child overlay).

    A translucent child overlay paints an opaque black slab on macOS and hides
    the MSA. Chrome is drawn in the host paintEvent after residues.

    Column live drag uses edge dirty only. Committed columns paint a wash band.
    Row boxes can overlay under a column band (row select, then column select).
    """

    def __init__(self, host: QWidget) -> None:
        self._host = host
        self._mode = ""  # "region" | "columns" | "rows"
        self._box = QRect()
        self._row_boxes: list[QRect] = []
        self._live = False  # True while dragging — outline only, no wash

    def _band_color(self) -> QColor:
        return _REGION_BAND_DARK if _is_dark(self._host) else _REGION_BAND_LIGHT

    def _edge_color(self) -> QColor:
        return QColor("#4ea1ff") if _is_dark(self._host) else QColor("#0b5cab")

    def _visible_clip(self) -> QRect:
        """Visible content rect — prefer the scroll viewport when available.

        ``visibleRegion()`` on a tall content widget can briefly cover the full
        canvas during updates, which made column-band dirties span every row.
        """
        host = self._host
        clamp = getattr(host, "scroll_viewport_rect", None)
        if callable(clamp):
            vis = clamp()
            if not vis.isEmpty():
                return vis
        vis = host.visibleRegion().boundingRect()
        if not vis.isEmpty():
            return vis.intersected(host.rect())
        return host.rect()

    def _request_update(self, dirty: QRegion) -> None:
        if dirty.isEmpty():
            return
        self._host.update(dirty)

    def _border_dirty_region(self, *rects: QRect) -> QRegion:
        """Dirty only outline strips — never the selection bbox interior."""
        region = QRegion()
        pad = 2
        for rect in rects:
            if rect is None or rect.isNull() or rect.width() <= 0 or rect.height() <= 0:
                continue
            r = rect.normalized()
            strips = (
                QRect(r.left() - pad, r.top() - pad, r.width() + 2 * pad, 2 * pad + 1),
                QRect(
                    r.left() - pad,
                    r.bottom() - pad,
                    r.width() + 2 * pad,
                    2 * pad + 1,
                ),
                QRect(r.left() - pad, r.top() - pad, 2 * pad + 1, r.height() + 2 * pad),
                QRect(
                    r.right() - pad,
                    r.top() - pad,
                    2 * pad + 1,
                    r.height() + 2 * pad,
                ),
            )
            for strip in strips:
                region = region.united(QRegion(strip))
        clip = self._visible_clip()
        if clip.isNull() or region.isEmpty():
            return QRegion()
        return region.intersected(QRegion(clip))

    def _area_dirty_region(self, *rects: QRect) -> QRegion:
        dirty = QRect()
        for rect in rects:
            if rect is None or rect.isNull():
                continue
            dirty = rect if dirty.isNull() else dirty.united(rect)
        if dirty.isNull():
            return QRegion()
        # Pad past the dashed pen so right/bottom edges never leave a hairline.
        area = dirty.adjusted(-4, -4, 4, 4).intersected(self._visible_clip())
        return QRegion() if area.isNull() else QRegion(area)

    def _column_edge_dirty_region(self, *rects: QRect) -> QRegion:
        """Dirty only left/right guide strips (never the span between them).

        Filling old→new travel made wide column drags recompute the whole band —
        the lag vs regional. Guides are 1px lines; erasing ``old±pad`` and
        painting ``new±pad`` is enough (same cost as a small regional box).
        """
        region = QRegion()
        pad = 4
        for rect in rects:
            if rect is None or rect.isNull() or rect.width() <= 0:
                continue
            r = rect.normalized()
            top = r.top() - pad
            height = r.height() + 2 * pad
            region = region.united(
                QRegion(QRect(r.left() - pad, top, 2 * pad + 1, height))
            )
            region = region.united(
                QRegion(QRect(r.right() - pad, top, 2 * pad + 1, height))
            )
        clip = self._visible_clip()
        if clip.isNull() or region.isEmpty():
            return QRegion()
        return region.intersected(QRegion(clip))

    def _erase_dirty_for_current(self) -> QRegion:
        """Cheapest dirty region that still erases the current chrome correctly."""
        if not self._mode and not self._row_boxes:
            return QRegion()
        dirty = QRegion()
        if self._mode == "columns":
            # Live guides = edges; committed column wash = band area.
            if self._live:
                dirty = self._column_edge_dirty_region(self._box)
            else:
                dirty = self._area_dirty_region(self._box)
        elif self._mode == "rows":
            if self._live:
                dirty = self._border_dirty_region(*self._row_boxes)
            else:
                dirty = self._area_dirty_region(*self._row_boxes)
        elif self._mode == "region":
            dirty = self._area_dirty_region(self._box)
        # Row wash may sit under a column band.
        if self._mode == "columns" and self._row_boxes:
            dirty = dirty.united(self._area_dirty_region(*self._row_boxes))
        return dirty

    def clear(self) -> None:
        if not self._mode and not self._row_boxes:
            return
        dirty = self._erase_dirty_for_current()
        self._mode = ""
        self._box = QRect()
        self._row_boxes = []
        self._live = False
        self._request_update(dirty)

    def set_region_box(self, box: QRect | None, *, live: bool = False) -> None:
        if box is None or box.isNull() or box.width() <= 0 or box.height() <= 0:
            self.clear()
            return
        new = QRect(box).normalized()
        if self._mode == "region" and self._box == new and self._live == live:
            return
        # Live regional moves must dirty old∪new. Border-only strips are filled
        # white then auto-clipped; at fit zoom MSA cells do not cover every
        # pixel in those strips, which left a white line on the box's right edge.
        dirty = self._erase_dirty_for_current()
        dirty = dirty.united(self._area_dirty_region(new))
        self._mode = "region"
        self._box = new
        self._row_boxes = []
        self._live = bool(live)
        self._request_update(dirty)

    def set_column_band(self, band: QRect | None, *, live: bool = False) -> None:
        if band is None or band.isNull() or band.width() <= 0:
            # Drop column chrome only; keep any row wash.
            if self._mode != "columns":
                return
            dirty = self._erase_dirty_for_current()
            self._mode = "rows" if self._row_boxes else ""
            self._box = QRect()
            self._live = False
            if self._mode == "rows":
                dirty = dirty.united(self._area_dirty_region(*self._row_boxes))
            self._request_update(dirty)
            return
        new = QRect(band).normalized()
        if self._mode == "columns" and self._box == new and self._live == live:
            return
        # Keep existing row boxes (row select then column select).
        kept_rows = list(self._row_boxes)
        dirty = self._erase_dirty_for_current()
        if live:
            dirty = dirty.united(self._column_edge_dirty_region(new))
        else:
            # Committed: wash the MSA band (one-shot area dirty is fine).
            dirty = dirty.united(self._area_dirty_region(new))
        if kept_rows:
            dirty = dirty.united(self._area_dirty_region(*kept_rows))
        self._mode = "columns"
        self._box = new
        self._row_boxes = kept_rows
        self._live = bool(live)
        self._request_update(dirty)

    def set_row_boxes(self, boxes: list[QRect], *, live: bool = False) -> None:
        valid = [
            QRect(b).normalized()
            for b in boxes
            if b is not None and not b.isNull() and b.width() > 0 and b.height() > 0
        ]
        if not valid:
            if self._mode == "columns":
                # Clear row overlay only; keep column band.
                if not self._row_boxes:
                    return
                dirty = self._area_dirty_region(*self._row_boxes)
                self._row_boxes = []
                self._request_update(dirty)
                return
            self.clear()
            return
        # Under an active column band: update row wash without dropping columns.
        if self._mode == "columns" and not self._box.isNull():
            if self._row_boxes == valid and not live:
                return
            dirty = self._area_dirty_region(*self._row_boxes) if self._row_boxes else QRegion()
            dirty = dirty.united(self._area_dirty_region(*valid))
            self._row_boxes = valid
            self._request_update(dirty)
            return
        if self._mode == "rows" and self._row_boxes == valid and self._live == live:
            return
        dirty = self._erase_dirty_for_current()
        if live:
            dirty = dirty.united(self._border_dirty_region(*valid))
        else:
            dirty = dirty.united(self._area_dirty_region(*valid))
        self._mode = "rows"
        self._box = QRect()
        self._row_boxes = valid
        self._live = bool(live)
        self._request_update(dirty)

    def paint(self, painter: QPainter, clip: QRect) -> None:
        if not self._mode and not self._row_boxes:
            return
        edge = self._edge_color()
        pen = QPen(edge, 1, Qt.PenStyle.DashLine)
        band = self._band_color()
        painter.setBrush(Qt.BrushStyle.NoBrush)

        def _paint_row_boxes(*, wash: bool) -> None:
            for box in self._row_boxes:
                if not box.intersects(clip):
                    continue
                if wash:
                    painter.fillRect(box, band)
                painter.setPen(pen)
                painter.drawRect(
                    box.left(),
                    box.top(),
                    max(1, box.width() - 1),
                    max(1, box.height() - 1),
                )

        if self._mode == "rows":
            _paint_row_boxes(wash=not self._live)
            return
        if self._mode == "columns":
            # Row select wash stays visible under column select.
            _paint_row_boxes(wash=True)
            box = self._box
            if box.intersects(clip):
                if not self._live:
                    painter.fillRect(box, band)
                painter.setPen(pen)
                painter.drawLine(box.left(), box.top(), box.left(), box.bottom())
                painter.drawLine(box.right(), box.top(), box.right(), box.bottom())
            return
        if self._mode == "region":
            box = self._box
            if box.intersects(clip):
                if not self._live:
                    painter.fillRect(box, band)
                painter.setPen(pen)
                painter.drawRect(
                    box.left(),
                    box.top(),
                    max(1, box.width() - 1),
                    max(1, box.height() - 1),
                )


class _BarsCanvas(QWidget):
    """Length bars + feature arrows, or aligned MSA residue rows (no name column)."""

    selectionChanged = Signal()
    featureActivated = Signal(int, int)
    featureResized = Signal(str, int, object)  # rid, store_index, SequenceFeature
    zoomChanged = Signal(float)
    regionZoomRequested = Signal(float, float)  # kept for Fit/Select zoom helpers
    regionSelectionChanged = Signal(object)  # tuple[int,int] | None (1-based inclusive)
    contextMenuRequested = Signal(QPoint)  # global pos
    # Prior chrome dismissed at drag start — clear scale band without resetting Position.
    selectionVisualsDismissed = Signal()
    # Live regional drag: (start, end, rows) while dragging; None when cleared.
    regionDragPreview = Signal(object)
    rowsChanged = Signal()  # layout / selection updated — refresh names column
    deleteRequested = Signal()  # Delete / Backspace while tracks focused

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._record_ids: list[str] = []
        self._lengths: dict[str, int] = {}
        self._features: dict[str, list[SequenceFeature]] = {}
        self._tracks: dict[str, list[list[SequenceFeature]]] = {}
        self._msa_seqs: list[str] = []
        self._msa_mode = False
        self._msa_protein = False
        self._msa_colors: dict[str, tuple[str, str]] = {}
        self._msa_qcolors: dict[str, tuple[QColor, QColor]] = {}
        self._msa_residue_height = _MSA_ROW_HEIGHT
        self._msa_font_pt = _MSA_LETTER_SIZE_PT
        self._track_height_level = _DEFAULT_TRACK_HEIGHT
        self._msa_color_mode = MSA_COLOR_ALL_BACKGROUNDS
        self._msa_divergence_threshold = DEFAULT_MSA_DIVERGENCE_THRESHOLD
        self._msa_col_freqs: list[dict[str, float]] = []
        # Explicit Color vs standard pin by record id (None → first sequence).
        self._msa_diff_standard_id: str | None = None
        # Apply default track height → overview band pixels.
        self._apply_track_height_to_residue_band()
        self._sequence_row_gap = 0  # pixels between sequences (Others → Sequence distance)
        self._feat_height = _FEAT_HEIGHT
        self._feat_gap = _FEAT_GAP
        self._feat_label_point_size = _FEAT_LABEL_POINT_SIZE
        self._column_distance = 0  # pixels between residue backgrounds
        self._color_unaligned_background = True
        self._show_annotation_track = True
        # None = no type filter; frozenset = only these annotation_visibility_key values.
        self._visible_annotation_kinds: frozenset[str] | None = None
        # Stable load order for “Original order” sorting.
        self._original_index: list[int] = []
        self._region_sel: tuple[int, int] | None = None
        # When set, column band is limited to these rows (MSA regional drag).
        # ``None`` means full-height column selection (scale / feature).
        self._regional_rows: set[int] | None = None
        # Rows highlighted in the name column while a regional drag is active.
        self._drag_preview_rows: set[int] | None = None
        self._drag_preview_columns: tuple[int, int] | None = None
        self._drag_preview_pending: tuple[int, int, list[int]] | None = None
        # AliView-style temp selection in matrix coords: (c0, c1, r0, r1).
        # Cols are 1-based inclusive; rows are 0-based inclusive.
        self._temp_matrix_sel: tuple[int, int, int, int] | None = None
        self._drag_preview_timer = QTimer(self)
        self._drag_preview_timer.setSingleShot(True)
        self._drag_preview_timer.setInterval(50)  # Position + name highlight (~20 Hz)
        self._drag_preview_timer.timeout.connect(self._flush_region_drag_preview)
        self._row_tops: list[int] = []
        self._row_heights: list[int] = []
        self._max_length = 1
        self._selected_rows: set[int] = set()
        # Selected annotation features as (row, flat_feature_index).
        self._selected_feats: set[tuple[int, int]] = set()
        # Hovered feature for overflow label tip: (row, flat_feature_index).
        self._hover_feat: tuple[int, int] | None = None
        self._show_similarity_database_name = True
        self._ppb = 0.05
        self._bar_pad = 8
        self._min_viewport = QSize(120, 40)
        self._drag_origin: QPointF | None = None
        self._drag_kind: str | None = None  # "region" | "features" | "resize"
        self._press_feat: tuple[int, int] | None = None
        self._press_edge: str | None = None  # "start" | "end"
        self._resize_active = False
        self._resize_rid = ""
        self._resize_store_index = -1
        self._resize_flat = -1
        self._resize_row = -1
        self._resize_edge: str | None = None
        self._resize_orig_start = 0
        self._resize_orig_end = 0
        self._rubber = QRect()
        self._rubber_active = False  # True only after drag exceeds threshold
        # Live selection previews (scale columns / name rows).
        self._preview_columns: tuple[int, int] | None = None
        self._preview_rows: set[int] | None = None
        self._live_row_snap: set[int] | None = None
        self._pending_scroll_bp: float | None = None
        self._pending_anchor: tuple[float, float] | None = None
        # Shift+Left / Shift+Right cycle: 0 = outer side, 1 = include feature edge.
        self._edge_expand_phase: dict[str, int] = {"left": 0, "right": 0}
        # Anchor span for Shift+arrow when expanding a column selection.
        self._edge_expand_anchor: tuple[int, int] | None = None
        # Edit > Find: hit ranges per record id, and the hit being shown.
        self._find_hits: dict[str, list[tuple[int, int]]] = {}
        self._find_current: tuple[str, int, int] | None = None
        self.setMouseTracking(True)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.setSizePolicy(QSizePolicy.Policy.MinimumExpanding, QSizePolicy.Policy.Minimum)
        self.setAttribute(Qt.WidgetAttribute.WA_OpaquePaintEvent, True)
        # Selection chrome is painted in paintEvent after residues (no overlay).
        self._chrome = _SelectionChrome(self)
        # Sticky overscrolled residue+feature snapshot:
        # - While the viewport sits inside the cover with edge slack → pure blit
        # - Near an edge → shift + paint only the newly exposed strip
        # Selection never invalidates this layer — chrome blits on top.
        self._msa_layer_cache: QPixmap | None = None
        self._msa_layer_cache_scratch: QPixmap | None = None
        self._msa_layer_cache_key: tuple | None = None
        self._msa_layer_cache_origin = QPoint(0, 0)
        self._msa_layer_cache_rect = QRect()
        self._msa_layer_serial = 0

    def _invalidate_msa_layer_cache(self) -> None:
        """Drop the MSA residue pixmap (zoom/data — not selection)."""
        self._msa_layer_cache = None
        self._msa_layer_cache_scratch = None
        self._msa_layer_cache_key = None
        self._msa_layer_cache_origin = QPoint(0, 0)
        self._msa_layer_cache_rect = QRect()
        self._msa_layer_serial = int(self._msa_layer_serial) + 1

    def _msa_layer_cache_identity_key(self, vis: QRect) -> tuple | None:
        """Identity that forces a full rebuild (not scroll position)."""
        if not self._msa_mode:
            return None
        if vis.isEmpty() or vis.width() < 2 or vis.height() < 2:
            return None
        dpr = max(1.0, float(self.devicePixelRatioF()))
        return (
            round(dpr, 3),
            int(vis.width()),
            int(vis.height()),
            round(float(self._ppb), 8),
            int(self._msa_layer_serial),
            int(self.width()),
            int(self.height()),
        )

    def _msa_layer_cover_size(self, vis: QRect) -> tuple[int, int]:
        """Logical cache size for this viewport (viewport + pad on each side)."""
        # One full viewport of pad → many scroll steps are pure blit.
        pad_x = max(64, int(vis.width()))
        pad_y = max(64, int(vis.height()))
        cw = min(self.width(), int(vis.width()) + 2 * pad_x)
        ch = min(self.height(), int(vis.height()) + 2 * pad_y)
        return max(1, cw), max(1, ch)

    def _msa_layer_overscroll_rect(self, vis: QRect) -> QRect:
        """Fixed-size cover centered on the viewport (clamped to content)."""
        cw, ch = self._msa_layer_cover_size(vis)
        if cw < 1 or ch < 1:
            return QRect()
        x0 = int(vis.center().x() - cw / 2)
        y0 = int(vis.center().y() - ch / 2)
        x0 = max(0, min(x0, max(0, self.width() - cw)))
        y0 = max(0, min(y0, max(0, self.height() - ch)))
        return QRect(x0, y0, cw, ch)

    def _msa_layer_cache_covers(self, vis: QRect) -> bool:
        cached = self._msa_layer_cache_rect
        if cached.isEmpty() or vis.isEmpty():
            return False
        return (
            cached.left() <= vis.left()
            and cached.top() <= vis.top()
            and cached.right() >= vis.right()
            and cached.bottom() >= vis.bottom()
        )

    def _msa_layer_edge_margin(self, vis: QRect) -> int:
        """Prefetch when the viewport comes this close to a cache edge."""
        return max(48, int(vis.height() * 0.2), int(vis.width() * 0.2))

    def _msa_layer_desired_cover(self, vis: QRect) -> QRect:
        """Keep the current cover while the viewport has comfortable slack."""
        cached = self._msa_layer_cache_rect
        ideal = self._msa_layer_overscroll_rect(vis)
        if ideal.isEmpty():
            return QRect()
        if (
            not cached.isEmpty()
            and cached.size() == ideal.size()
            and self._msa_layer_cache_covers(vis)
        ):
            margin = self._msa_layer_edge_margin(vis)
            slack_l = vis.left() - cached.left()
            slack_t = vis.top() - cached.top()
            slack_r = cached.right() - vis.right()
            slack_b = cached.bottom() - vis.bottom()
            if (
                slack_l >= margin
                and slack_t >= margin
                and slack_r >= margin
                and slack_b >= margin
            ):
                return QRect(cached)
        return ideal

    def _alloc_msa_layer_pixmap(self, cw: int, ch: int, dpr: float) -> QPixmap:
        pix = QPixmap(int(math.ceil(cw * dpr)), int(math.ceil(ch * dpr)))
        pix.setDevicePixelRatio(float(dpr))
        return pix

    def _msa_layer_pixmap_matches(
        self, pix: QPixmap | None, *, cw: int, ch: int, dpr: float
    ) -> bool:
        if pix is None or pix.isNull():
            return False
        return (
            abs(pix.width() / max(dpr, 1e-6) - cw) <= 1.5
            and abs(pix.height() / max(dpr, 1e-6) - ch) <= 1.5
            and abs(float(pix.devicePixelRatio()) - dpr) <= 0.01
        )

    def _rebuild_msa_layer_cache(self, identity: tuple, cover: QRect) -> bool:
        """Paint the full overscroll cover from scratch."""
        dpr = float(identity[0])
        cw = int(cover.width())
        ch = int(cover.height())
        if cw < 1 or ch < 1:
            return False
        if self._msa_layer_pixmap_matches(
            self._msa_layer_cache, cw=cw, ch=ch, dpr=dpr
        ):
            cache = self._msa_layer_cache
            assert cache is not None
        else:
            cache = self._alloc_msa_layer_pixmap(cw, ch, dpr)
        cache.fill(_msa_canvas_color(self))
        painter = QPainter(cache)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, False)
        painter.translate(-int(cover.x()), -int(cover.y()))
        self._paint_msa_rows(painter, cover, paint_feature_selection=False)
        painter.end()
        self._msa_layer_cache = cache
        self._msa_layer_cache_key = identity
        self._msa_layer_cache_origin = QPoint(int(cover.x()), int(cover.y()))
        self._msa_layer_cache_rect = QRect(cover)
        return True

    def _try_scroll_msa_layer_cache(self, new_cover: QRect) -> bool:
        """Shift retained pixels and paint only newly exposed strips.

        Returns False when a full rebuild is cheaper / required.
        """
        old = self._msa_layer_cache_rect
        cache = self._msa_layer_cache
        if (
            cache is None
            or cache.isNull()
            or old.isEmpty()
            or new_cover.isEmpty()
            or old.size() != new_cover.size()
        ):
            return False
        overlap = old.intersected(new_cover)
        if overlap.isEmpty():
            return False
        # Large jumps: full rebuild is clearer than copying a tiny sliver.
        old_area = max(1, old.width() * old.height())
        if overlap.width() * overlap.height() < int(0.35 * old_area):
            return False

        dx = int(old.x() - new_cover.x())
        dy = int(old.y() - new_cover.y())
        if dx == 0 and dy == 0:
            return True

        dpr = float(cache.devicePixelRatio())
        cw = int(new_cover.width())
        ch = int(new_cover.height())
        bg = _msa_canvas_color(self)
        if self._msa_layer_pixmap_matches(
            self._msa_layer_cache_scratch, cw=cw, ch=ch, dpr=dpr
        ):
            scratch = self._msa_layer_cache_scratch
            assert scratch is not None
        else:
            scratch = self._alloc_msa_layer_pixmap(cw, ch, dpr)
            self._msa_layer_cache_scratch = scratch
        painter = QPainter(scratch)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, False)
        # No full-clear: every pixel is either shifted from ``cache`` or in
        # the exposed strip painted below (avoids a Retina-wide fill hitch).
        painter.drawPixmap(QPoint(dx, dy), cache)
        exposed = QRegion(new_cover).subtracted(QRegion(old))
        painter.translate(-int(new_cover.x()), -int(new_cover.y()))
        for rect in exposed:
            if rect.isNull() or rect.width() <= 0 or rect.height() <= 0:
                continue
            painter.fillRect(rect, bg)
            self._paint_msa_rows(painter, rect, paint_feature_selection=False)
        painter.end()

        # Ping-pong: scratch becomes the live cache; old cache is reusable.
        self._msa_layer_cache_scratch = cache
        self._msa_layer_cache = scratch
        self._msa_layer_cache_origin = QPoint(int(new_cover.x()), int(new_cover.y()))
        self._msa_layer_cache_rect = QRect(new_cover)
        return True

    def _ensure_msa_layer_cache(self) -> bool:
        """Lazy sticky overscrolled residue snapshot without selection chrome."""
        vis = self.scroll_viewport_rect()
        identity = self._msa_layer_cache_identity_key(vis)
        if identity is None:
            return False
        have_cache = (
            self._msa_layer_cache is not None
            and not self._msa_layer_cache.isNull()
            and self._msa_layer_cache_key == identity
        )
        if have_cache and self._msa_layer_cache_covers(vis):
            desired = self._msa_layer_desired_cover(vis)
            if desired.isEmpty() or desired == self._msa_layer_cache_rect:
                # Comfortable slack — blit only, no residue work this frame.
                return True
            # Near an edge: prefetch the recentered cover while still covered.
            if self._try_scroll_msa_layer_cache(desired):
                return True
            return True
        cover = self._msa_layer_desired_cover(vis)
        if cover.isEmpty():
            cover = self._msa_layer_overscroll_rect(vis)
        if cover.isEmpty():
            return False
        if have_cache and self._try_scroll_msa_layer_cache(cover):
            if self._msa_layer_cache_covers(vis):
                return True
        return self._rebuild_msa_layer_cache(identity, cover)

    def scroll_viewport_rect(self) -> QRect:
        """This canvas's content rect currently shown in the scroll viewport."""
        parent = self.parent()
        if parent is None:
            vis = self.visibleRegion().boundingRect()
            return vis if not vis.isEmpty() else self.rect()
        area = parent.parent()
        if not isinstance(area, QScrollArea):
            vis = self.visibleRegion().boundingRect()
            return vis if not vis.isEmpty() else self.rect()
        vp = area.viewport()
        if vp is None:
            return self.rect()
        top_left = self.mapFrom(vp, QPoint(0, 0))
        return QRect(top_left, vp.size()).intersected(self.rect())

    def is_msa_mode(self) -> bool:
        return bool(self._msa_mode)

    def msa_sequences(self) -> list[str]:
        return list(self._msa_seqs)

    def set_data(
        self,
        record_ids: list[str],
        lengths: dict[str, int],
        features: dict[str, list[SequenceFeature]],
        *,
        sequences: list[str] | None = None,
        fit_bars_width: int | None = None,
    ) -> None:
        self._msa_mode = False
        self._msa_seqs = list(sequences or [])
        self._msa_protein = _guess_msa_is_protein(self._msa_seqs)
        self._msa_colors = (
            amino_acid_color_map()
            if self._msa_protein
            else nucleotide_base_color_map("multi")
        )
        self._msa_qcolors = {
            key: (QColor(bg), QColor(fg))
            for key, (bg, fg) in self._msa_colors.items()
        }
        self._msa_col_freqs = []
        self._clear_drag_paint_cache()
        self._record_ids = list(record_ids)
        self._original_index = list(range(len(self._record_ids)))
        self._lengths = dict(lengths)
        self._features = dict(features)
        self._rebuild_packed_tracks()
        self._max_length = max((lengths.get(rid, 1) for rid in record_ids), default=1)
        self._selected_rows = set()
        self._msa_diff_standard_id = None
        self._selected_feats = set()
        self._hover_feat = None
        self._region_sel = None
        self._regional_rows = None
        self._drag_preview_rows = None
        self._rubber = QRect()
        self._rubber_active = False
        self._drag_origin = None
        self._drag_kind = None
        self._press_feat = None
        self._chrome.clear()
        self._reset_edge_expand_phase()
        if fit_bars_width and fit_bars_width > 40 and self._max_length > 0:
            self._ppb = self.fit_pixels_per_base(fit_bars_width)
        fit_ppb, max_ppb = self.zoom_bounds()
        self._ppb = max(fit_ppb, min(max_ppb, self._ppb))
        self._recompute_layout()
        self.update()
        self.zoomChanged.emit(self._ppb)
        self.regionSelectionChanged.emit(None)
        self.rowsChanged.emit()

    def set_msa_data(
        self,
        record_ids: list[str],
        sequences: list[str],
        *,
        features: dict[str, list[SequenceFeature]] | None = None,
        fit_bars_width: int | None = None,
    ) -> None:
        """Show aligned residue rows (same length), optionally with feature tracks."""
        if len(record_ids) != len(sequences):
            raise ValueError("record_ids and sequences length mismatch")
        self._msa_mode = True
        self._msa_seqs = list(sequences)
        self._msa_protein = _guess_msa_is_protein(self._msa_seqs)
        self._msa_colors = (
            amino_acid_color_map()
            if self._msa_protein
            else nucleotide_base_color_map("msa")
        )
        self._rebuild_msa_qcolor_cache()
        self._rebuild_msa_conservation()
        self._clear_drag_paint_cache()
        self._record_ids = list(record_ids)
        self._original_index = list(range(len(self._record_ids)))
        width = len(sequences[0]) if sequences else 1
        self._lengths = {rid: width for rid in record_ids}
        feat_map = dict(features or {})
        self._features = {
            rid: list(feat_map.get(rid) or []) for rid in record_ids
        }
        self._rebuild_packed_tracks()
        self._max_length = max(1, width)
        self._selected_rows = set()
        self._msa_diff_standard_id = None
        self._selected_feats = set()
        self._hover_feat = None
        self._region_sel = None
        self._regional_rows = None
        self._drag_preview_rows = None
        self._rubber = QRect()
        self._rubber_active = False
        self._drag_origin = None
        self._drag_kind = None
        self._press_feat = None
        self._chrome.clear()
        self._reset_edge_expand_phase()
        if fit_bars_width and fit_bars_width > 40 and self._max_length > 0:
            self._ppb = self.fit_pixels_per_base(fit_bars_width)
        fit_ppb, max_ppb = self.zoom_bounds()
        self._ppb = max(fit_ppb, min(max_ppb, self._ppb))
        self._recompute_layout()
        self.update()
        self.zoomChanged.emit(self._ppb)
        self.regionSelectionChanged.emit(None)
        self.rowsChanged.emit()

    def set_features_in_place(
        self, features: dict[str, list[SequenceFeature]]
    ) -> None:
        """Update features without leaving MSA/bar mode or clearing selection."""
        self._features = {
            rid: list(features.get(rid) or []) for rid in self._record_ids
        }
        self._rebuild_packed_tracks()
        self._hover_feat = None
        self._clear_drag_paint_cache()
        self._recompute_layout()
        self.update()

    def _features_for_display(self, rid: str) -> list[SequenceFeature]:
        return filter_features_by_annotation_visibility(
            list(self._features.get(rid) or []),
            self._visible_annotation_kinds,
        )

    def _rebuild_packed_tracks(self) -> None:
        self._tracks = {
            rid: _pack_feature_tracks(self._features_for_display(rid))
            for rid in self._record_ids
        }

    def set_show_annotation_track(self, enabled: bool) -> None:
        self._show_annotation_track = bool(enabled)
        self._clear_drag_paint_cache()
        self._recompute_layout()
        self.update()

    def show_annotation_track(self) -> bool:
        return bool(self._show_annotation_track)

    def set_visible_annotation_kinds(
        self, keys: frozenset[str] | set[str] | None
    ) -> None:
        """Limit painted tracks to these visibility keys (None = all kinds)."""
        if keys is None:
            new_keys: frozenset[str] | None = None
        else:
            new_keys = frozenset(keys)
        if new_keys == self._visible_annotation_kinds:
            return
        self._visible_annotation_kinds = new_keys
        self._selected_feats = set()
        self._hover_feat = None
        self._rebuild_packed_tracks()
        self._clear_drag_paint_cache()
        self._recompute_layout()
        self.update()

    def visible_annotation_kinds(self) -> frozenset[str] | None:
        return self._visible_annotation_kinds

    def set_msa_residue_size(self, point_size: int) -> None:
        """Legacy hook: MSA letter size is fixed at ``_MSA_LETTER_SIZE_PT``.

        Overview vertical sizing uses ``set_track_height`` (1–10). Kept so older
        callers/tests do not break; values other than 10 are ignored for layout
        except as a temporary override of the letter font floor.
        """
        pt = max(2, min(36, int(point_size)))
        self._msa_font_pt = pt
        self._apply_track_height_to_residue_band()
        if self._msa_mode or self._msa_seqs:
            self._clear_drag_paint_cache()
            self._recompute_layout()
            self.update()

    def set_track_height(self, level: int) -> None:
        """Others → Track height (1–10).

        Level 10 matches the squeezed-row height at DNA/protein size 10
        (zoom-slider 3/4). Level 1 is the thinnest overview track.
        """
        new = max(1, min(10, int(level)))
        if new == self._track_height_level:
            return
        self._track_height_level = new
        self._msa_font_pt = _MSA_LETTER_SIZE_PT
        self._apply_track_height_to_residue_band()
        if self._msa_mode or self._msa_seqs or self._record_ids:
            self._clear_drag_paint_cache()
            self._recompute_layout()
            self.update()

    def track_height(self) -> int:
        return int(self._track_height_level)

    def letter_size_anchor_metrics(self) -> tuple[int, float]:
        """Band px and ppb for DNA/protein size 10 (slider 3/4 target)."""
        return _msa_letter_size_metrics(_MSA_LETTER_SIZE_PT)

    def _apply_track_height_to_residue_band(self) -> None:
        full_h, _ = _msa_letter_size_metrics(_MSA_LETTER_SIZE_PT)
        lo = max(2, min(4, max(2, full_h // 10)))
        level = max(1, min(10, int(self._track_height_level)))
        self._msa_residue_height = max(
            lo, int(round(lo + (full_h - lo) * (level - 1) / 9.0))
        )

    def msa_residue_height(self) -> int:
        return int(self._msa_residue_height)

    def _msa_shows_letters(self) -> bool:
        ppb = max(self._ppb, 1e-12)
        return bool(self._msa_mode) and ppb >= _MSA_LETTER_MIN_PPB and ppb >= 6

    def _bar_shows_letters(self) -> bool:
        ppb = max(self._ppb, 1e-12)
        has_rows = bool(not self._msa_mode and self._msa_seqs)
        return has_rows and ppb >= _MSA_LETTER_MIN_PPB and ppb >= 6

    def _row_letter_band_height(self) -> int:
        """Sequence band height: Track height when squeezed; square up to size-10."""
        if not (self._msa_shows_letters() or self._bar_shows_letters()):
            return max(1, int(self._msa_residue_height))
        letter_h, _ = _msa_letter_size_metrics(_MSA_LETTER_SIZE_PT)
        # Cap at size-10 band so further zoom only widens horizontal gaps.
        return max(1, min(letter_h, int(round(self._ppb))))

    def _msa_letter_band_height(self) -> int:
        """Residue-band height; square cells up to size-10, then fixed height."""
        if not self._msa_shows_letters():
            return max(1, int(self._msa_residue_height))
        letter_h, _ = _msa_letter_size_metrics(_MSA_LETTER_SIZE_PT)
        return max(1, min(letter_h, int(round(self._ppb))))

    def _msa_letter_font_for_cell(self, cell_w: float, cell_h: float) -> QFont:
        """Font sized to the letter background square (not the full cell)."""
        font = QFont(self.font())
        font.setFamily("Menlo")
        font.setStyleHint(QFont.StyleHint.Monospace)
        # Fit inside the square fill (~82% of bg side).
        bg = min(float(cell_w), float(cell_h)) * _MSA_LETTER_BG_SCALE
        target = max(2.0, bg * 0.82)
        lo = max(2, int(self._msa_font_pt), _MSA_LETTER_SIZE_PT)
        hi = min(96, max(lo, int(math.ceil(target)) + 2))
        best = lo
        while lo <= hi:
            mid = (lo + hi) // 2
            font.setPointSize(mid)
            fm = QFontMetrics(font)
            glyph = float(max(fm.height(), fm.ascent() + fm.descent()))
            if glyph <= target:
                best = mid
                lo = mid + 1
            else:
                hi = mid - 1
        font.setPointSize(best)
        return font

    def effective_letter_point_size(self) -> int:
        """Letter font at current zoom; floor is DNA/protein size 10 for MSA.

        At Fit (squeezed), returns Track height so names can use height+1.
        """
        base = max(2, int(self._msa_font_pt), _MSA_LETTER_SIZE_PT)
        if not (self._msa_shows_letters() or self._bar_shows_letters()):
            return max(2, int(self._track_height_level))
        cell_w = max(float(self._ppb), 1e-12)
        cell_h = float(
            self._msa_letter_band_height()
            if self._msa_mode
            else self._row_letter_band_height()
        )
        return max(base, int(self._msa_letter_font_for_cell(cell_w, cell_h).pointSize()))

    def set_sequence_row_gap(self, gap: int) -> None:
        """Vertical gap between sequences (Others → Sequence distance)."""
        self._sequence_row_gap = max(0, min(40, int(gap)))
        self._recompute_layout()
        self.update()

    def sequence_row_gap(self) -> int:
        return int(self._sequence_row_gap)

    def set_feature_bar_height(self, height: int) -> None:
        """Vertical height of annotation feature bars (Others → Annotation thickness)."""
        new = max(2, min(40, int(height)))
        if new == self._feat_height:
            return
        self._feat_height = new
        self._recompute_layout()
        self.update()

    def feature_bar_height(self) -> int:
        return int(self._feat_height)

    def set_feature_gap(self, gap: int) -> None:
        """Vertical gap between annotation rows (Others → Annotation distance)."""
        new = max(0, min(40, int(gap)))
        if new == self._feat_gap:
            return
        self._feat_gap = new
        self._recompute_layout()
        self.update()

    def feature_gap(self) -> int:
        return int(self._feat_gap)

    def set_feature_label_point_size(self, point_size: int) -> None:
        new = max(0, min(24, int(point_size)))
        if new == self._feat_label_point_size:
            return
        self._feat_label_point_size = new
        self.update()

    def feature_label_point_size(self) -> int:
        return int(self._feat_label_point_size)

    def set_column_distance(self, gap: int) -> None:
        """Horizontal gap between MSA residue backgrounds (Others → Column distance)."""
        new = max(0, min(40, int(gap)))
        if new == self._column_distance:
            return
        self._column_distance = new
        self._clear_drag_paint_cache()
        self.update()

    def column_distance(self) -> int:
        return int(self._column_distance)

    def set_color_unaligned_background(self, enabled: bool) -> None:
        enabled = bool(enabled)
        if enabled == self._color_unaligned_background:
            return
        self._color_unaligned_background = enabled
        self._clear_drag_paint_cache()
        self.update()

    def color_unaligned_background(self) -> bool:
        return bool(self._color_unaligned_background)

    def region_selection(self) -> tuple[int, int] | None:
        return self._region_sel

    def regional_rows(self) -> set[int] | None:
        """Rows limited by an MSA regional drag, or ``None`` for full columns."""
        if self._regional_rows is None:
            return None
        return set(self._regional_rows)

    def drag_preview_rows(self) -> set[int] | None:
        """Rows under an in-progress regional drag (name-column live highlight)."""
        if self._drag_preview_rows is None:
            return None
        return set(self._drag_preview_rows)

    def name_selection_preview_rows(self) -> set[int] | None:
        """Rows under a live name-column drag (overlay chrome)."""
        if self._preview_rows is None:
            return None
        return set(self._preview_rows)

    def is_regional_selection(self) -> bool:
        return self._region_sel is not None and self._regional_rows is not None

    def set_region_selection(
        self,
        start: int | None,
        end: int | None = None,
        *,
        keep_expand_anchor: bool = False,
    ) -> None:
        """Full-height column selection (scale / feature).

        Keeps name-selected sequences so both the name column and MSA row wash
        remain under the column band. Clears a prior regional matrix box only.
        """
        if start is None:
            if self._region_sel is None and self._regional_rows is None:
                return
            self._region_sel = None
            self._regional_rows = None
            if not keep_expand_anchor:
                self._reset_edge_expand_phase()
            self._sync_selection_chrome()
            self.regionSelectionChanged.emit(None)
            return
        maximum = max(1, int(self._max_length))
        a = max(1, min(maximum, int(start)))
        b = max(1, min(maximum, int(end if end is not None else start)))
        if b < a:
            a, b = b, a
        new = (a, b)
        same_band = new == self._region_sel and self._regional_rows is None
        if same_band:
            return
        had_regional = self._regional_rows is not None
        self._region_sel = new
        self._regional_rows = None
        # Keep ``_selected_rows`` so sequence names stay highlighted.
        if not keep_expand_anchor:
            self._reset_edge_expand_phase()
        self._sync_selection_chrome()
        self.regionSelectionChanged.emit(new)
        if had_regional:
            # Names were driven by regional_rows; refresh to whole-row select.
            self.rowsChanged.emit()

    def set_regional_selection(
        self,
        start: int,
        end: int,
        rows: list[int] | set[int],
    ) -> None:
        """Keep selection as the dragged MSA box (columns × rows).

        Name-column and scale highlights follow ``rows`` and ``[start, end]``.
        Whole-sequence ``selected_rows`` stay clear so tools treat this as a
        regional matrix selection, not a full-row select.
        """
        maximum = max(1, int(self._max_length))
        a = max(1, min(maximum, int(start)))
        b = max(1, min(maximum, int(end)))
        if b < a:
            a, b = b, a
        row_set = {
            int(r) for r in rows if 0 <= int(r) < len(self._record_ids)
        }
        if not row_set:
            self.set_region_selection(a, b)
            return
        new = (a, b)
        if (
            new == self._region_sel
            and row_set == self._regional_rows
            and not self._selected_rows
        ):
            return
        self._region_sel = new
        self._regional_rows = row_set
        self._selected_rows = set()
        self._reset_edge_expand_phase()
        self._sync_selection_chrome()
        self.regionSelectionChanged.emit(new)
        self.selectionChanged.emit()
        self.rowsChanged.emit()

    def select_all_rows(self) -> None:
        n = len(self._record_ids)
        self._selected_rows = set(range(n)) if n else set()
        self._selected_feats = set()
        self._reset_edge_expand_phase()
        if self._region_sel is not None or self._regional_rows is not None:
            self._region_sel = None
            self._regional_rows = None
            self.regionSelectionChanged.emit(None)
        self._sync_selection_chrome()
        self.selectionChanged.emit()
        self.rowsChanged.emit()

    def permute_rows(self, order: list[int]) -> None:
        """Reorder all row-parallel state by ``order`` (permutation of 0..n-1)."""
        n = len(self._record_ids)
        if n <= 1:
            return
        if sorted(order) != list(range(n)):
            return
        if order == list(range(n)):
            return

        def _perm(items: list):
            return [items[i] for i in order]

        self._record_ids = _perm(self._record_ids)
        if self._msa_seqs:
            self._msa_seqs = _perm(self._msa_seqs)
        if len(self._original_index) == n:
            self._original_index = _perm(self._original_index)
        else:
            self._original_index = list(range(n))
        inv = {old: new for new, old in enumerate(order)}
        self._selected_rows = {inv[r] for r in self._selected_rows if r in inv}
        # Standard is pinned by record id — survives row reordering.
        if self._regional_rows is not None:
            self._regional_rows = {inv[r] for r in self._regional_rows if r in inv}
        self._selected_feats = {
            (inv[r], f) for r, f in self._selected_feats if r in inv
        }
        self._find_hits = {
            rid: list(hits)
            for rid, hits in self._find_hits.items()
            if rid in set(self._record_ids)
        }
        self._hover_feat = None
        self._clear_drag_paint_cache()
        self._recompute_layout()
        self._sync_selection_chrome()
        self.update()
        self.rowsChanged.emit()
        self.selectionChanged.emit()

    def sort_rows(self, mode: str) -> None:
        """Sort display rows: ``original``, ``name``, or ``selected_columns``."""
        n = len(self._record_ids)
        if n <= 1:
            return
        key_mode = str(mode or "original").strip().lower()
        indices = list(range(n))
        if key_mode == "name":
            indices.sort(key=lambda i: self._record_ids[i].casefold())
        elif key_mode == "selected_columns":
            if not self._msa_mode or not self._msa_seqs:
                return
            region = self._region_sel
            if region is None:
                return
            start, end = region
            c0 = max(0, int(start) - 1)
            c1 = max(c0, int(end))

            def _col_key(i: int) -> tuple:
                seq = self._msa_seqs[i] if i < len(self._msa_seqs) else ""
                chunk = seq[c0:c1]
                # Gaps sort after letters so conserved blocks group first.
                return tuple(ch if ch not in "-." else "~" for ch in chunk)

            indices.sort(key=_col_key)
        else:
            # Original load order.
            if len(self._original_index) == n:
                indices.sort(key=lambda i: self._original_index[i])
            else:
                return
        self.permute_rows(indices)

    def _nucleotide_color_context(self) -> str:
        return "msa" if self._msa_mode else "multi"

    def refresh_msa_colors(self) -> None:
        """Reload residue palettes after Display color settings change."""
        ctx = self._nucleotide_color_context()
        invalidate_residue_color_maps(ctx)
        self._msa_colors = (
            amino_acid_color_map()
            if self._msa_protein
            else nucleotide_base_color_map(ctx)
        )
        self._rebuild_msa_qcolor_cache()
        self._clear_drag_paint_cache()
        # Nucleotide/protein palette changes must redraw the cached MSA layer;
        # otherwise colors only update after toggling MSA color scheme.
        self._invalidate_msa_layer_cache()
        self.update()

    def set_msa_color_mode(self, mode: str) -> None:
        """Display → Color MSA painting mode (nucleotide and protein)."""
        new = normalize_msa_color_mode(mode)
        if new == self._msa_color_mode:
            return
        self._msa_color_mode = new
        # Keep an explicit standard pin across mode switches for this MSA.
        if self._msa_mode:
            self._clear_drag_paint_cache()
            self._invalidate_msa_layer_cache()
            self.update()

    def set_msa_diff_standard_row(self, row: int | None) -> None:
        """Pin the Color vs standard reference by row (``None`` clears the pin)."""
        if row is None:
            self.set_msa_diff_standard_id(None)
            return
        r = int(row)
        if not (0 <= r < len(self._record_ids)):
            return
        self.set_msa_diff_standard_id(self._record_ids[r])

    def set_msa_diff_standard_id(self, record_id: str | None) -> None:
        """Pin the Color vs standard reference by sequence id."""
        new = str(record_id).strip() if record_id else None
        if new == "":
            new = None
        if new is not None and new not in self._record_ids:
            return
        if new == self._msa_diff_standard_id:
            return
        self._msa_diff_standard_id = new
        if self._msa_mode:
            self._clear_drag_paint_cache()
            self._invalidate_msa_layer_cache()
            self.update()

    def msa_color_mode(self) -> str:
        return str(self._msa_color_mode)

    def set_msa_divergence_threshold(self, value: float) -> None:
        """Fraction (0–1) of a column a letter must occupy to count as conserved."""
        new = normalize_msa_divergence_threshold(value)
        if abs(new - self._msa_divergence_threshold) < 1e-9:
            return
        self._msa_divergence_threshold = new
        if self._msa_mode:
            self._clear_drag_paint_cache()
            self._invalidate_msa_layer_cache()
            self.update()

    def msa_divergence_threshold(self) -> float:
        return float(self._msa_divergence_threshold)

    def _rebuild_msa_conservation(self) -> None:
        if not self._msa_seqs:
            self._msa_col_freqs = []
            return
        try:
            self._msa_col_freqs = column_residue_frequencies(self._msa_seqs)
        except ValueError:
            self._msa_col_freqs = [{} for _ in range(len(self._msa_seqs[0]))]

    def _rebuild_msa_qcolor_cache(self) -> None:
        cache: dict[str, tuple[QColor, QColor]] = {}
        for key, pair in self._msa_colors.items():
            bg = QColor(pair[0])
            fg = QColor(pair[1])
            if bg.isValid() and fg.isValid():
                cache[key] = (bg, fg)
                cache[key.upper()] = (bg, fg)
                cache[key.lower()] = (bg, fg)
        self._msa_qcolors = cache

    def msa_diff_standard_id(self) -> str | None:
        """Pinned Color vs standard record id, or ``None`` if unset."""
        rid = self._msa_diff_standard_id
        if rid and rid in self._record_ids:
            return str(rid)
        return None

    def msa_diff_standard_row(self) -> int | None:
        """Active Color vs standard row (pinned, else first sequence)."""
        return self._msa_standard_row()

    def msa_diff_standard_is_pinned(self) -> bool:
        """True when the user set an explicit standard for this MSA."""
        return self.msa_diff_standard_id() is not None

    def _msa_standard_row(self) -> int | None:
        """Pinned Color vs standard row, else first sequence while in that mode."""
        from repeatviewer.backends.msa_coordinates import MSA_COLOR_DIFF_FROM_SELECTED

        if self._msa_color_mode != MSA_COLOR_DIFF_FROM_SELECTED:
            return None
        if not self._msa_seqs or not self._record_ids:
            return None
        rid = self._msa_diff_standard_id
        if rid:
            try:
                return int(self._record_ids.index(rid))
            except ValueError:
                self._msa_diff_standard_id = None
        # No explicit standard yet → always the first sequence in this MSA.
        return 0

    def _msa_standard_residue(self, col: int) -> str | None:
        std = self._msa_standard_row()
        if std is None:
            return None
        seq = self._msa_seqs[std]
        if not (0 <= col < len(seq)):
            return None
        return seq[col]

    def _msa_should_fill_background(
        self, col: int, ch: str, *, row: int | None = None
    ) -> bool:
        """Whether the residue uses the scheme colour (vs overview light-grey).

        Applies at overview and letter zoom. Gaps fill only in Color vs
        selected when they differ from the standard.
        """
        from repeatviewer.backends.msa_coordinates import MSA_COLOR_DIFF_FROM_SELECTED

        if ch in "-." and self._msa_color_mode != MSA_COLOR_DIFF_FROM_SELECTED:
            return False
        freqs = (
            self._msa_col_freqs[col]
            if 0 <= col < len(self._msa_col_freqs)
            else None
        )
        std = self._msa_standard_row()
        is_ref = row is not None and std is not None and row == std
        reference = self._msa_standard_residue(col) if std is not None else None
        return msa_residue_fill_background(
            self._msa_color_mode,
            residue=ch,
            column_freqs=freqs,
            threshold=self._msa_divergence_threshold,
            reference_residue=reference,
            is_reference_row=is_ref,
        )

    def _msa_residue_paint_pair(self, col: int, ch: str) -> tuple[QColor, QColor] | None:
        """Background/foreground colours for a non-gap residue cell."""
        if ch in "-.":
            return None
        return self._msa_qcolors.get(ch)

    def _msa_overview_fill_color(
        self, col: int, ch: str, *, row: int | None = None
    ) -> QColor | None:
        """Solid fill for a zoomed-out MSA residue, or ``None`` for a gap.

        Scheme-coloured cells use the residue palette (or conservation
        greyscale). Cells the scheme leaves uncoloured use the same light
        grey as multi-sequence overview bars. Gap-vs-standard diffs use a
        dedicated accent fill.
        """
        if ch in "-.":
            if self._msa_should_fill_background(col, ch, row=row):
                return _msa_gap_diff_color(self)
            return None
        if self._msa_should_fill_background(col, ch, row=row):
            pair = self._msa_residue_paint_pair(col, ch)
            if pair is not None:
                return pair[0]
            return _msa_plain_bar_color(self)
        return _msa_plain_bar_color(self)

    def _clear_drag_paint_cache(self) -> None:
        """Drop live selection preview state (name kept for existing call sites)."""
        self._preview_columns = None
        self._preview_rows = None

    def _column_selection_x_range(self, start: int, end: int) -> tuple[float, float]:
        """Content-x span for columns ``start``–``end`` (1-based inclusive).

        Point features (start==end) and tiny zooms expand to at least
        ``_MIN_COLUMN_BAND_PX`` so dashed guides stay visible.
        """
        a, b = int(start), int(end)
        if b < a:
            a, b = b, a
        x0 = float(self.bp_to_x(a - 1))
        x1 = float(self.bp_to_x(b))
        if x1 < x0:
            x0, x1 = x1, x0
        min_w = max(_MIN_COLUMN_BAND_PX, float(self._ppb))
        if (x1 - x0) < min_w - 1e-9:
            mid = 0.5 * (x0 + x1)
            x0 = mid - 0.5 * min_w
            x1 = mid + 0.5 * min_w
        return x0, x1

    def _column_band_rect(self, start: int, end: int) -> QRect:
        x0, x1 = self._column_selection_x_range(start, end)
        left = int(math.floor(x0))
        right = int(math.ceil(x1))
        return QRect(left, 0, max(1, right - left), max(1, self.height()))

    def _rows_band_rect(self, rows: set[int]) -> QRect:
        """Bounding box of sequence bands for ``rows`` (excludes annotation tracks)."""
        boxes = self._row_selection_boxes(rows)
        if not boxes:
            return QRect()
        union = QRect(boxes[0])
        for b in boxes[1:]:
            union = union.united(b)
        return union

    def _row_selection_boxes(self, rows: set[int] | list[int]) -> list[QRect]:
        """One full-width rect per contiguous name-selected sequence block."""
        sorted_rows = sorted(
            r for r in rows if 0 <= int(r) < len(self._row_tops)
        )
        if not sorted_rows:
            return []
        width = max(1, self.width())
        boxes: list[QRect] = []
        i = 0
        while i < len(sorted_rows):
            start = sorted_rows[i]
            end = start
            while i + 1 < len(sorted_rows) and sorted_rows[i + 1] == end + 1:
                i += 1
                end = sorted_rows[i]
            y0 = int(self._row_tops[start])
            y1 = int(self._sequence_band_bottom(end))
            boxes.append(QRect(0, y0, width, max(1, y1 - y0)))
            i += 1
        return boxes

    def _committed_regional_box(
        self, start: int, end: int, rows: set[int]
    ) -> QRect:
        sorted_rows = sorted(r for r in rows if 0 <= r < len(self._row_tops))
        if not sorted_rows:
            return QRect()
        x0 = int(self.bp_to_x(start - 1))
        x1 = int(math.ceil(self.bp_to_x(end)))
        if x1 < x0:
            x0, x1 = x1, x0
        y0 = int(self._row_tops[sorted_rows[0]])
        y1 = int(self._sequence_band_bottom(sorted_rows[-1]))
        return QRect(x0, y0, max(1, x1 - x0), max(1, y1 - y0))

    def _sync_selection_chrome(self) -> None:
        """Committed selection chrome (column wash + optional row wash overlay)."""
        if (
            self._preview_columns is not None
            or self._preview_rows is not None
            or self._rubber_active
        ):
            return
        if self._region_sel is not None:
            start, end = self._region_sel
            if self._regional_rows is not None:
                box = self._committed_regional_box(start, end, self._regional_rows)
                if box.isNull():
                    self._chrome.clear()
                else:
                    self._chrome.set_region_box(box)
            else:
                self._chrome.set_column_band(self._column_band_rect(start, end))
                # Keep sequence-row highlighting on the MSA under the column band.
                row_boxes = self._row_selection_boxes(self._selected_rows)
                if row_boxes:
                    self._chrome.set_row_boxes(row_boxes)
                else:
                    self._chrome.set_row_boxes([])
            return
        boxes = self._row_selection_boxes(self._selected_rows)
        if not boxes:
            self._chrome.clear()
        else:
            self._chrome.set_row_boxes(boxes)

    def _sync_row_selection_chrome(self) -> None:
        self._sync_selection_chrome()

    def _dismiss_committed_selection_visuals(self) -> None:
        """Drop prior selection at drag start — chrome clear only."""
        had_rows = bool(self._selected_rows)
        had_region = self._region_sel is not None or self._regional_rows is not None
        self._region_sel = None
        self._regional_rows = None
        self._selected_rows = set()
        self._selected_feats = set()
        self._temp_matrix_sel = None
        self._reset_edge_expand_phase()
        self._chrome.clear()
        # Clear scale band immediately; do not emit regionSelectionChanged(None)
        # (that would also reset Position spins mid-drag).
        self.selectionVisualsDismissed.emit()
        # Keep name column in sync when starting a different selection kind.
        if had_rows or had_region:
            self.rowsChanged.emit()

    def _begin_column_selection_preview(self) -> None:
        """Start live column guides without clearing name-selected sequences."""
        had_regional = self._regional_rows is not None
        self._region_sel = None
        self._regional_rows = None
        self._selected_feats = set()
        self._temp_matrix_sel = None
        self._preview_rows = None
        self._reset_edge_expand_phase()
        # Keep sequence-row wash on the MSA while column guides follow the drag.
        row_boxes = self._row_selection_boxes(self._selected_rows)
        self._chrome.clear()
        if row_boxes:
            self._chrome.set_row_boxes(row_boxes)
        self.selectionVisualsDismissed.emit()
        if had_regional:
            self.rowsChanged.emit()

    def set_column_selection_preview(self, start: int, end: int) -> None:
        """Committed-style column preview from 1-based columns (snapped)."""
        maximum = max(1, int(self._max_length))
        a = max(1, min(maximum, int(start)))
        b = max(1, min(maximum, int(end)))
        if b < a:
            a, b = b, a
        new = (a, b)
        if new == self._preview_columns:
            return
        if self._preview_columns is None:
            self._begin_column_selection_preview()
        self._preview_columns = new
        self._preview_rows = None
        self._chrome.set_column_band(self._column_band_rect(a, b), live=True)

    def set_column_pixel_preview(self, x0: float, x1: float) -> None:
        """Live column guides; cell-snapped at nucleotide zoom, else pixel-follow."""
        starting = self._preview_columns is None and not self._chrome._live
        if starting:
            self._begin_column_selection_preview()
        self._preview_rows = None
        if self._letter_level_selection():
            start_i, end_i = self._columns_from_rubber(
                QRect(int(min(x0, x1)), 0, max(1, int(abs(x1 - x0))), 1)
            )
            new = (start_i, end_i)
            if new == self._preview_columns and self._chrome._mode == "columns":
                return
            self._preview_columns = new
            self._chrome.set_column_band(
                self._column_band_rect(start_i, end_i), live=True
            )
            return
        a = int(min(x0, x1))
        b = int(max(x0, x1))
        if b <= a:
            b = a + 1
        # Viewport only — never the full tall content height.
        vis = self.scroll_viewport_rect()
        if vis.isEmpty():
            y0, h = 0, max(1, self.height())
        else:
            y0, h = int(vis.top()), max(1, int(vis.height()))
        band = QRect(a, y0, b - a, h)
        if (
            self._preview_columns is not None
            and self._chrome._mode == "columns"
            and self._chrome._box == band
            and self._chrome._live
        ):
            return
        # Non-None sentinel so finish/dismiss know a column preview is active.
        self._preview_columns = self._preview_columns or (1, 1)
        self._chrome.set_column_band(band, live=True)

    def set_row_selection_preview(self, rows: list[int] | set[int]) -> None:
        """Live sequence-row outline while dragging names (overlay only)."""
        row_set = {
            int(r) for r in rows if 0 <= int(r) < len(self._record_ids)
        }
        if row_set == self._preview_rows:
            return
        if self._preview_rows is None:
            self._dismiss_committed_selection_visuals()
        self._preview_rows = set(row_set)
        self._preview_columns = None
        boxes = self._row_selection_boxes(row_set)
        if not boxes:
            self._chrome.clear()
        else:
            self._chrome.set_row_boxes(boxes, live=True)

    def set_row_pixel_preview(self, y0: float, y1: float) -> None:
        """Live row outline: cell-snapped at nucleotide zoom, else pixel-follow."""
        starting = self._preview_rows is None and not self._chrome._live
        rows = self._rows_in_y_range(float(y0), float(y1))
        row_set = set(rows)
        if starting:
            had_rows = bool(self._selected_rows)
            had_region = (
                self._region_sel is not None or self._regional_rows is not None
            )
            self._region_sel = None
            self._regional_rows = None
            self._selected_rows = set()
            self._selected_feats = set()
            self._temp_matrix_sel = None
            self._reset_edge_expand_phase()
            self._chrome.clear()
            self._preview_columns = None
            self._preview_rows = row_set
            if had_rows or had_region:
                self.rowsChanged.emit()
        else:
            self._preview_columns = None
            self._preview_rows = row_set

        if self._letter_level_selection():
            if (
                not starting
                and row_set == self._live_row_snap
                and self._chrome._mode == "rows"
            ):
                return
            self._live_row_snap = set(row_set)
            boxes = self._row_selection_boxes(row_set)
            if not boxes:
                self._chrome.clear()
            else:
                self._chrome.set_row_boxes(boxes, live=True)
            return

        # Fit / overview: pixel-follow outline (no cell snap).
        self._live_row_snap = None
        a = int(min(y0, y1))
        b = int(max(y0, y1))
        if b <= a:
            b = a + 1
        box = QRect(0, a, max(1, self.width()), b - a)
        if (
            not starting
            and self._chrome._mode == "rows"
            and self._chrome._row_boxes == [box]
            and self._chrome._live
        ):
            return
        self._chrome.set_row_boxes([box], live=True)

    def finish_live_selection_preview(self) -> None:
        """End live preview flags before commit.

        Leaves chrome in place so commit can upgrade live→committed in one
        paint (clear-then-rebuild was the column/row click lag vs regional).
        Cancel paths that need a blank canvas should call ``_chrome.clear()``.
        """
        self._preview_columns = None
        self._preview_rows = None
        self._live_row_snap = None
        self._clear_drag_paint_cache()

    def _sequence_band_height(self, row: int) -> int:
        """Height of the residue/bar band (excludes annotation tracks)."""
        if not (0 <= row < len(self._row_heights)):
            return _MIN_ROW_HEIGHT
        # MSA and multi-seq overview both size the band from Track height.
        if self._msa_mode or self._msa_seqs or self._bar_shows_letters():
            return int(self._row_letter_band_height())
        if not self._show_annotation_track:
            return int(self._row_heights[row])
        return _BAR_V_PAD + _BAR_HEIGHT + _BAR_V_PAD

    def _sequence_bar_y(self, row: int) -> int:
        """Top Y of the sequence bar within its band."""
        top = self._row_tops[row] if 0 <= row < len(self._row_tops) else 0
        if self._msa_mode or self._bar_shows_letters():
            return int(top)
        band_h = self._sequence_band_height(row)
        bar_h = self._multi_seq_bar_height(row)
        # Center only when the band is taller than the bar (legacy pad layout).
        return int(top + max(0, (band_h - bar_h) // 2))

    def _multi_seq_bar_height(self, row: int) -> int:
        """Painted bar thickness for multi-seq overview (not MSA letters).

        Matches Track height so Sequence distance gaps stay visible. The old
        fixed ``_BAR_HEIGHT`` (7px) overflowed thin tracks and erased gaps.
        """
        band_h = max(1, int(self._sequence_band_height(row)))
        # When annotations use the classic padded band, keep the historic bar.
        if not (self._msa_mode or self._msa_seqs):
            return min(_BAR_HEIGHT, band_h)
        return band_h

    def _sequence_band_bottom(self, row: int) -> int:
        """Y of the residue/bar band bottom (excludes annotation tracks)."""
        if not (0 <= row < len(self._row_tops)):
            return 0
        return int(self._row_tops[row] + self._sequence_band_height(row))

    def _nearest_row_at_y(self, y: float) -> int:
        """Row for ``y`` when outside a hit band (inter-row gap or past ends).

        Gaps between sequences are not part of any row. Mapping them to the
        last row made regional drags suddenly select every sequence below.
        """
        tops = self._row_tops
        heights = self._row_heights
        n = len(tops)
        if n <= 0:
            return -1
        if y < tops[0]:
            return 0
        last = n - 1
        last_bottom = tops[last] + heights[last]
        if y >= last_bottom:
            return last
        # Binary search the gap / band containing y.
        lo, hi = 0, last
        while lo <= hi:
            mid = (lo + hi) // 2
            top = tops[mid]
            bottom = top + heights[mid]
            if y < top:
                hi = mid - 1
            elif y >= bottom:
                # In gap after ``mid``, or further down.
                if mid < last and y < tops[mid + 1]:
                    mid_dist = y - bottom
                    next_dist = tops[mid + 1] - y
                    return mid if mid_dist <= next_dist else mid + 1
                lo = mid + 1
            else:
                return mid
        return max(0, min(last, lo))

    def _pane_to_matrix(self, x: float, y: float) -> tuple[int, int]:
        """Map pane pixels → ``(column, row)`` like AliView (cell-snapped).

        Column is 1-based inclusive; row is 0-based (clamped for drag).
        """
        maximum = max(1, int(self._max_length))
        col = int(math.floor(self.x_to_bp(x))) + 1
        col = max(1, min(maximum, col))
        row = self._row_at_y(y)
        if row < 0:
            row = self._nearest_row_at_y(y)
        return col, max(0, row)

    def _matrix_selection_box(
        self, c0: int, c1: int, r0: int, r1: int
    ) -> QRect:
        """Pixel rect covering matrix selection ``[c0..c1] × [r0..r1]``."""
        if c1 < c0:
            c0, c1 = c1, c0
        if r1 < r0:
            r0, r1 = r1, r0
        n = len(self._row_tops)
        if n <= 0:
            return QRect()
        r0 = max(0, min(n - 1, r0))
        r1 = max(0, min(n - 1, r1))
        x0 = int(self.bp_to_x(c0 - 1))
        x1 = int(math.ceil(self.bp_to_x(c1)))
        y0 = int(self._row_tops[r0])
        y1 = int(self._sequence_band_bottom(r1))
        return QRect(x0, y0, max(1, x1 - x0), max(1, y1 - y0))

    def _letter_level_selection(self) -> bool:
        """True when residues are drawn as letters — snap selection to cell edges."""
        return bool(self._msa_shows_letters() or self._bar_shows_letters())

    def _schedule_region_drag_preview(self, c0: int, c1: int, r0: int, r1: int) -> None:
        """Throttle Position spins + name-column highlight during regional drag."""
        if r1 < r0:
            r0, r1 = r1, r0
        rows = list(range(int(r0), int(r1) + 1))
        pending = (int(c0), int(c1), rows)
        if self._drag_preview_pending == pending:
            return
        if (
            self._drag_preview_pending is None
            and self._drag_preview_columns == (pending[0], pending[1])
            and self._drag_preview_rows == set(rows)
        ):
            return
        self._drag_preview_pending = pending
        if not self._drag_preview_timer.isActive():
            self._flush_region_drag_preview()
            self._drag_preview_timer.start()

    def _set_live_region_rubber(self, rubber: QRect) -> None:
        """Live rubber snapped to MSA cells so the box matches selected rows."""
        rect = QRect(rubber).normalized()
        if rect.isNull() or rect.width() <= 0 or rect.height() <= 0:
            return
        # Matrix selection for commit and chrome (overview + letter zoom).
        c0, r0 = self._pane_to_matrix(rect.left(), rect.top())
        c1, r1 = self._pane_to_matrix(rect.right(), rect.bottom())
        if c1 < c0:
            c0, c1 = c1, c0
        if r1 < r0:
            r0, r1 = r1, r0
        self._temp_matrix_sel = (c0, c1, r0, r1)
        # Live Position + sequence names while dragging.
        self._schedule_region_drag_preview(c0, c1, r0, r1)
        box = self._matrix_selection_box(c0, c1, r0, r1)
        if box.isNull():
            return
        if (
            self._rubber_active
            and self._chrome._mode == "region"
            and self._chrome._box == box
            and self._chrome._live
        ):
            return
        self._rubber_active = True
        self._rubber = box
        # Outline-only + border dirty — avoids MSA recomposite under a wash.
        self._chrome.set_region_box(box, live=True)

    def _ensure_drag_point_visible(self, point: QPoint) -> None:
        """Scroll only when the pointer is outside the viewport (cheap no-op inside)."""
        vp = self.parent()
        if vp is None:
            return
        area = vp.parent()
        if not isinstance(area, QScrollArea):
            return
        hbar = area.horizontalScrollBar()
        vbar = area.verticalScrollBar()
        vx = int(point.x()) - hbar.value()
        vy = int(point.y()) - vbar.value()
        vw = max(1, vp.width())
        vh = max(1, vp.height())
        if 0 <= vx < vw and 0 <= vy < vh:
            return
        margin = 28
        if vx < 0:
            hbar.setValue(max(hbar.minimum(), hbar.value() + vx - margin))
        elif vx >= vw:
            hbar.setValue(min(hbar.maximum(), hbar.value() + (vx - vw) + margin))
        if vy < 0:
            vbar.setValue(max(vbar.minimum(), vbar.value() + vy - margin))
        elif vy >= vh:
            vbar.setValue(min(vbar.maximum(), vbar.value() + (vy - vh) + margin))

    def _rubber_paint_box(self, rubber: QRect) -> QRect:
        """Pixel box for feature/annotation rubber (not cell-snapped)."""
        rect = QRect(rubber).normalized()
        if rect.isNull() or rect.isEmpty():
            return QRect()
        return rect

    def _update_rubber_overlay(self, rubber: QRect) -> None:
        """Feature-drag chrome — outline-only while dragging."""
        self._rubber = QRect(rubber)
        self._rubber_active = True
        self._chrome.set_region_box(self._rubber_paint_box(rubber), live=True)

    def pixels_per_base(self) -> float:
        return self._ppb

    def max_length(self) -> int:
        return self._max_length

    def bar_pad(self) -> int:
        return self._bar_pad

    def row_tops(self) -> list[int]:
        return list(self._row_tops)

    def row_heights(self) -> list[int]:
        return list(self._row_heights)

    def content_height(self) -> int:
        if not self._row_tops:
            return _MIN_ROW_HEIGHT
        return self._row_tops[-1] + self._row_heights[-1] + 4

    def record_ids(self) -> list[str]:
        return list(self._record_ids)

    def selected_rows(self) -> list[int]:
        return sorted(self._selected_rows)

    def selected_row_set(self) -> set[int]:
        return set(self._selected_rows)

    def selected_feature(self) -> tuple[int, int] | None:
        """Primary ``(row, flat_feature_index)`` when exactly one feature is selected."""
        if len(self._selected_feats) == 1:
            return next(iter(self._selected_feats))
        if not self._selected_feats:
            return None
        return sorted(self._selected_feats)[0]

    def selected_features(self) -> list[tuple[int, int]]:
        return sorted(self._selected_feats)

    def clear_feature_selection(self, *, repaint: bool = True) -> None:
        if not self._selected_feats:
            return
        self._selected_feats = set()
        self._reset_edge_expand_phase()
        if repaint:
            self.update()

    def _reset_edge_expand_phase(self) -> None:
        self._edge_expand_phase = {"left": 0, "right": 0}
        self._edge_expand_anchor = None

    def set_viewport_min_size(self, width: int, height: int) -> None:
        # Ignore collapsed/hidden viewports (e.g. PDF tab on top): adopting a
        # ~40px min size leaves the MSA in a corner with white around it.
        w = max(0, int(width))
        h = max(0, int(height))
        if w < 80 or h < 40:
            return
        new = QSize(w, h)
        if new == self._min_viewport:
            fit_ppb, max_ppb = self.zoom_bounds()
            self._ppb = max(fit_ppb, min(max_ppb, self._ppb))
            return
        self._min_viewport = new
        fit_ppb, max_ppb = self.zoom_bounds()
        self._ppb = max(fit_ppb, min(max_ppb, self._ppb))
        self._recompute_layout()

    def fit_pixels_per_base(self, bars_width: int | None = None) -> float:
        """Scale that fits the longest sequence inside the track viewport."""
        width = max(
            40,
            int(self._min_viewport.width() if bars_width is None else bars_width),
        )
        usable_width = max(1, width - self._bar_pad - _BAR_RIGHT_PAD)
        return max(_ABS_MIN_PPB, usable_width / max(1, self._max_length))

    def zoom_bounds(self, bars_width: int | None = None) -> tuple[float, float]:
        """Return data-dependent zoom limits that remain safe for QWidget."""
        fit_ppb = self.fit_pixels_per_base(bars_width)
        geometry_limit = max(
            _ABS_MIN_PPB,
            (
                _QT_MAX_WIDGET_SIZE
                - self._bar_pad
                - _BAR_RIGHT_PAD
            )
            / max(1, self._max_length),
        )
        if self._msa_mode or self._msa_seqs:
            # Cap at 3× size-10 so max zoom fills columns flush (no gutters).
            _band, letter_ppb = _msa_letter_size_metrics(_MSA_LETTER_SIZE_PT)
            ceil = max(
                _MSA_MAX_PPB,
                float(letter_ppb) * _MSA_LETTER_BG_MAX_WIDTH_SCALE,
            )
            desired_max = ceil
        else:
            ceil = _MAX_PPB
            desired_max = max(ceil, fit_ppb * _MAX_ZOOM_OVER_FIT)
        return fit_ppb, max(fit_ppb, min(desired_max, geometry_limit))

    def set_pixels_per_base(self, ppb: float) -> None:
        fit_ppb, max_ppb = self.zoom_bounds()
        self._ppb = max(fit_ppb, min(max_ppb, float(ppb)))
        self._clear_drag_paint_cache()
        self._invalidate_msa_layer_cache()
        self._recompute_layout()
        self.update()
        self.zoomChanged.emit(self._ppb)

    def zoom_by(self, factor: float, *, anchor_x: float | None = None) -> None:
        old = self._ppb
        fit_ppb, max_ppb = self.zoom_bounds()
        new = max(fit_ppb, min(max_ppb, old * factor))
        if abs(new - old) < 1e-12:
            return
        if anchor_x is not None and old > 0:
            self._pending_anchor = (
                (float(anchor_x) - self._bar_pad) / old,
                float(anchor_x),
            )
        else:
            self._pending_anchor = None
        self._pending_scroll_bp = None
        self._ppb = new
        self._clear_drag_paint_cache()
        self._invalidate_msa_layer_cache()
        self._recompute_layout()
        self.update()
        self.zoomChanged.emit(self._ppb)

    def fit_to_bars_width(self, bars_width: int) -> None:
        if self._max_length <= 0:
            return
        self._pending_scroll_bp = None
        self._pending_anchor = None
        self.set_pixels_per_base(self.fit_pixels_per_base(bars_width))

    def zoom_to_bp_range(
        self, start_bp: float, end_bp: float, bars_width: int
    ) -> None:
        span = max(1.0, float(end_bp) - float(start_bp))
        width = max(40, int(bars_width))
        self._pending_scroll_bp = float(start_bp)
        self._pending_anchor = None
        usable_width = max(1, width - self._bar_pad - _BAR_RIGHT_PAD)
        self.set_pixels_per_base(usable_width / span)

    def take_pending_scroll_bp(self) -> float | None:
        value = self._pending_scroll_bp
        self._pending_scroll_bp = None
        return value

    def take_pending_anchor(self) -> tuple[float, float] | None:
        value = self._pending_anchor
        self._pending_anchor = None
        return value

    def set_find_hits(
        self,
        hits_by_id: dict[str, list[tuple[int, int]]],
        current: tuple[str, int, int] | None = None,
    ) -> None:
        """Edit > Find hits drawn on the bars ({} clears)."""
        self._find_hits = {rid: list(ranges) for rid, ranges in (hits_by_id or {}).items()}
        self._find_current = current
        self.update()

    def set_selected_rows(self, rows: list[int] | set[int]) -> None:
        new = {int(r) for r in rows if 0 <= int(r) < len(self._record_ids)}
        # Sequence-row select and column/regional select are mutually exclusive.
        cleared_region = False
        if new and (
            self._region_sel is not None or self._regional_rows is not None
        ):
            self._region_sel = None
            self._regional_rows = None
            self._selected_feats = set()
            self._reset_edge_expand_phase()
            cleared_region = True
        if new == self._selected_rows and not cleared_region:
            self._sync_row_selection_chrome()
            return
        self._selected_rows = new
        # Overlay chrome only — Color vs standard is pinned separately.
        self._sync_row_selection_chrome()
        if cleared_region:
            self.regionSelectionChanged.emit(None)
        self.selectionChanged.emit()
        self.rowsChanged.emit()

    def _msa_track_block_height(self, rid: str) -> int:
        if not self._show_annotation_track:
            return 0
        n_tracks = len(self._tracks.get(rid) or [])
        if n_tracks <= 0:
            return 0
        # Lead with Annotation distance so sequence→feature matches feature→feature.
        gap = max(0, int(self._feat_gap))
        return (
            gap
            + n_tracks * self._feat_height
            + (n_tracks - 1) * gap
        )

    def _recompute_layout(self) -> None:
        # Layout/zoom/size change → residue snapshot must be rebuilt.
        self._invalidate_msa_layer_cache()
        self._row_tops = []
        self._row_heights = []
        y = 0
        n_rows = len(self._record_ids)
        gap = max(0, int(self._sequence_row_gap))
        for i, rid in enumerate(self._record_ids):
            if self._msa_mode or self._msa_seqs:
                # MSA / multi-seq: Track height when squeezed; square band up to
                # size-10 (slider ¾). Past that, height stays fixed.
                h = int(self._row_letter_band_height()) + self._msa_track_block_height(
                    rid
                )
            else:
                n_tracks = (
                    len(self._tracks.get(rid) or [])
                    if self._show_annotation_track
                    else 0
                )
                h = _row_height(
                    n_tracks,
                    feat_height=self._feat_height,
                    feat_gap=self._feat_gap,
                )
            self._row_tops.append(y)
            self._row_heights.append(h)
            y += h
            # Empty space between sequences (not part of the row hit-target).
            if i < n_rows - 1:
                y += gap
        content_w = math.ceil(
            self._bar_pad + self._max_length * self._ppb + _BAR_RIGHT_PAD
        )
        # Extra bottom pad so the last row clears the horizontal scrollbar edge.
        y = y + 4 if y else _MIN_ROW_HEIGHT
        w = min(
            _QT_MAX_WIDGET_SIZE,
            max(120, content_w, self._min_viewport.width()),
        )
        h = min(
            _QT_MAX_WIDGET_SIZE,
            max(_MIN_ROW_HEIGHT, y, self._min_viewport.height()),
        )
        self.setMinimumSize(w, h)
        self.resize(w, h)
        self.updateGeometry()
        # Row chrome geometry tracks sequence bands; refresh after relayout.
        self._sync_row_selection_chrome()
        self.rowsChanged.emit()

    def x_to_bp(self, x: float) -> float:
        return max(0.0, (float(x) - self._bar_pad) / max(self._ppb, 1e-12))

    def bp_to_x(self, bp: float) -> float:
        return self._bar_pad + float(bp) * self._ppb

    def _row_at_y(self, y: float) -> int:
        tops = self._row_tops
        heights = self._row_heights
        n = len(tops)
        if n <= 0:
            return -1
        lo, hi = 0, n - 1
        while lo <= hi:
            mid = (lo + hi) // 2
            top = tops[mid]
            bottom = top + heights[mid]
            if y < top:
                hi = mid - 1
            elif y >= bottom:
                lo = mid + 1
            else:
                return mid
        return -1

    def _rows_in_y_range(self, y0: float, y1: float) -> list[int]:
        tops = self._row_tops
        heights = self._row_heights
        n = len(tops)
        if n <= 0:
            return []
        lo_y, hi_y = (y0, y1) if y0 <= y1 else (y1, y0)
        lo, hi = 0, n - 1
        first = n
        while lo <= hi:
            mid = (lo + hi) // 2
            if tops[mid] + heights[mid] > lo_y:
                first = mid
                hi = mid - 1
            else:
                lo = mid + 1
        lo, hi = 0, n - 1
        last = -1
        while lo <= hi:
            mid = (lo + hi) // 2
            if tops[mid] < hi_y:
                last = mid
                lo = mid + 1
            else:
                hi = mid - 1
        if first > last:
            return []
        return list(range(first, last + 1))

    def _feature_band_y0(self, row: int) -> int:
        """Top Y of the first annotation track (below sequence band + gap).

        Gap matches Others → Annotation distance (same as between feature rows).
        """
        gap = max(0, int(self._feat_gap)) if self._show_annotation_track else 0
        return self._sequence_band_bottom(row) + gap

    def _y_in_annotation_band(self, row: int, y: float) -> bool:
        if not self._show_annotation_track or row < 0 or row >= len(self._row_tops):
            return False
        top = self._row_tops[row]
        height = self._row_heights[row] if row < len(self._row_heights) else 0
        feat_y0 = self._feature_band_y0(row)
        return feat_y0 <= y < top + height

    def _hit_feature(self, row: int, x: float, y: float) -> int | None:
        if not self._show_annotation_track:
            return None
        if not (0 <= row < len(self._record_ids)):
            return None
        rid = self._record_ids[row]
        length = max(1, int(self._lengths.get(rid, 1)))
        feat_y0 = self._feature_band_y0(row)
        flat_index = 0
        for track_i, track in enumerate(self._tracks.get(rid) or []):
            ty = feat_y0 + track_i * (self._feat_height + self._feat_gap)
            for feat in track:
                start = max(1, min(int(feat.start), length))
                end = max(1, min(int(feat.end), length))
                if start > end:
                    start, end = end, start
                x0 = self._bar_pad + (start - 1) * self._ppb
                x1 = self._bar_pad + end * self._ppb
                if x0 <= x <= max(x0 + 2, x1) and ty <= y <= ty + self._feat_height:
                    return flat_index
                flat_index += 1
        return None

    def _feature_by_flat_index(
        self, row: int, flat_index: int
    ) -> SequenceFeature | None:
        if not (0 <= row < len(self._record_ids)):
            return None
        rid = self._record_ids[row]
        i = 0
        for track in self._tracks.get(rid) or []:
            for feat in track:
                if i == flat_index:
                    return feat
                i += 1
        return None

    def flat_indices_for_annotation(
        self,
        row: int,
        start: int,
        end: int,
        label: str = "",
    ) -> list[int]:
        """Flat feature indices on ``row`` matching a Find annotation hit."""
        if not (0 <= row < len(self._record_ids)):
            return []
        rid = self._record_ids[row]
        a, b = sorted((int(start), int(end)))
        needle = str(label or "").strip()
        matches: list[int] = []
        flat_i = 0
        for track in self._tracks.get(rid) or []:
            for feat in track:
                fa, fb = sorted((int(feat.start), int(feat.end)))
                if fa == a and fb == b:
                    name = str(feat.name or "")
                    ftype = str(feat.feature_type or "")
                    if (
                        not needle
                        or needle == name
                        or needle == ftype
                        or needle.lower() in name.lower()
                        or needle.lower() in ftype.lower()
                    ):
                        matches.append(flat_i)
                flat_i += 1
        return matches

    def _feature_rect(self, row: int, flat_index: int) -> QRectF | None:
        if not (0 <= row < len(self._record_ids)):
            return None
        rid = self._record_ids[row]
        length = max(1, int(self._lengths.get(rid, 1)))
        feat_y0 = self._feature_band_y0(row)
        i = 0
        for track_i, track in enumerate(self._tracks.get(rid) or []):
            ty = feat_y0 + track_i * (self._feat_height + self._feat_gap)
            for feat in track:
                if i == flat_index:
                    start = max(1, min(int(feat.start), length))
                    end = max(1, min(int(feat.end), length))
                    if start > end:
                        start, end = end, start
                    x0 = self._bar_pad + (start - 1) * self._ppb
                    x1 = self._bar_pad + end * self._ppb
                    return QRectF(x0, ty, max(2.0, x1 - x0), float(self._feat_height))
                i += 1
        return None

    def _store_index_for_flat(
        self, row: int, flat_index: int
    ) -> tuple[str, int, SequenceFeature] | None:
        """Map packed flat index → (record_id, store_index, feature).

        Prefer object identity so stacked duplicates (same coords/name, e.g.
        multiple TSD flanks) each resolve to a distinct store entry.
        """
        feat = self._feature_by_flat_index(row, flat_index)
        if feat is None or not (0 <= row < len(self._record_ids)):
            return None
        rid = self._record_ids[row]
        store = list(self._features.get(rid) or [])
        for i, candidate in enumerate(store):
            if candidate is feat:
                return rid, i, candidate
        # Fallback when tracks were rebuilt from copies of the same fields.
        for i, candidate in enumerate(store):
            if (
                int(candidate.start) == int(feat.start)
                and int(candidate.end) == int(feat.end)
                and candidate.name == feat.name
                and candidate.feature_type == feat.feature_type
                and candidate.strand == feat.strand
            ):
                return rid, i, candidate
        return None

    def _edge_hit(self, row: int, x: float, y: float) -> tuple[int, str] | None:
        """Return ``(flat_index, 'start'|'end')`` when near a feature edge."""
        flat = self._hit_feature(row, x, y)
        if flat is None:
            return None
        rect = self._feature_rect(row, flat)
        if rect is None:
            return None
        threshold = min(
            max(4.0, _FEAT_EDGE_HIT_PX),
            max(2.0, rect.width() * 0.35),
        )
        dist_left = abs(float(x) - float(rect.left()))
        dist_right = abs(float(x) - float(rect.right()))
        if dist_left <= threshold and dist_left <= dist_right:
            return flat, "start"
        if dist_right <= threshold:
            return flat, "end"
        return None

    def _x_to_feature_bp(self, x: float, *, length: int) -> int:
        bp = int(math.floor(self.x_to_bp(x))) + 1
        return max(1, min(max(1, length), bp))

    def _apply_resize_at(self, x: float) -> None:
        if not self._resize_active or self._resize_edge is None:
            return
        rid = self._resize_rid
        length = max(1, int(self._lengths.get(rid, self._max_length)))
        bp = self._x_to_feature_bp(x, length=length)
        start, end = self._resize_orig_start, self._resize_orig_end
        if self._resize_edge == "start":
            start = min(bp, end)
        else:
            end = max(bp, start)
        store = list(self._features.get(rid) or [])
        if not (0 <= self._resize_store_index < len(store)):
            return
        current = store[self._resize_store_index]
        if int(current.start) == start and int(current.end) == end:
            return
        updated = _feature_with_bounds(current, start, end)
        store[self._resize_store_index] = updated
        self._features[rid] = store
        # Keep packed tracks in sync for live paint without a full re-layout.
        flat_i = 0
        for track in self._tracks.get(rid) or []:
            for j, feat in enumerate(track):
                if flat_i == self._resize_flat:
                    track[j] = updated
                    self.update()
                    return
                flat_i += 1

    def _feature_full_label(self, feat: SequenceFeature) -> str:
        label = _feature_display_label(
            feat, show_database_name=self._show_similarity_database_name
        )
        extras: list[str] = []
        if float(getattr(feat, "pident", 0) or 0) > 0 and "(" not in label:
            extras.append(f"{float(feat.pident):.1f}% id")
        subject = str(getattr(feat, "subject_id", "") or "").strip()
        if subject and subject not in label:
            extras.append(subject)
        if extras:
            return f"{label} ({', '.join(extras)})" if label else ", ".join(extras)
        return label

    def _feature_label_overflows(self, feat: SequenceFeature, rect: QRectF) -> bool:
        label = self._feature_full_label(feat)
        if not label or self._feat_label_point_size <= 0:
            return False
        if rect.width() <= _FEAT_LABEL_MIN_WIDTH:
            return True
        label_font = QFont(self.font())
        label_font.setPointSize(self._feat_label_point_size)
        fm = QFontMetrics(label_font)
        available = max(0.0, rect.width() - 10.0)
        return fm.horizontalAdvance(label) > available

    def set_show_similarity_database_name(self, show: bool) -> None:
        show = bool(show)
        if show == self._show_similarity_database_name:
            return
        self._show_similarity_database_name = show
        self.update()

    def _features_intersecting_rect(self, rubber: QRect) -> set[tuple[int, int]]:
        """Features whose arrows intersect the rubber-band (partial overlap OK)."""
        hits: set[tuple[int, int]] = set()
        if rubber.isNull() or not self._show_annotation_track:
            return hits
        rx0, ry0 = float(rubber.left()), float(rubber.top())
        rx1, ry1 = float(rubber.right()), float(rubber.bottom())
        for row, rid in enumerate(self._record_ids):
            if row >= len(self._row_tops):
                continue
            top = self._row_tops[row]
            height = self._row_heights[row]
            if top + height < ry0 or top > ry1:
                continue
            length = max(1, int(self._lengths.get(rid, 1)))
            feat_y0 = self._feature_band_y0(row)
            flat_index = 0
            for track_i, track in enumerate(self._tracks.get(rid) or []):
                ty = feat_y0 + track_i * (self._feat_height + self._feat_gap)
                ty1 = ty + self._feat_height
                for feat in track:
                    start = max(1, min(int(feat.start), length))
                    end = max(1, min(int(feat.end), length))
                    if start > end:
                        start, end = end, start
                    x0 = self._bar_pad + (start - 1) * self._ppb
                    x1 = self._bar_pad + end * self._ppb
                    if x1 >= rx0 and x0 <= rx1 and ty1 >= ry0 and ty <= ry1:
                        hits.add((row, flat_index))
                    flat_index += 1
        return hits

    def selected_feature_column_span(self) -> tuple[int, int] | None:
        """Leftmost/rightmost 1-based columns of selected features."""
        if not self._selected_feats:
            return None
        left: int | None = None
        right: int | None = None
        for row, flat_i in self._selected_feats:
            feat = self._feature_by_flat_index(row, flat_i)
            if feat is None:
                continue
            a, b = int(feat.start), int(feat.end)
            if b < a:
                a, b = b, a
            left = a if left is None else min(left, a)
            right = b if right is None else max(right, b)
        if left is None or right is None:
            return None
        return left, right

    def _set_feature_selection(
        self,
        feats: set[tuple[int, int]],
        *,
        set_columns: bool = False,
    ) -> None:
        self._selected_feats = set(feats)
        self._reset_edge_expand_phase()
        if set_columns:
            span = self.selected_feature_column_span()
            if span is not None:
                self.set_region_selection(span[0], span[1])
            else:
                self.set_region_selection(None)
        self.update()

    def _select_row_click(self, row: int, modifiers: Qt.KeyboardModifier) -> None:
        """Select sequences from the name column (not from MSA residue clicks)."""
        if row < 0:
            return
        had_band = (
            self._region_sel is not None or self._regional_rows is not None
        )
        if had_band:
            self._region_sel = None
            self._regional_rows = None
            self._reset_edge_expand_phase()
            self.regionSelectionChanged.emit(None)
        if modifiers & (
            Qt.KeyboardModifier.ControlModifier | Qt.KeyboardModifier.MetaModifier
        ):
            if row in self._selected_rows:
                self._selected_rows.discard(row)
            else:
                self._selected_rows.add(row)
        elif modifiers & Qt.KeyboardModifier.ShiftModifier and self._selected_rows:
            anchor = min(self._selected_rows)
            lo, hi = (anchor, row) if anchor <= row else (row, anchor)
            self._selected_rows = set(range(lo, hi + 1))
        elif row in self._selected_rows:
            self._selected_rows = set()
        else:
            self._selected_rows = {row}
        self._sync_selection_chrome()
        self.selectionChanged.emit()
        self.rowsChanged.emit()

    def _columns_from_rubber(self, rubber: QRect) -> tuple[int, int]:
        start_bp = self.x_to_bp(rubber.left())
        end_bp = self.x_to_bp(rubber.right())
        if end_bp < start_bp:
            start_bp, end_bp = end_bp, start_bp
        start_i = max(1, int(math.floor(start_bp)) + 1)
        end_i = max(start_i, int(math.ceil(end_bp)))
        return start_i, end_i

    def apply_shift_edge_expand(self, direction: str) -> bool:
        """Shift+Left / Shift+Right from a feature or column selection.

        Uses the feature span when features are selected; otherwise the current
        column / regional column range. Left/right cycle:

          Left:  ``[1, L-1]`` then ``[1, R]``
          Right: ``[R+1, end]`` then ``[L, end]``
        """
        span = self.selected_feature_column_span()
        if span is None and self._region_sel is not None:
            # Stick to the original edges across the two-step cycle.
            if self._edge_expand_anchor is None:
                self._edge_expand_anchor = self._region_sel
            span = self._edge_expand_anchor
        if span is None:
            return False
        left, right = int(span[0]), int(span[1])
        if right < left:
            left, right = right, left
        maximum = max(1, int(self._max_length))
        key = "left" if direction == "left" else "right"
        phase = int(self._edge_expand_phase.get(key, 0)) % 2
        if key == "left":
            if phase == 0:
                if left <= 1:
                    return False
                self.set_region_selection(1, left - 1, keep_expand_anchor=True)
            else:
                self.set_region_selection(1, right, keep_expand_anchor=True)
        else:
            if phase == 0:
                if right >= maximum:
                    return False
                self.set_region_selection(right + 1, maximum, keep_expand_anchor=True)
            else:
                self.set_region_selection(left, maximum, keep_expand_anchor=True)
        self._edge_expand_phase[key] = 1 - phase
        other = "right" if key == "left" else "left"
        self._edge_expand_phase[other] = 0
        return True

    def mousePressEvent(self, event) -> None:  # noqa: N802
        if event.button() == Qt.MouseButton.LeftButton:
            pos = event.position()
            self._drag_origin = QPointF(pos)
            self._rubber = QRect()
            self._rubber_active = False
            self._resize_active = False
            self._press_edge = None
            # Drop any live column/row preview before a regional drag.
            if self._preview_columns is not None or self._preview_rows is not None:
                self.finish_live_selection_preview()
            row = self._row_at_y(pos.y())
            edge = self._edge_hit(row, pos.x(), pos.y()) if row >= 0 else None
            if edge is not None:
                flat_i, edge_kind = edge
                self._press_feat = (row, flat_i)
                self._press_edge = edge_kind
                self._drag_kind = "resize"
            else:
                feat_i = self._hit_feature(row, pos.x(), pos.y()) if row >= 0 else None
                self._press_feat = (row, feat_i) if feat_i is not None else None
                if row >= 0 and self._y_in_annotation_band(row, pos.y()):
                    self._drag_kind = "features"
                else:
                    self._drag_kind = "region"
                    # Clear MSA chrome immediately; scale/Position sync on commit.
                    self._dismiss_committed_selection_visuals()
            self.grabMouse()
            # Sequence selection is via the name column; MSA press only starts a drag.
            event.accept()
            return
        super().mousePressEvent(event)

    def contextMenuEvent(self, event: QContextMenuEvent) -> None:  # noqa: N802
        self.contextMenuRequested.emit(event.globalPos())
        event.accept()

    def mouseMoveEvent(self, event) -> None:  # noqa: N802
        if self._drag_origin is not None and event.buttons() & Qt.MouseButton.LeftButton:
            if self._hover_feat is not None:
                self._hover_feat = None
            pos = event.position()
            if self._drag_kind == "resize" and self._press_feat is not None:
                origin = self._drag_origin
                delta = math.hypot(pos.x() - origin.x(), pos.y() - origin.y())
                if not self._resize_active and delta >= _RUBBER_MIN_PX:
                    row, flat_i = self._press_feat
                    resolved = self._store_index_for_flat(row, flat_i)
                    if resolved is None:
                        event.accept()
                        return
                    rid, store_i, feat = resolved
                    self._resize_active = True
                    self._resize_rid = rid
                    self._resize_store_index = store_i
                    self._resize_flat = flat_i
                    self._resize_row = row
                    self._resize_edge = self._press_edge
                    self._resize_orig_start = int(feat.start)
                    self._resize_orig_end = int(feat.end)
                    self._set_feature_selection({(row, flat_i)}, set_columns=False)
                    self.setCursor(Qt.CursorShape.SizeHorCursor)
                if self._resize_active:
                    self._apply_resize_at(pos.x())
                event.accept()
                return
            origin = self._drag_origin.toPoint()
            # Coalesce backlog: always track the latest cursor, not a stale event.
            current = self.mapFromGlobal(QCursor.pos())
            self._ensure_drag_point_visible(current)
            rubber = QRect(origin, current).normalized()
            if rubber.width() < _RUBBER_MIN_PX and rubber.height() < _RUBBER_MIN_PX:
                event.accept()
                return
            if not self._rubber_active:
                self._preview_columns = None
                self._preview_rows = None
            if self._drag_kind == "region":
                # Pixel-follow outline; Position + names throttle via drag preview.
                self._set_live_region_rubber(rubber)
            else:
                if not self._rubber_active:
                    self._dismiss_committed_selection_visuals()
                self._update_rubber_overlay(rubber)
            event.accept()
            return
        pos = event.position()
        row = self._row_at_y(pos.y())
        edge = self._edge_hit(row, pos.x(), pos.y()) if row >= 0 else None
        feat_i = (
            edge[0]
            if edge is not None
            else (self._hit_feature(row, pos.x(), pos.y()) if row >= 0 else None)
        )
        hover = (row, feat_i) if feat_i is not None else None
        if hover != self._hover_feat:
            self._hover_feat = hover
            self.update()
        if edge is not None:
            self.setCursor(Qt.CursorShape.SizeHorCursor)
        else:
            self.unsetCursor()
        super().mouseMoveEvent(event)

    def _flush_region_drag_preview(self) -> None:
        pending = self._drag_preview_pending
        if pending is None:
            return
        start_i, end_i, rows = pending
        self._drag_preview_pending = None
        row_set = {int(r) for r in rows}
        if (
            self._drag_preview_columns == (start_i, end_i)
            and self._drag_preview_rows == row_set
        ):
            return
        self._drag_preview_columns = (start_i, end_i)
        self._drag_preview_rows = row_set
        self.regionDragPreview.emit((start_i, end_i, rows))

    def _clear_region_drag_preview(self) -> None:
        self._drag_preview_timer.stop()
        self._drag_preview_pending = None
        self._drag_preview_rows = None
        self._drag_preview_columns = None
        self._temp_matrix_sel = None

    def leaveEvent(self, event) -> None:  # noqa: N802
        if self._hover_feat is not None:
            self._hover_feat = None
            self.unsetCursor()
            self.update()
        super().leaveEvent(event)

    def mouseReleaseEvent(self, event) -> None:  # noqa: N802
        if event.button() == Qt.MouseButton.LeftButton and self._drag_origin is not None:
            if QWidget.mouseGrabber() is self:
                self.releaseMouse()
            rubber = QRect(self._rubber)
            active = self._rubber_active
            temp_sel = self._temp_matrix_sel
            press_feat = self._press_feat
            drag_kind = self._drag_kind
            resize_active = self._resize_active
            resize_rid = self._resize_rid
            resize_store = self._resize_store_index
            resize_orig_start = self._resize_orig_start
            resize_orig_end = self._resize_orig_end
            self._drag_origin = None
            self._press_feat = None
            self._press_edge = None
            self._drag_kind = None
            self._rubber = QRect()
            self._rubber_active = False
            self._resize_active = False
            self._clear_region_drag_preview()
            self._chrome.clear()
            self._clear_drag_paint_cache()
            self.unsetCursor()
            if resize_active and resize_store >= 0 and resize_rid:
                store = list(self._features.get(resize_rid) or [])
                if 0 <= resize_store < len(store):
                    feat = store[resize_store]
                    if (
                        int(feat.start) != resize_orig_start
                        or int(feat.end) != resize_orig_end
                    ):
                        self.featureResized.emit(resize_rid, resize_store, feat)
                event.accept()
                return
            if active and (
                temp_sel is not None
                or rubber.width() >= _RUBBER_MIN_PX
                or rubber.height() >= _RUBBER_MIN_PX
            ):
                feats: set[tuple[int, int]] = set()
                if drag_kind == "features":
                    feats = self._features_intersecting_rect(rubber)
                if drag_kind == "features" and feats:
                    # Partial overlap selects the feature → column span.
                    self._set_feature_selection(feats, set_columns=True)
                else:
                    # Sequence/regional select — including when the drag started
                    # on the annotation band but hit no features.
                    if temp_sel is not None:
                        start_i, end_i, r0, r1 = temp_sel
                        rows = list(range(r0, r1 + 1))
                    else:
                        rows = self._rows_in_y_range(rubber.top(), rubber.bottom())
                        start_i, end_i = self._columns_from_rubber(rubber)
                    self._selected_feats = set()
                    self.set_regional_selection(start_i, end_i, rows)
                    # Position spins once at the end (never during drag).
                    self.regionDragPreview.emit((start_i, end_i, rows))
                    self.setFocus(Qt.FocusReason.MouseFocusReason)
            elif press_feat is not None:
                # Click a feature → select it and its columns (incl. point features).
                mods = event.modifiers()
                if mods & (
                    Qt.KeyboardModifier.ControlModifier
                    | Qt.KeyboardModifier.MetaModifier
                ):
                    feats = set(self._selected_feats)
                    if press_feat in feats:
                        feats.discard(press_feat)
                    else:
                        feats.add(press_feat)
                    self._set_feature_selection(feats, set_columns=True)
                else:
                    self._set_feature_selection({press_feat}, set_columns=True)
                self.setFocus(Qt.FocusReason.MouseFocusReason)
            else:
                # Plain click on residues/bars clears column, feature, and
                # name-selected sequence highlights (selected or not).
                self._selected_feats = set()
                self._reset_edge_expand_phase()
                rows_cleared = bool(self._selected_rows)
                if rows_cleared:
                    self._selected_rows = set()
                    self._chrome.clear()
                if self._region_sel is not None:
                    self.set_region_selection(None)
                if rows_cleared:
                    self.selectionChanged.emit()
                    self.rowsChanged.emit()
            event.accept()
            return
        super().mouseReleaseEvent(event)

    def keyPressEvent(self, event: QKeyEvent) -> None:  # noqa: N802
        from PySide6.QtGui import QKeySequence

        if event.key() == Qt.Key.Key_Escape:
            cleared = False
            if self._region_sel is not None:
                self.set_region_selection(None)
                cleared = True
            if self._selected_rows:
                self._selected_rows = set()
                self._chrome.clear()
                self.selectionChanged.emit()
                self.rowsChanged.emit()
                cleared = True
            if cleared:
                event.accept()
                return
        mods = event.modifiers()
        if mods & Qt.KeyboardModifier.ShiftModifier:
            if event.key() == Qt.Key.Key_Left and self.apply_shift_edge_expand("left"):
                event.accept()
                return
            if event.key() == Qt.Key.Key_Right and self.apply_shift_edge_expand("right"):
                event.accept()
                return
        # Plain Delete/Backspace; on macOS ⌘⌫ also matches StandardKey.Delete.
        if event.matches(QKeySequence.StandardKey.Delete) or event.key() in (
            Qt.Key.Key_Delete,
            Qt.Key.Key_Backspace,
        ):
            self.deleteRequested.emit()
            event.accept()
            return
        if event.matches(QKeySequence.StandardKey.Undo) or (
            event.key() == Qt.Key.Key_Z
            and event.modifiers()
            & (Qt.KeyboardModifier.ControlModifier | Qt.KeyboardModifier.MetaModifier)
            and not (event.modifiers() & Qt.KeyboardModifier.ShiftModifier)
        ):
            # Forward to panel undo (MSA / feature Restore).
            node = self.parent()
            while node is not None:
                restore = getattr(node, "_restore_features", None)
                if callable(restore):
                    restore(True)
                    event.accept()
                    return
                node = node.parent()
        if event.matches(QKeySequence.StandardKey.SelectAll) or (
            event.key() == Qt.Key.Key_A
            and event.modifiers()
            & (Qt.KeyboardModifier.ControlModifier | Qt.KeyboardModifier.MetaModifier)
        ):
            self.select_all_rows()
            event.accept()
            return
        super().keyPressEvent(event)

    def mouseDoubleClickEvent(self, event) -> None:  # noqa: N802
        if event.button() == Qt.MouseButton.LeftButton:
            pos = event.position()
            row = self._row_at_y(pos.y())
            if row >= 0:
                feat_i = self._hit_feature(row, pos.x(), pos.y())
                if feat_i is not None:
                    self.featureActivated.emit(row, feat_i)
                    event.accept()
                    return
        super().mouseDoubleClickEvent(event)

    def wheelEvent(self, event: QWheelEvent) -> None:  # noqa: N802
        mods = event.modifiers() | QGuiApplication.keyboardModifiers()
        if mods & Qt.KeyboardModifier.ShiftModifier:
            delta = event.angleDelta().y()
            if delta == 0:
                delta = event.angleDelta().x()
            pixel = event.pixelDelta()
            if delta == 0 and not pixel.isNull():
                delta = int(pixel.y()) or int(pixel.x())
            if delta != 0:
                factor = _ZOOM_STEP if delta > 0 else 1 / _ZOOM_STEP
                self.zoom_by(factor, anchor_x=event.position().x())
            event.accept()
            return
        super().wheelEvent(event)

    def paintEvent(self, event) -> None:  # noqa: N802
        painter = QPainter(self)
        base = (
            _msa_canvas_color(self)
            if self._msa_mode
            else self.palette().color(QPalette.ColorRole.Base)
        )
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, not self._msa_mode)
        # Paint each dirty rect — not event.rect() (bbox of all strips). Column /
        # row edge dirties must not recompose the whole selection span.
        region = event.region()
        rects = list(region) if not region.isEmpty() else [event.rect()]
        use_cache = self._msa_mode and self._ensure_msa_layer_cache()
        cache = self._msa_layer_cache if use_cache else None
        cache_origin = self._msa_layer_cache_origin if use_cache else QPoint()
        cache_rect = self._msa_layer_cache_rect if use_cache else QRect()
        for clip in rects:
            if clip.isNull() or clip.width() <= 0 or clip.height() <= 0:
                continue
            if use_cache and cache is not None and not cache_rect.isEmpty():
                # Always clear the dirty rect first so a stale/small cache cannot
                # leave uncleared white slabs after a tab switch.
                painter.fillRect(clip, base)
                # IMPORTANT: do not pass a source QRect into drawPixmap on a
                # HiDPI QPixmap — Qt/PySide treats that source in *device*
                # pixels, so only the top-left half of the MSA appeared (same
                # bug as SequenceMapCanvas letter cache). Draw the full logical
                # overscrolled pixmap at its content origin and clip to dirty.
                painter.save()
                painter.setClipRect(clip)
                painter.drawPixmap(cache_origin, cache)
                painter.restore()
                # Outside the overscrolled cache (rare): fall back to live paint.
                outside = QRegion(clip).subtracted(QRegion(cache_rect))
                if not outside.isEmpty():
                    for out in outside:
                        if out.isNull() or out.width() <= 0 or out.height() <= 0:
                            continue
                        painter.fillRect(out, base)
                        self._paint_msa_rows(painter, out)
                if self._selected_feats:
                    self._paint_feature_selection_outlines(painter, clip)
            else:
                # Neighbor columns so sub-pixel MSA cells cover strip edges.
                paint_clip = clip.adjusted(-2, 0, 2, 0).intersected(self.rect())
                painter.fillRect(clip, base)
                if self._msa_mode:
                    self._paint_msa_rows(
                        painter, paint_clip if not paint_clip.isNull() else clip
                    )
                else:
                    self._paint_bar_rows(painter, clip)
            self._chrome.paint(painter, clip)
        if not self._rubber_active:
            self._paint_feature_hover_tip(painter)
        painter.end()

    def _paint_feature_selection_outlines(
        self, painter: QPainter, clip: QRect
    ) -> None:
        """Selection borders only — kept out of the residue cache."""
        if not self._selected_feats:
            return
        pen = QPen(QColor("#4ea1ff" if _is_dark(self) else "#0b5cab"), 2)
        painter.setBrush(Qt.BrushStyle.NoBrush)
        painter.setPen(pen)
        for row, flat_i in self._selected_feats:
            rect = self._feature_rect(row, flat_i)
            if rect is None or not rect.intersects(QRectF(clip)):
                continue
            painter.drawRect(rect.adjusted(-1, -1, 1, 1))

    def _paint_feature_tracks(
        self,
        painter: QPainter,
        *,
        row: int,
        rid: str,
        length: int,
        feat_y0: int,
        paint_selection: bool = True,
        clip: QRect | None = None,
    ) -> None:
        if not self._show_annotation_track:
            return
        clip_left = float(clip.left()) - 2.0 if clip is not None else None
        clip_right = float(clip.right()) + 2.0 if clip is not None else None
        flat_index = 0
        for track_i, track in enumerate(self._tracks.get(rid) or []):
            ty = feat_y0 + track_i * (self._feat_height + self._feat_gap)
            for feat in track:
                start = max(1, min(int(feat.start), length))
                end = max(1, min(int(feat.end), length))
                if start > end:
                    start, end = end, start
                x0 = self._bar_pad + (start - 1) * self._ppb
                x1 = self._bar_pad + end * self._ppb
                if clip_left is not None and clip_right is not None:
                    if x1 < clip_left or x0 > clip_right:
                        flat_index += 1
                        continue
                rect = QRectF(x0, ty, max(2.0, x1 - x0), float(self._feat_height))
                try:
                    fc = QColor(feat.color) if feat.color else _DEFAULT_FEAT
                    if not fc.isValid():
                        fc = _DEFAULT_FEAT
                except Exception:
                    fc = _DEFAULT_FEAT
                label = ""
                if (
                    self._feat_label_point_size > 0
                    and rect.width() > _FEAT_LABEL_MIN_WIDTH
                ):
                    label = _feature_display_label(
                        feat,
                        show_database_name=self._show_similarity_database_name,
                    )
                _draw_arrow_feature(
                    painter,
                    rect,
                    fc,
                    forward=feat.is_forward(),
                    label=label,
                    label_point_size=self._feat_label_point_size,
                )
                if paint_selection and (row, flat_index) in self._selected_feats:
                    painter.setBrush(Qt.BrushStyle.NoBrush)
                    painter.setPen(
                        QPen(QColor("#4ea1ff" if _is_dark(self) else "#0b5cab"), 2)
                    )
                    painter.drawRect(rect.adjusted(-1, -1, 1, 1))
                flat_index += 1

    def _paint_feature_hover_tip(self, painter: QPainter) -> None:
        hover = self._hover_feat
        if hover is None or self._rubber_active:
            return
        row, flat_index = hover
        feat = self._feature_by_flat_index(row, flat_index)
        rect = self._feature_rect(row, flat_index)
        if feat is None or rect is None:
            return
        full_label = self._feature_full_label(feat)
        if not full_label or not self._feature_label_overflows(feat, rect):
            return
        tip_font = QFont(painter.font())
        tip_font.setPointSize(8)
        painter.setFont(tip_font)
        fm = QFontMetrics(tip_font)
        tip_w = min(
            max(fm.horizontalAdvance(full_label) + 10, 40),
            max(80, self.width() - 8),
        )
        tip_h = fm.height() + 6
        tip_x = min(
            max(4.0, rect.left()),
            max(4.0, float(self.width()) - tip_w - 4.0),
        )
        tip_y = rect.bottom() + 3.0
        if tip_y + tip_h > self.height() - 2:
            tip_y = max(2.0, rect.top() - tip_h - 3.0)
        tip_rect = QRectF(tip_x, tip_y, tip_w, tip_h)
        painter.setPen(QPen(QColor("#555555"), 1))
        painter.setBrush(QColor(255, 255, 255, 235))
        painter.drawRoundedRect(tip_rect, 3, 3)
        painter.setPen(QColor("#222222"))
        painter.drawText(
            tip_rect.adjusted(5, 0, -5, 0),
            Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft,
            full_label,
        )

    def _paint_msa_rows_pixel_sample(
        self,
        painter: QPainter,
        clip: QRect,
        *,
        residue_h: int,
        paint_feature_selection: bool,
    ) -> None:
        """Fit/overview: one colour sample per screen pixel (not per column).

        Matches the single-sequence map's ``_draw_sequence_bar_viewport`` path
        and AliView's zoomed-out rendering: cost is O(rows × viewport_width),
        not O(rows × alignment_length). Critical for ~15k-column MSAs at Fit.
        """
        ppb = max(self._ppb, 1e-12)
        max_len = max(1, int(self._max_length))
        px_left = max(0, int(math.floor(clip.left())))
        px_right = min(self.width(), int(math.ceil(clip.right())) + 1)
        if px_left >= px_right:
            return
        painter.setPen(Qt.PenStyle.NoPen)
        gap_line = QColor("#e0e0e0") if _is_dark(self) else QColor("#000000")
        for index in self._rows_in_y_range(float(clip.top()), float(clip.bottom()) + 1.0):
            rid = self._record_ids[index]
            top = self._row_tops[index]
            seq = self._msa_seqs[index] if index < len(self._msa_seqs) else ""
            y = float(top)
            y_mid = y + float(residue_h) * 0.5
            px = px_left
            while px < px_right:
                content_x = float(px) - float(self._bar_pad)
                if content_x < 0.0:
                    px += 1
                    continue
                col = int(content_x / ppb)
                if col >= max_len:
                    break
                ch = seq[col] if col < len(seq) else "-"
                fill = self._msa_overview_fill_color(col, ch, row=index)
                run_start = px
                px += 1
                while px < px_right:
                    content_x = float(px) - float(self._bar_pad)
                    if content_x < 0.0:
                        break
                    col2 = int(content_x / ppb)
                    if col2 >= max_len:
                        break
                    ch2 = seq[col2] if col2 < len(seq) else "-"
                    fill2 = self._msa_overview_fill_color(col2, ch2, row=index)
                    if fill is None or fill2 is None:
                        if fill is not fill2:
                            break
                        if (ch in "-.") != (ch2 in "-."):
                            break
                    elif fill.rgba() != fill2.rgba():
                        break
                    px += 1
                run_w = float(px - run_start)
                if fill is not None:
                    painter.fillRect(QRectF(float(run_start), y, run_w, float(residue_h)), fill)
                elif ch in "-.":
                    pen = QPen(gap_line, 1)
                    pen.setCosmetic(True)
                    painter.setPen(pen)
                    painter.drawLine(
                        QPointF(float(run_start), y_mid),
                        QPointF(float(px), y_mid),
                    )
                    painter.setPen(Qt.PenStyle.NoPen)
            length = max(1, int(self._lengths.get(rid, 1)))
            self._paint_feature_tracks(
                painter,
                row=index,
                rid=rid,
                length=length,
                feat_y0=self._feature_band_y0(index),
                paint_selection=paint_feature_selection,
                clip=clip,
            )

    def _paint_msa_rows(
        self,
        painter: QPainter,
        clip: QRect,
        *,
        paint_feature_selection: bool = True,
    ) -> None:
        """Column-aligned residue cells (visible rows/columns only)."""
        ppb = max(self._ppb, 1e-12)
        show_letters = self._msa_shows_letters()
        # Letter zoom: square band up to size-10; further zoom widens ppb only.
        residue_h = max(1, int(self._msa_letter_band_height()))
        if not show_letters and ppb < _MSA_PIXEL_SAMPLE_PPB:
            self._paint_msa_rows_pixel_sample(
                painter,
                clip,
                residue_h=residue_h,
                paint_feature_selection=paint_feature_selection,
            )
            return
        c0 = max(0, int((clip.left() - self._bar_pad) / ppb) - 1)
        c1 = min(self._max_length, int((clip.right() - self._bar_pad) / ppb) + 2)
        qcolors = self._msa_qcolors
        letter_side = float(_msa_letter_size_metrics(_MSA_LETTER_SIZE_PT)[0])
        residue_font = (
            self._msa_letter_font_for_cell(ppb, float(residue_h))
            if show_letters
            else None
        )
        # Overview: Column distance insets. Letter mode: square→wide fill.
        _bg_margin = max(0.0, float(self._column_distance) * 0.5)
        painter.setPen(Qt.PenStyle.NoPen)
        for index in self._rows_in_y_range(float(clip.top()), float(clip.bottom()) + 1.0):
            rid = self._record_ids[index]
            top = self._row_tops[index]
            # No selection wash on residues — dashed outline / name highlight only.
            seq = self._msa_seqs[index] if index < len(self._msa_seqs) else ""
            y = top
            # Feature tracks change the painter font; restore the residue font
            # for every row so later sequences do not inherit annotation text size.
            if residue_font is not None:
                painter.setFont(residue_font)
            # Coalesce same-residue runs when zoomed out (far fewer draw calls).
            j = c0
            while j < c1:
                x = self._bar_pad + j * ppb
                if x + ppb < clip.left():
                    j += 1
                    continue
                if x > clip.right():
                    break
                ch = seq[j] if j < len(seq) else "-"
                fill_bg = self._msa_should_fill_background(j, ch, row=index)
                paint_pair = self._msa_residue_paint_pair(j, ch)
                overview_fill = (
                    None
                    if show_letters
                    else self._msa_overview_fill_color(j, ch, row=index)
                )
                run_end = j + 1
                # Coalesce only when columns share one fill/colour (distance 0).
                if not show_letters and self._column_distance <= 0:
                    while run_end < c1:
                        nxt = seq[run_end] if run_end < len(seq) else "-"
                        nxt_fill = self._msa_overview_fill_color(
                            run_end, nxt, row=index
                        )
                        if overview_fill is None or nxt_fill is None:
                            if overview_fill is not nxt_fill:
                                break
                        elif overview_fill.rgba() != nxt_fill.rgba():
                            break
                        run_end += 1
                span = run_end - j
                # Cover [x, next_x). When zoomed out, span at least 1 logical px so
                # partial dirty strips (selection edges) cannot leave white hairlines.
                x_next = self._bar_pad + run_end * ppb
                width = max(x_next - x, 1.0 if not show_letters else max(0.5, ppb * span))
                cell = QRectF(x, y, width, float(residue_h))
                if ch in "-.":
                    # Gaps stay on the white canvas unless Color vs standard
                    # marks them as differing from the standard.
                    if show_letters:
                        letter_cell = QRectF(
                            x, y, max(0.5, ppb), float(residue_h)
                        )
                        if fill_bg:
                            gap_fill = _msa_gap_diff_color(self)
                            bg = _msa_letter_bg_rect(
                                x,
                                y,
                                max(0.5, ppb),
                                float(residue_h),
                                letter_side=letter_side,
                            )
                            painter.fillRect(bg, gap_fill)
                            from repeatviewer.widgets.residue_color_grid import (
                                ResidueColorGrid,
                            )

                            painter.setPen(
                                QColor(
                                    ResidueColorGrid.contrast_text_color(gap_fill)
                                )
                            )
                        else:
                            gap_fg = (
                                QColor("#d0d0d0")
                                if _is_dark(self)
                                else QColor("#666666")
                            )
                            painter.setPen(gap_fg)
                        painter.drawText(
                            letter_cell.toRect(),
                            int(Qt.AlignmentFlag.AlignCenter),
                            "-",
                        )
                        painter.setPen(Qt.PenStyle.NoPen)
                    elif overview_fill is not None:
                        m = _bg_margin
                        if m > 0 and (
                            cell.width() > 2.0 * m + 0.5
                            and cell.height() > 2.0 * m + 0.5
                        ):
                            painter.fillRect(
                                cell.adjusted(m, m, -m, -m), overview_fill
                            )
                        else:
                            painter.fillRect(cell, overview_fill)
                    else:
                        # Overview: thin backbone through gaps between residues.
                        y_mid = float(y) + float(residue_h) * 0.5
                        gap_line = (
                            QColor("#e0e0e0")
                            if _is_dark(self)
                            else QColor("#000000")
                        )
                        pen = QPen(gap_line, 1)
                        pen.setCosmetic(True)
                        painter.setPen(pen)
                        painter.drawLine(
                            QPointF(x, y_mid), QPointF(x + width, y_mid)
                        )
                        painter.setPen(Qt.PenStyle.NoPen)
                elif not show_letters:
                    # Overview: scheme colour or multi-seq light-grey backbone.
                    if overview_fill is not None:
                        m = _bg_margin
                        if m > 0 and (
                            cell.width() > 2.0 * m + 0.5
                            and cell.height() > 2.0 * m + 0.5
                        ):
                            painter.fillRect(
                                cell.adjusted(m, m, -m, -m), overview_fill
                            )
                        else:
                            painter.fillRect(cell, overview_fill)
                else:
                    pair = paint_pair if paint_pair is not None else qcolors.get(ch)
                    if pair is not None:
                        # One column per letter; square at ¾, wider toward max.
                        letter_cell = QRectF(
                            x, y, max(0.5, ppb), float(residue_h)
                        )
                        if fill_bg:
                            bg = _msa_letter_bg_rect(
                                x,
                                y,
                                max(0.5, ppb),
                                float(residue_h),
                                letter_side=letter_side,
                            )
                            painter.fillRect(bg, pair[0])
                            # Black or white from background luminance for readability.
                            from repeatviewer.widgets.residue_color_grid import (
                                ResidueColorGrid,
                            )

                            painter.setPen(
                                QColor(
                                    ResidueColorGrid.contrast_text_color(pair[0])
                                )
                            )
                        else:
                            # No fill → residue colour as the letter.
                            painter.setPen(pair[0])
                        painter.drawText(
                            letter_cell.toRect(),
                            int(Qt.AlignmentFlag.AlignCenter),
                            ch,
                        )
                        painter.setPen(Qt.PenStyle.NoPen)
                j = run_end
            hits = self._find_hits.get(rid)
            if hits:
                painter.setPen(QPen(_FIND_HIT_BORDER, 1))
                current = self._find_current
                for hit_start, hit_end in hits:
                    x0 = self._bar_pad + (hit_start - 1) * ppb
                    x1 = self._bar_pad + hit_end * ppb
                    mark = QRectF(x0, y, max(2.0, x1 - x0), float(residue_h))
                    is_current = (
                        current is not None and current == (rid, hit_start, hit_end)
                    )
                    painter.fillRect(
                        mark, _FIND_CURRENT_FILL if is_current else _FIND_HIT_FILL
                    )
                    painter.drawRect(mark)
                painter.setPen(Qt.PenStyle.NoPen)
            length = max(1, int(self._lengths.get(rid, 1)))
            self._paint_feature_tracks(
                painter,
                row=index,
                rid=rid,
                length=length,
                feat_y0=self._feature_band_y0(index),
                paint_selection=paint_feature_selection,
                clip=clip,
            )
        # Selection chrome lives entirely on ``_chrome`` (see _sync_selection_chrome).

    def _paint_bar_rows(self, painter: QPainter, clip: QRect) -> None:
        show_letters = self._bar_shows_letters()
        ppb = max(self._ppb, 1e-12)
        residue_h = max(1, int(self._row_letter_band_height()))
        letter_side = float(_msa_letter_size_metrics(_MSA_LETTER_SIZE_PT)[0])
        residue_font = (
            self._msa_letter_font_for_cell(ppb, float(residue_h))
            if show_letters
            else None
        )

        for index in self._rows_in_y_range(float(clip.top()), float(clip.bottom()) + 1.0):
            rid = self._record_ids[index]
            top = self._row_tops[index]
            # No selection wash on bars — dashed outline / name highlight only.
            length = max(1, int(self._lengths.get(rid, 1)))
            bar_y = top if show_letters else self._sequence_bar_y(index)
            if show_letters:
                # Feature tracks change the painter font; restore it per row.
                if residue_font is not None:
                    painter.setFont(residue_font)
                seq = self._msa_seqs[index] if index < len(self._msa_seqs) else ""
                c0 = max(0, int((clip.left() - self._bar_pad) / ppb) - 1)
                c1 = min(length, int((clip.right() - self._bar_pad) / ppb) + 2)
                for column in range(c0, c1):
                    ch = seq[column] if column < len(seq) else ""
                    if not ch:
                        continue
                    cell_w = max(0.5, ppb)
                    cell = QRectF(
                        self._bar_pad + column * ppb,
                        bar_y,
                        cell_w,
                        float(residue_h),
                    )
                    pair = self._msa_qcolors.get(ch)
                    if pair is not None and self._color_unaligned_background:
                        painter.fillRect(
                            _msa_letter_bg_rect(
                                cell.x(),
                                cell.y(),
                                cell_w,
                                float(residue_h),
                                letter_side=letter_side,
                            ),
                            pair[0],
                        )
                    if pair is None:
                        letter_color = self.palette().text().color()
                    elif self._color_unaligned_background:
                        letter_color = pair[1]
                    else:
                        letter_color = pair[0]
                    painter.setPen(letter_color)
                    painter.drawText(cell.toRect(), int(Qt.AlignmentFlag.AlignCenter), ch)
                painter.setPen(Qt.PenStyle.NoPen)
            else:
                bar_w = max(2.0, length * self._ppb)
                bar_h = float(self._multi_seq_bar_height(index))
                painter.fillRect(
                    QRectF(self._bar_pad, float(bar_y), bar_w, bar_h),
                    _bar_color(self),
                )
            hits = self._find_hits.get(rid)
            if hits:
                current = self._find_current
                painter.setPen(QPen(_FIND_HIT_BORDER, 1))
                for hit_start, hit_end in hits:
                    x0 = self._bar_pad + (hit_start - 1) * self._ppb
                    x1 = self._bar_pad + hit_end * self._ppb
                    if show_letters:
                        mark_h = residue_h
                        mark_y = bar_y
                    else:
                        bar_h = self._multi_seq_bar_height(index)
                        mark_h = max(2, bar_h + 2)
                        mark_y = max(top, bar_y - 1)
                    mark = QRectF(x0, mark_y, max(3.0, x1 - x0), mark_h)
                    is_current = current is not None and current == (rid, hit_start, hit_end)
                    painter.fillRect(mark, _FIND_CURRENT_FILL if is_current else _FIND_HIT_FILL)
                    painter.drawRect(mark)
                painter.setPen(Qt.PenStyle.NoPen)

            self._paint_feature_tracks(
                painter,
                row=index,
                rid=rid,
                length=length,
                feat_y0=self._feature_band_y0(index),
            )
        # Selection chrome lives entirely on ``_chrome``.


class _ScaleRuler(QWidget):
    """Horizontal bp scale under the tracks only.

    Click selects one column; drag selects a column range; double-click zooms
    to maximum at that position. The selected column band is painted on the
    scale to match the MSA overlay.
    """

    columnRegionSelected = Signal(int, int)  # 1-based inclusive start, end
    columnRegionPreview = Signal(object)  # tuple[int,int] | None (legacy / tests)
    columnPixelPreview = Signal(float, float)  # scale-local x0,x1 while dragging
    # Content-x under cursor + zoom-in flag (True=in, False=out).
    zoomAtRequested = Signal(float, bool)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._ppb = 0.05
        self._max_length = 1
        self._bar_pad = 8
        self._scroll_x = 0
        self._font_point_size = 8
        self._region_sel: tuple[int, int] | None = None
        self._drag_origin_x: float | None = None
        self._drag_current_x: float | None = None
        self._rubber_active = False
        self._suppress_click_select = False
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.setCursor(Qt.CursorShape.SizeHorCursor)
        self.setToolTip(
            "Click to select a column · "
            "Double-click to zoom to maximum at this position · "
            "Drag to select columns · "
            "Right-click or Shift+click to fit"
        )
        # Focus so Shift+Left/Right (column edge expand) works after scale select.
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self._sync_height_to_font()

    def set_font_point_size(self, point_size: int) -> None:
        self._font_point_size = max(4, min(24, int(point_size)))
        self._sync_height_to_font()
        self.update()

    def font_point_size(self) -> int:
        return int(self._font_point_size)

    def _sync_height_to_font(self) -> None:
        font = QFont(self.font())
        font.setFamily("Menlo")
        font.setStyleHint(QFont.StyleHint.Monospace)
        font.setPointSize(max(4, min(24, int(self._font_point_size))))
        fm = QFontMetrics(font)
        # Label ascent + gap + tick + bottom pad.
        height = max(_SCALE_HEIGHT, fm.ascent() + 8)
        self.setFixedHeight(height)

    def set_scale(self, *, ppb: float, max_length: int, bar_pad: int) -> None:
        self._ppb = max(_ABS_MIN_PPB, float(ppb))
        self._max_length = max(1, int(max_length))
        self._bar_pad = int(bar_pad)
        self.update()

    def set_scroll_x(self, value: int) -> None:
        self._scroll_x = int(value)
        self.update()

    def set_region_selection(self, start: int | None, end: int | None = None) -> None:
        """Mirror the MSA column band onto the scale."""
        if start is None:
            if self._region_sel is None:
                return
            self._region_sel = None
            self.update()
            return
        maximum = max(1, int(self._max_length))
        a = max(1, min(maximum, int(start)))
        b = max(1, min(maximum, int(end if end is not None else start)))
        if b < a:
            a, b = b, a
        new = (a, b)
        if new == self._region_sel:
            return
        self._region_sel = new
        self.update()

    def _x_to_bp(self, x: float) -> float:
        return (float(x) + self._scroll_x - self._bar_pad) / max(self._ppb, 1e-12)

    def _bp_to_x(self, bp: float) -> float:
        return self._bar_pad + float(bp) * max(self._ppb, 1e-12) - self._scroll_x

    def _columns_from_x_range(self, x0: float, x1: float) -> tuple[int, int]:
        a = self._x_to_bp(x0)
        b = self._x_to_bp(x1)
        if b < a:
            a, b = b, a
        start_i = max(1, int(math.floor(a)) + 1)
        end_i = max(start_i, min(self._max_length, int(math.ceil(b))))
        return start_i, end_i

    def _emit_zoom_at(self, scale_x: float, *, zoom_in: bool) -> None:
        """Zoom centered on the content x under ``scale_x`` (viewport-local)."""
        content_x = float(scale_x) + float(self._scroll_x)
        self.zoomAtRequested.emit(content_x, bool(zoom_in))

    def mousePressEvent(self, event) -> None:  # noqa: N802
        if event.button() == Qt.MouseButton.LeftButton:
            self._suppress_click_select = False
            self._drag_origin_x = float(event.position().x())
            self._drag_current_x = self._drag_origin_x
            self._rubber_active = False
            self.grabMouse()
            event.accept()
            return
        if event.button() == Qt.MouseButton.RightButton:
            # Right-click → zoom out at this position (no column select).
            self._emit_zoom_at(float(event.position().x()), zoom_in=False)
            event.accept()
            return
        super().mousePressEvent(event)

    def mouseMoveEvent(self, event) -> None:  # noqa: N802
        if (
            self._drag_origin_x is not None
            and event.buttons() & Qt.MouseButton.LeftButton
        ):
            # Coalesce backlog: follow the latest cursor, not a stale event.
            self._drag_current_x = float(self.mapFromGlobal(QCursor.pos()).x())
            if abs(self._drag_current_x - self._drag_origin_x) >= _RUBBER_MIN_PX:
                if not self._rubber_active:
                    # Drop the prior scale highlight as soon as the new drag starts.
                    self._region_sel = None
                self._rubber_active = True
                # Pixel-follow MSA chrome; Position spins update via the panel.
                self.columnPixelPreview.emit(self._drag_origin_x, self._drag_current_x)
                # Cheap local scale guides only (widget is short).
                self.update()
            event.accept()
            return
        super().mouseMoveEvent(event)

    def mouseDoubleClickEvent(self, event) -> None:  # noqa: N802
        if event.button() == Qt.MouseButton.LeftButton:
            # Suppress the trailing release so it does not re-select a column.
            self._suppress_click_select = True
            self._drag_origin_x = None
            self._drag_current_x = None
            self._rubber_active = False
            if QWidget.mouseGrabber() is self:
                self.releaseMouse()
            self._emit_zoom_at(float(event.position().x()), zoom_in=True)
            event.accept()
            return
        super().mouseDoubleClickEvent(event)

    def mouseReleaseEvent(self, event) -> None:  # noqa: N802
        if event.button() == Qt.MouseButton.LeftButton and self._drag_origin_x is not None:
            if QWidget.mouseGrabber() is self:
                self.releaseMouse()
            x0 = self._drag_origin_x
            x1 = float(event.position().x())
            active = self._rubber_active
            mods = event.modifiers()
            self._drag_origin_x = None
            self._drag_current_x = None
            self._rubber_active = False
            if self._suppress_click_select:
                self._suppress_click_select = False
                self.update()
                event.accept()
                return
            if active and abs(x1 - x0) >= _RUBBER_MIN_PX:
                start_i, end_i = self._columns_from_x_range(x0, x1)
                self._region_sel = (start_i, end_i)
                self.columnRegionSelected.emit(start_i, end_i)
            elif mods & Qt.KeyboardModifier.ShiftModifier:
                # Shift+click → fit (zoom out) at this position.
                self._emit_zoom_at(x1, zoom_in=False)
            else:
                # Click (no drag) → select the column under the cursor.
                start_i, end_i = self._columns_from_x_range(x1, x1)
                self._region_sel = (start_i, end_i)
                self.columnRegionSelected.emit(start_i, end_i)
            self.setFocus(Qt.FocusReason.MouseFocusReason)
            self.update()
            event.accept()
            return
        if (
            event.button() == Qt.MouseButton.LeftButton
            and self._suppress_click_select
        ):
            self._suppress_click_select = False
            event.accept()
            return
        super().mouseReleaseEvent(event)

    def _bars_canvas(self):
        """Owning overview panel's MSA canvas, if present."""
        node = self.parent()
        while node is not None:
            bars = getattr(node, "_bars", None)
            if bars is not None:
                return bars
            node = node.parent()
        return None

    def keyPressEvent(self, event: QKeyEvent) -> None:  # noqa: N802
        from PySide6.QtGui import QKeySequence

        bars = self._bars_canvas()
        # Same Shift+Left / Shift+Right cycle as feature / MSA column selection.
        if event.modifiers() & Qt.KeyboardModifier.ShiftModifier and event.key() in (
            Qt.Key.Key_Left,
            Qt.Key.Key_Right,
        ):
            apply = getattr(bars, "apply_shift_edge_expand", None) if bars else None
            if callable(apply):
                direction = "left" if event.key() == Qt.Key.Key_Left else "right"
                if apply(direction):
                    event.accept()
                    return
        if event.matches(QKeySequence.StandardKey.Delete) or event.key() in (
            Qt.Key.Key_Delete,
            Qt.Key.Key_Backspace,
        ):
            delete = getattr(bars, "deleteRequested", None) if bars else None
            if delete is not None:
                delete.emit()
                event.accept()
                return
        if event.key() == Qt.Key.Key_Escape:
            clear = getattr(bars, "set_region_selection", None) if bars else None
            if callable(clear) and getattr(bars, "region_selection", lambda: None)():
                clear(None)
                event.accept()
                return
        super().keyPressEvent(event)

    def _paint_column_selection(self, painter: QPainter) -> None:
        """Tint + dashed column guides for the selected region on the scale."""
        if (
            self._rubber_active
            and self._drag_origin_x is not None
            and self._drag_current_x is not None
        ):
            letter_level = max(self._ppb, 1e-12) >= _MSA_LETTER_MIN_PPB
            if letter_level:
                # Nucleotide zoom: guides snap to column edges (match MSA chrome).
                start, end = self._columns_from_x_range(
                    self._drag_origin_x, self._drag_current_x
                )
                x0 = self._bp_to_x(start - 1)
                x1 = self._bp_to_x(end)
            else:
                x0 = float(min(self._drag_origin_x, self._drag_current_x))
                x1 = float(max(self._drag_origin_x, self._drag_current_x))
            edge = QColor("#4ea1ff") if _is_dark(self) else QColor("#0b5cab")
            painter.setPen(QPen(edge, 1.5, Qt.PenStyle.DashLine))
            painter.drawLine(int(x0), 0, int(x0), self.height())
            painter.drawLine(int(x1), 0, int(x1), self.height())
            painter.setPen(Qt.PenStyle.NoPen)
            return
        sel = self._region_sel
        if sel is None:
            return
        start, end = sel
        x0 = self._bp_to_x(start - 1)
        x1 = self._bp_to_x(end)
        if x1 < x0:
            x0, x1 = x1, x0
        # Match MSA chrome: keep a visible band for point (start==end) columns.
        min_w = max(_MIN_COLUMN_BAND_PX, float(self._ppb))
        if (x1 - x0) < min_w - 1e-9:
            mid = 0.5 * (x0 + x1)
            x0 = mid - 0.5 * min_w
            x1 = mid + 0.5 * min_w
        band = _REGION_BAND_DARK if _is_dark(self) else _REGION_BAND_LIGHT
        painter.fillRect(
            QRectF(x0, 0, max(1.0, x1 - x0), float(self.height())),
            band,
        )
        edge = QColor("#4ea1ff") if _is_dark(self) else QColor("#0b5cab")
        painter.setPen(QPen(edge, 1.5, Qt.PenStyle.DashLine))
        painter.drawLine(int(round(x0)), 0, int(round(x0)), self.height())
        painter.drawLine(int(round(x1)), 0, int(round(x1)), self.height())
        painter.setPen(Qt.PenStyle.NoPen)

    def paintEvent(self, _event) -> None:  # noqa: N802
        """Ruler: short ticks centered on nucleotides; numbers above the ticks."""
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.TextAntialiasing, True)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        base = self.palette().color(QPalette.ColorRole.Base)
        painter.fillRect(self.rect(), base)
        self._paint_column_selection(painter)
        fg = self.palette().color(QPalette.ColorRole.Text)
        if _is_dark(self):
            fg = QColor("#c8c8c8")
        font = QFont(self.font())
        font.setFamily("Menlo")
        font.setStyleHint(QFont.StyleHint.Monospace)
        font.setPointSize(max(4, min(24, int(self._font_point_size))))
        painter.setFont(font)
        fm = QFontMetrics(font)

        width = self.width()
        bars_left = 0.0
        bars_right = float(width)
        bars_width = max(1.0, float(width))
        origin = float(self._bar_pad)
        ppb = max(self._ppb, 1e-12)
        scroll = float(self._scroll_x)

        # Visible 1-based inclusive range (position 1 = first residue / letter).
        start_bp = max(1.0, (scroll - origin) / ppb + 1.0)
        end_bp = min(
            float(self._max_length),
            (scroll - origin + bars_width) / ppb + 1.0,
        )
        span = max(1.0, end_bp - start_bp)
        sample = f"{int(end_bp)}"
        min_label_px = fm.horizontalAdvance(sample) + 12
        step = _tick_step_for_span(span, bars_width, float(min_label_px))

        # Number (top) → short tick → nucleotide. Tick sits just above the track.
        tick_len = 3
        tick_bottom = self.height() - 2
        tick_top = tick_bottom - tick_len
        label_baseline = max(fm.ascent() + 1, tick_top - 2)

        def _center_x(pos: float) -> float:
            return origin + (pos - 0.5) * ppb - scroll

        label_positions: list[int] = []
        if start_bp <= 1.0 <= end_bp:
            label_positions.append(1)
        first = int(math.ceil(start_bp / step) * step)
        if first < step:
            first = step
        pos_i = first
        guard = 0
        while pos_i <= int(end_bp) + step and guard < 500:
            guard += 1
            if 1 <= pos_i <= self._max_length and pos_i != 1:
                label_positions.append(pos_i)
            pos_i += step

        # Short tick above every residue when zoomed in; otherwise at label steps.
        tick_start = max(1, int(math.floor(start_bp)))
        tick_end = min(self._max_length, int(math.ceil(end_bp)))
        show_every_base = ppb >= 3.5
        tick_positions = (
            range(tick_start, tick_end + 1)
            if show_every_base
            else label_positions
        )
        painter.setPen(QPen(fg, 1))
        for pos in tick_positions:
            cx = _center_x(float(pos))
            if cx < -2 or cx > width + 2:
                continue
            painter.drawLine(int(cx), tick_top, int(cx), tick_bottom)

        last_label_right = -1.0e9
        for pos in label_positions:
            cx = _center_x(float(pos))
            if cx < -2 or cx > width + 2:
                continue
            label = f"{int(pos)}"
            label_w = fm.horizontalAdvance(label)
            label_left = cx - label_w / 2
            if label_left < bars_left:
                label_left = bars_left
            label_right = label_left + label_w
            if label_right > bars_right:
                label_left = bars_right - label_w
                label_right = label_left + label_w
            if label_left >= last_label_right + 4 and label_right <= bars_right + 0.5:
                painter.drawText(int(label_left), int(label_baseline), label)
                last_label_right = label_right
        painter.end()


class _OverviewToolbarSlot(QWidget):
    """Map toolbar host that never forces a wide minimum and clips paint overflow.

    In split panes the −/+/Position/Select row can prefer a width wider than the
    map column. Without a zero minimum and a clip mask, controls spill onto the
    scale below and under the Display rail instead of being hidden by overflow.
    """

    def minimumSizeHint(self) -> QSize:  # noqa: N802
        return QSize(0, self.minimumHeight() or ANNOTATION_TOOLBAR_RENDERED_HEIGHT)

    def sizeHint(self) -> QSize:  # noqa: N802
        return QSize(
            max(0, self.width()),
            self.height() or ANNOTATION_TOOLBAR_RENDERED_HEIGHT,
        )

    def resizeEvent(self, event) -> None:  # noqa: N802
        super().resizeEvent(event)
        # Clip children to this rect so a half-shown button cannot paint over
        # the scale ruler or the Display rail beside the map.
        self.setMask(QRegion(self.rect()))


class MultiSequenceOverviewPanel(QWidget):
    """Annotation overview: names (left) | tracks + scale (right)."""

    statusMessage = Signal(str)
    annotationsSaved = Signal(object)
    sequenceExtracted = Signal(object)
    # Open cleaned / derived FASTA as an extra Preview tab (not replace).
    openInAdditionalPreviewTab = Signal(object)
    # Open Dotplot (Self) for this record id (top selected when multi-select).
    showSelfDotplotRequested = Signal(str)
    # Display row order changed (Sorting) — Dotplot Sequence list should refresh.
    sequenceOrderChanged = Signal()

    # Split/duplicate windows must shrink below rail+toolbar natural width;
    # overflow and clipping handle the rest.
    _MIN_PANEL_WIDTH = 160
    _MIN_PANEL_HEIGHT = 120

    def __init__(self, fasta_path: Path, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.setMinimumWidth(0)
        self._path = fasta_path.expanduser().resolve()
        self._record_ids: list[str] = []
        self._lengths: dict[str, int] = {}
        self._features: dict[str, list[SequenceFeature]] = {}
        self._annotations_dirty = False
        self._did_initial_fit = False
        self._refit_pending_on_show = False
        self._name_width = _NAME_COLUMN_DEFAULT
        self._syncing_scroll = False
        self._search_needle = ""
        self._search_hits: list[int] = []
        self._search_hit_index = -1
        self._updating_zoom_slider = False
        # Debounce Fit-zoom reflow when the rail, splitter, or window changes width.
        self._viewport_refit_timer = QTimer(self)
        self._viewport_refit_timer.setSingleShot(True)
        self._viewport_refit_timer.setInterval(0)
        self._viewport_refit_timer.timeout.connect(self._refit_after_viewport_change)

        self._names = _NamesColumn(self)
        self._names_scroll = QScrollArea(self)
        self._names_scroll.setWidget(self._names)
        self._names_scroll.setWidgetResizable(False)
        self._names_scroll.setFrameShape(QScrollArea.Shape.NoFrame)
        self._names_scroll.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOff
        )
        # Tracks own the vertical scrollbar; names stay pixel-synced.
        self._names_scroll.setVerticalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOff
        )
        self._names_scroll.setMinimumWidth(_NAME_COLUMN_MIN)

        self._bars = _BarsCanvas(self)
        self._bars_scroll = QScrollArea(self)
        self._bars_scroll.setWidget(self._bars)
        self._bars_scroll.setWidgetResizable(False)
        self._bars_scroll.setFrameShape(QScrollArea.Shape.NoFrame)
        # Always reserve the horizontal scrollbar so viewport height matches names.
        self._bars_scroll.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOn
        )
        self._bars_scroll.setVerticalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAsNeeded
        )

        self._scale = _ScaleRuler(self)

        # Spacer under names matching the tracks h-scrollbar so last rows stay aligned.
        self._hbar_shim = QWidget(self)
        self._hbar_shim.setFixedHeight(
            max(12, self._bars_scroll.horizontalScrollBar().sizeHint().height())
        )
        self._scale_shim = QWidget(self)
        self._scale_shim.setFixedHeight(self._scale.height())

        base = self.palette().color(QPalette.ColorRole.Base)
        for area in (self._names_scroll, self._bars_scroll):
            area.setStyleSheet(
                f"QScrollArea {{ background: {base.name()}; border: none; }}"
            )
            area.viewport().setAutoFillBackground(True)
            vp_pal = area.viewport().palette()
            vp_pal.setColor(QPalette.ColorRole.Window, base)
            area.viewport().setPalette(vp_pal)

        self._names_scroll.verticalScrollBar().valueChanged.connect(
            self._on_names_vscroll
        )
        self._bars_scroll.verticalScrollBar().valueChanged.connect(
            self._on_bars_vscroll
        )
        self._bars_scroll.horizontalScrollBar().valueChanged.connect(self._on_hscroll)
        self._bars.rowsChanged.connect(self._sync_names_column)
        self._bars.featureActivated.connect(self._on_feature_activated)
        self._bars.featureResized.connect(self._on_feature_resized)
        self._bars.zoomChanged.connect(self._on_zoom_changed)
        self._bars.regionZoomRequested.connect(self._on_region_zoom)
        self._bars.regionSelectionChanged.connect(self._on_region_selection_changed)
        self._bars.selectionVisualsDismissed.connect(
            self._on_selection_visuals_dismissed
        )
        self._bars.regionDragPreview.connect(self._on_region_drag_preview)
        self._scale.columnRegionSelected.connect(self._on_scale_column_region)
        self._scale.columnRegionPreview.connect(self._on_scale_column_preview)
        self._scale.columnPixelPreview.connect(self._on_scale_column_pixel_preview)
        self._scale.zoomAtRequested.connect(self._on_scale_zoom_at)
        self._names.rowClicked.connect(self._on_name_row_clicked)
        self._names.rowsDragPreview.connect(self._on_name_rows_preview)
        self._names.rowsPixelPreview.connect(self._on_name_rows_pixel_preview)
        self._names.rowsDragSelected.connect(self._on_name_rows_dragged)
        self._names.selectAllRequested.connect(self._select_all_sequences)
        self._names.deleteRequested.connect(self._on_delete_clicked)
        self._names.contextMenuRequested.connect(self._on_names_context_menu)
        self._bars.contextMenuRequested.connect(self._on_bars_context_menu)
        # Bars owns the vertical scrollbar; names stays synced. Wheel over either
        # pane uses pixelDelta (trackpad) or a gentle mouse-wheel step.
        self._names_scroll.viewport().installEventFilter(self)
        self._bars_scroll.viewport().installEventFilter(self)
        self._bars_scroll.verticalScrollBar().setSingleStep(18)
        self._bars_scroll.verticalScrollBar().setPageStep(120)

        # Zoom / Position / Select / Fit live on AnnotationZoomBar (same row as
        # Text / Annotation). Keep hidden locals until bind, and for standalone use.
        self._header_zoom_bar = None
        self._zoom_bar_header_host: QWidget | None = None
        # Throttle Position spin updates while dragging the scale (~20 Hz).
        self._scale_preview_pending: tuple[int, int] | None = None
        self._scale_preview_timer = QTimer(self)
        self._scale_preview_timer.setSingleShot(True)
        self._scale_preview_timer.setInterval(50)
        self._scale_preview_timer.timeout.connect(self._flush_scale_preview_positions)
        from repeatviewer.widgets.spin_boxes import NoWheelSlider

        self._zoom_slider = NoWheelSlider(Qt.Orientation.Horizontal, self)
        self._zoom_slider.setObjectName("annotationZoomSlider")
        self._zoom_slider.setRange(0, 100)
        self._zoom_slider.setValue(0)
        self._zoom_slider.hide()
        self._zoom_slider.valueChanged.connect(self._on_zoom_slider_changed)
        self._name_filter = QLineEdit(self)
        self._name_filter.setObjectName("FileFilter")
        self._name_filter.setPlaceholderText("Filter by name…")
        self._name_filter.hide()
        self._name_filter.textChanged.connect(self._on_name_filter_edited)
        self._name_filter.returnPressed.connect(self._on_name_filter_next)

        # Host for AnnotationZoomBar — above the complete overview map so its
        # controls start at the same left edge as the sequence names.
        self._toolbar_slot = _OverviewToolbarSlot(self)
        self._toolbar_slot.setObjectName("OverviewMapToolbarSlot")
        self._toolbar_slot.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        self._toolbar_slot.setMinimumWidth(0)
        self._toolbar_slot.setFixedHeight(ANNOTATION_TOOLBAR_RENDERED_HEIGHT)
        self._toolbar_slot_layout = QHBoxLayout(self._toolbar_slot)
        # Match the single-sequence AnnotationMapToolbar exactly: its first
        # control begins at x=6 with no extra host inset.
        self._toolbar_slot_layout.setContentsMargins(0, 0, 0, 0)
        self._toolbar_slot_layout.setSpacing(0)
        self._toolbar_slot_layout.setSizeConstraint(
            QHBoxLayout.SizeConstraint.SetNoConstraint
        )
        self._toolbar_slot.hide()

        names_col = QWidget(self)
        names_layout = QVBoxLayout(names_col)
        names_layout.setContentsMargins(0, 0, 0, 0)
        names_layout.setSpacing(0)
        names_layout.addWidget(self._scale_shim, 0)
        names_layout.addWidget(self._names_scroll, 1)
        names_layout.addWidget(self._hbar_shim)

        tracks_col = QWidget(self)
        tracks_layout = QVBoxLayout(tracks_col)
        tracks_layout.setContentsMargins(0, 0, 0, 0)
        tracks_layout.setSpacing(0)
        tracks_layout.addWidget(self._scale, 0)
        tracks_layout.addWidget(self._bars_scroll, 1)

        names_col.setMinimumWidth(_NAME_COLUMN_MIN)
        tracks_col.setMinimumWidth(40)
        self._splitter = QSplitter(Qt.Orientation.Horizontal, self)
        self._splitter.addWidget(names_col)
        self._splitter.addWidget(tracks_col)
        self._splitter.setStretchFactor(0, 0)
        self._splitter.setStretchFactor(1, 1)
        # Allow narrow panes (split Preview); names floor is still _NAME_COLUMN_MIN.
        self._splitter.setChildrenCollapsible(True)
        self._splitter.splitterMoved.connect(self._on_splitter_moved)

        # Toolbar spans both splitter columns.  The scale and sequence tracks
        # retain their existing name-column offset below it.
        map_host = QWidget(self)
        map_layout = QVBoxLayout(map_host)
        map_layout.setContentsMargins(0, 0, 0, 0)
        map_layout.setSpacing(0)
        map_layout.addWidget(self._toolbar_slot, 0)
        map_layout.addWidget(self._splitter, 1)
        self._map_host = map_host

        self._auto_ann_settings = load_auto_annotation_settings()
        self._similarity_ann_settings = load_similarity_annotation_settings()
        self._bars.set_show_similarity_database_name(
            self._similarity_ann_settings.show_database_name_on_track
        )
        self._ann_cancel = JobCancellation()
        self._ann_thread: QThread | None = None
        self._ann_worker: Worker | None = None
        self._ann_job_kind = ""
        self._ann_progress_total = 0
        self._ann_last_sequence_status = ""
        self._ann_reveal_kinds: tuple[str, ...] = ()
        self._auto_all_updating = False
        self._similarity_all_updating = False
        self._similarity_db_rows: list[tuple] = []
        self._feature_undo_stack: list[_OverviewUndoSnapshot] = []

        self._rail = SequenceAnnotationRail(self)
        self._rail.set_overview_display_mode(True)
        self._rail.residueColorsChanged.connect(self._on_residue_colors_changed)
        self._rail.msaColoringModeChanged.connect(self._on_msa_coloring_mode_changed)
        self._rail.msaDivergenceChanged.connect(self._on_msa_divergence_changed)
        self._rail.colorsClicked.connect(self._on_colors_setting_clicked)
        self._transfer_cancel = JobCancellation()
        self._transfer_thread: QThread | None = None
        self._transfer_worker: Worker | None = None
        self._wire_overview_annotation_rail()
        self._bars.set_show_annotation_track(
            bool(self._rail.show_annotation_track.isChecked())
        )
        self._bars.set_visible_annotation_kinds(
            self._rail.visible_annotation_kinds()
        )

        # Same 1px map↔rail divider as single-seq AnnotationMapWindow.
        from repeatviewer.widgets.sequence_map_viewer import _MAP_RAIL_DIVIDER

        map_rail_divider = QFrame(self)
        map_rail_divider.setObjectName("annotationMapRailDivider")
        map_rail_divider.setFrameShape(QFrame.Shape.VLine)
        map_rail_divider.setFrameShadow(QFrame.Shadow.Plain)
        map_rail_divider.setFixedWidth(1)
        map_rail_divider.setAutoFillBackground(True)
        divider_palette = map_rail_divider.palette()
        divider_palette.setColor(QPalette.ColorRole.Window, _MAP_RAIL_DIVIDER)
        divider_palette.setColor(QPalette.ColorRole.Dark, _MAP_RAIL_DIVIDER)
        map_rail_divider.setPalette(divider_palette)
        map_rail_divider.setStyleSheet(
            f"QFrame#annotationMapRailDivider {{ background-color: {_MAP_RAIL_DIVIDER.name()}; "
            f"border: none; max-width: 1px; }}"
        )

        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        layout.setSizeConstraint(QHBoxLayout.SizeConstraint.SetNoConstraint)
        layout.addWidget(self._map_host, 1)
        layout.addWidget(map_rail_divider, 0)
        layout.addWidget(self._rail, 0)
        self._map_host.setMinimumWidth(0)

    def minimumSizeHint(self) -> QSize:  # noqa: N802
        return QSize(self._MIN_PANEL_WIDTH, self._MIN_PANEL_HEIGHT)

    def _compute_name_width(self, viewport_width: int | None = None) -> int:
        """Keep names compact so the MSA + right control panel get more width."""
        font = QFont(self.font())
        font.setPointSize(max(2, int(self._names.label_point_size())))
        fm = QFontMetrics(font)
        floor = _name_column_min_width(fm, len(self._record_ids))
        widest = 0
        for index, rid in enumerate(self._record_ids[:500]):
            widest = max(
                widest, fm.horizontalAdvance(f"{index + 1}. {rid}")
            )
        natural = widest + 2 * _NAME_PAD + 4
        cap = _NAME_COLUMN_MAX
        if viewport_width and viewport_width > 100:
            # At most ~22% of the overview content (names | tracks), not the rail.
            cap = min(cap, max(floor, int(viewport_width * 0.22)))
        return max(floor, min(cap, natural))

    def _current_name_width(self) -> int:
        sizes = self._splitter.sizes()
        floor = _NAME_COLUMN_MIN
        if sizes and sizes[0] >= floor:
            return int(sizes[0])
        return self._name_width

    def _on_splitter_moved(self, _pos: int = 0, _index: int = 0) -> None:
        self._name_width = self._current_name_width()
        # Keep the names widget flush with the pane so widening the splitter
        # does not leave a blank strip between labels and tracks.
        self._sync_names_column_width()
        self._schedule_viewport_refit()

    def _on_rail_panel_collapsed_changed(self, _collapsed: bool) -> None:
        """Rail hide/show changes track width without resizing this widget."""
        self._schedule_viewport_refit()
        bar = self._header_zoom_bar
        hider = getattr(bar, "_overview_overflow", None) if bar is not None else None
        if hider is not None:
            # Map column width changes when the Display panel opens/closes.
            def _deferred_overflow_sync(overflow: object = hider) -> None:
                sync = getattr(overflow, "sync", None)
                if callable(sync):
                    sync()

            QTimer.singleShot(0, _deferred_overflow_sync)

    def _schedule_viewport_refit(self) -> None:
        self._viewport_refit_timer.start()

    def _is_near_fit_zoom(self) -> bool:
        """True when the view is at (or past) whole-alignment Fit zoom.

        Use actual ``pixels_per_base`` — not the zoom slider alone. The slider
        can briefly lag after wheel/Select zoom; trusting ``value() <= 1`` would
        call Fit during name-column multi-select and jump the MSA to bp 1.
        """
        if not self._record_ids:
            return False
        fit_ppb, _ = self._slider_ppb_bounds()
        return float(self._bars.pixels_per_base()) <= fit_ppb * 1.02

    def _refit_after_viewport_change(self) -> None:
        """Keep Fit-zoom tracks edge-to-edge after rail / splitter / window resize."""
        if not self._overview_panel_alive() or not self._record_ids:
            return
        if not self.isVisible():
            # Defer until show (PDF → MSA). Do not adopt a collapsed viewport size.
            self._refit_pending_on_show = True
            return
        self._refit_pending_on_show = False
        if not self._did_initial_fit:
            self._try_initial_fit()
            return
        near_fit = self._is_near_fit_zoom()
        if near_fit:
            self._fit_bars()
            return
        self._sync_viewport_fill()
        self._sync_names_column()
        self._update_scale()
        self._sync_zoom_slider_from_ppb()
        self._bars.update()

    def _sync_hbar_shim(self) -> None:
        h = max(12, self._bars_scroll.horizontalScrollBar().sizeHint().height())
        if self._hbar_shim.height() != h:
            self._hbar_shim.setFixedHeight(h)

    def _sync_name_label_size_to_letters(self) -> None:
        """Grow sequence-name labels with zoomed-in nucleotide letters."""
        letter_pt = self._bars.effective_letter_point_size()
        desired = max(2, min(96, letter_pt + 1))
        if self._names.label_point_size() != desired:
            self._names.set_label_point_size(letter_pt)

    def _names_viewport_width(self) -> int:
        font = QFont(self.font())
        font.setPointSize(max(2, int(self._names.label_point_size())))
        floor = _name_column_min_width(QFontMetrics(font), len(self._record_ids))
        return max(
            floor,
            self._names_scroll.viewport().width() or self._name_width,
        )

    def _sync_names_column_width(self) -> None:
        """Resize the names widget to the splitter pane (no blank gap to tracks)."""
        font = QFont(self.font())
        font.setPointSize(max(2, int(self._names.label_point_size())))
        floor = _name_column_min_width(QFontMetrics(font), len(self._record_ids))
        if self._names_scroll.minimumWidth() != floor:
            self._names_scroll.setMinimumWidth(floor)
        width = self._names_viewport_width()
        if self._names.width() == width:
            return
        height = max(_MIN_ROW_HEIGHT, self._names.height())
        self._names.resize(width, height)
        self._names.update()

    def _overview_panel_alive(self) -> bool:
        """False once Qt has destroyed this panel (deferred timers must no-op)."""
        try:
            from shiboken6 import isValid

            return bool(isValid(self))
        except Exception:
            try:
                self.objectName()
                return True
            except RuntimeError:
                return False

    def _sync_sequence_count_label(self) -> None:
        """Update Overview toolbar sequence count (right of Cons)."""
        if not self._overview_panel_alive():
            return
        bar = self._header_zoom_bar
        if bar is None or not hasattr(bar, "set_sequence_count"):
            return
        try:
            from shiboken6 import isValid

            if not isValid(bar):
                return
        except Exception:
            pass
        bar.set_sequence_count(len(self._bars.record_ids()))

    def _sync_names_column(self) -> None:
        if not self._overview_panel_alive():
            return
        self._sync_hbar_shim()
        self._sync_name_label_size_to_letters()
        self._sync_sequence_count_label()
        width = self._names_viewport_width()
        # Match bars content; both viewports share the same height thanks to shims.
        content_h = max(
            self._bars.content_height(),
            self._bars_scroll.viewport().height(),
        )
        # Name highlight: whole-row select, regional rows, or live drag preview.
        highlight = self._bars.selected_row_set()
        regional = self._bars.regional_rows()
        if regional:
            highlight = highlight | regional
        drag_rows = self._bars.drag_preview_rows()
        if drag_rows:
            highlight = highlight | drag_rows
        name_preview = self._bars.name_selection_preview_rows()
        if name_preview:
            highlight = highlight | name_preview
        record_ids = self._bars.record_ids()
        tops = self._bars.row_tops()
        heights = self._bars.row_heights()
        # During rubber drag only the highlight set changes — skip geometry rebuild
        # unless the splitter pane width changed (would leave a blank gap).
        target_h = max(_MIN_ROW_HEIGHT, int(content_h))
        if (
            record_ids == self._names._record_ids
            and tops == self._names._row_tops
            and heights == self._names._row_heights
            and self._names.minimumHeight() == target_h
            and self._names.width() == width
        ):
            self._names.set_selected_rows(highlight)
            return
        self._names.set_rows(
            record_ids,
            tops,
            heights,
            highlight,
            width=width,
            content_height=content_h,
        )

    @staticmethod
    def _safe_disconnect(signal, slot) -> None:
        """Disconnect without raising when the slot/receiver is already gone."""
        import warnings

        try:
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    "ignore",
                    message=r".*Failed to disconnect.*",
                    category=RuntimeWarning,
                )
                signal.disconnect(slot)
        except (RuntimeError, TypeError, SystemError):
            pass

    def _restore_zoom_bar_to_header(self, bar: AnnotationZoomBar) -> None:
        """Return AnnotationZoomBar to PreviewHeaderBar after leaving Overview."""
        host = getattr(self, "_zoom_bar_header_host", None)
        bar.set_overview_mode(False)
        bar.hide()
        slot = getattr(self, "_toolbar_slot", None)
        if slot is not None:
            try:
                from shiboken6 import isValid

                if isValid(slot):
                    slot.hide()
            except Exception:
                try:
                    slot.hide()
                except RuntimeError:
                    pass
        if host is None:
            return
        try:
            from shiboken6 import isValid

            if not isValid(host):
                return
        except Exception:
            pass
        layout = host.layout()
        bar.setParent(host)
        if layout is not None and layout.indexOf(bar) < 0:
            layout.insertWidget(0, bar, 0, Qt.AlignmentFlag.AlignVCenter)

    def _mount_zoom_bar_in_panel(self, bar: AnnotationZoomBar) -> None:
        """Place −/+/Position/Select/Fit above names and sequence tracks."""
        from repeatviewer.widgets.sequence_map_viewer import (
            _app_chrome_background,
            _apply_map_paper_background,
        )

        # Remember header host so we can put the bar back on unbind.
        if bar.parent() is not None and bar.parent() is not self._toolbar_slot:
            self._zoom_bar_header_host = bar.parent()
        layout = self._toolbar_slot_layout
        while layout.count():
            item = layout.takeAt(0)
            w = item.widget()
            if w is not None and w is not bar:
                w.setParent(None)
        bar.setParent(self._toolbar_slot)
        if layout.indexOf(bar) < 0:
            layout.addWidget(bar, 1)
        bar.set_overview_mode(True)
        bar.show()
        self._toolbar_slot.show()
        # Overflow must measure the map column, not the preferred toolbar width.
        # Drop any stale sources from a previous Overview page first.
        hider = getattr(bar, "_overview_overflow", None)
        if hider is not None:
            clear = getattr(hider, "clear_width_sources", None)
            if callable(clear):
                clear()
            hider.add_width_source(self._toolbar_slot)
            hider.add_width_source(self._map_host)
            # If this panel is destroyed without unbind, drop our slot immediately.
            slot = self._toolbar_slot

            def _drop_slot_width_source(*_args, overflow=hider, source=slot) -> None:
                remove = getattr(overflow, "remove_width_source", None)
                if callable(remove):
                    remove(source)

            try:
                slot.destroyed.connect(_drop_slot_width_source)
            except (RuntimeError, TypeError):
                pass
            hider.sync()

            def _deferred_overflow_sync(
                overflow: object = hider,
            ) -> None:
                sync = getattr(overflow, "sync", None)
                if callable(sync):
                    sync()

            QTimer.singleShot(0, _deferred_overflow_sync)
        # Reparent leaves the header controls host empty — hide it so it does
        # not keep a blank gray strip under Text/Overview in duplicate windows.
        host = getattr(self, "_zoom_bar_header_host", None)
        if host is not None and host is not self._toolbar_slot:
            host.hide()
            header = host.parent()
            update_h = getattr(header, "_update_header_height", None)
            if callable(update_h):
                update_h()
        # Match single-seq AnnotationMapToolbar chrome (Select/Fit button fill).
        chrome = _app_chrome_background()
        _apply_map_paper_background(self._toolbar_slot, chrome)
        _apply_map_paper_background(bar, chrome)

    def _unbind_header_zoom_bar(self) -> None:
        """Disconnect this panel from the shared header zoom bar (no-op if unbound)."""
        bar = self._header_zoom_bar
        if bar is None:
            return
        # Overflow hider on the shared bar must not keep pointers to this
        # panel's toolbar slot / map host — those are destroyed with the page
        # and otherwise cause ``_OverviewToolbarSlot already deleted``.
        hider = getattr(bar, "_overview_overflow", None)
        if hider is not None:
            clear = getattr(hider, "clear_width_sources", None)
            if callable(clear):
                clear()
            else:
                remove = getattr(hider, "remove_width_source", None)
                if callable(remove):
                    remove(getattr(self, "_toolbar_slot", None))
                    remove(getattr(self, "_map_host", None))
        self._safe_disconnect(bar.zoom_slider.valueChanged, self._on_zoom_slider_changed)
        self._safe_disconnect(bar.fit_requested, self._fit_bars)
        self._safe_disconnect(bar.select_requested, self._select_by_position)
        self._safe_disconnect(bar.pos_start.valueChanged, self._on_position_spin_changed)
        self._safe_disconnect(bar.pos_end.valueChanged, self._on_position_spin_changed)
        self._safe_disconnect(bar.clean_requested, self._open_msa_clean_dialog)
        self._safe_disconnect(bar.consensus_requested, self._open_consensus_dialog)
        self._safe_disconnect(bar.extract_requested, self._extract_selected)
        self._safe_disconnect(bar.filter_text_changed, self._on_name_filter_edited)
        self._safe_disconnect(bar.filter_return_pressed, self._on_name_filter_next)
        if getattr(bar, "_bound_overview", None) is self:
            bar._bound_overview = None
        self._restore_zoom_bar_to_header(bar)
        self._header_zoom_bar = None

    def bind_annotation_zoom_bar(self, bar) -> None:
        """Attach −/+/slider/Position/Select/Fit beside the rail (single-seq top)."""
        if not isinstance(bar, AnnotationZoomBar):
            return

        prev_owner = getattr(bar, "_bound_overview", None)
        if prev_owner is self and self._header_zoom_bar is bar:
            self._mount_zoom_bar_in_panel(bar)
            bar.set_msa_mode(self.is_msa_mode())
            self._sync_sequence_count_label()
            self._sync_header_position_spins()
            return
        if prev_owner is not None and prev_owner is not self:
            try:
                from shiboken6 import isValid

                if isValid(prev_owner):
                    prev_owner._unbind_header_zoom_bar()
                else:
                    bar._bound_overview = None
            except (RuntimeError, TypeError, AttributeError, SystemError):
                bar._bound_overview = None
        if self._header_zoom_bar is not None:
            self._unbind_header_zoom_bar()

        prev_value = self._zoom_slider.value()
        prev_filter = self._name_filter.text()

        if self._zoom_slider is not bar.zoom_slider:
            self._safe_disconnect(
                self._zoom_slider.valueChanged, self._on_zoom_slider_changed
            )
        if self._name_filter is not bar.name_filter:
            self._safe_disconnect(
                self._name_filter.textChanged, self._on_name_filter_edited
            )
            self._safe_disconnect(
                self._name_filter.returnPressed, self._on_name_filter_next
            )

        self._header_zoom_bar = bar
        bar._bound_overview = self
        bar.set_msa_mode(self.is_msa_mode())
        self._sync_sequence_count_label()
        self._mount_zoom_bar_in_panel(bar)
        self._zoom_slider = bar.zoom_slider
        self._name_filter = bar.name_filter
        self._updating_zoom_slider = True
        try:
            self._zoom_slider.setValue(prev_value)
        finally:
            self._updating_zoom_slider = False
        self._name_filter.blockSignals(True)
        self._name_filter.setText(prev_filter)
        self._name_filter.blockSignals(False)
        self._zoom_slider.valueChanged.connect(self._on_zoom_slider_changed)
        bar.fit_requested.connect(self._fit_bars)
        bar.select_requested.connect(self._select_by_position)
        bar.pos_start.valueChanged.connect(self._on_position_spin_changed)
        bar.pos_end.valueChanged.connect(self._on_position_spin_changed)
        bar.clean_requested.connect(self._open_msa_clean_dialog)
        bar.consensus_requested.connect(self._open_consensus_dialog)
        bar.extract_requested.connect(self._extract_selected)
        bar.filter_text_changed.connect(self._on_name_filter_edited)
        bar.filter_return_pressed.connect(self._on_name_filter_next)
        self._sync_header_position_spins()

    def _sync_mode_specific_controls(self) -> None:
        msa = self.is_msa_mode()
        self._rail.set_msa_display_mode(msa)
        self._rail.set_nucleotide_color_context("msa" if msa else "multi")
        self._rail.sync_color_scheme_combos()
        settings = self._overview_display_settings()
        try:
            key = "sequence_distance" if msa else "multi_sequence_distance"
            default = (
                _DEFAULT_MSA_SEQUENCE_DISTANCE
                if msa
                else _DEFAULT_MULTI_SEQUENCE_DISTANCE
            )
            distance = max(0, min(40, int(settings.value(key, default))))
            color_background = settings.value(
                "multi_sequence_color_background", True, type=bool
            )
        finally:
            settings.endGroup()
        self._rail.sequence_distance.blockSignals(True)
        try:
            self._rail.sequence_distance.setValue(distance)
        finally:
            self._rail.sequence_distance.blockSignals(False)
        self._bars.set_sequence_row_gap(distance)
        self._rail.color_multi_sequence_background.blockSignals(True)
        try:
            self._rail.color_multi_sequence_background.setChecked(color_background)
        finally:
            self._rail.color_multi_sequence_background.blockSignals(False)
        self._bars.set_color_unaligned_background(color_background)
        if self._header_zoom_bar is not None:
            self._header_zoom_bar.set_msa_mode(msa)

    def is_msa_mode(self) -> bool:
        return self._bars.is_msa_mode()

    def _wire_overview_annotation_rail(self) -> None:
        """Connect Annotation / Display actions for multi-seq overview."""
        rail = self._rail
        self._auto_ann_checks: dict[AutoAnnotationMode, QCheckBox] = {}
        for label, mode in _AUTO_ANNOTATION_MODES:
            check = rail.auto_ann_checks.get(label)
            if check is None:
                continue
            self._auto_ann_checks[mode] = check
            check.toggled.connect(self._on_auto_annotation_check_toggled)
        rail.auto_ann_all_check.toggled.connect(self._on_auto_annotation_all_toggled)
        rail.similarity_ann_all_check.toggled.connect(
            self._on_similarity_annotation_all_toggled
        )
        rail.featureDetectClicked.connect(self._on_auto_run_clicked)
        rail.featureDetectSettingsClicked.connect(self._open_auto_annotation_settings)
        rail.similarityAnnotateClicked.connect(self._on_similarity_run_clicked)
        rail.similaritySettingsClicked.connect(
            self._open_similarity_annotation_settings
        )
        rail.saveClicked.connect(self._save_overview_annotations)
        rail.deleteClicked.connect(self._on_delete_clicked)
        rail.restoreClicked.connect(self._restore_features)
        rail.exportClicked.connect(self._export_overview_view)
        self._bars.deleteRequested.connect(self._on_delete_clicked)
        rail.addFeatureClicked.connect(self._on_add_feature_clicked)
        rail.addFeatureSettingsClicked.connect(
            self._open_add_feature_buttons_settings
        )
        rail.transferRequested.connect(self._on_transfer_annotations_clicked)
        rail.trackHeightChanged.connect(self._on_track_height_changed)
        rail.annotationThicknessChanged.connect(self._on_annotation_thickness_changed)
        rail.annotationDistanceChanged.connect(self._on_annotation_distance_changed)
        rail.annotationLabelSizeChanged.connect(self._on_annotation_label_size_changed)
        rail.baseNumberSizeChanged.connect(self._on_scale_label_size_changed)
        rail.sequenceDistanceChanged.connect(self._on_sequence_row_gap_changed)
        rail.columnDistanceChanged.connect(self._on_column_distance_changed)
        rail.multiSequenceBackgroundChanged.connect(
            self._on_multi_sequence_background_changed
        )
        rail.showAnnotationTrackChanged.connect(self._on_show_annotation_track_changed)
        rail.annotationTrackKindsChanged.connect(
            self._on_annotation_track_kinds_changed
        )
        rail.sequenceSortChanged.connect(self._on_sequence_sort_changed)
        rail.findRequested.connect(self._on_rail_find_requested)
        rail.panelCollapsedChanged.connect(self._on_rail_panel_collapsed_changed)
        self._find_items: list = []
        self._find_index: int = -1
        self._find_key: object | None = None
        self._load_overview_display_settings()
        track_h = int(rail.track_height.value())
        self._bars.set_track_height(track_h)
        self._names.set_label_point_size(track_h)  # names store height+1
        self._bars.set_feature_bar_height(int(rail.annotation_thickness.value()))
        self._bars.set_feature_gap(int(rail.annotation_distance.value()))
        self._bars.set_feature_label_point_size(int(rail.annotation_label_size.value()))
        self._bars.set_sequence_row_gap(int(rail.sequence_distance.value()))
        self._scale.set_font_point_size(int(rail.base_number_size.value()))
        self._scale_shim.setFixedHeight(self._scale.height())
        # Column spacing is no longer a user-facing MSA option.
        self._bars.set_column_distance(0)
        self._bars.set_color_unaligned_background(
            rail.color_multi_sequence_background.isChecked()
        )
        self._bars.set_msa_color_mode(rail.msa_coloring_mode_id())
        self._bars.set_msa_divergence_threshold(rail.msa_divergence_threshold())
        self._sync_auto_annotation_all_check()
        self._refresh_similarity_database_rows()

    def _overview_display_settings(self) -> QSettings:
        s = QSettings("RepeatViewer", "RepeatViewer")
        s.beginGroup(_OVERVIEW_DISPLAY_SETTINGS_GROUP)
        return s

    def _load_overview_display_settings(self) -> None:
        """MSA/overview defaults: Track height 4, MSA sequence distance 1."""
        rail = self._rail
        s = self._overview_display_settings()
        try:
            # Prefer track_height; migrate old bases_zoomed_in into 1–10 if needed.
            if s.contains("track_height"):
                track_h = int(s.value("track_height", _DEFAULT_TRACK_HEIGHT))
            else:
                legacy = int(s.value("bases_zoomed_in", _DEFAULT_TRACK_HEIGHT))
                track_h = max(1, min(10, legacy if legacy <= 10 else _DEFAULT_TRACK_HEIGHT))
            distance = int(
                s.value("sequence_distance", _DEFAULT_MSA_SEQUENCE_DISTANCE)
            )
            ann_thickness = int(s.value("annotation_thickness", _FEAT_HEIGHT))
            ann_distance = int(s.value("annotation_distance", _FEAT_GAP))
            ann_labels = int(s.value("annotation_label_size", _FEAT_LABEL_POINT_SIZE))
            scale_labels = int(s.value("base_number_size", 8))
            col_distance = int(
                s.value("column_distance", _DEFAULT_MSA_COLUMN_DISTANCE)
            )
            msa_mode = normalize_msa_color_mode(
                str(s.value("msa_coloring_mode", MSA_COLOR_ALL_BACKGROUNDS) or "")
            )
            msa_div = normalize_msa_divergence_threshold(
                s.value(
                    "msa_divergence_threshold",
                    DEFAULT_MSA_DIVERGENCE_THRESHOLD,
                )
            )
            sort_mode = str(s.value("sequence_sort_mode", "original") or "original")
        finally:
            s.endGroup()
        track_h = max(1, min(10, track_h))
        distance = max(0, min(40, distance))
        ann_thickness = max(2, min(40, ann_thickness))
        ann_distance = max(0, min(40, ann_distance))
        ann_labels = max(0, min(24, ann_labels))
        scale_labels = max(4, min(24, scale_labels))
        col_distance = max(0, min(40, col_distance))
        for widget in (
            rail.track_height,
            rail.sequence_distance,
            rail.annotation_thickness,
            rail.annotation_distance,
            rail.annotation_label_size,
            rail.base_number_size,
            rail.column_distance,
        ):
            widget.blockSignals(True)
        try:
            rail.track_height.setValue(track_h)
            rail.sequence_distance.setValue(distance)
            rail.annotation_thickness.setValue(ann_thickness)
            rail.annotation_distance.setValue(ann_distance)
            rail.annotation_label_size.setValue(ann_labels)
            rail.base_number_size.setValue(scale_labels)
            rail.column_distance.setValue(col_distance)
        finally:
            for widget in (
                rail.track_height,
                rail.sequence_distance,
                rail.annotation_thickness,
                rail.annotation_distance,
                rail.annotation_label_size,
                rail.base_number_size,
                rail.column_distance,
            ):
                widget.blockSignals(False)
        rail.set_msa_coloring_mode(msa_mode)
        rail.set_msa_divergence_threshold(msa_div)
        self._bars.set_msa_color_mode(msa_mode)
        self._bars.set_msa_divergence_threshold(msa_div)
        rail.set_sequence_sort_mode(sort_mode)

    def _save_overview_display_settings(self) -> None:
        rail = self._rail
        s = self._overview_display_settings()
        try:
            s.setValue("track_height", int(rail.track_height.value()))
            distance_key = (
                "sequence_distance"
                if self.is_msa_mode()
                else "multi_sequence_distance"
            )
            s.setValue(distance_key, int(rail.sequence_distance.value()))
            s.setValue("annotation_thickness", int(rail.annotation_thickness.value()))
            s.setValue("annotation_distance", int(rail.annotation_distance.value()))
            s.setValue("annotation_label_size", int(rail.annotation_label_size.value()))
            s.setValue("base_number_size", int(rail.base_number_size.value()))
            s.setValue("column_distance", int(rail.column_distance.value()))
            s.setValue(
                "multi_sequence_color_background",
                rail.color_multi_sequence_background.isChecked(),
            )
            s.setValue("msa_coloring_mode", rail.msa_coloring_mode_id())
            s.setValue(
                "msa_divergence_threshold",
                float(rail.msa_divergence_threshold()),
            )
            s.setValue("sequence_sort_mode", rail.sequence_sort_mode())
        finally:
            s.endGroup()

    def _on_track_height_changed(self, level: int) -> None:
        lvl = max(1, min(10, int(level)))
        self._bars.set_track_height(lvl)
        self._names.set_label_point_size(lvl)  # names store height+1
        self._save_overview_display_settings()
        self._sync_names_column()

    def _on_annotation_thickness_changed(self, height: int) -> None:
        self._bars.set_feature_bar_height(int(height))
        self._save_overview_display_settings()
        self._sync_names_column()

    def _on_annotation_distance_changed(self, gap: int) -> None:
        self._bars.set_feature_gap(int(gap))
        self._save_overview_display_settings()
        self._sync_names_column()

    def _on_annotation_label_size_changed(self, point_size: int) -> None:
        self._bars.set_feature_label_point_size(int(point_size))
        self._save_overview_display_settings()

    def _on_scale_label_size_changed(self, point_size: int) -> None:
        self._scale.set_font_point_size(int(point_size))
        self._scale_shim.setFixedHeight(self._scale.height())
        self._save_overview_display_settings()

    def _on_sequence_row_gap_changed(self, gap: int) -> None:
        self._bars.set_sequence_row_gap(int(gap))
        self._save_overview_display_settings()
        self._sync_names_column()

    def _on_column_distance_changed(self, gap: int) -> None:
        self._bars.set_column_distance(int(gap))
        self._save_overview_display_settings()

    def _on_multi_sequence_background_changed(self, enabled: bool) -> None:
        self._bars.set_color_unaligned_background(bool(enabled))
        self._save_overview_display_settings()

    def _on_show_annotation_track_changed(self, enabled: bool) -> None:
        self._bars.set_show_annotation_track(bool(enabled))
        self._sync_names_column()

    def _on_annotation_track_kinds_changed(self, keys: object) -> None:
        if keys is None:
            visible: frozenset[str] | None = None
        else:
            visible = frozenset(str(k) for k in keys)
        self._bars.set_visible_annotation_kinds(visible)
        self._sync_names_column()

    def _sync_annotation_kind_filters(self) -> None:
        """Refresh type checkboxes from loaded features and apply to the canvas."""
        all_feats: list[SequenceFeature] = []
        for rid in self._record_ids:
            all_feats.extend(self._features.get(rid) or [])
        self._rail.set_available_annotation_kinds(all_feats)
        self._bars.set_visible_annotation_kinds(
            self._rail.visible_annotation_kinds()
        )

    def _on_sequence_sort_changed(self, mode: str) -> None:
        self._bars.sort_rows(str(mode or "original"))
        # Bars own the display order; keep panel ids aligned for Dotplot / tools.
        self._record_ids = list(self._bars.record_ids())
        self._save_overview_display_settings()
        self._sync_names_column()
        self.sequenceOrderChanged.emit()

    def _select_all_sequences(self) -> None:
        self._bars.select_all_rows()

    def _on_selection_visuals_dismissed(self) -> None:
        """Regional/column/row drag started — drop scale highlight only."""
        self._scale.set_region_selection(None)

    def _on_region_selection_changed(self, region: object) -> None:
        # Scale band follows column span for full-height and regional selects.
        if region is None:
            self._scale.set_region_selection(None)
        else:
            try:
                self._scale.set_region_selection(int(region[0]), int(region[1]))  # type: ignore[index]
            except (TypeError, ValueError, IndexError):
                self._scale.set_region_selection(None)
        # Selected-column sort is one-shot from the Sorting combo — do not
        # re-order rows every time the column selection changes.
        bar = self._header_zoom_bar
        if bar is None or not hasattr(bar, "set_position_range"):
            return
        maximum = max(1, int(self._bars.max_length()))
        if region is None:
            # Cleared by click outside / Escape — reset Position spins to full span.
            bar.set_position_range(maximum, start=1, end=maximum)
            return
        try:
            start, end = int(region[0]), int(region[1])  # type: ignore[index]
        except (TypeError, ValueError, IndexError):
            return
        bar.set_position_range(maximum, start=start, end=end)

    def _on_region_drag_preview(self, payload: object) -> None:
        """Live Position spins + name/scale highlights while dragging a regional box."""
        if payload is None:
            return
        try:
            start, end, _rows = payload  # type: ignore[misc]
            start_i, end_i = int(start), int(end)
        except (TypeError, ValueError):
            return
        bar = self._header_zoom_bar
        if bar is not None and hasattr(bar, "set_position_range"):
            maximum = max(1, int(self._bars.max_length()))
            bar.set_position_range(maximum, start=start_i, end=end_i)
        # Live scale band for the column span (matches committed regional select).
        self._scale.set_region_selection(start_i, end_i)
        self._sync_names_column()

    def select_all_sequences(self) -> None:
        """Public hook for Edit → Select All."""
        self._select_all_sequences()

    def _on_residue_colors_changed(self) -> None:
        self._bars.refresh_msa_colors()

    def _on_msa_coloring_mode_changed(self, mode: str) -> None:
        from repeatviewer.backends.msa_coordinates import MSA_COLOR_DIFF_FROM_SELECTED

        self._bars.set_msa_color_mode(mode)
        self._save_overview_display_settings()
        if normalize_msa_color_mode(mode) == MSA_COLOR_DIFF_FROM_SELECTED:
            std = self._bars.msa_diff_standard_row()
            if std is None:
                self.statusMessage.emit(
                    "Color vs standard: no sequences loaded."
                )
            else:
                rid = self._bars.record_ids()[std]
                note = (
                    ""
                    if self._bars.msa_diff_standard_is_pinned()
                    else " (first sequence)"
                )
                self.statusMessage.emit(
                    f"Color vs standard: {rid}{note}."
                )

    def _msa_diff_standard_path_key(self) -> str:
        return str(Path(self._path).expanduser().resolve())

    def _persist_msa_diff_standard_id(self, record_id: str | None) -> None:
        """Remember Color vs standard per MSA path (survives reopen / duplicate)."""
        s = self._overview_display_settings()
        try:
            s.beginGroup(_MSA_DIFF_STANDARD_BY_PATH_GROUP)
            key = self._msa_diff_standard_path_key()
            text = str(record_id or "").strip()
            if text:
                s.setValue(key, text)
            else:
                s.remove(key)
        finally:
            s.endGroup()
            s.endGroup()
        s.sync()

    def _restore_msa_diff_standard_from_settings(self) -> None:
        """Apply a saved Color vs standard id after MSA rows are loaded."""
        if not self.is_msa_mode() or not self._record_ids:
            return
        s = self._overview_display_settings()
        try:
            s.beginGroup(_MSA_DIFF_STANDARD_BY_PATH_GROUP)
            raw = s.value(self._msa_diff_standard_path_key(), "")
        finally:
            s.endGroup()
            s.endGroup()
        rid = str(raw or "").strip()
        if not rid:
            return
        if rid in self._bars.record_ids():
            self._bars.set_msa_diff_standard_id(rid)
        else:
            # Sequence gone from this FASTA — drop the stale pin.
            self._persist_msa_diff_standard_id(None)

    def _set_selected_as_msa_standard(self) -> None:
        """Right-click → Set as standard for Color vs standard mode."""
        rows = self._bars.selected_rows()
        if len(rows) != 1:
            return
        row = int(rows[0])
        self._bars.set_msa_diff_standard_row(row)
        rid = self._bars.msa_diff_standard_id()
        self._persist_msa_diff_standard_id(rid)
        if rid:
            self.statusMessage.emit(f"Color vs standard: {rid}.")

    def _on_msa_divergence_changed(self, value: float) -> None:
        self._bars.set_msa_divergence_threshold(value)
        self._save_overview_display_settings()

    def _on_colors_setting_clicked(self) -> None:
        from repeatviewer.widgets.annotation_colors_dialog import (
            AnnotationElementColorsDialog,
        )

        dialog = AnnotationElementColorsDialog(
            parent=self,
            nucleotide_color_context=self._rail.nucleotide_color_context(),
        )
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        self._rail.sync_color_scheme_combos()
        self._rail.refresh_feature_detect_check_colors()
        self._rail.refresh_manual_feature_check_colors()
        self._on_residue_colors_changed()

    def _open_add_feature_buttons_settings(self) -> None:
        from repeatviewer.widgets.annotation_colors_dialog import (
            AddFeatureButtonsDialog,
        )

        dialog = AddFeatureButtonsDialog(parent=self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        self._rail.populate_manual_feature_types(dialog.buttons())
        self.statusMessage.emit("Feature types updated.")

    def _annotation_running(self) -> bool:
        thread = self._ann_thread
        return thread is not None and thread.isRunning()

    def _set_annotation_job_status(self, text: str) -> None:
        self._rail.set_annotation_job_status(text)
        if text:
            self.statusMessage.emit(text)

    def _clear_annotation_job_status(self) -> None:
        self._rail.set_annotation_job_status("")

    def _job_progress_hosts(self) -> list:
        from repeatviewer.widgets.task_progress import iter_job_progress_hosts

        return iter_job_progress_hosts(self)

    def _main_window_for_job_progress(self):
        """Prefer MainWindow (tool lock + status bar); else the local window."""
        hosts = self._job_progress_hosts()
        for win in hosts:
            if hasattr(win, "_set_job_actions_enabled"):
                return win
        return hosts[0] if hosts else None

    def _begin_window_job_progress(self, title: str, total_steps: int) -> None:
        steps = max(1, int(total_steps))
        for win in self._job_progress_hosts():
            win._show_job_progress(title, steps)
            if hasattr(win, "_set_job_actions_enabled"):
                win._set_job_actions_enabled(False)

    def _update_window_job_progress(self, completed: int, detail: str = "") -> None:
        value = max(0, int(completed))
        text = str(detail or "")
        for win in self._job_progress_hosts():
            if hasattr(win, "_on_job_progress"):
                win._on_job_progress(value, text)

    def _end_window_job_progress(self) -> None:
        for win in self._job_progress_hosts():
            if getattr(win, "_job_running", False) and hasattr(win, "_task_progress"):
                bar = win._task_progress
                if bar is not None:
                    total = max(1, int(getattr(bar, "_total", 1)))
                    bar.update(total, "Done")
            if hasattr(win, "_hide_job_progress"):
                win._hide_job_progress()
            # Unlock BLAST/TE-Aid unless a MainWindow QThread job is still active.
            # Do not call ``_any_background_job_running`` here — it includes this
            # panel and can keep the toolbar disabled while we are finishing.
            if hasattr(win, "_set_job_actions_enabled"):
                worker_busy = (
                    getattr(win, "_worker_thread", None) is not None
                    or getattr(win, "_genome_map_thread", None) is not None
                    or getattr(win, "_benchmark_thread", None) is not None
                )
                if not worker_busy:
                    win._set_job_actions_enabled(True)

    def _sync_auto_annotation_all_check(self) -> None:
        # "All" covers LTR/TIR/Tandem/ORF; TSD is optional (MSA-only) and independent.
        checks = [
            check
            for mode, check in self._auto_ann_checks.items()
            if mode != AutoAnnotationMode.TSD
        ]
        if not checks:
            return
        all_checked = all(check.isChecked() for check in checks)
        self._auto_all_updating = True
        self._rail.auto_ann_all_check.setChecked(all_checked)
        self._auto_all_updating = False

    def _on_auto_annotation_check_toggled(self, _checked: bool) -> None:
        if self._auto_all_updating:
            return
        self._sync_auto_annotation_all_check()

    def _on_auto_annotation_all_toggled(self, checked: bool) -> None:
        if self._auto_all_updating:
            return
        self._auto_all_updating = True
        try:
            for mode, check in self._auto_ann_checks.items():
                if mode == AutoAnnotationMode.TSD:
                    continue
                check.setChecked(checked)
        finally:
            self._auto_all_updating = False

    def _sync_similarity_annotation_all_check(self) -> None:
        checks = [check for _entry, check in self._similarity_db_rows]
        if not checks:
            self._similarity_all_updating = True
            self._rail.similarity_ann_all_check.setChecked(False)
            self._similarity_all_updating = False
            return
        all_checked = all(check.isChecked() for check in checks)
        self._similarity_all_updating = True
        self._rail.similarity_ann_all_check.setChecked(all_checked)
        self._similarity_all_updating = False

    def _on_similarity_annotation_all_toggled(self, checked: bool) -> None:
        if self._similarity_all_updating:
            return
        self._similarity_all_updating = True
        try:
            for _entry, check in self._similarity_db_rows:
                check.setChecked(checked)
        finally:
            self._similarity_all_updating = False

    def _on_similarity_db_check_toggled(self, _checked: bool) -> None:
        if self._similarity_all_updating:
            return
        self._sync_similarity_annotation_all_check()

    def _make_similarity_db_check(self, entry, index: int) -> QCheckBox:
        color = similarity_database_color(entry.name, index)
        check = SimilarityDatabaseCheck(
            entry.name, index, color, self._rail.similarity_db_host
        )
        readiness = database_search_readiness(entry)
        tooltip = f"{index + 1}. {entry.name}\n{entry.display_type_label}\n{entry.path}"
        if readiness is not None:
            tooltip += f"\n\nNot ready: {readiness}"
            check.setChecked(False)
        else:
            check.setChecked(True)
        check.setToolTip(tooltip)
        check.toggled.connect(self._on_similarity_db_check_toggled)
        return check

    def _refresh_similarity_database_rows(self) -> None:
        self._rail.clear_similarity_db_checks()
        self._similarity_db_rows = []
        databases = load_similarity_annotation_databases()
        if not databases:
            empty = QLabel(
                "No similarity databases configured.", self._rail.similarity_db_host
            )
            empty.setStyleSheet("color: palette(mid);")
            empty.setWordWrap(False)
            self._rail.add_similarity_db_widget(empty)
            self._sync_similarity_annotation_all_check()
            return
        for index, entry in enumerate(databases):
            check = self._make_similarity_db_check(entry, index)
            self._rail.add_similarity_db_check(check)
            self._similarity_db_rows.append((entry, check))
        self._sync_similarity_annotation_all_check()

    def _selected_auto_annotation_modes(self) -> list[AutoAnnotationMode]:
        return [
            mode for mode, check in self._auto_ann_checks.items() if check.isChecked()
        ]

    def _selected_similarity_databases(self) -> list[SimilarityDatabaseRun]:
        selected: list[SimilarityDatabaseRun] = []
        for index, (entry, check) in enumerate(self._similarity_db_rows):
            if check.isChecked():
                selected.append(
                    SimilarityDatabaseRun(
                        database=entry,
                        color=similarity_database_color(entry.name, index),
                    )
                )
        return selected

    def _set_annotation_controls_enabled(self, enabled: bool) -> None:
        rail = self._rail
        rail.auto_ann_settings_btn.setEnabled(enabled)
        rail.auto_ann_all_check.setEnabled(enabled)
        for check in self._auto_ann_checks.values():
            check.setEnabled(enabled)
        rail.similarity_ann_settings_btn.setEnabled(enabled)
        rail.similarity_ann_all_check.setEnabled(enabled)
        for _entry, check in self._similarity_db_rows:
            check.setEnabled(enabled)
        running = not enabled
        rail.auto_ann_run_btn.setEnabled(True)
        rail.similarity_ann_run_btn.setEnabled(True)
        rail.auto_ann_run_btn.setText("Stop" if running and self._ann_job_kind == "auto" else "Run")
        rail.similarity_ann_run_btn.setText(
            "Stop" if running and self._ann_job_kind == "similarity" else "Run"
        )

    def _open_auto_annotation_settings(self) -> None:
        from repeatviewer.job_options import save_tsd_search_options
        from repeatviewer.widgets.sequence_map_viewer import AutoAnnotationSettingsDialog

        self._auto_ann_settings = load_auto_annotation_settings()
        dialog = AutoAnnotationSettingsDialog(self._auto_ann_settings, parent=self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        self._auto_ann_settings = dialog.settings()
        save_auto_annotation_settings(self._auto_ann_settings)
        save_tsd_search_options(dialog.tsd_options())
        self.statusMessage.emit("Auto annotation settings saved.")

    def _open_similarity_annotation_settings(self) -> None:
        from repeatviewer.widgets.sequence_map_viewer import (
            SimilarityAnnotationSettingsDialog,
        )

        self._similarity_ann_settings = load_similarity_annotation_settings()
        dialog = SimilarityAnnotationSettingsDialog(
            self._similarity_ann_settings, parent=self
        )
        dialog.manageDatabaseClicked.connect(
            lambda: self._on_manage_database_clicked(host=dialog)
        )
        if dialog.exec() != QDialog.DialogCode.Accepted:
            self._refresh_similarity_database_rows()
            return
        self._similarity_ann_settings = dialog.settings()
        save_similarity_annotation_settings(self._similarity_ann_settings)
        self._bars.set_show_similarity_database_name(
            self._similarity_ann_settings.show_database_name_on_track
        )
        self._refresh_similarity_database_rows()
        self.statusMessage.emit("Similarity annotation settings saved.")

    def _on_manage_database_clicked(self, host: QWidget | None = None) -> None:
        window = self.window()
        opener = getattr(window, "_open_manage_references", None)
        if callable(opener):
            opener(role_filter="blast", host=host)
        self._refresh_similarity_database_rows()

    def _require_selected_record_ids(self, title: str) -> list[str] | None:
        ids = self._selected_record_ids()
        if ids:
            return ids
        QMessageBox.information(
            self,
            title,
            "Select one or more sequences first.\n\n"
            "Tip: Ctrl+A (Cmd+A on Mac) selects all sequences.",
        )
        return None

    def _target_record_ids_for_msa_job(self) -> list[str]:
        """Selected sequences, or every sequence in the open MSA when none are selected.

        Name-column selection and regional MSA row boxes both count as a selection.
        """
        ids = self._selected_record_ids()
        if ids:
            return ids
        regional = self._bars.regional_rows()
        if regional:
            return [
                self._record_ids[row]
                for row in sorted(regional)
                if 0 <= row < len(self._record_ids)
            ]
        return list(self._record_ids)

    def _target_record_ids_for_manual_add(self) -> list[str] | None:
        """Name-selected sequences, or rows inside a regional MSA selection box.

        Unlike Feature detection / Similarity, manual Add never broadcasts to the
        whole MSA — the user must pick sequence(s) or drag a box on them.
        """
        ids = self._selected_record_ids()
        if ids:
            return ids
        regional = self._bars.regional_rows()
        if not regional:
            return None
        return [
            self._record_ids[row]
            for row in sorted(regional)
            if 0 <= row < len(self._record_ids)
        ]

    def _on_auto_run_clicked(self) -> None:
        if self._annotation_running():
            if self._ann_job_kind == "auto":
                self._stop_annotation_job()
            return
        modes = self._selected_auto_annotation_modes()
        if not modes:
            QMessageBox.information(
                self,
                "Feature detection",
                "Select at least one detector to run.",
            )
            return
        target_ids = self._target_record_ids_for_msa_job()
        if not target_ids:
            QMessageBox.information(
                self,
                "Feature detection",
                "No sequences are loaded in this MSA.",
            )
            return
        self._start_batch_annotation(
            TeAnnotationOptions(auto_modes=tuple(modes)),
            kind="auto",
            title="Feature detection",
            total_steps=max(1, len(target_ids)),
            record_ids=target_ids,
        )

    def _on_similarity_run_clicked(self) -> None:
        if self._annotation_running():
            if self._ann_job_kind == "similarity":
                self._stop_annotation_job()
            return
        databases = self._selected_similarity_databases()
        if not databases:
            QMessageBox.information(
                self,
                "Similarity annotation",
                "Select at least one database to search.",
            )
            return
        issues = []
        for entry in databases:
            readiness = database_search_readiness(entry.database)
            if readiness is not None:
                issues.append(f"{entry.database.name}: {readiness}")
        if issues:
            QMessageBox.warning(
                self,
                "Similarity annotation",
                "These databases are not ready to search:\n\n" + "\n\n".join(issues),
            )
            return
        target_ids = self._target_record_ids_for_msa_job()
        if not target_ids:
            QMessageBox.information(
                self,
                "Similarity annotation",
                "No sequences are loaded in this MSA.",
            )
            return
        self._start_batch_annotation(
            TeAnnotationOptions(similarity_databases=tuple(databases)),
            kind="similarity",
            title="Similarity annotation",
            total_steps=max(1, len(target_ids)),
            record_ids=target_ids,
        )

    def _stop_annotation_job(self) -> None:
        if not self._annotation_running():
            return
        self._ann_cancel.request()
        self._set_annotation_job_status("Stopping annotation…")
        for win in self._job_progress_hosts():
            if hasattr(win, "_task_progress") and win._task_progress:
                win._task_progress.update(
                    getattr(win, "_job_work_steps", 0), "Stopping…"
                )
            # Keep Run as Stop until the worker exits; disable only the status Stop.
            if getattr(win, "_stop_btn", None) is not None:
                win._stop_btn.setEnabled(False)

    def cancel_auto_annotation(self) -> bool:
        """Request cancellation of a running Overview feature-detection job."""
        if not self._annotation_running() or self._ann_job_kind != "auto":
            return False
        self._stop_annotation_job()
        return True

    def cancel_similarity_annotation(self) -> bool:
        """Request cancellation of a running Overview similarity job."""
        if not self._annotation_running() or self._ann_job_kind != "similarity":
            return False
        self._stop_annotation_job()
        return True

    def _start_batch_annotation(
        self,
        options: TeAnnotationOptions,
        *,
        kind: str,
        title: str,
        total_steps: int,
        record_ids: list[str] | None = None,
    ) -> None:
        if self._annotation_running():
            return
        if not self._path.is_file():
            QMessageBox.warning(self, title, "No FASTA file loaded.")
            return
        target_ids = list(record_ids or self._selected_record_ids())
        if not target_ids:
            return
        # Refuse to start while MainWindow is already running BLAST / TE-Aid / etc.
        main = self._main_window_for_job_progress()
        if main is not None and hasattr(main, "_any_background_job_running"):
            if main._any_background_job_running():
                QMessageBox.information(
                    self,
                    "Job already running",
                    "Wait for the current job to finish before starting another. "
                    "You can still browse folders and preview files while it runs.",
                )
                return

        self._ann_cancel.reset()
        self._ann_job_kind = kind
        self._ann_progress_total = max(1, int(total_steps))
        self._ann_last_sequence_status = ""
        # Remember which Displayed-annotations kinds this run should turn on.
        if kind == "similarity":
            self._ann_reveal_kinds = tuple(
                f"domain::{run.database.name}"
                for run in options.similarity_databases
                if str(run.database.name or "").strip()
            )
        else:
            self._ann_reveal_kinds = tuple(
                str(mode.value) for mode in options.auto_modes
            )
        self._set_annotation_controls_enabled(False)
        n = len(target_ids)
        scope = (
            f"{n:,} selected sequence{'s' if n != 1 else ''}"
            if self._selected_record_ids() or self._bars.regional_rows()
            else f"all {n:,} sequence{'s' if n != 1 else ''}"
        )
        self._set_annotation_job_status(f"{title}: running on {scope}…")
        self.statusMessage.emit(
            f"Running {title.lower()} on {scope} in {self._path.name}…"
        )
        self._begin_window_job_progress(title, total_steps)

        path = self._path
        auto_cfg = self._auto_ann_settings
        sim_cfg = self._similarity_ann_settings
        cancel = self._ann_cancel
        job_ids = list(target_ids)

        # Snapshot current (possibly unsaved) features so the worker does not
        # discard Preview edits by reloading only the on-disk JSON.
        base_features = {
            rid: list(feats) for rid, feats in self._features.items()
        }

        def _job(*, progress=None):
            def _prog(completed: int, total: int, label: str) -> None:
                if progress is None:
                    return
                if completed < 0:
                    # Nested tool chatter (e.g. "0/1 LTR") — prefix so it cannot
                    # steal the sequence-level N/M progress channel.
                    progress(f"detail:{label}")
                else:
                    progress(f"{max(0, int(completed))}/{max(1, int(total))} {label}")

            return run_te_annotation_on_fasta(
                path,
                options,
                auto_settings=auto_cfg,
                similarity_settings=sim_cfg,
                use_system_blast=load_use_system_blast_from_settings(),
                cancel=cancel,
                progress=_prog,
                record_ids=job_ids,
                persist=False,
                base_features_by_id=base_features,
            )

        worker = Worker(_job, progress=True)
        from repeatviewer.widgets.process_monitor_dialog import (
            attach_annotation_worker_process_monitor,
        )

        if kind == "similarity":
            cfg = self._similarity_ann_settings
            mon_title = "Similarity annotation"
            tool_key = "similarity_annotation"
        else:
            cfg = self._auto_ann_settings
            mon_title = "Feature detection"
            tool_key = "feature_detection"
        attach_annotation_worker_process_monitor(
            self,
            worker,
            title=mon_title,
            enabled=bool(getattr(cfg, "show_process_window", True)),
            close_when_finished=bool(
                getattr(cfg, "close_process_window_when_finished", True)
            ),
            tool_key=tool_key,
        )
        thread = QThread(self)
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        worker.signals.progress.connect(self._on_batch_annotation_progress)
        worker.signals.finished.connect(self._on_batch_annotation_finished)
        worker.signals.error.connect(self._on_batch_annotation_error)
        worker.signals.finished.connect(thread.quit)
        worker.signals.error.connect(thread.quit)
        thread.finished.connect(worker.deleteLater)
        thread.finished.connect(thread.deleteLater)
        thread.finished.connect(self._clear_annotation_thread)
        self._ann_worker = worker
        self._ann_thread = thread
        thread.start()

    @Slot(str)
    def _on_batch_annotation_progress(self, message: str) -> None:
        text = str(message or "").strip()
        if text.startswith("detail:"):
            # Nested detector updates — keep sequence progress stable (no rail flash).
            return
        match = re.match(r"^(\d+)\s*/\s*(\d+)\s+(.+)$", text)
        if match is None:
            self._set_annotation_job_status(text)
            return
        completed, total, detail = (
            int(match.group(1)),
            int(match.group(2)),
            match.group(3).strip(),
        )
        expected = max(1, int(getattr(self, "_ann_progress_total", total) or total))
        # Ignore nested mode counters (e.g. "0/1 LTR") that are not per-sequence.
        if total != expected:
            return
        kind = "Feature detection" if self._ann_job_kind == "auto" else "Similarity"
        status = f"{kind} {completed}/{total}: {detail}"
        if status != getattr(self, "_ann_last_sequence_status", ""):
            self._ann_last_sequence_status = status
            self._set_annotation_job_status(status)
        self._update_window_job_progress(completed, detail)

    @Slot(object)
    def _on_batch_annotation_finished(self, result: object) -> None:
        from repeatviewer.backends.te_annotation_batch import TeAnnotationBatchResult

        self._clear_annotation_job_status()
        kind = "Feature detection" if self._ann_job_kind == "auto" else "Similarity"
        if not isinstance(result, TeAnnotationBatchResult):
            self.statusMessage.emit(f"{kind} finished with an unexpected result.")
            return
        self._push_feature_undo_snapshot()
        # Apply in memory only — same as single-sequence until the user clicks Save.
        for rid in self._record_ids:
            self._features[rid] = list(result.features_by_id.get(rid) or [])
        self._annotations_dirty = True
        reveal_kinds = tuple(self._ann_reveal_kinds)
        self._refresh_features_after_edit()
        # Show the kinds that were just predicted even if Tracks had them off.
        if reveal_kinds:
            self._rail.ensure_annotation_kinds_visible(reveal_kinds)
            self._bars.set_visible_annotation_kinds(
                self._rail.visible_annotation_kinds()
            )
        skipped = list(getattr(result, "skipped_sequences", ()) or ())
        if skipped:
            headline, body = result.skipped_reminder_parts(kind=kind)
            box = QMessageBox(self)
            box.setIcon(QMessageBox.Icon.Warning)
            box.setWindowTitle(f"{kind} — some sequences skipped")
            # macOS: setText is the bold key line; informative text stays regular.
            box.setText(headline)
            box.setInformativeText(body)
            box.setStandardButtons(QMessageBox.StandardButton.Ok)
            box.exec()
            self.statusMessage.emit(
                f"{kind} finished with {len(skipped)} sequence(s) skipped. "
                "Click Save to write changes to disk."
            )
        else:
            self.statusMessage.emit(
                f"{kind} finished for selected sequences. "
                "Click Save to write changes to disk."
            )

    @Slot(str)
    def _on_batch_annotation_error(self, message: str) -> None:
        self._clear_annotation_job_status()
        lower = str(message or "").lower()
        if "stopped by user" in lower or "job stopped" in lower:
            self.statusMessage.emit("Annotation stopped.")
            return
        QMessageBox.warning(self, "Annotation failed", str(message))
        self.statusMessage.emit(f"Annotation failed: {message}")

    @Slot()
    def _clear_annotation_thread(self) -> None:
        self._ann_thread = None
        self._ann_worker = None
        self._ann_job_kind = ""
        self._ann_progress_total = 0
        self._ann_last_sequence_status = ""
        self._ann_reveal_kinds = ()
        self._set_annotation_controls_enabled(True)
        self._end_window_job_progress()

    def _push_feature_undo_snapshot(self) -> None:
        seqs = self._bars.msa_sequences() if self.is_msa_mode() else None
        lengths = dict(self._lengths) if seqs is not None else None
        ids = list(self._record_ids) if seqs is not None else None
        self._feature_undo_stack.append(
            _OverviewUndoSnapshot(
                features=_clone_features_by_id(self._features),
                sequences=list(seqs) if seqs is not None else None,
                lengths=lengths,
                record_ids=ids,
            )
        )
        if len(self._feature_undo_stack) > _FEATURE_UNDO_LIMIT:
            self._feature_undo_stack.pop(0)

    def _refresh_features_after_edit(self) -> None:
        """Reload canvas features while preserving MSA mode, zoom, and selection."""
        keep_ppb = self._bars.pixels_per_base()
        keep_sel = self._bars.selected_rows()
        keep_region = self._bars.region_selection()
        if self.is_msa_mode():
            self._bars.set_features_in_place(self._features)
        else:
            self._bars.set_data(
                self._record_ids,
                self._lengths,
                self._features,
                fit_bars_width=None,
            )
        self._bars.set_pixels_per_base(keep_ppb)
        self._bars.set_selected_rows(keep_sel)
        if keep_region is not None:
            self._bars.set_region_selection(keep_region[0], keep_region[1])
        self._sync_annotation_kind_filters()
        self._sync_names_column()
        self._update_scale()

    def _persist_overview_features(self) -> Path | None:
        try:
            _json, gff_path = save_annotations_multi(
                self._path,
                features_by_id=self._features,
                lengths_by_id=self._lengths,
                record_ids=list(self._record_ids),
            )
        except (OSError, ValueError) as exc:
            QMessageBox.critical(self, "Save failed", str(exc))
            return None
        if gff_path is not None:
            self.annotationsSaved.emit(gff_path)
        return gff_path

    def _on_delete_clicked(self) -> None:
        """Delete selection: features, MSA region/columns, or whole sequences."""
        # Shift+Delete fills with '-'; plain Delete splices/removes (MSA only).
        fill_gaps = bool(
            QGuiApplication.keyboardModifiers() & Qt.KeyboardModifier.ShiftModifier
        )
        # Feature picks also set a column span for Shift+arrows — prefer
        # deleting the selected annotations, not that whole column region.
        if self._bars.selected_features():
            self._delete_selected_features()
            return
        if self.is_msa_mode() and self._bars.region_selection() is not None:
            self._delete_selected_msa_edit(hard_delete=not fill_gaps)
            return
        # Whole-sequence selection (name column / row pick), MSA or multi-FASTA.
        if self._bars.selected_rows() and not self._bars.is_regional_selection():
            if self.is_msa_mode():
                self._delete_selected_msa_rows(hard_delete=not fill_gaps)
            else:
                self._delete_selected_multi_rows()
            return
        if self.is_msa_mode():
            QMessageBox.information(
                self,
                "Delete",
                "Select sequences, a region, or columns, then Delete (remove) "
                "or Shift+Delete (→ gaps). Ctrl+Z undoes. "
                "Or select a feature and Delete.",
            )
            return
        QMessageBox.information(
            self,
            "Delete",
            "Select sequence name(s) or a feature, then Delete.",
        )

    def _delete_selected_features(self) -> None:
        """Delete every selected annotation (not the column span they cover)."""
        selected = self._bars.selected_features()
        if not selected:
            return
        # Resolve via canvas packed tracks (respects type visibility filters).
        # Use store indices — not field identity — so stacked duplicates
        # (identical TSD flanks) are each removed when selected.
        indices_by_rid: dict[str, set[int]] = {}
        for row, feat_index in selected:
            resolved = self._bars._store_index_for_flat(row, feat_index)
            if resolved is None:
                continue
            rid, store_i, _feat = resolved
            indices_by_rid.setdefault(rid, set()).add(store_i)
        if not indices_by_rid:
            return

        self._push_feature_undo_snapshot()
        deleted = 0
        for rid, remove_at in indices_by_rid.items():
            store = list(self._features.get(rid) or [])
            kept = [feat for i, feat in enumerate(store) if i not in remove_at]
            deleted += len(store) - len(kept)
            self._features[rid] = kept

        self._bars.clear_feature_selection()
        # Clear the column span that feature selection mirrored onto the scale.
        if self._bars.region_selection() is not None:
            self._bars.set_region_selection(None)
        self._refresh_features_after_edit()
        self._persist_overview_features()
        if deleted == 1:
            rid = next(iter(indices_by_rid))
            self.statusMessage.emit(f"Deleted feature on {rid}.")
        else:
            self.statusMessage.emit(f"Deleted {deleted} features.")

    def _delete_selected_feature(self) -> None:
        """Back-compat alias for single/multi feature delete."""
        self._delete_selected_features()

    def _delete_selected_multi_rows(self) -> None:
        """Remove name-selected sequences from a multi-FASTA (non-MSA) file."""
        if self.is_msa_mode():
            return
        rows = sorted(self._bars.selected_rows())
        if not rows:
            return
        n = len(self._record_ids)
        row_set = {r for r in rows if 0 <= r < n}
        if not row_set:
            return
        if len(row_set) >= n:
            QMessageBox.warning(
                self,
                "Delete",
                "Cannot delete every sequence. Leave at least one.",
            )
            return
        try:
            _ids, sequences = load_fasta_rows(self._path)
        except (OSError, ValueError) as exc:
            QMessageBox.critical(self, "Delete", str(exc))
            return
        if len(sequences) != n:
            QMessageBox.warning(
                self,
                "Delete",
                "FASTA record count no longer matches the overview.",
            )
            return
        keep = [i for i in range(n) if i not in row_set]
        new_ids = [self._record_ids[i] for i in keep]
        new_seqs = [sequences[i] for i in keep]
        new_features = {
            rid: list(self._features.get(rid) or []) for rid in new_ids
        }
        new_lengths = {
            rid: len(seq) for rid, seq in zip(new_ids, new_seqs)
        }
        self._push_feature_undo_snapshot()
        try:
            write_fasta_rows(self._path, new_ids, new_seqs)
        except (OSError, ValueError) as exc:
            self._feature_undo_stack.pop()
            QMessageBox.critical(self, "Delete", str(exc))
            return
        self._record_ids = new_ids
        self._features = new_features
        self._lengths = new_lengths
        keep_ppb = self._bars.pixels_per_base()
        self._bars.set_data(
            new_ids,
            new_lengths,
            new_features,
            sequences=new_seqs,
            fit_bars_width=None,
        )
        self._bars.set_pixels_per_base(keep_ppb)
        self._bars.set_selected_rows([])
        self._bars.set_region_selection(None)
        self._sync_names_column()
        self._update_scale()
        self._sync_header_position_spins()
        self._persist_overview_features()
        self.statusMessage.emit(
            f"Removed {len(row_set)} sequence(s) · "
            f"{len(new_ids)} remain. Ctrl+Z to undo."
        )

    def _delete_selected_msa_rows(self, *, hard_delete: bool) -> None:
        """Edit name-selected whole sequences.

        * Delete — remove the selected sequence(s) from the MSA.
        * Shift+Delete — fill each selected row with gaps.
        """
        if not self.is_msa_mode():
            return
        rows = sorted(self._bars.selected_rows())
        if not rows:
            return
        seqs = self._bars.msa_sequences()
        if not seqs:
            return
        n = len(seqs)
        row_set = {r for r in rows if 0 <= r < n}
        if not row_set:
            return
        title = "Delete" if hard_delete else "Gap fill"
        if hard_delete and len(row_set) >= n:
            QMessageBox.warning(
                self,
                title,
                "Cannot delete every sequence. Leave at least one.",
            )
            return

        self._push_feature_undo_snapshot()
        if hard_delete:
            keep = [i for i in range(n) if i not in row_set]
            new_ids = [self._record_ids[i] for i in keep]
            new_seqs = [seqs[i] for i in keep]
            new_features = {
                rid: list(self._features.get(rid) or []) for rid in new_ids
            }
            new_width = len(new_seqs[0]) if new_seqs else 0
            if new_width <= 0 or not new_ids:
                self._feature_undo_stack.pop()
                QMessageBox.warning(
                    self,
                    title,
                    "Cannot delete every sequence. Leave at least one.",
                )
                return
            try:
                write_aligned_fasta_rows(self._path, new_ids, new_seqs)
            except (OSError, ValueError) as exc:
                self._feature_undo_stack.pop()
                QMessageBox.critical(self, title, str(exc))
                return
            self._record_ids = new_ids
            self._features = new_features
            self._lengths = {rid: new_width for rid in new_ids}
            keep_ppb = self._bars.pixels_per_base()
            self._bars.set_msa_data(
                new_ids,
                new_seqs,
                features=new_features,
                fit_bars_width=None,
            )
            self._bars.set_pixels_per_base(keep_ppb)
            self._bars.set_selected_rows([])
            self._bars.set_region_selection(None)
            # Re-pin if the standard survived deletion; else drop a stale save.
            self._restore_msa_diff_standard_from_settings()
            self._sync_names_column()
            self._update_scale()
            self._sync_header_position_spins()
            self._persist_overview_features()
            self.statusMessage.emit(
                f"Removed {len(row_set)} sequence(s) · "
                f"{len(new_ids)} remain. Ctrl+Z to undo."
            )
            return

        # Shift+Delete: gap-fill entire selected rows.
        width = len(seqs[0])
        try:
            new_seqs = mask_rows_region(seqs, 1, width, row_set)
        except ValueError as exc:
            self._feature_undo_stack.pop()
            QMessageBox.warning(self, title, str(exc))
            return
        edited_ids = {
            self._record_ids[r] for r in row_set if 0 <= r < len(self._record_ids)
        }
        new_features: dict[str, list[SequenceFeature]] = {}
        for rid in self._record_ids:
            feats = list(self._features.get(rid) or [])
            if rid in edited_ids:
                feats = drop_features_overlapping_region(feats, 1, width)
            new_features[rid] = feats
        try:
            write_aligned_fasta_rows(self._path, self._record_ids, new_seqs)
        except (OSError, ValueError) as exc:
            self._feature_undo_stack.pop()
            QMessageBox.critical(self, title, str(exc))
            return
        self._features = new_features
        keep_ppb = self._bars.pixels_per_base()
        keep_sel = list(row_set)
        self._bars.set_msa_data(
            self._record_ids,
            new_seqs,
            features=self._features,
            fit_bars_width=None,
        )
        self._bars.set_pixels_per_base(keep_ppb)
        self._bars.set_selected_rows(keep_sel)
        self._bars.set_region_selection(None)
        self._sync_names_column()
        self._update_scale()
        self._persist_overview_features()
        self.statusMessage.emit(
            f"Filled {len(row_set)} sequence(s) with gaps. Ctrl+Z to undo."
        )

    def _delete_selected_msa_edit(self, *, hard_delete: bool) -> None:
        """Edit the selected MSA region.

        * Delete — remove residues/columns (splice).
        * Shift+Delete — replace residues with ``-`` (keep columns).
        Features overlapping the edited span on changed rows are dropped.
        """
        region = self._bars.region_selection()
        if region is None or not self.is_msa_mode():
            return
        start, end = int(region[0]), int(region[1])
        if end < start:
            start, end = end, start
        seqs = self._bars.msa_sequences()
        if not seqs:
            return
        width = len(seqs[0])
        start = max(1, min(width, start))
        end = max(start, min(width, end))
        n_del = end - start + 1
        regional_rows = self._bars.regional_rows()
        is_regional = regional_rows is not None and len(regional_rows) > 0
        if is_regional:
            rows = set(regional_rows)
            edited_ids = {
                self._record_ids[r]
                for r in rows
                if 0 <= r < len(self._record_ids)
            }
        else:
            rows = set(range(len(seqs)))
            edited_ids = set(self._record_ids)

        title = "Delete" if hard_delete else "Gap fill"
        if hard_delete and n_del >= width and not is_regional:
            QMessageBox.warning(
                self,
                title,
                "Cannot delete the entire alignment. Leave at least one column.",
            )
            return
        if hard_delete and is_regional and n_del >= width:
            QMessageBox.warning(
                self,
                title,
                "Cannot delete the entire sequence. Leave at least one column.",
            )
            return

        self._push_feature_undo_snapshot()
        try:
            if hard_delete:
                if is_regional:
                    new_seqs = splice_rows_region(seqs, start, end, rows)
                else:
                    new_seqs = remove_alignment_columns(seqs, start, end)
            else:
                new_seqs = mask_rows_region(
                    seqs,
                    start,
                    end,
                    rows if is_regional else None,
                )
        except ValueError as exc:
            self._feature_undo_stack.pop()
            QMessageBox.warning(self, title, str(exc))
            return

        new_features: dict[str, list[SequenceFeature]] = {}
        for rid in self._record_ids:
            feats = list(self._features.get(rid) or [])
            if rid in edited_ids:
                # Drop annotations that overlap the edited span.
                feats = drop_features_overlapping_region(feats, start, end)
                if hard_delete:
                    feats = remap_features_after_column_delete(feats, start, end)
            new_features[rid] = feats

        new_width = len(new_seqs[0]) if new_seqs else 0
        if new_width <= 0:
            self._feature_undo_stack.pop()
            QMessageBox.warning(
                self,
                title,
                "Cannot delete the entire alignment. Leave at least one column.",
            )
            return
        try:
            write_aligned_fasta_rows(self._path, self._record_ids, new_seqs)
        except (OSError, ValueError) as exc:
            self._feature_undo_stack.pop()
            QMessageBox.critical(self, title, str(exc))
            return

        self._features = new_features
        self._lengths = {rid: new_width for rid in self._record_ids}
        keep_ppb = self._bars.pixels_per_base()
        self._bars.set_msa_data(
            self._record_ids,
            new_seqs,
            features=self._features,
            fit_bars_width=None,
        )
        self._bars.set_pixels_per_base(keep_ppb)
        self._bars.set_region_selection(None)
        self._sync_names_column()
        self._update_scale()
        self._sync_header_position_spins()
        self._persist_overview_features()
        if hard_delete:
            scope = (
                f"{len(rows)} sequence(s)"
                if is_regional
                else "all sequences"
            )
            self.statusMessage.emit(
                f"Removed columns {start}–{end} from {scope} · "
                f"{new_width} columns remain. Ctrl+Z to undo."
            )
        else:
            scope = (
                f"{len(rows)} sequence(s)"
                if is_regional
                else "all sequences"
            )
            self.statusMessage.emit(
                f"Filled columns {start}–{end} with gaps on {scope}. "
                "Ctrl+Z to undo."
            )

    def _restore_features(self, quiet: bool = False) -> None:
        if not self._feature_undo_stack:
            if not quiet:
                QMessageBox.information(
                    self,
                    "Restore",
                    "No changes to restore.",
                )
            return
        snap = self._feature_undo_stack.pop()
        self._features = snap.features
        self._bars.clear_feature_selection()
        if snap.sequences is not None and self.is_msa_mode():
            new_seqs = list(snap.sequences)
            new_ids = (
                list(snap.record_ids)
                if snap.record_ids is not None
                else list(self._record_ids)
            )
            if len(new_ids) != len(new_seqs):
                QMessageBox.critical(
                    self,
                    "Restore",
                    "Undo snapshot is inconsistent (ids vs sequences).",
                )
                return
            new_width = len(new_seqs[0]) if new_seqs else 1
            self._record_ids = new_ids
            self._lengths = (
                dict(snap.lengths)
                if snap.lengths is not None
                else {rid: new_width for rid in new_ids}
            )
            try:
                write_aligned_fasta_rows(self._path, new_ids, new_seqs)
            except (OSError, ValueError) as exc:
                QMessageBox.critical(self, "Restore", str(exc))
                return
            keep_ppb = self._bars.pixels_per_base()
            self._bars.set_msa_data(
                new_ids,
                new_seqs,
                features=self._features,
                fit_bars_width=None,
            )
            self._bars.set_pixels_per_base(keep_ppb)
            self._bars.set_selected_rows([])
            self._bars.set_region_selection(None)
            self._sync_names_column()
            self._update_scale()
            self._sync_header_position_spins()
            self._persist_overview_features()
            self.statusMessage.emit("Restored previous MSA and features.")
            return
        self._refresh_features_after_edit()
        self._persist_overview_features()
        self.statusMessage.emit("Restored previous feature state.")

    def _open_exported_file(self, path: Path) -> None:
        """Open a successful export in Preview (or the system viewer as fallback)."""
        resolved = path.expanduser().resolve()
        if not resolved.is_file():
            return
        window = self.window()
        opener = getattr(window, "_auto_open_outputs", None)
        if callable(opener):
            try:
                opener([resolved])
                return
            except Exception:
                pass
        from PySide6.QtCore import QUrl
        from PySide6.QtGui import QDesktopServices

        QDesktopServices.openUrl(QUrl.fromLocalFile(str(resolved)))

    def _export_overview_view(self) -> None:
        """Export the names + tracks overview (same formats as single-seq Export)."""
        from repeatviewer.views.export_annotation_view_dialog import (
            ExportAnnotationViewDialog,
        )

        dialog = ExportAnnotationViewDialog(
            source_path=self._path,
            title="Export overview",
            name_tag="overview",
            default_format="pdf",
            parent=self,
        )
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        out_path = dialog.chosen_path()
        if out_path is None:
            return
        dpi = max(50, min(600, int(dialog.dpi())))
        scale = max(1.0, float(dpi) / 96.0)

        # Tight figure: scale + sequence content (not the tall empty scroll viewport).
        scale_h = max(1, int(self._scale.height() or self._scale.sizeHint().height()))
        body_h = max(1, int(self._bars.content_height()))
        names_w = max(1, int(self._names_scroll.width() or self._name_width))
        tracks_w = max(1, int(self._bars_scroll.width() or 200))
        handle_w = max(0, int(self._splitter.handleWidth()))
        width = max(1, names_w + handle_w + tracks_w)
        height = scale_h + body_h

        def _paint_overview(painter: QPainter) -> None:
            from PySide6.QtGui import QRegion

            # Names column: scale shim + names (clipped to content height).
            self._scale_shim.render(painter, QPoint(0, 0))
            self._names.render(
                painter,
                QPoint(0, scale_h),
                QRegion(QRect(0, 0, names_w, body_h)),
            )
            x = names_w + handle_w
            self._scale.render(painter, QPoint(x, 0))
            self._bars.render(
                painter,
                QPoint(x, scale_h),
                QRegion(QRect(0, 0, tracks_w, body_h)),
            )

        suffix = out_path.suffix.lower()
        if suffix == ".svg":
            try:
                from PySide6.QtSvg import QSvgGenerator  # type: ignore
            except Exception:  # noqa: BLE001
                QMessageBox.warning(
                    self,
                    "Export",
                    "SVG export needs PySide6 QtSvg module. Install it and retry.",
                )
                return
            gen = QSvgGenerator()
            gen.setFileName(str(out_path))
            gen.setSize(QSize(width, height))
            gen.setViewBox(QRect(0, 0, width, height))
            gen.setResolution(dpi)
            painter = QPainter(gen)
            try:
                _paint_overview(painter)
            finally:
                painter.end()
            self.statusMessage.emit(f"Exported overview → {out_path.name}")
            self._open_exported_file(out_path)
            return

        pixel = QSize(
            max(1, int(round(width * scale))),
            max(1, int(round(height * scale))),
        )
        img = QImage(pixel, QImage.Format.Format_ARGB32)
        img.fill(Qt.GlobalColor.transparent)
        painter = QPainter(img)
        try:
            painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
            painter.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform, True)
            painter.scale(scale, scale)
            _paint_overview(painter)
        finally:
            painter.end()

        if suffix == ".png":
            if not img.save(str(out_path), "PNG"):
                QMessageBox.warning(self, "Export", "Could not save PNG.")
                return
            self.statusMessage.emit(f"Exported overview → {out_path.name}")
            self._open_exported_file(out_path)
            return

        if suffix == ".pdf":
            from repeatviewer.views.export_annotation_view_dialog import (
                write_view_image_pdf,
            )

            write_view_image_pdf(
                out_path,
                img,
                dpi=dpi,
                content_width=width,
                content_height=height,
            )
            self.statusMessage.emit(f"Exported overview → {out_path.name}")
            self._open_exported_file(out_path)
            return

        QMessageBox.warning(
            self,
            "Export",
            f"Unsupported export format: {suffix!r}. Use .svg, .png, or .pdf.",
        )

    def is_dirty(self) -> bool:
        """True when feature detection / similarity results are not yet saved."""
        return bool(self._annotations_dirty)

    def _save_overview_annotations(self) -> None:
        if not self._record_ids:
            return
        try:
            json_path, _gff = save_annotations_multi(
                self._path,
                features_by_id=self._features,
                lengths_by_id=self._lengths,
                record_ids=self._record_ids,
            )
        except (OSError, ValueError) as exc:
            QMessageBox.warning(self, "Save annotations", str(exc))
            return
        self._annotations_dirty = False
        self.statusMessage.emit(f"Saved annotations · {json_path.name}")
        self.annotationsSaved.emit(self._path)

    def load_path(self, path: Path, *, msa_mode: bool | None = None) -> None:
        """Load multi-FASTA overview.

        When ``msa_mode`` is True (or auto-detected via equal lengths), paint
        aligned residue rows instead of length bars.
        """
        self._path = path.expanduser().resolve()
        self.clear_find_hits()  # hits belong to the previous file
        self._fasta_seq_cache = None
        use_msa = bool(msa_mode) if msa_mode is not None else is_aligned_fasta(self._path)
        if use_msa:
            self._load_msa_path()
            return
        try:
            record_ids, sequences = load_fasta_rows(self._path)
            lengths = {
                rid: len(seq) for rid, seq in zip(record_ids, sequences)
            }
        except (OSError, ValueError) as exc:
            self._record_ids = []
            self._lengths = {}
            self._features = {}
            self._annotations_dirty = False
            self._bars.set_data([], {}, {})
            self._sync_mode_specific_controls()
            self._scale.set_scale(ppb=_ABS_MIN_PPB, max_length=1, bar_pad=8)
            self.statusMessage.emit(str(exc))
            return

        # Prefer JSON sidecar (canonical store); fall back to companion GFF3.
        json_by_id = load_features_by_id_from_json(self._path, record_ids)
        gff = resolve_features_gff_path(self._path)
        sidecar = resolve_sidecar_path(self._path)
        if json_by_id and any(json_by_id.values()):
            features = {rid: list(json_by_id.get(rid) or []) for rid in record_ids}
            source_label = sidecar.name if sidecar is not None else "JSON"
        elif gff is not None:
            features = load_all_features_by_id(
                gff, record_ids, lengths_by_id=lengths
            )
            if json_by_id:
                for rid, json_feats in json_by_id.items():
                    if not json_feats:
                        continue
                    current = features.get(rid) or []
                    features[rid] = json_feats if not current else current
            source_label = gff.name
        else:
            features = {rid: [] for rid in record_ids}
            source_label = None

        self._record_ids = record_ids
        self._lengths = lengths
        self._features = features
        self._annotations_dirty = False
        self._did_initial_fit = False
        self._search_needle = ""
        self._search_hits = []
        self._search_hit_index = -1
        total_w = max(400, self.width() or 800)
        self._name_width = self._compute_name_width(total_w)
        tracks_w = max(200, total_w - self._name_width)
        self._splitter.setSizes([self._name_width, tracks_w])
        self._sync_viewport_fill()
        self._bars.set_data(
            record_ids,
            lengths,
            features,
            sequences=sequences,
            fit_bars_width=None,
        )
        self._sync_mode_specific_controls()
        self._bars.sort_rows(self._rail.sequence_sort_mode())
        self._record_ids = list(self._bars.record_ids())
        self._sync_annotation_kind_filters()
        self._sync_names_column()
        self._update_scale()
        self._sync_header_position_spins()
        # New file: always start at top-left (do not keep prior zoom/scroll).
        self._reset_overview_to_start()
        QTimer.singleShot(0, self._try_initial_fit)
        QTimer.singleShot(0, self._reset_scrollbars_to_origin)
        QTimer.singleShot(50, self._try_initial_fit)
        QTimer.singleShot(200, self._try_initial_fit)

        annotated = sum(1 for rid in record_ids if features.get(rid))
        feat_total = sum(len(features.get(rid) or []) for rid in record_ids)
        if source_label is not None and feat_total:
            self.statusMessage.emit(
                f"{len(record_ids):,} sequences · {annotated:,} with annotations · "
                f"{feat_total:,} features · {source_label}"
            )
        elif source_label is not None:
            self.statusMessage.emit(
                f"{len(record_ids):,} sequences · {source_label} (no features yet)"
            )
        else:
            self.statusMessage.emit(
                f"{len(record_ids):,} sequences · no annotation JSON / .RV.gff3 yet"
            )

    def _load_msa_path(self) -> None:
        """Load equal-length FASTA as colored residue rows with annotations."""
        try:
            record_ids, sequences = load_aligned_fasta_rows(self._path)
        except (OSError, ValueError) as exc:
            self._record_ids = []
            self._lengths = {}
            self._features = {}
            self._annotations_dirty = False
            self._bars.set_data([], {}, {})
            self._sync_mode_specific_controls()
            self._scale.set_scale(ppb=_ABS_MIN_PPB, max_length=1, bar_pad=8)
            self.statusMessage.emit(str(exc))
            return
        width = len(sequences[0]) if sequences else 1
        self._record_ids = record_ids
        self._lengths = {rid: width for rid in record_ids}
        # Prefer JSON sidecar only when it actually has features; an empty
        # sidecar must not hide companion GFF annotations (same as non-MSA).
        json_by_id = load_features_by_id_from_json(self._path, record_ids)
        gff = resolve_features_gff_path(self._path)
        if json_by_id and any(json_by_id.values()):
            self._features = {
                rid: list(json_by_id.get(rid) or []) for rid in record_ids
            }
        elif gff is not None:
            self._features = load_all_features_by_id(
                gff, record_ids, lengths_by_id=self._lengths
            )
            if json_by_id:
                for rid, json_feats in json_by_id.items():
                    if not json_feats:
                        continue
                    current = self._features.get(rid) or []
                    self._features[rid] = json_feats if not current else current
        else:
            self._features = {rid: [] for rid in record_ids}
        self._annotations_dirty = False
        self._did_initial_fit = False
        self._search_needle = ""
        self._search_hits = []
        self._search_hit_index = -1
        total_w = max(400, self.width() or 800)
        self._name_width = self._compute_name_width(total_w)
        tracks_w = max(200, total_w - self._name_width)
        self._splitter.setSizes([self._name_width, tracks_w])
        self._sync_viewport_fill()
        self._bars.set_msa_data(
            record_ids,
            sequences,
            features=self._features,
            fit_bars_width=None,
        )
        self._sync_mode_specific_controls()
        # Apply saved sort after rows are loaded (original order until then).
        self._bars.sort_rows(self._rail.sequence_sort_mode())
        self._record_ids = list(self._bars.record_ids())
        self._restore_msa_diff_standard_from_settings()
        self._sync_annotation_kind_filters()
        self._sync_names_column()
        self._update_scale()
        self._sync_header_position_spins()
        # New MSA: always start at top-left (do not keep prior zoom/scroll).
        self._reset_overview_to_start()
        QTimer.singleShot(0, self._try_initial_fit)
        QTimer.singleShot(0, self._reset_scrollbars_to_origin)
        QTimer.singleShot(50, self._try_initial_fit)
        QTimer.singleShot(200, self._try_initial_fit)
        kind = "protein" if self._bars._msa_protein else "nucleotide"
        n_feat = sum(len(v) for v in self._features.values())
        self.statusMessage.emit(
            f"{len(record_ids):,} aligned sequences · {width:,} columns · "
            f"{kind} MSA · {n_feat:,} features"
        )

    def _sync_viewport_fill(self) -> None:
        if not self.isVisible():
            return
        vp = self._bars_scroll.viewport()
        if vp.width() < 80 or vp.height() < 40:
            return
        self._bars.set_viewport_min_size(vp.width(), vp.height())

    def _try_initial_fit(self) -> None:
        if not self._overview_panel_alive():
            return
        if self._did_initial_fit or not self._record_ids:
            return
        if not self.isVisible():
            self._refit_pending_on_show = True
            return
        if self._bars_scroll.viewport().width() < 120:
            return
        self._fit_bars()
        self._did_initial_fit = True

    def _after_panel_shown(self) -> None:
        """Restore full MSA layout after the panel was hidden (e.g. PDF tab)."""
        # Tab replace can destroy this panel while a showEvent singleShot is queued.
        if not self._overview_panel_alive():
            return
        if not self.isVisible() or not self._record_ids:
            return
        self._bars._invalidate_msa_layer_cache()
        self._refit_pending_on_show = False
        if not self._did_initial_fit:
            self._try_initial_fit()
            return
        # Always re-measure the scroll viewport — it may have grown while
        # this panel was covered by another Preview page.
        self._refit_after_viewport_change()
        self._bars.update()

    def _bars_viewport_width(self) -> int:
        return max(80, self._bars_scroll.viewport().width())

    def _sync_header_position_spins(self) -> None:
        bar = self._header_zoom_bar
        if bar is None or not hasattr(bar, "set_position_range"):
            return
        maximum = max(1, int(self._bars.max_length()))
        region = self._bars.region_selection()
        if region is not None:
            start, end = int(region[0]), int(region[1])
            bar.set_position_range(maximum, start=start, end=end)
        else:
            bar.set_position_range(maximum, start=1, end=maximum)

    def _on_position_spin_changed(self, _value: int = 0) -> None:
        """Keep the MSA selection band in sync when Position spins change."""
        bar = self._header_zoom_bar
        if bar is None:
            return
        maximum = max(1, int(self._bars.max_length()))
        start = int(bar.pos_start.value())
        end = int(bar.pos_end.value())
        if end < start:
            start, end = end, start
        start = max(1, min(maximum, start))
        end = max(start, min(maximum, end))
        self._bars.set_region_selection(start, end)

    def _select_by_position(self) -> None:
        """Zoom/scroll to the current region or Position spins (Select).

        Regional (columns × rows) selections keep their row set and are
        centered both horizontally and vertically after zoom.
        """
        bar = self._header_zoom_bar
        if bar is None:
            return
        maximum = max(1, int(self._bars.max_length()))
        # Prefer the blue MSA band when present — spins can lag or mismatch.
        region = self._bars.region_selection()
        regional_rows = self._bars.regional_rows()
        if region is not None:
            start, end = int(region[0]), int(region[1])
        else:
            start = int(bar.pos_start.value())
            end = int(bar.pos_end.value())
        if end < start:
            start, end = end, start
        start = max(1, min(maximum, start))
        end = max(start, min(maximum, end))
        bar.set_position_range(maximum, start=start, end=end)
        # Keep a regional matrix box; only promote spins → full-height columns
        # when there was no row-limited selection.
        if regional_rows is not None:
            self._bars.set_regional_selection(start, end, regional_rows)
        else:
            self._bars.set_region_selection(start, end)
        # Inclusive columns → half-open window for zoom_to_bp_range.
        # Point features (start==end) use a 1-column window (span floored at 1).
        self._on_region_zoom(float(start - 1), float(end))
        self._sync_zoom_slider_from_ppb()
        # After zoom, scroll so the band is centered (pending scroll may be start-aligned).
        ppb = max(self._bars.pixels_per_base(), 1e-12)
        mid_x = self._bars.bar_pad() + ((start - 1 + end) / 2.0) * ppb
        vp_w = max(1, self._bars_scroll.viewport().width())
        target = int(round(mid_x - vp_w / 2.0))
        hbar = self._bars_scroll.horizontalScrollBar()
        hbar.setValue(max(0, min(hbar.maximum(), target)))
        # Center regional rows / selected feature in the viewport.
        self._center_selection_vertically(rows=regional_rows)
        self._update_scale()
        self._sync_names_column()

    def _reset_scrollbars_to_origin(self) -> None:
        """Pin Overview to bp 1 and the first sequence row."""
        self._bars_scroll.horizontalScrollBar().setValue(0)
        self._bars_scroll.verticalScrollBar().setValue(0)
        self._names_scroll.verticalScrollBar().setValue(0)

    def _reset_overview_to_start(self) -> None:
        """After loading a different FASTA/MSA, fit and show the top-left.

        Previous file zoom (``_ppb``) and scrollbar positions must not carry
        over when switching alignments in the same Overview panel.
        """
        if not self._overview_panel_alive():
            return
        self._name_width = self._current_name_width()
        self._sync_viewport_fill()
        # Drop any pending mid-alignment scroll from the previous file.
        self._bars.fit_to_bars_width(self._bars_viewport_width())
        self._reset_scrollbars_to_origin()
        self._sync_names_column()
        self._update_scale()
        self._sync_zoom_slider_from_ppb()
        self._sync_header_position_spins()

    def _fit_bars(self) -> None:
        # Keep the user's splitter ratio; only zoom tracks to the track viewport.
        if not self._overview_panel_alive():
            return
        self._reset_overview_to_start()

    def _open_msa_clean_dialog(self) -> None:
        """Open MSA clean parameters for the alignment shown in this overview."""
        from repeatviewer.job_options import (
            load_msa_clean_options,
            save_msa_clean_options,
        )
        from repeatviewer.views.tool_dialogs import MsaCleanParamsDialog

        if not self._path.is_file():
            QMessageBox.warning(
                self,
                "MSA Clean",
                f"Sequence file not found:\n{self._path}",
            )
            return
        opts = load_msa_clean_options()
        dialog = MsaCleanParamsDialog(initial=opts, parent=self)
        if not dialog.exec():
            return
        opts = dialog.options()
        save_msa_clean_options(opts)
        if not dialog.should_run():
            self.statusMessage.emit("MSA clean settings saved.")
            return
        self._run_msa_clean(dialog.selected_method(), opts)

    def _run_msa_clean(self, method: str, opts=None) -> None:
        """Clean the open MSA and open the result as an additional Preview tab."""
        from repeatviewer.backends.curation import run_on_sequence
        from repeatviewer.job_options import MsaCleanOptions, load_msa_clean_options

        msa_opts = opts if isinstance(opts, MsaCleanOptions) else load_msa_clean_options()
        method_id = str(method or "msa_clean").strip() or "msa_clean"
        project = None
        win = self.window()
        if win is not None:
            project = getattr(win, "_project", None)
        if project is None:
            QMessageBox.warning(
                self,
                "MSA Clean",
                "Could not access the project to run MSA clean.",
            )
            return
        self.statusMessage.emit(f"Cleaning {self._path.name}…")
        try:
            created = run_on_sequence(
                method_id,
                self._path,
                project,
                msa_clean_options=msa_opts,
            )
        except Exception as exc:  # noqa: BLE001 — surface tool failures
            QMessageBox.critical(self, "MSA Clean failed", str(exc))
            return
        paths = [
            Path(p).expanduser().resolve()
            for p in (created or [])
            if Path(p).is_file()
        ]
        if not paths:
            QMessageBox.warning(
                self,
                "MSA Clean",
                "MSA clean finished but no output file was created.",
            )
            return
        for path in paths:
            self.sequenceExtracted.emit(path)
            self.openInAdditionalPreviewTab.emit(path)
        names = ", ".join(p.name for p in paths[:3])
        if len(paths) > 3:
            names += ", …"
        self.statusMessage.emit(f"MSA clean finished · {names}")

    def _open_consensus_dialog(self) -> None:
        """Open Consensus parameters for the FASTA shown in this overview."""
        from repeatviewer.job_options import (
            load_consensus_options,
            save_consensus_options,
        )
        from repeatviewer.views.tool_dialogs import ConsensusParamsDialog

        if not self._path.is_file():
            QMessageBox.warning(
                self,
                "Consensus",
                f"Sequence file not found:\n{self._path}",
            )
            return
        opts = load_consensus_options()
        dialog = ConsensusParamsDialog(initial=opts, parent=self)
        if not dialog.exec():
            return
        opts = dialog.options()
        save_consensus_options(opts)
        if not dialog.should_run():
            self.statusMessage.emit("Consensus settings saved.")
            return
        self._run_consensus(opts)

    def _run_consensus(self, opts=None) -> None:
        """Write a consensus FASTA beside the current alignment."""
        from repeatviewer.backends.teaid_plotter import TeaidInputError, con_generater
        from repeatviewer.job_options import ConsensusOptions, load_consensus_options

        cons_opts = opts if opts is not None else load_consensus_options()
        if not isinstance(cons_opts, ConsensusOptions):
            cons_opts = load_consensus_options()
        try:
            out_path, length = con_generater(
                str(self._path),
                str(self._path.parent),
                threshold=cons_opts.threshold,
                ambiguous=cons_opts.ambiguous,
                sequence_type=cons_opts.sequence_type,
            )
        except TeaidInputError as exc:
            QMessageBox.warning(self, "Consensus", str(exc))
            return
        except Exception as exc:  # noqa: BLE001 — surface tool failures
            QMessageBox.critical(self, "Consensus failed", str(exc))
            return
        result = Path(out_path).expanduser().resolve()
        self.sequenceExtracted.emit(result)
        self.statusMessage.emit(
            f"Consensus written · {result.name} ({length:,} residues)"
        )

    def _slider_ppb_bounds(self) -> tuple[float, float]:
        return self._bars.zoom_bounds(self._bars_viewport_width())

    def _letter_anchor_ppb(self) -> float:
        """ppb at DNA/protein size 10 — zoom-slider 3/4 position."""
        fit_ppb, max_ppb = self._slider_ppb_bounds()
        _band, target = self._bars.letter_size_anchor_metrics()
        return max(fit_ppb, min(max_ppb, float(target)))

    @staticmethod
    def _interp_ppb(lo: float, hi: float, u: float) -> float:
        u = max(0.0, min(1.0, float(u)))
        if hi <= lo + 1e-15:
            return lo
        if hi / lo > 50:
            return math.exp(math.log(lo) + u * (math.log(hi) - math.log(lo)))
        return lo + (hi - lo) * u

    @staticmethod
    def _ppb_unit(lo: float, hi: float, ppb: float) -> float:
        ppb = max(lo, min(hi, float(ppb)))
        if hi <= lo + 1e-15:
            return 0.0
        if hi / lo > 50:
            return (math.log(ppb) - math.log(lo)) / (math.log(hi) - math.log(lo))
        return (ppb - lo) / (hi - lo)

    def _slider_value_to_ppb(self, value: int) -> float:
        """Map slider 0→fit, 75→DNA/protein size 10, 100→max (wider nucleotide gaps)."""
        fit_ppb, max_ppb = self._slider_ppb_bounds()
        if max_ppb <= fit_ppb + 1e-15:
            return fit_ppb
        letter = self._letter_anchor_ppb()
        t = max(0, min(100, int(value))) / 100.0
        pivot = _MSA_ZOOM_LETTER_SLIDER / 100.0
        if t <= pivot:
            return self._interp_ppb(fit_ppb, letter, t / pivot)
        return self._interp_ppb(letter, max_ppb, (t - pivot) / (1.0 - pivot))

    def _ppb_to_slider_value(self, ppb: float) -> int:
        fit_ppb, max_ppb = self._slider_ppb_bounds()
        ppb = max(fit_ppb, min(max_ppb, float(ppb)))
        if max_ppb <= fit_ppb + 1e-15:
            return 0
        letter = self._letter_anchor_ppb()
        pivot = _MSA_ZOOM_LETTER_SLIDER / 100.0
        if ppb <= letter + 1e-12:
            u = self._ppb_unit(fit_ppb, letter, ppb)
            return int(round(max(0.0, min(1.0, u)) * _MSA_ZOOM_LETTER_SLIDER))
        u = self._ppb_unit(letter, max_ppb, ppb)
        return int(
            round(
                _MSA_ZOOM_LETTER_SLIDER
                + max(0.0, min(1.0, u)) * (100 - _MSA_ZOOM_LETTER_SLIDER)
            )
        )

    def _sync_zoom_slider_from_ppb(self) -> None:
        if self._updating_zoom_slider:
            return
        self._updating_zoom_slider = True
        try:
            self._zoom_slider.setValue(
                self._ppb_to_slider_value(self._bars.pixels_per_base())
            )
        finally:
            self._updating_zoom_slider = False

    def _nudge_zoom_slider(self, delta: int) -> None:
        self._zoom_slider.setValue(
            max(0, min(100, self._zoom_slider.value() + int(delta)))
        )

    def _selection_zoom_anchor(self) -> tuple[float, float] | None:
        """``(mid_bp_0based, viewport_x)`` for the selected region/feature center.

        Uses the column band when present; otherwise the selected feature span
        (including point features with start==end).
        """
        region = self._bars.region_selection()
        if region is None:
            span = self._bars.selected_feature_column_span()
            if span is None:
                return None
            start, end = int(span[0]), int(span[1])
        else:
            start, end = int(region[0]), int(region[1])
        if end < start:
            start, end = end, start
        # Single column: mid of that base; multi: mid of the inclusive span.
        mid_bp = float(start - 1 + end) / 2.0
        mid_content_x = self._bars.bp_to_x(mid_bp)
        hbar = self._bars_scroll.horizontalScrollBar()
        view_x = float(mid_content_x) - float(hbar.value())
        return mid_bp, view_x

    def _scroll_to_keep_bp_at_viewport_x(self, bp: float, viewport_x: float) -> None:
        """Scroll so sequence coordinate ``bp`` (0-based) sits at ``viewport_x``."""
        new_x = self._bars.bp_to_x(float(bp))
        hbar = self._bars_scroll.horizontalScrollBar()
        target = int(round(new_x - float(viewport_x)))
        hbar.setValue(max(0, min(hbar.maximum(), target)))

    def _center_selection_vertically(
        self, *, rows: set[int] | list[int] | None = None
    ) -> None:
        """Scroll so the selected feature or regional rows sit mid-viewport."""
        tops = self._bars.row_tops()
        heights = self._bars.row_heights()
        vbar = self._bars_scroll.verticalScrollBar()
        vp_h = max(1, self._bars_scroll.viewport().height())

        def _center_y(mid_y: float) -> None:
            target = int(round(float(mid_y) - vp_h / 2.0))
            vbar.setValue(max(0, min(vbar.maximum(), target)))

        feat = self._bars.selected_feature()
        if feat is not None and rows is None:
            row, flat_i = feat
            rect = self._bars._feature_rect(row, flat_i)
            if rect is not None:
                _center_y(rect.center().y())
                return
            if 0 <= row < len(tops):
                _center_y(float(tops[row]) + float(heights[row]) * 0.5)
            return
        regional = (
            {int(r) for r in rows}
            if rows is not None
            else self._bars.regional_rows()
        )
        if not regional:
            return
        row_list = sorted(r for r in regional if 0 <= int(r) < len(tops))
        if not row_list:
            return
        y0 = float(tops[row_list[0]])
        y1 = float(
            self._bars._sequence_band_bottom(row_list[-1])
        )
        _center_y(0.5 * (y0 + y1))

    def _on_zoom_slider_changed(self, value: int) -> None:
        if self._updating_zoom_slider:
            return
        fit_ppb, _ = self._slider_ppb_bounds()
        ppb = self._slider_value_to_ppb(value)
        reset_scroll = value <= 0 or ppb <= fit_ppb * 1.001
        anchor = None if reset_scroll else self._selection_zoom_anchor()
        self._sync_viewport_fill()
        self._bars.set_pixels_per_base(ppb)
        if reset_scroll:
            self._bars_scroll.horizontalScrollBar().setValue(0)
        elif anchor is not None:
            mid_bp, view_x = anchor
            self._scroll_to_keep_bp_at_viewport_x(mid_bp, view_x)
            self._center_selection_vertically()
        self._sync_names_column()
        self._update_scale()

    def _scroll_to_row(self, row: int) -> None:
        tops = self._bars.row_tops()
        heights = self._bars.row_heights()
        if not (0 <= row < len(tops)):
            return
        y = tops[row]
        h = heights[row]
        vbar = self._bars_scroll.verticalScrollBar()
        vp_h = max(1, self._bars_scroll.viewport().height())
        if y < vbar.value():
            vbar.setValue(y)
        elif y + h > vbar.value() + vp_h:
            vbar.setValue(max(0, y + h - vp_h))

    def _on_name_filter_edited(self, text: str) -> None:
        needle = text.strip()
        if not needle:
            self._search_needle = ""
            self._search_hits = []
            self._search_hit_index = -1
            return
        self._jump_to_name_filter(needle, advance=False)

    def _on_name_filter_next(self) -> None:
        needle = self._name_filter.text().strip()
        if not needle:
            return
        self._jump_to_name_filter(needle, advance=True)

    def _jump_to_name_filter(self, needle: str, *, advance: bool) -> None:
        lowered = needle.lower()
        hits = [
            i
            for i, rid in enumerate(self._record_ids)
            if lowered in rid.lower()
        ]
        if not hits:
            self._search_needle = needle
            self._search_hits = []
            self._search_hit_index = -1
            self.statusMessage.emit(f"No sequence name contains “{needle}”")
            return

        if (
            advance
            and needle.lower() == self._search_needle.lower()
            and hits == self._search_hits
            and self._search_hit_index >= 0
        ):
            self._search_hit_index = (self._search_hit_index + 1) % len(hits)
        else:
            self._search_needle = needle
            self._search_hits = hits
            self._search_hit_index = 0

        row = hits[self._search_hit_index]
        self._bars.set_selected_rows([row])
        self._scroll_to_row(row)
        self.statusMessage.emit(
            f"Filter “{needle}”: {self._search_hit_index + 1}/{len(hits)} — "
            f"{self._record_ids[row]}"
        )

    def _zoom_step(self, factor: float) -> None:
        """Zoom from −/+ controls; prefer selected feature/region center."""
        vp = self._bars_scroll.viewport()
        hbar = self._bars_scroll.horizontalScrollBar()
        anchor = self._selection_zoom_anchor()
        if anchor is not None:
            mid_bp, _view_x = anchor
            # Content x of selection mid at the current zoom (before zoom_by).
            anchor_x = self._bars.bp_to_x(float(mid_bp))
        else:
            anchor_x = hbar.value() + vp.width() / 2
        self._bars.zoom_by(factor, anchor_x=anchor_x)
        if anchor is not None:
            self._center_selection_vertically()

    def _on_scale_zoom_at(self, content_x: float, zoom_in: bool) -> None:
        """Scale double-click: DNA/protein size 10 at slider 3/4; right-click → Fit.

        Zooming further past 3/4 (slider / −+) only increases horizontal
        nucleotide spacing (``ppb``), keeping the size-10 letter floor.
        """
        slider = self._zoom_slider
        old_ppb = max(self._bars.pixels_per_base(), 1e-12)
        # Double-click → size-10 / 3/4 slider; right-click → fit.
        new_val = _MSA_ZOOM_LETTER_SLIDER if zoom_in else 0
        if zoom_in:
            self._bars.set_msa_residue_size(_MSA_LETTER_SIZE_PT)
        if slider is not None:
            if int(slider.value()) == new_val:
                # Already at target zoom — still re-center on the click.
                anchor_bp = (float(content_x) - self._bars.bar_pad()) / old_ppb
                new_x = self._bars.bp_to_x(anchor_bp)
                hbar = self._bars_scroll.horizontalScrollBar()
                view_x = float(content_x) - hbar.value()
                hbar.setValue(max(0, int(new_x - view_x)))
                return
            new_ppb = self._slider_value_to_ppb(new_val)
        else:
            fit_ppb, _max_ppb = self._bars.zoom_bounds(self._bars_viewport_width())
            new_ppb = self._letter_anchor_ppb() if zoom_in else fit_ppb
        factor = float(new_ppb) / old_ppb
        self._bars.zoom_by(factor, anchor_x=float(content_x))

    def _update_scale(self) -> None:
        self._scale.set_scale(
            ppb=self._bars.pixels_per_base(),
            max_length=self._bars.max_length(),
            bar_pad=self._bars.bar_pad(),
        )
        self._scale.set_scroll_x(self._bars_scroll.horizontalScrollBar().value())

    def _on_hscroll(self, value: int) -> None:
        self._scale.set_scroll_x(value)

    def _on_names_vscroll(self, value: int) -> None:
        # Names bar is AlwaysOff; ignore user-driven changes (should not happen).
        if self._syncing_scroll:
            return
        self._syncing_scroll = True
        self._bars_scroll.verticalScrollBar().setValue(value)
        self._syncing_scroll = False

    def _on_bars_vscroll(self, value: int) -> None:
        if self._syncing_scroll:
            return
        self._syncing_scroll = True
        self._names_scroll.verticalScrollBar().setValue(value)
        self._syncing_scroll = False

    def _nudge_vertical_scroll(self, delta_y: int) -> None:
        if not delta_y:
            return
        bar = self._bars_scroll.verticalScrollBar()
        bar.setValue(bar.value() - int(delta_y))

    def _nudge_horizontal_scroll(self, delta_x: int) -> None:
        if not delta_x:
            return
        bar = self._bars_scroll.horizontalScrollBar()
        if bar.maximum() <= 0:
            return
        bar.setValue(bar.value() - int(delta_x))

    def eventFilter(self, obj, event) -> bool:  # noqa: N802
        bars_vp = self._bars_scroll.viewport()
        # Track viewport width changes from rail collapse (parent size unchanged).
        if obj is bars_vp and event.type() == QEvent.Type.Resize:
            self._schedule_viewport_refit()
            return super().eventFilter(obj, event)
        if event.type() != QEvent.Type.Wheel or not isinstance(event, QWheelEvent):
            return super().eventFilter(obj, event)
        names_vp = self._names_scroll.viewport()
        if obj not in (names_vp, bars_vp):
            return super().eventFilter(obj, event)
        mods = event.modifiers() | QGuiApplication.keyboardModifiers()
        if mods & Qt.KeyboardModifier.ShiftModifier:
            # Handle zoom here — wheel is delivered to the viewport, not the canvas.
            delta = event.angleDelta().y()
            if delta == 0:
                delta = event.angleDelta().x()
            pixel = event.pixelDelta()
            if delta == 0 and not pixel.isNull():
                delta = int(pixel.y()) or int(pixel.x())
            if delta != 0:
                factor = _ZOOM_STEP if delta > 0 else 1.0 / _ZOOM_STEP
                # Prefer selected-region center (single-seq parity); else cursor.
                sel_anchor = self._selection_zoom_anchor()
                if sel_anchor is not None:
                    mid_bp, _view_x = sel_anchor
                    anchor_x = self._bars.bp_to_x(float(mid_bp))
                else:
                    anchor_x = float(event.position().x()) + float(
                        self._bars_scroll.horizontalScrollBar().value()
                    )
                self._bars.zoom_by(factor, anchor_x=anchor_x)
                if sel_anchor is not None:
                    self._center_selection_vertically()
            event.accept()
            return True
        # Trackpad: apply both axes. Prefer the dominant axis when both move.
        dx = _wheel_delta_x(
            event, scroll_step=self._bars_scroll.horizontalScrollBar().singleStep() * 3
        )
        dy = _wheel_delta_y(event)
        if dx and dy and abs(dx) >= abs(dy):
            dy = 0
        elif dx and dy and abs(dy) > abs(dx):
            dx = 0
        handled = False
        if dx:
            self._nudge_horizontal_scroll(dx)
            handled = True
        if dy:
            self._nudge_vertical_scroll(dy)
            handled = True
        if handled:
            event.accept()
            return True
        return False

    def _on_zoom_changed(self, _ppb: float) -> None:
        self._sync_viewport_fill()
        self._sync_names_column()
        hbar = self._bars_scroll.horizontalScrollBar()
        scroll_to = self._bars.take_pending_scroll_bp()
        anchor = self._bars.take_pending_anchor()
        if scroll_to is not None:
            target = int(self._bars.bp_to_x(float(scroll_to)) - 4)
            hbar.setValue(max(0, target))
        elif anchor is not None:
            anchor_bp, anchor_x = anchor
            new_x = self._bars.bp_to_x(float(anchor_bp))
            view_x = float(anchor_x) - hbar.value()
            hbar.setValue(max(0, int(new_x - view_x)))
            self._center_selection_vertically()
        elif self._selection_zoom_anchor() is not None:
            # Slider/ppb set without pending anchor — still keep feature in view.
            self._center_selection_vertically()
        self._update_scale()
        self._sync_zoom_slider_from_ppb()

    # ------------------------------------------------ Edit > Find / Go To

    def _on_rail_find_requested(self, direction: int) -> None:
        """Find tab: next / previous / find-all across multi-seq or MSA overview."""
        from PySide6.QtWidgets import QApplication

        from repeatviewer.backends.find_in_file import find_in_fasta
        from repeatviewer.backends.motif_search import parse_position, query_problem

        query = self._rail.find_query.text().strip()
        if not query:
            self._rail.set_find_status("Type something to find.")
            self.statusMessage.emit("Type something to find.")
            self._rail.find_query.setFocus()
            return
        position = parse_position(query)
        if position is not None and (
            position[0] is not None or not any(c.isalpha() for c in query)
        ):
            name, start, end = position
            problem = self.go_to_find_position(name, start, end)
            where = f"{start:,}" if start == end else f"{start:,}-{end:,}"
            msg = problem or (
                f"Went to position {where}"
                + (f" of {name}" if name else "")
                + "."
            )
            self._rail.set_find_status(msg)
            self.statusMessage.emit(msg)
            return

        request = self._rail.find_request()
        if request.in_sequence:
            problem = query_problem(query, request.options)
            if problem and not request.in_names and not request.in_annotations:
                self._rail.set_find_status(problem)
                self.statusMessage.emit(problem)
                return
        try:
            path = Path(self._path).expanduser().resolve()
            stat = path.stat()
        except OSError:
            self._rail.set_find_status("The file is no longer there.")
            self.statusMessage.emit("The file is no longer there.")
            return
        key = (str(path), stat.st_mtime_ns, request)
        if self._find_key != key or not self._find_items:
            features = self.features_by_id() if request.in_annotations else None
            try:
                QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
                result = find_in_fasta(path, request, features_by_id=features)
            except (OSError, ValueError) as exc:
                msg = f"Could not search: {exc}"
                self._rail.set_find_status(msg)
                self.statusMessage.emit(msg)
                return
            finally:
                QApplication.restoreOverrideCursor()
            self._find_items = list(result.items)
            self._find_key = key
            self._find_index = -1
            if result.sequence_skipped and not self._find_items:
                self._rail.set_find_status(result.sequence_skipped)
                self.statusMessage.emit(result.sequence_skipped)
                self.clear_find_hits()
                return
            if not self._find_items:
                self._rail.set_find_status("No matches.")
                self.statusMessage.emit("No matches.")
                self.clear_find_hits()
                self._find_key = key
                return

        if direction == 0:
            self._find_index = 0
        else:
            step = 1 if direction >= 0 else -1
            if self._find_index < 0:
                self._find_index = 0 if step > 0 else len(self._find_items) - 1
            else:
                self._find_index = (self._find_index + step) % len(self._find_items)

        item = self._find_items[self._find_index]
        if item.kind == "annotation":
            # Select every matching annotation (so Delete can remove them all).
            # Do not paint sequence find-hit underlays or a column band.
            ann_hits = [h for h in self._find_items if h.kind == "annotation"]
            ok = self.show_find_hit(
                item.seq_id,
                item.start,
                item.end,
                {},
                select_annotations=ann_hits,
            )
        else:
            hits_by_id: dict[str, list[tuple[int, int]]] = {}
            for hit in self._find_items:
                if hit.kind == "sequence":
                    hits_by_id.setdefault(hit.seq_id, []).append((hit.start, hit.end))
            ok = self.show_find_hit(
                item.seq_id,
                item.start,
                item.end,
                hits_by_id or {item.seq_id: [(item.start, item.end)]},
                whole_sequence=(item.kind == "name"),
            )
        if not ok:
            self._rail.set_find_status("Could not show that hit.")
            self.statusMessage.emit("Could not show that hit.")
            return
        kind = item.kind
        label = f" ({item.label})" if item.label else ""
        msg = (
            f"{self._find_index + 1} of {len(self._find_items)} — "
            f"{item.seq_id}:{item.start:,}-{item.end:,} [{kind}]{label}"
        )
        self._rail.set_find_status(msg)
        self.statusMessage.emit(msg)

    def go_to_find_position(
        self, seq_id: str | None, start: int, end: int
    ) -> str | None:
        """Jump to a typed position (optional ``name:start-end``)."""
        rid = seq_id or self.selected_record_id()
        if not rid:
            return (
                "This file has several sequences: type the name too, "
                "e.g. scaffold_5:1000-2000."
            )
        row = self._row_for_id(rid)
        if row < 0:
            return f'There is no sequence "{rid}" in {Path(self._path).name}.'
        rid = self._record_ids[row]
        length = self.record_length(rid)
        if length and start > length:
            return f"Position {start:,} is past the end of {rid} ({length:,} bp)."
        end = min(end, length) if length else end
        self.show_find_hit(rid, start, end, self._bars._find_hits)
        return None

    def record_ids(self) -> list[str]:
        return list(self._record_ids)

    def record_length(self, seq_id: str) -> int:
        return int(self._lengths.get(seq_id, 0))

    def selected_record_id(self) -> str | None:
        ids = self._selected_record_ids()
        if ids:
            return ids[0]
        # Regional box (e.g. Dotplot → Overview): topmost highlighted row.
        regional = self._bars.regional_rows()
        if regional:
            row = min(regional)
            if 0 <= row < len(self._record_ids):
                return self._record_ids[row]
        return None

    def record_ids_for_dotplot(self) -> list[str]:
        """Sequence ids in Overview display order (current Sorting)."""
        return list(self._bars.record_ids())

    def sequence_text_for_dotplot(self, seq_id: str) -> str:
        """Load one sequence for Dotplot (Self); MSA rows stay gapped."""
        return self._sequence_text_for_id(seq_id)

    def sequences_for_dotplot(self) -> list[tuple[str, str]]:
        """``(id, sequence)`` rows for Dotplot (Self); MSA sequences stay gapped.

        Prefer ``record_ids_for_dotplot`` + ``sequence_text_for_dotplot`` for
        large files — this eager helper loads every sequence.
        """
        out: list[tuple[str, str]] = []
        for rid in self.record_ids_for_dotplot():
            try:
                out.append((rid, self._sequence_text_for_id(rid)))
            except (OSError, ValueError):
                continue
        return out

    def select_by_position(
        self,
        start: int,
        end: int,
        *,
        sequence_id: str | None = None,
    ) -> None:
        """Highlight a 1-based inclusive range.

        With ``sequence_id``, highlight only that sequence's span (regional box),
        not a full-height column across every track. Without it, keep a column
        band (scale / feature Selection).
        """
        maximum = max(1, int(self._bars.max_length()))
        lo, hi = (int(start), int(end)) if start <= end else (int(end), int(start))
        lo = max(1, min(maximum, lo))
        hi = max(lo, min(maximum, hi))
        row = self._row_for_id(sequence_id) if sequence_id else -1
        if row >= 0:
            self._bars.set_regional_selection(lo, hi, [row])
            self._scroll_to_row(row)
        else:
            self._bars.set_region_selection(lo, hi)
        self._sync_header_position_spins()
        # Inclusive columns → half-open window (same as Select / Find).
        self._on_region_zoom(float(lo - 1), float(hi))

    def features_by_id(self) -> dict[str, list[SequenceFeature]]:
        return {rid: list(feats) for rid, feats in self._features.items()}

    def _row_for_id(self, seq_id: str) -> int:
        wanted = str(seq_id or "").strip()
        if wanted in self._record_ids:
            return self._record_ids.index(wanted)
        from repeatviewer.backends.gff_features import _seqid_matches

        for row, rid in enumerate(self._record_ids):
            if _seqid_matches(rid, wanted):
                return row
        return -1

    def show_find_hit(
        self,
        seq_id: str,
        start: int,
        end: int,
        hits_by_id: dict[str, list[tuple[int, int]]],
        *,
        whole_sequence: bool = False,
        select_annotations: list | None = None,
    ) -> bool:
        """Scroll/zoom to a Find hit. Annotation hits select features only."""
        row = self._row_for_id(seq_id)
        if row < 0:
            return False
        rid = self._record_ids[row]

        if select_annotations is not None:
            feats: set[tuple[int, int]] = set()
            for hit in select_annotations:
                hit_row = self._row_for_id(getattr(hit, "seq_id", "") or "")
                if hit_row < 0:
                    continue
                for flat_i in self._bars.flat_indices_for_annotation(
                    hit_row,
                    int(getattr(hit, "start", start)),
                    int(getattr(hit, "end", end)),
                    str(getattr(hit, "label", "") or ""),
                ):
                    feats.add((hit_row, flat_i))
            # Annotation search: select features only — no sequence underlays,
            # no column band, no sequence-row selection.
            self._bars.set_find_hits({}, None)
            self._bars.set_selected_rows([])
            if self._bars.region_selection() is not None:
                self._bars.set_region_selection(None)
            self._bars._set_feature_selection(feats, set_columns=False)
            self._sync_names_column()
            self._scroll_to_row(row)
            span = max(1, int(end) - int(start) + 1)
            context = max(3 * span, 60)
            center = (int(start) + int(end)) / 2.0
            self._on_region_zoom(max(0.0, center - context / 2.0), center + context / 2.0)
            return True

        self._bars.set_find_hits(hits_by_id, (rid, int(start), int(end)))
        self._bars.clear_feature_selection(repaint=False)
        self._bars.set_selected_rows([row])
        self._sync_names_column()
        self._scroll_to_row(row)
        if whole_sequence:
            self._on_region_zoom(0.0, float(max(1, self._lengths.get(rid, end))))
        else:
            span = max(1, int(end) - int(start) + 1)
            context = max(3 * span, 60)
            center = (int(start) + int(end)) / 2.0
            self._on_region_zoom(max(0.0, center - context / 2.0), center + context / 2.0)
        return True

    def clear_find_hits(self) -> None:
        self._bars.set_find_hits({}, None)
        self._find_items = []
        self._find_index = -1
        self._find_key = None

    def add_features_and_save(self, new_by_id: dict[str, list[SequenceFeature]]) -> bool:
        """Append features per record and save the annotations (JSON sidecar)."""
        added = {
            rid: feats
            for rid, feats in new_by_id.items()
            if feats and rid in self._lengths
        }
        if not added:
            return False
        self._push_feature_undo_snapshot()
        for rid, feats in added.items():
            self._features[rid] = list(self._features.get(rid) or []) + list(feats)
        try:
            json_path, _gff = save_annotations_multi(
                self._path,
                features_by_id=self._features,
                lengths_by_id=self._lengths,
                record_ids=list(self._record_ids),
            )
        except OSError as exc:
            QMessageBox.critical(self, "Save failed", str(exc))
            return False
        keep_ppb = self._bars.pixels_per_base()
        keep_sel = self._bars.selected_rows()
        keep_region = self._bars.region_selection()
        was_msa = self._bars.is_msa_mode()
        if was_msa:
            self._bars.set_features_in_place(self._features)
        else:
            self._bars.set_data(
                self._record_ids, self._lengths, self._features, fit_bars_width=None
            )
            self._bars.set_pixels_per_base(keep_ppb)
        self._bars.set_selected_rows(keep_sel)
        if keep_region is not None:
            self._bars.set_region_selection(keep_region[0], keep_region[1])
        self._sync_annotation_kind_filters()
        self._sync_names_column()
        self._update_scale()
        self.annotationsSaved.emit(json_path)
        return True

    def _current_region_for_add(self) -> tuple[int, int] | None:
        """Region from canvas selection or header Position spins."""
        region = self._bars.region_selection()
        if region is not None:
            return region
        bar = self._header_zoom_bar
        if bar is None:
            return None
        start = int(bar.pos_start.value())
        end = int(bar.pos_end.value())
        if end < start:
            start, end = end, start
        return start, end

    def _on_add_feature_clicked(self) -> None:
        from repeatviewer.widgets.sequence_map_viewer import AddFeatureDialog

        target_ids = self._target_record_ids_for_manual_add()
        if not target_ids:
            QMessageBox.information(
                self,
                "Add feature",
                "Select one or more sequences (click names, or drag a box "
                "on the tracks), then Add.\n\n"
                "Tip: Ctrl+A (Cmd+A on Mac) selects all sequences.",
            )
            return
        region = self._current_region_for_add()
        if region is None:
            QMessageBox.information(
                self,
                "Add feature",
                "Drag on the tracks to select a region (or set Position), "
                "then Add.",
            )
            return
        start, end = region
        label, feature_type = self._rail.selected_manual_feature()
        existing = 0
        for rid in target_ids:
            existing += sum(
                1
                for feat in (self._features.get(rid) or [])
                if feat.feature_type == feature_type
            )
        if feature_type == "ORF":
            default_name = f"ORF_{existing + 1}"
        else:
            default_name = f"{(feature_type or label or 'feature').replace(' ', '_')}{existing + 1}"
        seq_len = max(
            (int(self._lengths.get(rid, 1)) for rid in target_ids),
            default=max(1, int(self._bars.max_length())),
        )
        dialog = AddFeatureDialog(
            start=start,
            end=end,
            seq_length=seq_len,
            parent=self,
            feature_type=feature_type,
            feature_label=label,
            default_name=default_name,
        )
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        feat = dialog.feature()
        if feat is None:
            QMessageBox.warning(self, "Add feature", "Enter a feature name.")
            return
        new_by_id: dict[str, list[SequenceFeature]] = {}
        for rid in target_ids:
            rid_len = max(1, int(self._lengths.get(rid, 1)))
            a = max(1, min(rid_len, int(feat.start)))
            b = max(1, min(rid_len, int(feat.end)))
            if b < a:
                a, b = b, a
            new_by_id[rid] = [
                SequenceFeature(
                    name=feat.name,
                    feature_type=feat.feature_type,
                    start=a,
                    end=b,
                    strand=feat.strand,
                    color=feat.color,
                    track=feat.track,
                    source=feat.source,
                    parent=feat.parent,
                )
            ]
        if self.add_features_and_save(new_by_id):
            n = len(target_ids)
            self.statusMessage.emit(
                f"Added “{feat.name}” to {n:,} sequence{'s' if n != 1 else ''}."
            )

    def _sequence_text_for_id(self, seq_id: str) -> str:
        """Return sequence string for a record (MSA cache or disk)."""
        wanted = str(seq_id or "").strip()
        if self._bars.is_msa_mode():
            ids = self._bars.record_ids()
            try:
                row = ids.index(wanted)
            except ValueError:
                row = -1
            if 0 <= row < len(self._bars._msa_seqs):
                return str(self._bars._msa_seqs[row])
        # Multi-FASTA: cache the file once so Dotplot switches are not O(n) parses.
        cache = getattr(self, "_fasta_seq_cache", None)
        if cache is None:
            from Bio import SeqIO

            cache = {}
            for rec in SeqIO.parse(str(self._path), "fasta"):
                cache[str(rec.id)] = str(rec.seq)
            self._fasta_seq_cache = cache
        if wanted in cache:
            return cache[wanted]
        short = wanted.split("#", 1)[0]
        for rid, seq in cache.items():
            if rid == wanted or rid.split("#", 1)[0] == short:
                return seq
        raise ValueError(f"Sequence “{seq_id}” was not found in {self._path.name}.")

    def _on_transfer_annotations_clicked(self) -> None:
        if self._transfer_thread is not None and self._transfer_thread.isRunning():
            self._transfer_cancel.request()
            self._rail.transfer_run_btn.setText("Stopping…")
            return
        try:
            self._start_overview_transfer_worker()
        except Exception as exc:  # noqa: BLE001
            QMessageBox.warning(self, "Transfer", str(exc))

    def _start_overview_transfer_worker(self) -> None:
        from repeatviewer.backends.annotation_transfer import (
            AnnotationTransferOptions,
            discover_annotated_transfer_sources,
            load_use_system_mafft_from_settings,
            transfer_annotations,
        )
        from repeatviewer.backends.similarity_annotation import (
            load_use_system_blast_from_settings,
        )

        target_ids = self._require_selected_record_ids("Transfer")
        if not target_ids:
            return
        source_path = None
        if hasattr(self._rail, "transfer_source_row"):
            source_path = self._rail.transfer_source_row.path()
        if source_path is None:
            QMessageBox.information(
                self,
                "Transfer",
                "Choose a source FASTA file or a folder of annotated FASTAs.",
            )
            return
        sources = discover_annotated_transfer_sources(source_path)
        if not sources:
            kind = "folder" if source_path.is_dir() else "file"
            QMessageBox.warning(
                self,
                "Transfer",
                f"No annotated sequences found in this {kind}.\n\n"
                "Only FASTA files with a RepeatViewer annotation JSON "
                "(*.seqmap_RV.json) and at least one feature are used.\n\n"
                f"Source:\n{source_path}",
            )
            return

        source_label = source_path.name
        if source_path.is_dir():
            source_label = f"{source_path.name}/ ({len(sources)} annotated)"
        ann_source = f"Transfer ({source_label})"
        opts = AnnotationTransferOptions(
            min_identity=float(self._rail.transfer_min_identity_value()),
            min_coverage=float(self._rail.transfer_min_coverage_value()),
            allow_reverse=bool(self._rail.transfer_allow_reverse.isChecked()),
            annotation_source=ann_source,
            use_system_blast=load_use_system_blast_from_settings(),
            use_system_mafft=load_use_system_mafft_from_settings(),
        )

        self._transfer_cancel = JobCancellation()
        self._rail.transfer_run_btn.setText("Stop")
        self._rail.transfer_run_btn.setFixedWidth(52)
        n_targets = len(target_ids)
        self.statusMessage.emit(
            f"Transferring annotations onto {n_targets:,} selected sequence(s)…"
        )
        self._begin_window_job_progress(
            "Transfer annotation", max(1, n_targets * max(1, len(sources)))
        )

        job_sources = list(sources)
        job_targets = list(target_ids)
        job_path = self._path
        job_opts = opts
        cancel = self._transfer_cancel
        # Preload target sequences on the UI thread (worker must not touch widgets).
        target_seqs: dict[str, str] = {}
        for rid in job_targets:
            try:
                target_seqs[rid] = self._sequence_text_for_id(rid)
            except ValueError as exc:
                QMessageBox.warning(self, "Transfer", str(exc))
                return

        def _job(*, progress=None):
            from Bio import SeqIO

            source_seqs: dict[tuple[str, str], tuple[str, list]] = {}
            for entry in job_sources:
                key = (str(entry.fasta_path), str(entry.sequence_id))
                if key in source_seqs:
                    continue
                feats = []
                seq = ""
                try:
                    from repeatviewer.backends.sequence_features import (
                        load_features_by_id_from_json as _load_json,
                    )

                    by_id = _load_json(entry.fasta_path, [entry.sequence_id])
                    feats = list(by_id.get(entry.sequence_id) or [])
                    for rec in SeqIO.parse(str(entry.fasta_path), "fasta"):
                        if str(rec.id) == entry.sequence_id:
                            seq = str(rec.seq)
                            break
                except Exception:
                    continue
                source_seqs[key] = (seq, feats)

            new_by_id: dict[str, list] = {rid: [] for rid in job_targets}
            step = 0
            total = max(1, len(job_targets) * max(1, len(job_sources)))
            for tgt_id in job_targets:
                if cancel.requested:
                    from repeatviewer.backends.job_cancel import JobCancelledError

                    raise JobCancelledError("Annotation transfer cancelled.")
                tgt_seq = target_seqs.get(tgt_id) or ""
                if not tgt_seq:
                    continue
                for entry in job_sources:
                    step += 1
                    if progress is not None:
                        progress(f"{step}/{total} {tgt_id} ← {entry.sequence_id}")
                    key = (str(entry.fasta_path), str(entry.sequence_id))
                    src_seq, src_feats = source_seqs.get(key, ("", []))
                    if not src_seq or not src_feats:
                        continue
                    try:
                        same_file = (
                            entry.fasta_path.resolve() == job_path.resolve()
                            and entry.sequence_id == tgt_id
                        )
                    except OSError:
                        same_file = False
                    if same_file:
                        continue
                    # MSA rows: strip gaps for transfer, project hits back to columns.
                    from repeatviewer.backends.msa_coordinates import (
                        project_features_to_gapped,
                        project_features_to_ungapped,
                        ungap_sequence,
                    )

                    src_u, src_keep = ungap_sequence(src_seq)
                    tgt_u, tgt_keep = ungap_sequence(tgt_seq)
                    src_feats_u = project_features_to_ungapped(
                        list(src_feats), src_keep
                    )
                    if not src_u or not tgt_u or not src_feats_u:
                        continue
                    result = transfer_annotations(
                        source_id=entry.sequence_id,
                        source_sequence=src_u,
                        target_id=tgt_id,
                        target_sequence=tgt_u,
                        features=list(src_feats_u),
                        opts=job_opts,
                        cancel=cancel,
                        source_fasta=entry.fasta_path,
                    )
                    for item in getattr(result, "lifted", []) or []:
                        feat = getattr(item, "feature", None)
                        if feat is None:
                            continue
                        projected = project_features_to_gapped([feat], tgt_keep)
                        new_by_id[tgt_id].extend(projected)
            return new_by_id

        worker = Worker(_job, progress=True)
        thread = QThread(self)
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        worker.signals.progress.connect(self._on_transfer_progress)
        worker.signals.finished.connect(self._on_transfer_finished)
        worker.signals.error.connect(self._on_transfer_error)
        worker.signals.finished.connect(thread.quit)
        worker.signals.error.connect(thread.quit)
        thread.finished.connect(worker.deleteLater)
        thread.finished.connect(thread.deleteLater)
        thread.finished.connect(self._clear_transfer_thread)
        self._transfer_worker = worker
        self._transfer_thread = thread
        thread.start()

    @Slot(str)
    def _on_transfer_progress(self, message: str) -> None:
        text = str(message or "").strip()
        match = re.match(r"^(\d+)\s*/\s*(\d+)\s+(.+)$", text)
        if match is None:
            self.statusMessage.emit(text)
            return
        completed, _total, detail = (
            int(match.group(1)),
            int(match.group(2)),
            match.group(3).strip(),
        )
        self._update_window_job_progress(completed, detail)
        self.statusMessage.emit(f"Transfer {completed}: {detail}")

    @Slot(object)
    def _on_transfer_finished(self, result: object) -> None:
        self._rail.transfer_run_btn.setText("Run")
        self._rail.transfer_run_btn.setFixedWidth(52)
        self._end_window_job_progress()
        if not isinstance(result, dict):
            return
        new_by_id = {
            rid: list(feats)
            for rid, feats in result.items()
            if feats and rid in self._lengths
        }
        if not new_by_id:
            self.statusMessage.emit("Transfer finished — no features placed.")
            return
        if self.add_features_and_save(new_by_id):
            n = sum(len(v) for v in new_by_id.values())
            self.statusMessage.emit(
                f"Transfer saved {n:,} feature(s) on {len(new_by_id):,} sequence(s)."
            )

    @Slot(str)
    def _on_transfer_error(self, message: str) -> None:
        self._rail.transfer_run_btn.setText("Run")
        self._rail.transfer_run_btn.setFixedWidth(52)
        self._end_window_job_progress()
        lower = str(message or "").lower()
        if "cancelled" in lower or "stopped" in lower:
            self.statusMessage.emit("Transfer stopped.")
            return
        QMessageBox.warning(self, "Transfer failed", str(message))
        self.statusMessage.emit(f"Transfer failed: {message}")

    @Slot()
    def _clear_transfer_thread(self) -> None:
        self._transfer_thread = None
        self._transfer_worker = None
        self._rail.transfer_run_btn.setText("Run")
        self._rail.transfer_run_btn.setFixedWidth(52)

    def _on_region_zoom(self, start_bp: float, end_bp: float) -> None:
        self._sync_viewport_fill()
        self._bars.zoom_to_bp_range(start_bp, end_bp, self._bars_viewport_width())

    def _on_name_row_clicked(self, row: int, modifiers: object) -> None:
        self._bars.finish_live_selection_preview()
        self._bars._select_row_click(row, modifiers)  # type: ignore[arg-type]

    def _on_name_rows_preview(self, rows: object) -> None:
        """Snapped row preview (legacy path); prefer pixel preview while dragging."""
        row_list = [int(r) for r in (rows or [])]
        self._bars.set_row_selection_preview(row_list)

    def _on_name_rows_pixel_preview(self, y0: float, y1: float) -> None:
        """Live row outline while dragging (names highlight is throttled separately)."""
        # Clear scale band once per drag — not on every move.
        if not getattr(self, "_name_row_drag_cleared_scale", False):
            self._scale.set_region_selection(None)
            self._name_row_drag_cleared_scale = True
        mid_y = int((float(y0) + float(y1)) * 0.5)
        # Keep the current horizontal zoom window. Passing x=0 made
        # ``_ensure_drag_point_visible`` treat the drag as left of the viewport
        # and jump the MSA back to bp 1 whenever several name rows were selected.
        hbar = self._bars_scroll.horizontalScrollBar()
        mid_x = int(hbar.value()) + max(1, self._bars_scroll.viewport().width()) // 2
        self._bars._ensure_drag_point_visible(QPoint(mid_x, mid_y))  # type: ignore[attr-defined]
        self._bars.set_row_pixel_preview(float(y0), float(y1))

    def _on_name_rows_dragged(self, rows: object) -> None:
        self._name_row_drag_cleared_scale = False
        row_list = [int(r) for r in (rows or [])]
        self._bars.finish_live_selection_preview()
        self._bars.set_selected_rows(row_list)

    def _on_scale_column_region(self, start: int, end: int) -> None:
        """Column selection from dragging the bottom scale (no row change)."""
        timer = getattr(self, "_scale_preview_timer", None)
        if timer is not None:
            timer.stop()
        self._scale_preview_pending = None
        self._bars.finish_live_selection_preview()
        # Clear feature picks without a full-canvas repaint — column chrome
        # update below covers the visible change.
        self._bars.clear_feature_selection(repaint=False)
        self._bars.set_region_selection(int(start), int(end))
        # Keep name-column highlights from a prior sequence-row select.
        self._sync_names_column()
        self._bars.setFocus(Qt.FocusReason.MouseFocusReason)
        # Position spins update via regionSelectionChanged (avoid double work).

    def _flush_scale_preview_positions(self) -> None:
        pending = self._scale_preview_pending
        if pending is None:
            return
        start, end = pending
        bar = self._header_zoom_bar
        if bar is not None and hasattr(bar, "set_position_range"):
            maximum = max(1, int(self._bars.max_length()))
            bar.set_position_range(maximum, start=start, end=end)

    def _on_scale_column_pixel_preview(self, x0: float, x1: float) -> None:
        """Pixel-follow column guides + live Position spins while dragging the scale."""
        scroll = int(self._bars_scroll.horizontalScrollBar().value())
        # Scale x is viewport-local; bars chrome uses content coordinates.
        cx0 = float(x0) + scroll
        cx1 = float(x1) + scroll
        self._bars.set_column_pixel_preview(cx0, cx1)
        # Autoscroll when the drag edge leaves the visible MSA viewport.
        self._bars._ensure_drag_point_visible(  # type: ignore[attr-defined]
            QPoint(int(cx1), 0)
        )
        start_i, end_i = self._scale._columns_from_x_range(float(x0), float(x1))
        if self._scale_preview_pending == (start_i, end_i):
            return
        self._scale_preview_pending = (start_i, end_i)
        # Immediate first update, then coalesce further moves (~20 Hz).
        if not self._scale_preview_timer.isActive():
            self._flush_scale_preview_positions()
            self._scale_preview_timer.start()

    def _on_scale_column_preview(self, region: object) -> None:
        """Snapped column preview (legacy); prefer pixel preview while dragging."""
        if region is None:
            return
        try:
            start, end = int(region[0]), int(region[1])  # type: ignore[index]
        except (TypeError, ValueError, IndexError):
            return
        self._bars.set_column_selection_preview(start, end)

    def _on_feature_activated(self, row: int, feat_index: int) -> None:
        self._modify_feature_at(row, feat_index)

    def _selected_record_ids(self) -> list[str]:
        ids: list[str] = []
        for row in self._bars.selected_rows():
            if 0 <= row < len(self._record_ids):
                ids.append(self._record_ids[row])
        return ids

    def _on_names_context_menu(self, global_pos: QPoint, row: int) -> None:
        if row >= 0 and row not in self._bars.selected_row_set():
            self._bars.set_selected_rows([row])
            self._sync_names_column()
        self._show_overview_context_menu(global_pos)

    def _on_bars_context_menu(self, global_pos: QPoint) -> None:
        self._show_overview_context_menu(global_pos)

    def _build_overview_context_menu(self) -> QMenu:
        """Right-click menu aligned with single-sequence map (plus Rename)."""
        menu = QMenu(self)
        menu.setObjectName("OverviewSequenceMenu")
        region = self._bars.region_selection()
        selected_rows = self._bars.selected_rows()
        # Name-column / whole-sequence selection (not a regional matrix box).
        sequence_sel = bool(selected_rows) and not self._bars.is_regional_selection()

        add_action = menu.addAction("Add feature…")
        add_action.triggered.connect(self._on_add_feature_clicked)

        if sequence_sel:
            menu.addSeparator()
            rename = menu.addAction("Rename sequence")
            rename.setEnabled(len(selected_rows) == 1)
            rename.triggered.connect(self._rename_selected_sequence)
            if self.is_msa_mode():
                set_std = menu.addAction("Set as standard")
                set_std.setEnabled(len(selected_rows) == 1)
                set_std.triggered.connect(self._set_selected_as_msa_standard)

        menu.addSeparator()
        copy_seq = menu.addAction("Copy Sequence")
        copy_ungapped = None
        if self.is_msa_mode():
            copy_ungapped = menu.addAction("Copy Sequence without gap")
        extract = menu.addAction("Extract Sequence")
        extract_ungapped = None
        if self.is_msa_mode():
            extract_ungapped = menu.addAction("Extract Sequence without gap")
        if region is not None:
            start, end = region
            copy_seq.triggered.connect(
                lambda _c=False, a=start, b=end: self._copy_region_sequence(a, b)
            )
            extract.triggered.connect(
                lambda _c=False, a=start, b=end: self._extract_region_sequence(a, b)
            )
            if copy_ungapped is not None:
                copy_ungapped.triggered.connect(
                    lambda _c=False, a=start, b=end: self._copy_region_sequence(
                        a, b, without_gaps=True
                    )
                )
            if extract_ungapped is not None:
                extract_ungapped.triggered.connect(
                    lambda _c=False, a=start, b=end: self._extract_region_sequence(
                        a, b, without_gaps=True
                    )
                )
        elif sequence_sel:
            copy_seq.triggered.connect(self._copy_selected_sequences)
            extract.triggered.connect(self._extract_selected)
            if copy_ungapped is not None:
                copy_ungapped.triggered.connect(
                    lambda _c=False: self._copy_selected_sequences(without_gaps=True)
                )
            if extract_ungapped is not None:
                extract_ungapped.triggered.connect(
                    lambda _c=False: self._extract_selected(without_gaps=True)
                )
        else:
            copy_seq.setEnabled(False)
            extract.setEnabled(False)
            if copy_ungapped is not None:
                copy_ungapped.setEnabled(False)
            if extract_ungapped is not None:
                extract_ungapped.setEnabled(False)

        show_dot = menu.addAction("Show self-dotplot")
        if selected_rows:
            show_dot.triggered.connect(self._show_self_dotplot_for_selection)
            if len(selected_rows) > 1:
                show_dot.setToolTip(
                    "Opens Dotplot (Self) for the topmost selected sequence."
                )
        else:
            show_dot.setEnabled(False)

        menu.addSeparator()
        delete = menu.addAction("Delete")
        can_delete = bool(
            self._bars.selected_features()
            or (self.is_msa_mode() and region is not None)
            or sequence_sel
        )
        delete.setEnabled(can_delete)
        if self.is_msa_mode():
            delete.setToolTip(
                "Remove the selection (sequences, region, or feature). "
                "Shift+Delete fills MSA residues with gaps."
            )
        else:
            delete.setToolTip("Remove selected sequence(s) or feature(s).")
        delete.triggered.connect(self._on_delete_clicked)
        return menu

    def _show_overview_context_menu(self, global_pos: QPoint) -> None:
        self._build_overview_context_menu().exec(global_pos)

    def _show_self_dotplot_for_selection(self) -> None:
        """Open Dotplot (Self) for the topmost selected sequence."""
        rows = self._bars.selected_rows()
        if not rows:
            return
        row = int(rows[0])  # selected_rows() is sorted → visually topmost
        if not (0 <= row < len(self._record_ids)):
            return
        self.showSelfDotplotRequested.emit(self._record_ids[row])

    def _rename_selected_sequence(self) -> None:
        """Rename one selected sequence in the FASTA and annotation sidecar."""
        from repeatviewer.backends.sequence_features import fasta_id_token
        from repeatviewer.widgets.file_conflict import prompt_rename_text

        rows = self._bars.selected_rows()
        if len(rows) != 1:
            return
        row = int(rows[0])
        if not (0 <= row < len(self._record_ids)):
            return
        old_id = self._record_ids[row]
        new_text, ok = prompt_rename_text(
            self,
            title="Rename sequence",
            label="New sequence name:",
            initial=old_id,
        )
        if not ok:
            return
        new_id = fasta_id_token(new_text)
        if not new_id:
            QMessageBox.warning(self, "Rename sequence", "Enter a sequence name.")
            return
        if new_id == old_id:
            return
        if new_id in self._record_ids:
            QMessageBox.warning(
                self,
                "Rename sequence",
                f'“{new_id}” is already used by another sequence.',
            )
            return
        try:
            sequences = list(self._bars.msa_sequences())
            if len(sequences) != len(self._record_ids):
                _ids, sequences = load_fasta_rows(self._path)
            if len(sequences) != len(self._record_ids):
                raise ValueError(
                    "FASTA record count no longer matches the overview."
                )
            new_ids = list(self._record_ids)
            new_ids[row] = new_id
            if self._bars.is_msa_mode():
                write_aligned_fasta_rows(self._path, new_ids, sequences)
            else:
                write_fasta_rows(self._path, new_ids, sequences)
        except (OSError, ValueError) as exc:
            QMessageBox.warning(self, "Rename sequence", str(exc))
            return

        feats = list(self._features.pop(old_id, []) or [])
        self._features[new_id] = feats
        length = int(self._lengths.pop(old_id, self._lengths.get(new_id, 0)))
        self._lengths[new_id] = length
        self._record_ids = new_ids
        was_standard = self._bars.msa_diff_standard_id() == old_id
        keep_ppb = self._bars.pixels_per_base()
        keep_region = self._bars.region_selection()
        if self._bars.is_msa_mode():
            self._bars.set_msa_data(
                new_ids,
                sequences,
                features=self._features,
            )
            self._bars.set_pixels_per_base(keep_ppb)
            if was_standard:
                self._bars.set_msa_diff_standard_id(new_id)
                self._persist_msa_diff_standard_id(new_id)
            else:
                self._restore_msa_diff_standard_from_settings()
        else:
            self._bars.set_data(
                new_ids,
                self._lengths,
                self._features,
                sequences=sequences,
                fit_bars_width=None,
            )
            self._bars.set_pixels_per_base(keep_ppb)
        self._bars.set_selected_rows([row])
        if keep_region is not None:
            self._bars.set_region_selection(keep_region[0], keep_region[1])
        self._sync_names_column()
        self._update_scale()
        try:
            json_path, _gff = save_annotations_multi(
                self._path,
                features_by_id=self._features,
                lengths_by_id=self._lengths,
                record_ids=list(self._record_ids),
            )
            self.annotationsSaved.emit(json_path)
        except OSError as exc:
            QMessageBox.warning(
                self,
                "Rename sequence",
                f"Sequence renamed, but annotations could not be saved:\n{exc}",
            )
            return
        self.statusMessage.emit(f"Renamed “{old_id}” → “{new_id}”")

    @staticmethod
    def _strip_alignment_gaps(sequence: str) -> str:
        """Remove MSA gap characters (``-`` / ``.``) from sequence text."""
        from repeatviewer.backends.msa_coordinates import ungap_sequence

        return ungap_sequence(sequence)[0]

    def _copy_region_sequence(
        self, start: int, end: int, *, without_gaps: bool = False
    ) -> None:
        """Copy the selected region from target sequence(s) to the clipboard."""
        title = "Copy Sequence without gap" if without_gaps else "Copy Sequence"
        targets = self._target_record_ids_for_manual_add()
        if not targets:
            QMessageBox.information(
                self,
                title,
                "Select one or more sequences (click names, or drag a box "
                f"on the tracks), then {title}.",
            )
            return
        lo, hi = (start, end) if start <= end else (end, start)
        chunks: list[str] = []
        total = 0
        try:
            for rid in targets:
                seq = self._sequence_text_for_id(rid)
                frag = seq[lo - 1 : hi]
                if without_gaps:
                    frag = self._strip_alignment_gaps(frag)
                if not frag:
                    continue
                total += len(frag)
                if len(targets) == 1:
                    chunks.append(frag)
                else:
                    chunks.append(f">{rid}\n{frag}")
        except (OSError, ValueError) as exc:
            QMessageBox.warning(self, title, str(exc))
            return
        if not chunks:
            QMessageBox.warning(self, title, "The selection is empty.")
            return
        QApplication.clipboard().setText("\n".join(chunks))
        unit = "aa" if self._bars.is_msa_mode() and self._bars._msa_protein else "bp"
        gap_note = " without gaps" if without_gaps else ""
        self.statusMessage.emit(
            f"Copied {total:,} {unit}{gap_note} ({lo:,}-{hi:,}) from "
            f"{len(targets):,} sequence{'s' if len(targets) != 1 else ''}"
        )

    def _copy_selected_sequences(self, *, without_gaps: bool = False) -> None:
        """Copy whole selected sequence(s) to the clipboard as FASTA text."""
        title = "Copy Sequence without gap" if without_gaps else "Copy Sequence"
        ids = self._selected_record_ids()
        if not ids:
            return
        chunks: list[str] = []
        try:
            for rid in ids:
                seq = self._sequence_text_for_id(rid)
                if without_gaps:
                    seq = self._strip_alignment_gaps(seq)
                chunks.append(f">{rid}\n{seq}")
        except (OSError, ValueError) as exc:
            QMessageBox.warning(self, title, str(exc))
            return
        QApplication.clipboard().setText("\n".join(chunks))
        gap_note = " without gaps" if without_gaps else ""
        self.statusMessage.emit(
            f"Copied {len(ids):,} sequence{'s' if len(ids) != 1 else ''}"
            f"{gap_note} to the clipboard"
        )

    def _extract_region_sequence(
        self, start: int, end: int, *, without_gaps: bool = False
    ) -> None:
        """Extract the selected region from target sequence(s) to a new FASTA."""
        from repeatviewer.widgets.sequence_map_viewer import _ExtractSequenceDialog

        title = "Extract Sequence without gap" if without_gaps else "Extract Sequence"
        targets = self._target_record_ids_for_manual_add()
        if not targets:
            QMessageBox.information(
                self,
                title,
                "Select one or more sequences (click names, or drag a box "
                f"on the tracks), then {title}.",
            )
            return
        lo, hi = (start, end) if start <= end else (end, start)
        primary = targets[0]
        try:
            primary_seq = self._sequence_text_for_id(primary)
        except (OSError, ValueError) as exc:
            QMessageBox.warning(self, title, str(exc))
            return
        seq_len = max(1, len(primary_seq))
        dialog = _ExtractSequenceDialog(
            sequence_id=primary if len(targets) == 1 else f"{len(targets)}_sequences",
            sequence_length=seq_len,
            start=lo,
            end=min(hi, seq_len),
            default_name=f"{primary}_{lo}_{hi}",
            reverse=False,
            parent=self,
        )
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        chosen_start, chosen_end, name, reverse = dialog.values()
        from repeatviewer.backends.dotplot import reverse_complement

        safe = "".join(ch if ch.isalnum() or ch in "-_." else "_" for ch in name)
        out_path = self._path.parent / f"{safe}.fasta"
        if out_path.exists():
            answer = QMessageBox.question(
                self,
                title,
                f"{out_path.name} already exists. Overwrite it?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No,
            )
            if answer != QMessageBox.StandardButton.Yes:
                return
        lines: list[str] = []
        try:
            for rid in targets:
                seq = self._sequence_text_for_id(rid)
                a = max(1, min(int(chosen_start), len(seq)))
                b = max(1, min(int(chosen_end), len(seq)))
                if b < a:
                    a, b = b, a
                frag = seq[a - 1 : b]
                if without_gaps:
                    frag = self._strip_alignment_gaps(frag)
                if reverse:
                    frag = reverse_complement(frag)
                strand = "-" if reverse else "+"
                header = (
                    f">{name}"
                    if len(targets) == 1
                    else f"{name}_{rid}"
                )
                lines.append(f">{header} {rid}:{a}-{b}({strand})")
                for offset in range(0, len(frag), 60):
                    lines.append(frag[offset : offset + 60])
            out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        except (OSError, ValueError) as exc:
            QMessageBox.warning(self, title, str(exc))
            return
        self.sequenceExtracted.emit(out_path)
        gap_note = " without gaps" if without_gaps else ""
        self.statusMessage.emit(
            f"Extracted {chosen_start:,}-{chosen_end:,}{gap_note} from "
            f"{len(targets):,} sequence{'s' if len(targets) != 1 else ''} → "
            f"{out_path.name}"
        )

    def _extract_selected(self, *, without_gaps: bool = False) -> None:
        from repeatviewer.backends.sequence_compose import (
            SequenceComposeError,
            extract_records,
            proposed_extract_path,
        )
        from repeatviewer.widgets.file_conflict import confirm_existing_file_path

        title = "Extract Sequence without gap" if without_gaps else "Extract"
        ids = self._selected_record_ids()
        if not ids:
            QMessageBox.information(
                self,
                title,
                "Select one or more sequences, then Extract.",
            )
            return

        proposed = proposed_extract_path(
            self._path, ids, self._path.parent, single_file=True
        )
        confirmed = confirm_existing_file_path(
            self,
            proposed,
            title=title,
            informative=(
                "Overwrite the existing file, or save under a different name?"
            ),
        )
        if confirmed is None:
            return

        if without_gaps:
            try:
                lines: list[str] = []
                for rid in ids:
                    seq = self._strip_alignment_gaps(self._sequence_text_for_id(rid))
                    lines.append(f">{rid}")
                    for offset in range(0, len(seq), 60):
                        lines.append(seq[offset : offset + 60])
                confirmed.write_text("\n".join(lines) + "\n", encoding="utf-8")
                written = [confirmed]
            except (OSError, ValueError) as exc:
                QMessageBox.warning(self, title, str(exc))
                return
        else:
            try:
                written = extract_records(
                    self._path,
                    set(ids),
                    self._path.parent,
                    single_file=True,
                    output_path=confirmed,
                )
            except SequenceComposeError as exc:
                QMessageBox.warning(self, title, str(exc))
                return

        for path in written:
            self.sequenceExtracted.emit(path)
        names = ", ".join(p.name for p in written[:4])
        if len(written) > 4:
            names += ", …"
        gap_note = " without gaps" if without_gaps else ""
        self.statusMessage.emit(
            f"Extracted {len(written)} file(s){gap_note}: {names}"
        )

    def _on_feature_resized(
        self, rid: str, store_index: int, feat: object
    ) -> None:
        """Commit an edge-drag resize from the overview canvas."""
        if not isinstance(feat, SequenceFeature):
            return
        store = list(self._features.get(rid) or [])
        if not (0 <= int(store_index) < len(store)):
            return
        self._push_feature_undo_snapshot()
        store[int(store_index)] = feat
        self._features[rid] = store
        self._refresh_features_after_edit()
        gff_path = self._persist_overview_features()
        label = feat.name or feat.feature_type
        if gff_path is not None:
            self.statusMessage.emit(
                f"Resized “{label}” on {rid} to {feat.start}–{feat.end} "
                f"→ {gff_path.name}"
            )
        else:
            self.statusMessage.emit(
                f"Resized “{label}” on {rid} to {feat.start}–{feat.end}."
            )

    def _modify_feature_at(self, row: int, feat_index: int) -> None:
        from repeatviewer.widgets.sequence_map_viewer import ModifyFeatureDialog

        resolved = self._bars._store_index_for_flat(row, feat_index)
        if resolved is None:
            return
        rid, store_index, _current = resolved
        store = list(self._features.get(rid) or [])
        if not (0 <= store_index < len(store)):
            return
        seq_length = max(1, int(self._lengths.get(rid, 1)))
        parent_candidates = [
            feat.gff_id()
            for i, feat in enumerate(store)
            if i != store_index and feat.gff_id()
        ]
        dialog = ModifyFeatureDialog(
            feature=store[store_index],
            seq_length=seq_length,
            parent_candidates=parent_candidates,
            parent=self,
        )
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        edited = dialog.feature()
        if edited is None:
            QMessageBox.warning(self, "Modify feature", "Enter a feature name and type.")
            return
        if edited.parent and edited.parent == edited.gff_id():
            QMessageBox.warning(
                self,
                "Modify feature",
                "A feature cannot be its own GFF3 Parent.",
            )
            return
        self._push_feature_undo_snapshot()
        old_id = store[store_index].gff_id()
        store[store_index] = edited
        store = remount_parent_ids(store, old_id=old_id, new_id=edited.gff_id())
        self._features[rid] = store
        self._refresh_features_after_edit()
        gff_path = self._persist_overview_features()
        if gff_path is not None:
            self.statusMessage.emit(f"Updated feature on {rid} → {gff_path.name}")
        else:
            self.statusMessage.emit(f"Updated feature on {rid}.")

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        # PDF (or any other Preview page) can hide this panel while the dock
        # resizes; refit + drop the residue cache so MSA fills the viewport.
        QTimer.singleShot(0, self._after_panel_shown)
        QTimer.singleShot(50, self._after_panel_shown)

    def resizeEvent(self, event) -> None:  # noqa: N802
        super().resizeEvent(event)
        self._sync_hbar_shim()
        self._schedule_viewport_refit()
        # Keep −/+/Position/… overflow in sync when the split pane or rail moves.
        bar = self._header_zoom_bar
        hider = getattr(bar, "_overview_overflow", None) if bar is not None else None
        if hider is not None:
            hider.sync()
