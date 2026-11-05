"""Side-panel controls for the genome-browser tracks: pick a track, then set how it
is drawn (expanded / squished, rows, colour by, MAPQ shading). The same settings
are in the track's right-click menu; this is the always-visible version."""

from __future__ import annotations

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QFormLayout,
    QLabel,
    QLineEdit,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends import alignment_style as style
from repeatviewer.backends import genome_tracks as gt

ROW_CHOICES = (4, 8, 20, 50, 100)


class TrackInspector(QWidget):
    settingChanged = Signal(int, str, object)    # track index, key, value

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._tracks: list = []
        self._sig = None
        self._busy = False
        self._picked = None        # the track the user chose, kept across changes
        col = QVBoxLayout(self)
        col.setContentsMargins(0, 0, 0, 0)
        col.setSpacing(6)
        self.track = QComboBox(self)
        self.track.setToolTip("Track to change (right-click a track name for the same options)")
        self.track.currentIndexChanged.connect(self._track_picked)
        col.addWidget(self.track)

        self.form = QFormLayout()
        self.form.setContentsMargins(0, 0, 0, 0)
        self.form.setHorizontalSpacing(8)
        self.form.setVerticalSpacing(6)
        self.form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow)
        self.mode = QComboBox(self)
        self.mode.addItem("Expanded", "expanded")
        self.mode.addItem("Squished", "squished")
        self.mode.addItem("Collapsed", "collapsed")
        self.mode.currentIndexChanged.connect(
            lambda _i: self._emit("mode", self.mode.currentData()))
        self.rows = QComboBox(self)
        for n in ROW_CHOICES:
            self.rows.addItem(str(n), n)
        self.rows.currentIndexChanged.connect(
            lambda _i: self._emit("max_rows", self.rows.currentData()))
        self.color = QComboBox(self)
        for value, label in style.COLOR_MODES:
            self.color.addItem(label.rstrip("…"), value)
        self.color.currentIndexChanged.connect(self._color_changed)
        self.tag = QLineEdit(self)
        self.tag.setPlaceholderText("e.g. HP")
        self.tag.setMaxLength(2)
        self.tag.setToolTip("Two-letter SAM tag whose values give the colours (HP, CB, XS…)")
        self.tag.editingFinished.connect(self._tag_changed)
        self.group = QComboBox(self)
        for value, label in style.GROUP_MODES:
            self.group.addItem(label.rstrip("…"), value)
        self.group.currentIndexChanged.connect(self._group_changed)
        self.group_tag = QLineEdit(self)
        self.group_tag.setPlaceholderText("e.g. HP")
        self.group_tag.setMaxLength(2)
        self.group_tag.editingFinished.connect(self._group_tag_changed)
        self.sort = QComboBox(self)
        self.sort.addItem("Not sorted", "")
        for value, label in style.SORT_MODES:
            if value != style.COLOR_TAG:
                self.sort.addItem(label, value)
        self.sort.setToolTip("Sort the reads over the middle of the view (right-click a "
                             "track to sort at another position)")
        self.sort.activated.connect(
            lambda _i: self._emit("sort_by", self.sort.currentData() or None))
        self.pairs = QCheckBox("View as pairs", self)
        self.pairs.setToolTip("Both reads of a pair in one row, joined by a line")
        self.pairs.toggled.connect(lambda on: self._emit("pairs", bool(on)))
        self.min_mapq = QSpinBox(self)
        self.min_mapq.setRange(0, 255)
        self.min_mapq.setSpecialValueText("show all")
        self.min_mapq.setToolTip("Hide reads with a mapping quality below this")
        self.min_mapq.editingFinished.connect(
            lambda: self._emit("min_mapq", int(self.min_mapq.value())))
        self.hide: dict[str, QCheckBox] = {}
        for key, label, tip in (
                ("hide_duplicates", "Duplicates", "PCR / optical duplicates (flag 0x400)"),
                ("hide_qcfail", "QC-failed", "Reads failing vendor quality checks (0x200)"),
                ("hide_secondary", "Secondary", "Other places of a multi-mapped read (0x100)"),
                ("hide_supplementary", "Supplementary",
                 "Split parts of a chimeric read (0x800)")):
            box = QCheckBox(label, self)
            box.setToolTip("Hide " + tip[0].lower() + tip[1:])
            box.toggled.connect(lambda on, k=key: self._emit(k, bool(on)))
            self.hide[key] = box
        self.shade = QCheckBox("Shade MAPQ 0 reads", self)
        self.shade.setToolTip("Draw reads with mapping quality 0 (placed equally well "
                              "elsewhere) pale, as IGV does")
        self.shade.toggled.connect(lambda on: self._emit("shade_mapq", bool(on)))
        self.quality = QCheckBox("Fade low-quality bases", self)
        self.quality.setToolTip("Mismatches at bases with quality below 20 are drawn "
                                "paler (below 5: almost clear)")
        self.quality.toggled.connect(lambda on: self._emit("shade_quality", bool(on)))
        self.clips = QCheckBox("Mark soft clips ≥ 20 bp", self)
        self.clips.setToolTip("A red edge where a read is clipped: piles of clipped reads "
                              "point at insertion breakpoints and structural variants")
        self.clips.toggled.connect(lambda on: self._emit("flag_clips", bool(on)))
        self._labels = {}
        for key, label, w in (("mode", "Display", self.mode), ("rows", "Rows", self.rows),
                              ("color", "Colour by", self.color), ("tag", "Tag", self.tag),
                              ("group", "Group by", self.group),
                              ("group_tag", "Tag", self.group_tag),
                              ("sort", "Sort by", self.sort),
                              ("min_mapq", "Min MAPQ", self.min_mapq)):
            lab = QLabel(label, self)
            self._labels[key] = lab
            self.form.addRow(lab, w)
        col.addLayout(self.form)
        for w in (self.pairs, self.shade, self.quality, self.clips):
            col.addWidget(w)
        from PySide6.QtWidgets import QGridLayout

        self.hide_label = QLabel("Hide", self)
        self.hide_label.setEnabled(False)
        col.addWidget(self.hide_label)
        grid = QGridLayout()
        grid.setContentsMargins(0, 0, 0, 0)
        grid.setHorizontalSpacing(10)
        grid.setVerticalSpacing(4)
        for k, box in enumerate(self.hide.values()):
            grid.addWidget(box, k // 2, k % 2)
        col.addLayout(grid)
        self.empty = QLabel("No tracks yet: Tracks ▾ > Add track…", self)
        self.empty.setEnabled(False)
        self.empty.setWordWrap(True)
        col.addWidget(self.empty)
        for combo in (self.track, self.mode, self.rows, self.color, self.group, self.sort):
            combo.setSizeAdjustPolicy(
                QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
            combo.setMinimumContentsLength(6)
        self.set_tracks([])

    def label_widgets(self) -> list[QLabel]:
        return list(self._labels.values())

    # ---------------------------------------------------------------- state

    @staticmethod
    def _signature(tracks) -> tuple:
        return tuple((id(t), t.name, tuple(sorted((k, str(v)) for k, v in t.settings.items())))
                     for t in tracks)

    def set_tracks(self, tracks: list) -> None:
        """Show ``tracks``; keeps the chosen track when it is still there."""
        sig = self._signature(tracks)
        if sig == self._sig:
            return
        self._sig = sig
        current = self._picked
        self._tracks = list(tracks)
        self._busy = True
        self.track.clear()
        for t in self._tracks:
            self.track.addItem(t.name)
        if current is not None and current in self._tracks:
            self.track.setCurrentIndex(self._tracks.index(current))
        elif self._tracks:
            # alignment tracks have the most settings: start there
            first = next((k for k, t in enumerate(self._tracks)
                          if t.kind != gt.KIND_ANNOTATION), 0)
            self.track.setCurrentIndex(first)
        self._busy = False
        self._show_track()

    def _track_picked(self, i: int) -> None:
        if not self._busy and 0 <= i < len(self._tracks):
            self._picked = self._tracks[i]
        self._show_track()

    def current(self) -> int:
        return self.track.currentIndex()

    def _show_track(self, *_args) -> None:
        has = bool(self._tracks)
        self.track.setVisible(has)
        self.empty.setVisible(not has)
        i = self.track.currentIndex()
        t = self._tracks[i] if has and 0 <= i < len(self._tracks) else None
        aln = t is not None and t.kind != gt.KIND_ANNOTATION
        for key, w in (("mode", self.mode), ("rows", self.rows)):
            w.setVisible(t is not None)
            self._labels[key].setVisible(t is not None)
        for key, w in (("color", self.color), ("group", self.group), ("sort", self.sort),
                       ("min_mapq", self.min_mapq)):
            w.setVisible(aln)
            self._labels[key].setVisible(aln)
        for w in (self.pairs, self.shade, self.quality, self.clips, self.hide_label,
                  *self.hide.values()):
            w.setVisible(aln)
        if t is None:
            for key, w in (("tag", self.tag), ("group_tag", self.group_tag)):
                w.setVisible(False)
                self._labels[key].setVisible(False)
            return
        self._busy = True
        s = t.settings
        self.mode.setCurrentIndex(max(0, self.mode.findData(s.get("mode", "expanded"))))
        default_rows = 8 if t.kind == gt.KIND_ANNOTATION else 20
        rows = int(s.get("max_rows", default_rows))
        k = self.rows.findData(rows)
        if k < 0:
            self.rows.addItem(str(rows), rows)
            k = self.rows.count() - 1
        self.rows.setCurrentIndex(k)
        mode = s.get("color_by", style.COLOR_STRAND)
        self.color.setCurrentIndex(max(0, self.color.findData(mode)))
        self.tag.setText(s.get("color_tag", ""))
        self.shade.setChecked(bool(s.get("shade_mapq", True)))
        self.pairs.setChecked(bool(s.get("pairs", False)))
        for key, box in self.hide.items():
            box.setChecked(bool(s.get(key, style.FILTER_DEFAULTS[key])))
        self.min_mapq.setValue(int(s.get("min_mapq", 0) or 0))
        # "Collapsed" is for reads only
        coll = self.mode.findData("collapsed")
        self.mode.model().item(coll).setEnabled(aln)
        group = s.get("group_by", style.GROUP_NONE)
        self.group.setCurrentIndex(max(0, self.group.findData(group)))
        self.group_tag.setText(s.get("group_tag", ""))
        self.sort.setCurrentIndex(max(0, self.sort.findData((s.get("sort") or {}).get("by", ""))))
        self.quality.setChecked(bool(s.get("shade_quality", True)))
        self.clips.setChecked(bool(s.get("flag_clips", True)))
        self._busy = False
        show_tag = aln and mode == style.COLOR_TAG
        self.tag.setVisible(show_tag)
        self._labels["tag"].setVisible(show_tag)
        show_gtag = aln and group == style.COLOR_TAG
        self.group_tag.setVisible(show_gtag)
        self._labels["group_tag"].setVisible(show_gtag)

    # ---------------------------------------------------------------- edits

    def _emit(self, key: str, value) -> None:
        if self._busy or not (0 <= self.current() < len(self._tracks)):
            return
        self.settingChanged.emit(self.current(), key, value)

    def _color_changed(self, _i: int) -> None:
        mode = self.color.currentData()
        show_tag = mode == style.COLOR_TAG
        self.tag.setVisible(show_tag)
        self._labels["tag"].setVisible(show_tag)
        if show_tag and not self.tag.text().strip():
            self.tag.setText("HP")
            self._emit("color_tag", "HP")
        self._emit("color_by", mode)
        if show_tag and not self._busy:
            self.tag.setFocus(Qt.FocusReason.OtherFocusReason)
            self.tag.selectAll()

    def _group_changed(self, _i: int) -> None:
        mode = self.group.currentData()
        show_tag = mode == style.COLOR_TAG
        self.group_tag.setVisible(show_tag)
        self._labels["group_tag"].setVisible(show_tag)
        if show_tag and not self.group_tag.text().strip():
            self.group_tag.setText("HP")
            self._emit("group_tag", "HP")
        self._emit("group_by", mode)

    def _group_tag_changed(self) -> None:
        name = self.group_tag.text().strip()
        if len(name) == 2:
            self._emit("group_tag", name)

    def _tag_changed(self) -> None:
        name = self.tag.text().strip()
        if len(name) == 2:
            self._emit("color_tag", name)
