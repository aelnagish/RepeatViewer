"""Primers page of the Annotation view's side rail (Preview > Annotation > Primers).

The page lists the primers kept for the open FASTA, shows where each binds on
the shown sequence (arrows in the annotation track, mismatches in red on the
bases) and the numbers of the selected primer. It talks to its host, the
:class:`SequenceMapPanel`, through a few ``primers_host_*`` methods.
"""

from __future__ import annotations

import html
import math
from pathlib import Path

from PySide6.QtCore import Qt, QTimer, Signal
from PySide6.QtGui import QColor, QFont, QFontDatabase
from PySide6.QtWidgets import (
    QAbstractItemView,
    QApplication,
    QDoubleSpinBox,
    QGridLayout,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QMenu,
    QMessageBox,
    QPushButton,
    QSizePolicy,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends.primers import on_sequence as osq
from repeatviewer.backends.primers import thermo as th
from repeatviewer.backends.primers.oligos import Primer
from repeatviewer.widgets.help_button import make_help_button

FORWARD_COLOR = osq.COLOR_FORWARD
REVERSE_COLOR = osq.COLOR_REVERSE
MISMATCH_COLOR = "#d0342c"
OK_COLOR = "#1a7f4b"
CHECK_COLOR = "#b26b00"

PRIMERS_HELP = (
    "<b>Primers</b> lists the primers kept with this FASTA. Each one is drawn as an arrow "
    "where it binds (blue: points right, orange: points left); mismatches are red on the "
    "bases.<br>"
    "<b>Tick</b> two primers and press <b>Run PCR</b> to see their products.<br>"
    "<b>Click</b> a primer to see its numbers and jump to it; <b>double-click</b> opens "
    "Edit Primer.<br>"
    "A <b>tail</b> is the 5′ part that does not bind (lower case in the sequence); it is "
    "copied into the product from cycle 3 on, so two Tm values are shown: annealed and full.<br>"
    "Primers are saved in a hidden file next to the FASTA (<b>.primers_RV.json</b>)."
)


def mono_font(point_delta: int = 0) -> QFont:
    font = QFontDatabase.systemFont(QFontDatabase.SystemFont.FixedFont)
    if point_delta:
        font.setPointSize(max(8, font.pointSize() + point_delta))
    return font


def fmt_tm(value: float) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "–"
    return f"{value:.1f} °C"


def primer_rich_text(seq: str, placement: osq.Placement | None = None) -> str:
    """5'→3' oligo as rich text: tail blue underlined, mismatches red."""
    tail, bind = th.split_tail(seq)
    mm_idx: set[int] = set()
    if placement is not None and placement.mismatch_positions:
        n = len(bind)
        for pos in placement.mismatch_positions:
            if placement.strand == "+":
                mm_idx.add(pos - placement.start)
            else:
                mm_idx.add(placement.end - pos)
        mm_idx = {i for i in mm_idx if 0 <= i < n}
    out = []
    if tail:
        out.append(f'<span style="color:{FORWARD_COLOR};text-decoration:underline">'
                   f"{html.escape(tail.upper())}</span>")
    for i, ch in enumerate(bind):
        if i in mm_idx:
            out.append(f'<span style="color:#ffffff;background:{MISMATCH_COLOR}">{ch}</span>')
        else:
            out.append(ch)
    return "".join(out)


def status_html(status: str) -> str:
    color = {"ok": OK_COLOR, "check": CHECK_COLOR}.get(status, MISMATCH_COLOR)
    return f'<span style="color:{color}">{status}</span>'


def find_main_window():
    """The RepeatViewer main window (previews can live in their own windows)."""
    app = QApplication.instance()
    if app is None:
        return None
    for w in app.topLevelWidgets():
        if hasattr(w, "_open_primer_design_dialog"):
            return w
    return None


class _Tile(QWidget):
    """Small caption + value box (Tm, GC …)."""

    def __init__(self, caption: str, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("primerTile")
        self.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)
        self.setStyleSheet(
            "QWidget#primerTile { background: palette(base); border: 1px solid palette(mid);"
            " border-radius: 6px; }")
        lay = QVBoxLayout(self)
        lay.setContentsMargins(8, 4, 8, 5)
        lay.setSpacing(0)
        cap = QLabel(caption, self)
        f = cap.font()
        f.setPointSize(max(8, f.pointSize() - 2))
        cap.setFont(f)
        cap.setStyleSheet("color: palette(placeholder-text); border: none;")
        self.value = QLabel("–", self)
        vf = self.value.font()
        vf.setPointSize(vf.pointSize() + 2)
        vf.setBold(True)
        self.value.setFont(vf)
        self.value.setStyleSheet("border: none;")
        lay.addWidget(cap)
        lay.addWidget(self.value)


class PrimersPanel(QWidget):
    """Rail page body. ``set_host`` must be called before use."""

    primersChanged = Signal()

    COL_TICK, COL_NAME, COL_TM, COL_SITE = range(4)

    def __init__(self, parent: QWidget | None = None, *, width: int = 260) -> None:
        super().__init__(parent)
        self._host = None
        self._primers: list[Primer] = []
        self._placements: dict[str, osq.Placement | None] = {}
        self._ticked: set[str] = set()
        self._selected: str = ""
        self._path: Path | None = None
        self._seq_id = ""
        self._template = ""
        self._loading = False
        self._report_timer = QTimer(self)
        self._report_timer.setSingleShot(True)
        self._report_timer.setInterval(120)
        self._report_timer.timeout.connect(self._refresh_inspector)

        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(8)

        grid = QGridLayout()
        grid.setContentsMargins(0, 0, 0, 0)
        grid.setHorizontalSpacing(6)
        grid.setVerticalSpacing(6)
        self.add_button = QPushButton("Add Primer…", self)
        self.database_button = QPushButton("From Database…", self)
        self.design_button = QPushButton("Design…", self)
        self.pcr_button = QPushButton("Run PCR", self)
        for i, b in enumerate((self.add_button, self.database_button,
                               self.design_button, self.pcr_button)):
            b.setMinimumWidth(0)
            b.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Fixed)
            grid.addWidget(b, i // 2, i % 2)
        self.pcr_button.setDefault(False)
        self.pcr_button.setToolTip("PCR with the ticked primers on this file")
        lay.addLayout(grid)

        self.table = QTableWidget(0, 4, self)
        self.table.setHorizontalHeaderLabels(["", "Primer", "Tm", "Site"])
        self.table.verticalHeader().hide()
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setShowGrid(False)
        self.table.setAlternatingRowColors(True)
        self.table.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        head = self.table.horizontalHeader()
        head.setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
        head.setDefaultAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        head.setStretchLastSection(True)
        self.table.setColumnWidth(self.COL_TICK, 24)
        self.table.setColumnWidth(self.COL_NAME, 96)
        self.table.setColumnWidth(self.COL_TM, 40)
        self.table.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        lay.addWidget(self.table)

        self.empty_label = QLabel(
            "No primers yet. Add one, pick saved primers from the Database, or design new ones.",
            self)
        self.empty_label.setWordWrap(True)
        self.empty_label.setStyleSheet("color: palette(placeholder-text);")
        lay.addWidget(self.empty_label)

        # ---- inspector
        self.inspector = QWidget(self)
        il = QVBoxLayout(self.inspector)
        il.setContentsMargins(0, 4, 0, 0)
        il.setSpacing(6)
        head_row = QHBoxLayout()
        head_row.setContentsMargins(0, 0, 0, 0)
        self.name_label = QLabel(self.inspector)
        nf = self.name_label.font()
        nf.setBold(True)
        nf.setPointSize(nf.pointSize() + 1)
        self.name_label.setFont(nf)
        head_row.addWidget(self.name_label, 1)
        head_row.addWidget(make_help_button(self.inspector, "Primers", PRIMERS_HELP))
        il.addLayout(head_row)
        self.seq_label = QLabel(self.inspector)
        self.seq_label.setFont(mono_font(-1))
        self.seq_label.setWordWrap(True)
        self.seq_label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        self.seq_label.setTextFormat(Qt.TextFormat.RichText)
        self.seq_label.setStyleSheet(
            "QLabel { background: palette(base); border: 1px solid palette(mid);"
            " border-radius: 6px; padding: 6px; }")
        il.addWidget(self.seq_label)
        tiles = QGridLayout()
        tiles.setContentsMargins(0, 0, 0, 0)
        tiles.setSpacing(6)
        self.tile_tm = _Tile("Tm annealed", self.inspector)
        self.tile_tm_full = _Tile("Tm full", self.inspector)
        self.tile_gc = _Tile("GC", self.inspector)
        self.tile_end = _Tile("3′ end exact", self.inspector)
        for i, t in enumerate((self.tile_tm, self.tile_tm_full, self.tile_gc, self.tile_end)):
            tiles.addWidget(t, i // 2, i % 2)
        il.addLayout(tiles)
        self.checks_label = QLabel(self.inspector)
        self.checks_label.setTextFormat(Qt.TextFormat.RichText)
        self.checks_label.setWordWrap(True)
        il.addWidget(self.checks_label)
        self.where_label = QLabel(self.inspector)
        self.where_label.setTextFormat(Qt.TextFormat.RichText)
        self.where_label.setWordWrap(True)
        il.addWidget(self.where_label)

        snap_row = QHBoxLayout()
        snap_row.setContentsMargins(0, 0, 0, 0)
        snap_row.setSpacing(6)
        self.snap_button = QPushButton("Snap to Tm", self.inspector)
        self.snap_tm = QDoubleSpinBox(self.inspector)
        self.snap_tm.setRange(40.0, 80.0)
        self.snap_tm.setDecimals(0)
        self.snap_tm.setValue(60.0)
        self.snap_tm.setSuffix(" °C")
        snap_row.addWidget(self.snap_button, 1)
        snap_row.addWidget(self.snap_tm)
        il.addLayout(snap_row)
        self.fix_button = QPushButton("Fix Mismatches to Match Template", self.inspector)
        self.edit_button = QPushButton("Edit Primer…", self.inspector)
        self.genome_button = QPushButton("Sites in Genome…", self.inspector)
        for b in (self.fix_button, self.edit_button, self.genome_button):
            b.setMinimumWidth(0)
            b.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Fixed)
            il.addWidget(b)
        lay.addWidget(self.inspector)
        lay.addStretch(1)

        self.add_button.clicked.connect(self._on_add)
        self.database_button.clicked.connect(self._on_from_database)
        self.design_button.clicked.connect(self._on_design)
        self.pcr_button.clicked.connect(self._on_pcr)
        self.edit_button.clicked.connect(lambda: self.edit_primer(self._selected))
        self.fix_button.clicked.connect(self._on_fix)
        self.snap_button.clicked.connect(self._on_snap)
        self.genome_button.clicked.connect(self._on_genome)
        self.table.itemSelectionChanged.connect(self._on_selection)
        self.table.itemChanged.connect(self._on_item_changed)
        self.table.cellDoubleClicked.connect(self._on_double_click)
        self.table.customContextMenuRequested.connect(self._on_context_menu)
        self._sync_visible()

    # ------------------------------------------------------------------ host

    def set_host(self, host) -> None:
        self._host = host

    def primers(self) -> list[Primer]:
        return list(self._primers)

    def placement(self, name: str) -> osq.Placement | None:
        return self._placements.get(name)

    def ticked_primers(self) -> list[Primer]:
        return [p for p in self._primers if p.name in self._ticked]

    def fasta_path(self) -> Path | None:
        return self._path

    def reload(self) -> None:
        """Read the primers of the host's file and place them on the shown sequence."""
        host = self._host
        if host is None:
            return
        try:
            path = Path(host.primers_host_path())
            seq_id, template = host.primers_host_sequence()
            features = list(host.primers_host_features())
        except Exception:  # noqa: BLE001 - no data loaded yet
            return
        same_file = self._path is not None and path == self._path
        self._path = path
        self._seq_id = seq_id
        self._template = template or ""
        stored = osq.load_sequence_primers(path)
        legacy = osq.primers_from_features(features, self._template, stored)
        self._primers = stored + legacy
        if not same_file:
            self._ticked = set()
            self._selected = ""
        self._place_all()
        self._fill_table()
        self._push_marks()

    def _place_all(self) -> None:
        self._placements = {}
        for p in self._primers:
            try:
                self._placements[p.name] = osq.place_primer(p.sequence, self._template,
                                                            name=p.name)
            except Exception:  # noqa: BLE001 - odd oligo: not placed
                self._placements[p.name] = None

    # ----------------------------------------------------------------- table

    def _fill_table(self) -> None:
        self._loading = True
        try:
            self.table.setRowCount(len(self._primers))
            for row, p in enumerate(self._primers):
                place = self._placements.get(p.name)
                tick = QTableWidgetItem()
                tick.setFlags(Qt.ItemFlag.ItemIsUserCheckable | Qt.ItemFlag.ItemIsEnabled
                              | Qt.ItemFlag.ItemIsSelectable)
                tick.setCheckState(Qt.CheckState.Checked if p.name in self._ticked
                                   else Qt.CheckState.Unchecked)
                self.table.setItem(row, self.COL_TICK, tick)
                name = QTableWidgetItem(p.name)
                name.setData(Qt.ItemDataRole.UserRole, p.name)
                if place is not None:
                    name.setForeground(QColor(FORWARD_COLOR if place.strand == "+"
                                              else REVERSE_COLOR))
                    f = name.font()
                    f.setBold(True)
                    name.setFont(f)
                name.setToolTip(p.sequence)
                self.table.setItem(row, self.COL_NAME, name)
                try:
                    tm = th.calc_tm(p.binding)
                except Exception:  # noqa: BLE001
                    tm = float("nan")
                self.table.setItem(row, self.COL_TM, QTableWidgetItem(
                    "–" if math.isnan(tm) else f"{tm:.1f}"))
                if place is None:
                    site = "not on this sequence"
                else:
                    arrow = "→" if place.strand == "+" else "←"
                    site = f"{arrow} {place.start:,}"
                    if place.mismatches:
                        site += f" · {place.mismatches} mm"
                self.table.setItem(row, self.COL_SITE, QTableWidgetItem(site))
            target = self._selected if self._selected in {p.name for p in self._primers} else ""
            if not target and self._primers:
                target = self._primers[0].name
            self._select_row_of(target)
        finally:
            self._loading = False
        rows = max(2, min(8, len(self._primers)))
        row_h = self.table.verticalHeader().defaultSectionSize()
        self.table.setFixedHeight(self.table.horizontalHeader().sizeHint().height()
                                  + rows * row_h + 4)
        self._selected = self._current_name()
        self._sync_visible()
        self._report_timer.start()

    def _select_row_of(self, name: str) -> None:
        for row in range(self.table.rowCount()):
            item = self.table.item(row, self.COL_NAME)
            if item is not None and item.data(Qt.ItemDataRole.UserRole) == name:
                self.table.selectRow(row)
                return

    def _current_name(self) -> str:
        rows = sorted({i.row() for i in self.table.selectedIndexes()})
        if not rows:
            return ""
        item = self.table.item(rows[0], self.COL_NAME)
        return str(item.data(Qt.ItemDataRole.UserRole)) if item else ""

    def _selected_names(self) -> list[str]:
        rows = sorted({i.row() for i in self.table.selectedIndexes()})
        out = []
        for r in rows:
            item = self.table.item(r, self.COL_NAME)
            if item is not None:
                out.append(str(item.data(Qt.ItemDataRole.UserRole)))
        return out

    def _sync_visible(self) -> None:
        has = bool(self._primers)
        self.table.setVisible(has)
        self.empty_label.setVisible(not has)
        self.inspector.setVisible(has and bool(self._selected))
        self.pcr_button.setEnabled(len(self._ticked) >= 1)

    def _on_item_changed(self, item: QTableWidgetItem) -> None:
        if self._loading or item.column() != self.COL_TICK:
            return
        name_item = self.table.item(item.row(), self.COL_NAME)
        if name_item is None:
            return
        name = str(name_item.data(Qt.ItemDataRole.UserRole))
        if item.checkState() == Qt.CheckState.Checked:
            self._ticked.add(name)
        else:
            self._ticked.discard(name)
        self._sync_visible()
        self._report_timer.start()

    def _on_selection(self) -> None:
        if self._loading:
            return
        self._selected = self._current_name()
        self._sync_visible()
        self._report_timer.start()
        place = self._placements.get(self._selected)
        if place is not None and self._host is not None:
            try:
                self._host.primers_host_show(place.start, place.end)
            except Exception:  # noqa: BLE001
                pass

    def _on_double_click(self, row: int, col: int) -> None:
        if col == self.COL_TICK:
            return
        item = self.table.item(row, self.COL_NAME)
        if item is not None:
            self.edit_primer(str(item.data(Qt.ItemDataRole.UserRole)))

    def _on_context_menu(self, pos) -> None:
        names = self._selected_names()
        if not names:
            return
        menu = QMenu(self)
        menu.addAction("Edit Primer…", lambda: self.edit_primer(names[0]))
        menu.addAction("Copy Sequence", lambda: self._copy_sequences(names))
        menu.addAction("Add to Primer Database", lambda: self._add_to_database(names))
        menu.addSeparator()
        label = "Remove Primer" if len(names) == 1 else f"Remove {len(names)} Primers"
        menu.addAction(label, lambda: self.remove_primers(names))
        menu.exec(self.table.viewport().mapToGlobal(pos))

    # -------------------------------------------------------------- inspector

    def _primer(self, name: str) -> Primer | None:
        for p in self._primers:
            if p.name == name:
                return p
        return None

    def _refresh_inspector(self) -> None:
        p = self._primer(self._selected)
        if p is None:
            self._sync_visible()
            return
        others = [o for o in self._primers if o.name in self._ticked and o.name != p.name]
        rep = osq.analyse_primer(p, self._template, others=others)
        place = rep.placement
        color = FORWARD_COLOR if place is None or place.strand == "+" else REVERSE_COLOR
        self.name_label.setText(f'<span style="color:{color}">{html.escape(p.name)}</span>')
        self.seq_label.setText(
            '<span style="color:gray">5′ </span>' + primer_rich_text(p.sequence, place)
            + '<span style="color:gray"> 3′</span>')
        self.tile_tm.value.setText(fmt_tm(rep.tm_annealed))
        self.tile_tm_full.value.setText(fmt_tm(rep.tm_full) if rep.tail else "same")
        self.tile_gc.value.setText(f"{rep.gc:.0f} %")
        if place is None:
            self.tile_end.value.setText("–")
        else:
            n = min(place.three_exact, 99)
            col = OK_COLOR if n >= 5 else (CHECK_COLOR if n >= 3 else MISMATCH_COLOR)
            self.tile_end.value.setText(f'<span style="color:{col}">{n}</span>')
        rows = []
        hp = rep.hairpin
        rows.append(("Hairpin", hp.dg, osq.dg_status(hp.dg, limit=th.WARN_HAIRPIN_DG)))
        sd = rep.self_dimer
        rows.append(("Self-dimer", sd.dg, osq.dg_status(sd.dg)))
        for other, d in rep.cross.items():
            rows.append((f"Dimer with {other}", d.dg, osq.dg_status(d.dg)))
        lines = []
        for label, dg, st in rows:
            val = f"ΔG {dg:.1f}" if dg else "none"
            lines.append(f"{html.escape(label)}: {val} {status_html(st)}")
        self.checks_label.setText("<br>".join(lines))
        if place is None:
            where = "Does not bind this sequence (up to 4 mismatches)."
        else:
            mm = ("no mismatch" if not place.mismatches else
                  f"{place.mismatches} mismatch{'es' if place.mismatches > 1 else ''}")
            where = (f"<b>Binds</b> {place.start:,}..{place.end:,} "
                     f"({'→' if place.strand == '+' else '←'}), {mm}, "
                     f"{place.annealed} annealed")
            if place.tail:
                where += f" + {place.tail} tail"
            if place.other_sites:
                where += (f'<br><span style="color:{CHECK_COLOR}">{place.other_sites} more '
                          "site(s) here with ≤ 2 mismatches</span>")
        self.where_label.setText(where)
        self.fix_button.setEnabled(place is not None and place.mismatches > 0)
        self.snap_button.setEnabled(place is not None)
        self._push_marks()

    def _push_marks(self) -> None:
        if self._host is None:
            return
        marks: list[int] = []
        for p in self._primers:
            place = self._placements.get(p.name)
            if place is not None:
                marks.extend(place.mismatch_positions)
        try:
            self._host.primers_host_set_marks(sorted(set(marks)))
        except Exception:  # noqa: BLE001
            pass

    # ---------------------------------------------------------------- changes

    def set_primers(self, primers: list[Primer], *, select: str = "") -> None:
        """Replace the list, save it and redraw the arrows."""
        self._primers = list(primers)
        names = {p.name for p in self._primers}
        self._ticked &= names
        if select:
            self._selected = select
        self._place_all()
        self._save()
        self._fill_table()
        self._push_marks()
        self.primersChanged.emit()

    def add_primers(self, primers: list[Primer], *, tick: bool = False) -> list[Primer]:
        out = list(self._primers)
        added = []
        for p in primers:
            name = osq.new_primer_name(p.name, out)
            seq = p.sequence
            if self._template and not th.split_tail(seq)[0]:
                try:
                    seq = osq.guess_tail(seq, self._template)
                except Exception:  # noqa: BLE001
                    pass
            np_ = Primer(name=name, sequence=seq, kind=p.kind, note=p.note, source=p.source,
                         extra=dict(p.extra))
            out.append(np_)
            added.append(np_)
            if tick:
                self._ticked.add(name)
        self.set_primers(out, select=added[-1].name if added else "")
        return added

    def replace_primer(self, old_name: str, new: Primer) -> None:
        out = [new if p.name == old_name else p for p in self._primers]
        if old_name in self._ticked:
            self._ticked.discard(old_name)
            self._ticked.add(new.name)
        self.set_primers(out, select=new.name)

    def remove_primers(self, names: list[str]) -> None:
        drop = set(names)
        self.set_primers([p for p in self._primers if p.name not in drop])

    def _save(self) -> None:
        if self._path is None:
            return
        try:
            osq.save_sequence_primers(self._path, self._primers)
        except OSError as exc:
            QMessageBox.warning(self, "Primers", f"Could not save the primers:\n{exc}")
            return
        if self._host is None:
            return
        try:
            features = list(self._host.primers_host_features())
            merged = osq.merge_primer_features(features, self._primers, self._placements)
            self._host.primers_host_apply_features(merged)
        except Exception as exc:  # noqa: BLE001
            QMessageBox.warning(self, "Primers", f"Could not update the arrows:\n{exc}")

    # ---------------------------------------------------------------- actions

    def _on_add(self) -> None:
        from repeatviewer.views.primer_edit_dialog import EditPrimerDialog

        name = osq.new_primer_name("primer_F", self._primers)
        dlg = EditPrimerDialog(Primer(name=name, sequence="N"), template=self._template,
                               others=self._primers, new=True, parent=self.window())
        if dlg.exec():
            self.add_primers([dlg.primer()])
            if dlg.save_to_database():
                self._add_to_database([dlg.primer().name])

    def edit_primer(self, name: str) -> None:
        from repeatviewer.views.primer_edit_dialog import EditPrimerDialog

        p = self._primer(name)
        if p is None:
            return
        dlg = EditPrimerDialog(p, template=self._template,
                               others=[o for o in self._primers if o.name != name],
                               parent=self.window())
        if dlg.exec():
            new = dlg.primer()
            if new.name != name and new.name in {o.name for o in self._primers}:
                new = Primer(name=osq.new_primer_name(new.name, self._primers),
                             sequence=new.sequence, kind=new.kind, note=new.note,
                             source=new.source, extra=new.extra)
            self.replace_primer(name, new)
            if dlg.save_to_database():
                self._add_to_database([new.name])

    def _on_from_database(self) -> None:
        from repeatviewer.views.primer_edit_dialog import PickPrimersDialog

        dlg = PickPrimersDialog(self._database_path(), template=self._template,
                                parent=self.window())
        if dlg.exec():
            picked = dlg.primers()
            if picked:
                self.add_primers(picked, tick=len(picked) <= 2)

    def _on_design(self) -> None:
        mw = find_main_window()
        if mw is None or self._path is None:
            return
        try:
            mw._open_primer_design_dialog(files=[self._path])
        except TypeError:
            mw._open_primer_design_dialog()

    def _on_pcr(self) -> None:
        mw = find_main_window()
        if mw is None or self._path is None:
            return
        opener = getattr(mw, "_open_pcr_window", None)
        if opener is not None:
            opener(templates=[self._path], primers=self.ticked_primers())

    def _on_genome(self) -> None:
        mw = find_main_window()
        p = self._primer(self._selected)
        if mw is None or p is None:
            return
        prims = self.ticked_primers() or [p]
        if p not in prims:
            prims = [p] + prims
        opener = getattr(mw, "_open_primer_genome_view", None)
        if opener is not None:
            opener(prims)

    def _on_fix(self) -> None:
        p = self._primer(self._selected)
        place = self._placements.get(self._selected)
        if p is None or place is None:
            return
        self.replace_primer(p.name, Primer(name=p.name, sequence=osq.fix_mismatches(p.sequence, place),
                                           kind=p.kind, note=p.note, source=p.source,
                                           extra=p.extra))

    def _on_snap(self) -> None:
        p = self._primer(self._selected)
        place = self._placements.get(self._selected)
        if p is None or place is None:
            return
        seq = osq.snap_to_tm(p.sequence, place, self._template, self.snap_tm.value())
        self.replace_primer(p.name, Primer(name=p.name, sequence=seq, kind=p.kind, note=p.note,
                                           source=p.source, extra=p.extra))

    def _copy_sequences(self, names: list[str]) -> None:
        text = "\n".join(f">{p.name}\n{p.sequence}" for p in self._primers if p.name in names)
        QApplication.clipboard().setText(text)

    def _database_path(self) -> Path | None:
        try:
            from repeatviewer.job_options import load_primer_options

            raw = load_primer_options().database_path
            return Path(raw) if raw else None
        except Exception:  # noqa: BLE001
            return None

    def _add_to_database(self, names: list[str]) -> None:
        try:
            from repeatviewer.backends.primers.inventory import open_database

            db = open_database(self._database_path())
            items = [Primer(name=p.name, sequence=p.sequence, kind=p.kind, note=p.note,
                            source=p.source or (self._path.name if self._path else ""),
                            extra={"status": "designed"})
                     for p in self._primers if p.name in names]
            added = db.add(items)
            db.save()
        except Exception as exc:  # noqa: BLE001
            QMessageBox.warning(self, "Primer Database", f"Could not add the primers:\n{exc}")
            return
        mw = find_main_window()
        msg = (f"Added {len(added)} primer(s) to the Primer Database" if added
               else "These primers are already in the Primer Database")
        if mw is not None and hasattr(mw, "_set_status_message"):
            mw._set_status_message(msg, 6000)
