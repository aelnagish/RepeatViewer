"""PCR tab of the primer workbench (Primers rail): live in-silico PCR on the open file.

Pick the forward and reverse primer of this sequence; products update by
themselves (no Run) and the chosen one is drawn on the map. Status box green
when the expected product fits the Design range, amber when it does not
(allowed: the tube still makes it). Bench Card, Open in a Gel and Save
Products use the same code as the PCR window.

Developer: Abdullah Yusuf
"""

from __future__ import annotations

import threading
from pathlib import Path

from PySide6.QtCore import QObject, QSettings, Qt, QTimer, Signal
from PySide6.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QFrame,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends.primers import pcr
from repeatviewer.backends.primers import thermo as th
from repeatviewer.widgets.primer_design_pane import (
    FIELD_W,
    GAP,
    LABEL_W,
    _row,
    _rule,
    _section,
)

_SETTINGS_KEY = "primers/workbench_pcr"
# The PCR rows have longer labels: wider label column, the same right edge.
# Older name for the shared label column, kept so nothing outside breaks.
PCR_LABEL_W = LABEL_W

RUN_THIS = "this"
RUN_FILES = "files"
RUN_ON_HELP = (
    "<b>This one</b>: the open sequence, live, as you change anything.<br>"
    "<b>Selected files</b>: every FASTA ticked in the Files panel, so one pair can be tested "
    "on many satellite families at once. The products of all of them are listed together, "
    "each with the file it came from."
)
PCR_HELP = (
    "<b>Forward</b>, <b>Reverse</b>: primers of this sequence (Primer page).<br>"
    "<b>Kind</b>: Standard lists every product; Satellite ladder groups products that step "
    "by a monomer (1×, 2×, 3×).<br>"
    "<b>Mismatches</b> and <b>3′ end exact</b> decide which sites count. "
    "<b>Polymerase</b> sets how fast long products fall off.<br>"
    "Products update by themselves when anything changes. The status box is amber when the "
    "expected product is outside the product range of the Design tab: allowed, the tube "
    "still makes it."
)
PRODUCTS_HELP = (
    "Every product the two primers make on this file, best first. Click one to see it on "
    "the sequence. <b>mm F / R</b>: mismatches of the forward and reverse primer at that site."
)


class _Relay(QObject):
    done = Signal(int, object, object)   # generation, PcrResult | None, error | None


def _combo(items, parent: QWidget, current=None) -> QComboBox:
    box = QComboBox(parent)
    for label, data in items:
        box.addItem(label, data)
    if current is not None:
        i = box.findData(current)
        box.setCurrentIndex(max(0, i))
    box.setFixedWidth(FIELD_W)
    return box


class PcrPane(QWidget):
    """``panel`` is the PrimersPanel (primers, host, Design settings)."""

    def __init__(self, panel, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._panel = panel
        self._result: pcr.PcrResult | None = None
        self._generation = 0
        self.active = False
        self._relay = _Relay(self)
        self._relay.done.connect(self._on_done)
        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(300)
        self._timer.timeout.connect(self._run)
        s = self._load_settings()

        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(10)

        sec, sl, _ = _section("Run on", self, "Run on", RUN_ON_HELP)
        from PySide6.QtWidgets import QButtonGroup

        self.run_group = QButtonGroup(self)
        self.run_group.setExclusive(True)
        buttons = []
        for i, (key, text) in enumerate(((RUN_THIS, "This one"), (RUN_FILES, "Selected files"))):
            b = QPushButton(text, sec)
            b.setCheckable(True)
            b.setProperty("run", key)
            b.setFixedWidth((FIELD_W - GAP) // 2)
            self.run_group.addButton(b, i)
            buttons.append(b)
        buttons[0].setChecked(True)
        sl.addLayout(_row("Sequence", *buttons, parent=sec, label_w=PCR_LABEL_W))
        self.run_note = QLabel("", sec)
        self.run_note.setWordWrap(True)
        self.run_note.setIndent(PCR_LABEL_W + GAP)
        self.run_note.setStyleSheet("color: palette(placeholder-text);")
        sl.addWidget(self.run_note)
        lay.addWidget(sec)
        lay.addWidget(_rule(self))

        sec, sl, _ = _section("Primers", self, "PCR", PCR_HELP)
        self.forward = _combo([], sec)
        self.reverse = _combo([], sec)
        sl.addLayout(_row("Forward", self.forward, parent=sec, label_w=PCR_LABEL_W))
        sl.addLayout(_row("Reverse", self.reverse, parent=sec, label_w=PCR_LABEL_W))
        self.kind = _combo([(pcr.MODE_LABELS[m], m) for m in pcr.MODES], sec, s.mode)
        sl.addLayout(_row("Kind", self.kind, parent=sec, label_w=PCR_LABEL_W))
        self.mismatches = _combo([("none" if n == 0 else f"up to {n}", n) for n in range(5)],
                                 sec, s.max_mismatches)
        sl.addLayout(_row("Mismatches", self.mismatches, parent=sec, label_w=PCR_LABEL_W))
        self.three = _combo([("no rule" if n == 0 else f"last {n} bases", n) for n in (0, 2, 3, 4, 5)],
                            sec, s.three_prime_exact)
        sl.addLayout(_row("3′ end exact", self.three, parent=sec, label_w=PCR_LABEL_W))
        self.polymerase = _combo([(k, k) for k in pcr.POLYMERASES], sec, s.polymerase)
        sl.addLayout(_row("Polymerase", self.polymerase, parent=sec, label_w=PCR_LABEL_W))
        self.status = QFrame(sec)
        self.status.setObjectName("workbenchPcrStatus")
        stl = QVBoxLayout(self.status)
        stl.setContentsMargins(8, 6, 8, 6)
        stl.setSpacing(6)
        self.status_label = QLabel(self.status)
        self.status_label.setWordWrap(True)
        self.status_label.setTextFormat(Qt.TextFormat.RichText)
        stl.addWidget(self.status_label)
        self.fit_button = QPushButton(self.status)
        fr = QHBoxLayout()
        fr.addStretch(1)
        fr.addWidget(self.fit_button)
        stl.addLayout(fr)
        sl.addWidget(self.status)
        lay.addWidget(sec)
        lay.addWidget(_rule(self))

        sec, sl, self.products_title = _section("Products", self, "Products", PRODUCTS_HELP)
        self.table = QTableWidget(0, 4, sec)
        self.table.setHorizontalHeaderLabels(["", "Range", "Size", "mm F / R"])
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
        for c, w in enumerate((30, 96, 70)):
            self.table.setColumnWidth(c, w)
        self.table.setMinimumHeight(130)
        self.table.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        sl.addWidget(self.table)
        lay.addWidget(sec)
        lay.addWidget(_rule(self))

        sec, sl, self.gel_title = _section("Expected gel", self, "Expected gel",
                                           "The products of this PCR next to a size ladder. "
                                           "<b>Open in a Gel</b> makes a gel document you can "
                                           "change (lanes, agarose, ladder).")
        from repeatviewer.widgets.gel_view import GelView

        self.gel = GelView(sec, compact=True)
        self.gel.setMinimumHeight(170)
        sl.addWidget(self.gel)
        self.bench_button = QPushButton("Bench Card", sec)
        self.gel_button = QPushButton("Open in a Gel", sec)
        self.save_button = QPushButton("Save Products", sec)
        br = QHBoxLayout()
        br.addStretch(1)
        for b in (self.bench_button, self.gel_button, self.save_button):
            br.addWidget(b)
        sl.addLayout(br)
        lay.addWidget(sec)
        lay.addStretch(1)

        for box in (self.forward, self.reverse, self.kind, self.mismatches, self.three,
                    self.polymerase):
            box.currentIndexChanged.connect(self.schedule)
        self.table.itemSelectionChanged.connect(self.push_overlay)
        self.run_group.idClicked.connect(lambda _i: self._on_run_where())
        self.fit_button.clicked.connect(self._on_fit)
        self.bench_button.clicked.connect(self._on_bench)
        self.gel_button.clicked.connect(self._on_gel)
        self.save_button.clicked.connect(self._on_save)
        self._show_status("", None)
        self._sync_buttons()

    # ------------------------------------------------------------------ inputs

    def _load_settings(self) -> pcr.PcrSettings:
        import json

        raw = QSettings().value(_SETTINGS_KEY, "")
        try:
            return pcr.PcrSettings.from_dict(json.loads(raw) if raw else None)
        except (TypeError, ValueError):
            return pcr.PcrSettings()

    def settings(self) -> pcr.PcrSettings:
        return pcr.PcrSettings(mode=str(self.kind.currentData()),
                               max_mismatches=int(self.mismatches.currentData()),
                               three_prime_exact=int(self.three.currentData()),
                               polymerase=str(self.polymerase.currentData()))

    def refresh_primers(self) -> None:
        """Fill Forward / Reverse from the sequence's primers (keeps the choice if it still exists)."""
        primers = self._panel.primers()
        places = {p.name: self._panel.placement(p.name) for p in primers}
        names = [p.name for p in primers]
        for box, strand in ((self.forward, "+"), (self.reverse, "-")):
            old = box.currentData()
            box.blockSignals(True)
            box.clear()
            for n in names:
                box.addItem(n, n)
            if old in names:
                box.setCurrentIndex(names.index(old))
            else:
                ticked = [n for n in names if n in self._panel._ticked]
                pick = next((n for n in ticked if places.get(n) and places[n].strand == strand),
                            None) or next((n for n in names
                                           if places.get(n) and places[n].strand == strand), None)
                box.setCurrentIndex(names.index(pick) if pick else -1)
            box.blockSignals(False)
        self.schedule()

    def chosen(self):
        by_name = {p.name: p for p in self._panel.primers()}
        f = by_name.get(self.forward.currentData())
        r = by_name.get(self.reverse.currentData())
        return f, r

    # ------------------------------------------------------------------ live run

    def _on_run_where(self) -> None:
        """Switching between one file and many: new columns, new run."""
        many = self.run_where() == RUN_FILES
        self.table.setColumnCount(5 if many else 4)
        self.table.setHorizontalHeaderLabels(
            ["", "File", "Range", "Size", "mm F / R"] if many
            else ["", "Range", "Size", "mm F / R"])
        for c, w in enumerate((30, 110, 96, 70) if many else (30, 96, 70)):
            self.table.setColumnWidth(c, w)
        self._update_run_note()
        self.push_overlay()
        self.schedule()

    def _update_run_note(self, files=None) -> None:
        files = list(files if files is not None else self.templates())
        if self.run_where() == RUN_THIS:
            self.run_note.setText("The open sequence, live.")
            return
        if not files:
            self.run_note.setText("No file ticked in Files.")
            return
        names = ", ".join(p.name for p in files[:2])
        more = f" and {len(files) - 2} more" if len(files) > 2 else ""
        self.run_note.setText(f"{len(files)} file{'s' if len(files) != 1 else ''}: "
                              f"{names}{more}")

    def schedule(self, *_a) -> None:
        if self.active:
            self._timer.start()

    def run_where(self) -> str:
        button = self.run_group.checkedButton()
        return str(button.property("run")) if button is not None else RUN_THIS

    def templates(self) -> list[Path]:
        """The file (or files) the PCR runs on."""
        if self.run_where() == RUN_FILES:
            from repeatviewer.widgets.primers_panel import find_main_window

            mw = find_main_window()
            chosen = getattr(mw, "_pcr_template_candidates", None)
            if chosen is not None:
                try:
                    return [Path(p) for p in chosen()]
                except Exception:  # noqa: BLE001 - nothing ticked
                    return []
            return []
        path = self._panel.fasta_path()
        return [Path(path)] if path else []

    def _run(self) -> None:
        f, r = self.chosen()
        files = self.templates()
        self._generation += 1
        gen = self._generation
        self._update_run_note(files)
        if f is None or r is None or not files:
            self._result = None
            self._fill()
            if files:
                note = "Choose a forward and a reverse primer."
            elif self.run_where() == RUN_FILES:
                note = "Tick the FASTA files to test in the Files panel."
            else:
                note = "Open a sequence to run the PCR on."
            self._show_status(note, None)
            return
        s = self.settings()
        QSettings().setValue(_SETTINGS_KEY, __import__("json").dumps(s.to_dict()))
        primers = [f] if f.name == r.name else [f, r]
        relay = self._relay

        def work() -> None:
            try:
                res = pcr.run_pcr(files, primers, s)
            except Exception as exc:  # noqa: BLE001 - shown in the pane
                relay.done.emit(gen, None, str(exc))
                return
            relay.done.emit(gen, res, None)

        self.status_label.setText("Running…")
        threading.Thread(target=work, daemon=True).start()

    def _on_done(self, gen: int, result, error) -> None:
        if gen != self._generation:
            return   # a newer run is on its way
        if error:
            self._result = None
            self._fill()
            self._show_status(f"PCR failed: {error}", None)
            return
        self._result = result
        self._fill()

    # ------------------------------------------------------------------ results

    def products(self) -> list[pcr.PcrProduct]:
        if self._result is None:
            return []
        real = [p for p in self._result.products if p.kind != "dimer"]
        order = {"expected": 0, "ladder": 1, "off-target": 2}
        return sorted(real, key=lambda p: (order.get(p.kind, 3), -p.yield_, p.size))

    def _fill(self) -> None:
        prods = self.products()
        self.table.blockSignals(True)
        self.table.setRowCount(len(prods))
        for r, p in enumerate(prods):
            vals = [p.label or str(r + 1), f"{p.start:,}–{p.end:,}", f"{p.size:,} bp",
                    f"{p.mm_forward} / {p.mm_reverse}"]
            if self.run_where() == RUN_FILES:
                vals.insert(1, Path(getattr(p, "template", "") or "").name
                            or str(getattr(p, "record", "")))
            for c, v in enumerate(vals):
                self.table.setItem(r, c, QTableWidgetItem(v))
        self.table.blockSignals(False)
        n = len(prods)
        self.products_title.setText(f"<b>Products</b> &nbsp;<span style='color:gray'>{n}</span>")
        if prods:
            self.table.selectRow(0)
        try:
            from repeatviewer.backends.primers import gel_doc as gd

            if self._result is not None and prods:
                self.gel.set_doc(gd.gel_from_pcr(self._result))
        except Exception:  # noqa: BLE001 - the gel is a convenience
            pass
        self._update_status()
        self.push_overlay()
        self._sync_buttons()

    def _update_status(self) -> None:
        prods = self.products()
        if self._result is None:
            return
        if not prods:
            self._show_status("<b>No product.</b> The primers do not face each other on this "
                              "file within the mismatches allowed.", "amber")
            return
        best = prods[0]
        more = len(prods) - 1
        extra = f"; {more} more product{'s' if more != 1 else ''} possible (other sites)." if more else "."
        d = self._panel.design_pane.settings()
        if d.product_min <= best.size <= d.product_max:
            self._show_status(f"<b>Ready for PCR.</b> Expected product {best.size:,} bp{extra}",
                              "green")
        else:
            side = "longer" if best.size > d.product_max else "shorter"
            self._show_status(
                f"<b>Ready, with a warning.</b> Expected product {best.size:,} bp is {side} than "
                f"the product range ({d.product_min:,} to {d.product_max:,} bp) set in Design. "
                "The tube still makes it.", "amber", fit=best.size)

    def _show_status(self, text: str, tone, fit: int | None = None) -> None:
        colours = {"green": ("#E7F6EC", "#9BD3AE", "#14532D"),
                   "amber": ("#FFF4E5", "#E8A33D", "#5A3700")}
        bg, border, fg = colours.get(tone, ("transparent", "transparent", ""))
        self.status.setStyleSheet(
            f"#workbenchPcrStatus {{ background: {bg}; border: 1px solid {border};"
            f" border-radius: 6px; }}" + (f" #workbenchPcrStatus QLabel {{ color: {fg}; }}" if fg else ""))
        self.status_label.setText(text)
        self._fit_size = fit
        self.fit_button.setVisible(fit is not None)
        if fit is not None:
            self.fit_button.setText(f"Set Range to {self._range_for(fit)[1]:,} bp")
        self.status.setVisible(bool(text))

    @staticmethod
    def _range_for(size: int) -> tuple[int, int]:
        step = 100 if size < 2000 else 500
        hi = ((size * 104 // 100) // step + 1) * step
        return max(10, size // 2), hi

    def _on_fit(self) -> None:
        if not getattr(self, "_fit_size", None):
            return
        d = self._panel.design_pane
        lo, hi = self._range_for(self._fit_size)
        if self._fit_size < d.prod_min.value():
            d.prod_min.setValue(lo)
        if self._fit_size > d.prod_max.value():
            d.prod_max.setValue(hi)
        self._update_status()

    # ------------------------------------------------------------------ drawing

    def push_overlay(self) -> None:
        host = self._panel._host
        if host is None or not hasattr(host, "primers_host_set_workbench") or not self.active:
            return
        f, r = self.chosen()
        flen = len(th.binding_part(f.sequence)) if f is not None else 20
        rlen = len(th.binding_part(r.sequence)) if r is not None else 20
        seq_id = self._panel._seq_id
        cands, selected = [], -1
        rows = self.table.selectionModel().selectedRows() if self.table.selectionModel() else []
        chosen_row = rows[0].row() if rows else -1
        here = str(self._panel.fasta_path() or "")
        for i, p in enumerate(self.products()):
            if p.record and seq_id and p.record != seq_id:
                continue
            if here and p.template and str(p.template) != here:
                continue   # a product of another ticked file: not on this map
            cands.append((p.start, p.start + flen - 1, p.end - rlen + 1, p.end))
            if i == chosen_row:
                selected = len(cands) - 1
        try:
            host.primers_host_set_workbench(None, cands, selected)
        except Exception:  # noqa: BLE001
            pass
        if 0 <= selected < len(cands):
            show = getattr(host, "primers_host_show_range", None)
            if show is not None:
                try:
                    show(cands[selected][0], cands[selected][3])
                except Exception:  # noqa: BLE001
                    pass

    # ------------------------------------------------------------------ actions

    def _sync_buttons(self) -> None:
        has = self._result is not None and bool(self.products())
        for b in (self.bench_button, self.gel_button, self.save_button):
            b.setEnabled(has)

    def _main_window(self):
        from repeatviewer.widgets.primers_panel import find_main_window

        return find_main_window()

    def _on_bench(self) -> None:
        mw = self._main_window()
        if mw is not None and self._result is not None:
            mw._open_bench_card(self._result)

    def _on_gel(self) -> None:
        mw = self._main_window()
        if mw is not None and self._result is not None:
            mw._open_gel_for_pcr(self._result)

    def _on_save(self) -> None:
        if self._result is None:
            return
        try:
            out = pcr.write_outputs(self._result)
        except Exception as exc:  # noqa: BLE001
            self._show_status(f"Could not save the products: {exc}", "amber")
            return
        fasta = out.get("fasta")
        mw = self._main_window()
        if fasta and mw is not None:
            try:
                mw._on_preview_open_requested(Path(fasta))
            except Exception:  # noqa: BLE001
                pass


__all__ = ["PcrPane"]
