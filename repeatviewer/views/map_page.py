"""Map as a tab in the middle: where the primers bind, and what they amplify.

Left: one row per sequence with the binding sites drawn on it, forward blue,
reverse orange, and a line between the two ends of a product. Right: the
primers, the sequences to test, how many mismatches are allowed, and what
counts as a product. The report is written only when you ask for it.

Developer: Abdullah Yusuf
"""

from __future__ import annotations

import html
import json
import threading
from pathlib import Path

from PySide6.QtCore import QObject, QRectF, QSettings, Qt, Signal
from PySide6.QtGui import QColor, QPainter, QPolygonF
from PySide6.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QComboBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends.primers import jobs as pj
from repeatviewer.widgets.help_button import make_help_button
from repeatviewer.widgets.primer_design_pane import (
    BOX2_W,
    FIELD_W,
    LABEL_W,
    _row,
    _rule,
    _section,
    _spin,
)

_SETTINGS_KEY = "primers/map_page"
PANEL_W = 330
FORWARD, REVERSE, PROBE = "#1F57C3", "#C2410C", "#7A3E9D"
MAX_DRAWN = 40

PRIMERS_HELP = (
    "Which primers to test. <b>Primer Database</b> uses the saved ones; pick a forward and "
    "a reverse, or leave one empty to look for single binding sites.<br>"
    "<b>Pairs only</b> keeps the primers that give a product with their partner."
)
SEQUENCES_HELP = (
    "The FASTA files ticked in the Files panel. <b>Consensus</b> tests the primers on one "
    "consensus of an alignment instead of each row."
)
MISMATCH_HELP = (
    "<b>Max mismatches</b>: how different a site may be and still count.<br>"
    "<b>3′ bases exact</b>: the last bases of the primer that must match, because a "
    "mismatch there stops the polymerase. For a probe this is 0."
)
PRODUCT_HELP = (
    "Two sites facing each other within this size range count as a product. "
    "<b>Count sites in the genome</b> also searches the project genome, which takes longer."
)
PICTURE_HELP = (
    "One row per sequence. Forward sites are blue, reverse orange, probes purple; the line "
    "between two sites is a product, with its size. Only the first "
    f"{MAX_DRAWN} sequences are drawn."
)


class _Relay(QObject):
    done = Signal(object, object)       # payload | None, error | None


class BindingPicture(QWidget):
    """Where every primer binds on every sequence, one row per sequence."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._rows: list[dict] = []
        self.setMinimumHeight(80)

    def set_payload(self, payload: dict | None) -> None:
        self._rows = []
        if payload:
            lengths = {t["id"]: int(t.get("length") or 0) for t in payload.get("templates", [])}
            by_record: dict[str, dict] = {}
            for site in payload.get("sites", []):
                row = by_record.setdefault(site["record"], {"sites": [], "products": []})
                row["sites"].append(site)
            for product in payload.get("products", []):
                by_record.setdefault(product["record"], {"sites": [], "products": []})
                by_record[product["record"]]["products"].append(product)
            for rid, row in list(by_record.items())[:MAX_DRAWN]:
                row["id"] = rid
                row["length"] = lengths.get(rid) or max(
                    [s["end"] for s in row["sites"]] + [p["end"] for p in row["products"]] + [1])
                self._rows.append(row)
        self.setMinimumHeight(24 + 46 * max(1, len(self._rows)))
        self.update()

    def paintEvent(self, event) -> None:  # noqa: N802
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        if not self._rows:
            painter.setPen(QColor("#6b7280"))
            painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter,
                             "Press Map to see where the primers bind.")
            painter.end()
            return
        font = painter.font()
        font.setPointSizeF(max(7.5, font.pointSizeF() * 0.82))
        painter.setFont(font)
        name_w, size_w = 150.0, 70.0
        left = name_w + 8
        right = max(left + 60, self.width() - size_w - 8)
        span = right - left
        y = 18.0
        for row in self._rows:
            length = max(1, int(row["length"]))

            def x_of(pos: int, length=length) -> float:
                return left + span * min(max(0, pos - 1), length) / length

            painter.setPen(QColor("#3A3A3C"))
            painter.drawText(QRectF(0, y - 8, name_w, 16),
                             Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter,
                             str(row["id"])[:28])
            painter.setPen(QColor("#5A5D63"))
            painter.drawText(QRectF(right + 6, y - 8, size_w, 16),
                             Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter,
                             f"{length:,} bp")
            painter.fillRect(QRectF(left, y - 4, span, 9), QColor("#D9D9DE"))
            for site in row["sites"]:
                x = x_of(int(site["start"]))
                strand = site.get("strand", "+")
                colour = QColor(FORWARD if strand == "+" else REVERSE)
                if site.get("mismatches"):
                    colour.setAlpha(150)
                tip = 12.0 if strand == "+" else -12.0
                poly = QPolygonF([
                    QRectF(x, y - 4, 1, 1).topLeft(),
                    QRectF(x + tip, y - 4, 1, 1).topLeft(),
                    QRectF(x + tip + (5 if strand == "+" else -5), y + 0.5, 1, 1).topLeft(),
                    QRectF(x + tip, y + 5, 1, 1).topLeft(),
                    QRectF(x, y + 5, 1, 1).topLeft(),
                ])
                painter.setBrush(colour)
                painter.setPen(Qt.PenStyle.NoPen)
                painter.drawPolygon(poly)
            painter.setBrush(Qt.BrushStyle.NoBrush)
            last_label_right = -1e9
            for product in row["products"][:12]:
                a, b = x_of(int(product["start"])), x_of(int(product["end"]))
                painter.setPen(QColor("#1C1C1E"))
                painter.drawLine(QRectF(a, y + 14, 1, 1).topLeft(),
                                 QRectF(b, y + 14, 1, 1).topLeft())
                text = f"{int(product['size']):,} bp"
                width = painter.fontMetrics().horizontalAdvance(text) + 8
                middle = (a + b) / 2
                if middle - width / 2 > last_label_right:     # no label on top of another
                    painter.drawText(QRectF(middle - width / 2, y + 16, width, 14),
                                     Qt.AlignmentFlag.AlignCenter, text)
                    last_label_right = middle + width / 2
            if len(row["products"]) > 1:
                painter.setPen(QColor("#5A5D63"))
                painter.drawText(QRectF(right + 6, y + 14, size_w, 14),
                                 Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter,
                                 f"{len(row['products'])} products")
            y += 46
        painter.end()


class MapPage(QWidget):
    """Test primers on the selected sequences and show where they bind."""

    def __init__(self, files=(), *, options=None, database_path: str = "",
                 parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.path = None
        self._files = [Path(p) for p in files]
        self._database_path = database_path
        self._options = options
        self._payload: dict | None = None
        self._busy = False
        self._relay = _Relay(self)
        self._relay.done.connect(self._on_done)
        saved = self._load()
        base = dict(getattr(options, "map", None) or {})
        base.update(saved)

        outer = QHBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)
        middle = QWidget(self)
        lay = QVBoxLayout(middle)
        lay.setContentsMargins(14, 12, 14, 12)
        lay.setSpacing(10)
        self.verdict = QLabel(middle)
        self.verdict.setWordWrap(True)
        self.verdict.setObjectName("mapVerdict")
        lay.addWidget(self.verdict)

        self.table = QTableWidget(0, 6, middle)
        self.table.setHorizontalHeaderLabels(
            ["Sequence", "Forward", "Reverse", "Sites", "Product", "Mismatches"])
        self.table.verticalHeader().hide()
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setShowGrid(False)
        self.table.setAlternatingRowColors(True)
        head = self.table.horizontalHeader()
        head.setDefaultAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        head.setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        for c, w in ((1, 90), (2, 90), (3, 70), (4, 150), (5, 110)):
            self.table.setColumnWidth(c, w)
        self.table.setMinimumHeight(120)
        lay.addWidget(self.table)

        pic_head = QHBoxLayout()
        pic_head.addWidget(QLabel("<b>Where they bind</b>", middle), 1)
        pic_head.addWidget(make_help_button(middle, "Where they bind", PICTURE_HELP))
        lay.addLayout(pic_head)
        self.picture = BindingPicture(middle)
        lay.addWidget(self.picture, 1)
        outer.addWidget(middle, 1)

        panel = QWidget(self)
        panel.setFixedWidth(PANEL_W)
        pl = QVBoxLayout(panel)
        pl.setContentsMargins(10, 12, 10, 12)
        pl.setSpacing(10)
        title = QHBoxLayout()
        title.addWidget(QLabel("<b>Map</b>", panel), 1)
        title.addWidget(make_help_button(panel, "Map", PRIMERS_HELP))
        pl.addLayout(title)
        pl.addWidget(_rule(panel))

        sec, sl, _ = _section("Primers", panel, "Primers", PRIMERS_HELP)
        self.source = QComboBox(sec)
        self.source.addItem("Primer Database", pj.MAP_DATABASE)
        self.source.addItem("Typed in", pj.MAP_SPECIFIC)
        self.source.setFixedWidth(FIELD_W)
        i = self.source.findData(base.get("mode", pj.MAP_DATABASE))
        self.source.setCurrentIndex(max(0, i))
        sl.addLayout(_row("From", self.source, parent=sec))
        self.forward = QComboBox(sec)
        self.reverse = QComboBox(sec)
        for box in (self.forward, self.reverse):
            box.setEditable(True)
            box.setFixedWidth(FIELD_W)
        sl.addLayout(_row("Forward", self.forward, parent=sec))
        sl.addLayout(_row("Reverse", self.reverse, parent=sec))
        self.pairs_only = QCheckBox("Pairs only", sec)
        self.pairs_only.setChecked(bool(base.get("pairs_only", False)))
        sl.addLayout(_row("", self.pairs_only, parent=sec))
        pl.addWidget(sec)
        pl.addWidget(_rule(panel))

        sec, sl, _ = _section("Sequences", panel, "Sequences", SEQUENCES_HELP)
        self.test_on = QComboBox(sec)
        self.test_on.addItem("The selected files", pj.TEST_ON_SEQUENCES)
        self.test_on.addItem("Consensus of an alignment", pj.TEST_ON_CONSENSUS)
        self.test_on.setFixedWidth(FIELD_W)
        i = self.test_on.findData(base.get("test_on", pj.TEST_ON_SEQUENCES))
        self.test_on.setCurrentIndex(max(0, i))
        sl.addLayout(_row("Test on", self.test_on, parent=sec))
        self.files_note = QLabel(sec)
        self.files_note.setWordWrap(True)
        self.files_note.setIndent(LABEL_W + 2)
        self.files_note.setStyleSheet("color: palette(placeholder-text);")
        sl.addWidget(self.files_note)
        pl.addWidget(sec)
        pl.addWidget(_rule(panel))

        sec, sl, _ = _section("Mismatches", panel, "Mismatches", MISMATCH_HELP)
        self.mismatches = _spin(0, 10, int(base.get("max_mismatches", 2)), parent=sec)
        sl.addLayout(_row("Max", self.mismatches, parent=sec, unit="nt"))
        self.exact = _spin(0, 20, int(base.get("three_prime_exact", 3)), parent=sec)
        sl.addLayout(_row("3′ exact", self.exact, parent=sec, unit="nt"))
        pl.addWidget(sec)
        pl.addWidget(_rule(panel))

        sec, sl, _ = _section("Products", panel, "Products", PRODUCT_HELP)
        self.size_min = _spin(1, 100_000, int(base.get("product_min", 1)), parent=sec,
                              width=BOX2_W)
        self.size_max = _spin(1, 1_000_000, int(base.get("product_max", 5000)), parent=sec,
                              width=BOX2_W)
        for box in (self.size_min, self.size_max):
            box.setGroupSeparatorShown(True)
        sl.addLayout(_row("Size", self.size_min, self.size_max, parent=sec, joins=("to",)))
        self.genome_check = QCheckBox("Count sites in the genome", sec)
        self.genome_check.setChecked(bool(base.get("check_genome", False)))
        sl.addLayout(_row("", self.genome_check, parent=sec))
        pl.addWidget(sec)
        pl.addWidget(_rule(panel))

        sec, sl, _ = _section("Options", panel, "Options", PRIMERS_HELP)
        self.annotate_check = QCheckBox("Add the sites to the sequence", sec)
        self.annotate_check.setChecked(bool(base.get("annotate", False)))
        self.report_check = QCheckBox("Write a report beside the files", sec)
        self.report_check.setChecked(bool(saved.get("report", False)))
        for box in (self.annotate_check, self.report_check):
            sl.addLayout(_row("", box, parent=sec))
        self.defaults_button = QPushButton("Restore Defaults", sec)
        self.save_button = QPushButton("Save", sec)
        self.map_button = QPushButton("Map", sec)
        self.map_button.setDefault(True)
        br = QHBoxLayout()
        br.addWidget(self.defaults_button)
        br.addStretch(1)
        br.addWidget(self.save_button)
        br.addWidget(self.map_button)
        sl.addLayout(br)
        self.saved_note = QLabel("", sec)
        self.saved_note.setWordWrap(True)
        self.saved_note.setStyleSheet("color: palette(placeholder-text);")
        sl.addWidget(self.saved_note)
        pl.addWidget(sec)
        pl.addStretch(1)
        outer.addWidget(panel)

        self.map_button.clicked.connect(self.run_map)
        self.save_button.clicked.connect(self._on_save)
        self.defaults_button.clicked.connect(self._on_defaults)
        self.source.currentIndexChanged.connect(self._fill_primers)
        self._fill_primers()
        self._refresh_files()
        self._show_idle()

    # ------------------------------------------------------------------ settings

    def _load(self) -> dict:
        raw = QSettings().value(_SETTINGS_KEY, "")
        try:
            d = json.loads(raw) if raw else {}
        except (TypeError, ValueError):
            d = {}
        return d if isinstance(d, dict) else {}

    def settings(self) -> pj.MapSettings:
        base = dict(getattr(self._options, "map", None) or {})
        base.update(
            mode=str(self.source.currentData()),
            use_forward=bool(self.forward.currentText().strip()),
            forward=self.forward.currentText().strip(),
            use_reverse=bool(self.reverse.currentText().strip()),
            reverse=self.reverse.currentText().strip(),
            pairs_only=self.pairs_only.isChecked(),
            test_on=str(self.test_on.currentData()),
            use_max_mismatches=True,
            max_mismatches=int(self.mismatches.value()),
            three_prime_exact=int(self.exact.value()),
            use_product_range=True,
            product_min=int(self.size_min.value()),
            product_max=int(self.size_max.value()),
            check_genome=self.genome_check.isChecked(),
            annotate=self.annotate_check.isChecked(),
        )
        if str(self.source.currentData()) == pj.MAP_DATABASE:
            names = [n for n in (self.forward.currentText().strip(),
                                 self.reverse.currentText().strip()) if n]
            base["primer_names"] = names
        return pj.MapSettings.from_dict(base)

    def _on_save(self) -> None:
        s = self.settings().to_dict()
        s["report"] = self.report_check.isChecked()
        QSettings().setValue(_SETTINGS_KEY, json.dumps(s))
        self.saved_note.setText("Saved.")

    def _on_defaults(self) -> None:
        d = pj.MapSettings()
        self.mismatches.setValue(d.max_mismatches)
        self.exact.setValue(d.three_prime_exact)
        self.size_min.setValue(d.product_min)
        self.size_max.setValue(d.product_max)
        self.pairs_only.setChecked(d.pairs_only)
        self.genome_check.setChecked(d.check_genome)
        self.annotate_check.setChecked(d.annotate)
        self.report_check.setChecked(False)

    # ------------------------------------------------------------------ inputs

    def tab_name(self) -> str:
        names = [n for n in (self.forward.currentText().strip(),
                             self.reverse.currentText().strip()) if n]
        return "Map · " + (" / ".join(names) if names else "primers")

    def set_files(self, files) -> None:
        self._files = [Path(p) for p in files]
        self._refresh_files()

    def _refresh_files(self) -> None:
        if not self._files:
            self.files_note.setText("No file ticked in Files.")
            return
        names = ", ".join(p.name for p in self._files[:3])
        more = f" and {len(self._files) - 3} more" if len(self._files) > 3 else ""
        self.files_note.setText(f"{len(self._files)} file"
                                f"{'s' if len(self._files) != 1 else ''}: {names}{more}")

    def _fill_primers(self) -> None:
        """The Database names, or free text when the primers are typed in."""
        database = str(self.source.currentData()) == pj.MAP_DATABASE
        names: list[str] = []
        if database:
            try:
                from repeatviewer.backends.primers.inventory import open_database

                db = open_database(Path(self._database_path) if self._database_path else None)
                names = [p.name for p in db.primers]
            except Exception:  # noqa: BLE001 - no database yet
                names = []
        for box in (self.forward, self.reverse):
            text = box.currentText()
            box.blockSignals(True)
            box.clear()
            box.addItem("")
            for name in names:
                box.addItem(name)
            box.setEditable(not database)
            box.setCurrentText(text)
            box.blockSignals(False)

    # ------------------------------------------------------------------ run

    def run_map(self) -> None:
        if self._busy:
            return
        if not self._files:
            self._set_verdict("amber", "Tick the FASTA files to test in the Files panel.")
            return
        s = self.settings()
        if not (s.forward or s.reverse):
            self._set_verdict("amber", "Choose at least one primer.")
            return
        self._busy = True
        self.map_button.setEnabled(False)
        self._set_verdict("", "Looking for binding sites…")
        files = list(self._files)
        database = Path(self._database_path) if self._database_path else None
        relay = self._relay

        def work() -> None:
            try:
                payload, _sites, _products = pj.map_payload(files, s, database_path=database)
            except Exception as exc:  # noqa: BLE001 - shown in the page
                relay.done.emit(None, str(exc))
                return
            relay.done.emit(payload, None)

        threading.Thread(target=work, daemon=True).start()

    def _on_done(self, payload, error) -> None:
        self._busy = False
        self.map_button.setEnabled(True)
        if error:
            self._payload = None
            self._set_verdict("amber", f"Map failed: {error}")
            self.table.setRowCount(0)
            self.picture.set_payload(None)
            return
        self._payload = payload
        self._fill(payload)
        if self.report_check.isChecked():
            self._write_report()

    # ------------------------------------------------------------------ view

    def _set_verdict(self, tone: str, text: str) -> None:
        if tone == "green":
            bg, border, fg = "#E7F6EC", "#9BD3AE", "#14532D"
        elif tone == "amber":
            bg, border, fg = "#FFF4E5", "#E8A33D", "#5A3700"
        else:
            self.verdict.setStyleSheet("")
            self.verdict.setText(text)
            return
        self.verdict.setStyleSheet(
            f"#mapVerdict {{ background: {bg}; border: 1px solid {border};"
            f" border-radius: 6px; padding: 8px; color: {fg}; }}")
        self.verdict.setText(html.escape(text))

    def _show_idle(self) -> None:
        self._set_verdict("", "Choose the primers and press Map.")
        self.picture.set_payload(None)

    def _fill(self, payload: dict) -> None:
        counts = payload.get("counts", {})
        rows: dict[str, dict] = {}
        for t in payload.get("templates", []):
            rows[t["id"]] = {"forward": [], "reverse": [], "sites": 0, "products": [],
                             "mismatches": []}
        for site in payload.get("sites", []):
            row = rows.setdefault(site["record"], {"forward": [], "reverse": [], "sites": 0,
                                                   "products": [], "mismatches": []})
            row["sites"] += 1
            (row["forward"] if site.get("strand", "+") == "+" else row["reverse"]).append(
                int(site["start"]))
            row["mismatches"].append(int(site.get("mismatches", 0)))
        for product in payload.get("products", []):
            rows.setdefault(product["record"], {"forward": [], "reverse": [], "sites": 0,
                                                "products": [], "mismatches": []})
            rows[product["record"]]["products"].append(int(product["size"]))
        self.table.setRowCount(len(rows))
        amber = QColor("#B26B00")
        for r, (rid, row) in enumerate(rows.items()):
            sizes = ", ".join(f"{s:,} bp" for s in sorted(set(row["products"]))[:3]) or "none"
            mism = ("/".join(str(m) for m in sorted(set(row["mismatches"]))) or "–")
            cells = [rid,
                     f"{row['forward'][0]:,}" if row["forward"] else "no site",
                     f"{row['reverse'][0]:,}" if row["reverse"] else "no site",
                     str(row["sites"]), sizes, mism]
            for c, text in enumerate(cells):
                item = QTableWidgetItem(text)
                if text in ("no site", "none"):
                    item.setForeground(amber)
                self.table.setItem(r, c, item)
        self.picture.set_payload(payload)
        with_product = int(counts.get("sequences_with_product", 0))
        total = int(counts.get("sequences", len(rows)))
        if with_product and with_product == total:
            self._set_verdict("green", f"A product in all {total} sequence"
                                       f"{'s' if total != 1 else ''}: "
                                       f"{counts.get('products', 0)} in total.")
        elif with_product:
            self._set_verdict("amber", f"A product in {with_product} of {total} sequences. "
                                       "The others are in the table.")
        else:
            self._set_verdict("amber", f"No product. {counts.get('sites', 0)} binding site"
                                       f"{'s' if counts.get('sites', 0) != 1 else ''} found.")

    def _write_report(self) -> None:
        if not self._files or self._payload is None:
            return
        try:
            layout = pj.plan_outputs(self._files[0], pj.MAP_TAG, None)
            files = pj.write_result(layout, self._payload)
        except Exception as exc:  # noqa: BLE001 - shown in the page
            self.saved_note.setText(f"Could not write the report: {exc}")
            return
        self.report_path = layout.main
        self.saved_note.setText(f"Report: {layout.main.name}")
        main = self.window()
        show = getattr(main, "_show_new_output", None)
        if show is not None and files:
            show(layout.main)
