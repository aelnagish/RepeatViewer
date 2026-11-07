"""Binding Sites on Chromosomes: where primers bind across a genome.

Each chromosome is a strip drawn to scale; forward sites are ticks above it,
reverse sites below, coloured by mismatches (darker = fewer). Products are
orange boxes. The side panel counts sites per chromosome and gives the number
of templates the pair can copy (the expected value for a qPCR copy-number assay).
"""

from __future__ import annotations

import html
import threading
from pathlib import Path

from PySide6.QtCore import QObject, QRectF, Qt, Signal
from PySide6.QtGui import QColor, QPainter
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
    QScrollArea,
    QSizePolicy,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends.primers import pcr
from repeatviewer.backends.primers.oligos import Primer
from repeatviewer.widgets.help_button import make_help_button
from repeatviewer.widgets.path_browse import PathBrowseRow

MM_COLORS = ("#0f3a7a", "#4f86d6", "#a9c6ef", "#d6e3f5", "#e9eef7")
PRODUCT_COLOR = "#f2a65a"
MAX_STRIPS = 60

GENOME_HELP = (
    "<b>Genome</b>: the FASTA to search (the project genome is filled in; click a file in "
    "the Files panel to change it).<br>"
    "Each chromosome is drawn to scale. <b>Forward</b> sites are ticks above the bar, "
    "<b>reverse</b> sites below; darker blue means fewer mismatches. <b>Orange</b> boxes are "
    "products the primers can make there.<br>"
    "<b>Templates for qPCR</b> counts the products of the pair per haploid genome, the value "
    "to expect in a copy-number assay.<br>"
    "Only the 60 longest sequences are drawn; the table lists every sequence with a site."
)


class _Relay(QObject):
    done = Signal(object)
    failed = Signal(str)
    progress = Signal(str)


class GenomeStrips(QWidget):
    recordClicked = Signal(str)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._records: list[pcr.RecordSites] = []
        self._shown: set[str] = set()
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.setMinimumWidth(500)

    def set_records(self, records: list[pcr.RecordSites], shown_primers: set[str]) -> None:
        recs = sorted(records, key=lambda r: -r.length)[:MAX_STRIPS]
        # keep the FASTA order for the strips that are drawn
        order = {r.record: i for i, r in enumerate(records)}
        self._records = sorted(recs, key=lambda r: order[r.record])
        self._shown = set(shown_primers)
        self.setFixedHeight(20 + 46 * max(1, len(self._records)))
        self.update()

    def paintEvent(self, event) -> None:  # noqa: N802
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        text = self.palette().color(self.foregroundRole())
        if not self._records:
            p.setPen(text)
            p.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, "Choose a genome and press Find Sites")
            p.end()
            return
        longest = max(r.length for r in self._records) or 1
        name_w = 90.0
        left = name_w + 12
        avail = max(100.0, self.width() - left - 220)
        fm = p.fontMetrics()
        y = 12.0
        for rs in self._records:
            w = max(4.0, rs.length / longest * avail)
            p.setPen(text)
            p.drawText(QRectF(0, y + 4, name_w, 20), Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter,
                       fm.elidedText(rs.record, Qt.TextElideMode.ElideMiddle, int(name_w)))
            bar = QRectF(left, y + 9, w, 10)
            p.setPen(Qt.PenStyle.NoPen)
            p.setBrush(QColor("#e1e1e6"))
            p.drawRoundedRect(bar, 5, 5)
            # products
            for a, b, _size in rs.products[:200]:
                x1 = left + a / rs.length * w
                x2 = left + b / rs.length * w
                p.setBrush(Qt.BrushStyle.NoBrush)
                p.setPen(QColor(PRODUCT_COLOR))
                p.drawRect(QRectF(x1 - 2, y + 4, max(4.0, x2 - x1) + 4, 20))
            p.setPen(Qt.PenStyle.NoPen)
            for primer, strand, pos, mm in rs.sites:
                if self._shown and primer not in self._shown:
                    continue
                x = left + pos / max(1, rs.length) * w
                p.setBrush(QColor(MM_COLORS[min(mm, len(MM_COLORS) - 1)]))
                if strand == "+":
                    p.drawRect(QRectF(x - 1, y - 2, 2, 11))
                else:
                    p.drawRect(QRectF(x - 1, y + 19, 2, 11))
            counts = rs.counts()
            if counts:
                summary = ", ".join(f"{n} × {mm} mm" for mm, n in sorted(counts.items()))
                p.setPen(text)
                p.drawText(QRectF(left + w + 10, y + 4, 210, 20),
                           Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter, summary)
            y += 46
        p.end()


class PrimerGenomeDialog(QDialog):
    fileWritten = Signal(object)

    def __init__(self, primers: list[Primer], *, genome: Path | None = None,
                 parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Binding Sites on Chromosomes")
        self.setMinimumSize(1200, 760)
        self._primers = list(primers)
        self._records: list[pcr.RecordSites] = []
        self._relay = _Relay(self)
        self._relay.done.connect(self._on_done)
        self._relay.failed.connect(self._on_failed)
        self._relay.progress.connect(lambda t: self.status.setText(t))

        root = QVBoxLayout(self)
        root.setContentsMargins(16, 12, 16, 12)
        root.setSpacing(10)
        top = QHBoxLayout()
        top.addWidget(QLabel("Genome", self))
        self.genome_row = PathBrowseRow(self, initial=genome, placeholder="Genome FASTA",
                                        browse_title="Choose the genome FASTA",
                                        file_filter="FASTA (*.fa *.fasta *.fna *.fa.gz);;All files (*)")
        self.genome_row.setMaximumWidth(460)
        top.addWidget(self.genome_row)
        top.addSpacing(10)
        top.addWidget(QLabel("Mismatches", self))
        self.mm_combo = QComboBox(self)
        for n in range(0, 4):
            self.mm_combo.addItem("none" if n == 0 else f"up to {n}", n)
        self.mm_combo.setCurrentIndex(2)
        top.addWidget(self.mm_combo)
        self.find_button = QPushButton("Find Sites", self)
        self.find_button.setDefault(True)
        top.addWidget(self.find_button)
        top.addStretch(1)
        legend = QLabel(
            " ".join(f'<span style="color:{c}">■</span> {i}' for i, c in enumerate(MM_COLORS[:3]))
            + " mismatches &nbsp; "
            f'<span style="color:{PRODUCT_COLOR}">□</span> product', self)
        top.addWidget(legend)
        top.addWidget(make_help_button(self, "Binding Sites on Chromosomes", GENOME_HELP))
        root.addLayout(top)
        prim = QLabel("Primers: " + (", ".join(html.escape(p.name) for p in self._primers)
                                     or "none (open from the Primers tab of a sequence)"), self)
        prim.setStyleSheet("color: palette(placeholder-text);")
        root.addWidget(prim)

        body = QHBoxLayout()
        self.strips = GenomeStrips(self)
        scroll = QScrollArea(self)
        scroll.setWidgetResizable(True)
        scroll.setWidget(self.strips)
        body.addWidget(scroll, 1)
        side = QVBoxLayout()
        side.addWidget(QLabel("<b>Sites per sequence</b>", self))
        self.table = QTableWidget(0, 5, self)
        self.table.setHorizontalHeaderLabels(["Sequence", "0 mm", "1 mm", "2+ mm", "Products"])
        self.table.verticalHeader().hide()
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setShowGrid(False)
        h = self.table.horizontalHeader()
        h.setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
        h.setDefaultAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        h.setStretchLastSection(True)
        for c, w in enumerate((110, 46, 46, 54)):
            self.table.setColumnWidth(c, w)
        self.table.setFixedWidth(340)
        side.addWidget(self.table, 1)
        self.summary = QLabel(self)
        self.summary.setWordWrap(True)
        self.summary.setTextFormat(Qt.TextFormat.RichText)
        self.summary.setFixedWidth(340)
        self.summary.setStyleSheet("QLabel { border: 1px solid palette(mid); border-radius: 8px;"
                                   " padding: 8px; }")
        side.addWidget(self.summary)
        self.gff_button = QPushButton("Save Sites as GFF3…", self)
        self.png_button = QPushButton("Save Figure (PNG)…", self)
        side.addWidget(self.gff_button)
        side.addWidget(self.png_button)
        body.addLayout(side)
        root.addLayout(body, 1)
        self.status = QLabel(self)
        self.status.setStyleSheet("color: palette(placeholder-text);")
        root.addWidget(self.status)

        self.find_button.clicked.connect(self._find)
        self.gff_button.clicked.connect(self._save_gff)
        self.png_button.clicked.connect(self._save_png)
        self.gff_button.setEnabled(False)
        self.png_button.setEnabled(False)
        self.summary.hide()

    def apply_path_from_panel(self, path: Path) -> bool:
        self.genome_row.set_path(path)
        return True

    def _find(self) -> None:
        genome = self.genome_row.path()
        if genome is None or not Path(genome).is_file():
            QMessageBox.information(self, "Binding Sites", "Choose the genome FASTA first.")
            return
        if not self._primers:
            QMessageBox.information(self, "Binding Sites",
                                    "No primers. Open this window from the Primers tab of a sequence.")
            return
        self.find_button.setEnabled(False)
        self.status.setText(f"Searching {Path(genome).name}…")
        mm = int(self.mm_combo.currentData())
        relay = self._relay
        primers = list(self._primers)

        def work() -> None:
            try:
                recs = pcr.genome_sites(Path(genome), primers, max_mismatches=mm,
                                        progress=relay.progress.emit)
                relay.done.emit(recs)
            except Exception as exc:  # noqa: BLE001
                relay.failed.emit(str(exc))

        threading.Thread(target=work, daemon=True).start()

    def _on_failed(self, error: str) -> None:
        self.find_button.setEnabled(True)
        self.status.setText("")
        QMessageBox.warning(self, "Binding Sites", error)

    def _on_done(self, records) -> None:
        self.find_button.setEnabled(True)
        self._records = list(records)
        self.strips.set_records(self._records, {p.name for p in self._primers})
        with_sites = [r for r in self._records if r.sites]
        self.table.setRowCount(len(with_sites))
        for i, r in enumerate(with_sites):
            c = r.counts()
            vals = [r.record, c.get(0, 0), c.get(1, 0), sum(n for mm, n in c.items() if mm >= 2),
                    len(r.products)]
            for j, v in enumerate(vals):
                self.table.setItem(i, j, QTableWidgetItem(str(v)))
        total = sum(len(r.sites) for r in self._records)
        prods = sum(len(r.products) for r in self._records)
        pair_products = 0
        if len(self._primers) >= 2:
            a, b = self._primers[0].name, self._primers[1].name
            pair_products = prods  # all combinations, as in a tube
            note = f"{a} + {b}"
        else:
            note = self._primers[0].name if self._primers else ""
        self.summary.setText(
            f"<b>{total:,}</b> sites on <b>{len(with_sites)}</b> of {len(self._records):,} sequences."
            f"<br><b>{prods:,}</b> products up to 5 kb.<br>"
            f"Templates for qPCR ({html.escape(note)}): <b>{pair_products:,}</b> per haploid genome.")
        self.summary.show()
        self.gff_button.setEnabled(bool(total))
        self.png_button.setEnabled(True)
        drawn = min(len(self._records), MAX_STRIPS)
        self.status.setText(f"{len(self._records):,} sequences searched; {drawn} drawn.")

    def _save_gff(self) -> None:
        genome = self.genome_row.path()
        start = str(Path(genome).with_name(Path(genome).name + "_primer_sites.gff3")) if genome else ""
        fn, _ = QFileDialog.getSaveFileName(self, "Save Sites as GFF3", start, "GFF3 (*.gff3)")
        if not fn:
            return
        Path(fn).write_text(pcr.sites_gff3(self._records), encoding="utf-8")
        self.fileWritten.emit(Path(fn))

    def _save_png(self) -> None:
        from repeatviewer.widgets.gel_view import save_view_image

        genome = self.genome_row.path()
        start = str(Path(genome).with_name(Path(genome).name + "_primer_sites.png")) if genome else ""
        fn, _ = QFileDialog.getSaveFileName(self, "Save Figure", start, "PNG (*.png)")
        if not fn:
            return
        save_view_image(self.strips, Path(fn), "png")
        self.fileWritten.emit(Path(fn))
