"""Details card for a read or feature clicked in a genome-browser track.

A small floating card (closes when clicking elsewhere or pressing Esc) with the
record's fields in two columns and buttons to copy the name, the sequence or all
the details.
"""

from __future__ import annotations

from PySide6.QtCore import QPoint, QRectF, Qt
from PySide6.QtGui import QColor, QGuiApplication, QPainter, QPen
from PySide6.QtWidgets import (
    QFrame,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QScrollArea,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends import alignment_style as style

MAX_ROWS_H = 360
MAX_W = 460


class TrackItemPopup(QFrame):
    def __init__(self, item, *, track_name: str = "", chrom: str = "", fetch=None,
                 parent: QWidget | None = None) -> None:
        """``fetch(start, end)`` reads the reference: features (and alignments without
        their own bases) then copy their sequence from it."""
        super().__init__(parent, Qt.WindowType.Popup | Qt.WindowType.FramelessWindowHint)
        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)
        # only the rounded card is painted: no square window behind its corners
        self.setAttribute(Qt.WidgetAttribute.WA_TranslucentBackground)
        self.setObjectName("trackItemPopup")
        self.item = item
        self.rows = style.details(item, chrom=chrom)
        self.chrom = chrom
        pal = self.palette()
        base = pal.color(pal.ColorRole.Base)
        dark = base.lightness() < 128
        border = "#3a3f47" if dark else "#d5d9e0"
        self._card = (QColor(base), QColor(border))
        muted = "#9aa3ad" if dark else "#6b7480"
        self.setStyleSheet(f"""
            QLabel#title {{ font-weight: 600; font-size: 13px; }}
            QLabel#sub {{ color: {muted}; font-size: 11px; }}
            QLabel.key {{ color: {muted}; }}
            QPushButton {{ padding: 3px 10px; }}
        """)
        outer = QVBoxLayout(self)
        outer.setContentsMargins(14, 12, 14, 10)
        outer.setSpacing(6)

        title_text = self.rows[0][1] if self.rows else "Item"
        title = QLabel(title_text, self)
        title.setObjectName("title")
        title.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        title.setMaximumWidth(MAX_W)
        title.setWordWrap(True)
        outer.addWidget(title)
        sub_parts = [r[1] for r in self.rows[1:2]]
        if track_name:
            sub_parts.append(track_name)
        sub_text = " · ".join(sub_parts)
        sub = QLabel(self)
        sub.setObjectName("sub")
        # long file names are shortened in the middle so the card stays narrow
        sub.setText(sub.fontMetrics().elidedText(sub_text, Qt.TextElideMode.ElideMiddle,
                                                 MAX_W))
        if sub.text() != sub_text:
            sub.setToolTip(sub_text)
        self.sub = sub
        outer.addWidget(sub)

        body = QWidget()
        grid = QGridLayout(body)
        grid.setContentsMargins(0, 4, 0, 0)
        grid.setHorizontalSpacing(14)
        grid.setVerticalSpacing(3)
        for r, (k, v) in enumerate(self.rows[2:]):
            key = QLabel(k, body)
            key.setProperty("class", "key")
            key.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignTop)
            val = QLabel(v, body)
            val.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
            val.setWordWrap(True)
            val.setMaximumWidth(360)
            grid.addWidget(key, r, 0)
            grid.addWidget(val, r, 1)
        grid.setColumnStretch(1, 1)
        scroll = QScrollArea(self)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        scroll.setWidget(body)
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setStyleSheet("QScrollArea { background: transparent; }")
        body.setStyleSheet("background: transparent;")
        body.adjustSize()
        scroll.setMinimumWidth(min(480, body.sizeHint().width() + 16))
        scroll.setFixedHeight(min(MAX_ROWS_H, body.sizeHint().height() + 6))
        outer.addWidget(scroll)

        buttons = QHBoxLayout()
        buttons.setSpacing(6)
        self.copy_name = QPushButton("Copy name", self)
        self.copy_name.clicked.connect(lambda: self._copy(getattr(item, "name", "") or ""))
        buttons.addWidget(self.copy_name)
        seq = getattr(item, "seq", "") or ""
        self.copy_seq = QPushButton("Copy sequence", self)
        self.fetch = fetch
        if seq:
            self.copy_seq.setToolTip("The read's bases")
            self.copy_seq.clicked.connect(lambda: self._copy(seq))
        elif fetch is not None:
            choices = style.feature_sequences(item, fetch)
            strand = " (reverse complement)" if getattr(item, "strand", "+") == "-" else ""
            if len(choices) == 1:
                self.copy_seq.setToolTip(f"Reference bases of the feature{strand}")
                self.copy_seq.clicked.connect(lambda: self._copy_feature("genomic"))
            else:
                from PySide6.QtWidgets import QMenu

                self.copy_seq.setText("Copy sequence ▾")
                menu = QMenu(self.copy_seq)
                for key, label, n in choices:
                    menu.addAction(f"{label} ({n:,} bp){strand}",
                                   lambda k=key: self._copy_feature(k))
                self.copy_seq.clicked.connect(
                    lambda: menu.exec(self.copy_seq.mapToGlobal(
                        self.copy_seq.rect().bottomLeft())))
                self.seq_menu = menu
        else:
            self.copy_seq.setEnabled(False)
            self.copy_seq.setToolTip("No bases: the file has no sequence for it and the "
                                     "reference is not loaded")
        buttons.addWidget(self.copy_seq)
        self.copy_all = QPushButton("Copy details", self)
        self.copy_all.clicked.connect(lambda: self._copy(self.details_text()))
        buttons.addWidget(self.copy_all)
        buttons.addStretch(1)
        outer.addLayout(buttons)

    def details_text(self) -> str:
        return "\n".join(f"{k}\t{v}" for k, v in self.rows)

    def feature_fasta(self, kind: str = "genomic") -> str:
        seq = style.feature_sequence(self.item, self.fetch, kind)
        name = getattr(self.item, "name", "") or "feature"
        suffix = {"genomic": "", "spliced": " spliced", "cds": " CDS"}.get(kind, "")
        lines = [seq[i:i + 60] for i in range(0, len(seq), 60)]
        it = self.item
        where = f"{self.chrom}:{it.start + 1}-{it.end}({getattr(it, 'strand', '+') or '+'})"
        return f">{name}{suffix} {where}\n" + "\n".join(lines) + "\n"

    def _copy_feature(self, kind: str) -> None:
        self._copy(self.feature_fasta(kind))

    def _copy(self, text: str) -> None:
        QGuiApplication.clipboard().setText(text)
        self.close()

    def show_at(self, global_pos: QPoint) -> None:
        self.adjustSize()
        screen = QGuiApplication.screenAt(global_pos) or QGuiApplication.primaryScreen()
        area = screen.availableGeometry()
        x = min(global_pos.x() + 12, area.right() - self.width() - 4)
        y = global_pos.y() + 12
        if y + self.height() > area.bottom():
            y = max(area.top(), global_pos.y() - self.height() - 12)
        self.move(max(area.left(), x), y)
        self.show()

    def paintEvent(self, event) -> None:
        # the card itself: rounded, the window around its corners stays transparent
        fill, edge = self._card
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        p.setPen(QPen(edge, 1))
        p.setBrush(fill)
        p.drawRoundedRect(QRectF(self.rect()).adjusted(0.5, 0.5, -0.5, -0.5), 8, 8)
        p.end()

    def keyPressEvent(self, event) -> None:
        if event.key() == Qt.Key.Key_Escape:
            self.close()
            return
        super().keyPressEvent(event)
