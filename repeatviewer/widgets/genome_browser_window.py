"""Tools > Genome Browser: a reference FASTA with annotation and alignment tracks in a
window of its own (the contig view in browser mode, without the Preview panel)."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QSettings, Signal
from PySide6.QtGui import QGuiApplication, QKeySequence
from PySide6.QtWidgets import QMainWindow, QMenu, QToolButton, QWidget

from repeatviewer.widgets.contig_view import ContigViewPanel

_GEOMETRY_KEY = "genome_browser/geometry"
RECENT_KEY = "genome_browser/recent"      # session files, last used first
MAX_RECENT = 8


def _settings() -> QSettings:
    return QSettings()


def recent_sessions() -> list[Path]:
    """Genome browser sessions opened before, last used first (existing ones only)."""
    raw = _settings().value(RECENT_KEY) or []
    if isinstance(raw, str):
        raw = [raw]
    return [Path(p) for p in raw if p and Path(p).is_file()]


def remember_session(session: Path) -> None:
    keep = [str(session)] + [str(p) for p in recent_sessions() if p != Path(session)]
    _settings().setValue(RECENT_KEY, keep[:MAX_RECENT])


def session_reference(session: Path) -> Path | None:
    """The reference FASTA a session opens, when it still exists."""
    from repeatviewer.backends import genome_tracks as gt

    try:
        ref = gt.read_session(Path(session)).get("reference")
    except gt.TrackError:
        return None
    return Path(ref) if ref and Path(ref).is_file() else None


class GenomeBrowserWindow(QMainWindow):
    # another genome: a session file to open, or None to pick a FASTA
    genomeRequested = Signal(object)

    def __init__(self, session: Path, reference: Path, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.session = Path(session)
        self.setWindowTitle(f"Genome Browser — {Path(reference).name}")
        self.setObjectName("GenomeBrowserWindow")
        self.panel = ContigViewPanel(None, self)
        self.setCentralWidget(self.panel)
        self.openRequested = self.panel.openRequested
        self.blastRequested = self.panel.blastRequested
        self.fileWritten = self.panel.fileWritten
        self.reference = Path(reference)
        # Genome ▾ next to Tracks ▾, and the same in the menu bar
        self.genome_btn = QToolButton(self.panel)
        self.genome_btn.setText("Genome ▾")
        self.genome_btn.setAutoRaise(True)
        self.genome_btn.setToolTip("Open another genome, or one used before")
        self.genome_menu = QMenu(self.genome_btn)
        self.genome_menu.aboutToShow.connect(self._fill_genome_menu)
        self.genome_btn.clicked.connect(lambda: self.genome_menu.exec(
            self.genome_btn.mapToGlobal(self.genome_btn.rect().bottomLeft())))
        nr = self.panel._nav_right.layout()
        # after the search box: Genome ▾ · Tracks ▾ · Files ▾
        nr.insertWidget(nr.indexOf(self.panel._tracks_btn), self.genome_btn)
        bar_menu = self.menuBar().addMenu("&Genome")
        bar_menu.aboutToShow.connect(lambda: self._fill_genome_menu(bar_menu))
        self._fill_genome_menu(bar_menu)
        geometry = QSettings().value(_GEOMETRY_KEY)
        if geometry is None or not self.restoreGeometry(geometry):
            screen = QGuiApplication.primaryScreen()
            area = screen.availableGeometry() if screen is not None else None
            if area is not None:
                self.resize(int(area.width() * 0.85), int(area.height() * 0.85))
                self.move(area.center() - self.rect().center())
            else:
                self.resize(1400, 900)

    def _fill_genome_menu(self, menu: QMenu | None = None) -> None:
        menu = menu if isinstance(menu, QMenu) else self.genome_menu
        menu.clear()
        act = menu.addAction("Open genome…", lambda: self.genomeRequested.emit(None))
        act.setShortcut(QKeySequence.StandardKey.Open)
        recent = [p for p in recent_sessions() if session_reference(p) is not None]
        if recent:
            menu.addSeparator()
            head = menu.addAction("Recent genomes")
            head.setEnabled(False)
            for sess in recent:
                ref = session_reference(sess)
                a = menu.addAction(ref.name, lambda s_=sess: self.genomeRequested.emit(s_))
                a.setToolTip(str(ref))
                a.setCheckable(True)
                a.setChecked(sess == self.session)

    def open_session(self) -> None:
        """Show the window, then load (the first fit uses the real window width)."""
        self.show()
        self.raise_()
        self.activateWindow()
        self.panel.load(self.session)
        remember_session(self.session)

    def closeEvent(self, event) -> None:
        QSettings().setValue(_GEOMETRY_KEY, self.saveGeometry())
        super().closeEvent(event)
