"""Tools > Genome Browser: a reference FASTA with annotation and alignment tracks in a
window of its own (the contig view in browser mode, without the Preview panel)."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QSettings, Signal
from PySide6.QtGui import QGuiApplication, QKeySequence
from PySide6.QtWidgets import (
    QFileDialog,
    QMainWindow,
    QMenu,
    QMessageBox,
    QToolButton,
    QWidget,
)

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
        # the genome sits in the breadcrumb, after Back / Forward and before the
        # sequence: its name is the drop-down (the old Genome ▾ button stays hidden)
        self.genome_btn.setVisible(False)
        self.panel._crumbs.set_genome(self.reference.name, self.genome_menu)
        # the menu bar: File (genomes, the files) and Tracks; the toolbar has neither
        self.file_menu = self.menuBar().addMenu("&File")
        self.file_menu.aboutToShow.connect(self._fill_file_menu)
        self._fill_file_menu()
        self.tracks_menu = self.menuBar().addMenu("&Tracks")
        self.tracks_menu.aboutToShow.connect(
            lambda: self.panel.fill_tracks_menu(self.tracks_menu))
        self.panel.fill_tracks_menu(self.tracks_menu)
        # Bookmarks: add the view, import / export BED, then the bookmarks to go to
        self.bookmarks_menu = self.menuBar().addMenu("&Bookmarks")
        self.bookmarks_menu.aboutToShow.connect(self._fill_bookmarks_menu)
        self._fill_bookmarks_menu()
        geometry = QSettings().value(_GEOMETRY_KEY)
        if geometry is None or not self.restoreGeometry(geometry):
            screen = QGuiApplication.primaryScreen()
            area = screen.availableGeometry() if screen is not None else None
            if area is not None:
                self.resize(int(area.width() * 0.85), int(area.height() * 0.85))
                self.move(area.center() - self.rect().center())
            else:
                self.resize(1400, 900)

    def _fill_file_menu(self) -> None:
        """File: open a genome, the recent ones, then the files of this session."""
        menu = self.file_menu
        self._fill_genome_menu(menu)
        menu.addSeparator()
        menu.addAction("Open Session…", lambda: self.open_session_file())
        act = menu.addAction("Save Session…", lambda: self.save_session_as())
        act.setShortcut(QKeySequence.StandardKey.SaveAs)
        files = self.panel._files_menu
        files.setTitle("Files")
        if files.actions():
            menu.addSeparator()
            menu.addMenu(files)
        menu.addSeparator()
        menu.addAction("Close window", self.close, QKeySequence.StandardKey.Close)

    # ---------------------------------------------------------------- sessions

    _SESSION_FILTER = "Genome browser session (*.rvgb);;All files (*)"

    def open_session_file(self, path: str | None = None) -> None:
        """File > Open Session…: a saved .rvgb (genome, tracks, view, bookmarks)."""
        if path is None:
            path, _ = QFileDialog.getOpenFileName(
                self, "Open Session", str(self.session.parent), self._SESSION_FILTER)
        if not path:
            return
        if session_reference(Path(path)) is None:
            QMessageBox.warning(self, "Open Session",
                                f"{Path(path).name} is not a genome browser session, or "
                                "its genome file has moved.")
            return
        if Path(path) == self.session:
            return
        self.genomeRequested.emit(Path(path))

    def save_session_as(self, path: str | None = None) -> Path | None:
        """File > Save Session…: the genome, tracks (with their heights and colours),
        the place in view and the bookmarks, to a .rvgb file; this window then keeps
        saving there."""
        from repeatviewer.backends import genome_tracks as gt

        if path is None:
            start = self.session.parent / f"{self.reference.name}_session.rvgb"
            path, _ = QFileDialog.getSaveFileName(self, "Save Session", str(start),
                                                  self._SESSION_FILTER)
        if not path:
            return None
        out = Path(path)
        if out.suffix.lower() != gt.SESSION_SUFFIX:
            out = out.with_name(out.name + gt.SESSION_SUFFIX)
        panel = self.panel
        extras = {"bookmarks": [b.to_json() for b in panel._bookmarks],
                  "recent": [p_.to_json() for p_ in panel._recent]}
        try:
            gt.write_session(out, self.reference, [t.to_json() for t in panel.canvas.tracks],
                             panel._current_locus(), extras)
        except OSError as exc:
            QMessageBox.warning(self, "Save Session", f"Cannot save {out.name}: {exc}")
            return None
        self.session = out
        panel._path = out
        remember_session(out)
        return out

    def _fill_bookmarks_menu(self) -> None:
        menu = self.bookmarks_menu
        menu.clear()
        panel = self.panel
        menu.addAction("Add bookmark for the view (B)",
                       lambda: panel._add_bookmark(whole_view=True))
        menu.addSeparator()
        menu.addAction("Import BED…", lambda: panel._import_bookmarks())
        act = menu.addAction("Export BED…", lambda: panel._export_bookmarks())
        act.setEnabled(bool(panel._bookmarks))
        if panel._bookmarks:
            menu.addSeparator()
            for bm in panel._bookmarks[:30]:
                menu.addAction(f"{bm.name or bm.locus()}",
                               lambda b=bm: panel._nav_go(b.chrom, b.start, b.end, b.name))

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
