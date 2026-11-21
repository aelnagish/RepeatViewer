"""Tools > Primers and Probes > PCR (with gel, binding sites, bench card) for ``MainWindow``.

Developer: Abdullah Yusuf
"""

from __future__ import annotations

import json
from pathlib import Path

from PySide6.QtCore import QObject, QSettings, Qt, Signal, Slot
from PySide6.QtWidgets import QMessageBox

from repeatviewer.backends.job_cancel import CANCELLED_MESSAGE


def sanger_read_suffixes() -> tuple[str, ...]:
    from repeatviewer.backends.primers import sanger as sg

    return sg.AB1_SUFFIXES

PCR_SETTINGS_GROUP = "pcr"
PCR_JOB_STEPS = 4


def load_pcr_settings():
    from repeatviewer.backends.primers.pcr import PcrSettings

    s = QSettings()
    s.beginGroup(PCR_SETTINGS_GROUP)
    raw = s.value("settings_json", "")
    rc = s.value("right_click_runs", False)
    s.endGroup()
    try:
        data = json.loads(str(raw)) if raw else {}
    except ValueError:
        data = {}
    if isinstance(rc, str):
        rc = rc.strip().lower() in {"1", "true", "yes"}
    return PcrSettings.from_dict(data), bool(rc)


def save_pcr_settings(settings, right_click_runs: bool) -> None:
    s = QSettings()
    s.beginGroup(PCR_SETTINGS_GROUP)
    s.setValue("settings_json", json.dumps(settings.to_dict()))
    s.setValue("right_click_runs", bool(right_click_runs))
    s.endGroup()


class PcrJobWorker(QObject):
    """In-silico PCR off the UI thread; writes the products next to the first template."""

    progress = Signal(int, str)
    log_line = Signal(str)
    finished = Signal(object)
    failed = Signal(str)

    def __init__(self, templates: list[Path], primers: list, settings, *, cancel=None) -> None:
        super().__init__()
        self._seq_files = [Path(t) for t in templates]
        self._primers = list(primers)
        self._settings = settings
        self._cancel = cancel

    def _cancelled(self) -> bool:
        return bool(self._cancel is not None and getattr(self._cancel, "requested", False))

    @Slot()
    def run(self) -> None:
        from repeatviewer.backends.primers import pcr
        from repeatviewer.backends.primers.binding import PrimerCancelled
        from repeatviewer.backends.primers.thermo import PrimerError

        def _progress(message: str) -> None:
            self.progress.emit(1, message)
            self.log_line.emit(message)

        try:
            self.progress.emit(0, "Running PCR…")
            result = pcr.run_pcr(self._seq_files, self._primers, self._settings,
                                 progress=_progress, cancel=self._cancelled)
            self.progress.emit(2, "Writing products…")
            outputs = pcr.write_outputs(result)
            self.progress.emit(3, "Done")
            self.finished.emit({"result": result, "outputs": outputs})
        except PrimerCancelled:
            self.failed.emit(CANCELLED_MESSAGE)
        except PrimerError as exc:
            self.failed.emit(str(exc))
        except Exception as exc:  # noqa: BLE001 - surface to the UI
            self.failed.emit(f"{type(exc).__name__}: {exc}")


class PcrActionsMixin:
    _primer_genome_dialog = None
    _bench_card_dialog = None

    # ------------------------------------------------------------- menu

    def _populate_pcr_menu_items(self, menu) -> None:
        a = menu.addAction("In-silico PCR", self._open_pcr_in_panel)
        a.setToolTip("In-silico PCR: products of the primers on this sequence or the ticked files")
        a = menu.addAction("Restriction Digest", self._open_digest_in_panel)
        a.setToolTip("Restriction map and digest of the sequence, with your methylation data")
        a = menu.addAction("Gel", self._open_gel_document)
        a.setToolTip("A gel: open the selected one, or start an empty one")
        menu.addSeparator()
        a = menu.addAction("Sanger Check", lambda: self._open_sanger_dialog())
        a.setToolTip("Sanger reads (.ab1) of a clone against the sequence it should carry")
        menu.addSeparator()
        a = menu.addAction("Cloning", lambda: self._open_cloning_dialog())
        a.setToolTip("Put a PCR product into a vector: orientation, colony PCR, warnings")

    def _open_pcr_in_panel(self) -> None:
        """PCR page of the Primers panel, live on the open sequence."""
        self._open_primers_panel(tab="pcr")

    def _open_digest_in_panel(self) -> None:
        """Enzymes page of the Primers panel, on its Digest view."""
        self._open_primers_panel(tab="enzymes", digest=True)

    def _open_gel_document(self) -> None:
        """Open the selected gel, or make an empty one to fill in."""
        from repeatviewer.backends.primers import gel_doc as gd

        chosen = [p for p in self._selected_files()
                  if Path(p).is_file() and Path(p).suffix.lower() == gd.GEL_SUFFIX]
        if chosen:
            self._show_new_output(Path(chosen[0]))
            return
        self._new_gel_document()

    # ------------------------------------------------------------- cloning

    def _open_cloning_dialog(self, *, insert: str = "", insert_name: str = "",
                             path=None, features=()) -> None:
        """Cloning as a tab in the middle, with its settings on the right."""
        from repeatviewer.views.cloning_page import CloningPage
        from repeatviewer.views.sanger_dialog import _read_fasta

        if not insert:
            source = Path(path) if path else next(
                (p for p in (self._selected_files() if hasattr(self, "_selected_files") else [])
                 if Path(p).is_file()), None)
            if source is not None:
                records = _read_fasta(Path(source))
                if records:
                    insert_name, insert = records[0][0], records[0][1]
                    path = source
        folder = Path(path).parent if path else None
        page = CloningPage(insert, insert_name=insert_name, folder=folder, features=features)
        if getattr(self, "_preview_tabs", None) is None:
            return
        self._host_in_workspace("cloning", page)
        self._cloning_page = page

    def _open_sanger_dialog(self, *, expected: str = "", expected_name: str = "",
                            path=None) -> None:
        """Sanger Check as a tab in the middle, with its settings on the right."""
        from repeatviewer.views.sanger_dialog import _read_fasta
        from repeatviewer.views.sanger_page import SangerPage

        picked = [Path(p) for p in self._selected_files()] if hasattr(self, "_selected_files") \
            else []
        reads = [p for p in picked if p.suffix.lower() in sanger_read_suffixes()]
        if not expected:
            source = Path(path) if path else next(
                (p for p in picked if p not in reads and p.is_file()), None)
            if source is not None:
                records = _read_fasta(source)
                if records:
                    expected_name, expected = records[0][0], records[0][1]
                    path = source
        folder = Path(path).parent if path else (reads[0].parent if reads else None)
        page = SangerPage(expected, expected_name=expected_name, folder=folder)
        if reads:
            page.add_reads(reads)
        if getattr(self, "_preview_tabs", None) is None:
            return
        self._host_in_workspace("sanger", page)
        self._sanger_page = page

    # ------------------------------------------------------------- templates

    def _pcr_template_candidates(self) -> list[Path]:
        try:
            return self._primer_candidate_files()
        except Exception:  # noqa: BLE001
            return []

    def _primers_of_files(self, files: list[Path]) -> list:
        """Primers kept with the files (Primers tab); ticked ones are not known here."""
        from repeatviewer.backends.primers.on_sequence import load_sequence_primers

        out, names = [], set()
        for f in files:
            for p in load_sequence_primers(f):
                if p.name not in names:
                    out.append(p)
                    names.add(p.name)
        return out

    def _quick_pcr_right_click(self) -> None:
        settings, rc = load_pcr_settings()
        files = self._pcr_template_candidates()
        primers = self._primers_of_files(files)
        if rc and files and primers:
            self._start_pcr_job(files, primers, settings)
            return
        self._open_pcr_in_panel()

    def _sync_primer_genome_from_files(self, clicked_path: Path | None) -> None:
        """A click in the Files panel fills the path of the open Binding Sites window."""
        if clicked_path is None:
            return
        gdlg = self._primer_genome_dialog
        if gdlg is None:
            return
        try:
            if gdlg.isVisible():
                gdlg.apply_path_from_panel(clicked_path)
        except RuntimeError:
            self._primer_genome_dialog = None

    def _start_pcr_job(self, templates, primers, settings) -> None:
        files = [Path(t) for t in templates]
        worker = PcrJobWorker(files, list(primers), settings, cancel=self._job_cancel)
        title = f"PCR on {files[0].name}" if len(files) == 1 else f"PCR on {len(files)} files"
        self._last_pcr_sources = files
        self._start_background_job(
            worker,
            job_title=title,
            total_steps=PCR_JOB_STEPS,
            source_files=files,
            on_finished=self._on_pcr_finished,
            on_failed=self._on_pcr_failed,
        )

    def _on_pcr_finished(self, payload) -> None:
        result = payload["result"]
        outputs = payload["outputs"]
        created = [p for p in outputs.get("files", []) if Path(p).is_file()]
        fasta = outputs.get("fasta")
        if self._current_folder:
            self._refresh_subfolders_under(self._current_folder)
            self._populate_file_list(
                self._current_folder,
                new_outputs=[Path(p) for p in created],
                source_files=getattr(self, "_last_pcr_sources", []),
                open_preview=False,
            )
        if fasta:
            self._auto_open_outputs([Path(fasta)])
        from repeatviewer.backends.primers import pcr

        self._set_status_message(f"PCR finished: {pcr.describe(result)}", 10000)

    def _on_pcr_failed(self, error: str) -> None:
        if error.strip() == CANCELLED_MESSAGE:
            return
        QMessageBox.warning(self, "PCR", error.strip() or "Unknown error.")

    # ------------------------------------------------------------- gel

    def _open_gel_for_pcr(self, result) -> None:
        from repeatviewer.backends.primers import gel_doc as gd
        from repeatviewer.backends.primers import pcr

        first = Path(result.templates[0]) if result.templates else Path.cwd()
        base = pcr.output_base(first, result.primers)
        path = first.parent / f"{base}{gd.GEL_SUFFIX}"
        if path.exists():
            # Same run: reuse; another run with this base: new name.
            try:
                same = gd.GelDoc.load(path).title == base
            except Exception:  # noqa: BLE001
                same = False
            if not same:
                path = gd.unique_gel_path(first.parent, base)
        if not path.exists():
            doc = gd.gel_from_pcr(result)
            doc.title = base
            try:
                doc.save(path)
            except OSError as exc:
                QMessageBox.warning(self, "Gel", f"Could not write the gel:\n{exc}")
                return
        self._show_new_output(path)

    def _new_gel_document(self) -> None:
        from repeatviewer.backends.primers import gel_doc as gd

        folder = Path(self._current_folder) if self._current_folder else Path.home()
        # An empty gel is no file yet: it is written when the user saves it,
        # so clicking Gel does not leave Gel.rvgel files behind.
        ensure = getattr(self, "_ensure_primers_workspace", None)
        workspace = ensure() if ensure is not None else None
        tabs = getattr(self, "_preview_tabs", None)
        if workspace is not None and tabs is not None and hasattr(workspace, "show_new_gel"):
            workspace.show_new_gel(folder)
            tabs.setCurrentWidget(workspace)
            return
        path = gd.unique_gel_path(folder, "Gel")
        doc = gd.GelDoc()
        doc.lanes = [gd.ladder_lane(doc.ladder)]
        doc.title = path.stem
        try:
            doc.save(path)
        except OSError as exc:
            QMessageBox.warning(self, "Gel", f"Could not write the gel:\n{exc}")
            return
        self._show_new_output(path)

    def _show_new_output(self, path: Path) -> None:
        if self._current_folder and Path(path).parent == Path(self._current_folder):
            self._populate_file_list(self._current_folder, new_outputs=[path],
                                     source_files=[], open_preview=False)
        if Path(path).suffix.lower() == ".rvgel":
            # a gel belongs to the Gel section of the Primers and Probes workspace
            ensure = getattr(self, "_ensure_primers_workspace", None)
            workspace = ensure() if ensure is not None else None
            tabs = getattr(self, "_preview_tabs", None)
            if workspace is not None and tabs is not None:
                workspace.show_gel(Path(path))
                tabs.setCurrentWidget(workspace)
                self._ensure_preview_visible()
                return
        self._open_preview_tab(path, replace=False)

    # ------------------------------------------------------------- genome view, bench card

    def _open_primer_genome_view(self, primers: list, *, genome=None, parent=None) -> None:
        """Binding Sites window: from the PCR window (Sites in Genome…) or the Primers tab."""
        if not self._primers_ready():
            return
        from repeatviewer.views.primer_genome_dialog import PrimerGenomeDialog

        if genome is None:
            try:
                genome = self._genome_path_resolved()
            except Exception:  # noqa: BLE001
                genome = None
        if not primers:
            primers = self._primers_of_files(self._pcr_template_candidates())
        dlg = PrimerGenomeDialog(primers, genome=genome, parent=parent or self)
        dlg.setWindowFlag(Qt.WindowType.Tool, True)
        dlg.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, True)
        dlg.fileWritten.connect(lambda p: self._show_new_output(Path(p)))
        self._primer_genome_dialog = dlg
        dlg.show()
        dlg.raise_()

    def _open_bench_card(self, result) -> None:
        from repeatviewer.views.bench_card_dialog import BenchCardDialog

        primers = list(result.primers) if result is not None else \
            self._primers_of_files(self._pcr_template_candidates())
        db = self._primer_options.database_path or None
        dlg = BenchCardDialog(primers, result=result, database_path=Path(db) if db else None,
                              parent=self)
        dlg.setWindowFlag(Qt.WindowType.Tool, True)
        dlg.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, True)
        self._bench_card_dialog = dlg
        dlg.show()
        dlg.raise_()
