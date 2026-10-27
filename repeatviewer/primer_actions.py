"""Tools > Primers in the main window (Design Primers, Map Primers, Primer Database).

Kept out of ``main_window.py``: :class:`PrimerActionsMixin` is mixed into
``MainWindow`` and only uses the window's shared job / dialog helpers
(``_present_te_viz_dialog``, ``_start_background_job``,
``_attach_process_monitor``, ``_populate_file_list`` ...).
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

from PySide6.QtCore import QObject, Signal, Slot
from PySide6.QtWidgets import QMenu, QMessageBox

from repeatviewer.backends.job_cancel import CANCELLED_MESSAGE
from repeatviewer.job_options import (
    PRIMER_DIALOG_TAB_PRIMERS,
    PRIMER_RIGHT_CLICK_DIALOG,
    PrimerOptions,
    load_primer_options,
    persistable_dialog_tab,
    save_primer_options,
)

PRIMER_JOB_DESIGN = "design"
PRIMER_JOB_MAP = "map"
PRIMER_TOOL_KEY = "primers"
_FASTA_SUFFIXES = (".fa", ".fasta", ".fna", ".fas", ".fsa", ".ffn", ".seq", ".aln", ".mfa")

# Progress phases for the bottom bar (status text, bar slice).
_PHASES = (
    ("Reading", 0),
    ("Masking", 1),
    ("Designing", 1),
    ("Scoring", 1),
    ("Testing", 1),
    ("Searching", 2),
    ("Checking", 2),
    ("Writing", 3),
)
PRIMER_JOB_STEPS = 4

# Tools > Primers shortcuts that open Design Primers on one mode.
PRIMER_MODE_MENU: tuple[tuple[str, str], ...] = (
    ("qPCR Copy Number Primers", "qpcr"),
    ("Satellite Ladder PCR Primers", "ladder"),
    ("FISH Oligo Probes", "fish"),
    ("Primers Conserved in a TE Family", "conserved"),
    ("Subfamily-Specific Primers", "specific"),
    ("TE Insertion Genotyping Primers", "insertion"),
    ("Bisulfite Primers (BSP)", "bsp"),
    ("Methylation-Specific Primers (MSP)", "msp"),
)


def is_primer_fasta(path: Path) -> bool:
    name = Path(path).name.lower()
    if name.endswith(".gz"):
        name = name[:-3]
    return name.endswith(_FASTA_SUFFIXES)


def phase_step(message: str) -> int:
    for prefix, step in _PHASES:
        if message.startswith(prefix):
            return step
    return 1


class PrimerJobWorker(QObject):
    """Design or map primers off the UI thread."""

    progress = Signal(int, str)
    log_line = Signal(str)
    finished = Signal(object)
    failed = Signal(str)

    def __init__(self, kind: str, seq_files: list[Path], settings, *, report_path: Path | None,
                 annotate: bool = True, database_path: str = "",
                 save_to_database: bool = False, cancel=None) -> None:
        super().__init__()
        self._kind = kind
        self._seq_files = list(seq_files)
        self._settings = settings
        self._report_path = report_path
        self._annotate = annotate
        self._database_path = database_path
        self._save_to_database = save_to_database
        self._cancel = cancel

    def _cancelled(self) -> bool:
        return bool(self._cancel is not None and getattr(self._cancel, "requested", False))

    @Slot()
    def run(self) -> None:
        from repeatviewer.backends.primers import jobs as pj
        from repeatviewer.backends.primers.binding import PrimerCancelled
        from repeatviewer.backends.primers.thermo import PrimerError
        from repeatviewer.external_process import process_log_sink

        def _log(message: str) -> None:
            self.log_line.emit(str(message))

        def _progress(message: str) -> None:
            self.progress.emit(phase_step(message), message)
            self.log_line.emit(message)

        db_path = Path(self._database_path).expanduser() if self._database_path else None
        try:
            with process_log_sink(_log):
                self.progress.emit(0, "Reading templates…")
                if self._kind == PRIMER_JOB_MAP:
                    main, files, payload = pj.run_map(
                        self._seq_files, self._settings, database_path=db_path,
                        progress=_progress, cancel=self._cancelled,
                        report_path=self._report_path)
                else:
                    main, files, payload = pj.run_design(
                        self._seq_files, self._settings, progress=_progress,
                        cancel=self._cancelled, annotate=self._annotate,
                        report_path=self._report_path)
                saved = 0
                if self._kind == PRIMER_JOB_DESIGN and self._save_to_database and payload.get("pairs"):
                    saved = len(pj.save_pairs_to_database(payload, db_path))
            self.progress.emit(PRIMER_JOB_STEPS - 1, "Writing results…")
            self.finished.emit({
                "kind": self._kind,
                "report": str(main),
                "files": [str(f) for f in files],
                "annotated": list(payload.get("annotated", [])),
                "pairs": len(payload.get("pairs", [])),
                "products": int((payload.get("counts") or {}).get("products", 0)),
                "sites": int((payload.get("counts") or {}).get("sites", 0)),
                "warnings": list(payload.get("warnings", [])),
                "saved_to_database": saved,
            })
        except PrimerCancelled:
            self.failed.emit(CANCELLED_MESSAGE)
        except PrimerError as exc:
            self.failed.emit(str(exc))
        except Exception as exc:  # noqa: BLE001 - surface to the UI
            self.failed.emit(f"{type(exc).__name__}: {exc}")


class PrimerActionsMixin:
    """Tools > Primers for ``MainWindow``."""

    _primer_options: PrimerOptions
    _primer_database_dialog = None

    # ------------------------------------------------------------- state

    def _init_primer_actions(self) -> None:
        self._primer_options = load_primer_options()
        self._primer_database_dialog = None

    def _save_primer_options(self, opts: PrimerOptions) -> None:
        self._primer_options = opts
        save_primer_options(opts)

    # ------------------------------------------------------------- menu

    def _populate_primers_menu(self, menu: QMenu) -> None:
        menu.addAction("Design Primers with default parameters", self._run_primer_design_default)
        menu.addAction("Change Design Primers parameters",
                       lambda: self._open_primer_design_dialog())
        menu.addSeparator()
        menu.addAction("Map Primers to Sequences", lambda: self._open_map_primers_dialog())
        menu.addAction("Primer Database", self._open_primer_database)
        menu.addSeparator()
        for label, mode in PRIMER_MODE_MENU:
            menu.addAction(label, lambda _checked=False, m=mode: self._open_primer_design_dialog(mode=m))

    # ------------------------------------------------------------- inputs

    def _primer_candidate_files(self) -> list[Path]:
        return [p for p in self._selected_files() if p.is_file() and is_primer_fasta(p)]

    def _primers_ready(self) -> bool:
        from repeatviewer.views.primer_dialogs import primer3_missing_message

        return not primer3_missing_message(self)

    # ------------------------------------------------------------- design

    def _run_primer_design_default(self) -> None:
        if not self._primers_ready():
            return
        files = self._primer_candidate_files()
        if not files:
            QMessageBox.information(
                self, "Primer Design",
                "Select the FASTA file to design primers on in the Files panel first.",
            )
            return
        self._start_primer_job(PRIMER_JOB_DESIGN, files, self._primer_options)

    def _quick_primer_design_right_click(self) -> None:
        if self._primer_options.right_click_action == PRIMER_RIGHT_CLICK_DIALOG:
            self._open_primer_design_dialog()
            return
        self._run_primer_design_default()

    def _open_primer_design_dialog(self, *, mode: str | None = None) -> None:
        if not self._primers_ready():
            return
        from repeatviewer.views.primer_dialogs import PrimerDesignDialog

        opts = self._primer_options
        if mode:
            from repeatviewer.backends.primers import design as dz

            current = opts.design_settings()
            if current.mode != mode:
                data = current.to_dict()
                data.update(dz.mode_preset(mode))
                data["mode"] = mode
                opts = replace(opts, design=data, last_dialog_tab=PRIMER_DIALOG_TAB_PRIMERS)
        dialog = PrimerDesignDialog(
            initial=opts,
            files=self._primer_candidate_files(),
            project_genome=self._genome_path_resolved(),
            parent=self,
        )

        def _accepted(d) -> None:
            self._save_primer_options(d.options())
            if not d.should_run():
                self._set_status_message("Primer Design settings saved", 4000)
                return
            files = d.template_files() or self._primer_candidate_files()
            if not files:
                QMessageBox.information(
                    self, "Primer Design",
                    "Click the FASTA to design on in the Files panel, then Run.",
                )
                return
            self._start_primer_job(PRIMER_JOB_DESIGN, files, self._primer_options)

        self._present_te_viz_dialog(dialog, _accepted)

    # ------------------------------------------------------------- map

    def _open_map_primers_dialog(self) -> None:
        if not self._primers_ready():
            return
        from repeatviewer.views.primer_dialogs import MapPrimersDialog

        dialog = MapPrimersDialog(
            initial=self._primer_options,
            files=self._primer_candidate_files(),
            project_genome=self._genome_path_resolved(),
            parent=self,
        )

        def _accepted(d) -> None:
            self._save_primer_options(d.options())
            if not d.should_run():
                self._set_status_message("Map Primers settings saved", 4000)
                return
            files = d.run_files() or self._primer_candidate_files()
            if not files:
                QMessageBox.information(
                    self, "Map Primers",
                    "Click the sequences to test in the Files panel, then Run.",
                )
                return
            self._start_primer_job(PRIMER_JOB_MAP, files, self._primer_options)

        self._present_te_viz_dialog(dialog, _accepted)

    def _quick_map_primers_right_click(self) -> None:
        if self._primer_options.map_right_click_action == PRIMER_RIGHT_CLICK_DIALOG:
            self._open_map_primers_dialog()
            return
        if not self._primers_ready():
            return
        files = self._primer_candidate_files()
        if not files:
            self._open_map_primers_dialog()
            return
        self._start_primer_job(PRIMER_JOB_MAP, files, self._primer_options)

    # ------------------------------------------------------------- database

    def _open_primer_database(self) -> None:
        if not self._primers_ready():
            return
        from repeatviewer.views.primer_dialogs import PrimerDatabaseDialog

        existing = self._primer_database_dialog
        if existing is not None:
            try:
                if existing.isVisible():
                    existing.raise_()
                    existing.activateWindow()
                    return
            except RuntimeError:
                self._primer_database_dialog = None
        path = self._primer_options.database_path or None
        dialog = PrimerDatabaseDialog(
            Path(path) if path else None,
            thermo=self._primer_options.map_settings().thermo,
            parent=self,
        )
        self._primer_database_dialog = dialog

        def _closed(_code: int) -> None:
            if self._primer_database_dialog is dialog:
                self._primer_database_dialog = None
            dialog.deleteLater()

        dialog.finished.connect(_closed)
        dialog.show()
        dialog.raise_()
        dialog.activateWindow()

    # ------------------------------------------------------------- tab memory

    def _remember_primer_dialog_tab(self, dialog) -> bool:
        """Cancel keeps the open tab (like the other tool windows)."""
        from repeatviewer.views.primer_dialogs import MapPrimersDialog, PrimerDesignDialog

        if isinstance(dialog, PrimerDesignDialog):
            tab = persistable_dialog_tab(dialog.current_dialog_tab(),
                                         first_tab=PRIMER_DIALOG_TAB_PRIMERS)
            self._save_primer_options(replace(self._primer_options, last_dialog_tab=tab))
            return True
        return isinstance(dialog, MapPrimersDialog)

    # ------------------------------------------------------------- jobs

    def _confirm_primer_report_path(self, kind: str, files: list[Path]) -> Path | None | bool:
        """Rule 9: ask before an existing report is replaced. False = cancelled."""
        from repeatviewer.backends.primers import jobs as pj
        from repeatviewer.widgets.file_conflict import confirm_existing_file_path

        tag = pj.MAP_TAG if kind == PRIMER_JOB_MAP else pj.DESIGN_TAG
        proposed = pj.default_report_path(files[0], tag)
        title = "Map Primers" if kind == PRIMER_JOB_MAP else "Primer Design"
        chosen = confirm_existing_file_path(self, proposed, title=title,
                                            required_suffix=pj.REPORT_SUFFIX)
        if chosen is None:
            return False
        return chosen

    def _start_primer_job(self, kind: str, files: list[Path], opts: PrimerOptions) -> None:
        files = [Path(p) for p in files]
        report_path = self._confirm_primer_report_path(kind, files)
        if report_path is False:
            return
        if kind == PRIMER_JOB_MAP:
            settings = opts.map_settings()
            title = (f"Map Primers on {files[0].name}" if len(files) == 1
                     else f"Map Primers on {len(files)} files")
        else:
            from repeatviewer.backends.primers import design as dz

            settings = opts.design_settings()
            label = dz.MODE_LABELS.get(settings.mode, "Primers")
            title = f"Primer Design ({label}) on {files[0].name}"
        worker = PrimerJobWorker(
            kind, files, settings,
            report_path=report_path,
            annotate=(opts.annotate_design if kind == PRIMER_JOB_DESIGN
                      else bool(getattr(settings, "annotate", True))),
            database_path=opts.database_path,
            save_to_database=opts.save_design_to_database,
            cancel=self._job_cancel,
        )
        self._attach_process_monitor(
            worker,
            title=title,
            enabled=bool(opts.show_process_window),
            tool_key=PRIMER_TOOL_KEY,
            mirror_progress=False,
            close_when_finished=bool(opts.close_process_window_when_finished),
        )
        self._start_background_job(
            worker,
            job_title=title,
            total_steps=PRIMER_JOB_STEPS,
            source_files=files,
            on_finished=self._on_primer_job_finished,
            on_failed=self._on_primer_job_failed,
        )

    @Slot(object)
    def _on_primer_job_finished(self, result) -> None:
        report = Path(result["report"])
        created = [report] if report.is_file() else []
        if self._current_folder:
            self._refresh_subfolders_under(self._current_folder)
            self._populate_file_list(
                self._current_folder,
                new_outputs=created,
                source_files=getattr(self, "_last_job_source_files", []),
                open_preview=False,
            )
        if created:
            self._auto_open_outputs(created)
        if result.get("kind") == PRIMER_JOB_MAP:
            msg = (f"Map Primers finished: {result.get('sites', 0):,} binding sites, "
                   f"{result.get('products', 0):,} products ({report.name})")
        else:
            n = int(result.get("pairs", 0))
            msg = f"Primer Design finished: {n} result{'s' if n != 1 else ''} ({report.name})"
            if result.get("saved_to_database"):
                msg += f"; {result['saved_to_database']} primer(s) added to the Primer Database"
        self._set_status_message(msg, 12000, override_file_link=True)
        if result.get("annotated"):
            self._reload_annotated_previews(getattr(self, "_last_job_source_files", []))

    def _reload_annotated_previews(self, fastas: list[Path]) -> None:
        """Open previews of the annotated inputs show the new primer arrows.

        Pages with unsaved edits are left alone (``reload_from_disk`` skips them).
        """
        tabs = getattr(self, "_preview_tabs", None)
        if tabs is None or not fastas:
            return
        from repeatviewer.widgets.preview_tabs import FilePreviewPage

        keys = set()
        for p in fastas:
            try:
                keys.add(str(Path(p).expanduser().resolve()))
            except OSError:
                continue
        current = tabs.currentWidget()
        for index in range(tabs.count()):
            page = tabs.widget(index)
            # The active tab now shows the report; only background tabs that
            # still show an annotated input need the new arrows.
            if page is current or not isinstance(page, FilePreviewPage):
                continue
            if str(getattr(page, "path", "")) in keys:
                try:
                    page.reload_from_disk()
                except RuntimeError:
                    continue

    @Slot(str)
    def _on_primer_job_failed(self, error: str) -> None:
        if error.strip() == CANCELLED_MESSAGE:
            return
        QMessageBox.warning(self, "Primers", error.strip() or "Unknown error.")
