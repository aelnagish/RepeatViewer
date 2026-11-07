"""Tools > Primers in the main window (Design Primers, Map Primers, Primer Database).

Kept out of ``main_window.py``: :class:`PrimerActionsMixin` is mixed into
``MainWindow`` and only uses the window's shared job / dialog helpers
(``_present_te_viz_dialog``, ``_start_background_job``,
``_attach_process_monitor``, ``_populate_file_list`` ...).
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

from PySide6.QtCore import QObject, Qt, QTimer, Signal, Slot
from PySide6.QtWidgets import QMenu, QMessageBox

from repeatviewer.backends.job_cancel import CANCELLED_MESSAGE
from repeatviewer.job_options import (
    PRIMER_DIALOG_TAB_PRIMERS,
    PrimerOptions,
    load_primer_options,
    persistable_dialog_tab,
    save_primer_options,
)

PRIMER_JOB_DESIGN = "design"
PRIMER_JOB_MAP = "map"
PRIMER_JOB_PAINTING = "painting"
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
    ("Running Chorus2", 1),
    ("Selecting", 2),
    ("Writing", 3),
)
PRIMER_JOB_STEPS = 4

# Tools > Primers and Probes: three items; the design type is the Mode inside Design.
PRIMER_MENU_ITEMS: tuple[str, ...] = ("Design", "Map", "Database")
PRIMER_MANUAL = Path(__file__).resolve().parent / "assets" / "manuals" / "primers_and_probes.html"


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
                if self._kind == PRIMER_JOB_PAINTING:
                    from repeatviewer.backends.primers import painting as pt

                    main, files, payload = pt.run_painting(
                        self._settings, report_path=self._report_path,
                        cancel=self._cancel, progress=_progress)
                elif self._kind == PRIMER_JOB_MAP:
                    main, files, payload = pj.run_map(
                        self._seq_files, self._settings, database_path=db_path,
                        progress=_progress, cancel=self._cancelled,
                        report_path=self._report_path)
                else:
                    # Nothing is written onto the sequences here: the user picks
                    # the pairs to add afterwards (Designed Primers window).
                    main, files, payload = pj.run_design(
                        self._seq_files, self._settings, progress=_progress,
                        cancel=self._cancelled, annotate=False,
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
                "selected": int(payload.get("selected", 0)),
                "chorus_input": int(payload.get("input", 0)),
                "products": int((payload.get("counts") or {}).get("products", 0)),
                "sites": int((payload.get("counts") or {}).get("sites", 0)),
                "warnings": list(payload.get("warnings", [])),
                "saved_to_database": saved,
                "pair_list": list(payload.get("pairs", [])) if self._kind == PRIMER_JOB_DESIGN else [],
                "ask_pairs": bool(self._annotate) if self._kind == PRIMER_JOB_DESIGN else False,
                "inputs": [str(p) for p in self._seq_files],
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
    _primer_check_dialog = None

    # ------------------------------------------------------------- state

    def _init_primer_actions(self) -> None:
        self._primer_options = load_primer_options()
        self._primer_database_dialog = None
        self._primer_check_dialog = None
        try:
            from repeatviewer.views.primer_database_dialog import apply_lab_user_setting

            apply_lab_user_setting()
        except Exception:  # noqa: BLE001 - only the name in the history
            pass

    def _save_primer_options(self, opts: PrimerOptions) -> None:
        self._primer_options = opts
        save_primer_options(opts)

    # ------------------------------------------------------------- menu

    def _populate_primers_menu(self, menu: QMenu) -> None:
        design, map_, database = PRIMER_MENU_ITEMS
        a = menu.addAction(design, lambda: self._open_primer_design_dialog())
        a.setToolTip("Design PCR primers or FISH probes (choose the type in Mode)")
        a = menu.addAction(map_, lambda: self._open_map_primers_dialog())
        a.setToolTip("Test primers against sequences: binding sites and products")
        a = menu.addAction(database, self._open_primer_database)
        a.setToolTip("Your saved primers, with Tm, hairpin and dimer checks")
        populate = getattr(self, "_populate_pcr_menu_items", None)
        if populate is not None:
            populate(menu)
        menu.setToolTipsVisible(True)

    def _open_primers_manual(self) -> None:
        from PySide6.QtCore import QUrl
        from PySide6.QtGui import QDesktopServices

        if not PRIMER_MANUAL.is_file():
            QMessageBox.warning(self, "Primers and Probes Manual",
                                f"The manual is missing:\n{PRIMER_MANUAL}")
            return
        QDesktopServices.openUrl(QUrl.fromLocalFile(str(PRIMER_MANUAL)))

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
        from repeatviewer.backends.primers import design as dz

        if self._primer_options.design_settings().mode == dz.MODE_PAINTING:
            self._start_painting_job(self._primer_options, open_dialog_on_error=True)
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
        self._run_primer_design_default()

    def _open_primer_design_dialog(self, *, mode: str | None = None,
                                   files: list[Path] | None = None) -> None:
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
        start_files = list(files or [])
        dialog = PrimerDesignDialog(
            initial=opts,
            files=start_files or self._primer_candidate_files(),
            project_genome=self._genome_path_resolved(),
            parent=self,
        )

        def _accepted(d) -> None:
            self._save_primer_options(d.options())
            if not d.should_run():
                self._set_status_message("Primer Design settings saved", 4000)
                return
            if d.is_painting():
                self._start_painting_job(self._primer_options)
                return
            # What the window shows; with Use selected and nothing shown yet, the panel.
            files = d.template_files() or (
                (start_files or self._primer_candidate_files())
                if d.uses_selected() else [])
            if not files:
                QMessageBox.information(
                    self, "Primer Design",
                    "Choose the template: click a FASTA in the Files panel (or fill "
                    "Template with Define inputs), then Run.",
                )
                return
            missing = [p for p in files if not Path(p).is_file()]
            if missing:
                QMessageBox.warning(self, "Primer Design",
                                    f"The template does not exist:\n{missing[0]}")
                return
            self._start_primer_job(PRIMER_JOB_DESIGN, files, self._primer_options)

        self._present_te_viz_dialog(dialog, _accepted)

    # ------------------------------------------------------------- map

    @staticmethod
    def _selected_primer_slots(database_path: str, names: list[str]) -> dict | None:
        """Forward / Reverse / Probe for the primers selected in the Database, or None."""
        if not names:
            return None
        from repeatviewer.backends.primers import inventory as inv
        from repeatviewer.backends.primers import jobs as pj
        from repeatviewer.backends.primers.thermo import PrimerError

        try:
            db = inv.open_database(database_path or None)
        except (OSError, PrimerError):
            return None
        wanted = set(names)
        chosen = [p for p in db.primers if p.name in wanted]
        if len(chosen) != len(wanted):
            return None
        return pj.selection_slots(chosen)

    def _check_stock_map_data(self, primer_names: list[str], database_path: str,
                              strict: dict | None = None) -> dict:
        """Map Primers settings for "test these Database primers on a sequence"."""
        from repeatviewer.backends.primers import jobs as pj

        ms = self._primer_options.map_settings()
        data = ms.to_dict()
        data.update(mode=pj.MAP_DATABASE, primer_file="", search_forward=True,
                    search_reverse=True, search_probe=False, pairs_only=False,
                    use_max_mismatches=True, max_mismatches=max(2, ms.max_mismatches),
                    three_prime_exact=3, use_product_range=True, product_min=40,
                    product_max=max(3000, ms.product_max), test_on=pj.TEST_ON_SEQUENCES,
                    primer_names=list(primer_names or []))
        for key in ("product_min", "product_max", "max_mismatches", "three_prime_exact"):
            if strict and key in strict:
                data[key] = int(strict[key])
        slots = self._selected_primer_slots(database_path, primer_names or [])
        if slots:
            # A pair (or one primer, or a pair and a probe) goes straight into the
            # Forward / Reverse / Probe boxes, so the window shows what is tested.
            data.update(mode=pj.MAP_SPECIFIC, primer_names=[],
                        use_forward=pj.KIND_FORWARD in slots,
                        forward=slots.get(pj.KIND_FORWARD, ""),
                        use_reverse=pj.KIND_REVERSE in slots,
                        reverse=slots.get(pj.KIND_REVERSE, ""),
                        use_probe=pj.KIND_PROBE in slots,
                        probe=slots.get(pj.KIND_PROBE, ""))
        return data

    def _open_primer_check(self, database_path: str, primer_names: list[str]) -> None:
        """Database › Check a Sequence: the small window for the selected primers."""
        if not self._primers_ready():
            return
        from repeatviewer.backends.primers import inventory as inv
        from repeatviewer.backends.primers.thermo import PrimerError
        from repeatviewer.views.primer_check_dialog import CheckPrimersDialog

        old = self._primer_check_dialog
        if old is not None:
            try:
                old.close()
            except RuntimeError:
                pass
        try:
            db = inv.open_database(database_path or None)
        except (OSError, PrimerError) as exc:
            QMessageBox.warning(self, "Check Primers", str(exc))
            return
        owner = self._primer_database_dialog or self
        files = self._primer_candidate_files()
        dialog = CheckPrimersDialog(
            db, list(primer_names), thermo=self._primer_options.map_settings().thermo,
            sequence=files[0] if files else None, parent=owner)
        self._primer_check_dialog = dialog

        def _run(values: dict) -> None:
            data = self._check_stock_map_data(values["primer_names"], database_path, values)
            opts = replace(self._primer_options, map=data)
            self._start_primer_job(PRIMER_JOB_MAP, [Path(values["sequence"])], opts)

        def _all(values: dict) -> None:
            QTimer.singleShot(0, lambda: self._open_map_primers_dialog(
                check_stock=True, primer_names=values["primer_names"],
                database_path=database_path, strict=values,
                sequence=values.get("sequence") or ""))

        def _closed(_code: int) -> None:
            if self._primer_check_dialog is dialog:
                self._primer_check_dialog = None
            dialog.deleteLater()

        dialog.runRequested.connect(_run)
        dialog.allSettingsRequested.connect(_all)
        dialog.finished.connect(_closed)
        dialog.show()
        dialog.raise_()
        dialog.activateWindow()

    def _open_map_primers_dialog(self, *, check_stock: bool = False,
                                 primer_names: list[str] | None = None,
                                 database_path: str = "", strict: dict | None = None,
                                 sequence: str = "") -> None:
        """Map Primers; ``check_stock``: preset to test the Database primers (the
        selected ones, or all) on a sequence."""
        if not self._primers_ready():
            return
        from repeatviewer.views.primer_dialogs import MapPrimersDialog

        opts = self._primer_options
        if check_stock:
            data = self._check_stock_map_data(list(primer_names or []), database_path, strict)
            opts = replace(opts, map=data)
        elif opts.map.get("primer_names"):
            # A selection from the Database window is for that one check only.
            opts = replace(opts, map={**opts.map, "primer_names": []})
        dialog = MapPrimersDialog(
            initial=opts,
            files=self._primer_candidate_files(),
            project_genome=self._genome_path_resolved(),
            parent=self,
        )
        if sequence:
            dialog._input_define.setChecked(True)
            dialog._seq_row.set_path(Path(sequence))

        def _no_sequence(d) -> str:
            files = d.run_files() or (
                self._primer_candidate_files() if d.uses_selected() else [])
            return "" if files else (
                "Choose the sequences: click a FASTA in the Files panel (or fill "
                "Sequences with Define inputs), then Run.")

        # Checked before the window closes, so OK brings you back to it.
        dialog._run_check = _no_sequence

        def _forget_selection() -> None:
            if self._primer_options.map.get("primer_names"):
                self._save_primer_options(replace(
                    self._primer_options, map={**self._primer_options.map, "primer_names": []}))

        def _accepted(d) -> None:
            self._save_primer_options(d.options())
            if not d.should_run():
                _forget_selection()
                self._set_status_message("Map Primers settings saved", 4000)
                return
            files = d.run_files() or (
                self._primer_candidate_files() if d.uses_selected() else [])
            if not files:
                QMessageBox.information(
                    self, "Map Primers",
                    "Choose the sequences: click a FASTA in the Files panel (or fill "
                    "Sequences with Define inputs), then Run.",
                )
                return
            missing = [p for p in files if not Path(p).is_file()]
            if missing:
                QMessageBox.warning(self, "Map Primers",
                                    f"The file does not exist:\n{missing[0]}")
                return
            self._start_primer_job(PRIMER_JOB_MAP, files, self._primer_options)
            _forget_selection()   # the job already has its own copy of the settings

        self._present_te_viz_dialog(dialog, _accepted)

    def _quick_map_primers_right_click(self) -> None:
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

        def _path_changed(new_path: str) -> None:
            self._save_primer_options(replace(self._primer_options, database_path=new_path))
            self._set_status_message(
                f"Primer Database: {new_path or 'the database on this computer'}", 6000)

        dialog.databasePathChanged.connect(_path_changed)
        dialog.checkSequenceRequested.connect(
            lambda p, names: self._open_primer_check(str(p or ""), list(names)))

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

    def _start_painting_job(self, opts: PrimerOptions, *, open_dialog_on_error: bool = False) -> None:
        """Chromosome painting: run Chorus2 here or import its BED, then select."""
        from repeatviewer.backends.primers import painting as pt

        ps = opts.painting_settings()
        if ps.source == pt.SOURCE_IMPORT:
            anchor = ps.import_bed
            problem = ("Choose the Chorus2 BED to import (Design, Chorus2 tab)." if not anchor
                       else "" if Path(anchor).expanduser().is_file()
                       else f"The Chorus2 BED does not exist:\n{anchor}")
        else:
            anchor = ps.genome_fasta
            problem = ("Choose the genome FASTA (Design, Chorus2 tab)." if not anchor
                       else "" if Path(anchor).expanduser().is_file()
                       else f"The genome FASTA does not exist:\n{anchor}")
            if not problem:
                inst = pt.find_chorus2(ps.conda_env)
                if not inst.ok:
                    problem = (inst.message + "\n\n" + pt.install_help()
                               + "\n\nOr run it on a server and choose \u201cImport results "
                               "from the server\u201d.")
        if problem:
            if open_dialog_on_error:
                self._open_primer_design_dialog()
            QMessageBox.information(self, "Chromosome painting", problem)
            return
        self._start_primer_job(PRIMER_JOB_PAINTING, [Path(anchor).expanduser()], opts)

    def _confirm_primer_report_path(self, kind: str, files: list[Path]) -> Path | None | bool:
        """Rule 9: ask before an existing report is replaced. False = cancelled."""
        from repeatviewer.backends.primers import jobs as pj
        from repeatviewer.backends.primers import painting as pt
        from repeatviewer.widgets.file_conflict import confirm_existing_file_path

        tag = {PRIMER_JOB_MAP: pj.MAP_TAG, PRIMER_JOB_PAINTING: pt.PAINTING_TAG}.get(
            kind, pj.DESIGN_TAG)
        proposed = pj.default_report_path(files[0], tag)
        title = {PRIMER_JOB_MAP: "Map Primers", PRIMER_JOB_PAINTING: "Chromosome painting"}.get(
            kind, "Primer Design")
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
        if kind == PRIMER_JOB_PAINTING:
            settings = opts.painting_settings()
            title = f"Chromosome painting (Chorus2) on {files[0].name}"
        elif kind == PRIMER_JOB_MAP:
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
        if result.get("kind") == PRIMER_JOB_PAINTING:
            msg = (f"Chromosome painting finished: {result.get('selected', 0):,} probes selected "
                   f"from {result.get('chorus_input', 0):,} ({report.name})")
        elif result.get("kind") == PRIMER_JOB_MAP:
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
        if result.get("ask_pairs") and result.get("pair_list"):
            self._offer_designed_pairs(result["pair_list"], result.get("inputs") or [])

    def _offer_designed_pairs(self, pairs: list, inputs: list) -> None:
        from repeatviewer.views.primer_edit_dialog import DesignedPairsDialog

        dialog = DesignedPairsDialog(pairs, [Path(p) for p in inputs], parent=self)
        dialog.setWindowFlag(Qt.WindowType.Tool, True)

        def _done(code: int) -> None:
            if code:
                try:
                    changed = dialog.add_ticked()
                except Exception as exc:  # noqa: BLE001
                    QMessageBox.warning(self, "Primer Design", f"Could not add the primers:\n{exc}")
                    changed = []
                if changed:
                    self._reload_annotated_previews(changed, include_current=True)
                    self._set_status_message(
                        f"Added the ticked primers to {', '.join(p.name for p in changed)}", 8000)
            dialog.deleteLater()

        dialog.finished.connect(_done)
        dialog.show()
        dialog.raise_()

    def _reload_annotated_previews(self, fastas: list[Path], *, include_current: bool = False) -> None:
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
            if (page is current and not include_current) or not isinstance(page, FilePreviewPage):
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
