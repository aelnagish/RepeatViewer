"""Tools > Download from NCBI in the main window.

Kept out of ``main_window.py``: :class:`NcbiDownloadActionsMixin` is mixed
into ``MainWindow`` and only uses the window's shared job / dialog helpers
(``_present_te_viz_dialog``, ``_start_background_job``,
``_attach_process_monitor``, ``_populate_file_list`` ...).
"""

from __future__ import annotations

import threading
from dataclasses import replace
from pathlib import Path

from PySide6.QtCore import QObject, Signal, Slot
from PySide6.QtGui import QGuiApplication
from PySide6.QtWidgets import QMessageBox

from repeatviewer.backends import ncbi_download as nd
from repeatviewer.backends.job_cancel import CANCELLED_MESSAGE
from repeatviewer.job_options import (
    NCBI_RIGHT_CLICK_CLIPBOARD,
    NcbiDownloadOptions,
    load_ncbi_download_options,
    save_ncbi_download_options,
)

NCBI_TOOL_KEY = "ncbi_download"
NCBI_TITLE = "Download from NCBI"
# Above this many existing files, ask once instead of once per file.
_MAX_SINGLE_CONFLICT_QUESTIONS = 3


class NcbiDownloadWorker(QObject):
    """Download NCBI records off the UI thread."""

    progress = Signal(int, str)
    log_line = Signal(str)
    finished = Signal(object)
    failed = Signal(str)

    def __init__(self, *, database: str, records: list, fmt: str, one_file: bool,
                 file_name: str, out_dir: Path, targets: dict, email: str = "",
                 api_key: str = "", cancel=None) -> None:
        super().__init__()
        self._database = database
        self._records = list(records)
        self._fmt = fmt
        self._one_file = one_file
        self._file_name = file_name
        self._out_dir = Path(out_dir)
        self._targets = dict(targets)
        self._email = email
        self._api_key = api_key
        self._cancel = cancel

    def _cancelled(self) -> bool:
        return bool(self._cancel is not None and getattr(self._cancel, "requested", False))

    @Slot()
    def run(self) -> None:
        client = nd.NcbiClient(email=self._email, api_key=self._api_key, timeout=120.0)
        last = {"detail": ""}

        def _progress(done: int, total: int, detail: str) -> None:
            self.progress.emit(int(done), f"Downloading {detail}")
            if detail != last["detail"] and "downloaded" not in detail:
                self.log_line.emit(detail)
            last["detail"] = detail

        try:
            self.log_line.emit(
                f"{nd.DATABASE_LABELS.get(self._database, self._database)}: "
                f"{len(self._records)} record(s) → {self._out_dir}"
            )
            result = nd.download_records(
                client, self._database, self._records, self._out_dir,
                fmt=self._fmt, one_file=self._one_file, file_name=self._file_name,
                targets=self._targets, progress=_progress, cancel=self._cancelled,
            )
            features = 0
            for fasta in result.fasta_with_features:
                gff_text = result.feature_tables.get(str(fasta), "")
                features += self._import_features(fasta, gff_text, result.warnings)
            for w in result.warnings:
                self.log_line.emit(f"Warning: {w}")
            self.finished.emit({
                "files": [str(p) for p in result.files],
                "records": len(self._records),
                "features": features,
                "warnings": list(result.warnings),
            })
        except nd.NcbiCancelled:
            self.failed.emit(CANCELLED_MESSAGE)
        except nd.NcbiError as exc:
            self.failed.emit(str(exc))
        except Exception as exc:  # noqa: BLE001 - surface to the UI
            self.failed.emit(f"{type(exc).__name__}: {exc}")

    def _import_features(self, fasta: Path, gff_text: str, warnings: list[str]) -> int:
        """NCBI features → the FASTA's annotation sidecar (shown on the Annotation track)."""
        if gff_text.count("\n") <= 1:
            warnings.append(f"{fasta.name}: NCBI has no annotated features for these records.")
            return 0
        from repeatviewer.backends.sequence_features import import_annotation_as_rv_gff3

        tmp = fasta.with_name(f".{fasta.name}.ncbi_features.gff3")
        try:
            tmp.write_text(gff_text, encoding="utf-8")
            result = import_annotation_as_rv_gff3(fasta, tmp, merge_existing=False)
            count = int(getattr(result, "feature_count", 0))
            self.log_line.emit(f"{fasta.name}: {count} NCBI feature(s) added to the annotation")
            return count
        except Exception as exc:  # noqa: BLE001 - the sequences are still fine
            warnings.append(f"{fasta.name}: could not add the NCBI annotations ({exc}).")
            return 0
        finally:
            tmp.unlink(missing_ok=True)


class _LookupBridge(QObject):
    """Right-click: accession lookup result back on the GUI thread."""

    done = Signal(object)
    failed = Signal(str)


class NcbiDownloadActionsMixin:
    """Tools > Download from NCBI for ``MainWindow``."""

    _ncbi_options: NcbiDownloadOptions

    def _init_ncbi_download_actions(self) -> None:
        self._ncbi_options = load_ncbi_download_options()
        self._ncbi_lookup_bridge = _LookupBridge(self)
        self._ncbi_lookup_bridge.done.connect(self._on_ncbi_lookup_done)
        self._ncbi_lookup_bridge.failed.connect(self._on_ncbi_lookup_failed)

    def _save_ncbi_options(self, opts: NcbiDownloadOptions) -> None:
        self._ncbi_options = opts
        save_ncbi_download_options(opts)

    def _ncbi_default_folder(self) -> Path | None:
        folder = getattr(self, "_current_folder", None)
        return Path(folder) if folder else None

    # ------------------------------------------------------------- window

    def _open_ncbi_download_dialog(self, query_text: str = "") -> None:
        from repeatviewer.views.ncbi_download_dialog import NcbiDownloadDialog

        dialog = NcbiDownloadDialog(
            initial=self._ncbi_options,
            folder=self._ncbi_default_folder(),
            query_text=query_text,
            parent=self,
        )

        def _accepted(d) -> None:
            self._save_ncbi_options(d.options())
            if not d.should_run():
                self._set_status_message("Download from NCBI settings saved", 4000)
                return
            self._start_ncbi_download(
                database=d.database(),
                records=d.chosen_records(),
                fmt=d.output_format(),
                one_file=d.one_file(),
                file_name=d.file_name(),
                out_dir=d.out_dir(),
                query_text=d.query_text(),
            )

        self._present_te_viz_dialog(dialog, _accepted)

    def _remember_ncbi_dialog_tab(self, dialog) -> bool:
        from repeatviewer.views.ncbi_download_dialog import NcbiDownloadDialog

        if not isinstance(dialog, NcbiDownloadDialog):
            return False
        self._save_ncbi_options(replace(self._ncbi_options,
                                        last_dialog_tab=dialog.current_dialog_tab()))
        return True

    # ---------------------------------------------------------- right-click

    def _quick_ncbi_download_right_click(self) -> None:
        text = QGuiApplication.clipboard().text() or ""
        if self._ncbi_options.right_click_action != NCBI_RIGHT_CLICK_CLIPBOARD:
            self._open_ncbi_download_dialog()
            return
        accessions = nd.clipboard_accessions(text)
        if not accessions:
            QMessageBox.information(
                self, NCBI_TITLE,
                "The clipboard has no NCBI accessions. Copy one or more accessions "
                "(e.g. NC_045512.2) and right-click again.",
            )
            return
        if self._ncbi_default_folder() is None:
            self._open_ncbi_download_dialog(" ".join(accessions))
            return
        database = nd.accession_database(accessions[0]) or nd.DB_NUCLEOTIDE
        same = [a for a in accessions if nd.accession_database(a) == database]
        skipped = [a for a in accessions if a not in same]
        opts = self._ncbi_options
        query = nd.NcbiQuery(database=database, text=" ".join(same))
        client = nd.NcbiClient(email=opts.email, api_key=opts.api_key)
        bridge = self._ncbi_lookup_bridge
        self._set_status_message(f"Looking up {len(same)} accession(s) on NCBI…", 0)

        def _work() -> None:
            try:
                page = nd.search(client, query, retmax=max(len(same), 1))
                bridge.done.emit({"database": database, "records": page.records,
                                  "wanted": same, "skipped": skipped})
            except Exception as exc:  # noqa: BLE001
                bridge.failed.emit(str(exc))

        threading.Thread(target=_work, name="ncbi-lookup", daemon=True).start()

    @Slot(object)
    def _on_ncbi_lookup_done(self, payload) -> None:
        records = list(payload.get("records") or [])
        if not records:
            self._set_status_message("Ready", 1000)
            QMessageBox.information(self, NCBI_TITLE, "NCBI found none of the copied accessions.")
            return
        from repeatviewer.views.ncbi_download_dialog import missing_accessions

        missing = missing_accessions(list(payload.get("wanted") or []), records)
        skipped = list(payload.get("skipped") or [])
        if missing or skipped:
            notes = []
            if missing:
                notes.append("Not found: " + ", ".join(missing[:10]))
            if skipped:
                notes.append("Other database, not downloaded: " + ", ".join(skipped[:10]))
            self._set_status_message("; ".join(notes), 10000)
        database = payload["database"]
        opts = self._ncbi_options
        fmt = opts.format_for(database)
        wanted_text = " ".join(payload.get("wanted") or [])
        self._start_ncbi_download(
            database=database,
            records=records,
            fmt=fmt,
            one_file=opts.one_file,
            file_name=nd.default_file_name(database, fmt, records, wanted_text),
            out_dir=self._ncbi_default_folder(),
            query_text=wanted_text,
        )

    @Slot(str)
    def _on_ncbi_lookup_failed(self, message: str) -> None:
        self._set_status_message("Ready", 1000)
        QMessageBox.warning(self, NCBI_TITLE, message or "NCBI lookup failed.")

    # ---------------------------------------------------------------- job

    def _confirm_ncbi_outputs(self, planned: list[Path]) -> dict[str, Path] | None:
        """Rule 9 for every planned file. Returns planned → chosen path, None = cancelled.

        A planned path mapped to ``None`` means "skip this one".
        """
        from repeatviewer.widgets.file_conflict import confirm_existing_file_path

        existing = [p for p in planned if p.exists()]
        chosen: dict[str, Path] = {}
        if len(existing) > _MAX_SINGLE_CONFLICT_QUESTIONS:
            box = QMessageBox(self)
            box.setWindowTitle(NCBI_TITLE)
            box.setIcon(QMessageBox.Icon.Warning)
            box.setText(f"{len(existing)} of these files already exist in this folder.")
            box.setInformativeText("Overwrite them, or skip the records already downloaded?")
            overwrite = box.addButton("Overwrite", QMessageBox.ButtonRole.DestructiveRole)
            skip = box.addButton("Skip Existing", QMessageBox.ButtonRole.AcceptRole)
            box.addButton(QMessageBox.StandardButton.Cancel)
            box.exec()
            clicked = box.clickedButton()
            if clicked is overwrite:
                return {str(p): p for p in planned}
            if clicked is skip:
                return {str(p): (None if p.exists() else p) for p in planned}
            return None
        for path in planned:
            if path.exists():
                picked = confirm_existing_file_path(self, path, title=NCBI_TITLE,
                                                    required_suffix=path.suffix)
                if picked is None:
                    return None
                chosen[str(path)] = picked
            else:
                chosen[str(path)] = path
        return chosen

    def _start_ncbi_download(self, *, database: str, records: list, fmt: str, one_file: bool,
                             file_name: str, out_dir: Path | None, query_text: str) -> None:
        if not records:
            return
        if out_dir is None or not Path(out_dir).is_dir():
            QMessageBox.information(self, NCBI_TITLE,
                                    "Open the folder to save to in the Folders panel first.")
            return
        out_dir = Path(out_dir)
        if one_file and len(records) > 1 and not file_name:
            file_name = nd.default_file_name(database, fmt, records, query_text)
        planned = nd.planned_outputs(database, fmt, records, one_file=one_file,
                                     file_name=file_name, out_dir=out_dir)
        targets = self._confirm_ncbi_outputs(planned)
        if targets is None:
            return
        # "Skip Existing" drops records whose file is already there.
        if any(v is None for v in targets.values()):
            keep = []
            for rec in records:
                rec_paths = nd.planned_outputs(database, fmt, [rec], one_file=False,
                                               file_name="", out_dir=out_dir)
                if all(targets.get(str(p)) is not None for p in rec_paths):
                    keep.append(rec)
            records = keep
            targets = {k: v for k, v in targets.items() if v is not None}
            if not records:
                self._set_status_message("Nothing to download: every file already exists", 6000)
                return
        opts = self._ncbi_options
        noun = nd.DATABASE_LABELS.get(database, database)
        title = (f"{NCBI_TITLE}: {records[0].accession}" if len(records) == 1
                 else f"{NCBI_TITLE}: {len(records)} {noun} records")
        worker = NcbiDownloadWorker(
            database=database, records=records, fmt=fmt, one_file=one_file,
            file_name=file_name, out_dir=out_dir, targets=targets,
            email=opts.email, api_key=opts.api_key, cancel=self._job_cancel,
        )
        self._attach_process_monitor(
            worker,
            title=title,
            enabled=bool(opts.show_process_window),
            tool_key=NCBI_TOOL_KEY,
            mirror_progress=False,
            close_when_finished=bool(opts.close_process_window_when_finished),
        )
        self._start_background_job(
            worker,
            job_title=title,
            # One slot per record plus the last one, kept until the files show.
            total_steps=len(records) + 1,
            source_files=[],
            on_finished=self._on_ncbi_download_finished,
            on_failed=self._on_ncbi_download_failed,
        )

    @Slot(object)
    def _on_ncbi_download_finished(self, result) -> None:
        created = [Path(p).resolve() for p in result.get("files", []) if Path(p).is_file()]
        if created:
            self._refresh_folders_after_job(created, [])
        sequences = [p for p in created if p.suffix.lower() not in (".gff3", ".gff")]
        if sequences:
            self._auto_open_outputs(sequences[:1])
        n = len(created)
        msg = f"Download from NCBI finished: {n} file{'s' if n != 1 else ''}"
        if result.get("features"):
            msg += f", {result['features']:,} NCBI features added to the annotation"
        if result.get("warnings"):
            msg += f" ({len(result['warnings'])} warning(s), see the process window)"
        self._set_status_message(msg, 12000, override_file_link=True)

    @Slot(str)
    def _on_ncbi_download_failed(self, error: str) -> None:
        if error.strip() == CANCELLED_MESSAGE:
            return
        QMessageBox.warning(self, NCBI_TITLE, error.strip() or "Unknown error.")
