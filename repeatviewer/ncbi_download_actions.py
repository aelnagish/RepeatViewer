"""Tools > Download from NCBI in the main window.

Kept out of ``main_window.py``: :class:`NcbiDownloadActionsMixin` is mixed
into ``MainWindow``.

Downloads do NOT use the window's single job slot (``_start_background_job``):
a genome can take an hour, and that slot greys out every tool. Instead each
download runs on its own daemon thread with its own progress bar and Stop
button in the status bar, so browsing, previews and other tools keep working.
Further downloads wait in a queue and start one after another.

Developer: Abdullah Yusuf
"""

from __future__ import annotations

import re
import threading
import time
from dataclasses import replace
from pathlib import Path

from PySide6.QtCore import QObject, Signal, Slot
from PySide6.QtGui import QGuiApplication
from PySide6.QtWidgets import QMessageBox

from repeatviewer.backends import ncbi_download as nd
from repeatviewer.backends.job_cancel import CANCELLED_MESSAGE
from repeatviewer.job_options import (
    NcbiDownloadOptions,
    load_ncbi_download_options,
    save_ncbi_download_options,
)

NCBI_TOOL_KEY = "ncbi_download"
NCBI_TITLE = "Download from NCBI"
# Above this many existing files, ask once instead of once per file.
_MAX_SINGLE_CONFLICT_QUESTIONS = 3
# The genome FASTA is not opened in Preview above this size (it would stall the window).
_AUTO_OPEN_MAX_BYTES = 10 * 1024 * 1024
# At most this many progress updates per second reach the window.
_PROGRESS_INTERVAL = 0.25
_DETAIL_RE = re.compile(r"^(?P<acc>[^\s:]+):?\s*(?P<rest>.*?)\s*(\(\d+ done, \d+ left\))?$")


def split_progress_detail(detail: str) -> tuple[str, str]:
    """``"GCF_1: 161.0 MB downloaded (0 done, 3 left)"`` → ``("GCF_1", "161.0 MB")``."""
    m = _DETAIL_RE.match(str(detail or "").strip())
    if not m:
        return str(detail or ""), ""
    rest = m.group("rest") or ""
    rest = rest.replace(" downloaded", "")
    if rest.endswith("done, 0 left"):
        rest = ""
    return m.group("acc"), rest


class _CancelFlag:
    """Stop button of one download (independent of the main job's Stop)."""

    def __init__(self) -> None:
        self._event = threading.Event()

    def request(self) -> None:
        self._event.set()

    @property
    def requested(self) -> bool:
        return self._event.is_set()


class NcbiDownloadWorker(QObject):
    """Download NCBI records off the UI thread."""

    progress = Signal(int, int, str, str)
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
        last = {"detail": "", "done": -1, "at": 0.0}

        def _progress(done: int, total: int, detail: str) -> None:
            # Throttled: a fast connection reports every MB, the window needs a few per second.
            now = time.monotonic()
            if done == last["done"] and now - last["at"] < _PROGRESS_INTERVAL:
                return
            last["done"], last["at"] = done, now
            acc, extra = split_progress_detail(detail)
            self.progress.emit(int(done), int(total), acc, extra)
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
    """Right-click: accession lookup result back on the GUI thread.

    The slots live on this QObject (not on the mixin) so Qt queues them to the
    GUI thread; a mixin method has no receiver object and would run on the
    lookup thread.
    """

    done = Signal(object)
    failed = Signal(str)

    def __init__(self, window) -> None:
        super().__init__(window)
        self._window = window
        self.done.connect(self._on_done)
        self.failed.connect(self._on_failed)

    @Slot(object)
    def _on_done(self, payload) -> None:
        self._window._on_ncbi_lookup_done(payload)

    @Slot(str)
    def _on_failed(self, message: str) -> None:
        self._window._on_ncbi_lookup_failed(message)


class _JobRelay(QObject):
    """Delivers one download's signals to the window, on the GUI thread."""

    def __init__(self, window, worker: NcbiDownloadWorker) -> None:
        super().__init__(window)
        self._window = window
        self._worker = worker
        worker.progress.connect(self._on_progress)
        worker.finished.connect(self._on_finished)
        worker.failed.connect(self._on_failed)

    @Slot(int, int, str, str)
    def _on_progress(self, done: int, total: int, accession: str, extra: str) -> None:
        self._window._on_ncbi_progress(self._worker, done, total, accession, extra)

    @Slot(object)
    def _on_finished(self, result) -> None:
        self._window._on_ncbi_download_finished(self._worker, result)
        self.deleteLater()

    @Slot(str)
    def _on_failed(self, error: str) -> None:
        self._window._on_ncbi_download_failed(self._worker, error)
        self.deleteLater()


class NcbiDownloadActionsMixin:
    """Tools > Download from NCBI for ``MainWindow``."""

    _ncbi_options: NcbiDownloadOptions

    def _init_ncbi_download_actions(self) -> None:
        self._ncbi_options = load_ncbi_download_options()
        self._ncbi_worker: NcbiDownloadWorker | None = None
        self._ncbi_cancel: _CancelFlag | None = None
        self._ncbi_queue: list[dict] = []
        self._ncbi_progress = None
        self._ncbi_stop_btn = None
        self._ncbi_job_title = ""
        self._ncbi_lookup_bridge = _LookupBridge(self)

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
        """Right-click: download clipboard accessions with saved settings."""
        text = QGuiApplication.clipboard().text() or ""
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
        job = {"database": database, "records": records, "fmt": fmt, "one_file": one_file,
               "file_name": file_name, "out_dir": out_dir, "targets": targets}
        if self._ncbi_worker is not None:
            self._ncbi_queue.append(job)
            self._set_status_message(
                f"NCBI download queued; it starts when the current one finishes "
                f"({len(self._ncbi_queue)} waiting)", 6000)
            self._sync_ncbi_progress_label()
            return
        self._run_ncbi_job(job)

    def _run_ncbi_job(self, job: dict) -> None:
        """Start one download on its own thread; the window stays free."""
        opts = self._ncbi_options
        records = job["records"]
        database = job["database"]
        noun = nd.DATABASE_LABELS.get(database, database)
        title = (f"{NCBI_TITLE}: {records[0].accession}" if len(records) == 1
                 else f"{NCBI_TITLE}: {len(records)} {noun} records")
        cancel = _CancelFlag()
        worker = NcbiDownloadWorker(
            database=database, records=records, fmt=job["fmt"], one_file=job["one_file"],
            file_name=job["file_name"], out_dir=job["out_dir"], targets=job["targets"],
            email=opts.email, api_key=opts.api_key, cancel=cancel,
        )
        # The worker QObject stays on the GUI thread; signals emitted from the
        # download thread are queued to it automatically.
        _JobRelay(self, worker)
        if opts.show_process_window:
            self._attach_process_monitor(
                worker,
                title=title,
                enabled=True,
                tool_key=NCBI_TOOL_KEY,
                mirror_progress=False,
                close_when_finished=bool(opts.close_process_window_when_finished),
            )
        self._ncbi_worker = worker
        self._ncbi_cancel = cancel
        self._ncbi_job_title = title
        self._show_ncbi_progress(len(records))
        # Daemon thread: quitting RepeatViewer mid-download never waits for NCBI.
        threading.Thread(target=worker.run, name="ncbi-download", daemon=True).start()

    # ------------------------------------------------------------ progress

    def _ensure_ncbi_progress_widgets(self) -> None:
        if self._ncbi_progress is not None:
            return
        from PySide6.QtWidgets import QPushButton

        from repeatviewer.widgets.task_progress import TaskProgressPanel

        status = self.statusBar()
        panel = TaskProgressPanel()
        stop = QPushButton("Stop")
        stop.setToolTip("Stop the NCBI download")
        stop.setMaximumHeight(22)
        stop.clicked.connect(self._stop_ncbi_download)
        stop.hide()
        # Left of the main job's progress bar and Stop button (index 0 is the
        # status-text file link, a normal widget; permanent ones come after it).
        status.insertPermanentWidget(1, panel)
        status.insertPermanentWidget(2, stop)
        self._ncbi_progress = panel
        self._ncbi_stop_btn = stop

    def _show_ncbi_progress(self, total: int) -> None:
        self._ensure_ncbi_progress_widgets()
        self._ncbi_progress.start("Downloading", max(1, int(total)))
        self._ncbi_stop_btn.setEnabled(True)
        self._ncbi_stop_btn.show()
        self._sync_ncbi_progress_label()

    def _sync_ncbi_progress_label(self, accession: str = "", extra: str = "") -> None:
        panel = self._ncbi_progress
        if panel is None or self._ncbi_worker is None:
            return
        waiting = len(self._ncbi_queue)
        suffix = f"  {extra}" if extra else ""
        if waiting:
            suffix += f"  (+{waiting} queued)"
        panel.set_bar_suffix(suffix)
        if accession:
            panel.update(panel.completed(), accession)

    def _hide_ncbi_progress(self) -> None:
        if self._ncbi_progress is not None:
            self._ncbi_progress.reset()
        if self._ncbi_stop_btn is not None:
            self._ncbi_stop_btn.hide()

    def _on_ncbi_progress(self, worker, done: int, total: int, accession: str, extra: str) -> None:
        panel = self._ncbi_progress
        if panel is None or worker is not self._ncbi_worker:
            return
        panel.update(min(done, max(0, total - 1)) if total > 1 else 0, accession)
        self._sync_ncbi_progress_label(extra=extra)

    def _stop_ncbi_download(self) -> None:
        """Stop the running download and drop the queued ones."""
        self._ncbi_queue.clear()
        if self._ncbi_cancel is not None:
            self._ncbi_cancel.request()
        if self._ncbi_stop_btn is not None:
            self._ncbi_stop_btn.setEnabled(False)
        if self._ncbi_progress is not None:
            self._ncbi_progress.update(self._ncbi_progress.completed(), "stopping…")

    def _ncbi_job_done(self) -> None:
        self._ncbi_worker = None
        self._ncbi_cancel = None
        if self._ncbi_queue:
            self._run_ncbi_job(self._ncbi_queue.pop(0))
        else:
            self._hide_ncbi_progress()

    def _on_ncbi_download_finished(self, worker, result) -> None:
        if worker is not self._ncbi_worker:
            return
        self._ncbi_job_done()
        created = [Path(p).resolve() for p in result.get("files", []) if Path(p).is_file()]
        self._show_ncbi_outputs(created)
        n = len(created)
        msg = f"Download from NCBI finished: {n} file{'s' if n != 1 else ''}"
        if result.get("features"):
            msg += f", {result['features']:,} NCBI features added to the annotation"
        if result.get("warnings"):
            msg += f" ({len(result['warnings'])} warning(s), see the process window)"
        self._set_status_message(msg, 12000, override_file_link=True)

    def _show_ncbi_outputs(self, created: list[Path]) -> None:
        """New files: green in the Files panel; open the first one when it is small.

        While another tool's job runs, only the file list is refreshed so that
        job's own results and highlights are left alone.
        """
        if not created:
            return
        current = self._current_folder.resolve() if getattr(self, "_current_folder", None) else None
        if self._any_background_job_running():
            if current is not None and any(p.parent == current for p in created):
                self._populate_file_list(current, open_preview=False)
            return
        if current is not None and all(p.parent != current for p in created):
            # The user browsed elsewhere meanwhile: do not pull them back.
            self._refresh_subfolders_under(created[0].parent)
            return
        self._refresh_folders_after_job(created, [])
        sequences = [p for p in created if p.suffix.lower() not in (".gff3", ".gff")]
        small = [p for p in sequences if p.stat().st_size <= _AUTO_OPEN_MAX_BYTES]
        if small:
            self._auto_open_outputs(small[:1])

    def _on_ncbi_download_failed(self, worker, error: str) -> None:
        if worker is not self._ncbi_worker:
            return
        self._ncbi_job_done()
        if error.strip() == CANCELLED_MESSAGE:
            self._set_status_message("NCBI download stopped", 6000)
            return
        QMessageBox.warning(self, NCBI_TITLE, error.strip() or "Unknown error.")
