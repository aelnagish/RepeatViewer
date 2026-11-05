"""Tools > Cytogenetics in the main window."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import QMessageBox

from repeatviewer.backends.cytogenetics import workspace as ws

CYTO_TITLE = "Cytogenetics"


class CytogeneticsActionsMixin:
    """Tools > Cytogenetics for ``MainWindow``."""

    def _open_cytogenetics_workspace(self, *, tab: str = "", sequence: Path | None = None, image_id: str = ""):
        """Open the workspace in the Preview area: the sequence's folder, else the current folder.

        Returns the workspace panel (or None). With ``image_id`` that image opens in the Image tab.
        """
        folder = Path(sequence).expanduser().resolve().parent if sequence is not None else getattr(self, "_current_folder", None)
        if (not folder or not Path(folder).is_dir()) and sequence is None:
            from PySide6.QtCore import QSettings

            last = str(QSettings().value("cytogenetics_workspace/last_folder", "") or "")
            if last and ws.workspace_path(Path(last)).exists():
                folder = Path(last)  # no folder chosen: back to the workspace used last time
        if not folder or not Path(folder).is_dir():
            QMessageBox.information(
                self, CYTO_TITLE, "Select a folder in the Folders panel first."
            )
            return
        try:
            marker = ws.ensure_workspace(Path(folder))
        except OSError as exc:
            QMessageBox.warning(self, CYTO_TITLE, f"Could not create the workspace:\n{exc}")
            return
        self._open_preview_tab(marker, replace=False)
        page = self._preview_tabs.currentWidget() if self._preview_tabs is not None else None
        panel = page.cyto_workspace_panel() if hasattr(page, "cyto_workspace_panel") else None
        if panel is None:
            return
        if sequence is not None:
            panel.set_sequence_filter(Path(sequence).name)
        if tab:
            panel.show_tab(tab)
        if image_id:
            panel.open_image(image_id)
        return panel

    def _cytogenetics_library_action(self, method: str) -> None:
        self._open_cytogenetics_workspace(tab="Library")
        page = self._preview_tabs.currentWidget() if self._preview_tabs is not None else None
        panel = page.cyto_workspace_panel() if hasattr(page, "cyto_workspace_panel") else None
        if panel is not None:
            getattr(panel.library(), method)()

    def _cytogenetics_import_session(self) -> None:
        self._cytogenetics_library_action("start_import")

    def _cytogenetics_edit_lists(self) -> None:
        self._cytogenetics_library_action("edit_lists")

    def _populate_cytogenetics_menu(self, menu) -> None:
        menu.addAction("Open workspace", lambda: self._open_cytogenetics_workspace())
        menu.addAction("Import", self._cytogenetics_import_session)
        menu.addAction("Add image...", lambda: self._cytogenetics_library_action("add_image"))
        menu.addAction("Edit lists...", self._cytogenetics_edit_lists)
        menu.addSeparator()
        menu.addAction("Make Chromosome Sizes File", lambda: self._make_chromosome_sizes())

    def _sync_cytogenetics_windows_from_files(self, path: Path) -> None:
        """A click in the Files or Folders panel fills the path field of an open Cytogenetics window."""
        from repeatviewer.widgets.cyto_panel_targets import offer

        offer(path)

    def _make_chromosome_sizes(self, path: Path | None = None) -> None:
        """Write ``<name>_sizes.txt`` (sequence name and length) for the selected FASTA files, in the background."""
        from repeatviewer.backends.cytogenetics import karyotype as kt
        from repeatviewer.backends.job_cancel import JobCancelledError
        from repeatviewer.main_window import WorkflowWorker
        from repeatviewer.widgets.file_conflict import confirm_existing_file_path

        title = "Make Chromosome Sizes File"
        files = [Path(path)] if path else [p for p in self._selected_files() if kt.is_fasta(p)]
        if not files:
            QMessageBox.information(self, title, "Select one or more FASTA files in the Files panel.")
            return
        jobs = []
        for fasta in files:
            out = confirm_existing_file_path(self, kt.sizes_path_for(fasta), title=title)
            if out is None:
                return
            jobs.append((fasta.expanduser().resolve(), out))

        def _run(progress_cb, *, jobs, cancel):
            created = []
            n = len(jobs)
            for i, (src, out) in enumerate(jobs):
                def step(pct, i=i, src=src):
                    progress_cb(i * 100 + int(pct), f"{src.name}: {i} done, {n - i} left")

                kt.make_sizes_file(src, out, step, lambda: cancel.requested)
                cancel.check()
                created.append(str(out))
                progress_cb((i + 1) * 100, f"{i + 1} done, {n - i - 1} left")
            if cancel.requested:
                raise JobCancelledError("Job stopped by user")
            return {"func_name": "chromosome_sizes", "created_outputs": created,
                    "source_files": [str(s) for s, _ in jobs]}

        worker = WorkflowWorker(_run, {"jobs": jobs, "cancel": self._job_cancel})
        self._start_background_job(
            worker,
            job_title="chromosome_sizes",
            total_steps=100 * len(jobs),
            source_files=[s for s, _ in jobs],
            on_finished=self._on_worker_finished,
            on_failed=self._on_worker_failed,
        )
