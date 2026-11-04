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
        menu.addAction("Import slide session...", self._cytogenetics_import_session)
        menu.addAction("Add image...", lambda: self._cytogenetics_library_action("add_image"))
        menu.addAction("Edit lists...", self._cytogenetics_edit_lists)
