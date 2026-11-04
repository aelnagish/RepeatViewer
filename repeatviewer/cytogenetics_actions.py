"""Tools > Cytogenetics in the main window."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import QMessageBox

from repeatviewer.backends.cytogenetics import workspace as ws

CYTO_TITLE = "Cytogenetics"


class CytogeneticsActionsMixin:
    """Tools > Cytogenetics for ``MainWindow``."""

    def _open_cytogenetics_workspace(self, *, tab: str = "", sequence: Path | None = None) -> None:
        """Open the workspace of the current folder in the Preview area."""
        folder = getattr(self, "_current_folder", None)
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
