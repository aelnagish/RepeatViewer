"""Define-database Files-panel path fill."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import QApplication

from repeatviewer.project import RepeatViewerProject
from repeatviewer.views.define_database_dialog import (
    DefineDatabaseDialog,
    resolve_panel_database_path,
)


def test_resolve_panel_database_path_fasta_and_hmm(tmp_path: Path) -> None:
    fasta = tmp_path / "Pfam.fa"
    fasta.write_text(">a\nACGT\n", encoding="utf-8")
    hmm = tmp_path / "Dfam.hmm"
    hmm.write_text("HMMER3/f\n", encoding="utf-8")
    cdd = tmp_path / "Cdd.pn"
    cdd.write_text("x\n", encoding="utf-8")

    assert resolve_panel_database_path(fasta) == ("user", fasta.resolve())
    assert resolve_panel_database_path(hmm) == ("user", hmm.resolve())
    assert resolve_panel_database_path(cdd) == ("cdd", tmp_path.resolve())


def test_apply_file_from_panel_fills_empty_user_row_when_cdd_focused(
    qapp: QApplication, tmp_path: Path
) -> None:
    del qapp
    fasta = tmp_path / "REXdb.fa"
    fasta.write_text(">a\nACDE\n", encoding="utf-8")
    dialog = DefineDatabaseDialog(RepeatViewerProject(), None)
    dialog.show()
    # Default focus target is CDD; a user DB file must still fill Other databases.
    assert dialog.apply_file_from_panel(fasta) is True
    rows = dialog._user_databases._rows
    assert any(row._path.text() == str(fasta.resolve()) for row in rows)
    dialog.close()


def test_apply_folder_from_panel_without_live_focus(
    qapp: QApplication, tmp_path: Path
) -> None:
    del qapp
    dialog = DefineDatabaseDialog(RepeatViewerProject(), None)
    dialog.show()
    assert dialog._focused_path_edit is dialog._cdd_path
    assert dialog.apply_folder_from_panel(tmp_path) is True
    assert dialog._cdd_path.text() == str(tmp_path.resolve())
    dialog.close()
