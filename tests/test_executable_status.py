"""Options → Executables: PATH status lines and named tools (FastGA / HISAT2)."""

from __future__ import annotations

from PySide6.QtWidgets import QApplication, QLabel

from repeatviewer.job_options import AlignOptions, BlastOptions, MapOptions
from repeatviewer.project import RepeatViewerProject
from repeatviewer.views.tool_dialogs import (
    AlignParamsDialog,
    BlastParamsDialog,
    set_named_executable_status,
    set_path_checkbox_status,
)


def test_set_path_checkbox_status_checked_found(qapp: QApplication) -> None:
    del qapp
    label = QLabel()
    label.resize(400, 20)
    set_path_checkbox_status(label, "/usr/bin/mafft", checked=True, has_builtin=True)
    assert "mafft" in label.text()
    assert "mid" in label.styleSheet()
    assert label.toolTip() == "/usr/bin/mafft"


def test_set_path_checkbox_status_checked_missing_builtin(qapp: QApplication) -> None:
    del qapp
    label = QLabel()
    set_path_checkbox_status(label, None, checked=True, has_builtin=True)
    assert "built-in" in label.text()
    assert "mid" in label.styleSheet()
    assert label.isVisible()


def test_set_path_checkbox_status_unchecked_clears(qapp: QApplication) -> None:
    del qapp
    label = QLabel()
    set_path_checkbox_status(label, "/usr/bin/mafft", checked=True)
    set_path_checkbox_status(label, "/usr/bin/mafft", checked=False)
    assert label.text() == ""
    assert label.isHidden()


def test_set_path_checkbox_status_elides_long_path(qapp: QApplication) -> None:
    del qapp
    label = QLabel()
    label.resize(120, 20)
    long_path = "/very/long/prefix/" + ("dir/" * 20) + "mafft"
    set_path_checkbox_status(label, long_path, checked=True)
    assert label.toolTip() == long_path
    assert len(label.text()) < len(long_path)
    assert "…" in label.text() or "..." in label.text() or label.text() != long_path


def test_set_named_executable_status(qapp: QApplication) -> None:
    del qapp
    label = QLabel()
    label.resize(400, 20)
    set_named_executable_status(label, "FastGA", "/opt/bin/FastGA")
    assert "<b>FastGA</b>" in label.text()
    assert "FastGA" in label.text()
    assert "mid" in label.styleSheet()
    set_named_executable_status(label, "HISAT2", None, missing_text="not configured")
    assert "<b>HISAT2</b>" in label.text()
    assert "not configured" in label.text()
    assert "mid" in label.styleSheet()


def test_align_options_shows_fastga_status(qapp: QApplication, tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(
        "repeatviewer.views.tool_dialogs.AlignParamsDialog._remind_fastga_if_needed",
        lambda self: None,
    )
    monkeypatch.setattr(
        "repeatviewer.backends.fastga_dotplot.fastga_executable",
        lambda: None,
    )
    dlg = AlignParamsDialog(
        parent=None,
        initial=AlignOptions(),
        fastga_files=[tmp_path / "a.fa"],
    )
    dlg.show()
    dlg._aligner_tabs.setCurrentWidget(dlg._options_tab)
    qapp.processEvents()
    assert "<b>FastGA</b>" in dlg._fastga_path_status.text()
    assert "not configured" in dlg._fastga_path_status.text()

    dlg._use_system_mafft.setChecked(True)
    qapp.processEvents()
    # Either a real path or the grey built-in fallback.
    status = dlg._mafft_path_status.text()
    assert status
    assert status.startswith("/") or "built-in" in status


def test_blast_options_path_status(qapp: QApplication, monkeypatch) -> None:
    monkeypatch.setattr(
        "repeatviewer.views.tool_dialogs.resolve_system_executable",
        lambda *names: None,
    )
    dlg = BlastParamsDialog(
        parent=None,
        project=RepeatViewerProject(),
        initial=BlastOptions(),
    )
    dlg.show()
    dlg._tabs.setCurrentIndex(dlg._options_tab_index)
    qapp.processEvents()
    dlg._use_system_blast.setChecked(True)
    qapp.processEvents()
    assert "built-in" in dlg._blast_path_status.text()
    dlg._use_system_blast.setChecked(False)
    qapp.processEvents()
    assert dlg._blast_path_status.text() == ""


def test_mapping_options_hisat2_status(qapp: QApplication, tmp_path, monkeypatch) -> None:
    from repeatviewer.views.map_reference_dialog import MapReferenceDialog

    monkeypatch.setattr(
        "repeatviewer.views.map_reference_dialog.resolve_system_executable",
        lambda *names: (
            "/fake/hisat2" if "hisat2" in names else None
        ),
    )
    q = tmp_path / "q.fa"
    q.write_text(">s\nACGT\n")
    dlg = MapReferenceDialog(files=[q], initial=MapOptions())
    dlg.show()
    dlg._tab_bar.setCurrentIndex(dlg._TAB_OPTIONS)
    qapp.processEvents()
    assert "<b>HISAT2</b>" in dlg._hisat2_path_status.text()
    assert "/fake/hisat2" in dlg._hisat2_path_status.text()

    dlg._use_system_path.setChecked(True)
    qapp.processEvents()
    assert "built-in" in dlg._minimap2_path_status.text()
