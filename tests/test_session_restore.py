"""Session restore settings round-trip."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QSettings
from PySide6.QtWidgets import QApplication

from repeatviewer.backends.sequence_features import SequenceMapData
from repeatviewer.database_settings import SETTINGS_GROUP
from repeatviewer.session_restore import (
    load_annotation_views,
    load_folder_file_memory,
    lookup_annotation_view,
    save_annotation_view,
    save_folder_file_memory,
)


def test_folder_file_memory_round_trip(tmp_path: Path) -> None:
    folder = tmp_path / "dir"
    folder.mkdir()
    file_path = folder / "seq.fa"
    file_path.write_text(">a\nACGT\n")

    s = QSettings()
    s.beginGroup(SETTINGS_GROUP)
    s.remove("folder_last_files_json")
    s.remove("folder_last_scroll_json")
    s.endGroup()

    key = str(folder.resolve())
    save_folder_file_memory(
        {key: str(file_path.resolve())},
        {key: 42},
    )
    files, scrolls = load_folder_file_memory()
    assert files[key].endswith("seq.fa")
    assert scrolls[key] == 42


def test_annotation_view_round_trip(tmp_path: Path) -> None:
    fasta = tmp_path / "x.fa"
    fasta.write_text(">q\nACGT\n")

    s = QSettings()
    s.beginGroup(SETTINGS_GROUP)
    s.remove("annotation_views_by_path_json")
    s.endGroup()

    state = {
        "start": 100,
        "end": 500,
        "has_selection": True,
        "ppb": 2.5,
        "scroll_px": 80.0,
        "sequence_id": "q",
    }
    save_annotation_view(fasta, state, sequence_id="q")
    hit = lookup_annotation_view(fasta, sequence_id="q")
    assert hit is not None
    assert hit["start"] == 100
    assert hit["end"] == 500
    assert hit["has_selection"] is True
    views = load_annotation_views()
    key = None
    for candidate in views:
        if candidate.endswith("\0q") or candidate.endswith("x.fa"):
            key = candidate
            break
    assert key is not None
    assert views[key]["ppb"] == 2.5


def test_annotation_view_restore_keeps_viewport_not_selection(qapp: QApplication) -> None:
    """Reopen must restore zoom/scroll, not jump to a prior selection."""
    from repeatviewer.widgets.sequence_map_viewer import SequenceMapWidget

    seq = "ACGT" * 1500  # 6000 bp
    data = SequenceMapData(
        fasta_path=Path("/tmp/rv_view_restore.fa"),
        sequence_id="s1",
        sequence=seq,
        features=[],
        is_protein=False,
        record_ids=["s1"],
        features_by_id={"s1": []},
        lengths_by_id={"s1": len(seq)},
    )
    widget = SequenceMapWidget()
    widget.resize(900, 400)
    widget.show()
    qapp.processEvents()
    widget.load(data, fit_after_load=False)
    qapp.processEvents()

    # Select a small early region, then leave for a different zoom/scroll.
    widget._canvas.set_selection_range(100, 200)
    widget._sync_position_spins(100, 200)
    target_ppb = 8.0
    widget._apply_zoom_ppb(target_ppb, sync_slider=True, reset_scroll=False)
    widget._update_horizontal_scroll()
    target_scroll = min(1200.0, float(widget._hbar.maximum()))
    widget._hbar.setValue(int(target_scroll))
    widget._canvas.set_scroll_px(target_scroll)
    qapp.processEvents()

    state = widget.view_state_for_session()
    assert state is not None
    assert state["has_selection"] is True
    assert state["start"] == 100
    assert state["end"] == 200
    assert abs(float(state["ppb"]) - target_ppb) < 0.01
    assert abs(float(state["scroll_px"]) - target_scroll) < 2.0

    # Simulate leaving and returning: clear view, then restore snapshot.
    widget._canvas.clear_selection()
    widget.fit_window()
    qapp.processEvents()
    widget.restore_view_state_for_session(state)
    # Drain deferred restore retries (layout may not be ready on first call).
    import time

    for _ in range(45):
        qapp.processEvents()
        if (
            abs(widget._canvas.pixels_per_bp() - target_ppb) < 0.05
            and abs(widget._canvas.scroll_px() - target_scroll) < 3.0
        ):
            break
        time.sleep(0.05)

    assert abs(widget._canvas.pixels_per_bp() - target_ppb) < 0.05
    assert abs(widget._canvas.scroll_px() - target_scroll) < 3.0
    assert widget.selection_range() == (100, 200)
    widget.close()
