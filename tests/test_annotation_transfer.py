"""Unit tests for alignment-based annotation transfer."""

from __future__ import annotations

from pathlib import Path

import pytest

from repeatviewer.backends.annotation_transfer import (
    AnnotationTransferOptions,
    AnnotationTransferResult,
    transfer_annotations,
)
from repeatviewer.backends.map_reference import (
    STATUS_INTACT,
    MappedRead,
    compute_stats,
    lift_features,
)
from repeatviewer.backends.sequence_features import SequenceFeature


def _perfect_read(*, ref_len: int = 100, name: str = "tgt") -> MappedRead:
    row = "A" * ref_len
    ref_row = "A" * ref_len
    read = MappedRead(
        name=name,
        ref_id="src",
        start=0,
        end=ref_len,
        strand="+",
        row=row,
        insertions=[],
        ref_row=ref_row,
        mapq=60,
        kind="primary",
        query_length=ref_len,
        clipped=0,
        clip_left=0,
    )
    return compute_stats(read)


def test_lift_features_intact_span() -> None:
    features = [
        SequenceFeature(start=11, end=20, name="LTR1", feature_type="LTR", strand="+"),
    ]
    lifted = lift_features(_perfect_read(), features, reference="A" * 100)
    assert len(lifted) == 1
    assert lifted[0].status == STATUS_INTACT
    assert lifted[0].start == 11
    assert lifted[0].end == 20


def test_lift_features_partial_when_alignment_clips_feature() -> None:
    features = [
        SequenceFeature(start=1, end=40, name="ORF1", feature_type="ORF", strand="+"),
    ]
    # Alignment only covers bases 10..29 on the reference.
    row = "A" * 20
    ref_row = "A" * 20
    read = MappedRead(
        name="tgt",
        ref_id="src",
        start=10,
        end=30,
        strand="+",
        row=row,
        insertions=[],
        ref_row=ref_row,
        mapq=60,
        kind="primary",
        query_length=20,
        clipped=0,
        clip_left=0,
    )
    read = compute_stats(read)
    lifted = lift_features(read, features, reference="A" * 100)
    assert len(lifted) == 1
    assert lifted[0].status == "partial"
    assert lifted[0].start > 0
    assert lifted[0].end >= lifted[0].start


def test_transfer_annotations_identical_sequences(tmp_path: Path) -> None:
    """End-to-end transfer when source and target are identical (needs minimap2)."""
    from repeatviewer.backends.map_reference import MAP_ENGINE_MINIMAP2, engine_executable

    if not engine_executable(MAP_ENGINE_MINIMAP2):
        pytest.skip("minimap2 not available")

    seq = ("ACGT" * 80) + ("TTAA" * 20)  # 400 bp
    features = [
        SequenceFeature(start=21, end=60, name="f1", feature_type="LTR", strand="+"),
        SequenceFeature(start=200, end=240, name="f2", feature_type="ORF", strand="-"),
    ]
    result = transfer_annotations(
        source_id="src",
        source_sequence=seq,
        target_id="tgt",
        target_sequence=seq,
        features=features,
        opts=AnnotationTransferOptions(min_identity=0.5, min_coverage=0.5),
    )
    assert isinstance(result, AnnotationTransferResult)
    assert result.features
    assert len(result.features) == 2
    counts = result.counts()
    assert counts["transferred"] == 2
    assert counts["lost"] == 0


def test_sequence_map_panel_has_geneious_rail() -> None:
    from PySide6.QtWidgets import QApplication
    import sys

    app = QApplication.instance() or QApplication(sys.argv)
    from repeatviewer.widgets.sequence_map_viewer import SequenceMapPanel
    from repeatviewer.widgets.sequence_annotation_rail import SequenceAnnotationRail

    fasta = Path("/tmp/rv_annot_rail_smoke.fa")
    fasta.write_text(">s1\nACGTACGTACGTACGT\n>s2\nACGTACGTACGTACGT\n", encoding="utf-8")
    panel = SequenceMapPanel(fasta, autoload=False)
    try:
        assert isinstance(panel._rail, SequenceAnnotationRail)
        assert panel._rail.add_feature_btn.text() == "Add feature"
        assert panel._rail._tab_buttons[0]._label == "Display"
        assert panel._rail._tab_buttons[1]._label == "Annotation"
        assert panel._rail._tab_buttons[2]._label == "Find"
        assert hasattr(panel._rail, "find_query")
        assert panel._rail._stack.count() == 3
        assert panel._auto_ann_run_btn.text() == "Run feature detection"
        assert panel._similarity_ann_run_btn.text() == "Run similarity"
        assert hasattr(panel._rail, "show_complement")
        assert hasattr(panel._rail, "transfer_run_btn")
        assert panel.layout().indexOf(panel._map) >= 0
        assert panel.layout().indexOf(panel._rail) >= 0
        assert panel.layout().indexOf(panel._map) < panel.layout().indexOf(panel._rail)
    finally:
        panel.close()


def test_visible_frame_residues_match_full_translate() -> None:
    """Viewport translation must match full-frame AA placement (all 6 frames)."""
    from repeatviewer.widgets.sequence_map_viewer import (
        _iter_visible_frame_residues,
        _translate_dna_frame,
    )

    seq = "ATGAAACCCGGGTTT" * 20
    vs, ve = 40, 80
    for frame in (1, 2, 3, -1, -2, -3):
        full = _translate_dna_frame(seq, frame)
        got = list(_iter_visible_frame_residues(seq, frame, vs, ve))
        expected: list[tuple[int, int, str]] = []
        if frame > 0:
            for ai, aa in enumerate(full):
                s0 = (frame - 1) + ai * 3
                e0 = s0 + 3
                if e0 <= vs or s0 >= ve:
                    continue
                expected.append((s0, e0, aa))
        else:
            offset = abs(frame) - 1
            n = len(seq)
            for ai, aa in enumerate(full):
                c0 = offset + ai * 3
                fs, fe = n - (c0 + 3), n - c0
                if fe <= vs or fs >= ve:
                    continue
                expected.append((fs, fe, aa))
        assert got == expected


def test_extra_rows_viewport_work_stays_fast_on_long_sequence() -> None:
    """Complement + 6 frames must not reverse-complement / translate the genome."""
    import time

    from repeatviewer.widgets.sequence_map_viewer import (
        _complement_base,
        _iter_visible_frame_residues,
    )

    # 5 Mbp — old path would RC + translate the whole string each paint.
    big = ("ACGTACGTACGTACGT" * 312_500)  # 5_000_000 bp
    assert len(big) == 5_000_000
    vs, ve = 2_500_000, 2_500_120
    t0 = time.perf_counter()
    for i in range(vs, ve):
        _complement_base(big[i])
    for frame in (1, 2, 3, -1, -2, -3):
        list(_iter_visible_frame_residues(big, frame, vs, ve))
    elapsed = time.perf_counter() - t0
    assert elapsed < 0.75, f"viewport extra-row work too slow: {elapsed:.3f}s"


def test_complement_and_translation_grow_canvas_height(tmp_path: Path) -> None:
    from PySide6.QtWidgets import QApplication
    import sys

    app = QApplication.instance() or QApplication(sys.argv)
    from repeatviewer.widgets.sequence_map_viewer import SequenceMapPanel

    fasta = tmp_path / "disp.fa"
    fasta.write_text(">s1\n" + ("ATGAAACCCGGGTTT" * 40) + "\n", encoding="utf-8")
    panel = SequenceMapPanel(fasta, autoload=True)
    try:
        panel.resize(900, 500)
        panel.show()
        app.processEvents()
        panel._map._apply_zoom_ppb(20.0, sync_slider=True)
        app.processEvents()
        canvas = panel._map._canvas
        assert canvas.shows_letters()
        # Settings may restore complement/translation — start from a clean baseline.
        canvas.set_show_complement(False)
        canvas.set_translation_frames(frozenset())
        panel._map._sync_header_views()
        app.processEvents()
        base_h = canvas._content_height()
        canvas.set_show_complement(True)
        panel._map._sync_header_views()
        app.processEvents()
        with_comp = canvas._content_height()
        assert with_comp > base_h
        assert canvas.height() == with_comp
        assert panel._map._canvas_row.height() == with_comp
        canvas.set_translation_frames({1, 2})
        with_tr = canvas._content_height()
        assert with_tr > with_comp

        # Zoomed out: complement / translation must not add rows.
        canvas.set_show_complement(True)
        canvas.set_translation_frames({1, 2, 3})
        panel._map._apply_zoom_ppb(0.5, sync_slider=True)
        panel._map._sync_header_views()
        app.processEvents()
        assert not canvas.shows_letters()
        assert canvas._extra_display_rows() == []
        assert canvas._extra_display_height() == 0.0
    finally:
        panel.close()


def test_show_complement_paints_second_nucleotide_row(tmp_path: Path) -> None:
    """Checking Show complement must draw a visible second base row (not just reserve height)."""
    from PySide6.QtWidgets import QApplication
    import sys

    app = QApplication.instance() or QApplication(sys.argv)
    from repeatviewer.widgets.sequence_map_viewer import (
        SequenceMapPanel,
        _complement_base,
    )

    fasta = tmp_path / "comp_paint.fa"
    fasta.write_text(">s1\nTGAAGGAAGG\n", encoding="utf-8")
    panel = SequenceMapPanel(fasta, autoload=True)
    try:
        panel.resize(900, 500)
        panel.show()
        app.processEvents()
        panel._rail.show_complement.setChecked(False)
        panel._rail.show_translation.setChecked(False)
        app.processEvents()
        panel._map._apply_zoom_ppb(24.0, sync_slider=True)
        app.processEvents()

        panel._rail.show_complement.setChecked(True)
        app.processEvents()
        panel._map._sync_header_views()
        app.processEvents()

        canvas = panel._map._canvas
        assert canvas._show_complement is True
        assert "complement" in canvas._extra_display_rows()
        assert canvas.height() == canvas._content_height()
        assert panel._map._canvas_row.height() >= canvas.height()

        canvas._invalidate_seq_layer_cache()
        canvas._ensure_seq_layer_cache()
        pm = canvas._seq_layer_cache
        assert pm is not None and not pm.isNull()
        dpr = float(pm.devicePixelRatio())
        img = pm.toImage()
        _, origin_y, cell_w, cell_h = canvas._nucleotide_cell_geometry()
        # Mid-forward (T → pink) and mid-complement (A → green).
        fx = int((cell_w / 2.0) * dpr)
        fy = int((origin_y + cell_h / 2.0) * dpr)
        # Forward ends at origin_y+cell_h; gap=2 then complement cell (see draw path).
        cy = int((origin_y + cell_h + 2.0 + cell_h / 2.0) * dpr)
        fwd = img.pixelColor(fx, fy)
        comp = img.pixelColor(fx, cy)
        assert fwd.red() > fwd.green() + 20, f"expected pink T, got {fwd.name()}"
        assert comp.green() > comp.red() + 20, f"expected green A under T, got {comp.name()}"
        assert _complement_base("T") == "A"
    finally:
        panel.close()


def test_annotation_map_keeps_paper_fill_and_hbar_layout(tmp_path: Path) -> None:
    """Empty space below tracks uses paper fill; zoom-in reserves hbar height."""
    from PySide6.QtWidgets import QApplication
    import sys

    app = QApplication.instance() or QApplication(sys.argv)
    from repeatviewer.widgets.sequence_map_viewer import (
        SequenceMapPanel,
        _ANNOT_PAPER,
        _MAP_PAPER,
    )

    fasta = tmp_path / "paper.fa"
    fasta.write_text(">s1\n" + ("ACGT" * 800) + "\n", encoding="utf-8")
    panel = SequenceMapPanel(fasta, autoload=True)
    try:
        panel.resize(900, 500)
        panel.show()
        app.processEvents()
        panel._map.fit_window()
        app.processEvents()
        block = panel._map._sequence_block
        # Annotation region stays light paper even in dark mode.
        assert block.palette().color(block.backgroundRole()).name() in {
            _ANNOT_PAPER.name(),
            _MAP_PAPER.name(),
        }
        # Zoom in far enough that horizontal scrolling is required.
        panel._map._apply_zoom_ppb(8.0, sync_slider=True)
        app.processEvents()
        panel._map._sync_header_views()
        app.processEvents()
        assert panel._map._hbar.maximum() > 0
        assert panel._map._hbar.isVisible()
        # Nucleotide row is canvas-only; scrollbar sits below and must not clip letters.
        assert panel._map._canvas_row.height() == panel._map._canvas.height()
    finally:
        panel.close()
