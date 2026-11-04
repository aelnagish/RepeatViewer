"""Multi-sequence Dotplot Grid (paged thumbnails) + Single mode."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import QApplication

from repeatviewer.widgets.dotplot_gallery import DotplotGallery
from repeatviewer.widgets.dotplot_view import DotplotPanel


def test_gallery_pages_twenty_then_next(qapp: QApplication) -> None:
    seqs = {f"s{i}": "ACGTACGTACGTACGT" for i in range(45)}

    def loader(sid: str) -> str:
        return seqs[sid]

    gallery = DotplotGallery()
    gallery.set_sequences(
        list(seqs),
        sequence_loader=loader,
        word_size=4,
        include_reverse=False,
        page_size=20,
    )
    assert gallery.page_count() == 3
    assert gallery.page() == 0
    assert len(gallery._thumbs) == 20
    assert gallery.page_size() == 20
    assert DotplotGallery.max_page_size() == 100

    gallery.set_page(1)
    assert gallery.page() == 1
    assert len(gallery._thumbs) == 20

    gallery.set_page(2)
    assert gallery.page() == 2
    assert len(gallery._thumbs) == 5

    gallery.set_page_size(10)
    assert gallery.page_size() == 10
    assert gallery.page_count() == 5

    gallery.set_thumb_size(200)
    assert gallery.thumb_size() == 200
    # Plot face is square: both X and Y follow the computed side length.
    gallery.resize(900, 600)
    qapp.processEvents()
    gallery._rebuild_page()
    assert gallery._thumbs
    thumb = gallery._thumbs[0]
    side = thumb._plot_size
    assert thumb.width() == side + 8
    assert thumb.height() == side + 22 + 10
    gallery.close()


def test_dotplot_panel_defaults_to_grid_and_returns(
    qapp: QApplication, tmp_path: Path
) -> None:
    path = tmp_path / "multi.fa"
    lines = []
    for i in range(25):
        lines.append(f">seq{i}")
        lines.append("ACGTACGTACGTACGTACGT")
    path.write_text("\n".join(lines) + "\n")

    panel = DotplotPanel(path)
    panel.set_record_choices(
        [f"seq{i}" for i in range(25)],
        sequence_loader=lambda sid: "ACGTACGTACGTACGTACGT",
    )
    panel.show()
    qapp.processEvents()
    panel.ensure_loaded()
    qapp.processEvents()

    assert not panel._rail.view_mode_combo.isVisible()
    assert panel._view_mode == "grid"
    assert panel._plot_stack.currentIndex() == 1
    assert panel._gallery.page_size() == 20
    assert len(panel._gallery._thumbs) == 20

    panel._gallery.set_page(1)
    qapp.processEvents()
    assert panel._gallery.page() == 1

    panel._on_gallery_thumb_clicked("seq22")
    qapp.processEvents()
    assert panel._view_mode == "single"
    assert panel._plot_stack.currentIndex() == 0
    assert panel.sequence_id() == "seq22"
    assert panel._rail.return_grid_btn.isVisible()

    panel._return_to_gallery()
    qapp.processEvents()
    assert panel._view_mode == "grid"
    assert panel._plot_stack.currentIndex() == 1
    # seq22 is on page 1 with page_size 20 (indices 20-24).
    assert panel._gallery.page() == 1
    assert any(t.sequence_id() == "seq22" for t in panel._gallery._thumbs)

    panel.close()


def test_open_single_for_sequence_from_grid(
    qapp: QApplication, tmp_path: Path
) -> None:
    """Overview 'Show self-dotplot' must open Single, not stay on Grid."""
    path = tmp_path / "multi.fa"
    lines = []
    for i in range(10):
        lines.append(f">seq{i}")
        lines.append("ACGTACGTACGTACGTACGT")
    path.write_text("\n".join(lines) + "\n")

    panel = DotplotPanel(path)
    panel.set_record_choices(
        [f"seq{i}" for i in range(10)],
        sequence_loader=lambda sid: "ACGTACGTACGTACGTACGT",
    )
    panel.show()
    qapp.processEvents()
    panel.ensure_loaded()
    qapp.processEvents()
    assert panel._view_mode == "grid"

    panel.open_single_for_sequence("seq7")
    qapp.processEvents()
    assert panel._view_mode == "single"
    assert panel._plot_stack.currentIndex() == 0
    assert panel.sequence_id() == "seq7"
    assert panel._rail.return_grid_btn.isVisible()
    panel.close()


def test_gallery_word_size_change_recomputes(qapp: QApplication) -> None:
    """Changing Word size must not keep stale thumb results (unlike resize)."""
    seqs = {
        "s0": "ACGTACGTACGTACGTACGTACGTACGTACGT",
        "s1": "ACGTACGTACGTACGTACGTACGTACGTACGT",
    }

    def loader(sid: str) -> str:
        return seqs[sid]

    gallery = DotplotGallery()
    gallery.set_sequences(
        list(seqs),
        sequence_loader=loader,
        word_size=6,
        include_reverse=False,
        page_size=20,
    )
    # Drain background compute so thumbs have results for word=6.
    for _ in range(20):
        qapp.processEvents()
        if gallery._pending:
            gallery._compute_next()
        elif all(t._result is not None for t in gallery._thumbs):
            break
    assert gallery._thumbs
    assert all(t._result is not None for t in gallery._thumbs)
    assert all(int(t._result.word_size) == 6 for t in gallery._thumbs)

    gallery.set_compute_params(word_size=10, include_reverse=False)
    # Stale word=6 results must not be reused; thumbs should recompute.
    assert gallery._word == 10
    assert gallery._pending == ["s0", "s1"] or all(
        t._busy or (t._result is not None and int(t._result.word_size) == 10)
        for t in gallery._thumbs
    )
    for _ in range(20):
        qapp.processEvents()
        if gallery._pending:
            gallery._compute_next()
        elif all(
            t._result is not None and int(t._result.word_size) == 10
            for t in gallery._thumbs
        ):
            break
    assert all(
        t._result is not None and int(t._result.word_size) == 10
        for t in gallery._thumbs
    )
    gallery.close()


def test_gallery_line_width_updates_thumbs(qapp: QApplication) -> None:
    seqs = {f"s{i}": "ACGTACGTACGTACGT" for i in range(5)}

    def loader(sid: str) -> str:
        return seqs[sid]

    gallery = DotplotGallery()
    gallery.set_sequences(
        list(seqs),
        sequence_loader=loader,
        word_size=4,
        include_reverse=False,
        page_size=20,
    )
    assert gallery.line_width_step() == 8
    assert all(t.line_width_step() == 8 for t in gallery._thumbs)

    gallery.set_line_width_step(12)
    assert gallery.line_width_step() == 12
    assert all(t.line_width_step() == 12 for t in gallery._thumbs)
    gallery.close()


def test_gallery_set_selected_sequence_highlights_thumb(qapp: QApplication) -> None:
    seqs = {f"s{i}": "ACGTACGTACGTACGT" for i in range(6)}

    def loader(sid: str) -> str:
        return seqs[sid]

    gallery = DotplotGallery()
    gallery.set_sequences(
        list(seqs),
        sequence_loader=loader,
        word_size=4,
        include_reverse=False,
        page_size=20,
    )
    gallery.set_selected_sequence("s3")
    assert gallery.selected_sequence() == "s3"
    selected = [t for t in gallery._thumbs if t.is_selected()]
    assert len(selected) == 1
    assert selected[0].sequence_id() == "s3"
    gallery.set_selected_sequence("s1")
    assert [t.sequence_id() for t in gallery._thumbs if t.is_selected()] == ["s1"]
    gallery.close()


def test_gallery_page_survives_set_record_choices_resync(
    qapp: QApplication, tmp_path: Path
) -> None:
    """Tab switches re-call set_record_choices; Grid page must not reset."""
    path = tmp_path / "multi.fa"
    lines = []
    for i in range(45):
        lines.append(f">seq{i}")
        lines.append("ACGTACGTACGTACGTACGT")
    path.write_text("\n".join(lines) + "\n")
    ids = [f"seq{i}" for i in range(45)]

    panel = DotplotPanel(path)
    panel.set_record_choices(
        ids,
        sequence_loader=lambda sid: "ACGTACGTACGTACGTACGT",
    )
    panel.show()
    qapp.processEvents()
    panel.ensure_loaded()
    qapp.processEvents()

    panel._gallery.set_page(2)
    qapp.processEvents()
    assert panel._gallery.page() == 2

    # Simulate leaving Dotplot and returning (same id list).
    panel.set_record_choices(
        ids,
        preferred_id="seq40",
        sequence_loader=lambda sid: "ACGTACGTACGTACGTACGT",
    )
    qapp.processEvents()
    assert panel._view_mode == "grid"
    assert panel._gallery.page() == 2

    # Single mode must also survive a resync.
    panel._on_gallery_thumb_clicked("seq40")
    qapp.processEvents()
    assert panel._view_mode == "single"
    panel.set_record_choices(
        ids,
        preferred_id="seq40",
        sequence_loader=lambda sid: "ACGTACGTACGTACGTACGT",
    )
    qapp.processEvents()
    assert panel._view_mode == "single"
    assert panel.sequence_id() == "seq40"

    panel.close()
