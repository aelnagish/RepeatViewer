"""Annotation view export uses content height, not the stretched Preview pane."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import QApplication

from repeatviewer.backends.sequence_features import SequenceFeature


def test_annotation_export_height_excludes_viewport_blank(
    qapp: QApplication, tmp_path: Path
) -> None:
    from repeatviewer.widgets.sequence_map_viewer import SequenceMapPanel

    fasta = tmp_path / "short.fa"
    fasta.write_text(">s1\nATGAAACCCGGGTTTAAA\n", encoding="utf-8")
    panel = SequenceMapPanel(fasta, autoload=True)
    try:
        # Tall Preview pane leaves empty paper under a short annotation track.
        panel.resize(900, 700)
        panel.show()
        qapp.processEvents()
        panel._map._sync_header_views()
        qapp.processEvents()

        panel._map.set_features(
            [
                SequenceFeature(
                    start=1,
                    end=6,
                    name="ORF1",
                    feature_type="ORF",
                    strand="+",
                )
            ]
        )
        panel._map._sync_header_views()
        qapp.processEvents()

        track = panel._map._annot_track
        seq_block_h = int(panel._map._sequence_block.height())
        content_h = int(track.content_height())
        assert seq_block_h > content_h + 80  # stretched empty region present on screen

        _header, _canvas, _track, _w, height, header_h, canvas_h = (
            panel._annotation_view_export_layout()
        )
        assert height == header_h + canvas_h + content_h
        assert height < seq_block_h + header_h
    finally:
        panel.close()


def test_write_view_image_pdf_page_matches_image(tmp_path: Path) -> None:
    from PySide6.QtGui import QImage, QColor
    from repeatviewer.views.export_annotation_view_dialog import write_view_image_pdf

    img = QImage(120, 40, QImage.Format.Format_ARGB32)
    img.fill(QColor(255, 0, 0))
    out = tmp_path / "fig.pdf"
    write_view_image_pdf(
        out,
        img,
        dpi=100,
        content_width=9999,
        content_height=9999,
    )
    assert out.is_file() and out.stat().st_size > 0
