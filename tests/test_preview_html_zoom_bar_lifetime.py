"""Regression: HTML open must not delete the shared AnnotationZoomBar."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import QApplication, QVBoxLayout, QWidget
from shiboken6 import isValid

from repeatviewer.widgets.preview_controls import AnnotationZoomBar, PreviewHeaderBar, PreviewZoomBar
from repeatviewer.widgets.preview_tabs import PreviewTabWidget


def test_replace_fasta_overview_with_html_keeps_annotation_zoom_bar(
    qapp: QApplication, tmp_path: Path
) -> None:
    fasta = tmp_path / "aln.fa"
    fasta.write_text(">a\nACGT\n>b\nACGT\n", encoding="utf-8")
    html = tmp_path / "report.html"
    html.write_text("<html><body>ok</body></html>", encoding="utf-8")

    host = QWidget()
    layout = QVBoxLayout(host)
    preview_zoom = PreviewZoomBar()
    ann_zoom = AnnotationZoomBar()
    header = PreviewHeaderBar(preview_zoom, annotation_zoom_bar=ann_zoom)
    tabs = PreviewTabWidget(
        zoom_bar=preview_zoom,
        header_bar=header,
        annotation_zoom_bar=ann_zoom,
    )
    layout.addWidget(header)
    layout.addWidget(tabs)
    host.show()
    qapp.processEvents()

    assert tabs.open_file(fasta, replace=False) is True
    page = tabs.currentWidget()
    assert page is not None
    page.show_annotation_tab()
    qapp.processEvents()
    overview = getattr(page, "_annotation_overview", None)
    assert overview is not None
    overview.bind_annotation_zoom_bar(ann_zoom)
    qapp.processEvents()
    assert ann_zoom.parent() is overview._toolbar_slot

    assert tabs.open_file(html, replace=True) is True
    qapp.processEvents()

    assert isValid(ann_zoom)
    # Bar returned to header host (or at least not destroyed with the old page).
    assert ann_zoom.parent() is not None
    host.close()
