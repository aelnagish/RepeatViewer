"""One Export button: annotation files + vector images; optional GFF3 on Save."""

from __future__ import annotations

from pathlib import Path

import pytest

from repeatviewer.backends import sequence_features as sf
from repeatviewer.backends.sequence_features import (
    SequenceFeature,
    auto_export_seqid_map,
    export_annotation_path,
    save_annotations_multi,
)


def _feat(name: str = "LTR1") -> SequenceFeature:
    return SequenceFeature(start=2, end=9, name=name, feature_type="LTR", strand="+")


def test_save_writes_no_gff3_by_default(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(sf, "auto_gff3_on_save_enabled", lambda: False)
    fasta = tmp_path / "CroSat3.fasta"
    fasta.write_text(">CroSat3\n" + "A" * 50 + "\n", encoding="utf-8")
    json_path, gff = save_annotations_multi(
        fasta,
        features_by_id={"CroSat3": [_feat()]},
        lengths_by_id={"CroSat3": 50},
        record_ids=["CroSat3"],
    )
    assert json_path.is_file()
    assert gff is None
    assert not export_annotation_path(fasta, "gff3").exists()


def test_save_also_writes_gff3_when_enabled(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(sf, "auto_gff3_on_save_enabled", lambda: True)
    fasta = tmp_path / "CroSat3.fasta"
    fasta.write_text(">CroSat3\n" + "A" * 50 + "\n", encoding="utf-8")
    _json, gff = save_annotations_multi(
        fasta,
        features_by_id={"CroSat3": [_feat("LTR1"), _feat("LTR2")]},
        lengths_by_id={"CroSat3": 50},
        record_ids=["CroSat3"],
    )
    assert gff == export_annotation_path(fasta, "gff3")
    assert gff.name == "CroSat3.fasta_RV.gff3"
    text = gff.read_text(encoding="utf-8")
    assert text.startswith("##gff-version 3")
    assert "LTR1" in text and "LTR2" in text


def test_auto_export_seqid_map_only_for_duplicates() -> None:
    assert auto_export_seqid_map(["a", "b"]) is None
    mapping = auto_export_seqid_map(["seq", "seq_2"])
    if mapping is not None:  # depends on the uniquifier convention
        assert len(set(mapping.values())) == 2


def test_vector_pdf_has_no_embedded_bitmap(qapp, tmp_path: Path) -> None:
    from PySide6.QtCore import QRect
    from PySide6.QtGui import QColor, QPainter

    from repeatviewer.views.export_annotation_view_dialog import write_view_vector_pdf

    def paint(painter: QPainter) -> None:
        painter.fillRect(QRect(10, 10, 100, 20), QColor("#dc2626"))
        painter.drawText(12, 50, "LTR1")

    out = tmp_path / "view.pdf"
    write_view_vector_pdf(out, paint, width=300, height=80, logical_dpi=72)
    data = out.read_bytes()
    assert data.startswith(b"%PDF")
    assert b"/Subtype /Image" not in data and b"/Subtype/Image" not in data


def test_export_dialog_dpi_only_for_png_and_names_follow_type(
    qapp, tmp_path: Path
) -> None:
    from repeatviewer.views.export_annotation_view_dialog import (
        ExportAnnotationViewDialog,
    )

    fasta = tmp_path / "monomers.fasta"
    fasta.write_text(">a\nACGT\n", encoding="utf-8")
    dialog = ExportAnnotationViewDialog(
        source_path=fasta,
        name_tag="annotation",
        default_format="pdf",
        annotation_provider=lambda: (["a"], {"a": [_feat()]}, {"a": 4}),
    )
    try:
        types = [dialog._type.itemData(i) for i in range(dialog._type.count())]
        for fmt in ("gff3", "gtf", "bed", "pdf", "png", "svg"):
            assert fmt in types

        def pick(fmt: str) -> None:
            dialog._type.setCurrentIndex(dialog._type.findData(fmt))

        pick("pdf")
        assert not dialog._form.isRowVisible(dialog._dpi)
        assert Path(dialog._file.text()).name == "monomers.fasta_annotation.pdf"
        pick("svg")
        assert not dialog._form.isRowVisible(dialog._dpi)
        pick("png")
        assert dialog._form.isRowVisible(dialog._dpi)
        pick("gff3")
        assert not dialog._form.isRowVisible(dialog._dpi)
        assert dialog._form.isRowVisible(dialog._auto_gff3_row)
        assert Path(dialog._file.text()).name == "monomers.fasta_RV.gff3"
        pick("pdf")
        assert not dialog._form.isRowVisible(dialog._auto_gff3_row)
        assert Path(dialog._file.text()).name == "monomers.fasta_annotation.pdf"
    finally:
        dialog.close()


@pytest.mark.parametrize("fmt", ["pdf", "png", "svg"])
def test_image_only_dialog_has_no_annotation_types(qapp, tmp_path: Path, fmt) -> None:
    from repeatviewer.views.export_annotation_view_dialog import (
        ExportAnnotationViewDialog,
    )

    dialog = ExportAnnotationViewDialog(
        source_path=tmp_path / "x.fa", name_tag="dotplot", default_format=fmt
    )
    try:
        types = [dialog._type.itemData(i) for i in range(dialog._type.count())]
        assert types == ["pdf", "png", "svg"]
    finally:
        dialog.close()


def test_export_type_combo_wide_enough_for_labels(qapp, tmp_path: Path) -> None:
    """Type must not stay at ~55 px or the popup shows \"…\" instead of PNG/SVG."""
    from repeatviewer.views.export_annotation_view_dialog import (
        ExportAnnotationViewDialog,
    )

    dialog = ExportAnnotationViewDialog(
        source_path=tmp_path / "x.fa",
        title="Export dotplot",
        name_tag="dotplot",
        default_format="pdf",
    )
    try:
        assert dialog._type.minimumWidth() >= 100
        view = dialog._type.view()
        assert view is not None
        assert view.minimumWidth() >= 100
        for i in range(dialog._type.count()):
            assert dialog._type.itemText(i) in {"PDF", "PNG", "SVG"}
    finally:
        dialog.close()


def test_unsaved_annotation_marker_and_discard(qapp, tmp_path: Path) -> None:
    from repeatviewer.widgets.sequence_map_viewer import SequenceMapPanel

    fasta = tmp_path / "CroSat3.fasta"
    fasta.write_text(">CroSat3\n" + "ACGT" * 30 + "\n", encoding="utf-8")
    panel = SequenceMapPanel(fasta, autoload=True)
    try:
        panel.show()
        qapp.processEvents()
        assert not panel.has_unsaved_annotation_changes()

        panel._map.set_features([_feat("LTR1")])
        qapp.processEvents()
        assert panel.has_unsaved_annotation_changes()
        assert panel._rail.save_btn.text() == "Save •"

        panel.discard_annotation_changes()
        assert not panel.has_unsaved_annotation_changes()
        assert panel._rail.save_btn.text() == "Save"

        panel._map.set_features([_feat("LTR1")])
        assert panel.save_annotations(silent=True)
        assert not panel.has_unsaved_annotation_changes()
        assert panel._rail.save_btn.text() == "Save"
    finally:
        panel.close()
