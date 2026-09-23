"""Cytogenetic images UI: tick column, attach dialog, panel, markers, caption export."""

from __future__ import annotations

from pathlib import Path

import pytest
from PySide6.QtCore import QCoreApplication, QSettings
from PySide6.QtGui import QColor, QImage

from repeatviewer.backends.cyto_images import (
    FeatureRef,
    ImageTarget,
    add_image,
    images_for_sequence,
    load_images,
)


@pytest.fixture()
def isolated_settings(qapp, tmp_path, monkeypatch):
    """Keep QSettings writes out of the user's real preferences."""
    QSettings.setDefaultFormat(QSettings.Format.IniFormat)
    QSettings.setPath(QSettings.Format.IniFormat, QSettings.Scope.UserScope, str(tmp_path / "settings"))
    old_org = QCoreApplication.organizationName()
    QCoreApplication.setOrganizationName("RepeatViewerTests")
    yield
    QCoreApplication.setOrganizationName(old_org)


def _png(path: Path, color: str = "#ff0000", size: int = 40) -> Path:
    image = QImage(size, size, QImage.Format.Format_RGB32)
    image.fill(QColor(color))
    assert image.save(str(path))
    return path


def _fasta(folder: Path, name: str = "CroSat1.fasta") -> Path:
    path = folder / name
    path.write_text(">CroSat1\n" + "ACGT" * 50 + "\n", encoding="utf-8")
    return path


def test_file_list_tick_column(qapp, isolated_settings, tmp_path: Path) -> None:
    from repeatviewer.widgets.file_list_table import COL_IMAGES, FileListTableWidget

    fasta = _fasta(tmp_path)
    plain = _fasta(tmp_path, "Other.fasta")
    table = FileListTableWidget()
    table.add_file_row(fasta)
    table.add_file_row(plain)
    assert table.item(0, COL_IMAGES).text() == ""

    add_image(tmp_path, _png(tmp_path / "a.png"), targets=[ImageTarget(fasta.name)])
    add_image(tmp_path, _png(tmp_path / "b.png"), targets=[ImageTarget(fasta.name)])
    table.refresh_cyto_images()
    rows = {Path(table._path_for_row(r)).name: r for r in range(table.rowCount())}
    cell = table.item(rows[fasta.name], COL_IMAGES)
    assert cell.text() == "✓"
    assert "2 cytogenetic images" in cell.toolTip()
    assert table.item(rows[plain.name], COL_IMAGES).text() == ""

    opened: list[Path] = []
    table.cyto_images_requested.connect(opened.append)
    clicked: list[object] = []
    table.itemClicked.connect(clicked.append)
    table._on_cell_clicked(rows[fasta.name], COL_IMAGES)
    assert opened == [fasta.resolve()]
    assert clicked == []  # the tick opens images, not the normal row click
    table._on_cell_clicked(rows[plain.name], COL_IMAGES)
    assert len(clicked) == 1  # empty cell behaves like any other cell


def test_attach_dialog_save_as_is(qapp, isolated_settings, tmp_path: Path) -> None:
    from repeatviewer.widgets.cyto_images_ui import AttachImageDialog, cyto_notifier

    fasta = _fasta(tmp_path)
    source = _png(tmp_path / "fish.png")
    changed: list[str] = []
    cyto_notifier().changed.connect(changed.append)
    dialog = AttachImageDialog(
        tmp_path,
        source=source,
        targets=[ImageTarget(fasta.name, feature=FeatureRef("probeA", 5, 40))],
    )
    dialog._save()
    images = load_images(tmp_path)
    assert len(images) == 1
    target = images[0].targets[0]
    assert target.sequence == fasta.name
    assert target.feature == FeatureRef("probeA", 5, 40)
    # Only auto-filled date is stored; nothing else was typed.
    assert set(images[0].metadata) <= {"date"}
    assert changed and changed[-1] == str(tmp_path.resolve())


def test_attach_dialog_multi_probe_and_details(qapp, isolated_settings, tmp_path: Path) -> None:
    from repeatviewer.widgets.cyto_images_ui import AttachImageDialog

    fasta = _fasta(tmp_path)
    other = _fasta(tmp_path, "CroSat3.fasta")
    dialog = AttachImageDialog(tmp_path, source=_png(tmp_path / "fish.png"), targets=[ImageTarget(fasta.name)])
    dialog._probes.add_target(ImageTarget(other.name))
    dialog._probes.cellWidget(0, dialog._probes.COL_DYE).setCurrentText("Cy3")
    dialog._probes.cellWidget(1, dialog._probes.COL_DYE).setCurrentText("FITC")
    dialog._fields["species"].setText("Crocus sativus")
    dialog._fields["technique"].setCurrentText("FISH")
    dialog._save()
    image = images_for_sequence(other)[0]
    assert [t.sequence for t in image.targets] == [fasta.name, other.name]
    assert [t.dye for t in image.targets] == ["Cy3", "FITC"]
    assert image.targets[0].color  # dye picked a default colour
    assert image.metadata["species"] == "Crocus sativus"

    # Species is remembered for the next attach.
    again = AttachImageDialog(tmp_path, source=_png(tmp_path / "fish2.png"), targets=[ImageTarget(fasta.name)])
    assert again._fields["species"].text() == "Crocus sativus"


def test_panel_lists_and_removes(qapp, isolated_settings, tmp_path: Path, monkeypatch) -> None:
    from PySide6.QtWidgets import QMessageBox

    from repeatviewer.widgets.cyto_images_ui import CytogeneticsPanel

    fasta = _fasta(tmp_path)
    panel = CytogeneticsPanel(fasta)
    assert not panel._empty.isHidden()
    add_image(tmp_path, _png(tmp_path / "a.png"), targets=[ImageTarget(fasta.name, feature=FeatureRef("p", 1, 9))])
    from repeatviewer.widgets.cyto_images_ui import notify_changed

    notify_changed(tmp_path)
    assert panel._thumbs.count() == 1
    assert panel._empty.isHidden()
    sources: list[object] = []
    panel.showSourceRequested.connect(sources.append)
    panel._show_source()
    assert sources == [FeatureRef("p", 1, 9)]

    monkeypatch.setattr(QMessageBox, "question", lambda *a, **k: QMessageBox.StandardButton.Yes)
    panel._remove()
    assert load_images(tmp_path) == []
    assert panel._thumbs.count() == 0


def test_annotation_markers(qapp, isolated_settings, tmp_path: Path) -> None:
    from repeatviewer.backends.sequence_features import SequenceFeature, save_annotations
    from repeatviewer.widgets.sequence_map_viewer import SequenceMapPanel

    fasta = _fasta(tmp_path)
    features = [
        SequenceFeature(10, 60, "probeA", "others"),
        SequenceFeature(80, 120, "probeB", "others"),
        SequenceFeature(130, 150, "rep", "others"),
        SequenceFeature(160, 180, "rep", "others"),
    ]
    save_annotations(fasta, sequence_id="CroSat1", sequence_length=200, features=features)
    add_image(
        tmp_path,
        _png(tmp_path / "a.png"),
        targets=[
            ImageTarget(fasta.name, feature=FeatureRef("probeA", 10, 60)),
            # moved since attaching, unique name -> still marked
            ImageTarget(fasta.name, feature=FeatureRef("probeB", 85, 120)),
            # ambiguous name, coordinates gone -> not marked
            ImageTarget(fasta.name, feature=FeatureRef("rep", 1, 5)),
        ],
    )
    panel = SequenceMapPanel(fasta, autoload=True)
    data_features = panel._map.data().features
    marked = panel._marked_feature_indices(data_features)
    names = sorted(data_features[i].name for i in marked)
    assert names == ["probeA", "probeB"]
    assert panel.show_feature(FeatureRef("probeA", 10, 60)) is True


def test_render_with_caption_adds_strip(qapp) -> None:
    from repeatviewer.widgets.cyto_images_ui import render_with_caption

    image = QImage(300, 200, QImage.Format.Format_RGB32)
    image.fill(QColor("black"))
    out = render_with_caption(image, "FISH of CroSat1 (Cy3, red) on Crocus sativus chromosomes.")
    assert out.width() >= 300
    assert out.height() > 200


def test_viewer_and_preview_tab(qapp, isolated_settings, tmp_path: Path) -> None:
    from repeatviewer.widgets.cyto_images_ui import CytoImageViewer
    from repeatviewer.widgets.preview_tabs import FilePreviewPage

    fasta = _fasta(tmp_path)
    add_image(tmp_path, _png(tmp_path / "a.png"), targets=[ImageTarget(fasta.name, dye="Cy3", color="#ff3b30")])
    add_image(tmp_path, _png(tmp_path / "b.png", "#00ff00"), targets=[ImageTarget(fasta.name)])
    viewer = CytoImageViewer(tmp_path, images_for_sequence(fasta), 0)
    assert "CroSat1" in viewer._legend.text()
    viewer._show(1)
    assert viewer._counter.text() == "2 / 2"
    viewer.close()

    page = FilePreviewPage(fasta)
    try:
        assert page.is_single_fasta_preview()
        page.open_cyto_tab()
        assert page.is_cyto_preview_active()
        assert page._cyto_panel is not None and page._cyto_panel._thumbs.count() == 2
    finally:
        # Let the page's deferred loads run while it is alive, then dispose of it,
        # so no queued callback outlives it into later tests.
        from PySide6.QtCore import QEventLoop, QTimer

        loop = QEventLoop()
        QTimer.singleShot(200, loop.quit)
        loop.exec()
        page.deleteLater()
        qapp.processEvents()
