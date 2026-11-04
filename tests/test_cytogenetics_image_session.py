"""Image tab logic: display changes, channels, scale bar, apply to slide."""

from __future__ import annotations

import numpy as np
import pytest

from repeatviewer.backends import cyto_images as ci
from repeatviewer.backends.cytogenetics import image_session as sess
from repeatviewer.backends.cytogenetics import model as md
from repeatviewer.backends.cytogenetics.importer import ImportRequest, run_request

tifffile = pytest.importorskip("tifffile")


def _tiff(path, seed, shape=(40, 50)):
    rng = np.random.default_rng(seed)
    tifffile.imwrite(str(path), (rng.random(shape) * 3000).astype(np.uint16), resolution=(155000, 155000), metadata={"unit": "um"})


def _session_folder(tmp_path, cells=3):
    src = tmp_path / "scope"
    src.mkdir()
    for n in range(1, cells + 1):
        for dye in ("DAPI", "FITC"):
            _tiff(src / f"s12_c0{n}_{dye}.tif", n)
    work = tmp_path / "data"
    work.mkdir()
    slide, images, errors = run_request(ImportRequest(folder=work, source=src, slide_name="Slide 12", probes_by_dye={"FITC": "CroSat2"}))
    assert not errors and len(images) == cells
    return work, slide, images


def test_levels_visibility_and_save(tmp_path):
    work, _slide, images = _session_folder(tmp_path)
    s = sess.ImageSession(work, ci.load_images(work)[0])
    assert s.is_multichannel and len(s.channels) == 2
    before = s.merged()
    s.set_visible(1, False)
    assert not np.array_equal(before, s.merged())
    s.set_visible(1, True)
    s.set_levels(0, 100, 900)
    s.set_gamma(0, 1.5)
    s.set_color(0, "#ff0000")
    s.save()
    again = sess.ImageSession(work, ci.load_images(work)[0])
    d = again.channels[0].display
    assert (d.low, d.high, d.gamma, d.color) == (100, 900, 1.5, "#ff0000")
    assert any("gamma 1.5" in line for line in again.details.adjustment_log)
    assert s.image.extra["scale_bar"]["show"] is True


def test_auto_levels_and_range(tmp_path):
    work, _slide, _images = _session_folder(tmp_path)
    s = sess.ImageSession(work, ci.load_images(work)[0])
    s.set_levels(0, 5, 6)
    s.auto_levels(0)
    d = s.channels[0].display
    assert d.high > 100 and d.gamma == 1.0
    assert s.level_range(0) == (0.0, 65535.0)


def test_side_by_side_and_labels(tmp_path):
    work, _slide, _images = _session_folder(tmp_path)
    s = sess.ImageSession(work, ci.load_images(work)[0])
    view = s.side_by_side()
    assert view.ndim == 3 and view.shape[2] == 3 and view.shape[1] >= 2 * 50
    assert s.panel_labels() == ["DAPI", "FITC · CroSat2", "Merge"]


def test_add_and_remove_channel(tmp_path):
    work, _slide, _images = _session_folder(tmp_path)
    s = sess.ImageSession(work, ci.load_images(work)[0])
    extra = tmp_path / "x_Cy3.tif"
    _tiff(extra, 9)
    c = s.add_channel(extra, probe="45S rDNA")
    assert c.dye == "Cy3" and len(s.channels) == 3
    stored = md.channel_path(work, c)
    assert stored.exists()
    wrong = tmp_path / "wrong.tif"
    _tiff(wrong, 1, shape=(10, 10))
    with pytest.raises(ValueError):
        s.add_channel(wrong)
    s.remove_channel(2)
    assert len(s.channels) == 2 and not stored.exists()


def test_scale_bar():
    assert sess.nice_length(13.0) == 10.0 and sess.nice_length(0.1) == 0.5
    assert sess.parse_length("10 µm") == 10.0 and sess.parse_length("Auto") == 0.0


def test_scale_bar_session(tmp_path):
    work, _slide, _images = _session_folder(tmp_path)
    s = sess.ImageSession(work, ci.load_images(work)[0])
    s.set_pixel_size(0.0645)
    length, px = s.scale_bar(2000)
    assert length == 20.0 and px == pytest.approx(20 / 0.0645)
    s.scale_length_um = 5
    assert s.scale_bar(2000)[0] == 5
    s.set_pixel_size(None)
    assert s.scale_bar(2000) is None


def test_apply_to_slide_and_neighbours(tmp_path):
    work, _slide, _images = _session_folder(tmp_path, cells=3)
    first = ci.load_images(work)[0]
    s = sess.ImageSession(work, first)
    s.set_levels(0, 50, 777)
    s.save()
    assert s.slide_cell_count() == 3
    assert s.apply_to_slide() == 2
    for other in ci.load_images(work):
        d = md.image_details(other).channels[0].display
        assert (d.low, d.high) == (50, 777)
    pos, total, prev_id, next_id = sess.neighbours(work, first.id)
    assert (pos, total, prev_id) == (1, 3, "") and next_id
