"""Import without dye names: numbered files, multichannel TIFFs, RGB merges."""

from __future__ import annotations

import numpy as np
import pytest

from repeatviewer.backends import cyto_images as ci
from repeatviewer.backends.cytogenetics import channels as ch
from repeatviewer.backends.cytogenetics import lists as ls
from repeatviewer.backends.cytogenetics import model as md
from repeatviewer.backends.cytogenetics.importer import ImportRequest, channel_layers, plan_cells, run_request, scan_tiffs

tifffile = pytest.importorskip("tifffile")


def test_numbered_files_one_folder_one_cell(tmp_path):
    src = tmp_path / "66-644"
    for cell in ("cell1", "cell2"):
        (src / cell).mkdir(parents=True)
        for n in range(1, 4):
            tifffile.imwrite(str(src / cell / f"{n}.tif"), np.full((20, 30), n * 100, dtype=np.uint16))
        (src / cell / "merge.psd").write_bytes(b"8BPS")
    cells, notes = plan_cells(scan_tiffs(src, recursive=True), ls.load_lists(tmp_path))
    assert [c.key for c in cells] == ["cell1", "cell2"] and notes == []
    assert [p.name for p, _d, _i in cells[0].entries()] == ["1.tif", "2.tif", "3.tif"]
    work = tmp_path / "data"
    work.mkdir()
    slide, images, errors = run_request(ImportRequest(folder=work, source=src, recursive=True, slide_name="S"))
    assert len(images) == 2 and not errors
    d = md.image_details(ci.load_images(work)[0])
    assert len(d.channels) == 3 and d.channels[0].dye == "" and d.channels[0].display.color != d.channels[1].display.color
    assert d.adjustment_log[1].startswith("Channel 1:")


def test_multichannel_imagej_tiff(tmp_path):
    stack = np.stack([np.full((20, 30), v, dtype=np.uint16) for v in (100, 200, 300)])
    path = tmp_path / "metaphase.tif"
    tifffile.imwrite(str(path), stack, imagej=True, metadata={"axes": "CYX", "Labels": ["DAPI", "FITC", "Cy3"]})
    info = ch.inspect_tiff(path)
    assert info.kind == "channels" and info.count == 3
    assert int(ch.read_plane(path, 2).array.mean()) == 300
    cells, _n = plan_cells([path], ls.load_lists(tmp_path))
    assert len(cells) == 1 and cells[0].dyes == ["DAPI", "FITC", "Cy3"]
    work = tmp_path / "data"
    work.mkdir()
    src = tmp_path / "src"
    src.mkdir()
    path.rename(src / path.name)
    _slide, images, errors = run_request(ImportRequest(folder=work, source=src, slide_name="S"))
    assert len(images) == 1 and not errors
    d = md.image_details(ci.load_images(work)[0])
    assert [c.index for c in d.channels] == [0, 1, 2] and len({c.file for c in d.channels}) == 1
    layers = channel_layers(work, d, {})
    assert [int(a.mean()) for a, _ in layers] == [100, 200, 300]


def test_rgb_merge_alone_is_split(tmp_path):
    rgb = np.zeros((20, 30, 3), dtype=np.uint8)
    rgb[..., 0] = np.arange(30, dtype=np.uint8)[None, :] * 3
    rgb[..., 1] = 20
    rgb[..., 2] = np.arange(20, dtype=np.uint8)[:, None] * 5
    (tmp_path / "s").mkdir()
    path = tmp_path / "s" / "overlay.tif"
    tifffile.imwrite(str(path), rgb, photometric="rgb")
    assert ch.inspect_tiff(path).kind == "rgb"
    cells, _n = plan_cells([path], ls.load_lists(tmp_path))
    assert [i for _p, _d, i in cells[0].entries()] == [0, 1, 2]
    assert int(ch.read_plane(path, 1).array.mean()) == 20


def test_split_picture_and_add_channel_image(tmp_path):
    from repeatviewer.backends.cytogenetics.image_session import ImageSession
    from repeatviewer.backends.cytogenetics.importer import add_channel_image

    rgb = np.zeros((20, 30, 3), dtype=np.uint8)
    rgb[..., 2] = 200
    pic = tmp_path / "merge.tif"
    tifffile.imwrite(str(pic), rgb, photometric="rgb")
    work = tmp_path / "data"
    work.mkdir()
    image = ci.add_image(work, pic, targets=[], metadata={})
    s = ImageSession(work, image)
    s.split_picture()
    s.save()
    again = ci.load_images(work)[0]
    d = md.image_details(again)
    assert [c.index for c in d.channels] == [0, 1, 2] and again.file.endswith("_preview.png")
    assert ci.image_file_path(work, again).exists()
    assert int(channel_layers(work, d, {})[2][0].mean()) == 200
    stack = np.stack([np.full((20, 30), v, dtype=np.uint16) for v in (100, 200)])
    multi = tmp_path / "two.tif"
    tifffile.imwrite(str(multi), stack, imagej=True, metadata={"axes": "CYX"})
    img = add_channel_image(work, multi)
    assert len(md.image_details(img).channels) == 2


def test_tinted_channels_sizes_and_merges(tmp_path):
    import numpy as np
    import tifffile

    from repeatviewer.backends.cytogenetics import channels as chn
    from repeatviewer.backends.cytogenetics.importer import plan_cells
    from repeatviewer.backends.cytogenetics.lists import default_lists

    cell = tmp_path / "cell1"
    cell.mkdir()
    rng = np.random.default_rng(1)
    base = (rng.random((40, 50)) * 255).astype(np.uint8)
    zero = np.zeros_like(base)
    tifffile.imwrite(cell / "1.tif", np.dstack([zero, zero, base]), photometric="rgb")  # blue DAPI
    tifffile.imwrite(cell / "2.tif", np.dstack([zero, base, zero]), photometric="rgb")  # green
    tifffile.imwrite(cell / "3.tif", np.dstack([base, base // 1, zero]), photometric="rgb")  # yellow
    other = (rng.random((40, 50, 3)) * 255).astype(np.uint8)
    tifffile.imwrite(cell / "merge.tif", other, photometric="rgb")  # a merge: planes differ
    tifffile.imwrite(cell / "figure.tif", (rng.random((90, 120)) * 255).astype(np.uint8))  # another size
    assert chn.inspect_tiff(cell / "1.tif").kind == "colour" and chn.inspect_tiff(cell / "1.tif").colours == ["#0000ff"]
    assert chn.inspect_tiff(cell / "3.tif").colours == ["#ffff00"]
    assert chn.inspect_tiff(cell / "merge.tif").kind == "rgb"
    cells, notes = plan_cells(sorted(cell.iterdir()), default_lists())
    assert [(g.key, g.size) for g in cells] == [("cell1", 3)]
    assert cells[0].colours == ["#0000ff", "#00ff00", "#ffff00"]
    assert "left out" in notes[0]
    cells, _notes = plan_cells(sorted(cell.iterdir()), default_lists(), merges=True)
    assert sorted(g.key for g in cells) == ["cell1", "cell1 (120×90)", "cell1 merge"]


def _spots(h, w, centres, value):
    plane = np.zeros((h, w), dtype=np.float32)
    yy, xx = np.mgrid[:h, :w]
    for y, x in centres:
        plane = np.maximum(plane, value * np.exp(-((yy - y) ** 2 + (xx - x) ** 2) / 18.0))
    return plane


def test_mostly_black_merge_is_not_a_tinted_channel(tmp_path):
    h, w = 200, 200
    dapi = _spots(h, w, [(50, 50), (150, 150), (60, 140)], 230)
    fitc = _spots(h, w, [(52, 54)], 255)
    cy3 = _spots(h, w, [(148, 146)], 255)
    tint = np.stack([dapi * 0.0, dapi * 0.5, dapi * 1.0], axis=-1).round().astype(np.uint8)
    assert ch.off_tint_share(tint) == 0.0
    assert ch.colour_of(tint) == "#0080ff"
    merge = np.clip(np.stack([cy3, fitc, dapi], axis=-1), 0, 255).round().astype(np.uint8)
    assert ch.off_tint_share(merge) > ch.MERGE_OFF_TINT_SHARE
    assert ch.colour_of(merge) is None
    cell = tmp_path / "cell1"
    cell.mkdir()
    tifffile.imwrite(str(cell / "1.tif"), tint, photometric="rgb")
    tifffile.imwrite(str(cell / "2.tif"), np.stack([fitc * 0, fitc, fitc * 0], -1).round().astype(np.uint8), photometric="rgb")
    tifffile.imwrite(str(cell / "Abdullah-101.tif"), merge, photometric="rgb")
    cells, notes = plan_cells(scan_tiffs(tmp_path, recursive=True), ls.load_lists(tmp_path))
    assert [p.name for p, _d, _i in cells[0].entries()] == ["1.tif", "2.tif"]
    assert any("Abdullah-101.tif" in n for n in notes)
