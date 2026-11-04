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
    rgb[..., 0], rgb[..., 1], rgb[..., 2] = 10, 20, 30
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
