"""Figure builder logic: panels, layout, legend, methods, store."""

from __future__ import annotations

import numpy as np
import pytest

from repeatviewer.backends import cyto_images as ci
from repeatviewer.backends.cytogenetics import figure as fg
from repeatviewer.backends.cytogenetics import karyotype as kt
from repeatviewer.backends.cytogenetics import model as md
from repeatviewer.backends.cytogenetics.importer import ImportRequest, run_request

tifffile = pytest.importorskip("tifffile")


def _folder(tmp_path, cells=1):
    src = tmp_path / "scope"
    src.mkdir()
    for n in range(1, cells + 1):
        for dye in ("DAPI", "FITC", "Cy3"):
            arr = (np.random.default_rng(n).random((60, 80)) * 3000).astype(np.uint16)
            tifffile.imwrite(str(src / f"s12_c0{n}_{dye}.tif"), arr)
    work = tmp_path / "data"
    work.mkdir()
    run_request(ImportRequest(folder=work, source=src, slide_name="Slide 12",
                              slide_values={"species": "Crocus sativus", "technique": "FISH", "counterstain": "DAPI"},
                              stage="Mitotic metaphase", probes_by_dye={"FITC": "CroSat2", "Cy3": "45S rDNA"},
                              pixel_size_um=0.0645))
    k = kt.new_karyotype("Cs", 8, 3, species="Crocus sativus")
    kt.save_karyotype(work, k)
    return work


def test_panels_for_cell_and_labels(tmp_path):
    work = _folder(tmp_path)
    image = ci.load_images(work)[0]
    fig = fg.Figure()
    panels = fg.add_cell(fig, image)
    d = md.image_details(image)
    assert len(panels) == 4 and panels[-1].view == fg.VIEW_MERGE
    assert [fg.auto_label(p, d) for p in panels] == ["DAPI", "CroSat2", "45S rDNA", "Merge"]
    assert fg.label_colour(panels[0], d) != "#FFFFFF" and fg.label_colour(panels[3], d) == "#FFFFFF"
    assert fg.source_text(panels[1], image, "Slide 12") == "Slide 12 cell 01 · FITC"
    assert [fg.letter(i) for i in (0, 3, 25, 26)] == ["A", "D", "Z", "AA"]


def test_render_and_crop(tmp_path):
    work = _folder(tmp_path)
    image = ci.load_images(work)[0]
    fig = fg.Figure()
    panels = fg.add_cell(fig, image)
    rgb = fg.render_panel(work, image, panels[1], (0.0, 0.0, 0.5, 0.5))
    assert rgb.shape == (30, 40, 3)
    merged = fg.render_panel(work, image, panels[3], (0, 0, 1, 1))
    assert merged.shape == (60, 80, 3)
    assert fg.normalise_crop(0.5, 0.5, 0.1, 0.2) == (0.1, 0.2, 0.5, 0.5)
    assert fg.normalise_crop(0.1, 0.1, 0.105, 0.5) is None
    frac = fg.scale_bar_fraction(md.image_details(image), (0, 0, 0.5, 1), 80, 1.0)
    assert frac == pytest.approx(1.0 / (40 * 0.0645))


def test_layout(tmp_path):
    fig = fg.Figure(width_mm=174, spacing=2)
    fig.panels = [fg.Panel(id=str(i)) for i in range(4)]
    boxes, height = fg.layout(fig, {str(i): 0.75 for i in range(4)})
    assert len(boxes) == 4
    w = (174 - 1.0) / 2
    assert boxes[1][1].x == pytest.approx(w + 1.0) and boxes[2][1].y == pytest.approx(w * 0.75 + 1.0)
    assert height == pytest.approx(2 * w * 0.75 + 1.0)
    fig.layout = "4row"
    _boxes, h2 = fg.layout(fig, {str(i): 0.75 for i in range(4)})
    assert h2 < height
    fg.add_idiogram(fig)
    fig.layout = "3+idiogram"
    boxes, _h = fg.layout(fig, {str(i): 0.75 for i in range(4)})
    idio = boxes[-1][1]
    assert idio.w == pytest.approx(174) and idio.x == 0


def test_move_remove(tmp_path):
    fig = fg.Figure()
    fig.panels = [fg.Panel(id=x) for x in "abc"]
    fg.move_panel(fig, "c", -1)
    assert [p.id for p in fig.panels] == ["a", "c", "b"]
    fg.remove_panel(fig, "a")
    assert [p.id for p in fig.panels] == ["c", "b"]


def test_legend_and_methods(tmp_path):
    work = _folder(tmp_path)
    image = ci.load_images(work)[0]
    fig = fg.Figure(number="Figure 3.")
    fg.add_cell(fig, image)
    fig.panels[3].arrows.append(fg.Arrow(0.5, 0.5))
    html = fg.legend_html(work, fig)
    text = fg.html_to_text(html)
    assert text.startswith("Figure 3. FISH of CroSat2 (green, FITC) and 45S rDNA (red, Cy3) on mitotic metaphase chromosomes of Crocus sativus (2n = 3x = 24)")
    assert "counterstained with DAPI (blue)" in text
    assert "(A) DAPI, (B) CroSat2, (C) 45S rDNA, (D) merge; arrowheads" in text
    assert "adjusted linearly" in text and text.endswith("Scale bar 10 µm.")
    assert "<i>Crocus sativus</i>" in html
    methods = fg.methods_text(work, fig)
    assert "separate 16-bit image" in methods and "Raw image files were not changed" in methods
    assert "Pixel size 0.0645 µm" in methods


def test_store_roundtrip(tmp_path):
    fig = fg.Figure(name="Figure CroSat2", layout="1row", same_crop=False, crop=(0.1, 0.1, 0.9, 0.8))
    fig.panels = [fg.Panel(id="p1", image_id="i", view="c1", arrows=[fg.Arrow(0.2, 0.3, 200)], crop=(0, 0, 0.5, 0.5))]
    fg.save_figure(tmp_path, fig)
    [again] = fg.list_figures(tmp_path)
    assert again.to_dict() == fig.to_dict()
    assert fg.new_figure_name(tmp_path) == "Figure 1"
    fg.save_figure(tmp_path, fg.Figure(name="Figure 1"))
    assert fg.new_figure_name(tmp_path) == "Figure 2"
    fig.name = "Renamed"
    fg.save_figure(tmp_path, fig, old_name="Figure CroSat2")
    assert sorted(f.name for f in fg.list_figures(tmp_path)) == ["Figure 1", "Renamed"]
    assert fg.delete_figure(tmp_path, "Renamed")


def test_layout_choices_and_font():
    from repeatviewer.backends.cytogenetics import figure as fg2

    fig = fg2.Figure()
    fig.panels = [fg2.Panel(id=str(i), image_id="x") for i in range(5)]
    for layout, cols in (("1row", 5), ("1col", 1), ("2x3", 3), ("3x2", 2), ("3x3", 3)):
        fig.layout = layout
        assert fg2.columns_for(fig) == cols
    fig.layout, fig.rows, fig.cols = "custom", 1, 4
    assert fg2.columns_for(fig) == 4
    boxes, _h = fg2.layout(fig, {})
    assert len({round(b.x, 3) for _p, b in boxes}) == 4
    again = fg2.Figure.from_dict({"layout": "4row", "font_family": "Helvetica", "cols": 5})
    assert again.layout == "1row" and again.font_family == "Helvetica" and again.cols == 5
