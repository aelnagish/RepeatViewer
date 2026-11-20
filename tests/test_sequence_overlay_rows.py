"""Rows above the bases: nothing the workspace draws may cover a letter or another graphic."""

import random
from itertools import combinations

import pytest
from PySide6.QtCore import QCoreApplication, QSettings
from PySide6.QtGui import QColor


@pytest.fixture()
def isolated_settings(qapp, tmp_path):
    org, app = QCoreApplication.organizationName(), QCoreApplication.applicationName()
    QCoreApplication.setOrganizationName("RVOverlayTests")
    QCoreApplication.setApplicationName(f"t{id(tmp_path)}")
    QSettings().clear()
    yield
    QSettings().clear()
    QCoreApplication.setOrganizationName(org)
    QCoreApplication.setApplicationName(app)


@pytest.fixture()
def canvas(qapp, isolated_settings, tmp_path):
    from repeatviewer.widgets.sequence_map_viewer import SequenceMapPanel

    r = random.Random(4)
    fasta = tmp_path / "s.fa"
    fasta.write_text(">s\n" + "".join(r.choice("ACGT") for _ in range(900)) + "\n")
    panel = SequenceMapPanel(fasta, autoload=True)
    panel.resize(1200, 500)
    panel.show()
    qapp.processEvents()
    yield panel._map._canvas
    panel.close()


def _base_top(canvas):
    return canvas._nucleotide_cell_geometry()[1]


def _no_overlaps(canvas):
    """No two graphics of the rows intersect, and none reaches into the bases."""
    boxes = canvas.overlay_boxes()
    top = _base_top(canvas)
    for name, rect in boxes:
        assert rect.top() >= -0.01, name
        assert rect.bottom() <= top + 0.01, f"{name} reaches into the bases"
    for (a, ra), (b, rb) in combinations(boxes, 2):
        assert not ra.adjusted(0.3, 0.3, -0.3, -0.3).intersects(rb), f"{a} covers {b}"
    return boxes


PAIRS = [(100, 119, 480, 499), (110, 129, 470, 489), (95, 114, 500, 519)]


def test_no_rows_without_workspace_graphics(canvas):
    assert canvas._overlay_band_height() == 0 and canvas.overlay_boxes() == []
    assert _base_top(canvas) == 0                      # the strip sits where it always did


def test_each_candidate_has_a_row_above_the_bases(canvas, qapp):
    from repeatviewer.widgets import sequence_map_viewer as smv

    plain = canvas._content_height()
    canvas.set_workbench_overlay((200, 400, "amplify"), PAIRS, 0)
    band = canvas._overlay_band_height()
    assert band == 3 * (smv._OVERLAY_LANE_H + smv._OVERLAY_LANE_GAP) + smv._OVERLAY_PAD
    assert _base_top(canvas) == band and canvas._content_height() == plain + band
    assert canvas.height() == plain + band             # the widget made room
    boxes = _no_overlaps(canvas)
    names = [n for n, _ in boxes]
    for i in (1, 2, 3):
        assert f"pair {i} forward" in names and f"pair {i} reverse" in names
    tops = {n: r.top() for n, r in boxes if n.endswith("forward")}
    assert len(set(tops.values())) == 3                # three pairs, three rows
    canvas.grab()                                      # paints without an error

    canvas.set_workbench_overlay((200, 400, "amplify"), [], -1)   # a region alone needs no row
    assert canvas._overlay_band_height() == 0 and canvas.height() == plain
    canvas.grab()


def test_edit_primer_enzymes_and_candidates_never_cover_each_other(canvas, qapp):
    canvas.set_workbench_overlay((200, 400, "amplify"), PAIRS, 1)
    canvas.set_edit_primer((100, 119), "+", False)
    canvas.set_primer_marks([105])
    sites = [(300, 303, 301, "MspI", False, ""), (302, 305, 303, "AluI", False, ""),
             (304, 307, 305, "HpaII", True, ""), (306, 309, 307, "TaqI", False, ""),
             (700, 703, 701, "MspI", False, "")]
    canvas.set_enzyme_sites(sites)
    boxes = _no_overlaps(canvas)
    names = [n for n, _ in boxes]
    assert "edit primer" in names and "enzyme MspI" in names
    # the primer being edited sits in the row right above the bases
    edit = next(r for n, r in boxes if n == "edit primer")
    assert max(r.bottom() for _n, r in boxes) == edit.bottom()
    x = (edit.left() + edit.right()) / 2
    assert canvas._edit_primer_hit(x, edit.center().y()) == "move"
    assert canvas._edit_primer_hit(x, _base_top(canvas) + 2) == "move"   # its bases count too
    canvas.grab()

    canvas.set_enzyme_sites([])
    canvas.set_edit_primer(None)
    canvas.set_workbench_overlay(None, [], -1)
    assert canvas._overlay_band_height() == 0 and _base_top(canvas) == 0
    canvas.grab()


def test_more_candidates_than_rows(canvas):
    from repeatviewer.widgets import sequence_map_viewer as smv

    many = [(50 + 5 * i, 69 + 5 * i, 600 + 5 * i, 619 + 5 * i) for i in range(12)]
    canvas.set_workbench_overlay(None, many, 10)
    rows = smv._OVERLAY_MAX_CANDIDATES
    assert canvas._overlay_rows() == (0, rows, 0)
    names = {n for n, _ in _no_overlaps(canvas)}
    shown = sorted(int(n.split()[1]) for n in names if n.endswith("forward"))
    assert shown == list(range(1, rows)) + [11]        # the last row is the selected pair


def test_feature_bars_keep_their_colour_family(qapp):
    from repeatviewer.widgets.sequence_map_viewer import _feature_paint_colors

    base = QColor("#1f57c3")
    fill, outline, text = _feature_paint_colors(base)
    assert fill.lightness() > base.lightness() > outline.lightness()
    assert fill.hue() == pytest.approx(base.hue(), abs=4)
    assert text == QColor("#1c1c1e")                   # dark name on the light fill
    faded = QColor(base)
    faded.setAlpha(90)
    assert _feature_paint_colors(faded)[0].alpha() == 90
    assert _feature_paint_colors(QColor("#101010"))[2] == QColor("#ffffff")


def test_drawing_a_feature_bar(qapp):
    from PySide6.QtCore import QRectF
    from PySide6.QtGui import QImage, QPainter

    from repeatviewer.widgets.sequence_map_viewer import _draw_directed_feature

    image = QImage(200, 40, QImage.Format.Format_ARGB32)
    image.fill(QColor("#ffffff"))
    painter = QPainter(image)
    _draw_directed_feature(painter, QRectF(10, 10, 150, 16), QColor("#1f57c3"),
                           forward=True, label="monomer")
    _draw_directed_feature(painter, QRectF(170, 10, 3, 16), QColor("#1f57c3"), forward=False)
    painter.end()
    inside = image.pixelColor(100, 12)
    assert inside != QColor("#ffffff") and inside.lightness() > QColor("#1f57c3").lightness()
    assert image.pixelColor(5, 18) == QColor("#ffffff")           # nothing outside the bar
    assert image.pixelColor(159, 11) == QColor("#ffffff")         # the corner is cut: a pointed end
