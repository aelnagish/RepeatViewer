"""Genome-browser tracks in the contig view: an .rvmap with an extra track, and a
reference opened on its own (.rvgb session) with annotation, PAF and BAM tracks."""

from __future__ import annotations

import json
import time

import pytest
from PySide6.QtCore import QPointF

from repeatviewer.backends import genome_tracks as gt
from repeatviewer.widgets import contig_view as cv
from tests.test_contig_view import _payload


def _wait_tracks(qapp, panel, n, timeout=10.0):
    t0 = time.time()
    while (len(panel.canvas.tracks) < n or not panel._tracks_btn.isEnabled()) \
            and time.time() - t0 < timeout:
        qapp.processEvents()
    qapp.processEvents()


def test_rvmap_with_an_annotation_track(qapp, tmp_path):
    path = tmp_path / "x.fa_map.rvmap"
    path.write_text(json.dumps(_payload()))
    bed = tmp_path / "te.bed"
    bed.write_text("ref\t2\t12\tTE1\t0\t-\n")
    panel = cv.ContigViewPanel(path)
    panel.resize(900, 600)
    panel.show()
    qapp.processEvents()
    c = panel.canvas
    before = c.header_h()
    c.add_track(gt.open_track(bed))
    assert c.header_h() > before
    top, _h = c._extra_spans[0]
    c.show_base(8, zoom_to_letters=True)
    tip = c.track_tooltip(QPointF(c.x_of(6), top + 4))
    assert tip and "TE1" in tip
    # reads still sit below the tracks
    assert c._track_y()["reads"] > top
    panel.grab()
    c.remove_track(0)
    assert c.header_h() == before


@pytest.fixture
def genome(tmp_path):
    seq = "ACGT" * 2500
    fa = tmp_path / "g.fa"
    fa.write_text(">chr1\n" + "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\n")
    (tmp_path / "a.gff3").write_text(
        "chr1\t.\tgene\t101\t900\t.\t+\t.\tID=g1;Name=GeneA\n"
        "chr1\t.\tLTR_retrotransposon\t2001\t3000\t.\t-\t.\tID=te1;Name=Gypsy-7\n")
    (tmp_path / "r.paf").write_text(
        "q1\t100\t0\t100\t+\tchr1\t10000\t200\t300\t99\t100\t60\tcs:Z::50*ac:49\n")
    return fa


def test_session_opens_reference_with_tracks(qapp, genome):
    session = gt.session_path_for(genome)
    gt.write_session(session, genome, [{"path": str(genome.parent / "a.gff3")},
                                       {"path": str(genome.parent / "r.paf")}])
    panel = cv.ContigViewPanel(session)
    panel.resize(1000, 600)
    panel.show()
    _wait_tracks(qapp, panel, 2)
    m = panel.current_model()
    assert m.browser and m.length == 10000 and m.reference.startswith("ACGT")
    assert [t.name for t in panel.canvas.tracks] == ["a.gff3", "r.paf"]
    assert not panel._tabs.isTabVisible(panel._tabs.indexOf(panel._reads_view))
    panel._goto.setText("Gypsy")
    panel._go_to()
    s, e = panel.canvas.visible_span()
    assert s <= 2000 and e >= 3000
    panel.canvas.show_range(240, 260)
    assert panel.canvas.letters()
    panel.grab()
    # the session remembers changed tracks
    panel._remove_track(0)
    saved = gt.read_session(session)
    assert [t["path"] for t in saved["tracks"]] == [str(genome.parent / "r.paf")]


def test_bam_track_in_browser(qapp, genome, tmp_path):
    pysam = pytest.importorskip("pysam")
    sam = tmp_path / "r.sam"
    sam.write_text("@SQ\tSN:chr1\tLN:10000\n"
                   "r1\t0\tchr1\t101\t60\t20M\t*\t0\t0\tACGTACGTACGTACGTTCGT\t*\n")
    bam = tmp_path / "r.bam"
    pysam.sort("-o", str(bam), str(sam))
    pysam.index(str(bam))
    session = gt.session_path_for(genome)
    gt.write_session(session, genome, [{"path": str(bam)}])
    panel = cv.ContigViewPanel(session)
    panel.resize(1000, 600)
    panel.show()
    _wait_tracks(qapp, panel, 1)
    t = panel.canvas.tracks[0]
    assert isinstance(t, gt.IndexedBamTrack)
    panel.canvas.show_range(95, 130)
    assert gt.mismatches_against(t.fetch("chr1", 95, 130)[0][0],
                                 panel.current_model().reference, 95, 130) == [(116, "T")]
    panel.grab()


def _browser(qapp, genome):
    session = gt.session_path_for(genome)
    gt.write_session(session, genome, [{"path": str(genome.parent / "a.gff3")},
                                       {"path": str(genome.parent / "r.paf")}])
    panel = cv.ContigViewPanel(session)
    panel.resize(1000, 600)
    panel.show()
    _wait_tracks(qapp, panel, 2)
    return session, panel


def test_drag_track_name_reorders_and_is_saved(qapp, genome):
    from PySide6.QtCore import QPoint, Qt
    from PySide6.QtTest import QTest

    session, panel = _browser(qapp, genome)
    c = panel.canvas
    vp = c.viewport()
    (top0, _h0), (top1, h1) = c._extra_spans
    QTest.mousePress(vp, Qt.MouseButton.LeftButton, Qt.KeyboardModifier.NoModifier,
                     QPoint(20, top0 + 5))
    QTest.mouseMove(vp, QPoint(20, top1 + h1 - 3))
    QTest.mouseRelease(vp, Qt.MouseButton.LeftButton, Qt.KeyboardModifier.NoModifier,
                       QPoint(20, top1 + h1 - 3))
    assert [t.name for t in c.tracks] == ["r.paf", "a.gff3"]
    assert [p_["path"].rsplit("/", 1)[-1] for p_ in gt.read_session(session)["tracks"]] == \
        ["r.paf", "a.gff3"]
    c.move_track(1, 0)
    assert [t.name for t in c.tracks] == ["a.gff3", "r.paf"]


def test_drag_on_ruler_zooms_to_region(qapp, genome):
    from PySide6.QtCore import QPoint, Qt
    from PySide6.QtTest import QTest

    _session, panel = _browser(qapp, genome)
    c = panel.canvas
    c.fit()
    qapp.processEvents()
    vp = c.viewport()
    x0, x1 = int(c.x_of(2000)), int(c.x_of(4000))
    QTest.mousePress(vp, Qt.MouseButton.LeftButton, Qt.KeyboardModifier.NoModifier,
                     QPoint(x0, 10))
    QTest.mouseMove(vp, QPoint(x1, 10))
    QTest.mouseRelease(vp, Qt.MouseButton.LeftButton, Qt.KeyboardModifier.NoModifier,
                       QPoint(x1, 10))
    # releasing zooms straight there: no selection chip, no Go
    assert c.region_sel is None and (c._sel_chip is None or not c._sel_chip.isVisible())
    s, e = c.visible_span()
    assert abs(s - 2000) < 150 and abs(e - 4000) < 150
    # the breadcrumb shows the sequence and the region in view
    assert panel._crumbs.chrom.text() == "chr1" and panel._crumbs.range.isVisible()
    # Esc clears a selection
    c.set_region_selection(100, 200)
    QTest.keyClick(c, Qt.Key.Key_Escape)
    assert c.region_sel is None
    # a sideways step is a fixed share of the view, not a few bases
    assert c.horizontalScrollBar().singleStep() >= c._data_w() // 20 - 1


def _mouse(c, kind, x, y):
    from PySide6.QtCore import QEvent, Qt
    from PySide6.QtGui import QMouseEvent

    t = {"press": QEvent.Type.MouseButtonPress, "move": QEvent.Type.MouseMove,
         "release": QEvent.Type.MouseButtonRelease}[kind]
    pos = QPointF(x, y)
    buttons = Qt.MouseButton.NoButton if kind == "release" else Qt.MouseButton.LeftButton
    ev = QMouseEvent(t, pos, c.viewport().mapToGlobal(pos), Qt.MouseButton.LeftButton,
                     buttons, Qt.KeyboardModifier.NoModifier)
    {"press": c.mousePressEvent, "move": c.mouseMoveEvent,
     "release": c.mouseReleaseEvent}[kind](ev)


def test_click_read_opens_details_and_drag_pans(qapp, genome):
    from repeatviewer.widgets.track_item_popup import TrackItemPopup

    _session, panel = _browser(qapp, genome)
    c = panel.canvas
    c.show_range(150, 450)
    qapp.processEvents()
    panel.grab()
    i = [t.name for t in c.tracks].index("r.paf")
    top, _h = c._extra_spans[i]
    y = top + 28 + 4 - c.verticalScrollBar().value()      # first read row
    got = []
    c.trackItemClicked.connect(lambda k, it, _p: got.append((k, it.name)))
    _mouse(c, "press", c.x_of(250), y)
    _mouse(c, "release", c.x_of(250), y)
    assert got == [(i, "q1")]
    popups = [w for w in qapp.topLevelWidgets() if isinstance(w, TrackItemPopup)
              and w.isVisible()]
    assert popups and "q1" in popups[-1].details_text()
    popups[-1].close()
    # drag to pan: no click, the view moves
    h0 = c.horizontalScrollBar().value()
    _mouse(c, "press", c.x_of(250), y)
    _mouse(c, "move", c.x_of(250) - 120, y)
    _mouse(c, "release", c.x_of(250) - 120, y)
    assert len(got) == 1
    assert c.horizontalScrollBar().value() > h0


def test_colour_by_settings_are_saved_and_drawn(qapp, genome):
    from repeatviewer.backends import alignment_style as style

    session, panel = _browser(qapp, genome)
    c = panel.canvas
    i = [t.name for t in c.tracks].index("r.paf")
    c.show_range(150, 450)
    for value, _label in style.COLOR_MODES:
        panel._set_track_setting(i, "color_by", value)
        panel.grab()
    panel._set_track_setting(i, "shade_mapq", False)
    saved = gt.read_session(session)["tracks"]
    paf = next(s for s in saved if s["path"].endswith("r.paf"))
    assert paf["settings"]["color_by"] == style.COLOR_NONE
    assert paf["settings"]["shade_mapq"] is False
    from PySide6.QtWidgets import QMenu

    menu = QMenu()
    panel._fill_track_menu(menu, i)
    labels = [a.text() for a in menu.actions()]
    assert "Colour by" in labels and "Shade MAPQ 0 reads" in labels


def test_browser_has_the_navigator_rail_not_the_side_panel(qapp, genome):
    _session, panel = _browser(qapp, genome)
    assert panel._rail.isVisible() and not panel._side.isVisible()
    assert panel._crumbs.isVisible()
    assert [c for c, _n in panel._rail.map.chroms] == ["chr1"]
    # a click on a track name (no drag) selects the track, without a menu
    from PySide6.QtCore import QPoint, Qt
    from PySide6.QtTest import QTest

    c = panel.canvas
    asked = []
    c.trackContextMenu.disconnect()
    c.trackContextMenu.connect(lambda k, _pos: asked.append(k))
    top0, _h0 = c._extra_spans[0]
    QTest.mouseClick(c.viewport(), Qt.MouseButton.LeftButton, Qt.KeyboardModifier.NoModifier,
                     QPoint(20, top0 + 5))
    assert asked == [] and c.selected_track == 0
    QTest.keyClick(c, Qt.Key.Key_Escape)
    assert c.selected_track is None


def test_rail_folds_hides_and_picks_the_density_track(qapp, genome):
    _session, panel = _browser(qapp, genome)
    r = panel._rail
    # the density comes from the first annotation track unless another is picked
    combo = r.density_track
    assert [combo.itemText(i) for i in range(combo.count())] == ["None", "a.gff3"]
    assert combo.currentText() == "a.gff3" and r.map.density
    combo.setCurrentIndex(0)
    assert r.map.density == {}
    assert gt.read_session(_session)["density_track"] == ""
    # sections fold, and the panel icon hides the whole rail
    r.genome.set_open(False)
    assert not r.genome.body.isVisible() and r.bm_section.body.isVisible()
    r.genome.set_open(True)
    panel._side_btn.setChecked(False)
    assert not r.isVisible() and not panel._side.isVisible()
    panel._side_btn.setChecked(True)
    assert r.isVisible()


def test_sequence_picker_is_a_short_filtered_list(qapp):
    from repeatviewer.widgets.genome_navigator import ChromPicker

    chroms = [(f"scaffold_{i}", 1000 + i) for i in range(300)]
    picked = []
    pk = ChromPicker(chroms, "scaffold_150")
    pk.picked.connect(picked.append)
    pk.show()
    assert pk.list.height() < 400 and pk.list.count() == 300
    assert pk.list.currentItem().data(0x0100) == "scaffold_150"
    pk.filter.setText("_29")
    pk._pick_first()
    assert picked == ["scaffold_29"]


def test_back_forward_recent_and_gene_chip(qapp, genome):
    session, panel = _browser(qapp, genome)
    c = panel.canvas
    panel._nav_go("chr1", 1000, 2000, "first")
    panel._nav_go("chr1", 5000, 6000, "second")
    assert [p_.label for p_ in panel._recent][:2] == ["second", "first"]
    assert panel._crumbs.back.isEnabled() and not panel._crumbs.forward.isEnabled()
    panel._nav_back()
    s, _e = c.visible_span()
    assert abs(s - 1000) < 100
    assert panel._crumbs.forward.isEnabled()
    panel._nav_forward()
    assert abs(c.visible_span()[0] - 5000) < 100
    # the recent places are kept in the session
    saved = gt.read_session(session)
    assert [r["label"] for r in saved["recent"]][:2] == ["second", "first"]
    # typing a feature name goes there and names the gene chip
    panel._goto.setText("Gypsy")
    panel._go_to()
    assert panel._goto.text() == ""
    assert panel._recent[0].label.startswith("Gypsy")


def test_bookmarks_are_flagged_and_saved(qapp, genome, monkeypatch):
    from repeatviewer.widgets import genome_navigator as nav

    session, panel = _browser(qapp, genome)
    c = panel.canvas
    c.show_range(2000, 3000)
    monkeypatch.setattr(nav.BookmarkDialog, "exec", lambda self: 1)
    c.set_region_selection(2100, 2200)
    panel._add_bookmark()
    assert len(panel._bookmarks) == 1
    bm = panel._bookmarks[0]
    assert (bm.chrom, bm.start, bm.end) == ("chr1", 2100, 2200)
    assert c.bookmark_flags and c.bookmark_flags[0][:2] == (2100, 2200)
    c.viewport().grab()                       # the flag and its band draw
    saved = gt.read_session(session)
    assert saved["bookmarks"][0]["start"] == 2100
    # rewriting the tracks keeps the bookmarks
    panel._remember_tracks()
    assert gt.read_session(session)["bookmarks"][0]["start"] == 2100
    panel._delete_bookmark(0)
    assert gt.read_session(session)["bookmarks"] == [] and not c.bookmark_flags
    # B adds one at the middle of the view
    from PySide6.QtCore import Qt
    from PySide6.QtTest import QTest

    QTest.keyClick(c, Qt.Key.Key_B)
    assert len(panel._bookmarks) == 1 and panel._bookmarks[0].end - panel._bookmarks[0].start == 1


def test_mismatch_details_draw(qapp, genome, tmp_path):
    """Base-quality fading, clip marks, insertion lengths and allele-coloured coverage
    all draw (smoke test at base zoom)."""
    pysam = pytest.importorskip("pysam")
    seq = "ACGT" * 2500
    lines = ["@SQ\tSN:chr1\tLN:10000"]
    for i in range(12):
        p = 100 + i * 3
        read = list(seq[p:p + 60])
        read[150 - p] = "T" if seq[150] != "T" else "G"            # a SNP at 150
        read = "".join(read[:20]) + "AAAAAAA" + "".join(read[20:53])
        lines.append(f"r{i}\t0\tchr1\t{p + 1}\t60\t20M7I23M10S\t*\t0\t0\t{read}\t"
                     + "I" * 30 + "#" * 30)
    sam = tmp_path / "m.sam"
    sam.write_text("\n".join(lines) + "\n")
    bam = tmp_path / "m.bam"
    pysam.sort("-o", str(bam), str(sam))
    pysam.index(str(bam))
    session = gt.session_path_for(genome)
    gt.write_session(session, genome, [{"path": str(bam)}])
    panel = cv.ContigViewPanel(session)
    panel.resize(1000, 600)
    panel.show()
    _wait_tracks(qapp, panel, 1)
    c = panel.canvas
    c.show_range(100, 220)
    img = c.viewport().grab().toImage()
    top, _h = c._extra_spans[0]
    v = c.verticalScrollBar().value()
    y = int(top + 20 - v)
    from repeatviewer.widgets.genome_tracks_view import _COV_FILL

    def px(bp):
        return img.pixelColor(int(c.x_of(bp) + c.ppb / 2), y).name()

    # all 12 reads carry T/G at 150: that coverage bar takes the base's colour
    assert px(140) == _COV_FILL.name()
    assert px(150) != _COV_FILL.name()
    for key in ("shade_quality", "flag_clips"):
        panel._set_track_setting(0, key, False)
        c.viewport().grab()


def _snp_bam(tmp_path):
    pysam = pytest.importorskip("pysam")
    seq = "ACGT" * 2500
    lines = ["@SQ\tSN:chr1\tLN:10000", "@RG\tID:a\tSM:s1", "@RG\tID:b\tSM:s2"]
    for i in range(16):
        p = 100 + i * 2
        read = list(seq[p:p + 60])
        if i % 2:                                     # odd reads carry T/G at 150
            read[150 - p] = "T" if seq[150] != "T" else "G"
        lines.append(f"r{i}\t{16 if i % 3 == 0 else 0}\tchr1\t{p + 1}\t60\t60M\t*\t0\t0\t"
                     f"{''.join(read)}\t*\tRG:Z:{'ab'[i % 2]}")
    sam = tmp_path / "s.sam"
    sam.write_text("\n".join(lines) + "\n")
    bam = tmp_path / "s.bam"
    pysam.sort("-o", str(bam), str(sam))
    pysam.index(str(bam))
    return bam


def _open_bam(qapp, genome, bam):
    session = gt.session_path_for(genome)
    gt.write_session(session, genome, [{"path": str(bam)}])
    panel = cv.ContigViewPanel(session)
    panel.resize(1000, 700)
    panel.show()
    _wait_tracks(qapp, panel, 1)
    return session, panel


def test_sort_by_base_puts_the_variant_reads_on_top(qapp, genome, tmp_path):
    session, panel = _open_bam(qapp, genome, _snp_bam(tmp_path))
    c = panel.canvas
    c.show_range(100, 220)
    c.viewport().grab()
    t = c.tracks[0]
    panel.sort_track(0, "base", 150)
    c.viewport().grab()
    placed = sorted((yoff, it.name) for it, yoff in c._placed(t, 150, 151) if yoff is not None)
    names = [n for _y, n in placed]
    odd = [n for n in names if int(n[1:]) % 2]
    assert names[:len(odd)] == odd and len(odd) == 8
    assert gt.read_session(session)["tracks"][0]["settings"]["sort"]["by"] == "base"
    panel._set_track_setting(0, "sort", None)
    assert "sort" not in t.settings


def test_group_by_sample_draws_bands(qapp, genome, tmp_path):
    _session, panel = _open_bam(qapp, genome, _snp_bam(tmp_path))
    c = panel.canvas
    c.show_range(100, 220)
    c.viewport().grab()
    h0 = c._extra_track_h(c.tracks[0])
    panel._set_track_setting(0, "group_by", "sample")
    c.viewport().grab()
    qapp.processEvents()
    groups, _h = c._layout(c.tracks[0])
    assert [g for g, _top, _n in groups] == ["s1", "s2"]
    by_group = {}
    for it, yoff in c._placed(c.tracks[0], 100, 220):
        by_group.setdefault("ab"[int(it.name[1:]) % 2], []).append(yoff)
    # s1 (read group a) reads sit above every s2 read
    assert max(y for y in by_group["a"] if y is not None) < \
        min(y for y in by_group["b"] if y is not None)
    assert c._extra_track_h(c.tracks[0]) != h0
    from PySide6.QtWidgets import QMenu

    menu = QMenu()
    c.track_menu_bp = 150
    panel._fill_track_menu(menu, 0)
    labels = [a.text() for a in menu.actions()]
    assert "Group by" in labels and "Sort reads at 151 by" in labels


def test_collapsed_pairs_and_filters_in_the_view(qapp, genome, tmp_path):
    from repeatviewer.widgets.genome_tracks_view import MODE_COLLAPSED

    session, panel = _open_bam(qapp, genome, _snp_bam(tmp_path))
    c = panel.canvas
    c.show_range(100, 220)
    t = c.tracks[0]
    expanded = c._extra_track_h(t)
    panel._set_track_mode(0, MODE_COLLAPSED)
    assert c._extra_track_h(t) < expanded
    # collapsed: thin touching rows (2 px), not one row
    ys = {y for _it, y in c._placed(t, 100, 220) if y is not None}
    assert c._row_h(t) == 2 and len(ys) > 1 and all(y % 2 == 0 for y in ys)
    c.viewport().grab()
    panel._set_track_setting(0, "pairs", True)
    panel._set_track_setting(0, "min_mapq", 61)            # every read has MAPQ 60
    c.viewport().grab()
    assert c._placed(t, 100, 220) == []
    insp = panel._track_inspector
    assert insp.pairs.isChecked() and insp.min_mapq.value() == 61
    saved = gt.read_session(session)["tracks"][0]["settings"]
    assert saved["pairs"] is True and saved["min_mapq"] == 61 and saved["mode"] == "collapsed"


@pytest.fixture(autouse=True)
def _own_browser_settings(tmp_path, monkeypatch):
    """The recent-genomes list goes to a test file, never the user's settings."""
    from PySide6.QtCore import QSettings

    import repeatviewer.widgets.genome_browser_window as gbw

    ini = str(tmp_path / "browser_settings.ini")
    monkeypatch.setattr(gbw, "_settings", lambda: QSettings(ini, QSettings.Format.IniFormat))


def test_genome_browser_reopens_the_last_genome(qapp, genome, tmp_path, monkeypatch):
    import repeatviewer.widgets.genome_browser_window as gbw

    other = tmp_path / "other.fa"
    other.write_text(">c\n" + "ACGT" * 100 + "\n")
    s1 = gt.session_path_for(genome)
    gt.write_session(s1, genome, [])
    s2 = gt.session_path_for(other)
    gt.write_session(s2, other, [])
    assert gbw.recent_sessions() == []
    gbw.remember_session(s1)
    gbw.remember_session(s2)
    gbw.remember_session(s1)
    assert gbw.recent_sessions() == [s1, s2]
    assert gbw.session_reference(s2) == other

    # Tools > Genome Browser opens the last one, whatever is selected
    opened = []

    class _Main:
        _show_genome_browser = staticmethod(lambda s, r, w=None: opened.append((s, r)))
        _pick_genome_for_browser = staticmethod(lambda w=None: opened.append("ask"))

    from repeatviewer.main_window import MainWindow

    MainWindow._open_genome_browser(_Main())
    assert opened == [(s1, genome)]
    monkeypatch.setattr(gbw, "recent_sessions", list)
    MainWindow._open_genome_browser(_Main())
    assert opened[-1] == "ask"


def test_genome_menu_lists_recent_genomes(qapp, genome):
    from repeatviewer.widgets.genome_browser_window import GenomeBrowserWindow

    session = gt.session_path_for(genome)
    gt.write_session(session, genome, [])
    win = GenomeBrowserWindow(session, genome)
    win.open_session()
    asked = []
    win.genomeRequested.connect(asked.append)
    win._fill_genome_menu()
    acts = [a for a in win.genome_menu.actions() if not a.isSeparator()]
    assert acts[0].text() == "Open genome…"
    assert [a.text() for a in acts if a.isCheckable()] == ["g.fa"]
    assert acts[-1].isChecked()
    acts[0].trigger()
    acts[-1].trigger()
    assert asked == [None, session]
    win.close()


def test_genome_browser_window_has_no_title_or_tabs_and_fills_the_width(qapp, genome):
    from repeatviewer.widgets.genome_browser_window import GenomeBrowserWindow

    session = gt.session_path_for(genome)
    gt.write_session(session, genome, [{"path": str(genome.parent / "a.gff3")}])
    win = GenomeBrowserWindow(session, genome)
    win.resize(500, 400)
    win.open_session()
    _wait_tracks(qapp, win.panel, 1)
    assert win.windowTitle() == "Genome Browser — g.fa"
    assert not win.panel._title.isVisible()
    assert not win.panel._tabs.tabBar().isVisible()
    c = win.panel.canvas
    win.resize(1400, 800)
    qapp.processEvents()
    # the whole reference still spans the (now wider) view
    assert c.x_of(c.model.length) > c.viewport().width() - 20
    win.close()


def test_gene_label_leaves_the_intron_visible(qapp, tmp_path):
    from PySide6.QtGui import QColor

    seq = "ACGT" * 2500
    fa = tmp_path / "g.fa"
    fa.write_text(">chr1\n" + seq + "\n")
    (tmp_path / "m.gff3").write_text(
        "chr1\t.\tmRNA\t101\t3000\t.\t+\t.\tID=m1;Name=BLGH_00875-mRNA-1\n"
        "chr1\t.\texon\t101\t200\t.\t+\t.\tParent=m1\n"
        "chr1\t.\texon\t2901\t3000\t.\t+\t.\tParent=m1\n"
        "chr1\t.\tCDS\t101\t200\t.\t+\t0\tParent=m1\n"
        "chr1\t.\tCDS\t2901\t3000\t.\t+\t0\tParent=m1\n")
    session = gt.session_path_for(fa)
    gt.write_session(session, fa, [{"path": str(tmp_path / "m.gff3")}])
    panel = cv.ContigViewPanel(session)
    panel.resize(1000, 400)
    panel.show()
    _wait_tracks(qapp, panel, 1)
    c = panel.canvas
    c.show_range(50, 3050)
    t = c.tracks[0]
    it = t.fetch("chr1", 0, 4000)[0][0]
    fm_w = 90
    _rect, _pen, over_line, _fill = c._gene_label(it, QColor("#4a90d9"), 0, 11, fm_w,
                                          c.x_of(it.start), c.x_of(it.end))
    assert over_line          # exons too narrow at this zoom: written over the line
    c.show_range(100, 210)    # zoomed into the first exon: the name goes inside it
    _rect, _pen, over_line, _fill = c._gene_label(it, QColor("#4a90d9"), 0, 11, fm_w,
                                          c.x_of(it.start), c.x_of(it.end))
    assert not over_line
    panel.grab()


def test_feature_popup_copies_its_sequence_and_files_drop_in_as_tracks(qapp, genome):
    from PySide6.QtCore import QMimeData, QUrl
    from PySide6.QtGui import QDropEvent, QGuiApplication

    from repeatviewer.widgets.track_item_popup import TrackItemPopup

    session = gt.session_path_for(genome)
    gt.write_session(session, genome, [])
    panel = cv.ContigViewPanel(session)
    panel.resize(1000, 500)
    panel.show()
    qapp.processEvents()
    mime = QMimeData()
    mime.setUrls([QUrl.fromLocalFile(str(genome.parent / "a.gff3")),
                  QUrl.fromLocalFile(str(genome))])          # the FASTA is not a track
    assert panel._dropped_tracks(mime) == [genome.parent / "a.gff3"]
    from PySide6.QtCore import QPointF, Qt

    ev = QDropEvent(QPointF(200, 200), Qt.DropAction.CopyAction, mime,
                    Qt.MouseButton.LeftButton, Qt.KeyboardModifier.NoModifier)
    panel.dropEvent(ev)
    _wait_tracks(qapp, panel, 1)
    assert [t.name for t in panel.canvas.tracks] == ["a.gff3"]
    gene = next(it for it, _r in panel.canvas.tracks[0].fetch("chr1", 0, 10000)
                if it.name == "GeneA")
    pop = TrackItemPopup(gene, chrom="chr1", fetch=panel._reference_fetcher())
    assert pop.copy_seq.isEnabled()
    fasta = pop.feature_fasta()
    assert fasta.startswith(">GeneA chr1:101-900(+)\n")
    assert "".join(fasta.splitlines()[1:]) == ("ACGT" * 2500)[100:900]
    pop.copy_seq.click()
    assert QGuiApplication.clipboard().text() == fasta


def test_trackpad_scrolls_the_browser_faster(qapp, genome, tmp_path):
    from PySide6.QtCore import QPoint, Qt
    from PySide6.QtGui import QWheelEvent

    _session, panel = _open_bam(qapp, genome, _snp_bam(tmp_path))
    c = panel.canvas
    c.resize(600, 150)
    qapp.processEvents()
    v = c.verticalScrollBar()
    v.setValue(0)
    ev = QWheelEvent(QPointF(300, 100), QPointF(300, 100), QPoint(0, -10), QPoint(0, -30),
                     Qt.MouseButton.NoButton, Qt.KeyboardModifier.NoModifier,
                     Qt.ScrollPhase.ScrollUpdate, False)
    c.wheelEvent(ev)
    assert v.maximum() > 30 and v.value() == 30


def test_added_gene_track_goes_above_the_reads(qapp, genome, tmp_path):
    session, panel = _open_bam(qapp, genome, _snp_bam(tmp_path))
    panel._load_tracks([genome.parent / "a.gff3"])
    _wait_tracks(qapp, panel, 2)
    assert panel.canvas.tracks[0].kind == gt.KIND_ANNOTATION
    saved = [p_["path"].rsplit("/", 1)[-1] for p_ in gt.read_session(session)["tracks"]]
    assert saved[0] == "a.gff3"


def test_dragging_a_track_near_the_bottom_scrolls(qapp, genome, tmp_path):
    from PySide6.QtCore import QPoint, Qt
    from PySide6.QtTest import QTest

    _session, panel = _open_bam(qapp, genome, _snp_bam(tmp_path))
    panel._load_tracks([genome.parent / "a.gff3"])
    _wait_tracks(qapp, panel, 2)
    c = panel.canvas
    c.resize(600, 160)
    qapp.processEvents()
    v = c.verticalScrollBar()
    v.setValue(0)
    vp = c.viewport()
    top0, _h0 = c._extra_spans[0]
    QTest.mousePress(vp, Qt.MouseButton.LeftButton, Qt.KeyboardModifier.NoModifier,
                     QPoint(20, top0 + 3))
    QTest.mouseMove(vp, QPoint(20, vp.height() - 5))
    assert v.value() > 0
    QTest.mouseRelease(vp, Qt.MouseButton.LeftButton, Qt.KeyboardModifier.NoModifier,
                       QPoint(20, vp.height() - 5))


def test_track_is_highlighted_while_its_menu_is_open(qapp, genome, monkeypatch):
    from PySide6.QtCore import QPoint
    from PySide6.QtWidgets import QMenu

    _session, panel = _browser(qapp, genome)
    c = panel.canvas
    seen = []

    class _Menu(QMenu):
        def exec(self, *_a):
            seen.append(c.menu_track)

    monkeypatch.setattr(cv, "QMenu", _Menu)
    panel._track_menu(1, QPoint(0, 0))
    assert seen == [1] and c.menu_track is None
    c.set_menu_track(0)
    c.viewport().grab()                       # the highlight draws
    c.set_menu_track(None)


def test_frame_codons_on_both_strands():
    from repeatviewer.widgets.contig_view import frame_codons

    ref = "CCATGAAATAGCC"
    assert frame_codons(ref, "+", 2, 0, 13) == [(2, "M"), (5, "K"), (8, "*")]
    # reverse complement GGCTATTTCATGG read from the end: G Y F M
    assert frame_codons(ref, "−", 0, 0, 13) == [(1, "M"), (4, "F"), (7, "Y"), (10, "G")]


def test_reference_menu_adds_reverse_complement_and_translation(qapp, genome):
    _session, panel = _browser(qapp, genome)
    c = panel.canvas
    before = c._track_y()["reads"]
    panel._set_sequence_rows(revcomp=True, frames="all")
    ys = c._track_y()
    assert "revcomp" in ys and "frames" in ys and len(c.translation_frames()) == 6
    assert ys["reads"] > before
    c.show_range(100, 160)                    # letters and amino acids draw
    c.viewport().grab()
    c.show_range(0, 9000)                     # small boxes for starts / stops only
    c.viewport().grab()
    assert c.complement_seq()[:4] == c.model.reference[:4].translate(
        str.maketrans("ACGT", "TGCA"))
    panel._set_sequence_rows(revcomp=False, frames="off")
    assert c._track_y()["reads"] == before


def test_hover_text_stays_narrow():
    from repeatviewer.widgets import genome_tracks_view as gtv

    long = ",".join(f"InterPro:IPR{i:06d}" for i in range(40))
    info = f"geneID BLGH_02691, Dbxref {gtv._short(long, gtv.TIP_VALUE_CHARS)}, Ontology {long}"
    rows = gtv._wrap(info)
    assert len(rows) <= gtv.TIP_INFO_LINES
    assert all(len(r) <= gtv.TIP_LINE_CHARS for r in rows)
    assert rows[-1].endswith("…")
    assert gtv._short("short", 10) == "short"


def test_bookmarks_from_bed_and_back():
    from repeatviewer.widgets.genome_navigator import (
        Bookmark,
        bookmarks_from_bed,
        bookmarks_to_bed,
    )

    text = (
        "track name=x itemRgb=On\n"
        "# comment\n"
        "chr1\t99\t200\n"                                   # BED3: no name -> locus
        "chr1\t10\t20\tTE1\t0\t+\t10\t20\t0,128,255\n"       # BED9 colour
        "chr1\t30\t40\tLTR\tGypsy element near gene\n"         # description in col 5
        "chrX\t1\t5\tnope\n"                                  # not in this genome
        "chr1 500 20000 tail\n"                                 # spaces; end clipped
        "bad line\n")
    bms, skipped = bookmarks_from_bed(text, {"chr1": 10000})
    assert [b.name for b in bms] == ["chr1:100-200", "TE1", "LTR", "tail"]
    assert (bms[0].start, bms[0].end) == (99, 200)
    assert bms[1].color == "#0080ff"
    assert bms[2].note == "Gypsy element near gene"
    assert bms[3].end == 10000
    assert len(skipped) == 2 and "chrX" in skipped[0]
    out = bookmarks_to_bed(bms + [Bookmark("p", "chr1", 5, 6, note="a\tb")])
    again, sk = bookmarks_from_bed(out)
    assert not sk
    assert [(b.chrom, b.start, b.end, b.name) for b in again] == [
        (b.chrom, b.start, b.end, b.name) for b in bms] + [("chr1", 5, 6, "p")]
    assert again[1].color == "#0080ff" and again[2].note == "Gypsy element near gene"
    assert again[-1].note == "a b"


def test_import_bookmarks_from_bed_file(qapp, genome, tmp_path, monkeypatch):
    from PySide6.QtWidgets import QMessageBox

    monkeypatch.setattr(QMessageBox, "information", lambda *a, **k: None)
    _s, panel = _browser(qapp, genome)
    bed = tmp_path / "marks.bed"
    bed.write_text("chr1\t100\t300\tmyTE\t0\t-\t100\t300\t255,0,0\n")
    panel._import_bookmarks(str(bed))
    assert [b.name for b in panel._bookmarks] == ["myTE"]
    assert panel._rail.bookmarks[0].color == "#ff0000"
    panel._import_bookmarks(str(bed))           # the same line twice: not duplicated
    assert len(panel._bookmarks) == 1
    out = tmp_path / "out.bed"
    panel._export_bookmarks(str(out))
    assert "chr1\t100\t300\tmyTE" in out.read_text()



def test_no_read_hover_when_the_view_is_too_wide_for_reads(qapp, genome, tmp_path):
    from PySide6.QtCore import QPointF

    from repeatviewer.widgets import genome_tracks_view as gtv

    _session, panel = _open_bam(qapp, genome, _snp_bam(tmp_path))
    c = panel.canvas
    t = c.tracks[0]
    c.show_range(100, 220)
    c.viewport().grab()
    top, _h = c._extra_spans[0]
    yoff = next(y for _it, y in c._placed(t, 150, 151) if y is not None)
    y = top + gtv.ALN_COV_H + 2 + yoff + 1 - c.verticalScrollBar().value()
    pos = QPointF(c.x_of(150) + 1, y)
    assert c.track_items_at(pos)[1]
    # the view is wider than the read window: only "Zoom in" is drawn, no hover
    t.read_window = 50
    c.viewport().grab()
    assert c.track_items_at(pos)[1] == []
