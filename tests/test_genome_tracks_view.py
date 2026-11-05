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
    # the region is selected, not zoomed to: its Go button zooms there
    assert c.visible_span()[1] - c.visible_span()[0] > 9000
    a, b = c.region_sel
    assert abs(a - 2000) < 150 and abs(b - 4000) < 150
    chip = c._sel_chip
    assert chip.isVisible() and "chr1:" in chip.label.text()
    panel.grab()
    chip.go.click()
    assert c.region_sel is None and not chip.isVisible()
    s, e = c.visible_span()
    assert abs(s - 2000) < 150 and abs(e - 4000) < 150
    # the locus box shows the region in view
    assert panel._goto.text().startswith("chr1:")
    assert gt.parse_locus(panel._goto.text())[1] == s + 1
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


def test_side_panel_track_options(qapp, genome):
    from repeatviewer.backends import alignment_style as style

    session, panel = _browser(qapp, genome)
    insp = panel._track_inspector
    assert panel._sections["tracks"].isVisible()
    assert not panel._sections["display"].isVisible()
    assert insp.track.count() == 2
    assert panel.canvas.tracks[insp.current()].name == "r.paf"   # alignment track first
    insp.color.setCurrentIndex(insp.color.findData(style.COLOR_TAG))
    assert insp.tag.isVisible()
    t = panel.canvas.tracks[insp.current()]
    assert t.settings["color_by"] == style.COLOR_TAG and t.settings["color_tag"] == "HP"
    insp.rows.setCurrentIndex(insp.rows.findData(50))
    assert t.settings["max_rows"] == 50
    insp.track.setCurrentIndex(0)                                  # the annotation track
    assert not insp.color.isVisible() and insp.mode.isVisible()
    saved = gt.read_session(session)["tracks"]
    assert next(s for s in saved if s["path"].endswith("r.paf"))["settings"]["max_rows"] == 50


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
