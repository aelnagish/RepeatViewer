"""Light genome browser preview: reference matching, loading and drawing (offscreen)."""

from __future__ import annotations

import os
import random
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from repeatviewer.backends.light_browser import genome_tracks as gt  # noqa: E402
from repeatviewer.backends.light_browser import reference as refmod  # noqa: E402


def _fasta(path: Path, seqs: dict[str, str]) -> Path:
    with path.open("w") as fh:
        for name, s in seqs.items():
            fh.write(f">{name} description\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")
    return path


def _sample(tmp_path: Path, contig_prefix: str = "chr"):
    rng = random.Random(1)
    seqs = {f"{contig_prefix}{i}": "".join(rng.choice("ACGT") for _ in range(6000))
            for i in (1, 2)}
    ref = _fasta(tmp_path / "ref.fa", seqs)
    lines = [f"@SQ\tSN:{n}\tLN:{len(s)}" for n, s in seqs.items()]
    for k in range(200):
        n = f"{contig_prefix}{1 + k % 2}"
        pos = rng.randrange(0, 5800)
        read = list(seqs[n][pos:pos + 100])
        read[50] = "A" if read[50] != "A" else "C"
        flag = 16 if k % 3 == 0 else 0
        lines.append(f"r{k}\t{flag}\t{n}\t{pos + 1}\t60\t100M\t*\t0\t0\t{''.join(read)}\t*")
    sam = tmp_path / "reads.sam"
    sam.write_text("\n".join(lines) + "\n")
    gff = tmp_path / "genes.gff3"
    gff.write_text(
        "##gff-version 3\n"
        f"{contig_prefix}1\tsrc\tgene\t101\t2000\t.\t+\t.\tID=g1;Name=GeneOne\n"
        f"{contig_prefix}1\tsrc\tmRNA\t101\t2000\t.\t+\t.\tID=t1;Parent=g1\n"
        f"{contig_prefix}1\tsrc\texon\t101\t500\t.\t+\t.\tParent=t1\n"
        f"{contig_prefix}1\tsrc\texon\t1001\t2000\t.\t+\t.\tParent=t1\n"
        f"{contig_prefix}1\tsrc\tCDS\t301\t500\t.\t+\t0\tParent=t1\n"
        f"{contig_prefix}1\tsrc\tCDS\t1001\t1800\t.\t+\t1\tParent=t1\n"
        f"{contig_prefix}2\tsrc\tgene\t3001\t4000\t.\t-\t.\tID=g2;Name=GeneTwo\n")
    return ref, sam, gff, seqs


class _Project:
    def __init__(self, genome: Path | None) -> None:
        self.genome = genome
        self.managed_genomes = []

    def genome_targets(self):
        return [("", self.genome)] if self.genome else []


def test_name_matcher_accepts_chr_prefix_and_versions():
    find = refmod.NameMatcher(["1", "2", "NC_000.1"])
    assert find("chr1") == "1"
    assert find("2") == "2"
    assert find("NC_000") == "NC_000.1"
    assert find("chrX") is None


def test_project_genome_is_chosen_when_names_match(tmp_path):
    ref, sam, _gff, _ = _sample(tmp_path)
    track = gt.open_track(sam)
    names = refmod.track_sequence_names(track)
    cands = refmod.reference_candidates(_Project(ref))
    best, tried = refmod.choose_reference(names, cands)
    assert best is not None and best.path == ref and best.matched == 2
    lengths = refmod.sequence_lengths(track, best)
    assert lengths == {"chr1": 6000, "chr2": 6000}


def test_mismatched_project_genome_is_reported(tmp_path):
    other_dir = tmp_path / "other"
    other_dir.mkdir()
    other, _s, _g, _ = _sample(other_dir, contig_prefix="scaffold_x")
    _ref, sam, _gff, _ = _sample(tmp_path)
    track = gt.open_track(sam)
    best, tried = refmod.choose_reference(refmod.track_sequence_names(track),
                                          refmod.reference_candidates(_Project(other)))
    assert best is None and len(tried) == 1 and tried[0].matched == 0
    # lengths still come from the SAM header
    assert refmod.sequence_lengths(track, None)["chr1"] == 6000


def test_a_few_lenient_name_hits_are_not_a_match(tmp_path):
    # scaffold_1 ~ Chr1 after normalising: one hit of many names is a coincidence
    fa = _fasta(tmp_path / "plant.fa", {"Chr1": "ACGT" * 10, "Chr2": "ACGT" * 10})
    names = [f"scaffold_{i}" for i in range(1, 20)]
    m = refmod.match_reference(fa, names)
    assert m.matched == 2 and not m.ok
    best, _tried = refmod.choose_reference(names, [("Project genome", fa)])
    assert best is None


def test_sequences_keep_the_genome_order(tmp_path):
    # reference order chr10, chr2, chr1: never sorted by name
    fa = _fasta(tmp_path / "g.fa", {"chr10": "A" * 50, "chr2": "C" * 50, "chr1": "G" * 50})
    gff = tmp_path / "a.gff3"
    gff.write_text("##gff-version 3\n"
                   "chr1\ts\tgene\t1\t10\t.\t+\t.\tID=a\n"
                   "chr2\ts\tgene\t1\t10\t.\t+\t.\tID=b\n"
                   "chr10\ts\tgene\t1\t10\t.\t+\t.\tID=c\n")
    track = gt.open_track(gff)
    ref = refmod.match_reference(fa, refmod.track_sequence_names(track))
    assert refmod.genome_order(track, ref) == ["chr10", "chr2", "chr1"]
    # without a reference: the order of the file
    assert refmod.genome_order(track, None) == ["chr1", "chr2", "chr10"]


def test_remembered_reference_comes_first(tmp_path):
    ref, _sam, _gff, _ = _sample(tmp_path)
    other = _fasta(tmp_path / "b.fa", {"x": "ACGT"})
    cands = refmod.reference_candidates(_Project(other), remembered=ref)
    assert [p for _l, p in cands] == [ref, other]


def test_mismatch_helper_finds_the_changed_base(tmp_path):
    from repeatviewer.widgets.light_genome_browser import _mismatches

    ref, sam, _gff, seqs = _sample(tmp_path)
    track = gt.open_track(sam)
    (it, _row), *_ = track.fetch("chr1", 0, 6000)
    w0 = max(0, it.start - 10)
    window = seqs["chr1"][w0:it.end + 10]
    got = _mismatches(it, window, w0, it.start, it.end)
    assert [pos for pos, _b in got] == [it.start + 50]


@pytest.fixture()
def qapp():
    from PySide6.QtWidgets import QApplication

    return QApplication.instance() or QApplication([])


def _wait(qapp, cond, timeout=10.0):
    import time

    end = time.monotonic() + timeout
    while time.monotonic() < end:
        qapp.processEvents()
        if cond():
            return True
        time.sleep(0.01)
    return False


@pytest.mark.parametrize("which", ["sam", "gff", "bam"])
def test_panel_loads_and_draws(qapp, tmp_path, monkeypatch, which):
    from repeatviewer.widgets import light_genome_browser as lgb

    ref, sam, gff, _ = _sample(tmp_path)
    path = {"sam": sam, "gff": gff}.get(which)
    if which == "bam":
        pysam = pytest.importorskip("pysam")
        path = tmp_path / "reads.bam"
        pysam.sort("-o", str(path), str(sam))
        pysam.index(str(path))
    monkeypatch.setattr(lgb, "_load_project", lambda: _Project(ref))
    monkeypatch.setattr(lgb, "remembered_reference", lambda p: None)
    monkeypatch.setattr(lgb, "_text_mode_kinds", lambda: set())
    panel = lgb.LightGenomeBrowserPanel()
    panel.resize(900, 500)
    panel.show()
    panel.load(path)
    assert _wait(qapp, lambda: panel._result is not None)
    assert panel._result.error == ""
    assert panel._result.reference is not None
    assert panel.canvas.chrom == "chr1"
    assert not panel._banner.isVisible()
    # Alignments and annotations open at Fit (whole chromosome).
    assert panel.canvas.view_range() == (0, panel.canvas.length)
    panel.canvas.show_range(0, 1500)
    img = panel.canvas.grab().toImage()
    assert img.width() > 0
    out = os.environ.get("LIGHT_BROWSER_SHOTS")
    if out:
        _wait(qapp, lambda: False, 0.5)
        panel.grab().save(str(Path(out) / f"light_{which}.png"))
    # zoom to bases: the sequence lane draws letters
    panel.canvas.show_range(200, 260)
    panel.canvas.repaint()
    # ruler drag marks a region; Select zooms to it
    panel.canvas.set_selection((300, 400))
    panel.canvas.go_to_selection()
    assert panel.canvas.view_range() == (300, 400) and panel.canvas.sel is None
    panel._on_locus_entered("1,001-2,000")       # the position box: current sequence
    assert panel.canvas.view_range() == (1000, 2000)
    panel._set_mode(text=True, remember=False)
    assert panel._text.toPlainText()
    panel.shutdown()


def test_panel_banner_when_reference_does_not_match(qapp, tmp_path, monkeypatch):
    from repeatviewer.widgets import light_genome_browser as lgb

    other_dir = tmp_path / "other"
    other_dir.mkdir()
    other, _s, _g, _ = _sample(other_dir, contig_prefix="scaffold_x")
    _ref, sam, _gff, _ = _sample(tmp_path)
    monkeypatch.setattr(lgb, "_load_project", lambda: _Project(other))
    monkeypatch.setattr(lgb, "remembered_reference", lambda p: None)
    monkeypatch.setattr(lgb, "_text_mode_kinds", lambda: set())
    panel = lgb.LightGenomeBrowserPanel()
    panel.resize(900, 500)
    panel.show()
    panel.load(sam)
    assert _wait(qapp, lambda: panel._result is not None)
    assert panel._result.reference is None
    assert panel._banner.isVisible()
    assert panel._banner_label.text() == "Please check if the reference genome is correct."
    # Mismatched reference: warn and do not draw the track.
    assert panel._browser_stack.currentIndex() == 0
    assert panel.canvas.track is None
    panel.shutdown()


def test_one_per_row_names_and_menu(qapp, tmp_path, monkeypatch):
    from PySide6.QtCore import QPointF

    from repeatviewer.widgets import light_genome_browser as lgb

    ref, sam, _gff, _ = _sample(tmp_path)
    monkeypatch.setattr(lgb, "_load_project", lambda: _Project(ref))
    monkeypatch.setattr(lgb, "remembered_reference", lambda p: None)
    monkeypatch.setattr(lgb, "_text_mode_kinds", lambda: set())
    panel = lgb.LightGenomeBrowserPanel()
    panel.resize(900, 500)
    panel.show()
    panel.load(sam)
    assert _wait(qapp, lambda: panel._result is not None and panel.canvas.track is not None)
    panel.canvas.one_per_row = True
    panel._sync_name_column()
    panel.canvas.repaint()
    assert panel.canvas._names_on
    assert panel.canvas.view_range() == (0, panel.canvas.length)
    # one per row keeps every sequence on the scaffold, not only the window
    assert len(panel.canvas._named_in_view(0, 20)) == 100
    panel.canvas.show_range(0, 20)
    panel.canvas.repaint()
    assert panel.canvas._names_on
    assert len(panel.canvas._name_hits) > 5
    panel._sync_name_column()
    panel.canvas.repaint()
    assert panel.canvas._name_hits
    it = panel.canvas.name_at(panel.canvas._name_hits[0][0].center())
    assert it is not None and it.name
    # the name sits in the left column, not on the read body
    assert panel.canvas._name_hits[0][0].left() == 0
    assert panel.canvas._name_hits[0][0].right() < lgb.NAME_W
    texts = [a.text() for a in panel._read_name_menu(it).actions()]
    assert "Copy name" in texts
    assert "Copy sequence (FASTA)" in texts
    assert "Extract sequence to a file…" in texts
    assert "Zoom to this sequence" in texts
    assert "Hide sequence" in texts
    assert "BLAST sequence…" in texts
    before = len(panel.canvas._named_in_view(0, panel.canvas.length))
    panel.canvas.hide_read(it)
    assert not panel.canvas.read_passes(it)
    assert len(panel.canvas._named_in_view(0, panel.canvas.length)) == before - 1
    panel.canvas.show_hidden_reads()
    assert panel.canvas.read_passes(it)
    from PySide6.QtCore import QEvent, QPointF, Qt
    from PySide6.QtGui import QMouseEvent

    pos = panel.canvas._name_hits[0][0].center()
    panel.canvas.mousePressEvent(QMouseEvent(
        QEvent.Type.MouseButtonPress, QPointF(pos), QPointF(pos),
        Qt.MouseButton.LeftButton, Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier))
    assert panel.canvas.picked is it or panel.canvas.picked.name == it.name
    panel._sync_view_menu()
    labels = [a.text() for a in panel._view_menu.actions() if not a.isSeparator()]
    assert labels[:2] == ["Browser", "Text"]
    panel._set_view("Text")
    assert panel._view_btn.text().startswith("Text")
    panel._set_view("Browser")
    assert panel._stack.currentIndex() == 0
    # more mapped reads than the cap: the names stay off and One per row is grey
    monkeypatch.setattr(lgb, "ONE_PER_ROW_LIMIT", 10)
    panel._sync_one_per_row_choice()
    row_item = panel._run_panel._layout_mode.model().item(1)
    assert row_item is not None and not row_item.isEnabled()
    assert "10" in panel._run_panel._layout_mode.toolTip()
    panel.canvas.repaint()
    assert not panel.canvas._names_on
    assert panel.canvas._name_hits == []
    assert panel.canvas.name_at(QPointF(8, panel.canvas._rows_top() + 4)) is None
    panel.shutdown()


def test_sam_length_mismatch_is_not_ok(tmp_path):
    """Same chromosome names but different @SQ / FASTA lengths → not a match."""
    ref = _fasta(tmp_path / "ref.fa", {"chr1": "A" * 1000, "chr2": "C" * 1000})
    sam = tmp_path / "reads.sam"
    sam.write_text(
        "@SQ\tSN:chr1\tLN:5000\n"
        "@SQ\tSN:chr2\tLN:5000\n"
        "r1\t0\tchr1\t1\t60\t10M\t*\t0\t0\tAAAAAAAAAA\t*\n"
        "r2\t0\tchr2\t1\t60\t10M\t*\t0\t0\tCCCCCCCCCC\t*\n"
    )
    track = gt.open_track(sam)
    names = refmod.track_sequence_names(track)
    lengths, auth = refmod.track_sequence_lengths(track)
    assert set(names) == {"chr1", "chr2"}
    assert auth == frozenset({"chr1", "chr2"})
    m = refmod.match_reference(
        ref, names, track_lengths=lengths, authoritative_lengths=auth
    )
    assert m.matched == 2 and m.length_mismatches == 2 and not m.ok


def test_copy_and_blast_use_the_placed_sequence_when_the_alignment_has_none(qapp, tmp_path):
    """PAF records carry no sequence. Copy and BLAST still use the run's placed FASTA."""
    from repeatviewer.widgets import light_genome_browser as lgb

    details = tmp_path / "details"
    details.mkdir()
    (details / "_mapped.fa").write_text(
        ">Harbinger_reversed mapped to chr1:101-110 (-)\n"
        "TTGCAATGCA\n"
        ">Other mapped to chr1:1-4 (+)\n"
        "AAAA\n"
    )
    query = _fasta(tmp_path / "queries.fa", {"OnlyInQuery": "ACGTAAAA"})
    panel = lgb.LightGenomeBrowserPanel()
    panel.path = tmp_path / "reads.paf"
    panel.canvas.chrom = "chr1"
    panel._run = {
        "_folder": str(tmp_path),
        "details_dir": "details",
        "details": {"mapped": "_mapped.fa"},
        "queries": [str(query)],
        "contigs": [{
            "id": "chr1",
            "reads": [{
                "name": "FromRow", "start": 20, "end": 24, "strand": "+",
                "row": "A-C?", "ins": [[21, "TT"]],
            }],
        }],
    }
    placed = gt.TrackItem(100, 110, "-", "Harbinger", seq="")
    menu = panel._read_name_menu(placed)
    by_text = {a.text(): a for a in menu.actions()}
    assert by_text["Copy sequence (FASTA)"].isEnabled()
    assert by_text["BLAST sequence…"].isEnabled()
    assert panel._read_bases(placed) == "TTGCAATGCA"

    # No placed-FASTA record: the original query, turned onto the reference strand.
    minus = gt.TrackItem(0, 8, "-", "OnlyInQuery", seq="")
    assert panel._read_bases(minus) == "TTTTACGT"

    # No files for this name: the alignment row, with an insertion put back.
    row_item = gt.TrackItem(20, 24, "+", "FromRow", seq="")
    assert panel._read_bases(row_item) == "ATTC"
    panel.shutdown()


def test_preview_routes_genome_files_to_light_browser():
    from repeatviewer.widgets.preview_tabs import is_previewable_file

    for name in ("a.bam", "a.sam", "a.gff3", "a.gtf", "a.gff.gz", "a.cram"):
        assert is_previewable_file(Path(name)), name
