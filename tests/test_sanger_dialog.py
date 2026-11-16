"""Sanger Check window: verdict, differences, chromatogram, report."""

import random

from repeatviewer.backends.primers import thermo as th
from tests.test_primers_sanger import write_ab1


def _insert(n=400, seed=11):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


def _ab1(path, seq, quality, with_trace=True):
    locs = [12 * i + 8 for i in range(len(seq))]
    chans = None
    if with_trace:
        chans = {b: [0] * (len(seq) * 12 + 40) for b in "GATC"}
        for i, b in enumerate(seq):
            for k in range(-4, 5):
                chans[b][locs[i] + k] = 800 - abs(k) * 100
    return write_ab1(path, seq, quality, chans, locs)


def test_window_lists_differences_and_draws_the_trace(qapp, tmp_path):
    from repeatviewer.views.sanger_dialog import SangerDialog

    seq = _insert()
    read = list(seq[20:260])
    read[100] = "A" if read[100] != "A" else "G"
    fwd = _ab1(tmp_path / "cloneA_T7.ab1", "".join(read), [40] * len(read))
    rev = _ab1(tmp_path / "cloneA_SP6.ab1", th.reverse_complement(seq[150:]), [40] * 250)

    dlg = SangerDialog(seq, expected_name="pGEM insert", folder=tmp_path)
    dlg.add_reads([fwd, rev])
    assert "1 difference" in dlg.verdict.text() and "121" in dlg.verdict.text()
    texts = [dlg.table.item(r, 2).text() for r in range(dlg.table.rowCount())]
    assert texts[0].startswith("121:") and texts[1] == "none"
    assert dlg.table.item(0, 3).text() == "Q40"
    assert "forward" in dlg.table.item(0, 1).text()
    assert "reverse" in dlg.table.item(1, 1).text()

    dlg.table.selectRow(0)
    assert dlg.trace._centre == 100 and dlg.trace._expected == seq[120]
    dlg.trace.resize(600, 180)
    dlg.trace.grab()                      # paints without raising
    dlg.deleteLater()


def test_window_reports_gaps_and_saves_a_report(qapp, tmp_path):
    from repeatviewer.views.sanger_dialog import SangerDialog

    seq = _insert(seed=3)
    short = _ab1(tmp_path / "cloneB.ab1", seq[:150], [40] * 150)
    dlg = SangerDialog(seq, expected_name="pGEM insert", folder=tmp_path)
    dlg.add_reads([short])
    assert "not read" in dlg.verdict.text() and "151" in dlg.verdict.text()
    dlg.report_button.click()
    out = tmp_path / "pGEM_insert_sanger.html"
    assert out.is_file() and "cloneB" in out.read_text()
    dlg.deleteLater()


def test_window_reads_a_fasta_as_the_expected_sequence(qapp, tmp_path):
    from repeatviewer.views.sanger_dialog import SangerDialog, _read_fasta

    seq = _insert(seed=5)
    fa = tmp_path / "insert.fa"
    fa.write_text(f">insert_1 first\n{seq[:200]}\n>insert_2\n{seq[200:]}\n")
    records = _read_fasta(fa)
    assert [n for n, _s in records] == ["insert_1", "insert_2"]

    dlg = SangerDialog(folder=tmp_path)
    assert "No expected sequence" in dlg.expected_label.text()
    dlg.set_expected(records[0][1], records[0][0])
    read = _ab1(tmp_path / "c.ab1", seq[10:190], [40] * 180)
    dlg.add_reads([read])
    assert "insert_1" in dlg.expected_label.text()
    assert dlg.table.rowCount() == 1 and dlg.table.item(0, 2).text() == "none"
    dlg.deleteLater()


def test_menu_offers_sanger_check(qapp):
    from PySide6.QtWidgets import QMenu

    from repeatviewer.pcr_actions import PcrActionsMixin
    from repeatviewer.primer_actions import PrimerActionsMixin

    class _Host(PrimerActionsMixin, PcrActionsMixin):
        pass

    menu = QMenu()
    _Host()._populate_primers_menu(menu)
    texts = [a.text() for a in menu.actions() if not a.isSeparator()]
    assert texts[:5] == ["Design", "Map", "Database", "In-silico PCR", "Sanger Check"]
    assert all("…" not in t for t in texts)
