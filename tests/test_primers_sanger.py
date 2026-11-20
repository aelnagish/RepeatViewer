"""Sanger verification of a clone: trimming, differences, coverage, verdict."""

import struct

from repeatviewer.backends.primers import sanger as sg
from repeatviewer.backends.primers import thermo as th

INSERT = ("ACGTTGCATGCATTGACCATGGCTAGCTAGCATCGATCGGGATCCTTTTAAGGCCAATCGATCGATCGTTAGCAT"
          "GGCCTTAAGGCCATCGATTTACCGGTTAACCGGTTAAGGCCTTAACCGGAATTCCGGTTAACCGGTTAACCGGTT")

_HEAD = ">H4sI2H3I"
_DIR = ">4sI2H4I"


def write_ab1(path, bases, quality, channels=None, locations=None):
    """Smallest .ab1 Biopython will read: PBAS2, PCON2, PLOC2 and the four traces."""
    tags = [(b"PBAS", 2, 2, 1, bases.encode()),
            (b"PCON", 2, 2, 1, bytes(min(255, max(0, q)) for q in quality)),
            (b"FWO_", 1, 2, 1, b"GATC")]
    if locations:
        tags.append((b"PLOC", 2, 4, 2, struct.pack(f">{len(locations)}h", *locations)))
    for n, base in enumerate("GATC"):
        data = (channels or {}).get(base)
        if data:
            tags.append((f"DATA".encode(), 9 + n, 4, 2, struct.pack(f">{len(data)}h", *data)))
    blob, offsets = b"", []
    for _name, _num, _code, _size, data in tags:
        if len(data) <= 4:          # four bytes or fewer live inside the entry
            offsets.append(None)
            continue
        offsets.append(len(blob))
        blob += data + b"\0" * ((4 - len(data) % 4) % 4)
    head_len = 4 + struct.calcsize(_HEAD)
    data_start = head_len
    dir_offset = data_start + len(blob)
    entries = b""
    for (name, num, code, size, data), off in zip(tags, offsets):
        n_elem = len(data) // size
        where = (int.from_bytes(data.ljust(4, b"\0"), "big") if off is None
                 else data_start + off)
        entries += struct.pack(_DIR, name, num, code, size, n_elem, len(data), where, 0)
    head = b"ABIF" + struct.pack(_HEAD, 101, b"tdir", 1, 1023, 28, len(tags),
                                 len(entries), dir_offset)
    path.write_bytes(head + blob + entries)
    return path


def _trace(seq, q=40):
    return sg.Trace(name="r", bases=seq, quality=[q] * len(seq))


def test_mott_trim_keeps_the_good_middle():
    qual = [5] * 20 + [45] * 200 + [4] * 30
    start, end = sg.mott_trim(qual)
    assert 15 <= start <= 21 and 219 <= end <= 225
    assert sg.mott_trim([40] * 10) == (0, 10)          # shorter than the minimum: keep all
    assert sg.mott_trim([]) == (0, 0)


def test_read_matches_both_strands():
    fwd = sg.check_read(INSERT, _trace(INSERT[20:120]))
    assert (fwd.strand, fwd.start, fwd.end) == ("+", 21, 120)
    assert fwd.differences == [] and fwd.identity == 1.0
    rev = sg.check_read(INSERT, _trace(th.reverse_complement(INSERT[20:120])))
    assert (rev.strand, rev.start, rev.end) == ("-", 21, 120)
    assert rev.differences == []
    assert "no differences" in rev.summary() and "reverse" in rev.summary()


def test_mismatch_insertion_and_deletion_are_listed():
    read = INSERT[10:90]
    swapped = read[:30] + ("A" if read[30] != "A" else "C") + read[31:]
    check = sg.check_read(INSERT, _trace(swapped))
    assert [d.kind for d in check.differences] == ["mismatch"]
    assert check.differences[0].position == 41 and check.differences[0].quality == 40

    missing = sg.check_read(INSERT, _trace(read[:30] + read[31:]))
    assert [d.kind for d in missing.differences] == ["deletion"]
    extra = sg.check_read(INSERT, _trace(read[:30] + "G" + read[30:]))
    assert [d.kind for d in extra.differences] == ["insertion"]


def test_low_quality_difference_is_doubtful():
    read = list(INSERT[10:90])
    read[30] = "A" if read[30] != "A" else "C"
    qual = [40] * 80
    qual[30] = 12
    check = sg.check_read(INSERT, sg.Trace(name="r", bases="".join(read), quality=qual))
    assert check.differences[0].doubtful and check.real_differences == []
    assert sg.verdict(INSERT, [check])[0] == "amber"       # still a gap in coverage


def test_coverage_gaps_and_verdict():
    left = sg.check_read(INSERT, _trace(INSERT[:90]))
    right = sg.check_read(INSERT, _trace(INSERT[80:]))
    assert sg.gaps_in_coverage(len(INSERT), [left, right]) == []
    tone, line = sg.verdict(INSERT, [left, right])
    assert tone == "green" and "Insert confirmed" in line and f"{len(INSERT):,} bp" in line

    short = sg.check_read(INSERT, _trace(INSERT[:60]))
    assert sg.gaps_in_coverage(len(INSERT), [short]) == [(61, len(INSERT))]
    tone, line = sg.verdict(INSERT, [short])
    assert tone == "amber" and "not read" in line
    assert sg.verdict(INSERT, [])[1].startswith("No read")


def test_reads_an_ab1_with_traces(tmp_path):
    bases = INSERT[5:85]
    locs = [10 * i + 5 for i in range(len(bases))]
    channels = {b: [i % 50 for i in range(len(bases) * 10)] for b in "GATC"}
    path = write_ab1(tmp_path / "clone_A.ab1", bases, [38] * len(bases), channels, locs)
    trace = sg.read_trace(path)
    assert trace.bases == bases and trace.quality[:3] == [38, 38, 38]
    assert trace.has_trace and set(trace.channels) == set("GATC")
    assert trace.locations == locs and trace.name == "clone_A"

    checks = sg.check_reads(INSERT, [path, tmp_path / "missing.ab1"])
    assert checks[0].start == 6 and checks[0].differences == []
    assert "Could not read this file" in checks[1].note


def test_plain_read_and_report(tmp_path):
    fa = tmp_path / "clone_B.seq"
    fa.write_text(">clone_B\n" + INSERT[:70] + "\n")
    check = sg.check_read(INSERT, sg.read_trace(fa))
    assert check.name == "clone_B" and check.trace.quality == []
    assert check.differences == [] and "no quality" in check.note

    out = sg.write_report(tmp_path / "r.html", "pGEM-CroSat3", INSERT, [check])
    html = out.read_text()
    assert "Sanger check: pGEM-CroSat3" in html and "clone_B" in html
    assert "not read" in html.lower() and "\u2014" not in html
