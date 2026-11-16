"""Restriction tails: padding, second sites, enzymes cutting the insert."""

from repeatviewer.backends.primers import tails as tl

BIND = "ATGCTAGCTAGGTTCAGATCC"


def test_site_padding_and_build():
    assert tl.site_of("EcoRI") == "GAATTC" and tl.site_of("Nope") == ""
    pad = tl.pick_padding("EcoRI", BIND)
    assert len(pad) == tl.PAD_BASES and not tl.is_palindrome(pad)
    assert tl.tail_for("EcoRI", binding=BIND) == (pad + "GAATTC").lower()
    assert tl.tail_for("EcoRI", "AT") == "atgaattc"
    assert tl.tail_for("Nope") == ""
    assert tl.pick_padding("EcoRI", BIND, bases=0) == ""


def test_padding_avoids_making_a_second_site():
    # A primer starting with AATTC would turn a pad ending in G into a second site.
    pad = tl.pick_padding("EcoRI", "AATTCGGGGG")
    assert (pad + "GAATTC" + "AATTCGGGGG").count("GAATTC") == 1


def test_read_tail():
    info = tl.read_tail(tl.tail_for("KpnI", binding=BIND) + BIND)
    assert (info.enzyme, info.site, len(info.padding)) == ("KpnI", "GGTACC", 6)
    assert not tl.read_tail(BIND).has_site          # no tail at all
    assert not tl.read_tail("atgc" + BIND).has_site  # tail without a site


def test_warnings():
    good = tl.tail_for("EcoRI", binding=BIND) + BIND
    assert tl.tail_warnings(good) == []
    assert tl.tail_warnings(BIND) == []

    short = "atgaattc" + BIND
    words = " | ".join(tl.tail_warnings(short))
    assert "Only 2 bases" in words and "both strands" in words

    insert = "AAAA" + "GAATTC" + "TTTT"
    words = " | ".join(tl.tail_warnings(good, insert=insert))
    assert "EcoRI also cuts the product at 6" in words and "cut as well" in words
    words = " | ".join(tl.tail_warnings(good, insert=insert, insert_name="open sequence"))
    assert "cuts the open sequence" in words

    partner = tl.tail_for("EcoRI", binding=BIND) + "TTGACCTAGCTAGCATGCA"
    words = " | ".join(tl.tail_warnings(good, partner_seq=partner, partner_name="A_R"))
    assert "other primer (A_R) carries EcoRI" in words
    other = tl.tail_for("BamHI", binding=BIND) + "TTGACCTAGCTAGCATGCA"
    assert tl.tail_warnings(good, partner_seq=other) == []

    twice = ("gtcagtgaattc").replace("gtcagt", "gaattcgt") + BIND
    assert any("twice" in w for w in tl.tail_warnings(twice))


def test_free_enzymes():
    insert = "AAAA" + "GAATTC" + "TTTT" + "GGATCC" + "AAAA"
    free = tl.free_enzymes(insert)
    assert "EcoRI" not in free and "BamHI" not in free and len(free) == 8
    assert tl.cuts_in("EcoRI", insert) == [6]
