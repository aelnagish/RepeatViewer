"""Restriction tails: the extra 5' bases, and what can go wrong with them.

A restriction tail is ``padding + recognition site`` written in lower case in
front of the binding part. Two things decide whether it works at the bench:

* the enzyme needs DNA to hold on to, so there must be enough bases 5' of the
  site. NEB's rule of thumb is six bases, not palindromic and not forming a
  dimer (some enzymes cut with fewer, but an enzyme that was not tested close
  to the end may not cut at all);
* the enzyme must not cut inside the fragment it is meant to open, otherwise
  the insert falls apart in the same tube.

``pick_padding`` chooses the bases, ``tail_for`` builds the tail and
``tail_warnings`` lists, in plain words, what is wrong with a primer's tail.
Nothing here blocks anything: the warnings are shown in amber and the user
decides.

Pure Python, no Qt; unit-tested in ``tests/test_primers_tails.py``.

Developer: Abdullah Yusuf
"""

from __future__ import annotations

from dataclasses import dataclass

from repeatviewer.backends.primers import thermo as th

# NEB: six bases between the 5' end and the site cut well for enzymes that were
# not measured close to the end.
PAD_BASES = 6
# Pads to try, in order: none is a palindrome, none makes a dimer with itself.
PAD_POOL = ("GTCAGT", "AGTCAG", "TGACTG", "CAGTCA", "ACTGAC", "TCAGTG")

# Enzymes offered as a tail (the cloning work-horses plus the Type IIS pair).
TAIL_ENZYMES = (
    "EcoRI", "BamHI", "HindIII", "KpnI", "SacI", "XhoI", "NotI", "XbaI", "SalI", "PstI",
    "NcoI", "NdeI", "SpeI", "SmaI", "EcoRV", "BglII", "NheI", "ApaI", "ClaI", "SphI",
    "MluI", "AscI", "PacI", "AgeI", "AvrII", "MfeI", "NsiI", "SacII", "XmaI", "AflII",
    "HpaI", "PvuII", "StuI", "AatII", "BsaI", "BsmBI",
)

_PLAIN = set("ACGT")


def _enzyme(name: str):
    from Bio import Restriction

    return getattr(Restriction, str(name), None)


def site_of(name: str) -> str:
    """Recognition sequence, upper case ("" when the enzyme is unknown)."""
    enz = _enzyme(name)
    return str(enz.site).upper() if enz is not None else ""


def is_palindrome(seq: str) -> bool:
    s = str(seq or "").upper()
    return bool(s) and s == th.reverse_complement(s)


def _occurrences(needle: str, hay: str) -> int:
    """How often a plain site occurs in ``hay``, both strands (0 if degenerate)."""
    needle = str(needle or "").upper()
    hay = str(hay or "").upper()
    if not needle or not hay or not set(needle) <= _PLAIN:
        return 0
    wanted = {needle, th.reverse_complement(needle)}
    n = 0
    for w in wanted:
        start = hay.find(w)
        while start >= 0:
            n += 1
            start = hay.find(w, start + 1)
    return n


def cuts_in(name: str, sequence: str) -> list[int]:
    """1-based cut positions of the enzyme on a sequence (top strand)."""
    enz = _enzyme(name)
    if enz is None or not sequence:
        return []
    from Bio.Seq import Seq

    try:
        return [int(p) for p in enz.search(Seq(str(sequence).upper()), linear=True)]
    except Exception:  # noqa: BLE001 - enzymes without a defined cut
        return []


def pick_padding(enzyme: str, binding: str = "", *, bases: int = PAD_BASES) -> str:
    """Padding bases that are not a palindrome and make no second site.

    Tries the pool in order against ``padding + site + the first bases that
    bind``; falls back to the first pad when every one of them clashes.
    """
    if bases <= 0:
        return ""
    site = site_of(enzyme)
    head = str(binding or "").upper()[:8]
    for pad in PAD_POOL:
        pad = pad[:bases] if bases <= len(pad) else (pad * 3)[:bases]
        if is_palindrome(pad):
            continue
        if _occurrences(site, pad + site + head) > 1:
            continue
        return pad
    return PAD_POOL[0][:bases]


def tail_for(enzyme: str, padding: str | None = None, binding: str = "") -> str:
    """The lower-case 5' tail for this enzyme (padding chosen when not given)."""
    site = site_of(enzyme)
    if not site:
        return ""
    pad = pick_padding(enzyme, binding) if padding is None else str(padding).upper()
    return (pad + site).lower()


@dataclass
class TailInfo:
    """What the lower-case tail of a primer holds."""

    tail: str = ""
    enzyme: str = ""
    padding: str = ""
    site: str = ""

    @property
    def has_site(self) -> bool:
        return bool(self.enzyme)


def read_tail(primer_seq: str, candidates=TAIL_ENZYMES) -> TailInfo:
    """Find which enzyme's site sits in a primer's tail, and the bases before it."""
    tail, _bind = th.split_tail(primer_seq)
    info = TailInfo(tail=tail)
    if not tail:
        return info
    up = tail.upper()
    best = ""
    for name in candidates:
        site = site_of(name)
        if site and set(site) <= _PLAIN and site in up:
            # the longest site wins (GGTACC before GGTAC-like short ones)
            if len(site) > len(site_of(best)) if best else True:
                best = name
    if not best:
        return info
    site = site_of(best)
    info.enzyme, info.site = best, site
    info.padding = up[:up.index(site)]
    return info


def tail_warnings(primer_seq: str, *, insert: str = "", insert_name: str = "product",
                  partner_seq: str = "", partner_name: str = "") -> list[str]:
    """Plain-word warnings about a primer's restriction tail (empty: nothing wrong)."""
    info = read_tail(primer_seq)
    if not info.has_site:
        return []
    out: list[str] = []
    pad = info.padding
    if len(pad) < PAD_BASES:
        out.append(f"Only {len(pad)} base{'s' if len(pad) != 1 else ''} before the "
                   f"{info.enzyme} site; 6 cut reliably")
    if pad and is_palindrome(pad):
        out.append(f"The {len(pad)} bases before the site read the same on both strands")
    bind = th.binding_part(primer_seq)
    whole = info.tail.upper() + bind[:10]
    if _occurrences(info.site, whole) > 1:
        out.append(f"The tail holds the {info.enzyme} site twice")
    if insert:
        inside = cuts_in(info.enzyme, insert)
        if inside:
            where = ", ".join(f"{c:,}" for c in inside[:3])
            more = "…" if len(inside) > 3 else ""
            out.append(f"{info.enzyme} also cuts the {insert_name} at {where}{more}: "
                       f"it is cut as well")
    if partner_seq:
        other = read_tail(partner_seq)
        if other.enzyme and other.enzyme == info.enzyme:
            who = f" ({partner_name})" if partner_name else ""
            out.append(f"The other primer{who} carries {info.enzyme} too: the insert can go "
                       "in either direction")
    return out


def free_enzymes(insert: str, candidates=TAIL_ENZYMES, limit: int = 8) -> list[str]:
    """Enzymes of the tail list that do not cut the insert."""
    out = []
    for name in candidates:
        if not cuts_in(name, insert):
            out.append(name)
            if len(out) >= limit:
                break
    return out
