"""Basic cloning of a PCR product into a vector: TA, blunt and restriction.

What this does is the bookkeeping around a ligation, not a full cloning suite:
which ends the insert has, whether the chosen enzymes really open the insert
cleanly, how long the clone is in each orientation, what a colony PCR gives,
and how far a sequencing primer reaches into the insert.

Sizes that depend on the vector's own sequence (colony PCR, how much vector a
sequencing read sees) come from one number the user reads off their vector's
map or manual: the product the primer pair gives on the **empty** vector. The
clone's product is that number plus the insert. No vector sequence is shipped
with RepeatViewer, so nothing here is guessed from a half-remembered map.

Tandem arrays are the reason for the last warning: direct repeats recombine in
E. coli and come back shorter, so a satellite array needs a recombination-poor
strain and a slow, cool culture.

Pure Python, no Qt; unit-tested in ``tests/test_primers_cloning.py``.

Developer: Abdullah Yusuf
"""

from __future__ import annotations

from dataclasses import dataclass, field

from repeatviewer.backends.primers import tails as tl
from repeatviewer.backends.primers import thermo as th

METHOD_TA = "ta"
METHOD_BLUNT = "blunt"
METHOD_RE = "restriction"
METHOD_LABELS = {METHOD_TA: "TA (A overhangs)", METHOD_BLUNT: "Blunt",
                 METHOD_RE: "Restriction and ligation"}
METHODS = (METHOD_TA, METHOD_BLUNT, METHOD_RE)

# Polymerases that leave the 3' A a TA ligation needs.
A_TAILING = ("Taq", "Long-range")
BLUNT_POLYMERASES = ("Phusion", "Q5", "Pfu")
# Above this, a tandem array is likely to be lost in an ordinary strain.
ARRAY_COPIES_WARN = 3


@dataclass
class SeqPrimer:
    """A primer that reads out of the vector into the insert."""

    name: str
    note: str = ""


@dataclass
class Vector:
    """A cloning vector, as far as the bookkeeping needs it."""

    name: str
    length: int = 0
    methods: tuple[str, ...] = METHODS
    primers: tuple[SeqPrimer, ...] = ()
    note: str = ""

    def takes(self, method: str) -> bool:
        return method in self.methods


PRESETS: tuple[Vector, ...] = (
    Vector("pGEM-T Easy", 3015, (METHOD_TA,),
           (SeqPrimer("M13 forward (-20)"), SeqPrimer("M13 reverse"),
            SeqPrimer("T7 promoter"), SeqPrimer("SP6 promoter")),
           "Comes linear with 3′ T overhangs: the insert needs a 3′ A, so amplify with Taq "
           "or add A to a blunt product."),
    Vector("pJET1.2/blunt", 2974, (METHOD_BLUNT,),
           (SeqPrimer("pJET1.2 forward"), SeqPrimer("pJET1.2 reverse")),
           "Blunt ends. A product made with a proofreading polymerase goes straight in; a "
           "Taq product has to be blunted first."),
    Vector("Your own vector", 0, METHODS, (SeqPrimer("Forward"), SeqPrimer("Reverse")),
           "Give the length and the empty-vector product of your primer pair, from its map."),
)


def preset(name: str) -> Vector:
    for v in PRESETS:
        if v.name == name:
            return v
    return PRESETS[-1]


@dataclass
class Clone:
    """The ligated plasmid in one orientation."""

    orientation: str           # forward, reverse
    insert_length: int = 0
    total_length: int = 0
    colony_product: int = 0    # 0 when the empty-vector product is unknown

    def label(self) -> str:
        way = "insert as drawn" if self.orientation == "forward" else "insert turned around"
        return f"{self.orientation.title()} ({way})"


@dataclass
class CloneResult:
    insert_length: int = 0
    ends: str = ""                      # what the insert's ends look like
    clones: list[Clone] = field(default_factory=list)
    empty_product: int = 0
    warnings: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return bool(self.clones) and not self.warnings


def insert_from_primers(product: str, forward_tail: str = "", reverse_tail: str = "") -> str:
    """The product as it comes out of the tube: the tails are part of it."""
    return (str(forward_tail or "").upper() + str(product or "").upper()
            + th.reverse_complement(str(reverse_tail or "").upper()))


def ends_text(method: str, insert: str, *, enzymes: tuple[str, str] = ("", ""),
              polymerase: str = "") -> str:
    if method == METHOD_TA:
        return "3′ A overhangs (Taq)"
    if method == METHOD_BLUNT:
        return "blunt"
    left, right = enzymes
    return f"{left or '?'} on the left, {right or '?'} on the right"


def _array_warning(insert: str, features=()) -> str:
    from repeatviewer.backends.primers import array_digest as ad

    info = ad.array_info(insert, features)
    if info.period and info.copies >= ARRAY_COPIES_WARN:
        return (f"The insert is a tandem array ({info.copies} copies of {info.period:,} bp): "
                "direct repeats recombine in E. coli and clones often come back shorter. "
                "Use a recombination-poor strain (recA, sbcC) and grow at 30 °C, then check "
                "the length of every clone")
    return ""


def plan(insert: str, vector: Vector, method: str, *, enzymes: tuple[str, str] = ("", ""),
         polymerase: str = "", empty_product: int = 0, features=(),
         forward_primer: str = "", reverse_primer: str = "") -> CloneResult:
    """Work out both orientations, the colony PCR sizes and what to watch out for."""
    insert = str(insert or "").upper()
    res = CloneResult(insert_length=len(insert), empty_product=int(empty_product or 0))
    if not insert:
        res.warnings.append("No insert sequence")
        return res
    res.ends = ends_text(method, insert, enzymes=enzymes, polymerase=polymerase)
    if not vector.takes(method):
        res.warnings.append(f"{vector.name} does not take a {METHOD_LABELS[method].lower()} "
                            "insert")
    if method == METHOD_TA:
        if polymerase and not any(polymerase.startswith(p) for p in A_TAILING):
            res.warnings.append(f"{polymerase} leaves blunt ends: add A to the product first, "
                                "or use a blunt vector")
        res.notes.append("The 3′ A comes from Taq, so run the PCR with Taq or A-tail the "
                         "product for 20 min at 72 °C")
    elif method == METHOD_BLUNT:
        if polymerase and any(polymerase.startswith(p) for p in A_TAILING):
            res.warnings.append(f"{polymerase} leaves a 3′ A: blunt the product first "
                                "(the kit's blunting step), or use a TA vector")
    else:
        left, right = enzymes
        if not left or not right:
            res.warnings.append("Choose the enzyme at each end of the insert")
        for name, where in ((left, "left"), (right, "right")):
            if not name:
                continue
            inside = tl.cuts_in(name, insert)
            # the site in the tail itself is not "inside"
            inside = [c for c in inside if 20 < c < len(insert) - 20]
            if inside:
                at = ", ".join(f"{c:,}" for c in inside[:3])
                res.warnings.append(f"{name} ({where} end) also cuts inside the insert at {at}: "
                                    "the insert falls apart in the same tube")
        if left and right and left == right:
            res.notes.append(f"Both ends are {left}: the insert goes in either way, and the "
                             "vector can close on itself (dephosphorylate it)")
    for primer_seq, side in ((forward_primer, "forward"), (reverse_primer, "reverse")):
        if not primer_seq:
            continue
        info = tl.read_tail(primer_seq)
        if method == METHOD_RE and not info.has_site:
            res.warnings.append(f"The {side} primer carries no restriction site in its tail")
    array = _array_warning(insert, features)
    if array:
        res.warnings.append(array)
    total = (vector.length + len(insert)) if vector.length else 0
    for orientation in ("forward", "reverse"):
        res.clones.append(Clone(orientation=orientation, insert_length=len(insert),
                                total_length=total,
                                colony_product=(res.empty_product + len(insert))
                                if res.empty_product else 0))
    if not res.empty_product:
        res.notes.append("Give the product your primer pair makes on the empty vector (its "
                         "map or manual says so) and the colony PCR sizes appear")
    if not vector.length:
        res.notes.append("Give the vector's length for the size of the finished plasmid")
    return res


def reads_into_insert(empty_product: int, read_length: int = 700) -> int:
    """How far a sequencing read reaches into the insert, give or take.

    Half the empty-vector product sits between the primer and the insert, so a
    read of ``read_length`` good bases covers about that much of the insert.
    """
    flank = max(0, int(empty_product)) // 2
    return max(0, int(read_length) - flank)


def colony_gel(result: CloneResult, *, title: str = ""):
    """A gel of a colony PCR: empty vector beside a clone with the insert."""
    from repeatviewer.backends.primers import gel_doc as gd

    empty = result.empty_product
    clone = result.clones[0].colony_product if result.clones else 0
    if not empty or not clone:
        return None
    sizes = [empty, clone]
    ladder = gd.pick_ladder_name(sizes)
    doc = gd.GelDoc(ladder=ladder, agarose=gd.pick_agarose(sizes), title=title)
    doc.lanes = [gd.ladder_lane(ladder),
                 gd.GelLane(label="No insert", bands=[gd.GelBand(size=empty, label="empty")],
                            note="colony that took up re-closed vector"),
                 gd.GelLane(label="With insert",
                            bands=[gd.GelBand(size=clone, label="insert")],
                            note=f"empty {empty:,} bp + insert {result.insert_length:,} bp")]
    return doc


def summary(result: CloneResult, vector: Vector, method: str) -> tuple[str, str]:
    """(tone, sentence) for the window."""
    if result.warnings:
        return "amber", result.warnings[0] + "."
    if not result.clones:
        return "amber", "Nothing to ligate yet."
    clone = result.clones[0]
    bits = [f"{result.insert_length:,} bp insert, {result.ends}, into {vector.name}"]
    if clone.total_length:
        bits.append(f"plasmid {clone.total_length:,} bp")
    if clone.colony_product:
        bits.append(f"colony PCR {clone.colony_product:,} bp against {result.empty_product:,} bp "
                    "without insert")
    return "green", "; ".join(bits) + "."


def write_report(path, insert_name: str, result: CloneResult, vector: Vector,
                 method: str) -> "object":
    import html as _html
    from pathlib import Path

    tone, line = summary(result, vector, method)
    colour = "#14532D" if tone == "green" else "#8a4b00"
    rows = "".join(
        f"<tr><td>{_html.escape(c.label())}</td><td>{c.insert_length:,} bp</td>"
        f"<td>{c.total_length:,} bp</td>"
        f"<td>{c.colony_product:,} bp</td></tr>" if c.total_length else
        f"<tr><td>{_html.escape(c.label())}</td><td>{c.insert_length:,} bp</td>"
        f"<td>not known</td><td>{c.colony_product:,} bp</td></tr>"
        for c in result.clones)
    items = "".join(f"<li>{_html.escape(w)}.</li>" for w in result.warnings + result.notes)
    doc = f"""<!DOCTYPE html>
<html><head><meta charset="utf-8"><title>Cloning: {_html.escape(insert_name)}</title>
<style>body{{font-family:-apple-system,Segoe UI,sans-serif;margin:2em;max-width:60em}}
table{{border-collapse:collapse;width:100%}}td,th{{border-bottom:1px solid #ddd;padding:6px;
text-align:left}}</style></head><body>
<h1>Cloning: {_html.escape(insert_name)}</h1>
<p style="color:{colour}"><b>{_html.escape(line)}</b></p>
<p>{_html.escape(METHOD_LABELS.get(method, method))} into {_html.escape(vector.name)}.
{_html.escape(vector.note)}</p>
<table><tr><th>Orientation</th><th>Insert</th><th>Plasmid</th><th>Colony PCR</th></tr>
{rows}</table>
<ul>{items}</ul>
<p style="color:gray">Sequencing primers: {_html.escape(', '.join(p.name for p in vector.primers))}.
Made by RepeatViewer.</p></body></html>
"""
    out = Path(path)
    out.write_text(doc, encoding="utf-8")
    return out
