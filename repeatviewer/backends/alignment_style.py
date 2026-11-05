"""How alignment-track reads are coloured, grouped and described (genome browser).

Works on every alignment item: BAM records (``genome_tracks.BamItem``, read from
pysam), SAM lines and PAF / BLAST hits (``TrackItem``), reading whatever each one
carries (flags, mate, tags). No Qt; tested in ``tests/test_alignment_style.py``.

Colour modes
------------
* ``strand``: forward blue, reverse rose (the default).
* ``mapq``: amber (MAPQ 0) to blue (60).
* ``identity``: red (80 % or lower) to green (100 %), for PAF / BLAST.
* ``insert``: pairs at the expected distance stay grey; farther apart red,
  closer blue, mate on another sequence a colour per sequence, mate unmapped
  amber; wrong orientation (outward, or both on one strand) teal / green.
* ``first_strand``: strand of the first read of the pair (stranded RNA-seq).
* ``read_group`` / ``sample`` / ``library``: a colour per value (from @RG).
* ``tag``: a colour per value of any tag (``HP``, ``CB``…).
"""

from __future__ import annotations

import hashlib

COLOR_STRAND = "strand"
COLOR_MAPQ = "mapq"
COLOR_IDENTITY = "identity"
COLOR_INSERT = "insert"
COLOR_FIRST_STRAND = "first_strand"
COLOR_READ_GROUP = "read_group"
COLOR_SAMPLE = "sample"
COLOR_LIBRARY = "library"
COLOR_TAG = "tag"
COLOR_NONE = "none"

COLOR_MODES = [
    (COLOR_STRAND, "Strand"),
    (COLOR_FIRST_STRAND, "First-of-pair strand"),
    (COLOR_INSERT, "Insert size and pair orientation"),
    (COLOR_MAPQ, "Mapping quality"),
    (COLOR_IDENTITY, "Identity"),
    (COLOR_READ_GROUP, "Read group"),
    (COLOR_SAMPLE, "Sample"),
    (COLOR_LIBRARY, "Library"),
    (COLOR_TAG, "Tag…"),
    (COLOR_NONE, "No colour"),
]

FWD = "#b9cfe8"
REV = "#e8c3cf"
NEUTRAL = "#c9c9c9"
LARGER = "#e06a5c"        # mates farther apart than expected
SMALLER = "#5c8fe0"       # closer than expected
MATE_UNMAPPED = "#e0b43a"
OUTWARD = "#3fa89a"       # RL: reads point away from each other
SAME_STRAND = "#6aae4f"   # LL / RR
# soft, distinct colours for categories (read groups, samples, tag values…)
CATEGORY = ["#7aa6d6", "#f0a35e", "#7cc08a", "#e5807a", "#9fd3cd", "#e9cf6b",
            "#c59ac2", "#f4b3bd", "#b59a86", "#aab4bf", "#8fc1e8", "#d4a5e0"]

FLAG_PAIRED = 0x1
FLAG_PROPER = 0x2
FLAG_UNMAPPED = 0x4
FLAG_MATE_UNMAPPED = 0x8
FLAG_REVERSE = 0x10
FLAG_MATE_REVERSE = 0x20
FLAG_FIRST = 0x40
FLAG_SECOND = 0x80
FLAG_SECONDARY = 0x100
FLAG_QCFAIL = 0x200
FLAG_DUPLICATE = 0x400
FLAG_SUPPLEMENTARY = 0x800

_FLAG_NAMES = [(FLAG_PAIRED, "paired"), (FLAG_PROPER, "proper pair"),
               (FLAG_MATE_UNMAPPED, "mate unmapped"), (FLAG_REVERSE, "reverse"),
               (FLAG_MATE_REVERSE, "mate reverse"), (FLAG_FIRST, "first in pair"),
               (FLAG_SECOND, "second in pair"), (FLAG_SECONDARY, "secondary"),
               (FLAG_QCFAIL, "QC fail"), (FLAG_DUPLICATE, "duplicate"),
               (FLAG_SUPPLEMENTARY, "supplementary")]


# ---------------------------------------------------------------- accessors


def _seg(it):
    return getattr(it, "seg", None)


def flag(it) -> int | None:
    s = _seg(it)
    if s is not None:
        return int(s.flag)
    f = (getattr(it, "info", None) or {}).get("flag")
    return int(f) if f is not None else None


def tag(it, name: str):
    s = _seg(it)
    if s is not None:
        return s.get_tag(name) if s.has_tag(name) else None
    tags = _info_tags(it)
    raw = tags.get(name)
    if raw is None or not isinstance(raw, str):
        return raw
    return int(raw) if raw.lstrip("-").isdigit() else raw


def _info_tags(it) -> dict:
    info = getattr(it, "info", None) or {}
    if "_tags" in info:
        return info["_tags"] or {}
    raw = info.get("_rawtags")
    if raw:
        from repeatviewer.backends.genome_tracks import bam_tags

        return bam_tags(raw)
    return {}


def all_tags(it) -> list[tuple[str, object]]:
    s = _seg(it)
    if s is not None:
        return [(k, v) for k, v in s.get_tags()]
    return list(_info_tags(it).items())


def mate(it) -> tuple[str | None, int, int]:
    """(mate sequence or None, mate start (0-based) or -1, template length)."""
    s = _seg(it)
    if s is not None:
        name = s.next_reference_name if s.next_reference_id >= 0 else None
        return name, int(s.next_reference_start), int(s.template_length)
    info = getattr(it, "info", None) or {}
    rn = info.get("_rnext")
    return (None if rn in (None, "*") else rn), int(info.get("_pnext", -1)), \
        int(info.get("_tlen", 0))


def insertions(it) -> list[tuple[int, int]]:
    """(reference position, inserted length) of each insertion in the read."""
    s = _seg(it)
    if s is not None:
        from repeatviewer.backends.genome_tracks import _CIGAR_OPS, cigar_extras

        return cigar_extras(it.start, [(_CIGAR_OPS[op], n)
                                       for op, n in (s.cigartuples or [])])[0]
    info = getattr(it, "info", None) or {}
    if "_ins" in info:
        return info["_ins"]
    return [(pos, 1) for pos in (getattr(it, "ins", None) or [])]


def clips(it) -> tuple[int, int]:
    """Soft-clipped (or, for PAF, unaligned) bases at the read's left and right end on
    the reference."""
    s = _seg(it)
    if s is not None:
        cig = s.cigartuples or []
        left = cig[0][1] if cig and cig[0][0] == 4 else 0
        right = cig[-1][1] if len(cig) > 1 and cig[-1][0] == 4 else 0
        return left, right
    return tuple((getattr(it, "info", None) or {}).get("_clip", (0, 0)))


def qualities(it):
    """Base qualities (phred, indexable by query offset), or None."""
    s = _seg(it)
    if s is not None:
        return s.query_qualities
    q = (getattr(it, "info", None) or {}).get("_qual")
    if q is None:
        return None
    if isinstance(q, str):
        return [ord(c) - 33 for c in q]
    return q


def quality_alpha(q: int | None) -> float:
    """Opacity of a mismatch for its base quality, as IGV fades low-quality bases:
    faint below 5, full from 20."""
    if q is None or q >= 20:
        return 1.0
    if q < 5:
        return 0.15
    return 0.15 + 0.85 * (q - 5) / 15.0


def describe_flag(f: int | None) -> str:
    if f is None:
        return ""
    names = [n for bit, n in _FLAG_NAMES if f & bit]
    return f"{f} ({', '.join(names)})" if names else str(f)


# ---------------------------------------------------------------- insert sizes


def insert_bounds(items, low_pct: float = 0.5, high_pct: float = 99.5) -> tuple[int, int] | None:
    """Expected insert size range from proper pairs among ``items`` (percentiles)."""
    sizes = []
    for it in items:
        f = flag(it)
        if f is None or not f & FLAG_PAIRED or not f & FLAG_PROPER:
            continue
        _m, _p, tlen = mate(it)
        if tlen:
            sizes.append(abs(tlen))
        if len(sizes) >= 20000:
            break
    if len(sizes) < 20:
        return None
    sizes.sort()

    def pct(q: float) -> int:
        return sizes[min(len(sizes) - 1, int(len(sizes) * q / 100.0))]

    return pct(low_pct), pct(high_pct)


def pair_orientation(it) -> str | None:
    """'LR' (normal, inward), 'RL' (outward), 'LL' / 'RR' (same strand), or None."""
    f = flag(it)
    if f is None or not f & FLAG_PAIRED or f & FLAG_MATE_UNMAPPED:
        return None
    mchrom, mpos, _tlen = mate(it)
    if mchrom is None or mpos < 0:
        return None
    rev = bool(f & FLAG_REVERSE)
    mrev = bool(f & FLAG_MATE_REVERSE)
    if rev == mrev:
        return "RR" if rev else "LL"
    left_is_me = it.start <= mpos
    left_rev = rev if left_is_me else mrev
    return "RL" if left_rev else "LR"


# ---------------------------------------------------------------- colours


def category_color(value) -> str:
    if value is None or value == "":
        return NEUTRAL
    digest = hashlib.md5(str(value).encode()).digest()
    return CATEGORY[digest[0] % len(CATEGORY)]


def _mix(a: tuple[int, int, int], b: tuple[int, int, int], t: float) -> str:
    t = max(0.0, min(1.0, t))
    r, g, b_ = (round(x + (y - x) * t) for x, y in zip(a, b))
    return f"#{r:02x}{g:02x}{b_:02x}"


def group_colors(read_groups: dict | None, mode: str) -> dict[str, str]:
    """Read group / sample / library -> colour, in header order, so the first eight
    groups never share a colour (hashing is used for values not in the header)."""
    groups = read_groups or {}
    if mode == COLOR_READ_GROUP:
        values = list(groups)
    elif mode in (COLOR_SAMPLE, COLOR_LIBRARY):
        k = "SM" if mode == COLOR_SAMPLE else "LB"
        values = list(dict.fromkeys(str(g.get(k, rid)) for rid, g in groups.items()))
    else:
        return {}
    return {v: CATEGORY[n % len(CATEGORY)] for n, v in enumerate(values)}


def read_color(it, mode: str, *, chrom: str = "", tag_name: str = "",
               insert_range: tuple[int, int] | None = None,
               read_groups: dict | None = None,
               palette: dict[str, str] | None = None) -> str:
    """Fill colour (hex) of one read for ``mode``. ``palette`` (from
    ``group_colors``) fixes the colours of read groups, samples or libraries."""
    if mode == COLOR_NONE:
        return NEUTRAL
    if mode == COLOR_MAPQ:
        mq = getattr(it, "mapq", None)
        if mq is None or mq < 0:
            return NEUTRAL
        return _mix((224, 180, 58), (99, 146, 214), mq / 60.0)
    if mode == COLOR_IDENTITY:
        score = getattr(it, "score", None)
        if score is None:
            return NEUTRAL
        return _mix((208, 87, 74), (95, 174, 102), (score - 0.8) / 0.2)
    if mode == COLOR_INSERT:
        f = flag(it)
        if f is None or not f & FLAG_PAIRED:
            return NEUTRAL
        if f & FLAG_MATE_UNMAPPED:
            return MATE_UNMAPPED
        mchrom, _mpos, tlen = mate(it)
        if mchrom is not None and chrom and mchrom != chrom:
            return category_color(mchrom)
        orient = pair_orientation(it)
        if orient == "RL":
            return OUTWARD
        if orient in ("LL", "RR"):
            return SAME_STRAND
        if insert_range and tlen:
            lo, hi = insert_range
            if abs(tlen) > hi:
                return LARGER
            if abs(tlen) < lo:
                return SMALLER
        return NEUTRAL
    if mode == COLOR_FIRST_STRAND:
        f = flag(it)
        rev = it.strand == "-"
        if f is not None and f & FLAG_PAIRED and f & FLAG_SECOND:
            rev = not rev
        return REV if rev else FWD
    if mode in (COLOR_READ_GROUP, COLOR_SAMPLE, COLOR_LIBRARY):
        rg = tag(it, "RG")
        if rg is None:
            return NEUTRAL
        value = str(rg)
        if mode != COLOR_READ_GROUP:
            info = (read_groups or {}).get(value, {})
            value = str(info.get("SM" if mode == COLOR_SAMPLE else "LB", rg))
        if palette and value in palette:
            return palette[value]
        return category_color(value)
    if mode == COLOR_TAG:
        return category_color(tag(it, tag_name) if tag_name else None)
    return REV if getattr(it, "strand", "+") == "-" else FWD


def legend(mode: str, *, tag_name: str = "") -> list[tuple[str, str]]:
    """(colour, meaning) rows for the track's colour mode."""
    if mode == COLOR_INSERT:
        return [(NEUTRAL, "expected"), (LARGER, "larger insert"), (SMALLER, "smaller insert"),
                (OUTWARD, "outward (RL)"), (SAME_STRAND, "same strand (LL/RR)"),
                (MATE_UNMAPPED, "mate unmapped"), ("", "mate elsewhere: a colour per sequence")]
    if mode == COLOR_MAPQ:
        return [(_mix((224, 180, 58), (99, 146, 214), 0), "MAPQ 0"),
                (_mix((224, 180, 58), (99, 146, 214), 1), "MAPQ 60")]
    if mode == COLOR_IDENTITY:
        return [(_mix((208, 87, 74), (95, 174, 102), 0), "≤ 80 %"),
                (_mix((208, 87, 74), (95, 174, 102), 1), "100 %")]
    if mode in (COLOR_STRAND, COLOR_FIRST_STRAND):
        return [(FWD, "forward"), (REV, "reverse")]
    return []


def mode_label(mode: str, tag_name: str = "") -> str:
    if mode == COLOR_TAG:
        return f"tag {tag_name}" if tag_name else "tag"
    return dict(COLOR_MODES).get(mode, mode).lower()


# ---------------------------------------------------------------- read details


def details(it, *, chrom: str = "") -> list[tuple[str, str]]:
    """(label, value) rows for the read (or feature) popup."""
    if getattr(it, "kind", "alignment") not in ("alignment", "") and flag(it) is None:
        return _feature_details(it, chrom=chrom)
    rows: list[tuple[str, str]] = [("Read", it.name or "(unnamed)"),
                                   ("Position", (f"{chrom + ':' if chrom else ''}"
                                                 f"{it.start + 1:,}–{it.end:,} ({it.strand})")),
                                   ("Length on reference", f"{it.end - it.start:,} bp")]
    s = _seg(it)
    info = getattr(it, "info", None) or {}
    cigar = s.cigarstring if s is not None else info.get("cigar")
    if cigar:
        rows.append(("CIGAR", cigar if len(cigar) <= 80 else cigar[:77] + "…"))
    if s is not None and s.query_length:
        rows.append(("Read length", f"{s.query_length:,} bp"))
    left, right = clips(it)
    if left or right:
        rows.append(("Clipped", f"{left:,} bp left · {right:,} bp right"))
    big = [n for _p, n in insertions(it) if n >= 2]
    if big:
        rows.append(("Insertions", f"{len(insertions(it)):,} (largest {max(big):,} bp)"))
    mq = getattr(it, "mapq", None)
    if mq is not None:
        rows.append(("MAPQ", str(mq)))
    score = getattr(it, "score", None)
    if score is not None:
        rows.append(("Identity", f"{score:.1%}"))
    f = flag(it)
    if f is not None:
        rows.append(("Flags", describe_flag(f)))
        if f & FLAG_PAIRED:
            mchrom, mpos, tlen = mate(it)
            if f & FLAG_MATE_UNMAPPED:
                rows.append(("Mate", "unmapped"))
            elif mchrom is not None:
                mstrand = "-" if f & FLAG_MATE_REVERSE else "+"
                rows.append(("Mate", f"{mchrom}:{mpos + 1:,} ({mstrand})"))
            if tlen:
                rows.append(("Insert size", f"{abs(tlen):,} bp"))
            orient = pair_orientation(it)
            if orient:
                rows.append(("Pair orientation", orient))
    for k, v in info.items():
        if not k.startswith("_") and k not in ("flag", "cigar"):
            rows.append((k.capitalize(), str(v)))
    for k, v in all_tags(it):
        text = str(v)
        rows.append((k, text if len(text) <= 60 else text[:57] + "…"))
    return rows


def _feature_details(it, *, chrom: str = "") -> list[tuple[str, str]]:
    rows = [("Name", it.name or "(unnamed)"), ("Type", it.kind or "feature"),
            ("Position", f"{chrom + ':' if chrom else ''}{it.start + 1:,}–{it.end:,}"
                         + (f" ({it.strand})" if it.strand in "+-" and it.strand else "")),
            ("Length", f"{it.end - it.start:,} bp")]
    blocks = getattr(it, "blocks", None) or []
    if len(blocks) > 1:
        rows.append(("Exons", str(len(blocks))))
    if getattr(it, "score", None) is not None:
        rows.append(("Score", f"{it.score:g}"))
    for k, v in (getattr(it, "info", None) or {}).items():
        if not str(k).startswith("_"):
            text = str(v)
            rows.append((str(k), text if len(text) <= 80 else text[:77] + "…"))
    return rows


# ---------------------------------------------------------------- group and sort

GROUP_NONE = "none"
GROUP_MODES = [
    (GROUP_NONE, "No groups"),
    (COLOR_STRAND, "Strand"),
    (COLOR_FIRST_STRAND, "First-of-pair strand"),
    ("pair_orientation", "Pair orientation"),
    ("mate_chrom", "Mate sequence"),
    ("chimeric", "Chimeric (SA tag)"),
    ("supplementary", "Primary / secondary / supplementary"),
    (COLOR_MAPQ, "Mapping quality"),
    (COLOR_READ_GROUP, "Read group"),
    (COLOR_SAMPLE, "Sample"),
    (COLOR_LIBRARY, "Library"),
    (COLOR_TAG, "Tag…"),
]


def group_function(mode: str, *, tag_name: str = "", read_groups: dict | None = None,
                   chrom: str = ""):
    """item -> group label (str) for ``mode``, or None for no groups."""
    groups = read_groups or {}
    if mode in (None, "", GROUP_NONE):
        return None
    if mode == COLOR_STRAND:
        return lambda it: "forward" if it.strand != "-" else "reverse"
    if mode == COLOR_FIRST_STRAND:
        def first(it):
            f = flag(it)
            rev = it.strand == "-"
            if f is not None and f & FLAG_PAIRED and f & FLAG_SECOND:
                rev = not rev
            return "first-of-pair reverse" if rev else "first-of-pair forward"
        return first
    if mode == "pair_orientation":
        return lambda it: pair_orientation(it) or "unpaired / mate unmapped"
    if mode == "mate_chrom":
        def mate_seq(it):
            f = flag(it)
            if f is None or not f & FLAG_PAIRED:
                return "unpaired"
            if f & FLAG_MATE_UNMAPPED:
                return "mate unmapped"
            m = mate(it)[0]
            return "same sequence" if m is None or m == chrom else f"mate on {m}"
        return mate_seq
    if mode == "chimeric":
        return lambda it: "chimeric" if tag(it, "SA") is not None else "not chimeric"
    if mode == "supplementary":
        def kind(it):
            f = flag(it) or 0
            if f & FLAG_SUPPLEMENTARY:
                return "supplementary"
            if f & FLAG_SECONDARY or getattr(it, "secondary", False):
                return "secondary"
            return "primary"
        return kind
    if mode == COLOR_MAPQ:
        def mq(it):
            q = getattr(it, "mapq", None)
            if q is None:
                return "MAPQ unknown"
            return "MAPQ 0" if q == 0 else ("MAPQ 1–29" if q < 30 else "MAPQ ≥ 30")
        return mq
    if mode in (COLOR_READ_GROUP, COLOR_SAMPLE, COLOR_LIBRARY):
        key = {COLOR_SAMPLE: "SM", COLOR_LIBRARY: "LB"}.get(mode)

        def rg(it):
            v = tag(it, "RG")
            if v is None:
                return "no read group"
            if key is None:
                return str(v)
            return str(groups.get(str(v), {}).get(key, v))
        return rg
    if mode == COLOR_TAG:
        if not tag_name:
            return None
        return lambda it: (f"{tag_name}: {v}" if (v := tag(it, tag_name)) is not None
                           else f"no {tag_name}")
    return None


SORT_MODES = [
    ("base", "Base"),
    ("start", "Start"),
    (COLOR_STRAND, "Strand"),
    (COLOR_MAPQ, "Mapping quality"),
    (COLOR_INSERT, "Insert size"),
    (COLOR_READ_GROUP, "Read group"),
    (COLOR_SAMPLE, "Sample"),
    (COLOR_TAG, "Tag…"),
    ("name", "Read name"),
]


def base_at(it, pos: int) -> str:
    """Read base aligned to reference position ``pos``: a letter, '-' (deletion or
    intron) or '' (no bases kept)."""
    seq = getattr(it, "seq", "") or ""
    blocks = getattr(it, "blocks", None) or [(it.start, it.end, 0)]
    for s, e, q in blocks:
        if s <= pos < e:
            k = q + pos - s
            return seq[k].upper() if 0 <= k < len(seq) else ""
    return "-"


def sort_key_function(mode: str, pos: int, items=(), *, reference: str = "",
                      tag_name: str = "", read_groups: dict | None = None):
    """item -> sort key for the reads over ``pos`` (IGV's "Sort alignments by")."""
    if mode == "base":
        ref = reference[pos].upper() if 0 <= pos < len(reference) else ""
        counts: dict[str, int] = {}
        bases = {}
        for it in items:
            b = bases[id(it)] = base_at(it, pos)
            counts[b] = counts.get(b, 0) + 1

        def by_base(it):
            b = bases.get(id(it))
            if b is None:
                b = base_at(it, pos)
            # differing bases first, the most frequent first; the reference base last
            return (b == ref, -counts.get(b, 0), b)
        return by_base
    if mode == "start":
        return lambda it: (it.start, it.end)
    if mode == COLOR_STRAND:
        return lambda it: (it.strand == "-", it.start)
    if mode == COLOR_MAPQ:
        return lambda it: (-(getattr(it, "mapq", None) or 0), it.start)
    if mode == COLOR_INSERT:
        return lambda it: (-abs(mate(it)[2]), it.start)
    if mode in (COLOR_READ_GROUP, COLOR_SAMPLE):
        fn = group_function(mode, read_groups=read_groups)
        return lambda it: (fn(it), it.start)
    if mode == COLOR_TAG:
        def by_tag(it):
            v = tag(it, tag_name) if tag_name else None
            return (v is None, str(v) if v is not None else "", it.start)
        return by_tag
    if mode == "name":
        return lambda it: (it.name or "",)
    return lambda it: (it.start,)


# ---------------------------------------------------------------- filters

FILTER_DEFAULTS = {"hide_duplicates": True, "hide_qcfail": True, "hide_secondary": False,
                   "hide_supplementary": False, "min_mapq": 0}


def read_filter(settings: dict):
    """(key, item -> keep?) for the track's filter settings (IGV's defaults: duplicates
    and QC-failed reads hidden); (None, None) when nothing is hidden."""
    opts = {k: settings.get(k, v) for k, v in FILTER_DEFAULTS.items()}
    mask = 0
    if opts["hide_duplicates"]:
        mask |= FLAG_DUPLICATE
    if opts["hide_qcfail"]:
        mask |= FLAG_QCFAIL
    if opts["hide_secondary"]:
        mask |= FLAG_SECONDARY
    if opts["hide_supplementary"]:
        mask |= FLAG_SUPPLEMENTARY
    min_mapq = int(opts["min_mapq"] or 0)
    if not mask and min_mapq <= 0:
        return None, None
    hide_secondary_attr = bool(opts["hide_secondary"] and opts["hide_supplementary"])

    def keep(it) -> bool:
        f = flag(it)
        if f is not None:
            if f & mask:
                return False
        elif hide_secondary_attr and getattr(it, "secondary", False):
            return False                      # PAF / BLAST: no flags, only primary or not
        if min_mapq > 0:
            mq = getattr(it, "mapq", None)
            if mq is not None and mq < min_mapq:
                return False
        return True

    return (mask, min_mapq), keep


# ---------------------------------------------------------------- feature sequences

_COMPLEMENT = str.maketrans("ACGTRYKMBVDHNacgtrykmbvdhn", "TGCAYRMKVBHDNtgcayrmkvbhdn")


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def feature_sequences(it, fetch) -> list[tuple[str, str, int]]:
    """(key, label, length) of the sequences a feature can give: the genomic span,
    and for a gene model the spliced exons and the CDS. ``fetch(start, end)`` reads
    the reference."""
    out = [("genomic", "Genomic", it.end - it.start)]
    blocks = getattr(it, "blocks", None) or []
    if len(blocks) > 1:
        out.append(("spliced", "Spliced exons", sum(e - s for s, e, _q in blocks)))
    thick = getattr(it, "thick", None) or []
    if thick:
        out.append(("cds", "CDS", sum(e - s for s, e in thick)))
    return out


def feature_sequence(it, fetch, kind: str = "genomic") -> str:
    """Bases of a feature from the reference, 5'→3' on its strand (reverse complement
    on the minus strand)."""
    if kind == "spliced" and getattr(it, "blocks", None):
        parts = [(s, e) for s, e, _q in sorted(it.blocks)]
    elif kind == "cds" and getattr(it, "thick", None):
        parts = sorted(it.thick)
    else:
        parts = [(it.start, it.end)]
    seq = "".join(fetch(s, e) for s, e in parts)
    return reverse_complement(seq) if getattr(it, "strand", "+") == "-" else seq
