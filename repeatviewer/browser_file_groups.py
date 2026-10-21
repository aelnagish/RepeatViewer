"""Group related Files-panel rows (query + tool outputs) for visual banding."""

from __future__ import annotations

import re
from pathlib import Path

# Longest suffixes first so ``_bn.txt.fa`` wins over ``.fa`` / ``_bn.txt``.
_TOOL_OUTPUT_SUFFIXES: tuple[str, ...] = (
    # BLAST hit FASTA extracted from tables
    "_bn.txt.fa",
    "_bp.txt.fa",
    "_bx.txt.fa",
    "_b.txt.fa",
    "_blast.txt.fa",
    "_blastp.txt.fa",
    "_blastx.txt.fa",
    "_obn.txt.fa",
    "_obp.txt.fa",
    "_obx.txt.fa",
    "_obtn.txt.fa",
    "_obtx.txt.fa",
    "_ob.txt.fa",
    "_online.txt.fa",
    # BLAST / online tables + GFF
    "_bn.txt",
    "_bp.txt",
    "_bx.txt",
    "_b.txt",
    "_blast.txt",
    "_blastp.txt",
    "_blastx.txt",
    "_obn.txt",
    "_obp.txt",
    "_obx.txt",
    "_obtn.txt",
    "_obtx.txt",
    "_ob.txt",
    "_online.txt",
    "_bn.gff",
    "_bp.gff",
    "_bx.gff",
    "_b.gff",
    "_blast.gff",
    "_blastp.gff",
    "_blastx.gff",
    # Online HTML
    "_obn.html",
    "_obp.html",
    "_obx.html",
    "_obtn.html",
    "_obtx.html",
    "_ob.html",
    "_online.html",
    # Align
    "_cmafft.fa",
    "_cmuscle.fa",
    "_mafft.fa",
    "_muscle.fa",
    # Sequence ops
    "_translated.faa",
    "_orfs.fa",
    "_rc.fa",
    "_extr.fa",
    "_hmmer_orfs.faa",
    # MSA clean / consensus / extend (extend matched via regex below)
    "_crowdiv.fa",
    "_crowgap.fa",
    "_cgap.fa",
    "_rowdiv.fa",
    "_rowgap.fa",
    "_gap.fa",
    # Plots / reports
    "_msa_plot.png",
    "_msa_plot.pdf",
    "_den_plot_region.html",
    "_dens_plot_region.html",
    "_den_plot.html",
    "_dens_plot.html",
    "_div_plot.html",
    "_abund_plot.html",
    "_abundance_plot.html",
    # Pan-genome before single-genome: ``_pteaid.html`` also ends with ``_teaid.html``.
    "_pteaid.html",
    "_TEAid.html",
    "_teaid.html",
    # Flexidot
    "_multi_flexidot.png",
    "_multi_flexidot.html",
    "_flexidot_multi.png",  # legacy
    "_flexidot_multi.html",  # legacy
    "_flexidot.png",
    "_flexidot.html",
    # Built-in Dotplot PAF (self / multi)
    "_multi_dot.paf",
    "_self_dot.paf",
    # HMMER / TSD / annotation
    ".tsd_search.txt",
    ".seqmap_RV.json",
    "_hmmer.txt",
    # Map / tree primary (engine tags + legacy `_map` / `_tree`)
    "_minimap.rvmap",
    "_hisat.rvmap",
    "_map.rvmap",
    "_iqtree.treefile",
    "_fastme.treefile",
    "_fasttree.treefile",
    "_tree.treefile",
    "_iqtree_masked.treefile",
    "_fastme_masked.treefile",
    "_fasttree_masked.treefile",
    "_tree_masked.treefile",
    "_masked.treefile",
)

# Numbered / parameterized tails glued before a normal suffix.
_NUMBERED_ALIGN_RE = re.compile(
    r"_(?:c?mafft|c?muscle)(?:_\d+)?\.(?:fa|fasta)$",
    re.IGNORECASE,
)
_NUMBERED_HMMER_RE = re.compile(r"_hmmer\d+\.txt$", re.IGNORECASE)
_NUMBERED_TREE_RE = re.compile(
    r"_(?:iqtree|fastme|fasttree|tree)"
    r"(?:_rt|_rot|_cgchg|_cg)?(?:_\d+)?(?:_masked)?\.treefile$",
    re.IGNORECASE,
)
_NUMBERED_MAP_RE = re.compile(
    r"_(?:minimap|hisat|map)(?:_\d+)?\.rvmap$",
    re.IGNORECASE,
)
_NUMBERED_TEAID_RE = re.compile(r"_(?:p)?TEAid(?:_\d+)?\.html$", re.IGNORECASE)
_NUMBERED_FLEXIDOT_RE = re.compile(
    r"(?:_multi_flexidot|_flexidot_multi|_flexidot)(?:_\d+)?\.(?:png|html)$",
    re.IGNORECASE,
)
_NUMBERED_BUILTIN_DOTPLOT_RE = re.compile(
    r"_(?:self|multi)_dot(?:_\d+)?\.paf$",
    re.IGNORECASE,
)
_EXTEND_RE = re.compile(r"_L\d+_R\d+(?:_\d+)?\.fa$", re.IGNORECASE)
_CONSENSUS_RE = re.compile(r"_\d+(?:\.\d+)?_co(?:_\d+)?\.fa$", re.IGNORECASE)
_FASTGA_RE = re.compile(
    r"_FastGA(?:_vs_.+?)?(?:_\d+)?\.(?:paf|psl|1aln)$"
    r"|_multi_FastGA(?:_vs_.+?)?(?:_\d+)?\.(?:paf|psl|1aln)$",
    re.IGNORECASE,
)
_MSA_CLEAN_RE = re.compile(
    r"_(?:cgap|crowdiv|crowgap|gap|rowdiv|rowgap)(?:_\d+)?\.fa$",
    re.IGNORECASE,
)
_PLOT_TAG_RE = re.compile(
    r"_(?:bn|bp|bx|b|blast|RM|obn|obp|obx|ob|multi)_(?:den|dens|div|abund)_plot"
    r"(?:_region)?\.html$",
    re.IGNORECASE,
)
_COMBINE_RE = re.compile(r"_and_(?:\d+|.+)_comb(?:_\d+)?\.fa$", re.IGNORECASE)

_REGEX_STRIPPERS: tuple[re.Pattern[str], ...] = (
    _PLOT_TAG_RE,
    _FASTGA_RE,
    _MSA_CLEAN_RE,
    _NUMBERED_BUILTIN_DOTPLOT_RE,
    _NUMBERED_FLEXIDOT_RE,
    _NUMBERED_TEAID_RE,
    _NUMBERED_TREE_RE,
    _NUMBERED_MAP_RE,
    _NUMBERED_ALIGN_RE,
    _NUMBERED_HMMER_RE,
    _EXTEND_RE,
    _CONSENSUS_RE,
    _COMBINE_RE,
)

# Stronger than online-BLAST `#fafbfc` so name-family groups read clearly.
FAMILY_BAND_LIGHT = "#eef2f6"
FAMILY_BAND_DARK = "#32363c"
# Simple even/odd zebra when the Files list is not sorted by Name.
ROW_ZEBRA_LIGHT = "#f3f6fa"
ROW_ZEBRA_DARK = "#2c3036"


def file_family_key(name: str) -> str:
    """Return the source-family key for a filename (query + derived outputs).

    ``reads.fa``, ``reads.fa_bn.txt``, and ``reads.fa_bn.txt.fa`` share
    ``reads.fa``. Unknown names keep their own full name as the key.
    """
    current = name
    # Repeat until stable so stacked suffixes (``_bn.txt.fa``) peel cleanly.
    for _ in range(8):
        lowered = current.lower()
        nxt = current
        for pattern in _REGEX_STRIPPERS:
            match = pattern.search(current)
            if match and match.end() == len(current):
                nxt = current[: match.start()]
                break
        else:
            for suffix in _TOOL_OUTPUT_SUFFIXES:
                if lowered.endswith(suffix):
                    nxt = current[: -len(suffix)]
                    break
            else:
                # Numbered plain tool suffix: ``…_mafft_2.fa`` already covered;
                # also ``…_flexidot_2.png`` via regex. Generic ``_2`` before ext
                # after a known stem is handled by regexes above.
                break
        if not nxt or nxt == current:
            break
        current = nxt
    return current.lower() or name.lower()


def family_band_flags(names: list[str]) -> list[bool]:
    """For names in display order, True = tinted band, False = plain.

    Consecutive rows with the same :func:`file_family_key` share a band.
    Band parity alternates when the family key changes (online-BLAST style).
    """
    flags: list[bool] = []
    prev_key: str | None = None
    band = False
    for name in names:
        key = file_family_key(name)
        if prev_key is None:
            band = False
        elif key != prev_key:
            band = not band
        flags.append(band)
        prev_key = key
    return flags


def family_band_flags_for_paths(paths: list[Path]) -> list[bool]:
    return family_band_flags([p.name for p in paths])
