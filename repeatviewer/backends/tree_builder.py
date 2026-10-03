"""Phylogenetic tree building (Tools → Tree).

RepeatViewer prepares, checks and reports; established programs compute:

* **MAFFT** (bundled, or from PATH) aligns unaligned input, by default with
  L-INS-i and strand adjustment (TE / satellite copies from both strands).
* **IQ-TREE** (maximum likelihood; ModelFinder by BIC, independent searches,
  SH-aLRT + UFBoot2 with ``-bnni``) is the recommended engine.
* **FastME** (balanced minimum evolution, NNI + SPR, bootstrap) and
  **FastTree** (approximate ML, > 1,000 sequences) are the alternatives.

Nothing is silently downgraded: a missing engine is an error with an install
hint, never a quiet switch to a weaker method.

Like every RepeatViewer tool, the results go next to the input and keep the
full input name, so the file panel can show them green and open the report:
``<input>_tree_report.html``, ``<input>_tree.treefile``, ``<input>_tree.pdf`` /
``.svg`` / ``.png`` and ``<input>_tree.aln.fa``. Supporting files (IQ-TREE
``.iqtree`` / ``.log``, ``.names.tsv``, ``_methods.txt``, ``_commands.txt``,
``.excluded.fa``) go into ``<input>_tree_details/``. A second run on the same
input becomes ``<input>_tree_2…`` — nothing is overwritten.

Pure functions (checks, statistics, split handling, rooting, reports) have no
Qt dependency and are unit-tested in ``tests/test_tree_builder.py``.
"""

from __future__ import annotations

import functools
import html
import io
import re
import shutil
import statistics
import subprocess
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Callable, Iterable

from repeatviewer.job_options import (
    TREE_COVERAGE_LONGEST,
    TREE_ALIGN_ALWAYS,
    TREE_ALIGN_NEVER,
    TREE_ENGINE_FASTME,
    TREE_ENGINE_FASTTREE,
    TREE_ENGINE_IQTREE,
    TREE_MODEL_FINDER,
    TREE_ROOT_MIDPOINT,
    TREE_ROOT_OUTGROUP,
    TREE_SEQ_RT,
    TreeOptions,
)
from repeatviewer.tool_output_paths import kept_input_name

RUN_STAMP_FORMAT = "%Y-%m-%d_%H%M"
MIN_TAXA = 4  # UFBoot and any meaningful unrooted topology need >= 4 sequences

_IQTREE_NAMES = ("iqtree3", "iqtree2", "iqtree")
_FASTME_NAMES = ("fastme",)
_FASTTREE_NAMES = ("FastTreeMP", "FastTree", "fasttree")
_ENGINE_EXECUTABLES = {
    TREE_ENGINE_IQTREE: _IQTREE_NAMES,
    TREE_ENGINE_FASTME: _FASTME_NAMES,
    TREE_ENGINE_FASTTREE: _FASTTREE_NAMES,
}
ENGINE_NAMES = {  # for reports and methods text
    TREE_ENGINE_IQTREE: "Maximum likelihood (IQ-TREE)",
    TREE_ENGINE_FASTME: "Balanced minimum evolution (FastME)",
    TREE_ENGINE_FASTTREE: "Approximate maximum likelihood (FastTree)",
}
ENGINE_LABELS = {  # for the dialog
    TREE_ENGINE_IQTREE: "Maximum likelihood (IQ-TREE)",
    TREE_ENGINE_FASTME: "Minimum evolution (FastME)",
    TREE_ENGINE_FASTTREE: "Approximate maximum likelihood (FastTree)",
}
_CONDA_PACKAGE = {
    TREE_ENGINE_IQTREE: "iqtree",
    TREE_ENGINE_FASTME: "fastme",
    TREE_ENGINE_FASTTREE: "fasttree",
}

MAFFT_STRATEGY_LABELS = (
    ("L-INS-i (accurate, up to ~200 sequences)", "linsi"),
    ("G-INS-i (similar over the whole length)", "ginsi"),
    ("E-INS-i (long unalignable regions)", "einsi"),
    ("FFT-NS-2 (fast, large sets)", "fft_ns_2"),
    ("Auto (MAFFT decides)", "auto"),
)
MAFFT_CHOOSE = "choose"          # Tree window default: L-INS-i ≤ 200 sequences, else FFT-NS-2
MAFFT_CHOOSE_LIMIT = 200

_DNA_CHARS = set("ACGTUNRYSWKMBDHV-.?")
_NUCLEOTIDES = set("ACGTU")
_GAP_CHARS = set("-.?")


class TreeBuildError(RuntimeError):
    """A tree cannot be built; the message is shown to the user as is."""


# --------------------------------------------------------------------------
# Engines
# --------------------------------------------------------------------------


def engine_uses_system_path(opts: TreeOptions, engine: str | None = None) -> bool:
    """True when ``engine`` (default: the selected one) must come from PATH only.

    The old single ``use_system_path`` switch still counts for every engine."""
    engine = engine or opts.engine
    per_engine = {
        TREE_ENGINE_IQTREE: opts.use_system_iqtree,
        TREE_ENGINE_FASTME: opts.use_system_fastme,
        TREE_ENGINE_FASTTREE: opts.use_system_fasttree,
    }.get(engine, False)
    return bool(opts.use_system_path or per_engine)


def engine_executable(engine: str, *, system_path_only: bool = False) -> str | None:
    """Executable for ``engine``.

    ``system_path_only`` (the "Use from system PATH" box): PATH only. Otherwise the
    bundled program (IQ-TREE), then the conda env from Settings, then PATH."""
    names = _ENGINE_EXECUTABLES.get(engine, ())
    if system_path_only:
        for name in names:
            found = shutil.which(name)
            if found:
                return found
        return None
    from repeatviewer.backends.tree_portable_tools import bundled_engine_executable
    from repeatviewer.conda_tools import which_tool

    bundled = bundled_engine_executable(engine)
    if bundled:
        return bundled
    for name in names:
        found = which_tool(name) or shutil.which(name)
        if found:
            return found
    return None


def engine_status_executable(engine: str, *, system_path_only: bool = False) -> str | None:
    """Fast Options-tab lookup: bundled file or PATH — no conda probe."""
    names = _ENGINE_EXECUTABLES.get(engine, ())
    if system_path_only:
        for name in names:
            found = shutil.which(name)
            if found:
                return found
        return None
    from repeatviewer.backends.tree_portable_tools import bundled_engine_executable

    return bundled_engine_executable(engine)


def engine_install_hint(engine: str) -> str:
    from repeatviewer.backends.tree_portable_tools import (
        ENGINE_IQTREE,
        IQTREE_VERSION,
        bundled_engine_executable,
    )

    if engine == ENGINE_IQTREE:
        if bundled_engine_executable(engine):
            return (
                f"IQ-TREE {IQTREE_VERSION} is bundled with RepeatViewer. "
                "Leave \"Use IQ-TREE from the system PATH or defined Conda environment\" "
                "off to use it, or tick that box to use your own install."
            )
        return (
            f"IQ-TREE {IQTREE_VERSION} should be bundled with RepeatViewer but was not found. "
            "Reinstall the app, or tick "
            "\"Use IQ-TREE from the system PATH or defined Conda environment\" and install "
            "iqtree yourself (e.g. conda install -c bioconda iqtree)."
        )
    pkg = _CONDA_PACKAGE.get(engine, engine)
    return (
        f"Install it into the conda env set in File → Setting → External tools:\n"
        f"  conda install -c bioconda {pkg}\n"
        f"or put it on your system PATH and tick "
        f"\"Use … from the system PATH or defined Conda environment\"."
    )


def tool_version(executable: str, engine: str) -> str:
    """First informative line of the program's version output ('' if unknown)."""
    args = {
        TREE_ENGINE_IQTREE: ["--version"],
        TREE_ENGINE_FASTME: ["--version"],
        TREE_ENGINE_FASTTREE: ["-expert"],
    }.get(engine, ["--version"])
    try:
        from repeatviewer.conda_tools import prepare_command

        result = subprocess.run(
            prepare_command([executable, *args]), capture_output=True, text=True, timeout=20
        )
    except (OSError, subprocess.SubprocessError):
        return ""
    text = result.stdout + "\n" + result.stderr
    m = re.search(r"(?:version|Version|v)\s*(\d+\.\d+(?:\.\d+)*)", text) or re.search(
        r"\b(\d+\.\d+(?:\.\d+)+)\b", text
    )
    return m.group(1) if m else ""


# --------------------------------------------------------------------------
# Input handling and checks
# --------------------------------------------------------------------------


@dataclass
class SeqRecordLite:
    name: str       # name shown in outputs (Newick-safe, unique)
    original: str   # header as found in the input
    source: str     # file it came from
    seq: str        # upper case, as read (may contain gaps)


def newick_safe(name: str) -> str:
    """Make a sequence name safe for Newick / IQ-TREE while staying readable."""
    cleaned = re.sub(r"[^A-Za-z0-9_.|\-]+", "_", name.strip())
    return cleaned.strip("_") or "seq"


def read_sequences(files: Iterable[Path]) -> list[SeqRecordLite]:
    """Read FASTA records from one or more files; names made unique and safe."""
    from Bio import SeqIO

    records: list[SeqRecordLite] = []
    used: dict[str, int] = {}
    for path in files:
        path = Path(path)
        for rec in SeqIO.parse(str(path), "fasta"):
            base = newick_safe(rec.id)
            n = used.get(base, 0)
            used[base] = n + 1
            name = base if n == 0 else f"{base}_{n + 1}"
            records.append(
                SeqRecordLite(name=name, original=rec.description, source=path.name,
                              seq=str(rec.seq).upper())
            )
    return records


def looks_like_dna(records: list[SeqRecordLite]) -> bool:
    text = "".join(r.seq for r in records[:200])
    if not text:
        return True
    return sum(c in _DNA_CHARS for c in text) / len(text) >= 0.9


def is_aligned(records: list[SeqRecordLite]) -> bool:
    """Equal lengths = treated as an existing alignment."""
    return len({len(r.seq) for r in records}) == 1


def ungapped_length(seq: str) -> int:
    return sum(c not in _GAP_CHARS for c in seq)


def _reference_length(lengths: list[int], reference: str) -> float:
    if not lengths:
        return 0.0
    return float(max(lengths)) if reference == TREE_COVERAGE_LONGEST else float(
        statistics.median(lengths))


def coverage_filter(
    records: list[SeqRecordLite], min_coverage: float, reference: str = "median"
) -> tuple[list[SeqRecordLite], list[tuple[SeqRecordLite, float]]]:
    """Drop fragments shorter than ``min_coverage`` x the reference length.

    ``reference`` is the median length (default; robust when most sequences are
    complete) or the longest sequence (when most sequences may be fragments, or
    all should be about the same length, e.g. satellite monomers).
    Returns (kept, excluded-with-coverage). ``min_coverage`` 0 disables.
    Fragments (5'-truncated copies, partial monomers) otherwise land at random
    places in the tree and pull down support for everything else.
    """
    if min_coverage <= 0 or not records:
        return list(records), []
    median = _reference_length([ungapped_length(r.seq) for r in records], reference)
    if median <= 0:
        return list(records), []
    kept: list[SeqRecordLite] = []
    excluded: list[tuple[SeqRecordLite, float]] = []
    for r in records:
        cov = ungapped_length(r.seq) / median
        if cov >= min_coverage:
            kept.append(r)
        else:
            excluded.append((r, cov))
    return kept, excluded


@dataclass
class AlignmentStats:
    n_sequences: int
    length: int
    informative_sites: int
    constant_sites: int
    gap_fraction: float


def alignment_stats(seqs: list[str], *, dna: bool = True) -> AlignmentStats:
    """Parsimony-informative sites, constant sites and gap share of an MSA."""
    if not seqs:
        return AlignmentStats(0, 0, 0, 0, 0.0)
    length = len(seqs[0])
    informative = constant = gaps = 0
    valid = _NUCLEOTIDES if dna else None
    for i in range(length):
        counts: dict[str, int] = {}
        for s in seqs:
            c = s[i] if i < len(s) else "-"
            if c in _GAP_CHARS:
                gaps += 1
                continue
            if valid is not None and c not in valid:
                continue  # ambiguity codes / N are not states
            if dna and c == "U":
                c = "T"
            counts[c] = counts.get(c, 0) + 1
        if len(counts) <= 1:
            constant += 1
        elif sum(1 for v in counts.values() if v >= 2) >= 2:
            informative += 1
    total = length * len(seqs)
    return AlignmentStats(len(seqs), length, informative, constant,
                          gaps / total if total else 0.0)


@dataclass
class InputSummary:
    n_sequences: int
    names: list[str]
    dna: bool
    aligned: bool
    length: int | None           # alignment length when already aligned
    informative_sites: int | None
    fragments: int               # would be excluded at the given min_coverage
    error: str = ""
    lengths: list[int] = field(default_factory=list)   # ungapped, for live recounts


def fragment_count(lengths: list[int], min_coverage: float, reference: str = "median") -> int:
    """How many sequences ``coverage_filter`` would exclude."""
    if min_coverage <= 0 or not lengths:
        return 0
    median = _reference_length(lengths, reference)
    return 0 if median <= 0 else sum(1 for n in lengths if n / median < min_coverage)


def input_summary(files: list[Path], min_coverage: float) -> InputSummary:
    """Quick look at the selected input for the dialog's input-check strip."""
    try:
        records = read_sequences(files)
    except Exception as exc:  # noqa: BLE001 - shown in the dialog, never raised
        return InputSummary(0, [], True, False, None, None, 0, error=str(exc))
    if not records:
        return InputSummary(0, [], True, False, None, None, 0,
                            error="No FASTA sequences in the selection.")
    dna = looks_like_dna(records)
    aligned = is_aligned(records)
    _kept, excluded = coverage_filter(records, min_coverage)
    length = informative = None
    if aligned:
        st = alignment_stats([r.seq for r in records], dna=dna)
        length, informative = st.length, st.informative_sites
    return InputSummary(len(records), [r.name for r in records], dna, aligned,
                        length, informative, len(excluded),
                        lengths=[ungapped_length(r.seq) for r in records])


# --------------------------------------------------------------------------
# Trees: splits, support, rooting (Bio.Phylo)
# --------------------------------------------------------------------------


def read_newick(text: str):
    from Bio import Phylo

    return Phylo.read(io.StringIO(text.strip()), "newick")


def tree_to_newick(tree) -> str:
    from Bio import Phylo

    buf = io.StringIO()
    Phylo.write(tree, buf, "newick")
    return buf.getvalue().strip()


def splits(tree) -> set[frozenset[str]]:
    """Non-trivial bipartitions, each as the side NOT containing the first
    taxon (alphabetically), so rooted and unrooted drawings compare equal."""
    taxa = sorted(t.name for t in tree.get_terminals())
    if len(taxa) < 4:
        return set()
    ref = taxa[0]
    all_taxa = frozenset(taxa)
    out: set[frozenset[str]] = set()
    for clade in tree.get_nonterminals():
        side = frozenset(t.name for t in clade.get_terminals())
        if ref in side:
            side = all_taxa - side
        if 1 < len(side) < len(taxa) - 1:
            out.add(side)
    return out


def rf_distance(tree_a, tree_b) -> int:
    return len(splits(tree_a) ^ splits(tree_b))


def support_from_replicates(main_tree, replicate_trees: list) -> None:
    """Label each internal clade of ``main_tree`` with the % of replicates
    that contain its split (classic bootstrap proportion)."""
    if not replicate_trees:
        return
    counts: dict[frozenset[str], int] = {}
    for rep in replicate_trees:
        for s in splits(rep):
            counts[s] = counts.get(s, 0) + 1
    taxa = sorted(t.name for t in main_tree.get_terminals())
    ref, all_taxa = taxa[0], frozenset(taxa)
    for clade in main_tree.get_nonterminals():
        side = frozenset(t.name for t in clade.get_terminals())
        if ref in side:
            side = all_taxa - side
        if 1 < len(side) < len(taxa) - 1:
            clade.confidence = None
            clade.name = str(round(100 * counts.get(side, 0) / len(replicate_trees)))


def _split_labels(tree) -> dict[frozenset[str], object]:
    taxa = sorted(t.name for t in tree.get_terminals())
    ref, all_taxa = taxa[0], frozenset(taxa)
    labels = {}
    for clade in tree.get_nonterminals():
        side = frozenset(t.name for t in clade.get_terminals())
        if ref in side:
            side = all_taxa - side
        if 1 < len(side) < len(taxa) - 1:
            labels[side] = (clade.name, clade.confidence)
    return labels


def _restore_split_labels(tree, labels: dict[frozenset[str], object]) -> None:
    """Re-attach support labels by split after rerooting.

    Supports belong to branches (splits), but Newick stores them on nodes;
    rerooting reverses parent/child along the old-root→new-root path and
    would leave those labels on the wrong branches."""
    taxa = sorted(t.name for t in tree.get_terminals())
    ref, all_taxa = taxa[0], frozenset(taxa)
    for clade in tree.get_nonterminals():
        side = frozenset(t.name for t in clade.get_terminals())
        if ref in side:
            side = all_taxa - side
        name, conf = labels.get(side, (None, None))
        clade.name, clade.confidence = name, conf
    # A root inserted on a branch splits it in two; both halves carry the same
    # split, so show its support once (on the larger side).
    kids = tree.root.clades
    if len(kids) == 2 and not kids[0].is_terminal() and not kids[1].is_terminal():
        a, b = kids
        if (a.name, a.confidence) == (b.name, b.confidence):
            smaller = a if a.count_terminals() <= b.count_terminals() else b
            smaller.name, smaller.confidence = None, None


def apply_rooting(tree, mode: str, outgroup: list[str]) -> str:
    """Root ``tree`` in place, keeping supports on their branches.

    Returns a sentence for the report."""
    labels = _split_labels(tree) if len(tree.get_terminals()) >= 4 else {}
    try:
        return _apply_rooting(tree, mode, outgroup)
    finally:
        if labels:
            _restore_split_labels(tree, labels)


def _apply_rooting(tree, mode: str, outgroup: list[str]) -> str:
    if mode == TREE_ROOT_OUTGROUP and outgroup:
        present = [t for t in tree.get_terminals() if t.name in set(outgroup)]
        if not present:
            return "Unrooted (none of the outgroup sequences are in the tree)."
        # Root on the middle of the branch above the outgroup (clade). First move
        # the root into the ingroup so the outgroup's ancestor is not the root.
        wanted = {t.name for t in present}
        ingroup = [t for t in tree.get_terminals() if t.name not in wanted]
        if ingroup:
            tree.root_with_outgroup(ingroup[0])
        target = present[0] if len(present) == 1 else tree.common_ancestor(present)
        half = (target.branch_length or 0.0) / 2.0
        tree.root_with_outgroup(target, outgroup_branch_length=half)
        names = ", ".join(t.name for t in present)
        return f"Rooted with outgroup {names}."
    if mode == TREE_ROOT_MIDPOINT:
        tree.root_at_midpoint()
        return ("Midpoint-rooted (assumes roughly equal rates of evolution; "
                "an outgroup is preferable when one exists).")
    tree.rooted = False
    return "Unrooted (drawn with an arbitrary root; no rooting method applied)."


def rename_leaves(tree, mapping: dict[str, str]) -> None:
    for leaf in tree.get_terminals():
        if leaf.name in mapping:
            leaf.name = mapping[leaf.name]


# --------------------------------------------------------------------------
# Rendering
# --------------------------------------------------------------------------


def render_tree(tree, out_base: Path, *, title: str = "",
                formats: tuple[str, ...] = ("svg", "pdf", "png")) -> list[Path]:
    """Draw the tree (support values on branches) in ``formats``."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from Bio import Phylo

    n = max(len(tree.get_terminals()), 1)
    height = max(3.0, min(0.22 * n + 1.4, 180.0))
    plt.rcParams["svg.fonttype"] = "none"
    fig, ax = plt.subplots(figsize=(9.0, height))

    def _branch_label(clade):
        if clade.is_terminal():
            return None
        if clade.name:
            return clade.name
        if clade.confidence is not None:
            return f"{clade.confidence:g}"
        return None

    Phylo.draw(
        tree, axes=ax, do_show=False, show_confidence=False,
        label_func=lambda c: c.name if c.is_terminal() else None,
        branch_labels=_branch_label,
    )
    leaf_names = {t.name for t in tree.get_terminals()}
    for text in ax.texts:
        text.set_fontsize(7 if n > 60 else 8)
        if text.get_text().strip() not in leaf_names:
            # Support value: sit just above its branch instead of on the line.
            text.set_verticalalignment("bottom")
            text.set_horizontalalignment("center")
            text.set_fontsize(6 if n > 60 else 7)
            text.set_color("#555555")
    ax.set_ylabel("")
    ax.set_yticks([])
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.set_xlabel("substitutions per site", fontsize=8)
    if title:
        ax.set_title(title, fontsize=9)
    fig.tight_layout()
    written = []
    for ext in formats:
        path = out_base.with_name(out_base.name + f".{ext}")
        fig.savefig(str(path), dpi=200 if ext == "png" else None)
        written.append(path)
    plt.close(fig)
    return written


# --------------------------------------------------------------------------
# Output layout and reporting
# --------------------------------------------------------------------------


def run_stamp(now: datetime | None = None) -> str:
    return (now or datetime.now()).strftime(RUN_STAMP_FORMAT)


def output_base_name(inputs: list[Path]) -> str:
    """``CroSat3.fa`` → ``CroSat3.fa``; several files → ``<first>_multi``."""
    first = kept_input_name(inputs[0])
    return first if len(inputs) == 1 else f"{first}_multi"


# Only the tree files sit next to the input; the report, the alignment and
# everything else go into the ``<base>_details`` folder (decided 25 Sep 2026).
# Figures come from the tree viewer (Save ▾ → Figure as shown).
MAIN_SUFFIXES = (".treefile", "_masked.treefile")


@dataclass
class OutputLayout:
    parent: Path    # folder of the (first) input: main results go here
    base: str       # "<input>_tree" or "<input>_tree_2", …
    details: Path   # "<base>_details" folder for supporting files

    def path(self, suffix: str) -> Path:
        folder = self.parent if suffix in MAIN_SUFFIXES else self.details
        return folder / f"{self.base}{suffix}"


def run_tag(opts: TreeOptions, *, dna: bool = True) -> str:
    """Name tag for runs that differ from a plain tree, so they stand apart:
    ``_rot`` (monomers rotated), ``_cg`` / ``_cgchg`` (methylation check),
    ``_rt`` (RT domains; rotation and masking do not apply)."""
    if opts.sequence_mode == TREE_SEQ_RT:
        return "_rt"
    tag = ""
    if opts.rotate_monomers and dna and opts.align_mode != TREE_ALIGN_NEVER:
        tag += "_rot"
    if opts.methylation_mask and dna:
        tag += "_cgchg" if opts.methylation_mask == "CG_CHG" else "_cg"
    return tag


def plan_outputs(inputs: list[Path], tag: str = "") -> OutputLayout:
    """Pick ``<input>_tree<tag>`` (or ``…_2``, ``…_3`` when a previous run's
    results exist) and create its details folder. Never overwrites."""
    parent = Path(inputs[0]).expanduser().resolve().parent
    stem = f"{output_base_name(inputs)}_tree{tag}"
    n = 1
    while True:
        base = stem if n == 1 else f"{stem}_{n}"
        taken = any((parent / f"{base}{suf}").exists() for suf in MAIN_SUFFIXES)
        details = parent / f"{base}_details"
        if not taken:
            try:
                details.mkdir(exist_ok=False)
                return OutputLayout(parent, base, details)
            except FileExistsError:
                pass
        n += 1


def iqtree_model_from_report(text: str) -> str:
    m = re.search(r"Best-fit model according to BIC:\s*(\S+)", text)
    if m:
        return m.group(1)
    m = re.search(r"Model of substitution:\s*(\S+)", text)
    return m.group(1) if m else ""


def iqtree_loglik_from_report(text: str) -> float | None:
    m = re.search(r"Log-likelihood of the tree:\s*(-?[\d.]+)", text)
    return float(m.group(1)) if m else None


def parse_runtrees(text: str) -> list[tuple[float | None, str]]:
    """IQ-TREE ``.runtrees``: ``[ lh=-850.13 ](newick);`` per line."""
    out = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        m = re.match(r"\[\s*lh=(-?[\d.eE+-]+)\s*\]\s*(.*)$", line)
        if m:
            out.append((float(m.group(1)), m.group(2)))
        else:
            out.append((None, line))
    return out


@dataclass
class TreeRunReport:
    input_files: list[str]
    n_input: int
    n_used: int
    excluded: list[tuple[str, float]]
    dna: bool
    aligned_by: str                 # "" when the input was used as is
    reversed_by_mafft: list[str]
    stats: AlignmentStats | None
    engine: str
    engine_version: str
    mafft_version: str
    model: str
    loglik: float | None
    runs: int
    runs_same_topology: int | None
    support: str
    rooting: str
    seed: int
    warnings: list[str] = field(default_factory=list)
    commands: list[str] = field(default_factory=list)
    rotation: str = ""              # monomer rotation summary ("" = not run)
    rotation_reference: str = ""
    masking: str = ""               # sensitivity paragraph ("" = not run)
    masked_lost: list[list[str]] = field(default_factory=list)
    rt_summary: str = ""            # RT-domain mode summary ("" = not used)
    rt_methods: str = ""


def reference_phrase(opts: TreeOptions) -> str:
    return ("longest sequence" if opts.coverage_reference == TREE_COVERAGE_LONGEST
            else "median sequence length")


def _plural(n: int, word: str) -> str:
    return f"{n} {word}" + ("" if n == 1 else "s")


def _rooting_sentence(rooting: str) -> str:
    if rooting.startswith("Rooted with outgroup "):
        names = rooting.removeprefix("Rooted with outgroup ").rstrip(".")
        return f"The tree was rooted with {names} as outgroup."
    if rooting.startswith("Midpoint"):
        return "The tree was midpoint-rooted."
    return "The tree is unrooted."


def methods_paragraph(r: TreeRunReport, opts: TreeOptions) -> str:
    """A methods-section paragraph describing exactly what was run."""
    parts: list[str] = []
    seq_word = "nucleotide" if r.dna else "protein"
    if r.rt_methods:
        parts.append(r.rt_methods)
    if r.rotation:
        parts.append(
            "Before alignment, monomers were rotated to a common start and orientation "
            f"relative to the monomer of median length ({r.rotation_reference}), using "
            "shared k-mer diagonals (RepeatViewer)."
        )
    if r.aligned_by:
        direction = (" with automatic strand adjustment (--adjustdirectionaccurately)"
                     if opts.adjust_direction and r.dna else "")
        parts.append(
            f"{_plural(r.n_used, seq_word + ' sequence')} were aligned with MAFFT"
            f"{(' ' + r.mafft_version) if r.mafft_version else ''} "
            f"using the {r.aligned_by} strategy{direction}."
        )
    else:
        parts.append(f"The provided alignment of {_plural(r.n_used, seq_word + ' sequence')} "
                     "was used.")
    if r.excluded:
        n = len(r.excluded)
        parts.append(
            f"{_plural(n, 'sequence')} shorter than {opts.min_coverage:.0%} of the "
            f"{reference_phrase(opts)} {'was' if n == 1 else 'were'} excluded."
        )
    if r.stats:
        parts.append(
            f"The alignment comprised {r.stats.length} columns, "
            f"{r.stats.informative_sites} of which were parsimony-informative."
        )
    ver = f" {r.engine_version}" if r.engine_version else ""
    if r.engine == TREE_ENGINE_IQTREE:
        model = (f"the best-fit model {r.model}, selected by ModelFinder according to the "
                 "Bayesian information criterion"
                 if opts.model.upper() == TREE_MODEL_FINDER else f"the {r.model or opts.model} model")
        parts.append(f"A maximum-likelihood tree was inferred with IQ-TREE{ver} under {model}.")
        if r.runs > 1:
            parts.append(
                f"The tree with the highest likelihood from {r.runs} independent searches "
                f"(random seeds {r.seed}–{r.seed + r.runs - 1}) was retained."
            )
        else:
            parts.append(f"A single tree search was run (random seed {r.seed}).")
        measures = []
        if opts.alrt:
            measures.append(f"{opts.alrt} SH-aLRT replicates")
        if opts.ufboot:
            measures.append(f"{max(1000, opts.ufboot)} ultrafast bootstrap replicates "
                            f"(UFBoot2{', with the -bnni correction' if opts.bnni else ''})")
        if measures:
            parts.append("Branch support was assessed with " + " and ".join(measures) + ".")
    elif r.engine == TREE_ENGINE_FASTME:
        dist = "Tamura–Nei (TN93) distances" if r.dna else "LG protein distances"
        parts.append(
            f"A balanced minimum-evolution tree was inferred with FastME{ver} from {dist}, "
            "with NNI and SPR tree rearrangements."
        )
        if opts.distance_bootstrap:
            parts.append(f"Branch support was assessed with "
                         f"{_plural(opts.distance_bootstrap, 'bootstrap replicate')}.")
    else:
        model = "GTR" if r.dna else "LG"
        parts.append(
            f"An approximately maximum-likelihood tree was inferred with FastTree{ver} "
            f"under the {model}+CAT model with Gamma20 rescaling; branch support is given as "
            "SH-like local support values."
        )
    parts.append(_rooting_sentence(r.rooting))
    if r.masking:
        parts.append(r.masking)
    return " ".join(p for p in parts if p)


CITATIONS = {
    "mafft": "Katoh K, Standley DM (2013) MAFFT multiple sequence alignment software version 7. Mol Biol Evol 30:772-780.",
    TREE_ENGINE_IQTREE: "Minh BQ et al. (2020) IQ-TREE 2: new models and efficient methods for phylogenetic inference in the genomic era. Mol Biol Evol 37:1530-1534.",
    "iqtree3": "Wong T et al. (2026) IQ-TREE 3: phylogenomic inference software using complex evolutionary models. Mol Biol Evol 43:msag117.",
    "modelfinder": "Kalyaanamoorthy S et al. (2017) ModelFinder. Nat Methods 14:587-589.",
    "ufboot": "Hoang DT et al. (2018) UFBoot2: improving the ultrafast bootstrap approximation. Mol Biol Evol 35:518-522.",
    "alrt": "Guindon S et al. (2010) New algorithms and methods to estimate maximum-likelihood phylogenies: PhyML 3.0. Syst Biol 59:307-321.",
    TREE_ENGINE_FASTME: "Lefort V, Desper R, Gascuel O (2015) FastME 2.0. Mol Biol Evol 32:2798-2800.",
    TREE_ENGINE_FASTTREE: "Price MN, Dehal PS, Arkin AP (2010) FastTree 2. PLoS ONE 5:e9490.",
    "rexdb": "Neumann P, Novák P, Hoštáková N, Macas J (2019) Systematic survey of plant LTR-retrotransposons elucidates phylogenetic relationships of their polyprotein domains and provides a reference for element classification. Mobile DNA 10:1.",
    "hmmer": "Eddy SR (2011) Accelerated profile HMM searches. PLoS Comput Biol 7:e1002195.",
    "pyhmmer": "Larralde M, Zeller G (2023) PyHMMER: a Python library binding to HMMER for efficient sequence analysis. Bioinformatics 39:btad214.",
}


def citations_for(r: TreeRunReport, opts: TreeOptions) -> list[str]:
    out = []
    if r.rt_summary:
        out += [CITATIONS["rexdb"], CITATIONS["hmmer"], CITATIONS["pyhmmer"]]
    if r.aligned_by:
        out.append(CITATIONS["mafft"])
    if r.engine == TREE_ENGINE_IQTREE and r.engine_version.startswith("3"):
        out.append(CITATIONS["iqtree3"])
    else:
        out.append(CITATIONS[r.engine])
    if r.engine == TREE_ENGINE_IQTREE:
        if opts.model.upper() == TREE_MODEL_FINDER:
            out.append(CITATIONS["modelfinder"])
        if opts.ufboot:
            out.append(CITATIONS["ufboot"])
        if opts.alrt:
            out.append(CITATIONS["alrt"])
    return out


def report_html(r: TreeRunReport, opts: TreeOptions, svg_text: str, files: list[Path],
                title: str) -> str:
    e = html.escape
    rows = [
        ("Input", ", ".join(r.input_files)),
        ("Sequences used", f"{r.n_used} of {r.n_input}"),
        ("Alignment", r.aligned_by or "input used as is"),
    ]
    if r.stats:
        rows += [
            ("Alignment length", f"{r.stats.length} columns"),
            ("Parsimony-informative sites", str(r.stats.informative_sites)),
            ("Gaps", f"{r.stats.gap_fraction:.1%}"),
        ]
    rows += [
        ("Method", ENGINE_NAMES.get(r.engine, r.engine) + (f" {r.engine_version}" if r.engine_version else "")),
    ]
    if r.engine == TREE_ENGINE_IQTREE:
        rows.append(("Model", r.model + (" (ModelFinder, BIC)" if opts.model.upper() == TREE_MODEL_FINDER else "")))
        if r.loglik is not None:
            rows.append(("Log-likelihood", f"{r.loglik:.3f}"))
        if r.runs > 1 and r.runs_same_topology is not None:
            rows.append(("Independent searches", f"{r.runs_same_topology} of {r.runs} found the same topology"))
    support = r.support.removeprefix("Branch support: ")
    rows += [("Branch support", support[:1].upper() + support[1:]), ("Root", r.rooting),
             ("Seed", str(r.seed))]
    if r.rotation:
        rows.insert(3, ("Monomer rotation", r.rotation))
    if r.rt_summary:
        rows.insert(1, ("Sequences", r.rt_summary))
    table = "\n".join(f"<tr><th>{e(k)}</th><td>{e(v)}</td></tr>" for k, v in rows)
    warn = "".join(f"<li>{e(w)}</li>" for w in r.warnings)
    excl = "".join(f"<li>{e(n)} ({c:.0%} of the {e(reference_phrase(opts))})</li>"
                   for n, c in r.excluded)
    # Main results sit next to the input; everything else is in <base>_details/.
    main = [p for p in files if not p.parent.name.endswith("_details")]
    extra = [p for p in files if p not in main]
    flist = "".join(f"<li>{e(p.name)}</li>" for p in main)
    if extra:
        flist += (f"<li>{e(extra[0].parent.name)}/<ul>"
                  + "".join(f"<li>{e(p.name)}</li>" for p in extra) + "</ul></li>")
    cites = "".join(f"<li>{e(c)}</li>" for c in citations_for(r, opts))
    masking = ""
    if r.masking:
        lost = "".join(f"<li>{e(', '.join(names))}</li>" for names in r.masked_lost)
        masking = (f'<h2>Methylation-masked sensitivity run</h2><p class="methods">{e(r.masking)}</p>'
                   + (f"<p>Well-supported groupings of the full tree not recovered:</p><ul>{lost}</ul>"
                      if lost else ""))
    methods = e(methods_paragraph(r, opts))
    return f"""<!doctype html><html><head><meta charset="utf-8"><title>{e(title)}</title>
<style>
body{{font-family:-apple-system,Segoe UI,Helvetica,Arial,sans-serif;margin:18px;color:#1d1d1b;background:#fff}}
h1{{font-size:17px;font-weight:600;margin:0 0 12px}} h2{{font-size:14px;margin:18px 0 6px}}
table{{border-collapse:collapse;font-size:13px}} th{{text-align:left;color:#555;font-weight:500;padding:3px 14px 3px 0;vertical-align:top}}
td{{padding:3px 0}} .warn{{background:#fff4e0;border:1px solid #f0c36d;border-radius:6px;padding:8px 12px;font-size:13px}}
.tree svg{{max-width:100%;height:auto}} .methods{{font-size:13px;line-height:1.5;background:#f6f6f4;border-radius:6px;padding:10px 12px}}
ul{{font-size:13px;margin:4px 0 0 18px;padding:0}} small{{color:#666}}
</style></head><body>
<h1>{e(title)}</h1>
{f'<div class="warn"><b>Check before use</b><ul>{warn}</ul></div>' if warn else ''}
<div class="tree">{svg_text}</div>
<h2>Run report</h2><table>{table}</table>
{f'<h2>Excluded fragments</h2><ul>{excl}</ul>' if excl else ''}
{masking}
<h2>Methods text</h2><div class="methods">{methods}</div>
<h2>Cite</h2><ul>{cites}</ul>
<h2>Files written</h2><ul>{flist}</ul>
<p><small>Support values are shown as SH-aLRT / UFBoot2 for IQ-TREE (commonly
trusted at ≥ 80 / ≥ 95), bootstrap % for FastME and SH-like local support for FastTree.</small></p>
</body></html>"""


def iqtree_thread_args(threads: int) -> list[str]:
    """``-T AUTO`` capped at the user's thread count; ``-T 1`` for one thread.

    AUTO lets IQ-TREE use a single thread on short alignments, where more
    threads only slow it down (and trigger its 'too many threads' warning)."""
    threads = max(1, int(threads))
    return ["-T", "1"] if threads == 1 else ["-T", "AUTO", "--threads-max", str(threads)]


_STACK_LINE = re.compile(r"^ERROR:\s*(\d+\s|\*{3}|STACK TRACE|$)")


@functools.lru_cache(maxsize=8)
def fastme_supports_threads(executable: str) -> bool:
    """FastME builds without OpenMP reject --nb_threads; ask its --help first."""
    try:
        from repeatviewer.conda_tools import prepare_command

        result = subprocess.run(prepare_command([executable, "--help"]), capture_output=True,
                                text=True, timeout=20)
    except (OSError, subprocess.SubprocessError):
        return False
    return "nb_threads" in (result.stdout or "") + (result.stderr or "")


MIN_SHARED_SITES = 20     # fewer shared aligned positions: no distance can be measured
UNRELATED_P = 0.5         # median p-distance at which copies look unrelated


def pairwise_overlap(rows: list[str], *, dna: bool = True):
    """(shared, differing, residues): aligned positions shared and differing for
    every pair (n × n arrays), and each row's number of residues."""
    import numpy as np

    arr = np.frombuffer("".join(r.upper() for r in rows).encode("ascii", "replace"),
                        dtype=np.uint8).reshape(len(rows), -1)
    if dna:
        valid = np.isin(arr, np.frombuffer(b"ACGTU", dtype=np.uint8))
    else:
        valid = ~np.isin(arr, np.frombuffer(b"-.?X*", dtype=np.uint8))
    n = len(rows)
    shared = np.zeros((n, n), dtype=np.int64)
    diff = np.zeros((n, n), dtype=np.int64)
    for i in range(n):
        both = valid[i] & valid
        shared[i] = both.sum(axis=1)
        diff[i] = (both & (arr != arr[i])).sum(axis=1)
    return shared, diff, valid.sum(axis=1)


def overlap_problems(names: list[str], rows: list[str], *, dna: bool = True,
                     min_shared: int = MIN_SHARED_SITES) -> tuple[list[tuple[str, str, int]], float]:
    """Pairs that share too few aligned positions to be compared, and the median
    p-distance over the other pairs (nan when there are none).

    "Too few" is fewer than ``min_shared``, or fewer than half of the shorter
    sequence when that is less (so short repeats are not flagged)."""
    import numpy as np

    shared, diff, residues = pairwise_overlap(rows, dna=dna)
    need = np.minimum(min_shared, np.maximum(1, np.minimum.outer(residues, residues) // 2))
    iu = np.triu_indices(len(rows), k=1)
    s, d = shared[iu], diff[iu]
    too_few = [(names[i], names[j], int(shared[i, j]))
               for i, j in zip(*iu) if shared[i, j] < need[i, j]]
    ok = s >= need[iu]
    median_p = float(np.median(d[ok] / s[ok])) if ok.any() else float("nan")
    return too_few, median_p


def summarize_tool_error(tool: str, returncode: int, output: str) -> str:
    """A readable one-paragraph reason for a failed external program."""
    text = output or ""
    if "OMP: Error" in text or "unable to allocate necessary resources for OMP" in text:
        return (
            f"{tool} stopped because the system could not create more threads "
            f"(OpenMP). Lower 'Threads' in the Tree tab or close other heavy programs, "
            f"then run again."
        )
    if "bad_alloc" in text or "Not enough memory" in text or "out of memory" in text.lower():
        return f"{tool} ran out of memory. Reduce the number of sequences or threads."
    if "Unable to compute all distances" in text:
        return (
            f"{tool} could not measure the distance between some pairs of sequences: they "
            "share no aligned positions or are too different. The sequences are probably "
            "not all copies of one family."
        )
    errors = [line.strip() for line in text.splitlines()
              if line.strip().startswith("ERROR") and not _STACK_LINE.match(line.strip())]
    if errors:
        return f"{tool} failed (exit {returncode}): " + " ".join(errors[:3])
    tail = [line.strip() for line in text.splitlines() if line.strip()][-4:]
    return f"{tool} failed (exit {returncode})." + ("\n" + "\n".join(tail) if tail else "")


def _add_command(work: Path, commands: list[str], text: str) -> None:
    """Record a command line (also on disk, so a failure log can include it)."""
    commands.append(text)
    try:
        with (work / "commands.txt").open("a", encoding="utf-8") as fh:
            fh.write(text + "\n")
    except OSError:
        pass


def failure_log_path(inputs: list[Path]) -> Path:
    """``<input>_tree_failed.log`` next to the (first) input."""
    first = Path(inputs[0]).expanduser().resolve()
    return first.parent / f"{output_base_name(inputs)}_tree_failed.log"


def _write_failure_log(inputs: list[Path], run_dir: Path, error: str) -> Path | None:
    parts = [f"RepeatViewer tree builder — failed run ({datetime.now():%Y-%m-%d %H:%M})",
             "", "Reason:", error, ""]
    work = run_dir / ".work"
    for name in ("commands.txt",):
        f = work / name
        if f.is_file():
            parts += [f"== {name} ==", f.read_text(encoding="utf-8", errors="replace"), ""]
    for log in sorted(work.rglob("*.log")) if work.is_dir() else []:
        parts += [f"== {log.relative_to(work)} ==",
                  log.read_text(encoding="utf-8", errors="replace"), ""]
    path = failure_log_path(inputs)
    try:
        path.write_text("\n".join(parts), encoding="utf-8")
        return path
    except OSError:
        return None


# --------------------------------------------------------------------------
# The run
# --------------------------------------------------------------------------


def _write_fasta(path: Path, items: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8") as fh:
        for name, seq in items:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")


def _write_phylip(path: Path, items: list[tuple[str, str]]) -> None:
    """Sequential PHYLIP with <=10-char IDs (safe IDs are S00001-style)."""
    length = len(items[0][1])
    with path.open("w", encoding="utf-8") as fh:
        fh.write(f"{len(items)} {length}\n")
        for name, seq in items:
            fh.write(f"{name:<10} {seq.replace('.', '-')}\n")


def _read_fasta_pairs(path: Path) -> list[tuple[str, str]]:
    from Bio import SeqIO

    return [(rec.id, str(rec.seq).upper()) for rec in SeqIO.parse(str(path), "fasta")]


def build_tree_run(
    inputs: list[Path],
    opts: TreeOptions,
    *,
    progress: Callable[[str], None] | None = None,
    cancel=None,
    stamp: str | None = None,
) -> tuple[Path, list[Path], Path]:
    """Run the whole pipeline. Returns (details folder, created files, report HTML)."""
    from repeatviewer.external_process import run_logged_process

    say = progress or (lambda _m: None)
    inputs = [Path(p) for p in inputs]
    if not inputs:
        raise TreeBuildError("Select one or more FASTA files first.")

    # --- engine present? (fail early, never downgrade) ---------------------
    exe = engine_executable(opts.engine, system_path_only=engine_uses_system_path(opts))
    if not exe:
        raise TreeBuildError(
            f"{ENGINE_LABELS.get(opts.engine, opts.engine)} is selected, but the program "
            f"was not found.\n\n{engine_install_hint(opts.engine)}"
        )

    records = read_sequences(inputs)
    n_input = len(records)
    if n_input < MIN_TAXA:
        raise TreeBuildError(
            f"A tree needs at least {MIN_TAXA} sequences; the selection has {n_input}."
        )
    dna = looks_like_dna(records)
    warnings: list[str] = []
    rt = None
    if opts.sequence_mode == TREE_SEQ_RT:
        if not dna:
            raise TreeBuildError(
                "'RT domain from DNA' is selected, but the selection looks like protein. "
                "The RT domain is searched in the DNA; choose 'As they are' in the Alignment "
                "tab for protein sequences."
            )
        records, rt = rt_domain_records(records, opts, say, cancel, warnings)
        dna = False
    if rt is None:
        records, excluded = coverage_filter(records, opts.min_coverage, opts.coverage_reference)
    else:
        # RT types differ in domain length (Ty1 vs Ty3): fragments are judged
        # against sequences of the same RT type only.
        records, excluded = _coverage_filter_by_rt_type(records, rt, opts)
    if len(records) < MIN_TAXA:
        raise TreeBuildError(
            f"Only {len(records)} sequence(s) remain after excluding fragments shorter than "
            f"{opts.min_coverage:.0%} of the {reference_phrase(opts)}. Lower 'Minimum length' in the "
            "Alignment tab or check the input."
        )
    outgroup = [o for o in opts.outgroup if o]
    if opts.root_mode == TREE_ROOT_OUTGROUP:
        dropped = [r.name for r, _c in excluded if r.name in outgroup]
        if dropped:
            warnings.append("Outgroup sequence(s) excluded as fragments: " + ", ".join(dropped))
        present = {r.name for r in records} | set(dropped)
        missing = [o for o in outgroup if o not in present]
        if missing:
            # The outgroup is saved with the settings; it may belong to another file.
            warnings.append("Outgroup sequence(s) not in this input: " + ", ".join(missing)
                            + ". Choose the outgroup again in the Support and rooting tab.")
        if not outgroup:
            warnings.append("Outgroup rooting was selected, but no outgroup was chosen.")

    stamp = stamp or run_stamp()
    layout = plan_outputs(inputs, run_tag(opts, dna=dna))
    created: list[Path] = []
    from repeatviewer.backends.job_cancel import JobCancelledError

    def _cleanup() -> None:
        # Failed or stopped: no half-finished results are left behind.
        for path in created:
            if path.parent == layout.parent:
                path.unlink(missing_ok=True)
        shutil.rmtree(layout.details, ignore_errors=True)

    try:
        return _run_in_dir(layout, created, inputs, records, n_input, excluded, dna,
                           outgroup, warnings, opts, stamp, say, cancel, rt=rt)
    except JobCancelledError:
        _cleanup()
        raise
    except Exception as exc:  # noqa: BLE001 - reported with a saved log
        message = str(exc) if isinstance(exc, TreeBuildError) else f"{type(exc).__name__}: {exc}"
        log = _write_failure_log(inputs, layout.details, message)
        _cleanup()
        suffix = f"\n\nFull log: {log}" if log else ""
        raise TreeBuildError(message + suffix) from exc


def _check_cancel(cancel) -> None:
    if cancel is not None and getattr(cancel, "requested", False):
        from repeatviewer.backends.job_cancel import JobCancelledError

        raise JobCancelledError("Job stopped by user")


@dataclass
class _Inference:
    tree: object                    # Bio.Phylo tree, internal sequence IDs, supports attached
    model: str
    loglik: float | None
    runs_same: int | None
    support: str
    best_prefix: Path | None = None   # IQ-TREE: best search
    first_prefix: Path | None = None  # IQ-TREE: first search (model selection)


def _infer(exe, opts, safe_aln, aln_pairs, dna, safe_outgroup, work, commands, say, cancel,
           warnings, *, model_override: str = "", label: str = "") -> _Inference:
    """Build one tree from ``safe_aln`` with the engine in ``opts``.

    Used for the main tree and, with ``label`` and a fixed ``model_override``,
    for the methylation-masked sensitivity tree.
    """
    from repeatviewer.external_process import run_logged_process

    wd = work / label.strip().replace(" ", "_") if label else work
    wd.mkdir(parents=True, exist_ok=True)
    loglik = None
    model = ""
    runs_same = None
    support = ""
    best_prefix = None
    if opts.engine == TREE_ENGINE_IQTREE:
        # Independent searches run as separate IQ-TREE processes (seed, seed+1, …)
        # rather than one ``--runs N`` process: IQ-TREE 3 on macOS exhausts OpenMP
        # threads when many searches share one process, and a crash in one search
        # must not lose the others. ModelFinder runs once, in the first search;
        # later searches reuse the selected model.
        n_runs = max(1, opts.runs)
        results: list[tuple[float, str, Path]] = []  # (logL, newick, prefix)
        model = ""
        failures: list[str] = []
        for i in range(n_runs):
            _check_cancel(cancel)
            prefix = wd / f"run{i + 1}" / "iq"
            prefix.parent.mkdir(parents=True)
            if model_override:
                run_model = model_override
            else:
                run_model = (opts.model or TREE_MODEL_FINDER) if i == 0 or not model else model
            cmd = [exe, "-s", str(safe_aln), "-m", run_model, "--seed", str(opts.seed + i),
                   "--prefix", str(prefix), "-redo"]
            cmd += iqtree_thread_args(opts.threads)
            if getattr(opts, "iqtree_fast", False):
                cmd.append("--fast")
            if opts.ufboot:
                cmd += ["-B", str(max(1000, opts.ufboot))]
                if opts.bnni:
                    cmd.append("-bnni")
            if opts.alrt:
                cmd += ["--alrt", str(opts.alrt)]
            if opts.root_mode == TREE_ROOT_OUTGROUP and safe_outgroup:
                cmd += ["-o", ",".join(safe_outgroup)]
            _add_command(work, commands, " ".join(Path(c).name if k == 0 else c for k, c in enumerate(cmd))
                            .replace(str(work) + "/", ""))
            if n_runs > 1:
                say(f"Building the {label or ''}tree — search {i + 1} of {n_runs}…".replace("  ", " "))
            result = run_logged_process(cmd, tool="IQ-TREE", inputs=[safe_aln], cwd=work,
                                        check=False, cancel=cancel)
            _check_cancel(cancel)
            treefile = prefix.with_suffix(".treefile")
            report_file = prefix.with_suffix(".iqtree")
            if result.returncode != 0 or not treefile.is_file() or not report_file.is_file():
                failures.append(summarize_tool_error(
                    "IQ-TREE", result.returncode, (result.stdout or "") + "\n" + (result.stderr or "")
                ))
                if i == 0:
                    break  # the first search also selects the model; nothing to reuse
                continue
            report_text = report_file.read_text(encoding="utf-8", errors="replace")
            if i == 0:
                model = iqtree_model_from_report(report_text) or model_override or opts.model
            lnl = iqtree_loglik_from_report(report_text)
            results.append((lnl if lnl is not None else float("-inf"),
                            treefile.read_text(encoding="utf-8"), prefix))
        if not results:
            raise TreeBuildError(failures[0] if failures else "IQ-TREE produced no tree.")
        if failures:
            warnings.append(f"{len(failures)} of {n_runs} IQ-TREE {label}searches failed and were "
                            f"skipped: {failures[0].splitlines()[0]}")
        loglik, tree_text, best_prefix = max(results, key=lambda r: r[0])
        best = read_newick(tree_text)
        runs_same = sum(1 for _l, nwk, _p in results if rf_distance(best, read_newick(nwk)) == 0)
        parts = []
        if opts.alrt:
            parts.append(f"{opts.alrt} SH-aLRT replicates")
        if opts.ufboot:
            parts.append(f"{max(1000, opts.ufboot)} ultrafast bootstrap replicates (UFBoot2"
                         + (", -bnni" if opts.bnni else "") + ")")
        support = ("Branch support: " + " and ".join(parts) + " (shown as SH-aLRT/UFBoot)."
                   if parts else "No branch support computed.")
    elif opts.engine == TREE_ENGINE_FASTME:
        phy = wd / "aln.phy"
        _write_phylip(phy, aln_pairs)
        main_out = wd / "fastme.nwk"
        boot_out = wd / "fastme_boot.nwk"
        cmd = [exe, f"--input_data={phy}", f"--output_tree={main_out}",
               "--dna=T" if dna else "--protein=L", "--method=B", "--nni=B", "--spr",
               f"--seed={opts.seed}"]
        if fastme_supports_threads(exe):
            cmd.append(f"--nb_threads={max(1, opts.threads)}")
        if opts.distance_bootstrap:
            cmd += [f"--bootstrap={opts.distance_bootstrap}", "-B", str(boot_out)]
        _add_command(work, commands, " ".join(cmd).replace(str(work) + "/", ""))
        result = run_logged_process(cmd, tool="FastME", inputs=[phy], cwd=work,
                                    check=False, cancel=cancel)
        _check_cancel(cancel)
        if result.returncode != 0 or not main_out.is_file():
            raise TreeBuildError(summarize_tool_error(
                "FastME", result.returncode, (result.stdout or "") + "\n" + (result.stderr or "")))
        tree_text = main_out.read_text(encoding="utf-8").strip().splitlines()[0]
        model = "TN93" if dna else "LG"
        support = (f"Branch support: {opts.distance_bootstrap} nonparametric bootstrap replicates."
                   if opts.distance_bootstrap else "No branch support computed.")
    else:  # FastTree
        cmd = [exe, "-quiet", "-nopr", "-gamma", "-seed", str(opts.seed)]
        cmd += ["-nt", "-gtr"] if dna else ["-lg"]
        cmd.append(str(safe_aln))
        _add_command(work, commands, " ".join(Path(c).name if i == 0 else c for i, c in enumerate(cmd))
                        .replace(str(work) + "/", "") + " > tree.nwk")
        result = run_logged_process(cmd, tool="FastTree", inputs=[safe_aln], cwd=work,
                                    check=False, cancel=cancel, stdout_is_result=True)
        _check_cancel(cancel)
        tree_text = (result.stdout or "").strip()
        if result.returncode != 0 or not tree_text:
            raise TreeBuildError(summarize_tool_error("FastTree", result.returncode,
                                                      result.stderr or ""))
        model = "GTR+CAT+Gamma20" if dna else "LG+CAT+Gamma20"
        support = "Branch support: FastTree SH-like local support values (0-1)."

    tree = read_newick(tree_text)
    if opts.engine == TREE_ENGINE_FASTME and opts.distance_bootstrap and boot_out.is_file():
        reps = [read_newick(line) for line in boot_out.read_text(encoding="utf-8").splitlines()
                if line.strip()]
        support_from_replicates(tree, reps)
    for clade in tree.find_clades():
        if clade.branch_length is not None and clade.branch_length < 0:
            clade.branch_length = 0.0
    return _Inference(tree, model, loglik, runs_same, support, best_prefix,
                      wd / "run1" / "iq" if opts.engine == TREE_ENGINE_IQTREE else None)


def rt_domain_records(records: list[SeqRecordLite], opts: TreeOptions, say, cancel,
                      warnings: list[str]):
    """RT-domain mode: replace the DNA copies by their RT protein domains and add
    REXdb reference RTs. Returns (protein records, RTResult)."""
    from repeatviewer import app_preferences
    from repeatviewer.backends import rt_domains as rd
    from repeatviewer.backends.sequence_aligner import mafft_command

    if not opts.rexdb_fasta:
        raise TreeBuildError(
            "'RT domain from DNA' needs REXdb. Choose the REXdb protein FASTA in the "
            "Alignment tab (the .classification file must sit next to it)."
        )
    cmd = mafft_command(use_system_mafft=opts.use_system_mafft)
    if cmd is None:
        raise TreeBuildError("MAFFT was not found (bundled or on PATH).")
    _check_cancel(cancel)
    try:
        res = rd.rt_domains_for(
            [(r.name, r.seq.replace("-", "").replace(".", "")) for r in records],
            opts.rexdb_fasta, mafft=cmd,
            cache_dir=app_preferences.cache_dir() / "rexdb_profiles",
            per_lineage=opts.rexdb_references, threads=opts.threads, say=say,
        )
    except rd.RexdbError as exc:
        raise TreeBuildError(str(exc)) from exc
    _check_cancel(cancel)
    if len(res.hits) < MIN_TAXA:
        raise TreeBuildError(
            f"An RT domain was found in {len(res.hits)} of {len(records)} sequence(s); "
            f"a tree needs at least {MIN_TAXA}. The copies may be too short or too "
            "degraded, or not retrotransposons."
        )
    if res.missing:
        shown = ", ".join(res.missing[:12]) + ("…" if len(res.missing) > 12 else "")
        warnings.append(f"No RT domain was found in {len(res.missing)} sequence(s); "
                        f"they are not in the tree: {shown}")
    if opts.methylation_mask:
        warnings.append("The CG/CHG masking run does not apply to RT domains; it was skipped.")
    out = [SeqRecordLite(name=r.name, original=r.original, source=r.source,
                         seq=res.hits[r.name].protein)
           for r in records if r.name in res.hits]
    used = {r.name for r in out}
    for rid, lineage, protein in res.references:
        base = newick_safe(rd.reference_name(rid, lineage))
        name, n = base, 1
        while name in used:
            n += 1
            name = f"{base}_{n}"
        used.add(name)
        out.append(SeqRecordLite(name=name, original=rid, source="REXdb",
                                 seq=protein.upper()))
    return out, res


def _coverage_filter_by_rt_type(records: list[SeqRecordLite], rt, opts: TreeOptions):
    def rt_type(r: SeqRecordLite) -> str:
        if r.name in rt.hits:
            return rt.hits[r.name].rt_type
        return rt.rexdb.rt.get(r.original, ("", ""))[0] if rt.rexdb else ""

    groups: dict[str, list[SeqRecordLite]] = {}
    for r in records:
        groups.setdefault(rt_type(r), []).append(r)
    kept_names: set[str] = set()
    excluded = []
    for group in groups.values():
        k, e = coverage_filter(group, opts.min_coverage, opts.coverage_reference)
        kept_names |= {r.name for r in k}
        excluded += e
    return [r for r in records if r.name in kept_names], excluded


def _rt_texts(rt, opts: TreeOptions, n_copies: int) -> tuple[str, str]:
    """(report row, methods sentence) for RT-domain mode."""
    from repeatviewer.backends.rt_domains import DEFAULT_EVALUE

    db = rt.rexdb.name if rt.rexdb else "REXdb"
    types = ", ".join(f"{t} {n}" for t, n in rt.by_type.most_common())
    lineages = len({h.lineage for h in rt.hits.values() if h.lineage})
    row = (f"RT domain from DNA ({db}): found in {len(rt.hits)} of {n_copies} copies "
           f"({types}); {lineages} lineage(s); {len(rt.references)} REXdb reference RT(s) added")
    refs = (f", together with {len(rt.references)} reference RTs from {db} (up to "
            f"{opts.rexdb_references} per lineage of the RT types found)" if rt.references else "")
    methods = (
        "Reverse transcriptase (RT) domains were identified in six-frame translations of "
        f"{n_copies} copies with profile HMMs built from the {db} RT sequences (Neumann et al. "
        f"2019), searched with HMMER through pyhmmer (E-value ≤ {DEFAULT_EVALUE:g}; best domain per "
        "copy). Each domain was assigned to the REXdb lineage whose profile HMM scored it "
        f"highest. The RT protein sequences{refs} were used for the tree."
    )
    return row, methods


def _write_rt_traits(tree_path: Path, rt, records: list[SeqRecordLite]) -> Path:
    """Trait table so the viewer colours tips by lineage straight away."""
    from repeatviewer.backends import tree_traits as tt

    ref_lineage = {}
    for r in records:
        if r.source == "REXdb" and rt.rexdb is not None:
            ref_lineage[r.name] = (rt.rexdb.lineage(r.original), rt.rexdb.superfamily(r.original),
                                   rt.rexdb.rt.get(r.original, ("", ""))[0])
    # "lineage" colours the copies and the references of the lineages the copies
    # fall into (REXdb has ~40 lineages, more than the palette's 17 colours);
    # "rexdb_lineage" names the lineage of every tip.
    found = {h.lineage for h in rt.hits.values() if h.lineage}
    rows: dict[str, dict[str, str]] = {}
    for r in records:
        if r.name in rt.hits:
            h = rt.hits[r.name]
            rows[r.name] = {"lineage": h.lineage, "rexdb_lineage": h.lineage,
                            "superfamily": h.superfamily, "rt_type": h.rt_type,
                            "kind": "copy", "header": r.original}
        else:
            lin, sf, rtt = ref_lineage.get(r.name, ("", "", ""))
            rows[r.name] = {"lineage": lin if lin in found else "", "rexdb_lineage": lin,
                            "superfamily": sf, "rt_type": rtt, "kind": "REXdb",
                            "header": r.original}
    table = tt.TraitTable(["lineage", "rexdb_lineage", "superfamily", "rt_type", "kind",
                           "header"], rows)
    return tt.write_trait_table(tt.traits_path(tree_path), table)


def _run_in_dir(layout, created, inputs, records, n_input, excluded, dna, outgroup,
                warnings, opts, stamp, say, cancel, rt=None):
    exe = engine_executable(opts.engine, system_path_only=engine_uses_system_path(opts))
    out = layout.path
    work = layout.details / ".work"
    work.mkdir()
    commands: list[str] = []

    safe_to_name = {f"S{i:05d}": r.name for i, r in enumerate(records, 1)}
    name_to_safe = {v: k for k, v in safe_to_name.items()}
    with out(".names.tsv").open("w", encoding="utf-8") as fh:
        fh.write("tree_name\toriginal_header\tsource_file\n")
        for r in records:
            fh.write(f"{r.name}\t{r.original}\t{r.source}\n")
    created.append(out(".names.tsv"))
    if excluded:
        _write_fasta(out(".excluded.fa"), [(r.name, r.seq.replace("-", "")) for r, _c in excluded])
        created.append(out(".excluded.fa"))
    rt_summary = rt_methods = ""
    if rt is not None:
        out(".rt_domains.tsv").write_text(rt.table_tsv(), encoding="utf-8")
        created.append(out(".rt_domains.tsv"))
        _write_fasta(out(".rt_proteins.fa"), [(r.name, r.seq) for r in records])
        created.append(out(".rt_proteins.fa"))
        rt_summary, rt_methods = _rt_texts(rt, opts, len(rt.hits) + len(rt.missing))

    # --- monomer phase rotation (satellites) ----------------------------------
    rotation_summary = ""
    rotation_reference = ""
    rotate = bool(opts.rotate_monomers) and rt is None
    if rotate and not dna:
        warnings.append("Monomer rotation applies to DNA only; the protein input was not rotated.")
        rotate = False
    if rotate and opts.align_mode == TREE_ALIGN_NEVER:
        warnings.append("Monomer rotation was skipped: it needs a new alignment, but "
                        "'Use the input as an alignment' is selected.")
        rotate = False
    if rotate:
        from repeatviewer.backends.monomer_phase import rotate_to_common_start

        say("Rotating monomers to a common start…")
        rot = rotate_to_common_start([(r.name, r.seq) for r in records])
        rotated = dict(rot.sequences)
        for r in records:
            r.seq = rotated[r.name]
        out(".rotation.tsv").write_text(rot.table_tsv(), encoding="utf-8")
        created.append(out(".rotation.tsv"))
        rotation_reference = rot.reference
        rotation_summary = (f"{rot.n_rotated} of {len(records)} rotated, "
                            f"{rot.n_reversed} reverse-complemented (reference {rot.reference})")
        if rot.unmatched:
            shown = ", ".join(rot.unmatched[:12]) + ("…" if len(rot.unmatched) > 12 else "")
            warnings.append(
                f"{len(rot.unmatched)} monomer(s) shared too few k-mers with the reference "
                f"{rot.reference} to be rotated and were left as they are: {shown}"
            )

    # --- alignment -----------------------------------------------------------
    already_aligned = is_aligned(records)
    do_align = opts.align_mode == TREE_ALIGN_ALWAYS or rotate or (
        opts.align_mode != TREE_ALIGN_NEVER and not already_aligned
    )
    if opts.align_mode == TREE_ALIGN_NEVER and not is_aligned(records):
        raise TreeBuildError(
            "'Use the input as an alignment' is selected, but the sequences differ in "
            "length. Choose 'Align if not aligned' in the Alignment tab."
        )
    reversed_ids: list[str] = []
    aligned_by = ""
    mafft_version = ""
    safe_aln = work / "aln.fa"
    if do_align:
        from repeatviewer.backends.sequence_aligner import _run_mafft, mafft_command

        cmd = mafft_command(use_system_mafft=opts.use_system_mafft)
        if cmd is None:
            raise TreeBuildError("MAFFT was not found (bundled or on PATH).")
        command, env = cmd
        raw = work / "unaligned.fa"
        _write_fasta(raw, [(name_to_safe[r.name], r.seq.replace("-", "").replace(".", ""))
                           for r in records])
        strategy = opts.mafft_strategy
        if strategy == MAFFT_CHOOSE:
            # RepeatViewer picks: the accurate L-INS-i up to ~200 sequences, else FFT-NS-2.
            strategy = "linsi" if len(records) <= MAFFT_CHOOSE_LIMIT else "fft_ns_2"
        label = dict((v, k) for k, v in MAFFT_STRATEGY_LABELS).get(strategy, strategy)
        if strategy in ("linsi", "ginsi", "einsi") and len(records) > 200:
            warnings.append(
                f"{len(records)} sequences with {label.split(' (')[0]}: accurate but slow; "
                "FFT-NS-2 is the usual choice above ~200 sequences."
            )
        say("Aligning with MAFFT…")
        mafft_out = work / "aln_raw.fa"
        _run_mafft(
            raw, mafft_out, command=command, env=env, num_threads=opts.threads,
            strategy=strategy, sequence_type="nucleotide" if dna else "protein",
            retree=2, max_iterate=0, gap_open_penalty=opts.mafft_gap_open,
            offset=opts.mafft_offset, reorder=False,
            adjust_direction=bool(opts.adjust_direction and dna), cancel=cancel,
        )
        _check_cancel(cancel)
        _add_command(work, commands, 
            "mafft " + {"linsi": "--localpair --maxiterate 1000", "ginsi": "--globalpair --maxiterate 1000",
                        "einsi": "--genafpair --maxiterate 1000", "fft_ns_2": "--retree 2 --maxiterate 0",
                        "auto": "--auto"}.get(strategy, strategy)
            + f" --op {opts.mafft_gap_open:g} --ep {opts.mafft_offset:g}"
            + (" --adjustdirectionaccurately" if opts.adjust_direction and dna else "")
            + f" --thread {opts.threads} unaligned.fa > aln.fa"
        )
        pairs = []
        for sid, seq in _read_fasta_pairs(mafft_out):
            if sid.startswith("_R_"):
                sid = sid[3:]
                reversed_ids.append(safe_to_name[sid])
            pairs.append((sid, seq))
        _write_fasta(safe_aln, pairs)
        aligned_by = label.split(" (")[0]
        try:
            vres = subprocess.run([*command, "--version"], capture_output=True, text=True,
                                  env=env, timeout=20)
            mver = re.search(r"v\d+\.\d+", vres.stdout + vres.stderr)
            mafft_version = mver.group(0) if mver else ""
        except (OSError, subprocess.SubprocessError):
            pass
    else:
        _write_fasta(safe_aln, [(name_to_safe[r.name], r.seq.replace(".", "-")) for r in records])

    aln_pairs = _read_fasta_pairs(safe_aln)
    if len({len(s) for _n, s in aln_pairs}) != 1:
        raise TreeBuildError("The alignment has rows of different lengths.")
    _write_fasta(out(".aln.fa"), [(safe_to_name[n], s) for n, s in aln_pairs])
    created.append(out(".aln.fa"))
    stats = alignment_stats([s for _n, s in aln_pairs], dna=dna)
    if len(aln_pairs) <= 2000:
        too_few, median_p = overlap_problems([safe_to_name[n] for n, _s in aln_pairs],
                                             [s for _n, s in aln_pairs], dna=dna)
        if too_few:
            shown = "; ".join(f"{a} and {b} ({k})" for a, b, k in too_few[:6])
            more = f" and {len(too_few) - 6} more" if len(too_few) > 6 else ""
            message = (
                f"{len(too_few)} pair(s) of sequences share (almost) no aligned positions "
                f"({shown}{more}; shared positions in brackets), so the distance between "
                "them cannot be measured. The sequences are probably not all copies of one family (for "
                "example consensus sequences of different satellite families); a tree needs "
                "related copies."
            )
            if opts.engine == TREE_ENGINE_FASTME:
                raise TreeBuildError(message + " FastME cannot build a tree from this alignment.")
            warnings.append(message)
        elif median_p == median_p and median_p >= UNRELATED_P:
            warnings.append(
                f"The sequences differ at {median_p:.0%} of their shared aligned positions "
                "(median over all pairs), about what unrelated sequences give. Check that they "
                "are copies of one family."
            )
    if stats.informative_sites < 10:
        warnings.append(
            f"Only {stats.informative_sites} parsimony-informative sites: most branches "
            "cannot be resolved, whatever the method."
        )
    if reversed_ids:
        warnings.append(
            f"MAFFT reverse-complemented {len(reversed_ids)} sequence(s) to match the "
            "others' strand: " + ", ".join(reversed_ids[:12])
            + ("…" if len(reversed_ids) > 12 else "")
        )

    # --- inference -----------------------------------------------------------
    engine_version = tool_version(exe, opts.engine)
    safe_outgroup = [name_to_safe[o] for o in outgroup if o in name_to_safe]
    loglik = None
    model = ""
    runs_same = None
    support = ""
    say(f"Building the tree ({ENGINE_NAMES[opts.engine]})…")
    inf = _infer(exe, opts, safe_aln, aln_pairs, dna, safe_outgroup, work, commands, say,
                 cancel, warnings)
    tree, model, loglik, runs_same, support = (inf.tree, inf.model, inf.loglik,
                                               inf.runs_same, inf.support)
    if opts.engine == TREE_ENGINE_IQTREE:
        best_prefix = inf.best_prefix
        # Keep the report and log of the best search, plus the model selection
        # (first search) when that was a different one.
        first_prefix = work / "run1" / "iq"
        keep = [(best_prefix, "")]
        if best_prefix != first_prefix and first_prefix.with_suffix(".iqtree").is_file():
            keep.append((first_prefix, ".modelfinder"))
        for src_prefix, tag in keep:
            for ext in (".iqtree", ".log"):
                src = src_prefix.with_suffix(ext)
                if src.is_file():
                    text = src.read_text(encoding="utf-8", errors="replace")
                    for sid, name in safe_to_name.items():
                        text = text.replace(sid, name)
                    dest = out(f"{tag}{ext}")
                    dest.write_text(text, encoding="utf-8")
                    created.append(dest)
    rename_leaves(tree, safe_to_name)
    if opts.root_mode == TREE_ROOT_OUTGROUP and len(outgroup) > 1:
        present = [t for t in tree.get_terminals() if t.name in set(outgroup)]
        if len(present) > 1 and tree.is_monophyletic(present) is False:
            warnings.append(
                "The outgroup sequences do not form one clade in this tree, so the root "
                "position is uncertain. Use a single outgroup or midpoint rooting."
            )
    rooting = apply_rooting(tree, opts.root_mode, outgroup)
    tree.ladderize()
    out(".treefile").write_text(tree_to_newick(tree) + "\n", encoding="utf-8")
    created.append(out(".treefile"))
    if rt is not None:
        created.append(_write_rt_traits(out(".treefile"), rt, records))

    # --- methylation-masked sensitivity tree ----------------------------------
    masking_text = ""
    masking_lost: list[list[str]] = []
    if opts.methylation_mask and not dna and rt is None:
        warnings.append("The CG/CHG masking run applies to DNA only; it was skipped.")
    elif opts.methylation_mask and dna:
        from repeatviewer.backends import methylation_mask as mm

        _check_cancel(cancel)
        rows = [s for _n, s in aln_pairs]
        cols = mm.methylation_columns(rows, opts.methylation_mask)
        masked_rows = mm.remove_columns(rows, cols)
        masked_pairs = [(n, s) for (n, _s), s in zip(aln_pairs, masked_rows)]
        mstats = alignment_stats(masked_rows, dna=True) if masked_rows[0] else None
        if mstats is None or mstats.informative_sites < 10:
            warnings.append(
                f"After removing {len(cols)} {mm.MASK_LABELS[opts.methylation_mask]} columns, "
                f"only {mstats.informative_sites if mstats else 0} parsimony-informative sites "
                "remain; the masked tree was not built."
            )
        else:
            mdir = work / "masked"
            mdir.mkdir(exist_ok=True)
            masked_aln = mdir / "aln.fa"
            _write_fasta(masked_aln, masked_pairs)
            say("Building the methylation-masked tree…")
            fixed = model if opts.engine == TREE_ENGINE_IQTREE else ""
            minf = _infer(exe, opts, masked_aln, masked_pairs, dna, safe_outgroup, work,
                          commands, say, cancel, warnings, model_override=fixed,
                          label="masked ")
            mtree = minf.tree
            rename_leaves(mtree, safe_to_name)
            apply_rooting(mtree, opts.root_mode, outgroup)
            mtree.ladderize()
            out("_masked.treefile").write_text(tree_to_newick(mtree) + "\n", encoding="utf-8")
            created.append(out("_masked.treefile"))
            _write_fasta(out("_masked.aln.fa"),
                         [(safe_to_name[n], s) for n, s in masked_pairs])
            created.append(out("_masked.aln.fa"))
            cmp = mm.compare_trees(tree, mtree, opts.engine)
            masking_text = mm.comparison_sentence(cmp, len(cols), len(rows[0]),
                                                  opts.methylation_mask)
            if fixed:
                masking_text += (f" The masked tree used the model selected for the full "
                                 f"tree ({fixed}).")
            masking_lost = cmp.lost_supported
            if cmp.lost_supported:
                warnings.append(
                    f"{len(cmp.lost_supported)} well-supported grouping(s) of the full tree "
                    "are not recovered without CG/CHG sites; see the sensitivity section."
                )

    title = f"{output_base_name(inputs)} — {ENGINE_NAMES[opts.engine].split(' (')[0].lower()} tree"
    say("Drawing the tree…")
    # Drawing for the report only (the report embeds it); not kept as files.
    images = render_tree(tree, work / "report_tree", title="", formats=("svg",))

    report = TreeRunReport(
        input_files=[p.name for p in inputs], n_input=n_input, n_used=len(records),
        excluded=[(r.name, c) for r, c in excluded], dna=dna, aligned_by=aligned_by,
        reversed_by_mafft=reversed_ids, stats=stats, engine=opts.engine,
        engine_version=engine_version, mafft_version=mafft_version, model=model,
        loglik=loglik, runs=opts.runs if opts.engine == TREE_ENGINE_IQTREE else 1,
        runs_same_topology=runs_same, support=support, rooting=rooting, seed=opts.seed,
        warnings=warnings, commands=commands, rotation=rotation_summary,
        rotation_reference=rotation_reference, masking=masking_text, masked_lost=masking_lost,
        rt_summary=rt_summary, rt_methods=rt_methods,
    )
    methods = methods_paragraph(report, opts)
    out("_methods.txt").write_text(
        methods + "\n\nReferences\n" + "\n".join(citations_for(report, opts)) + "\n",
        encoding="utf-8",
    )
    created.append(out("_methods.txt"))
    versions = [f"RepeatViewer tree builder, run {stamp}"]
    if mafft_version:
        versions.append(f"MAFFT {mafft_version}")
    versions.append(f"{opts.engine}: {engine_version or 'version unknown'} ({exe})")
    out("_commands.txt").write_text(
        "\n".join(versions) + "\n\n" + "\n".join(commands)
        + "\n\nSequence IDs S00001… map to names in " + out(".names.tsv").name + "\n",
        encoding="utf-8",
    )
    created.append(out("_commands.txt"))
    svg_text = next(p for p in images if p.suffix == ".svg").read_text(encoding="utf-8")
    svg_text = svg_text[svg_text.find("<svg"):]
    report_path = out("_report.html")
    report_path.write_text(
        report_html(report, opts, svg_text, [*created, report_path], title), encoding="utf-8"
    )
    created.append(report_path)
    shutil.rmtree(work, ignore_errors=True)
    return layout.details, created, report_path
