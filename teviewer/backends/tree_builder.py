"""Phylogenetic tree building from a multiple sequence alignment.

Two engines, mirroring how alignment backends are handled in
``teviewer.backends.sequence_aligner``:

* **External (preferred when installed)** -- IQ-TREE or FastTree found on PATH.
  Maximum-likelihood trees, far faster on large sets, with branch support
  values included. Install with::

      conda install -c bioconda iqtree fasttree

  These are optional: like BLAST+ and system MAFFT they are detected at runtime
  and never required for the app to work.

* **Biopython (always available)** -- distance-based Neighbor-Joining or UPGMA.
  No extra installation, works on every platform, but slows down noticeably
  past a few hundred sequences.

Distance corrections (Jukes-Cantor, Kimura 2-parameter, Tamura-Nei) are
implemented here because Biopython's ``DistanceCalculator`` only ships
substitution *scoring matrices* (``blastn``, ``trans``, ...), not nucleotide
substitution models.

Note on HKY: unlike JC69/K2P/TN93 it has no closed-form distance estimator --
it is fitted by maximum likelihood. Use an external engine for HKY-style
modelling; the distance path offers TN93, which generalises it.
"""

from __future__ import annotations

import math
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO
from Bio.Phylo.BaseTree import Tree
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

# --------------------------------------------------------------------------
# Options (mirror the Tree dialog)
# --------------------------------------------------------------------------

ENGINE_AUTO = "auto"
ENGINE_BIOPYTHON = "biopython"
ENGINE_FASTTREE = "fasttree"
ENGINE_IQTREE = "iqtree"

MODEL_JC69 = "jc69"
MODEL_K2P = "k2p"
MODEL_TN93 = "tn93"
MODEL_P_DISTANCE = "p"

DISTANCE_MODEL_LABELS = (
    ("Jukes-Cantor", MODEL_JC69),
    ("Kimura 2-parameter", MODEL_K2P),
    ("Tamura-Nei", MODEL_TN93),
    ("Uncorrected p-distance", MODEL_P_DISTANCE),
)

METHOD_NJ = "nj"
METHOD_UPGMA = "upgma"

BUILD_METHOD_LABELS = (
    ("Neighbor-Joining", METHOD_NJ),
    ("UPGMA", METHOD_UPGMA),
)

NO_OUTGROUP = ""

# Pairs of purines / pyrimidines: a change within a pair is a transition.
_PURINES = frozenset("AG")
_PYRIMIDINES = frozenset("CT")
_VALID_BASES = frozenset("ACGT")

# Distance assigned when a pair is so divergent the model cannot be applied
# (the logarithm goes negative). Keeps the matrix usable instead of crashing.
_SATURATED_DISTANCE = 5.0


class TreeBuildError(RuntimeError):
    """Raised when a tree cannot be built (too few sequences, tool failed...)."""


@dataclass(frozen=True)
class TreeResult:
    newick_path: Path
    engine_used: str
    sequence_count: int
    notes: str = ""


# --------------------------------------------------------------------------
# Engine discovery
# --------------------------------------------------------------------------

_IQTREE_NAMES = ("iqtree3", "iqtree2", "iqtree")
_FASTTREE_NAMES = ("FastTreeMP", "FastTree", "fasttreeMP", "fasttree")


def _find_executable(names: tuple[str, ...]) -> str | None:
    for name in names:
        found = shutil.which(name)
        if found:
            return found
    return None


def iqtree_executable() -> str | None:
    return _find_executable(_IQTREE_NAMES)


def fasttree_executable() -> str | None:
    return _find_executable(_FASTTREE_NAMES)


def available_tree_engines() -> list[str]:
    """Engines usable right now, most preferred first.

    Biopython is always last but always present, so the list is never empty.
    """
    found: list[str] = []
    if iqtree_executable():
        found.append(ENGINE_IQTREE)
    if fasttree_executable():
        found.append(ENGINE_FASTTREE)
    found.append(ENGINE_BIOPYTHON)
    return found


def preferred_engine() -> str:
    """Best available engine: external when installed, else Biopython."""
    return available_tree_engines()[0]


def engine_display_name(engine: str) -> str:
    return {
        ENGINE_IQTREE: "IQ-TREE (maximum likelihood)",
        ENGINE_FASTTREE: "FastTree (approximate maximum likelihood)",
        ENGINE_BIOPYTHON: "Biopython (distance-based)",
        ENGINE_AUTO: "Automatic",
    }.get(engine, engine)


# --------------------------------------------------------------------------
# Alignment input
# --------------------------------------------------------------------------


def read_alignment(path: Path) -> list[tuple[str, str]]:
    """Return [(id, uppercased sequence)] from a FASTA alignment."""
    records = [
        (str(rec.id), str(rec.seq).upper())
        for rec in SeqIO.parse(str(path), "fasta")
    ]
    if not records:
        raise TreeBuildError(f"No sequences found in {path.name}.")
    return records


def is_aligned(records: list[tuple[str, str]]) -> bool:
    """True when every sequence has the same length (i.e. already aligned)."""
    lengths = {len(seq) for _, seq in records}
    return len(lengths) == 1


def check_tree_input(records: list[tuple[str, str]]) -> None:
    """Raise a clear error when a tree cannot be built from these sequences."""
    if len(records) < 3:
        raise TreeBuildError(
            f"A tree needs at least 3 sequences; {len(records)} provided.\n\n"
            "Select three or more sequences, or an alignment containing "
            "three or more sequences."
        )
    if not is_aligned(records):
        raise TreeBuildError(
            "Sequences are not aligned (they differ in length).\n\n"
            "Align them first, then build the tree from the alignment."
        )


# --------------------------------------------------------------------------
# Distance models
# --------------------------------------------------------------------------


def _pair_counts(seq_a: str, seq_b: str) -> tuple[int, int, int, int, int]:
    """Compare two aligned sequences, ignoring gaps and ambiguous bases.

    Returns (compared, transitions_ag, transitions_ct, transversions, matches).
    Positions where either sequence is not an unambiguous A/C/G/T are skipped
    (pairwise deletion), which is what most tree builders do by default.
    """
    compared = ti_ag = ti_ct = tv = matches = 0
    for base_a, base_b in zip(seq_a, seq_b):
        if base_a not in _VALID_BASES or base_b not in _VALID_BASES:
            continue
        compared += 1
        if base_a == base_b:
            matches += 1
            continue
        pair = {base_a, base_b}
        if pair == {"A", "G"}:
            ti_ag += 1
        elif pair == {"C", "T"}:
            ti_ct += 1
        else:
            tv += 1
    return compared, ti_ag, ti_ct, tv, matches


def _base_frequencies(sequences: list[str]) -> dict[str, float]:
    counts = {base: 0 for base in "ACGT"}
    total = 0
    for seq in sequences:
        for base in "ACGT":
            n = seq.count(base)
            counts[base] += n
            total += n
    if total == 0:
        return {base: 0.25 for base in "ACGT"}
    return {base: counts[base] / total for base in "ACGT"}


def _safe_log(value: float) -> float | None:
    """ln(value), or None when the model has broken down (value <= 0)."""
    if value <= 0.0:
        return None
    return math.log(value)


def pairwise_distance(
    seq_a: str,
    seq_b: str,
    model: str,
    freqs: dict[str, float] | None = None,
) -> float:
    """Corrected evolutionary distance between two aligned sequences."""
    compared, ti_ag, ti_ct, tv, matches = _pair_counts(seq_a, seq_b)
    if compared == 0:
        # No comparable positions (all gaps/ambiguous) -- treat as maximally distant.
        return _SATURATED_DISTANCE
    if matches == compared:
        return 0.0

    p_diff = (compared - matches) / compared

    if model == MODEL_P_DISTANCE:
        return p_diff

    if model == MODEL_JC69:
        ln = _safe_log(1.0 - (4.0 / 3.0) * p_diff)
        return _SATURATED_DISTANCE if ln is None else -0.75 * ln

    transitions = (ti_ag + ti_ct) / compared
    transversions = tv / compared

    if model == MODEL_K2P:
        ln1 = _safe_log(1.0 - 2.0 * transitions - transversions)
        ln2 = _safe_log(1.0 - 2.0 * transversions)
        if ln1 is None or ln2 is None:
            return _SATURATED_DISTANCE
        return -0.5 * ln1 - 0.25 * ln2

    if model == MODEL_TN93:
        f = freqs or {base: 0.25 for base in "ACGT"}
        g_a, g_c, g_g, g_t = f["A"], f["C"], f["G"], f["T"]
        g_r = g_a + g_g
        g_y = g_c + g_t
        if g_r <= 0 or g_y <= 0:
            return _SATURATED_DISTANCE
        p1 = ti_ag / compared  # A<->G transitions
        p2 = ti_ct / compared  # C<->T transitions
        q = transversions
        term1 = _safe_log(1.0 - (g_r * p1) / (2.0 * g_a * g_g) - q / (2.0 * g_r)) \
            if g_a > 0 and g_g > 0 else None
        term2 = _safe_log(1.0 - (g_y * p2) / (2.0 * g_c * g_t) - q / (2.0 * g_y)) \
            if g_c > 0 and g_t > 0 else None
        term3 = _safe_log(1.0 - q / (2.0 * g_r * g_y))
        if term1 is None or term2 is None or term3 is None:
            return _SATURATED_DISTANCE
        d = (
            -2.0 * (g_a * g_g / g_r) * term1
            - 2.0 * (g_c * g_t / g_y) * term2
            - 2.0 * (g_r * g_y - (g_a * g_g * g_y / g_r) - (g_c * g_t * g_r / g_y)) * term3
        )
        return max(d, 0.0)

    raise TreeBuildError(f"Unknown distance model: {model}")


def build_distance_matrix(records: list[tuple[str, str]], model: str) -> DistanceMatrix:
    """Lower-triangular distance matrix in the form Bio.Phylo expects."""
    names = [name for name, _ in records]
    seqs = [seq for _, seq in records]
    freqs = _base_frequencies(seqs) if model == MODEL_TN93 else None

    matrix: list[list[float]] = []
    for i in range(len(seqs)):
        row = [
            round(pairwise_distance(seqs[i], seqs[j], model, freqs), 8)
            for j in range(i)
        ]
        row.append(0.0)
        matrix.append(row)
    return DistanceMatrix(names=_unique_names(names), matrix=matrix)


def _unique_names(names: list[str]) -> list[str]:
    """Bio.Phylo requires unique taxon names; disambiguate duplicates."""
    seen: dict[str, int] = {}
    out: list[str] = []
    for name in names:
        if name in seen:
            seen[name] += 1
            out.append(f"{name}_{seen[name]}")
        else:
            seen[name] = 0
            out.append(name)
    return out


# --------------------------------------------------------------------------
# Tree construction
# --------------------------------------------------------------------------


def build_tree_biopython(
    records: list[tuple[str, str]],
    *,
    model: str = MODEL_JC69,
    method: str = METHOD_NJ,
    outgroup: str = NO_OUTGROUP,
) -> Tree:
    """Distance-based tree (no external tools required)."""
    matrix = build_distance_matrix(records, model)
    constructor = DistanceTreeConstructor()
    if method == METHOD_UPGMA:
        tree = constructor.upgma(matrix)
    elif method == METHOD_NJ:
        tree = constructor.nj(matrix)
    else:
        raise TreeBuildError(f"Unknown tree build method: {method}")

    # Biopython labels internal nodes Inner1, Inner2... which clutters the
    # rendered tree and confuses downstream viewers.
    for clade in tree.get_nonterminals():
        if clade.name and clade.name.startswith("Inner"):
            clade.name = None

    _clamp_negative_branches(tree)

    if outgroup and outgroup != NO_OUTGROUP:
        _root_with_outgroup(tree, outgroup)
    return tree


def _clamp_negative_branches(tree: Tree) -> None:
    """Set negative branch lengths to zero.

    Neighbor-Joining can produce negative branches when distances violate
    additivity. They are not biologically meaningful, render oddly, and upset
    some tree viewers, so they are clamped -- the same thing MEGA and PAUP do.
    """
    for clade in tree.find_clades():
        if clade.branch_length is not None and clade.branch_length < 0:
            clade.branch_length = 0.0


def _root_with_outgroup(tree: Tree, outgroup: str) -> None:
    try:
        target = next(
            leaf for leaf in tree.get_terminals() if leaf.name == outgroup
        )
    except StopIteration:
        return  # Outgroup not in this tree; leave the default rooting.
    tree.root_with_outgroup(target)


def _run_tool(cmd: list[str], *, cwd: Path | None = None) -> subprocess.CompletedProcess:
    try:
        return subprocess.run(
            cmd,
            cwd=str(cwd) if cwd else None,
            capture_output=True,
            text=True,
            check=False,
        )
    except OSError as exc:
        raise TreeBuildError(f"Could not run {cmd[0]}: {exc}") from exc


def build_tree_fasttree(alignment_path: Path, out_newick: Path) -> None:
    exe = fasttree_executable()
    if not exe:
        raise TreeBuildError("FastTree is not installed.")
    # -nt: nucleotide, -gtr: GTR model, -quiet/-nopr: keep stderr usable.
    result = _run_tool([exe, "-nt", "-gtr", "-quiet", "-nopr", str(alignment_path)])
    if result.returncode != 0 or not result.stdout.strip():
        raise TreeBuildError(
            f"FastTree failed (exit {result.returncode}).\n{result.stderr.strip()[:500]}"
        )
    out_newick.write_text(result.stdout, encoding="utf-8")


def build_tree_iqtree(alignment_path: Path, out_newick: Path, *, threads: str = "AUTO") -> None:
    exe = iqtree_executable()
    if not exe:
        raise TreeBuildError("IQ-TREE is not installed.")
    with tempfile.TemporaryDirectory(prefix="teviewer_iqtree_") as tmp:
        tmp_dir = Path(tmp)
        prefix = tmp_dir / "tree"
        result = _run_tool(
            [
                exe,
                "-s", str(alignment_path),
                "-m", "MFP",          # ModelFinder picks the substitution model
                "-B", "1000",         # ultrafast bootstrap -> branch supports
                "-T", threads,
                "--prefix", str(prefix),
                "--quiet",
                "-redo",
            ]
        )
        treefile = prefix.with_suffix(".treefile")
        if result.returncode != 0 or not treefile.exists():
            raise TreeBuildError(
                f"IQ-TREE failed (exit {result.returncode}).\n"
                f"{(result.stderr or result.stdout).strip()[:500]}"
            )
        out_newick.write_text(treefile.read_text(encoding="utf-8"), encoding="utf-8")


def build_tree(
    alignment_path: Path,
    out_newick: Path,
    *,
    engine: str = ENGINE_AUTO,
    model: str = MODEL_JC69,
    method: str = METHOD_NJ,
    outgroup: str = NO_OUTGROUP,
) -> TreeResult:
    """Build a tree from an aligned FASTA, writing Newick to ``out_newick``.

    ``engine=auto`` prefers an installed external tool and falls back to
    Biopython. UPGMA and explicit distance models are only available on the
    Biopython path, so requesting them forces that engine.
    """
    records = read_alignment(alignment_path)
    check_tree_input(records)

    resolved = engine
    notes = ""
    if engine == ENGINE_AUTO:
        resolved = preferred_engine()
    elif engine in {ENGINE_IQTREE, ENGINE_FASTTREE}:
        exe = iqtree_executable() if engine == ENGINE_IQTREE else fasttree_executable()
        if not exe:
            resolved = ENGINE_BIOPYTHON
            notes = (
                f"{engine} is not installed; used the built-in distance method "
                "instead. Install with: conda install -c bioconda iqtree fasttree"
            )

    # UPGMA has no equivalent in the ML tools.
    if method == METHOD_UPGMA and resolved != ENGINE_BIOPYTHON:
        resolved = ENGINE_BIOPYTHON
        notes = "UPGMA is only available with the built-in distance method."

    out_newick.parent.mkdir(parents=True, exist_ok=True)

    if resolved == ENGINE_IQTREE:
        build_tree_iqtree(alignment_path, out_newick)
    elif resolved == ENGINE_FASTTREE:
        build_tree_fasttree(alignment_path, out_newick)
        if outgroup:
            _reroot_newick_file(out_newick, outgroup)
    else:
        tree = build_tree_biopython(
            records, model=model, method=method, outgroup=outgroup
        )
        _write_newick(tree, out_newick)

    return TreeResult(
        newick_path=out_newick,
        engine_used=resolved,
        sequence_count=len(records),
        notes=notes,
    )


def _write_newick(tree: Tree, out_newick: Path) -> None:
    from Bio import Phylo

    with out_newick.open("w", encoding="utf-8") as handle:
        Phylo.write(tree, handle, "newick")


def _reroot_newick_file(newick_path: Path, outgroup: str) -> None:
    """Re-root an existing Newick file on ``outgroup`` (external engines)."""
    from Bio import Phylo

    try:
        tree = Phylo.read(str(newick_path), "newick")
    except (ValueError, OSError):
        return
    _root_with_outgroup(tree, outgroup)
    _write_newick(tree, newick_path)


def render_tree_image(
    newick_path: Path,
    out_path: Path,
    *,
    title: str | None = None,
    show_support: bool = True,
) -> Path:
    """Draw a Newick tree to PNG/SVG.

    Forces matplotlib's headless Agg backend before importing pyplot: this runs
    on a worker thread, and the interactive Qt backend would build QTimers off
    the GUI thread and crash.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from Bio import Phylo

    tree = Phylo.read(str(newick_path), "newick")
    leaf_count = max(len(tree.get_terminals()), 1)

    # Grow the canvas with the number of taxa so labels stay legible.
    height = max(3.0, min(0.25 * leaf_count + 1.5, 200.0))
    fig, axes = plt.subplots(figsize=(10.0, height))

    def _label(clade):
        return clade.name or ""

    try:
        Phylo.draw(
            tree,
            axes=axes,
            do_show=False,
            label_func=_label,
            show_confidence=show_support,
        )
    except Exception as exc:  # noqa: BLE001 - drawing must never kill the job
        plt.close(fig)
        raise TreeBuildError(f"Could not draw the tree: {exc}") from exc

    if title:
        axes.set_title(title, fontsize=10)
    axes.set_xlabel("branch length (substitutions per site)", fontsize=8)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(str(out_path), dpi=150)
    plt.close(fig)
    return out_path


def taxon_names(alignment_path: Path) -> list[str]:
    """Sequence ids, for populating the Outgroup dropdown."""
    try:
        return [name for name, _ in read_alignment(alignment_path)]
    except (TreeBuildError, OSError, ValueError):
        return []
