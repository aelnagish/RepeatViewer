"""CG/CHG masking sensitivity run for trees (Tools → Tree).

Methylated cytosines deaminate to thymine, so C→T (and G→A on the other
strand) changes pile up at CG sites and, in plants, CHG sites (H = A, C or T).
In repeats these hypermutable sites can create shared changes that look like
common ancestry. The sensitivity run removes every alignment column where any
sequence has a C or G in such a context, rebuilds the tree with the same
settings, and compares the two trees. If the groupings hold, the result does
not depend on methylation-driven mutation.

Contexts are read on the ungapped sequence, on both strands:

* **CG**: the C of ``CG`` and its G (the C of the other strand);
* **CHG**: the C of ``C[ACT]G``, and the G of ``C[AGT]G`` (that G is the C of
  a CHG on the other strand).

A column is removed when the context occurs in at least one sequence: a C→T
change turns ``CG`` into ``TG``, so the sequences that still carry ``CG`` are
what flags the column.

Pure Python, no Qt; unit-tested in ``tests/test_methylation_mask.py``.
"""

from __future__ import annotations

from dataclasses import dataclass, field

MASK_CG = "CG"
MASK_CG_CHG = "CG_CHG"
MASK_LABELS = {MASK_CG: "CG", MASK_CG_CHG: "CG and CHG"}

_GAPS = set("-.")


def methylation_columns(aligned: list[str], contexts: str = MASK_CG_CHG) -> list[int]:
    """Alignment columns holding a C or G in a CG (or CHG) context in any sequence."""
    chg = contexts == MASK_CG_CHG
    flagged: set[int] = set()
    for row in aligned:
        cols = [i for i, c in enumerate(row) if c not in _GAPS]
        seq = "".join(row[i] for i in cols).upper()
        n = len(seq)
        for p, base in enumerate(seq):
            if base == "C":
                nxt = seq[p + 1] if p + 1 < n else ""
                if nxt == "G" or (chg and nxt in ("A", "C", "T") and p + 2 < n and seq[p + 2] == "G"):
                    flagged.add(cols[p])
            elif base == "G":
                prv = seq[p - 1] if p >= 1 else ""
                if prv == "C" or (chg and prv in ("A", "G", "T") and p >= 2 and seq[p - 2] == "C"):
                    flagged.add(cols[p])
    return sorted(flagged)


def remove_columns(aligned: list[str], columns: list[int] | set[int]) -> list[str]:
    drop = set(columns)
    return ["".join(c for i, c in enumerate(row) if i not in drop) for row in aligned]


# --------------------------------------------------------------------------
# Comparing the full and the masked tree
# --------------------------------------------------------------------------


def _numbers(label: str) -> list[float]:
    out = []
    for part in str(label).split("/"):
        try:
            out.append(float(part))
        except ValueError:
            return []
    return out


def is_well_supported(label: str | None, engine: str) -> bool:
    """Support label passes the usual threshold for its method.

    IQ-TREE ``SH-aLRT/UFBoot``: 80 and 95. A single IQ-TREE value cannot be
    told apart (SH-aLRT or UFBoot), so the stricter 95 applies. FastME
    bootstrap: 70. FastTree local support: 0.9.
    """
    if not label:
        return False
    vals = _numbers(label)
    if not vals:
        return False
    if engine == "fasttree":
        return vals[0] >= 0.9
    if engine == "fastme":
        return vals[0] >= 70
    if len(vals) >= 2:
        return vals[0] >= 80 and vals[1] >= 95
    return vals[0] >= 95


def _labelled_splits(tree) -> dict[frozenset[str], str]:
    taxa = sorted(t.name for t in tree.get_terminals())
    if len(taxa) < 4:
        return {}
    ref, all_taxa = taxa[0], frozenset(taxa)
    out: dict[frozenset[str], str] = {}
    for clade in tree.get_nonterminals():
        side = frozenset(t.name for t in clade.get_terminals())
        if ref in side:
            side = all_taxa - side
        if 1 < len(side) < len(taxa) - 1:
            if clade.name:
                label = str(clade.name)
            elif clade.confidence is not None:
                label = f"{clade.confidence:g}"
            else:
                label = ""
            if side not in out or (label and not out[side]):
                out[side] = label
    return out


@dataclass
class TopologyComparison:
    n_taxa: int
    rf: int                      # Robinson–Foulds distance (splits in one tree only)
    max_rf: int                  # 2 × (n − 3) for binary unrooted trees
    supported_full: int          # well-supported splits in the full tree
    lost_supported: list[list[str]] = field(default_factory=list)  # absent after masking

    @property
    def normalized_rf(self) -> float:
        return self.rf / self.max_rf if self.max_rf else 0.0

    @property
    def kept_supported(self) -> int:
        return self.supported_full - len(self.lost_supported)


def compare_trees(full_tree, masked_tree, engine: str) -> TopologyComparison:
    """Robinson–Foulds distance and well-supported groupings lost after masking."""
    full = _labelled_splits(full_tree)
    masked = _labelled_splits(masked_tree)
    n = len(full_tree.get_terminals())
    rf = len(set(full) ^ set(masked))
    supported = [s for s, lab in full.items() if is_well_supported(lab, engine)]
    all_taxa = frozenset(t.name for t in full_tree.get_terminals())
    lost = []
    for s in supported:
        if s not in masked:
            other = all_taxa - s
            side = s if len(s) <= len(other) else other
            lost.append(sorted(side))
    lost.sort(key=lambda names: (len(names), names))
    return TopologyComparison(n, rf, max(0, 2 * (n - 3)), len(supported), lost)


def comparison_sentence(cmp: TopologyComparison, masked_cols: int, total_cols: int,
                        contexts: str) -> str:
    """One paragraph for the report and methods text."""
    what = MASK_LABELS.get(contexts, contexts)
    share = masked_cols / total_cols if total_cols else 0.0
    head = (f"As a sensitivity test, alignment columns with a C or G in a {what} context "
            f"in any sequence were removed ({masked_cols} of {total_cols} columns, {share:.0%}) "
            "and the tree was rebuilt with the same settings.")
    if cmp.max_rf == 0:
        return head + " The tree has too few sequences for a topology comparison."
    body = (f" The normalised Robinson–Foulds distance between the full and the masked tree "
            f"was {cmp.normalized_rf:.2f} ({cmp.rf} of {cmp.max_rf} splits differ).")
    if cmp.supported_full:
        body += (f" {cmp.kept_supported} of {cmp.supported_full} well-supported groupings "
                 "of the full tree were recovered")
        body += (", so the differences concern only branches without strong support."
                 if not cmp.lost_supported and cmp.rf else ".")
    return head + body
