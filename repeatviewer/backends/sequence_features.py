"""Load/save sequence annotations for the in-app sequence map viewer."""

from __future__ import annotations

import json
import re
from collections.abc import Iterable
from dataclasses import asdict, dataclass, field
from pathlib import Path

from Bio import SeqIO

# Trailing #N (N>=2) added by uniquify_record_ids for duplicate FASTA headers.
_UNIQUIFY_SUFFIX_RE = re.compile(r"^(.+)#([2-9]|[1-9]\d+)$")


def fasta_id_token(header_or_id: str) -> str:
    """Text after ``>`` until the first whitespace (FASTA / GFF seqid token)."""
    text = str(header_or_id or "").strip()
    if not text:
        return ""
    if text.startswith(">"):
        text = text[1:].strip()
    return text.split()[0] if text else ""


def strip_duplicate_uniquifier(
    seq_id: str,
    *,
    sibling_ids: Iterable[str] | None = None,
) -> str:
    """Remove trailing ``#N`` (N≥2) only when ``uniquify_record_ids`` added it.

    Literal ``#`` characters in a real FASTA id (e.g. ``Class_I#LINE``) are kept.
    When ``sibling_ids`` includes the base name, ``name#2`` → ``name``.
    """
    text = fasta_id_token(seq_id)
    match = _UNIQUIFY_SUFFIX_RE.fullmatch(text)
    if not match:
        return text
    base = match.group(1)
    if sibling_ids is None:
        return text
    siblings = {fasta_id_token(sid) for sid in sibling_ids}
    if base in siblings:
        return base
    return text


def export_seqid(seq_id: str, *, seqid_map: dict[str, str] | None = None) -> str:
    """Seqid written to GFF/GTF/BED: full FASTA token (before first space)."""
    if seqid_map and seq_id in seqid_map:
        return fasta_id_token(seqid_map[seq_id]) or seqid_map[seq_id]
    return fasta_id_token(seq_id) or str(seq_id)

FEATURE_TYPE_COLORS: dict[str, str] = {
    "LTR": "#c0392b",
    "TIR": "#8e44ad",
    "PBS": "#d35400",
    "PPT": "#e67e22",
    "ORF": "#27ae60",
    "tandem repeat": "#000000",
    "domain": "#2980b9",
    "others": "#7f8c8d",
    # Legacy labels from older annotations
    "Domain": "#2980b9",
    "Other": "#7f8c8d",
}

DEFAULT_FEATURE_TYPES = (
    "LTR",
    "TIR",
    "PBS",
    "PPT",
    "ORF",
    "tandem repeat",
    "domain",
    "others",
)


@dataclass
class SequenceFeature:
    """1-based inclusive interval on the sequence."""

    start: int
    end: int
    name: str
    feature_type: str = "others"
    color: str = ""
    strand: str = "+"
    track: int = 0
    # GFF3 Parent attribute (e.g. exon → gene ID). Empty when top-level.
    parent: str = ""
    # Where the feature came from (import filename, or repeatviewer_internal).
    source: str = ""
    # BLAST percent identity (0 = unset / legacy). Used for Online BLAST coloring.
    pident: float = 0.0
    # Normalized subject accession for multi-HSP grouping (Online BLAST).
    subject_id: str = ""

    def __post_init__(self) -> None:
        self.start = int(self.start)
        self.end = int(self.end)
        if self.start > self.end:
            self.start, self.end = self.end, self.start
        if not self.color:
            try:
                from repeatviewer.backends.annotation_colors import feature_type_color

                self.color = feature_type_color(self.feature_type)
            except ImportError:
                self.color = FEATURE_TYPE_COLORS.get(
                    self.feature_type,
                    FEATURE_TYPE_COLORS["others"],
                )
        self.strand = self._normalize_strand(self.strand)
        self.track = max(0, int(self.track))
        self.parent = str(self.parent or "").strip()
        self.source = str(self.source or "").strip()
        try:
            self.pident = float(self.pident or 0.0)
        except (TypeError, ValueError):
            self.pident = 0.0
        self.subject_id = str(self.subject_id or "").strip()

    @staticmethod
    def _normalize_strand(strand: str) -> str:
        s = str(strand).strip().lower()
        if s in {"-", "minus", "reverse", "rev", "c"}:
            return "-"
        return "+"

    def is_forward(self) -> bool:
        return self.strand != "-"

    def gff_id(self) -> str:
        """Stable GFF3 ID used for Parent links (Name/ID in our exports)."""
        return str(self.name or "").strip()

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict) -> SequenceFeature:
        try:
            pident = float(data.get("pident", 0.0) or 0.0)
        except (TypeError, ValueError):
            pident = 0.0
        return cls(
            start=int(data.get("start", 1)),
            end=int(data.get("end", 1)),
            name=str(data.get("name", "")),
            feature_type=str(data.get("feature_type", "others")),
            color=str(data.get("color", "")),
            strand=str(data.get("strand", "+")),
            track=int(data.get("track", 0)),
            parent=str(data.get("parent", "") or ""),
            source=str(data.get("source", "") or ""),
            pident=pident,
            subject_id=str(data.get("subject_id", "") or ""),
        )


def remount_parent_ids(
    features: list[SequenceFeature],
    *,
    old_id: str,
    new_id: str,
) -> list[SequenceFeature]:
    """Retarget ``Parent`` attributes when a feature's GFF ID (name) changes."""
    old = str(old_id or "").strip()
    new = str(new_id or "").strip()
    if not old or old == new:
        return features
    updated: list[SequenceFeature] = []
    for feat in features:
        if feat.parent != old:
            updated.append(feat)
            continue
        updated.append(
            SequenceFeature(
                start=feat.start,
                end=feat.end,
                name=feat.name,
                feature_type=feat.feature_type,
                color=feat.color,
                strand=feat.strand,
                track=feat.track,
                parent=new,
                source=feat.source,
                pident=feat.pident,
                subject_id=feat.subject_id,
            )
        )
    return updated


# Features from Sequence Annotation or manual curation share this source.
# Imported annotations use the source file name instead.
REPEATVIEWER_INTERNAL_SOURCE = "repeatviewer_internal"
# Back-compat alias for annotation pipelines.
SEQUENCE_ANNOTATION_SOURCE = REPEATVIEWER_INTERNAL_SOURCE


def feature_display_pack_key(
    feat: SequenceFeature,
    *,
    split_by_source: bool,
) -> tuple[str, int]:
    """Layout key for annotation rows: ``(source, track)`` or track-only."""
    source = str(feat.source or "") if split_by_source else ""
    return (source, int(feat.track))


def ordered_feature_pack_keys(
    features: list[SequenceFeature],
    *,
    split_by_source: bool,
) -> list[tuple[str, int]]:
    """Empty/legacy source first, then sources alphabetically; tracks ascending."""
    keys = {
        feature_display_pack_key(feat, split_by_source=split_by_source)
        for feat in features
    }
    return sorted(keys, key=lambda item: (0 if item[0] == "" else 1, item[0].lower(), item[1]))


def with_feature_source(
    feat: SequenceFeature,
    source: str,
) -> SequenceFeature:
    """Return a copy of ``feat`` with ``source`` set (no-op when unchanged)."""
    label = str(source or "").strip()
    if feat.source == label:
        return feat
    return SequenceFeature(
        start=feat.start,
        end=feat.end,
        name=feat.name,
        feature_type=feat.feature_type,
        color=feat.color,
        strand=feat.strand,
        track=feat.track,
        parent=feat.parent,
        source=label,
        pident=feat.pident,
        subject_id=feat.subject_id,
    )


def _pack_feature_row_group(features: list[SequenceFeature]) -> tuple[list[int], int]:
    """Assign stacked rows for one pack band, packing overlaps vertically."""
    if not features:
        return [], 0
    row_for_index = [0] * len(features)
    row_ends: list[int] = []
    for idx, feat in sorted(enumerate(features), key=lambda item: item[1].start):
        for row, last_end in enumerate(row_ends):
            if feat.start > last_end:
                row_ends[row] = feat.end
                row_for_index[idx] = row
                break
        else:
            row_for_index[idx] = len(row_ends)
            row_ends.append(feat.end)
    return row_for_index, max(0, len(row_ends))


def _pack_feature_display_rows_internal(
    features: list[SequenceFeature],
    *,
    split_by_source: bool = True,
) -> tuple[list[int], int]:
    """Assign display rows by ``(source, track)`` or ``track``, stacking bands."""
    if not features:
        return [], 0
    keys = ordered_feature_pack_keys(features, split_by_source=split_by_source)
    row_for_index = [0] * len(features)
    row_offset = 0
    for key in keys:
        indices = [
            idx
            for idx, feat in enumerate(features)
            if feature_display_pack_key(feat, split_by_source=split_by_source) == key
        ]
        group = [features[idx] for idx in indices]
        group_rows, group_count = _pack_feature_row_group(group)
        for local_idx, global_idx in enumerate(indices):
            row_for_index[global_idx] = group_rows[local_idx] + row_offset
        row_offset += group_count
    return row_for_index, row_offset


def _online_blast_subject_key(feat: SequenceFeature) -> str:
    sid = str(getattr(feat, "subject_id", "") or "").strip()
    if sid:
        return sid
    parent = str(getattr(feat, "parent", "") or "").strip()
    if parent:
        return parent
    return str(feat.name or feat.feature_type or "subject")


def pack_feature_display_rows(
    features: list[SequenceFeature],
    *,
    split_by_source: bool = True,
) -> tuple[list[int], int]:
    """Assign display rows; Online BLAST uses one row per subject (NCBI graphic style)."""
    if not features:
        return [], 0
    from repeatviewer.backends.ncbi_blast_online import ONLINE_BLAST_SOURCE

    non_online_idx = [
        i for i, f in enumerate(features) if str(f.source or "") != ONLINE_BLAST_SOURCE
    ]
    online_idx = [
        i for i, f in enumerate(features) if str(f.source or "") == ONLINE_BLAST_SOURCE
    ]
    row_for_index = [0] * len(features)
    row_offset = 0
    if non_online_idx:
        subset = [features[i] for i in non_online_idx]
        sub_rows, count = _pack_feature_display_rows_internal(
            subset, split_by_source=split_by_source
        )
        for local, global_idx in enumerate(non_online_idx):
            row_for_index[global_idx] = sub_rows[local]
        row_offset = count
    if online_idx:
        groups: dict[str, list[int]] = {}
        for idx in online_idx:
            key = _online_blast_subject_key(features[idx])
            groups.setdefault(key, []).append(idx)

        def _group_score(key: str) -> float:
            return max(float(features[i].pident or 0) for i in groups[key])

        # One subject band, but pack overlapping HSPs into sub-lanes so they
        # do not paint on top of each other (common for TE / multi-domain hits).
        for key in sorted(groups.keys(), key=_group_score, reverse=True):
            members = groups[key]
            member_feats = [features[i] for i in members]
            local_rows, local_count = _pack_feature_row_group(member_feats)
            for local_i, global_idx in enumerate(members):
                row_for_index[global_idx] = row_offset + local_rows[local_i]
            row_offset += max(1, local_count)
    return row_for_index, row_offset


@dataclass
class SequenceMapData:
    """Active view of one sequence from a (possibly multi-record) FASTA.

    ``record_ids`` / ``features_by_id`` hold the full map for switching;
    ``sequence_id`` / ``sequence`` / ``features`` are the currently displayed record.
    """

    fasta_path: Path
    sequence_id: str
    sequence: str
    features: list[SequenceFeature] = field(default_factory=list)
    is_protein: bool = False
    record_ids: list[str] = field(default_factory=list)
    features_by_id: dict[str, list[SequenceFeature]] = field(default_factory=dict)
    lengths_by_id: dict[str, int] = field(default_factory=dict)

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def is_multi_sequence(self) -> bool:
        return len(self.record_ids) > 1

    def sync_active_features(self, features: list[SequenceFeature]) -> None:
        """Write the canvas feature list back into ``features_by_id``."""
        self.features = list(features)
        if self.sequence_id:
            self.features_by_id[self.sequence_id] = list(features)

    def switch_to(self, seq_id: str) -> None:
        """Load another FASTA record into this view (lazy; mutates in place)."""
        rid, sequence = load_fasta_record(self.fasta_path, seq_id)
        self.sequence_id = rid
        self.sequence = sequence
        self.features = list(self.features_by_id.get(rid, []))
        self.is_protein = detect_sequence_is_protein(sequence, path=self.fasta_path)
        if rid not in self.features_by_id:
            self.features_by_id[rid] = []
        if rid not in self.lengths_by_id:
            self.lengths_by_id[rid] = len(sequence)
        else:
            self.lengths_by_id[rid] = len(sequence)


_NUCLEOTIDE_CHARS = frozenset("ACGTUN")
_AMINO_ACID_CHARS = frozenset("ACDEFGHIKLMNPQRSTVWYBXZ*")


def detect_sequence_is_protein(sequence: str, *, path: Path | None = None) -> bool:
    """Guess protein vs nucleotide from suffix and residue composition."""
    if path is not None and path.suffix.lower() == ".faa":
        return True
    upper = sequence.upper()
    residues = [ch for ch in upper if ch.isalpha() or ch == "*"]
    if not residues:
        return False
    sample = residues[: min(len(residues), 50_000)]
    n_nuc = sum(1 for ch in sample if ch in _NUCLEOTIDE_CHARS)
    n_aa = sum(1 for ch in sample if ch in _AMINO_ACID_CHARS)
    if n_nuc / max(1, len(sample)) >= 0.95:
        return False
    if n_aa / max(1, len(sample)) >= 0.8:
        return True
    return n_aa > n_nuc


def _resolved_fasta_path(fasta_path: Path) -> Path:
    return fasta_path.expanduser().resolve()


def sidecar_path(fasta_path: Path) -> Path:
    """JSON sidecar next to the FASTA: ``{stem}.seqmap_RV.json``."""
    path = _resolved_fasta_path(fasta_path)
    return path.parent / f"{path.stem}.seqmap_RV.json"


def legacy_stem_sidecar_path(fasta_path: Path) -> Path:
    """Previous JSON name ``{stem}.seqmap.json`` (pre-seqmap_RV)."""
    path = _resolved_fasta_path(fasta_path)
    return path.parent / f"{path.stem}.seqmap.json"


def legacy_sidecar_path(fasta_path: Path) -> Path:
    """Older JSON sidecar name that duplicated the ``.fa`` suffix."""
    path = _resolved_fasta_path(fasta_path)
    return path.parent / f"{path.name}.seqmap.json"


def features_gff_path(fasta_path: Path) -> Path:
    """Canonical GFF3 path for sequence-map annotations: ``{stem}.RV.gff3``."""
    path = _resolved_fasta_path(fasta_path)
    return path.parent / f"{path.stem}.RV.gff3"


def legacy_features_gff_path(fasta_path: Path) -> Path:
    """Older GFF3 name that duplicated the ``.fa`` suffix."""
    path = _resolved_fasta_path(fasta_path)
    return path.parent / f"{path.name}.features.gff3"


def legacy_stem_features_gff_path(fasta_path: Path) -> Path:
    """Previous canonical name ``{stem}.features.gff3`` (pre-.RV.gff3)."""
    path = _resolved_fasta_path(fasta_path)
    return path.parent / f"{path.stem}.features.gff3"


def resolve_sidecar_path(fasta_path: Path) -> Path | None:
    """Return an existing annotation JSON path next to ``fasta_path``, if any.

    Does not rename or write files in the source directory.
    """
    for candidate in (
        sidecar_path(fasta_path),
        legacy_stem_sidecar_path(fasta_path),
        legacy_sidecar_path(fasta_path),
    ):
        if candidate.is_file():
            return candidate
    return None


def resolve_features_gff_path(fasta_path: Path) -> Path | None:
    for candidate in (
        features_gff_path(fasta_path),
        legacy_stem_features_gff_path(fasta_path),
        legacy_features_gff_path(fasta_path),
    ):
        if candidate.is_file():
            return candidate
    return None


def _remove_sidecar_if_exists(path: Path) -> None:
    if path.is_file():
        path.unlink()


def _cleanup_legacy_sidecars(fasta_path: Path) -> None:
    """Drop pre-stem / pre-.RV sidecar names after writing the canonical files."""
    canonical_json = sidecar_path(fasta_path)
    canonical_gff = features_gff_path(fasta_path)
    for legacy_json in (
        legacy_stem_sidecar_path(fasta_path),
        legacy_sidecar_path(fasta_path),
    ):
        if legacy_json != canonical_json:
            _remove_sidecar_if_exists(legacy_json)
    for legacy_gff in (
        legacy_features_gff_path(fasta_path),
        legacy_stem_features_gff_path(fasta_path),
    ):
        if legacy_gff != canonical_gff:
            _remove_sidecar_if_exists(legacy_gff)


def is_sequence_map_hidden_file(path: Path) -> bool:
    """True for sequence-map JSON sidecars (hidden from the file list)."""
    name = path.name
    return name.endswith(".seqmap_RV.json") or name.endswith(".seqmap.json")


def companion_annotation_files_for_delete(path: Path) -> list[Path]:
    """Return paired annotation files that should be removed with ``path``.

    - Deleting a FASTA also removes its JSON sidecar(s) only.
    - Exported ``.RV.gff3`` / ``.features.gff3`` files are independent and are
      not paired for delete, copy, or rename.
    """
    resolved = path.expanduser().resolve()
    if not resolved.is_file():
        return []
    name = resolved.name
    companions: list[Path] = []
    seen: set[str] = set()

    def _add(candidate: Path) -> None:
        try:
            key = str(candidate.resolve())
        except OSError:
            key = str(candidate)
        if key in seen or key == str(resolved):
            return
        seen.add(key)
        if candidate.is_file():
            companions.append(candidate)

    # FASTA deleted → take JSON companions only (not exported GFF3).
    suffix = resolved.suffix.lower()
    if suffix in {".fa", ".fasta", ".fna", ".faa", ".fas"} or name.lower().endswith(
        ".fa.gz"
    ):
        for candidate in (
            sidecar_path(resolved),
            legacy_stem_sidecar_path(resolved),
            legacy_sidecar_path(resolved),
        ):
            _add(candidate)
        return companions

    return companions


def companion_annotation_paths_for_fasta(fasta_path: Path) -> list[Path]:
    """Existing JSON companions beside a FASTA (for move/copy/delete/rename)."""
    return companion_annotation_files_for_delete(fasta_path)


def companion_annotation_renames_for_fasta(
    old_fasta: Path,
    new_fasta: Path,
) -> list[tuple[Path, Path]]:
    """Map existing JSON sidecars of ``old_fasta`` onto ``new_fasta`` names.

    Exported ``.RV.gff3`` files are left in place; only ``*.seqmap_RV.json``
    (plus legacy JSON names) stay paired after a FASTA rename.
    """
    old = old_fasta.expanduser().resolve()
    new = new_fasta.expanduser().resolve()
    if old == new:
        return []

    planned: list[tuple[Path, Path]] = []
    seen_src: set[str] = set()
    for src, dst in (
        (sidecar_path(old), sidecar_path(new)),
        (legacy_stem_sidecar_path(old), sidecar_path(new)),
        (legacy_sidecar_path(old), sidecar_path(new)),
    ):
        if not src.is_file():
            continue
        try:
            src_resolved = src.resolve()
            dst_resolved = dst.expanduser().resolve()
        except OSError:
            continue
        src_key = str(src_resolved)
        if src_key in seen_src:
            continue
        if src_resolved == dst_resolved:
            continue
        if dst_resolved.exists():
            continue
        seen_src.add(src_key)
        planned.append((src_resolved, dst_resolved))
    return planned


def load_single_sequence_fasta(path: Path) -> tuple[str, str]:
    """Return (sequence_id, sequence) when the file contains exactly one record."""
    resolved = path.expanduser().resolve()
    if not resolved.is_file():
        raise FileNotFoundError(resolved)
    records = list(SeqIO.parse(str(resolved), "fasta"))
    if not records:
        raise ValueError(f"No sequences in {resolved.name}")
    if len(records) != 1:
        raise ValueError(
            f"{resolved.name} has {len(records)} sequences. "
            "Sequence map works on single-sequence FASTA files only."
        )
    record = records[0]
    seq = str(record.seq).upper()
    if not seq:
        raise ValueError(f"Empty sequence in {resolved.name}")
    return str(record.id), seq


def list_fasta_record_ids(path: Path) -> list[str]:
    """Return FASTA record ids in file order (sequences are not kept in memory).

    Duplicate header ids are preserved (same name may appear more than once).
    """
    resolved = path.expanduser().resolve()
    if not resolved.is_file():
        raise FileNotFoundError(resolved)
    ids: list[str] = []
    for record in SeqIO.parse(str(resolved), "fasta"):
        rid = str(record.id).strip()
        if rid:
            ids.append(rid)
    if not ids:
        raise ValueError(f"No sequences in {resolved.name}")
    return ids


def uniquify_record_ids(record_ids: list[str]) -> list[str]:
    """Make record ids unique for in-memory maps while preserving order.

    First occurrence keeps the original name; later duplicates become
    ``name#2``, ``name#3``, … so annotations for identical headers can coexist.
    """
    seen: dict[str, int] = {}
    unique: list[str] = []
    for rid in record_ids:
        base = str(rid or "").strip() or "sequence"
        count = seen.get(base, 0) + 1
        seen[base] = count
        unique.append(base if count == 1 else f"{base}#{count}")
    return unique


def duplicate_original_sequence_ids(record_ids: list[str]) -> list[str]:
    """Original ids that appear more than once (stable unique list)."""
    counts: dict[str, int] = {}
    for rid in record_ids:
        counts[rid] = counts.get(rid, 0) + 1
    return [rid for rid, n in counts.items() if n > 1]


def suggest_unique_export_ids(record_ids: list[str]) -> list[str]:
    """Suggested seqids for export when FASTA headers collide (``name``, ``name_2``, …)."""
    seen: dict[str, int] = {}
    out: list[str] = []
    siblings = list(record_ids)
    for rid in record_ids:
        original = strip_duplicate_uniquifier(rid, sibling_ids=siblings) or "sequence"
        count = seen.get(original, 0) + 1
        seen[original] = count
        out.append(original if count == 1 else f"{original}_{count}")
    return out


def fasta_record_lengths(path: Path) -> dict[str, int]:
    """Map each unique FASTA record id to its sequence length.

    Duplicate headers are uniquified (``name#2``, …) so later records are not
    dropped when building annotation maps.
    """
    lengths, _single_seq = _fasta_lengths_and_single_sequence(path)
    return lengths


def is_aligned_fasta(path: Path) -> bool:
    """True when the FASTA has ≥2 records that all share one positive length.

    Used to choose residue-row MSA view in Preview Annotation. Streams once;
    does not load sequences into memory.
    """
    resolved = Path(path).expanduser().resolve()
    if not resolved.is_file():
        return False
    count = 0
    length: int | None = None
    try:
        for record in SeqIO.parse(str(resolved), "fasta"):
            rid = str(record.id).strip()
            if not rid:
                continue
            seq_len = len(record.seq)
            if seq_len <= 0:
                return False
            if length is None:
                length = seq_len
            elif seq_len != length:
                return False
            count += 1
    except (OSError, ValueError):
        return False
    return count >= 2 and length is not None and length > 0


def load_aligned_fasta_rows(path: Path) -> tuple[list[str], list[str]]:
    """Load an aligned FASTA as ``(record_ids, uppercase_sequences)``.

    Raises ``ValueError`` when the file is empty or sequences differ in length.
    Record ids are uniquified like :func:`fasta_record_lengths`.
    """
    resolved = Path(path).expanduser().resolve()
    if not resolved.is_file():
        raise FileNotFoundError(resolved)
    raw_ids: list[str] = []
    seqs: list[str] = []
    for record in SeqIO.parse(str(resolved), "fasta"):
        rid = str(record.id).strip()
        if not rid:
            continue
        raw_ids.append(rid)
        seqs.append(str(record.seq).upper())
    if len(raw_ids) < 2:
        raise ValueError(f"Need at least two sequences in {resolved.name}")
    lengths = {len(s) for s in seqs}
    if len(lengths) != 1:
        raise ValueError(
            "Sequences have different lengths; provide an aligned FASTA (MSA)."
        )
    if lengths.pop() == 0:
        raise ValueError("Alignment has zero columns.")
    return uniquify_record_ids(raw_ids), seqs


def load_fasta_rows(path: Path) -> tuple[list[str], list[str]]:
    """Load FASTA records without requiring equal sequence lengths.

    Record ids are uniquified in the same way as the overview length loader.
    """
    resolved = Path(path).expanduser().resolve()
    if not resolved.is_file():
        raise FileNotFoundError(resolved)
    raw_ids: list[str] = []
    seqs: list[str] = []
    for record in SeqIO.parse(str(resolved), "fasta"):
        rid = str(record.id).strip()
        if not rid:
            continue
        raw_ids.append(rid)
        seqs.append(str(record.seq).upper())
    if not raw_ids:
        raise ValueError(f"No FASTA sequences found in {resolved.name}")
    return uniquify_record_ids(raw_ids), seqs


def write_fasta_rows(
    path: Path,
    record_ids: list[str],
    sequences: list[str],
    *,
    line_width: int = 60,
) -> None:
    """Overwrite ``path`` with FASTA rows (same order as ``record_ids``).

    Lengths may differ (multi-sequence). Duplicate-display ids (``name#2``, …)
    are written back as the original FASTA token when siblings still include
    the base name.
    """
    if len(record_ids) != len(sequences):
        raise ValueError("record_ids and sequences length mismatch")
    if not sequences:
        raise ValueError("No sequences to write")
    resolved = Path(path).expanduser().resolve()
    siblings = list(record_ids)
    wrap = max(0, int(line_width))
    chunks: list[str] = []
    for rid, seq in zip(record_ids, sequences, strict=True):
        out_id = strip_duplicate_uniquifier(rid, sibling_ids=siblings) or "sequence"
        chunks.append(f">{out_id}\n")
        if wrap <= 0:
            chunks.append(f"{seq}\n")
        else:
            for i in range(0, len(seq), wrap):
                chunks.append(f"{seq[i : i + wrap]}\n")
    resolved.write_text("".join(chunks), encoding="utf-8")


def write_aligned_fasta_rows(
    path: Path,
    record_ids: list[str],
    sequences: list[str],
    *,
    line_width: int = 60,
) -> None:
    """Overwrite ``path`` with aligned FASTA rows (equal lengths required)."""
    if not sequences:
        raise ValueError("No sequences to write")
    width = len(sequences[0])
    if width <= 0:
        raise ValueError("Alignment has zero columns.")
    if any(len(s) != width for s in sequences):
        raise ValueError("Sequences have different lengths; cannot write MSA.")
    write_fasta_rows(path, record_ids, sequences, line_width=line_width)


def _fasta_lengths_and_single_sequence(
    path: Path,
) -> tuple[dict[str, int], str | None]:
    """One SeqIO pass: lengths for every record; sequence only when there is one.

    Multi-record files return ``None`` for the sequence so callers load only the
    chosen record afterward (avoids holding every sequence in memory).
    """
    resolved = path.expanduser().resolve()
    if not resolved.is_file():
        raise FileNotFoundError(resolved)
    raw_ids: list[str] = []
    raw_lens: list[int] = []
    first_seq: str | None = None
    multi = False
    for record in SeqIO.parse(str(resolved), "fasta"):
        rid = str(record.id).strip()
        if not rid:
            continue
        seq = str(record.seq).upper()
        raw_ids.append(rid)
        raw_lens.append(len(seq))
        if len(raw_ids) == 1:
            first_seq = seq
        else:
            multi = True
            first_seq = None
    if not raw_ids:
        raise ValueError(f"No sequences in {resolved.name}")
    unique = uniquify_record_ids(raw_ids)
    lengths = {uid: length for uid, length in zip(unique, raw_lens, strict=True)}
    return lengths, None if multi else first_seq


def load_fasta_record(path: Path, seq_id: str) -> tuple[str, str]:
    """Load one record by id (exact or soft match). Returns (id, uppercase sequence)."""
    from repeatviewer.backends.gff_features import _seqid_matches

    resolved = path.expanduser().resolve()
    if not resolved.is_file():
        raise FileNotFoundError(resolved)
    wanted = str(seq_id or "").strip()
    if not wanted:
        raise ValueError("Sequence id is empty.")

    index = SeqIO.index(str(resolved), "fasta")
    try:
        if wanted in index:
            record = index[wanted]
            seq = str(record.seq).upper()
            if not seq:
                raise ValueError(f"Empty sequence “{wanted}” in {resolved.name}")
            return str(record.id), seq
        for key in index:
            if _seqid_matches(wanted, key) or _seqid_matches(key, wanted):
                record = index[key]
                seq = str(record.seq).upper()
                if not seq:
                    raise ValueError(f"Empty sequence “{key}” in {resolved.name}")
                return str(record.id), seq
    finally:
        index.close()
    raise ValueError(f"Sequence “{wanted}” not found in {resolved.name}")


@dataclass(frozen=True)
class SingleSequenceSummary:
    sequence_id: str
    length: int
    gc_percent: float | None = None


_DNA_LETTERS = frozenset("ACGTN")


_ASCII_UPPERCASE = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"


def _gc_percent_for_dna_sequence(sequence: str) -> float | None:
    """Return G+C percent on unambiguous A/C/G/T bases, or None if not DNA-like.

    Counts with ``str.count`` (C speed, no allocation). The previous version
    built Python lists holding one entry per base; on a 500 Mbp scaffold that
    was ~8 GB of list and ~24 s of work on the UI thread, which hung the app.

    Note: bases are counted as A-Z rather than by ``str.isalpha()``, so exotic
    non-ASCII "letters" no longer count toward the DNA-likeness ratio. For
    nucleotide data the result is identical.
    """
    upper = sequence.upper()
    counts = {letter: upper.count(letter) for letter in _ASCII_UPPERCASE}
    total_letters = sum(counts.values())
    if not total_letters:
        return None
    dna_like = sum(count for letter, count in counts.items() if letter in _DNA_LETTERS)
    if dna_like / total_letters < 0.85:
        return None
    acgt = counts["A"] + counts["C"] + counts["G"] + counts["T"]
    if not acgt:
        return None
    return 100.0 * (counts["G"] + counts["C"]) / acgt


def summarize_loaded_sequence(sequence_id: str, sequence: str) -> SingleSequenceSummary:
    """Summary for a sequence already held in memory (no re-read from disk)."""
    return SingleSequenceSummary(
        sequence_id=sequence_id,
        length=len(sequence),
        gc_percent=_gc_percent_for_dna_sequence(sequence),
    )


def summarize_single_sequence_fasta(path: Path) -> SingleSequenceSummary:
    """Length and optional GC content for a one-record FASTA file."""
    sequence_id, sequence = load_single_sequence_fasta(path)
    return summarize_loaded_sequence(sequence_id, sequence)


def format_single_sequence_summary_text(summary: SingleSequenceSummary) -> str:
    length_part, gc_part = format_single_sequence_summary_parts(summary)
    if gc_part is None:
        return length_part
    return f"{length_part} | {gc_part}"


def format_single_sequence_summary_parts(
    summary: SingleSequenceSummary,
) -> tuple[str, str | None]:
    length_text = f"{summary.length:,}"
    length_part = f"Length: {length_text} bp"
    if summary.gc_percent is None:
        return length_part, None
    return length_part, f"GC: {summary.gc_percent:.1f}%"


def format_multi_sequence_fasta_summary_text(sequence_count: int) -> str:
    return f"Sequence number: {sequence_count:,}"


def _infer_feature_color(name: str, feature_type: str) -> str:
    """Best-effort color when loading annotations without an explicit color."""
    if feature_type == "domain" and ": " in name:
        db_name, _, _ = name.partition(": ")
        db_name = db_name.strip()
        if db_name:
            try:
                from repeatviewer.backends.annotation_colors import similarity_database_color

                return similarity_database_color(db_name, 0)
            except Exception:
                pass
    try:
        from repeatviewer.backends.annotation_colors import feature_type_color

        return feature_type_color(feature_type)
    except Exception:
        return FEATURE_TYPE_COLORS.get(feature_type, FEATURE_TYPE_COLORS["others"])


def _feature_color_key(feat: SequenceFeature) -> tuple[int, int, str, str]:
    return (feat.start, feat.end, feat.name, feat.feature_type)


def _merge_feature_colors_from_json(
    features: list[SequenceFeature],
    json_features: list[SequenceFeature],
) -> list[SequenceFeature]:
    """Overlay colors from the JSON sidecar onto features loaded from GFF."""
    color_by_key = {
        _feature_color_key(feat): feat.color
        for feat in json_features
        if feat.color
    }
    if not color_by_key:
        return features
    merged: list[SequenceFeature] = []
    for feat in features:
        color = color_by_key.get(_feature_color_key(feat), feat.color)
        if color == feat.color:
            merged.append(feat)
            continue
        merged.append(
            SequenceFeature(
                start=feat.start,
                end=feat.end,
                name=feat.name,
                feature_type=feat.feature_type,
                color=color,
                strand=feat.strand,
                track=feat.track,
                parent=feat.parent,
                source=feat.source,
                pident=feat.pident,
                subject_id=feat.subject_id,
            )
        )
    return merged


def load_sequence_map_data(
    path: Path,
    *,
    active_id: str | None = None,
) -> SequenceMapData:
    """Load sequence-map data for a single- or multi-record FASTA.

    Only the active sequence string is held in memory. Annotations prefer the
    JSON sidecar (canonical store); a companion GFF3 is used when JSON is
    absent (import / legacy).

    Single-record FASTA files are parsed once (lengths + sequence together);
    multi-record files still load only the chosen sequence afterward.
    """
    resolved = path.expanduser().resolve()
    lengths, single_sequence = _fasta_lengths_and_single_sequence(resolved)
    record_ids = list(lengths.keys())

    json_by_id = load_features_by_id_from_json(resolved, record_ids)
    gff_path = resolve_features_gff_path(resolved)
    # JSON sidecar is canonical whenever it exists — even if every feature list
    # is empty (e.g. user deleted all annotations). Otherwise we would fall back
    # to a leftover `.RV.gff3` export and undelete everything.
    if resolve_sidecar_path(resolved) is not None:
        features_by_id = {rid: list(json_by_id.get(rid) or []) for rid in record_ids}
    elif gff_path is not None:
        features_by_id = load_all_features_by_id(
            gff_path, record_ids, lengths_by_id=lengths
        )
        # Optional color overlay from an empty-ish or partial JSON.
        if json_by_id:
            for rid, json_feats in json_by_id.items():
                if not json_feats:
                    continue
                current = features_by_id.get(rid) or []
                if current:
                    features_by_id[rid] = _merge_feature_colors_from_json(
                        current, json_feats
                    )
                else:
                    features_by_id[rid] = json_feats
    else:
        features_by_id = {rid: [] for rid in record_ids}

    chosen = str(active_id or "").strip()
    if chosen and chosen in lengths:
        seq_id = chosen
    else:
        seq_id = next(
            (rid for rid in record_ids if features_by_id.get(rid)),
            record_ids[0],
        )
    if single_sequence is not None and seq_id == record_ids[0]:
        sequence = single_sequence
    else:
        # load_fasta_record matches soft ids; strip duplicate #N uniquifier for lookup.
        lookup_id = strip_duplicate_uniquifier(seq_id, sibling_ids=record_ids)
        try:
            _loaded_id, sequence = load_fasta_record(resolved, lookup_id)
        except ValueError:
            # Ambiguous duplicates: load by index order via SeqIO.
            _loaded_id, sequence = _load_fasta_record_by_index(
                resolved, record_ids.index(seq_id) if seq_id in record_ids else 0
            )
    features = list(features_by_id.get(seq_id, []))
    return SequenceMapData(
        fasta_path=resolved,
        sequence_id=seq_id,
        sequence=sequence,
        features=features,
        is_protein=detect_sequence_is_protein(sequence, path=resolved),
        record_ids=record_ids,
        features_by_id=features_by_id,
        lengths_by_id=lengths,
    )


def _load_fasta_record_by_index(path: Path, index: int) -> tuple[str, str]:
    resolved = path.expanduser().resolve()
    for i, record in enumerate(SeqIO.parse(str(resolved), "fasta")):
        if i != index:
            continue
        seq = str(record.seq).upper()
        if not seq:
            raise ValueError(f"Empty sequence at index {index} in {resolved.name}")
        return str(record.id), seq
    raise ValueError(f"Sequence index {index} not found in {resolved.name}")


def load_features_by_id_from_json(
    path: Path,
    record_ids: list[str],
) -> dict[str, list[SequenceFeature]]:
    """Load JSON sidecar features keyed by sequence id (legacy or multi).

    Format version 2 stores a list of ``{record_index, sequence_id, features}``
    so two FASTA headers with the same name keep separate annotation sets.
    """
    sidecar = resolve_sidecar_path(path)
    if sidecar is None:
        return {}
    try:
        payload = json.loads(sidecar.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}

    result: dict[str, list[SequenceFeature]] = {rid: [] for rid in record_ids}

    sequences = payload.get("sequences")
    # v2: ordered list with record_index (supports duplicate sequence names).
    if isinstance(sequences, list):
        for block in sequences:
            if not isinstance(block, dict):
                continue
            raw = block.get("features")
            if not isinstance(raw, list):
                continue
            feats: list[SequenceFeature] = []
            for item in raw:
                if isinstance(item, dict):
                    try:
                        feats.append(SequenceFeature.from_dict(item))
                    except (TypeError, ValueError):
                        continue
            feats = sorted(feats, key=lambda f: (f.start, f.end))
            try:
                idx = int(block.get("record_index", -1))
            except (TypeError, ValueError):
                idx = -1
            if 0 <= idx < len(record_ids):
                result[record_ids[idx]] = feats
                continue
            sid = str(block.get("sequence_id", "") or "").strip()
            if sid and sid in result:
                result[sid] = feats
        return result

    if isinstance(sequences, dict):
        for rid in record_ids:
            # Match uniquified keys (name#2) to original name in legacy JSON.
            original = strip_duplicate_uniquifier(rid, sibling_ids=record_ids)
            block = sequences.get(rid) or sequences.get(original)
            if not isinstance(block, dict):
                continue
            raw = block.get("features")
            if not isinstance(raw, list):
                continue
            feats = []
            for item in raw:
                if isinstance(item, dict):
                    try:
                        feats.append(SequenceFeature.from_dict(item))
                    except (TypeError, ValueError):
                        continue
            result[rid] = sorted(feats, key=lambda f: (f.start, f.end))
        return result

    # Legacy single-id payload.
    expected = str(payload.get("sequence_id", "") or "")
    raw = payload.get("features")
    if not isinstance(raw, list):
        return result
    feats = []
    for item in raw:
        if isinstance(item, dict):
            try:
                feats.append(SequenceFeature.from_dict(item))
            except (TypeError, ValueError):
                continue
    feats = sorted(feats, key=lambda f: (f.start, f.end))
    if expected and expected in result:
        result[expected] = feats
    elif len(record_ids) == 1:
        result[record_ids[0]] = feats
    return result


def load_all_features_by_id(
    gff_path: Path,
    record_ids: list[str],
    *,
    lengths_by_id: dict[str, int] | None = None,
) -> dict[str, list[SequenceFeature]]:
    """Bucket GFF features by FASTA record id (no cross-id remapping)."""
    from repeatviewer.backends.gff_features import features_for_records, parse_gff

    resolved = gff_path.expanduser().resolve()
    result: dict[str, list[SequenceFeature]] = {rid: [] for rid in record_ids}
    if not resolved.is_file() or not record_ids:
        return result
    rows = parse_gff(resolved)
    by_rec = features_for_records(rows, record_ids)
    lengths = lengths_by_id or {}
    for rid, matched_rows in by_rec.items():
        if not matched_rows:
            continue
        seq_len = max(1, int(lengths.get(rid, 10**9)))
        result[rid] = _gff_rows_to_sequence_features(
            matched_rows,
            sequence_id=rid,
            sequence_length=seq_len,
            allow_remap=False,
        )
    return result


def load_features(path: Path, *, expected_id: str | None = None) -> list[SequenceFeature]:
    sidecar = resolve_sidecar_path(path)
    if sidecar is None:
        return []
    try:
        payload = json.loads(sidecar.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return []
    if expected_id and str(payload.get("sequence_id", "")) not in ("", expected_id):
        return []
    raw = payload.get("features")
    if not isinstance(raw, list):
        return []
    features: list[SequenceFeature] = []
    for item in raw:
        if isinstance(item, dict):
            try:
                features.append(SequenceFeature.from_dict(item))
            except (TypeError, ValueError):
                continue
    return sorted(features, key=lambda f: (f.start, f.end))


def load_features_from_gff(path: Path, *, expected_id: str) -> list[SequenceFeature]:
    """Load sequence-map features from a RepeatViewer or standard GFF3 file."""
    from repeatviewer.backends.gff_features import _seqid_matches, parse_gff

    resolved = path.expanduser().resolve()
    if not resolved.is_file():
        return []
    features: list[SequenceFeature] = []
    for row in parse_gff(resolved):
        if not _seqid_matches(expected_id, row.seqid):
            continue
        name = row.label()
        feat_type = row.attributes.get("type", row.feature_type)
        raw_color = str(row.attributes.get("color", "")).strip()
        color = raw_color or _infer_feature_color(name, feat_type)
        parent = str(row.attributes.get("Parent", "") or "").strip()
        try:
            pident = float(row.attributes.get("pident", 0) or 0)
        except (TypeError, ValueError):
            pident = 0.0
        features.append(
            SequenceFeature(
                start=row.start,
                end=row.end,
                name=name,
                feature_type=feat_type,
                color=color,
                strand=row.strand if row.strand in {"+", "-"} else "+",
                track=int(row.attributes.get("track", 0) or 0),
                parent=parent,
                source=str(row.attributes.get("source", "") or "").strip(),
                pident=pident,
                subject_id=str(row.attributes.get("subject_id", "") or "").strip(),
            )
        )
    return sorted(features, key=lambda f: (f.start, f.end))


def save_features(path: Path, *, sequence_id: str, features: list[SequenceFeature]) -> Path:
    sidecar = sidecar_path(path)
    payload = {
        "format_version": 2,
        "sequences": [
            {
                "record_index": 0,
                "sequence_id": fasta_id_token(sequence_id) or sequence_id,
                "features": [f.to_dict() for f in features],
            }
        ],
    }
    sidecar.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return sidecar


def save_features_multi(
    path: Path,
    features_by_id: dict[str, list[SequenceFeature]],
    *,
    record_ids: list[str] | None = None,
) -> Path:
    """Write multi-sequence JSON beside the FASTA (canonical annotation store).

    Uses format_version 2: a list of blocks with ``record_index`` so two records
    that share the same FASTA header name keep separate feature lists.
    """
    sidecar = sidecar_path(path)
    ordered = list(record_ids) if record_ids else list(features_by_id.keys())
    # Ensure every annotated id appears even if missing from record_ids.
    for sid in features_by_id:
        if sid not in ordered:
            ordered.append(sid)
    blocks = []
    for index, seq_id in enumerate(ordered):
        feats = features_by_id.get(seq_id) or []
        if not feats:
            continue
        original = strip_duplicate_uniquifier(seq_id, sibling_ids=ordered)
        blocks.append(
            {
                "record_index": index,
                "sequence_id": original,
                "features": [f.to_dict() for f in feats],
            }
        )
    payload = {"format_version": 2, "sequences": blocks}
    sidecar.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return sidecar


def export_features_gff(
    path: Path,
    *,
    sequence_id: str,
    sequence_length: int,
    features: list[SequenceFeature],
) -> Path:
    """Write a simple GFF3 next to the FASTA (single sequence id)."""
    return export_features_gff_multi(
        path,
        features_by_id={sequence_id: features},
        lengths_by_id={sequence_id: sequence_length},
    )


def _gff_feature_line(sequence_id: str, feat: SequenceFeature, sequence_length: int) -> str:
    start = max(1, min(feat.start, sequence_length))
    end = max(1, min(feat.end, sequence_length))
    if start > end:
        start, end = end, start
    safe_name = feat.name.replace("\t", " ").replace("\n", " ").replace(";", ",")
    safe_type = str(feat.feature_type).replace("\t", " ").replace("\n", " ").replace(";", ",")
    strand = "-" if not feat.is_forward() else "+"
    attrs = f"ID={safe_name};Name={safe_name};type={safe_type}"
    if feat.parent:
        safe_parent = feat.parent.replace("\t", " ").replace("\n", " ").replace(";", ",")
        attrs += f";Parent={safe_parent}"
    if feat.track:
        attrs += f";track={feat.track}"
    if feat.color:
        attrs += f";color={feat.color}"
    if feat.source:
        safe_source = feat.source.replace("\t", " ").replace("\n", " ").replace(";", ",")
        attrs += f";source={safe_source}"
    if feat.pident > 0:
        attrs += f";pident={feat.pident:.3f}"
    if feat.subject_id:
        safe_subject = feat.subject_id.replace("\t", " ").replace("\n", " ").replace(";", ",")
        attrs += f";subject_id={safe_subject}"
    return f"{sequence_id}\tRepeatViewer\t{safe_type}\t{start}\t{end}\t.\t{strand}\t.\t{attrs}"


def _ordered_export_ids(
    features_by_id: dict[str, list[SequenceFeature]],
    lengths_by_id: dict[str, int],
) -> list[str]:
    ordered_ids = [
        sid for sid in lengths_by_id if sid in features_by_id or sid in lengths_by_id
    ]
    seen = set(ordered_ids)
    for sid in features_by_id:
        if sid not in seen:
            ordered_ids.append(sid)
            seen.add(sid)
    return ordered_ids


def export_features_gff_multi(
    path: Path,
    *,
    features_by_id: dict[str, list[SequenceFeature]],
    lengths_by_id: dict[str, int],
    output_path: Path | None = None,
    seqid_map: dict[str, str] | None = None,
) -> Path:
    """Write GFF3 (default ``{stem}.RV.gff3``) with features for one or more ids."""
    resolved = path.expanduser().resolve()
    out = Path(output_path).expanduser().resolve() if output_path else features_gff_path(resolved)
    lines = ["##gff-version 3"]
    ordered_ids = _ordered_export_ids(features_by_id, lengths_by_id)
    rename = seqid_map or {}
    for seq_id in ordered_ids:
        export_id = export_seqid(seq_id, seqid_map=rename)
        length = max(1, int(lengths_by_id.get(seq_id, 1)))
        lines.append(f"##sequence-region {export_id} 1 {length}")
    for seq_id in ordered_ids:
        feats = features_by_id.get(seq_id) or []
        if not feats:
            continue
        export_id = export_seqid(seq_id, seqid_map=rename)
        length = max(1, int(lengths_by_id.get(seq_id, 10**9)))
        for feat in feats:
            lines.append(_gff_feature_line(export_id, feat, length))
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out


def export_features_gtf_multi(
    path: Path,
    *,
    features_by_id: dict[str, list[SequenceFeature]],
    lengths_by_id: dict[str, int],
    output_path: Path,
    seqid_map: dict[str, str] | None = None,
) -> Path:
    """Write a simple GTF next to / at ``output_path``."""
    out = Path(output_path).expanduser().resolve()
    rename = seqid_map or {}
    lines: list[str] = []
    for seq_id in _ordered_export_ids(features_by_id, lengths_by_id):
        feats = features_by_id.get(seq_id) or []
        if not feats:
            continue
        export_id = export_seqid(seq_id, seqid_map=rename)
        length = max(1, int(lengths_by_id.get(seq_id, 10**9)))
        for feat in feats:
            start = max(1, min(feat.start, length))
            end = max(1, min(feat.end, length))
            if start > end:
                start, end = end, start
            strand = "-" if not feat.is_forward() else "+"
            safe_name = feat.name.replace('"', "'")
            safe_type = str(feat.feature_type).replace('"', "'")
            attrs = (
                f'gene_id "{safe_name}"; transcript_id "{safe_name}"; '
                f'feature_type "{safe_type}";'
            )
            lines.append(
                f"{export_id}\tRepeatViewer\t{safe_type}\t{start}\t{end}\t.\t{strand}\t.\t{attrs}"
            )
    out.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")
    return out


def export_features_bed_multi(
    path: Path,
    *,
    features_by_id: dict[str, list[SequenceFeature]],
    lengths_by_id: dict[str, int],
    output_path: Path,
    seqid_map: dict[str, str] | None = None,
) -> Path:
    """Write BED6 (0-based start) to ``output_path``."""
    del path  # FASTA path reserved for callers / future relative naming
    out = Path(output_path).expanduser().resolve()
    rename = seqid_map or {}
    lines: list[str] = []
    for seq_id in _ordered_export_ids(features_by_id, lengths_by_id):
        feats = features_by_id.get(seq_id) or []
        if not feats:
            continue
        export_id = export_seqid(seq_id, seqid_map=rename)
        length = max(1, int(lengths_by_id.get(seq_id, 10**9)))
        for feat in feats:
            start = max(1, min(feat.start, length))
            end = max(1, min(feat.end, length))
            if start > end:
                start, end = end, start
            strand = "-" if not feat.is_forward() else "+"
            # BED is 0-based, half-open.
            lines.append(
                f"{export_id}\t{start - 1}\t{end}\t{feat.name}\t0\t{strand}"
            )
    out.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")
    return out


def save_annotations(
    path: Path,
    *,
    sequence_id: str,
    sequence_length: int,
    features: list[SequenceFeature],
) -> tuple[Path, Path | None]:
    """Persist JSON sidecar only (GFF3 is export-on-demand)."""
    return save_annotations_multi(
        path,
        features_by_id={sequence_id: features},
        lengths_by_id={sequence_id: sequence_length},
        record_ids=[sequence_id],
    )


def save_annotations_multi(
    path: Path,
    *,
    features_by_id: dict[str, list[SequenceFeature]],
    lengths_by_id: dict[str, int],
    record_ids: list[str] | None = None,
) -> tuple[Path, Path | None]:
    """Persist multi-sequence annotations to JSON only.

    Returns ``(json_path, None)``. Write GFF3/GTF/BED via the export helpers
    when the user explicitly requests them.
    """
    del lengths_by_id  # kept for call-site compatibility
    ordered = list(record_ids) if record_ids else list(features_by_id.keys())
    non_empty = {sid: list(feats) for sid, feats in features_by_id.items() if feats}
    if len(ordered) <= 1 and len(non_empty) <= 1:
        if non_empty:
            sid, feats = next(iter(non_empty.items()))
        else:
            sid = next(iter(ordered), "")
            feats = []
        json_path = save_features(path, sequence_id=sid, features=feats)
    else:
        json_path = save_features_multi(
            path, features_by_id, record_ids=ordered
        )
    _cleanup_legacy_sidecars(path)
    return json_path, None


def _gff_rows_to_sequence_features(
    rows: list,
    *,
    sequence_id: str,
    sequence_length: int,
    allow_remap: bool = False,
) -> list[SequenceFeature]:
    """Map parsed GFF/BED rows onto SequenceFeature objects for one FASTA id.

    Only rows whose seqid matches ``sequence_id`` are kept. ``allow_remap`` is
    retained for callers/tests but defaults to False so a GFF for sequence A
    cannot silently attach to FASTA B.
    """
    from repeatviewer.backends.gff_features import _seqid_matches

    matched = [row for row in rows if _seqid_matches(sequence_id, row.seqid)]
    if matched:
        use_rows = matched
    elif allow_remap:
        use_rows = list(rows)
    else:
        use_rows = []
    features: list[SequenceFeature] = []
    for row in use_rows:
        start = max(1, min(int(row.start), sequence_length))
        end = max(1, min(int(row.end), sequence_length))
        if start > end:
            start, end = end, start
        name = row.label()
        feat_type = row.attributes.get("type", row.feature_type)
        raw_color = str(row.attributes.get("color", "")).strip()
        color = raw_color or _infer_feature_color(name, feat_type)
        strand = row.strand if row.strand in {"+", "-"} else "+"
        try:
            track = int(row.attributes.get("track", 0) or 0)
        except (TypeError, ValueError):
            track = 0
        parent = str(row.attributes.get("Parent", "") or "").strip()
        try:
            pident = float(row.attributes.get("pident", 0) or 0)
        except (TypeError, ValueError):
            pident = 0.0
        features.append(
            SequenceFeature(
                start=start,
                end=end,
                name=name,
                feature_type=str(feat_type),
                color=color,
                strand=strand,
                track=track,
                parent=parent,
                source=str(row.attributes.get("source", "") or "").strip(),
                pident=pident,
                subject_id=str(row.attributes.get("subject_id", "") or "").strip(),
            )
        )
    return sorted(features, key=lambda f: (f.start, f.end))


@dataclass(frozen=True)
class AnnotationImportResult:
    json_path: Path
    feature_count: int
    sequence_count: int

    @property
    def gff_path(self) -> Path:
        """Deprecated alias — import writes JSON only; GFF is export-on-demand."""
        return self.json_path

    def __iter__(self):
        """Allow ``path, count = import_annotation_as_rv_gff3(...)``."""
        yield self.json_path
        yield self.feature_count


def import_annotation_as_rv_gff3(
    fasta_path: Path,
    annotation_path: Path,
    *,
    merge_existing: bool = True,
) -> AnnotationImportResult:
    """Import GFF/GFF3/GTF/BED into the JSON sidecar beside the FASTA.

    Does not write ``.RV.gff3``; use Export annotation for GFF/GTF/BED output.

    Features are matched to FASTA records by sequence id (FASTA header token
    before the first space ↔ GFF column 1). Unmatched rows are skipped; if
    nothing matches, import raises.
    """
    from repeatviewer.backends.gff_features import features_for_records, parse_annotation

    fasta = _resolved_fasta_path(fasta_path)
    source = annotation_path.expanduser().resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    suffix = source.suffix.lower()
    if source.name.lower().endswith(".gz"):
        raise ValueError("Gzipped annotation files are not supported for import.")
    if suffix not in {".gff", ".gff3", ".gtf", ".bed"}:
        raise ValueError(
            f"Unsupported annotation type ({suffix or 'unknown'}). "
            "Use .gff, .gff3, .gtf, or .bed."
        )

    lengths = fasta_record_lengths(fasta)
    record_ids = list(lengths.keys())
    rows = parse_annotation(source)
    if not rows:
        raise ValueError(f"No features found in {source.name}")

    features_by_id: dict[str, list[SequenceFeature]] = {rid: [] for rid in record_ids}
    by_rec = features_for_records(rows, record_ids)
    for rid, matched_rows in by_rec.items():
        if not matched_rows:
            continue
        features_by_id[rid] = _gff_rows_to_sequence_features(
            matched_rows,
            sequence_id=rid,
            sequence_length=lengths[rid],
            allow_remap=False,
        )

    imported_total = sum(len(feats) for feats in features_by_id.values())
    if imported_total == 0:
        gff_ids = sorted({str(row.seqid).strip() for row in rows if str(row.seqid).strip()})
        fasta_preview = ", ".join(record_ids[:5])
        if len(record_ids) > 5:
            fasta_preview += f", … ({len(record_ids)} total)"
        gff_preview = ", ".join(gff_ids[:8])
        if len(gff_ids) > 8:
            gff_preview += f", … ({len(gff_ids)} total)"
        raise ValueError(
            f"No features in {source.name} matched any sequence id in {fasta.name}.\n\n"
            f"FASTA id(s): {fasta_preview}\n"
            f"GFF column 1: {gff_preview}\n\n"
            "Column 1 of the annotation must match the FASTA header text after “>” "
            "and before the first space."
        )

    # Stamp every newly imported feature with the annotation file name.
    source_label = source.name
    for rid in list(features_by_id):
        features_by_id[rid] = [
            with_feature_source(feat, source_label) for feat in features_by_id[rid]
        ]

    if merge_existing:
        # Prefer JSON store; fall back to existing GFF3.
        json_existing = load_features_by_id_from_json(fasta, record_ids)
        if json_existing and any(json_existing.values()):
            existing_by_id = json_existing
        else:
            existing_path = resolve_features_gff_path(fasta)
            existing_by_id = (
                load_all_features_by_id(
                    existing_path,
                    record_ids,
                    lengths_by_id=lengths,
                )
                if existing_path is not None
                else {}
            )
        for rid in record_ids:
            imported = features_by_id.get(rid) or []
            existing = existing_by_id.get(rid) or []
            if not existing:
                continue
            seen = {_feature_color_key(f) for f in imported}
            merged = list(imported)
            for feat in existing:
                if _feature_color_key(feat) not in seen:
                    merged.append(feat)
            features_by_id[rid] = sorted(merged, key=lambda f: (f.start, f.end))

    json_path, _ = save_annotations_multi(
        fasta,
        features_by_id=features_by_id,
        lengths_by_id=lengths,
        record_ids=record_ids,
    )
    feature_count = sum(len(feats) for feats in features_by_id.values())
    sequence_count = sum(1 for feats in features_by_id.values() if feats)
    return AnnotationImportResult(
        json_path=json_path,
        feature_count=feature_count,
        sequence_count=sequence_count,
    )
