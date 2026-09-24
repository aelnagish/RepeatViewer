"""Load/save sequence annotations for the in-app sequence map viewer."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path

from Bio import SeqIO

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
        return cls(
            start=int(data.get("start", 1)),
            end=int(data.get("end", 1)),
            name=str(data.get("name", "")),
            feature_type=str(data.get("feature_type", "others")),
            color=str(data.get("color", "")),
            strand=str(data.get("strand", "+")),
            track=int(data.get("track", 0)),
            parent=str(data.get("parent", "") or ""),
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
            )
        )
    return updated


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
    canonical = sidecar_path(fasta_path)
    if canonical.is_file():
        return canonical
    for candidate in (
        legacy_stem_sidecar_path(fasta_path),
        legacy_sidecar_path(fasta_path),
    ):
        if not candidate.is_file():
            continue
        # Migrate legacy name to ``*.seqmap_RV.json`` when possible.
        try:
            candidate.replace(canonical)
            return canonical
        except OSError:
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
    """Return JSON sidecars that should be removed with an RV/features GFF3."""
    resolved = path.expanduser().resolve()
    if not resolved.is_file():
        return []
    name = resolved.name
    folder = resolved.parent
    stems: list[str] = []
    if name.endswith(".RV.gff3"):
        stems.append(name[: -len(".RV.gff3")])
    elif name.endswith(".features.gff3"):
        stem = name[: -len(".features.gff3")]
        stems.append(stem)
        # Legacy ``demo.fa.features.gff3`` → also try ``demo.seqmap*``.
        if stem.endswith(".fa") or stem.endswith(".fasta") or stem.endswith(".fna"):
            stems.append(Path(stem).stem)
    else:
        return []

    companions: list[Path] = []
    seen: set[str] = set()
    for stem in stems:
        for candidate in (
            folder / f"{stem}.seqmap_RV.json",
            folder / f"{stem}.seqmap.json",
        ):
            key = str(candidate)
            if key in seen:
                continue
            seen.add(key)
            if candidate.is_file() and candidate != resolved:
                companions.append(candidate)
    return companions


def companion_annotation_renames_for_fasta(
    old_fasta: Path,
    new_fasta: Path,
) -> list[tuple[Path, Path]]:
    """Map existing annotation sidecars of ``old_fasta`` onto ``new_fasta`` names.

    Renames ``{stem}.RV.gff3`` / legacy ``.features.gff3`` and ``*.seqmap_RV.json``
    (plus legacy JSON names) so they stay paired after a FASTA rename. Legacy
    companion names are migrated to the canonical ``.RV.gff3`` /
    ``.seqmap_RV.json`` forms when possible.
    """
    old = old_fasta.expanduser().resolve()
    new = new_fasta.expanduser().resolve()
    if old == new:
        return []

    planned: list[tuple[Path, Path]] = []
    seen_src: set[str] = set()
    for src, dst in (
        (features_gff_path(old), features_gff_path(new)),
        (legacy_stem_features_gff_path(old), features_gff_path(new)),
        (legacy_features_gff_path(old), features_gff_path(new)),
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
    """Return FASTA record ids in file order (sequences are not kept in memory)."""
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


def fasta_record_lengths(path: Path) -> dict[str, int]:
    """Map each FASTA record id to its sequence length."""
    resolved = path.expanduser().resolve()
    if not resolved.is_file():
        raise FileNotFoundError(resolved)
    lengths: dict[str, int] = {}
    for record in SeqIO.parse(str(resolved), "fasta"):
        rid = str(record.id).strip()
        if rid:
            lengths[rid] = len(record.seq)
    if not lengths:
        raise ValueError(f"No sequences in {resolved.name}")
    return lengths


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
            )
        )
    return merged


def load_sequence_map_data(
    path: Path,
    *,
    active_id: str | None = None,
) -> SequenceMapData:
    """Load sequence-map data for a single- or multi-record FASTA.

    Only the active sequence string is held in memory; annotations for all
    record ids are loaded from the companion GFF when present.
    """
    resolved = path.expanduser().resolve()
    lengths = fasta_record_lengths(resolved)
    record_ids = list(lengths.keys())

    gff_path = resolve_features_gff_path(resolved)
    features_by_id = (
        load_all_features_by_id(gff_path, record_ids, lengths_by_id=lengths)
        if gff_path is not None
        else {rid: [] for rid in record_ids}
    )

    # JSON color overlay (single-id legacy or multi-id payload).
    json_by_id = load_features_by_id_from_json(resolved, record_ids)
    if json_by_id:
        for rid, json_feats in json_by_id.items():
            if not json_feats:
                continue
            current = features_by_id.get(rid) or []
            if current:
                features_by_id[rid] = _merge_feature_colors_from_json(current, json_feats)
            else:
                features_by_id[rid] = json_feats

    chosen = str(active_id or "").strip()
    if chosen and chosen in lengths:
        seq_id = chosen
    else:
        seq_id = next(
            (rid for rid in record_ids if features_by_id.get(rid)),
            record_ids[0],
        )
    seq_id, sequence = load_fasta_record(resolved, seq_id)
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


def load_features_by_id_from_json(
    path: Path,
    record_ids: list[str],
) -> dict[str, list[SequenceFeature]]:
    """Load JSON sidecar features keyed by sequence id (legacy or multi)."""
    sidecar = resolve_sidecar_path(path)
    if sidecar is None:
        return {}
    try:
        payload = json.loads(sidecar.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}

    result: dict[str, list[SequenceFeature]] = {rid: [] for rid in record_ids}

    sequences = payload.get("sequences")
    if isinstance(sequences, dict):
        for rid in record_ids:
            block = sequences.get(rid)
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
            )
        )
    return sorted(features, key=lambda f: (f.start, f.end))


def save_features(path: Path, *, sequence_id: str, features: list[SequenceFeature]) -> Path:
    sidecar = sidecar_path(path)
    payload = {
        "sequence_id": sequence_id,
        "features": [f.to_dict() for f in features],
    }
    sidecar.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return sidecar


def save_features_multi(
    path: Path,
    features_by_id: dict[str, list[SequenceFeature]],
) -> Path:
    """Write a multi-sequence JSON color/edit cache beside the FASTA."""
    sidecar = sidecar_path(path)
    payload = {
        "sequences": {
            seq_id: {"features": [f.to_dict() for f in feats]}
            for seq_id, feats in features_by_id.items()
            if feats
        }
    }
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
    return f"{sequence_id}\tRepeatViewer\t{safe_type}\t{start}\t{end}\t.\t{strand}\t.\t{attrs}"


def export_features_gff_multi(
    path: Path,
    *,
    features_by_id: dict[str, list[SequenceFeature]],
    lengths_by_id: dict[str, int],
) -> Path:
    """Write ``{stem}.RV.gff3`` with features for one or more sequence ids."""
    resolved = path.expanduser().resolve()
    out = features_gff_path(resolved)
    lines = ["##gff-version 3"]
    # Prefer a stable order: lengths_by_id insertion order, then any extras.
    ordered_ids = [sid for sid in lengths_by_id if sid in features_by_id or sid in lengths_by_id]
    seen = set(ordered_ids)
    for sid in features_by_id:
        if sid not in seen:
            ordered_ids.append(sid)
            seen.add(sid)
    for seq_id in ordered_ids:
        length = max(1, int(lengths_by_id.get(seq_id, 1)))
        lines.append(f"##sequence-region {seq_id} 1 {length}")
    for seq_id in ordered_ids:
        feats = features_by_id.get(seq_id) or []
        if not feats:
            continue
        length = max(1, int(lengths_by_id.get(seq_id, 10**9)))
        for feat in feats:
            lines.append(_gff_feature_line(seq_id, feat, length))
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out


def save_annotations(
    path: Path,
    *,
    sequence_id: str,
    sequence_length: int,
    features: list[SequenceFeature],
) -> tuple[Path, Path | None]:
    """Persist JSON sidecar and companion GFF3 when features are present."""
    return save_annotations_multi(
        path,
        features_by_id={sequence_id: features},
        lengths_by_id={sequence_id: sequence_length},
    )


def save_annotations_multi(
    path: Path,
    *,
    features_by_id: dict[str, list[SequenceFeature]],
    lengths_by_id: dict[str, int],
) -> tuple[Path, Path | None]:
    """Persist multi-sequence annotations (JSON cache + ``.RV.gff3``)."""
    non_empty = {sid: list(feats) for sid, feats in features_by_id.items() if feats}
    if len(lengths_by_id) <= 1 and len(non_empty) <= 1:
        # Keep legacy single-id JSON shape for one-record FASTA.
        if non_empty:
            sid, feats = next(iter(non_empty.items()))
            json_path = save_features(path, sequence_id=sid, features=feats)
        else:
            sid = next(iter(lengths_by_id), "")
            json_path = save_features(path, sequence_id=sid, features=[])
    else:
        json_path = save_features_multi(path, non_empty)

    gff_path = features_gff_path(path)
    legacy_gff = legacy_features_gff_path(path)
    if non_empty:
        gff_path = export_features_gff_multi(
            path,
            features_by_id=non_empty,
            lengths_by_id=lengths_by_id,
        )
        _cleanup_legacy_sidecars(path)
    else:
        for candidate in (
            gff_path,
            legacy_gff,
            legacy_stem_features_gff_path(path),
        ):
            _remove_sidecar_if_exists(candidate)
        gff_path = None
        _cleanup_legacy_sidecars(path)
    return json_path, gff_path


def _gff_rows_to_sequence_features(
    rows: list,
    *,
    sequence_id: str,
    sequence_length: int,
    allow_remap: bool = True,
) -> list[SequenceFeature]:
    """Map parsed GFF/BED rows onto SequenceFeature objects for one FASTA id."""
    from repeatviewer.backends.gff_features import _seqid_matches

    matched = [row for row in rows if _seqid_matches(sequence_id, row.seqid)]
    # Single-consensus TE files often use a different seqid than the FASTA
    # header; remap all rows when nothing matched but rows exist.
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
            )
        )
    return sorted(features, key=lambda f: (f.start, f.end))


@dataclass(frozen=True)
class AnnotationImportResult:
    gff_path: Path
    feature_count: int
    sequence_count: int

    def __iter__(self):
        """Allow ``path, count = import_annotation_as_rv_gff3(...)``."""
        yield self.gff_path
        yield self.feature_count


def import_annotation_as_rv_gff3(
    fasta_path: Path,
    annotation_path: Path,
    *,
    merge_existing: bool = True,
) -> AnnotationImportResult:
    """Convert GFF/GFF3/GTF/BED into ``{stem}.RV.gff3`` beside the FASTA.

    Features are matched to FASTA records by sequence id. For a one-record
    FASTA, unmatched GFF seqids are remapped onto that record (legacy TE
    consensus behaviour). For multi-record FASTA, remapping is disabled.
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
    multi = len(record_ids) > 1
    rows = parse_annotation(source)
    if not rows:
        raise ValueError(f"No features found in {source.name}")

    features_by_id: dict[str, list[SequenceFeature]] = {rid: [] for rid in record_ids}
    if multi:
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
    else:
        rid = record_ids[0]
        features_by_id[rid] = _gff_rows_to_sequence_features(
            rows,
            sequence_id=rid,
            sequence_length=lengths[rid],
            allow_remap=True,
        )

    imported_total = sum(len(feats) for feats in features_by_id.values())
    if imported_total == 0:
        raise ValueError(
            f"No features in {source.name} matched any sequence id in {fasta.name}."
        )

    if merge_existing:
        existing_path = resolve_features_gff_path(fasta)
        if existing_path is not None:
            existing_by_id = load_all_features_by_id(
                existing_path,
                record_ids,
                lengths_by_id=lengths,
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

    _json_path, gff_path = save_annotations_multi(
        fasta,
        features_by_id=features_by_id,
        lengths_by_id=lengths,
    )
    if gff_path is None:
        raise RuntimeError("Failed to write RepeatViewer annotation GFF3.")
    feature_count = sum(len(feats) for feats in features_by_id.values())
    sequence_count = sum(1 for feats in features_by_id.values() if feats)
    return AnnotationImportResult(
        gff_path=gff_path,
        feature_count=feature_count,
        sequence_count=sequence_count,
    )
