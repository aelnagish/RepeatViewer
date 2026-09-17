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
    fasta_path: Path
    sequence_id: str
    sequence: str
    features: list[SequenceFeature] = field(default_factory=list)

    @property
    def length(self) -> int:
        return len(self.sequence)


def _resolved_fasta_path(fasta_path: Path) -> Path:
    return fasta_path.expanduser().resolve()


def sidecar_path(fasta_path: Path) -> Path:
    """JSON sidecar next to the FASTA: ``{stem}.seqmap.json``."""
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
    for candidate in (sidecar_path(fasta_path), legacy_sidecar_path(fasta_path)):
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
    legacy_json = legacy_sidecar_path(fasta_path)
    if legacy_json != canonical_json:
        _remove_sidecar_if_exists(legacy_json)
    for legacy_gff in (
        legacy_features_gff_path(fasta_path),
        legacy_stem_features_gff_path(fasta_path),
    ):
        if legacy_gff != canonical_gff:
            _remove_sidecar_if_exists(legacy_gff)


def is_sequence_map_hidden_file(path: Path) -> bool:
    """True for internal JSON sidecars (hidden from the file list)."""
    return path.name.endswith(".seqmap.json")


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
            from repeatviewer.backends.annotation_colors import similarity_database_color

            return similarity_database_color(db_name, 0)
    try:
        from repeatviewer.backends.annotation_colors import feature_type_color

        return feature_type_color(feature_type)
    except ImportError:
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


def load_sequence_map_data(path: Path) -> SequenceMapData:
    seq_id, sequence = load_single_sequence_fasta(path)
    resolved = path.expanduser().resolve()
    gff_path = resolve_features_gff_path(resolved)
    json_features = load_features(resolved, expected_id=seq_id)
    if gff_path is not None:
        features = load_features_from_gff(gff_path, expected_id=seq_id)
        if json_features:
            features = _merge_feature_colors_from_json(features, json_features)
    else:
        features = json_features
    return SequenceMapData(
        fasta_path=resolved,
        sequence_id=seq_id,
        sequence=sequence,
        features=features,
    )


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


def export_features_gff(
    path: Path,
    *,
    sequence_id: str,
    sequence_length: int,
    features: list[SequenceFeature],
) -> Path:
    """Write a simple GFF3 next to the FASTA."""
    resolved = path.expanduser().resolve()
    out = features_gff_path(resolved)
    lines = ["##gff-version 3"]
    for feat in features:
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
        lines.append(
            f"{sequence_id}\tRepeatViewer\t{safe_type}\t{start}\t{end}\t.\t{strand}\t.\t{attrs}"
        )
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
    json_path = save_features(path, sequence_id=sequence_id, features=features)
    gff_path = features_gff_path(path)
    legacy_gff = legacy_features_gff_path(path)
    if features:
        gff_path = export_features_gff(
            path,
            sequence_id=sequence_id,
            sequence_length=sequence_length,
            features=features,
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
) -> list[SequenceFeature]:
    """Map parsed GFF/BED rows onto SequenceFeature objects for one FASTA id."""
    from repeatviewer.backends.gff_features import _seqid_matches

    matched = [row for row in rows if _seqid_matches(sequence_id, row.seqid)]
    # Single-consensus TE files often use a different seqid than the FASTA
    # header; remap all rows when nothing matched but rows exist.
    use_rows = matched if matched else list(rows)
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


def import_annotation_as_rv_gff3(
    fasta_path: Path,
    annotation_path: Path,
    *,
    merge_existing: bool = True,
) -> tuple[Path, int]:
    """Convert GFF/GFF3/GTF/BED into ``{stem}.RV.gff3`` beside the FASTA.

    Returns ``(written_gff_path, feature_count)``.
    """
    from repeatviewer.backends.gff_features import parse_annotation

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

    sequence_id, sequence = load_single_sequence_fasta(fasta)
    rows = parse_annotation(source)
    if not rows:
        raise ValueError(f"No features found in {source.name}")

    imported = _gff_rows_to_sequence_features(
        rows,
        sequence_id=sequence_id,
        sequence_length=len(sequence),
    )
    if not imported:
        raise ValueError(
            f"No usable features for sequence “{sequence_id}” in {source.name}"
        )

    features = list(imported)
    if merge_existing:
        existing_path = resolve_features_gff_path(fasta)
        existing = (
            load_features_from_gff(existing_path, expected_id=sequence_id)
            if existing_path is not None
            else []
        )
        if existing:
            seen = {_feature_color_key(f) for f in imported}
            for feat in existing:
                if _feature_color_key(feat) not in seen:
                    features.append(feat)
            features = sorted(features, key=lambda f: (f.start, f.end))

    _json_path, gff_path = save_annotations(
        fasta,
        sequence_id=sequence_id,
        sequence_length=len(sequence),
        features=features,
    )
    if gff_path is None:
        raise RuntimeError("Failed to write RepeatViewer annotation GFF3.")
    return gff_path, len(features)
