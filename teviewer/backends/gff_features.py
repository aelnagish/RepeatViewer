"""Parse GFF/GTF and build a simple sequence feature map for preview."""

from __future__ import annotations

import html
from dataclasses import dataclass
from pathlib import Path

_ANNOTATION_SUFFIXES = (".gff", ".gff3", ".gtf", ".bed")
_GFF_SUFFIXES = (".gff", ".gff3", ".gtf")


def is_gene_annotation_file(path: Path) -> bool:
    """True for GFF/GTF/BED annotation files usable as genome-map feature input."""
    name = path.name.lower()
    return any(name.endswith(suffix) for suffix in _ANNOTATION_SUFFIXES)


@dataclass(frozen=True)
class GffFeature:
    seqid: str
    source: str
    feature_type: str
    start: int
    end: int
    strand: str
    score: str
    phase: str
    attributes: dict[str, str]

    @property
    def length(self) -> int:
        return max(0, self.end - self.start + 1)

    def label(self) -> str:
        for key in ("Name", "ID", "gene", "locus_tag", "product"):
            if key in self.attributes:
                return self.attributes[key]
        return self.feature_type


def find_companion_gff(sequence_path: Path) -> Path | None:
    """Look for a GFF/GTF next to the sequence file (same stem or common names)."""
    folder = sequence_path.parent
    stem = sequence_path.stem
    candidates = [
        folder / f"{stem}.features.gff3",
        folder / f"{sequence_path.name}.features.gff3",
        folder / f"{stem}.gff3",
        folder / f"{stem}.gff",
        folder / f"{stem}.gtf",
        folder / f"{stem}.annotations.gff3",
        folder / f"{stem}.annotations.gff",
    ]
    for path in candidates:
        if path.is_file():
            return path
    for suffix in _GFF_SUFFIXES:
        for path in sorted(folder.glob(f"*{suffix}")):
            if path.is_file():
                return path
    return None


def fasta_record_ids(fasta_path: Path, *, limit: int = 200) -> list[str]:
    ids: list[str] = []
    try:
        with fasta_path.open(encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if not line.startswith(">"):
                    continue
                header = line[1:].strip()
                if not header:
                    continue
                record_id = header.split()[0]
                ids.append(record_id)
                if len(ids) >= limit:
                    break
    except OSError:
        return []
    return ids


def _parse_attributes(raw: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    if not raw or raw == ".":
        return attrs
    for part in raw.split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            key, value = part.split("=", 1)
        else:
            key, value = part, ""
        attrs[key.strip()] = value.strip()
    return attrs


def parse_gff(path: Path) -> list[GffFeature]:
    features: list[GffFeature] = []
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return features
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        cols = line.split("\t")
        if len(cols) < 9:
            continue
        try:
            start = int(cols[3])
            end = int(cols[4])
        except ValueError:
            continue
        features.append(
            GffFeature(
                seqid=cols[0],
                source=cols[1],
                feature_type=cols[2],
                start=start,
                end=end,
                strand=cols[6] if len(cols) > 6 else ".",
                score=cols[5] if len(cols) > 5 else ".",
                phase=cols[7] if len(cols) > 7 else ".",
                attributes=_parse_attributes(cols[8]),
            )
        )
    return features


def parse_bed(path: Path) -> list[GffFeature]:
    """Parse BED intervals (0-based start, exclusive end) into GffFeature rows."""
    features: list[GffFeature] = []
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return features
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        lower = line.lower()
        if lower.startswith("track") or lower.startswith("browser"):
            continue
        cols = line.split("\t")
        if len(cols) < 3:
            continue
        try:
            start0 = int(cols[1])
            end0 = int(cols[2])
        except ValueError:
            continue
        if end0 <= start0:
            continue
        name = cols[3].strip() if len(cols) > 3 else "."
        score = cols[4].strip() if len(cols) > 4 else "."
        strand = cols[5].strip() if len(cols) > 5 else "."
        if strand not in {"+", "-"}:
            strand = "."
        attrs: dict[str, str] = {}
        if name and name != ".":
            attrs["Name"] = name
        features.append(
            GffFeature(
                seqid=cols[0],
                source="bed",
                feature_type="gene",
                start=start0 + 1,
                end=end0,
                strand=strand,
                score=score,
                phase=".",
                attributes=attrs,
            )
        )
    return features


def parse_annotation(path: Path) -> list[GffFeature]:
    """Parse GFF/GTF (all rows) or BED (interval rows) from an annotation file."""
    if path.name.lower().endswith(".bed"):
        return parse_bed(path)
    return parse_gff(path)


def list_feature_types_in_annotation(path: Path, *, max_lines: int = 200_000) -> list[str]:
    """Return sorted unique feature types in a GFF/GTF/BED annotation file."""
    path = path.expanduser().resolve()
    if not path.is_file():
        return []
    if path.name.lower().endswith(".bed"):
        return ["interval"]
    found: set[str] = set()
    try:
        with path.open(encoding="utf-8", errors="replace") as handle:
            for line_no, line in enumerate(handle):
                if line_no >= max_lines:
                    break
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                cols = line.split("\t")
                if len(cols) < 3:
                    continue
                feature_type = cols[2].strip()
                if feature_type and feature_type != ".":
                    found.add(feature_type)
    except OSError:
        return []
    return sorted(found, key=str.lower)


def _seqid_matches(record_id: str, gff_seqid: str) -> bool:
    if record_id == gff_seqid:
        return True
    if record_id.startswith(gff_seqid) or gff_seqid.startswith(record_id):
        return True
    base = record_id.split()[0]
    return base == gff_seqid or gff_seqid.startswith(base) or base.startswith(gff_seqid)


def features_for_records(
    features: list[GffFeature],
    record_ids: list[str],
) -> dict[str, list[GffFeature]]:
    by_seq: dict[str, list[GffFeature]] = {rid: [] for rid in record_ids}
    if not record_ids:
        return by_seq
    for feat in features:
        for rid in record_ids:
            if _seqid_matches(rid, feat.seqid):
                by_seq[rid].append(feat)
    return by_seq


_TYPE_COLORS = {
    "gene": "#4a90d9",
    "mrna": "#5cb85c",
    "exon": "#8bc34a",
    "cds": "#f0ad4e",
    "repeat_region": "#d9534f",
    "transposable_element": "#c9302c",
    "ltr": "#e67e22",
    "tir": "#9b59b6",
}


def _feature_color(feature_type: str) -> str:
    key = feature_type.lower().replace(" ", "_")
    if key in _TYPE_COLORS:
        return _TYPE_COLORS[key]
    for prefix, color in _TYPE_COLORS.items():
        if key.startswith(prefix):
            return color
    return "#7f8c8d"


def build_gff_annotation_html(
    fasta_path: Path,
    gff_path: Path,
    *,
    max_records: int = 12,
    max_features_per_record: int = 80,
) -> str:
    record_ids = fasta_record_ids(fasta_path, limit=max_records)
    all_features = parse_gff(gff_path)
    by_record = features_for_records(all_features, record_ids)

    sections: list[str] = []
    for rid in record_ids:
        feats = sorted(by_record.get(rid, []), key=lambda f: (f.start, f.end))
        if not feats:
            continue
        feats = feats[:max_features_per_record]
        seq_len = max(f.end for f in feats)
        if seq_len <= 0:
            continue
        bars: list[str] = []
        rows: list[str] = []
        for feat in feats:
            left = 100.0 * (feat.start - 1) / seq_len
            width = max(0.4, 100.0 * feat.length / seq_len)
            color = _feature_color(feat.feature_type)
            label = html.escape(feat.label())
            ftype = html.escape(feat.feature_type)
            bars.append(
                f'<div class="feat" style="left:{left:.2f}%;width:{width:.2f}%;'
                f'background:{color}" title="{ftype}: {feat.start}-{feat.end} ({feat.strand}) '
                f'{label}"></div>'
            )
            rows.append(
                "<tr>"
                f"<td>{ftype}</td><td>{feat.start}</td><td>{feat.end}</td>"
                f"<td>{html.escape(feat.strand)}</td><td>{label}</td>"
                "</tr>"
            )
        sections.append(
            f"""
            <section class="record">
              <h3>{html.escape(rid)} <span class="meta">({len(feats)} features, ~{seq_len:,} bp span)</span></h3>
              <div class="track">{"".join(bars)}</div>
              <table>
                <thead><tr><th>Type</th><th>Start</th><th>End</th><th>Strand</th><th>Label</th></tr></thead>
                <tbody>{"".join(rows)}</tbody>
              </table>
            </section>
            """
        )

    if not sections:
        unmatched = len({f.seqid for f in all_features})
        body = (
            f"<p>No GFF features matched FASTA record IDs in <code>{html.escape(fasta_path.name)}</code>.</p>"
            f"<p>GFF file <code>{html.escape(gff_path.name)}</code> has {len(all_features)} features "
            f"on {unmatched} sequence(s). Check that GFF <code>seqid</code> columns match FASTA headers.</p>"
        )
    else:
        body = "".join(sections)

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>GFF annotations — {html.escape(fasta_path.name)}</title>
  <style>
    body {{ font-family: system-ui, sans-serif; margin: 1rem 1.25rem; color: #222; }}
    h1 {{ font-size: 1.2rem; margin-bottom: 0.25rem; }}
    .sub {{ color: #666; font-size: 0.9rem; margin-bottom: 1.25rem; }}
    .record {{ margin-bottom: 2rem; }}
    .record h3 {{ font-size: 1rem; margin: 0 0 0.5rem; }}
    .meta {{ color: #888; font-weight: normal; font-size: 0.85rem; }}
    .track {{
      position: relative; height: 28px; background: #f4f4f4; border: 1px solid #ddd;
      border-radius: 4px; margin-bottom: 0.5rem;
    }}
    .feat {{
      position: absolute; top: 4px; height: 20px; border-radius: 2px; opacity: 0.92;
      min-width: 2px;
    }}
    table {{ border-collapse: collapse; width: 100%; font-size: 0.85rem; }}
    th, td {{ border: 1px solid #e0e0e0; padding: 4px 8px; text-align: left; }}
    th {{ background: #fafafa; }}
    code {{ background: #f0f0f0; padding: 0 4px; border-radius: 3px; }}
  </style>
</head>
<body>
  <h1>GFF feature map</h1>
  <p class="sub">Sequence: <code>{html.escape(fasta_path.name)}</code> ·
  Annotations: <code>{html.escape(gff_path.name)}</code></p>
  {body}
  <p class="sub">Edit the GFF on disk and reopen this view to refresh. For full editing use IGV, Apollo, or GenomeView.</p>
</body>
</html>
"""
