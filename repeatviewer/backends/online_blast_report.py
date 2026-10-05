"""Self-contained NCBI Online BLAST HTML report (BlastOutput XML → tabs)."""

from __future__ import annotations

import html
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path

from repeatviewer.backends.ncbi_blast_online import NcbiBlastOnlineError, _xml_text

_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
_NCBI_TOOL = "RepeatViewer"
_NCBI_EMAIL = "repeatviewer-noreply@local"


# NCBI Graphic Summary alignment-score color key (bit score).
_SCORE_COLOR_LT40 = "#000000"
_SCORE_COLOR_40_50 = "#0000FF"
_SCORE_COLOR_50_80 = "#00FF00"
_SCORE_COLOR_80_200 = "#FF00FF"
_SCORE_COLOR_GE200 = "#FF0000"

@dataclass
class BlastHsp:
    bit_score: float
    score: float
    evalue: str
    identity: int
    positive: int
    gaps: int
    align_len: int
    qstart: int
    qend: int
    sstart: int
    send: int
    sstrand: str
    pident: float
    qseq: str = ""
    hseq: str = ""
    midline: str = ""


@dataclass
class BlastHit:
    accession: str
    hit_id: str
    title: str
    length: int
    hsps: list[BlastHsp] = field(default_factory=list)

    @property
    def max_score(self) -> float:
        return max((h.bit_score for h in self.hsps), default=0.0)

    @property
    def total_score(self) -> float:
        return sum(h.bit_score for h in self.hsps)

    @property
    def best_evalue(self) -> str:
        if not self.hsps:
            return "N/A"
        return min(self.hsps, key=lambda h: _evalue_sort_key(h.evalue)).evalue

    @property
    def best_pident(self) -> float:
        return max((h.pident for h in self.hsps), default=0.0)


@dataclass
class BlastIteration:
    query_id: str
    query_def: str
    query_len: int
    hits: list[BlastHit] = field(default_factory=list)


@dataclass
class BlastOutputResult:
    program: str
    database: str
    version: str
    rid: str
    query_name: str
    iterations: list[BlastIteration] = field(default_factory=list)

    @property
    def primary(self) -> BlastIteration | None:
        return self.iterations[0] if self.iterations else None


def _evalue_sort_key(text: str) -> float:
    try:
        return float(text)
    except (TypeError, ValueError):
        return 1e9


def _short_subject_id(hit: ET.Element) -> str:
    accession = _xml_text(hit, "Hit_accession")
    if accession:
        return accession.split()[0]
    hit_id = _xml_text(hit, "Hit_id")
    if hit_id:
        parts = [p for p in hit_id.split("|") if p]
        if len(parts) >= 4 and parts[-2].lower() in {
            "ref",
            "gb",
            "emb",
            "dbj",
            "pdb",
            "sp",
            "tpg",
            "tpe",
            "tpd",
            "gnl",
        }:
            return parts[-1] or parts[-2]
        if len(parts) >= 2:
            return parts[-1]
        return hit_id.split()[0]
    return "subject"


def _hit_title(hit: ET.Element) -> str:
    title = _xml_text(hit, "Hit_def")
    if title:
        return " ".join(title.replace("\t", " ").split())
    return _short_subject_id(hit)


def parse_blast_output_xml(
    xml_text: str,
    *,
    rid: str = "",
    query_name: str = "",
    program: str = "",
    database: str = "",
) -> BlastOutputResult:
    """Parse classic NCBI BlastOutput XML into structured hits/HSPs."""
    text = (xml_text or "").strip()
    if not text:
        raise NcbiBlastOnlineError("Empty NCBI BLAST XML.")
    if "<BlastOutput" in text and not text.lstrip().startswith("<?xml"):
        start = text.find("<?xml")
        if start < 0:
            start = text.find("<BlastOutput")
        text = text[start:]
    try:
        root = ET.fromstring(text)
    except ET.ParseError as exc:
        raise NcbiBlastOnlineError(
            f"Could not parse NCBI BLAST XML results: {exc}\n"
            f"Response starts with:\n{xml_text[:400]}"
        ) from exc

    result = BlastOutputResult(
        program=program or _xml_text(root, "BlastOutput_program") or "blast",
        database=database or _xml_text(root, "BlastOutput_db") or "",
        version=_xml_text(root, "BlastOutput_version"),
        rid=rid,
        query_name=query_name,
    )
    for iteration in root.findall(".//Iteration"):
        qdef = _xml_text(iteration, "Iteration_query-def") or _xml_text(
            iteration, "Iteration_query-ID", "query"
        )
        qid = qdef.split()[0] if qdef else "query"
        try:
            qlen = int(float(_xml_text(iteration, "Iteration_query-len", "0") or 0))
        except (TypeError, ValueError):
            qlen = 0
        it = BlastIteration(query_id=qid, query_def=qdef, query_len=qlen)
        for hit_el in iteration.findall("Iteration_hits/Hit"):
            try:
                hit_len = int(float(_xml_text(hit_el, "Hit_len", "0") or 0))
            except (TypeError, ValueError):
                hit_len = 0
            hit = BlastHit(
                accession=_short_subject_id(hit_el),
                hit_id=_xml_text(hit_el, "Hit_id"),
                title=_hit_title(hit_el),
                length=hit_len,
            )
            for hsp_el in hit_el.findall("Hit_hsps/Hsp"):
                identity = int(float(_xml_text(hsp_el, "Hsp_identity", "0") or 0))
                positive = int(
                    float(
                        _xml_text(hsp_el, "Hsp_positive", str(identity)) or identity
                    )
                )
                align_len = int(float(_xml_text(hsp_el, "Hsp_align-len", "0") or 0))
                gaps = int(float(_xml_text(hsp_el, "Hsp_gaps", "0") or 0))
                qstart = int(float(_xml_text(hsp_el, "Hsp_query-from", "0") or 0))
                qend = int(float(_xml_text(hsp_el, "Hsp_query-to", "0") or 0))
                sstart = int(float(_xml_text(hsp_el, "Hsp_hit-from", "0") or 0))
                send = int(float(_xml_text(hsp_el, "Hsp_hit-to", "0") or 0))
                if sstart <= send:
                    sstrand = "plus"
                else:
                    sstrand = "minus"
                    sstart, send = send, sstart
                if qstart > qend:
                    qstart, qend = qend, qstart
                pident = (100.0 * identity / align_len) if align_len else 0.0
                hit.hsps.append(
                    BlastHsp(
                        bit_score=float(
                            _xml_text(hsp_el, "Hsp_bit-score", "0") or 0
                        ),
                        score=float(_xml_text(hsp_el, "Hsp_score", "0") or 0),
                        evalue=_xml_text(hsp_el, "Hsp_evalue", "0") or "0",
                        identity=identity,
                        positive=positive,
                        gaps=gaps,
                        align_len=align_len,
                        qstart=qstart,
                        qend=qend,
                        sstart=sstart,
                        send=send,
                        sstrand=sstrand,
                        pident=pident,
                        qseq=_xml_text(hsp_el, "Hsp_qseq"),
                        hseq=_xml_text(hsp_el, "Hsp_hseq"),
                        midline=_xml_text(hsp_el, "Hsp_midline"),
                    )
                )
            if hit.hsps:
                hit.hsps.sort(key=lambda h: h.bit_score, reverse=True)
                it.hits.append(hit)
        it.hits.sort(key=lambda h: h.max_score, reverse=True)
        result.iterations.append(it)
    return result


def _union_coverage(intervals: list[tuple[int, int]], query_len: int) -> float:
    if query_len <= 0 or not intervals:
        return 0.0
    cleaned = sorted((min(a, b), max(a, b)) for a, b in intervals if a > 0 and b > 0)
    if not cleaned:
        return 0.0
    merged: list[list[int]] = []
    for start, end in cleaned:
        if not merged or start > merged[-1][1] + 1:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    covered = sum(end - start + 1 for start, end in merged)
    return 100.0 * min(covered, query_len) / query_len


def _guess_organism(title: str) -> str:
    """Best-effort organism label from a Hit_def / title string.

    Returns a binomial (Genus species), optionally with one subsp./var. epithet.
    Does **not** consume following gene/product words (those belong in Description).
    """
    text = " ".join(str(title or "").split())
    if not text:
        return "Unknown"
    cleaned = re.sub(
        r"^(PREDICTED:|PROVISIONAL:|RecName:|SubName:)\s*",
        "",
        text,
        flags=re.IGNORECASE,
    )
    # Common NCBI pattern: "... [Genus species ...]"
    bracket = re.search(
        r"\[([A-Z][a-z]+(?:\s+(?:subsp\.|ssp\.|var\.|f\.|[a-z]+)){1,3})\]",
        cleaned,
    )
    if bracket:
        return _normalize_binomial(bracket.group(1))
    # Leading Genus + species, then optional subsp./var. epithet only.
    match = re.match(
        r"^([A-Z][a-z]+)\s+([a-z]+|sp\.|spp\.)"
        r"(?:\s+(subsp\.|ssp\.|var\.|f\.)\s+([a-z]+))?",
        cleaned,
    )
    if match:
        genus, species, rank, epi = match.groups()
        if rank and epi:
            return f"{genus} {species} {rank} {epi}"
        return f"{genus} {species}"
    return "Unknown"


def _normalize_binomial(text: str) -> str:
    """Keep Genus species (+ optional subsp/var) from a bracket organism string."""
    match = re.match(
        r"^([A-Z][a-z]+)\s+([a-z]+|sp\.|spp\.)"
        r"(?:\s+(subsp\.|ssp\.|var\.|f\.)\s+([a-z]+))?",
        " ".join(str(text or "").split()),
    )
    if not match:
        parts = " ".join(str(text or "").split()).split()
        return " ".join(parts[:2]) if len(parts) >= 2 else (parts[0] if parts else "Unknown")
    genus, species, rank, epi = match.groups()
    if rank and epi:
        return f"{genus} {species} {rank} {epi}"
    return f"{genus} {species}"


def _description_cluster_key(title: str) -> str:
    """Sort key so similar titles (e.g. VulMITE variants) sit together."""
    text = " ".join(str(title or "").split())
    text = re.sub(
        r"^(PREDICTED:|PROVISIONAL:|RecName:|SubName:)\s*",
        "",
        text,
        flags=re.IGNORECASE,
    )
    return text.casefold()


def _ncbi_accession_url(accession: str) -> str:
    acc = (accession or "").strip()
    if not acc:
        return "#"
    return f"https://www.ncbi.nlm.nih.gov/nuccore/{html.escape(acc, quote=True)}"


def _fmt_float(value: float, digits: int = 1) -> str:
    return f"{value:.{digits}f}"


def _fmt_pct(value: float) -> str:
    return f"{value:.2f}%"


def _fmt_score(value: float) -> str:
    if value <= 0:
        return ""
    if abs(value - round(value)) < 0.05:
        return str(int(round(value)))
    return _fmt_float(value)


def _organism_lineage_parts(organism: str) -> list[str]:
    """Build a shallow lineage path from a binomial / subspecies string."""
    org = " ".join(str(organism or "").split())
    if not org or org == "Unknown":
        return ["Unknown"]
    tokens = org.split()
    if len(tokens) == 1:
        return [tokens[0]]
    parts: list[str] = [tokens[0]]
    if tokens[1][:1].islower() or tokens[1] in {"sp.", "spp."}:
        parts.append(f"{tokens[0]} {tokens[1]}")
        i = 2
        while i < len(tokens):
            if tokens[i].lower() in {"subsp.", "ssp.", "var.", "f."} and i + 1 < len(
                tokens
            ):
                parts.append(" ".join(tokens[: i + 2]))
                i += 2
                continue
            break
        return parts
    return [org]


# Common NCBI “Blast Name” labels for high-level taxa (taxdump name class).
_KNOWN_BLAST_NAMES: dict[str, str] = {
    "Eukaryota": "eukaryotes",
    "Bacteria": "bacteria",
    "Archaea": "archaea",
    "Viruses": "viruses",
    "Fungi": "fungi",
    "Metazoa": "metazoans",
    "Viridiplantae": "green plants",
    "Streptophyta": "streptophytes",
    "Embryophyta": "land plants",
    "Tracheophyta": "vascular plants",
    "Spermatophyta": "seed plants",
    "Magnoliopsida": "flowering plants",
    "Mesangiospermae": "flowering plants",
    "eudicotyledons": "eudicots",
    "Pentapetalae": "eudicots",
    "Liliopsida": "monocots",
    "monocotyledons": "monocots",
    "Amniota": "amniotes",
    "Mammalia": "mammals",
    "Primates": "primates",
    "Rodentia": "rodents",
    "Aves": "birds",
    "Insecta": "insects",
    "Arthropoda": "arthropods",
    "Nematoda": "nematodes",
    "Proteobacteria": "proteobacteria",
    "Firmicutes": "firmicutes",
    "Actinobacteria": "high GC Gram+",
}


def _blast_name_along_path(path: list[str], index: int) -> str:
    """Nearest known Blast Name at or above this lineage node."""
    for name in reversed(path[: index + 1]):
        mapped = _KNOWN_BLAST_NAMES.get(name)
        if mapped:
            return mapped
    return ""


@dataclass
class _TaxonRecord:
    tax_id: str
    scientific_name: str
    lineage_names: list[str]  # ancestors + self, excluding cellular organisms


def _eutils_get(url: str, *, timeout: float = 45.0) -> str:
    req = urllib.request.Request(
        url,
        headers={"User-Agent": f"{_NCBI_TOOL}/1.0 ({_NCBI_EMAIL})"},
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read().decode("utf-8", errors="replace")


def _taxonomy_search_taxid(organism: str) -> str | None:
    term = urllib.parse.quote(organism)
    url = (
        f"{_EUTILS}/esearch.fcgi?db=taxonomy&term={term}"
        f"&retmode=json&retmax=1&tool={_NCBI_TOOL}&email={_NCBI_EMAIL}"
    )
    raw = _eutils_get(url)
    try:
        data = json.loads(raw)
        ids = data.get("esearchresult", {}).get("idlist") or []
        return str(ids[0]) if ids else None
    except (TypeError, ValueError, AttributeError):
        return None


def _taxonomy_parse_efetch(xml_text: str) -> dict[str, _TaxonRecord]:
    out: dict[str, _TaxonRecord] = {}
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        return out
    for taxon in root.findall("Taxon"):
        tax_id = (taxon.findtext("TaxId") or "").strip()
        sci = (taxon.findtext("ScientificName") or "").strip()
        if not tax_id or not sci:
            continue
        names: list[str] = []
        lineage_ex = taxon.find("LineageEx")
        if lineage_ex is not None:
            for anc in lineage_ex.findall("Taxon"):
                name = (anc.findtext("ScientificName") or "").strip()
                if not name or name.lower() in {"cellular organisms"}:
                    continue
                names.append(name)
        if not names or names[-1] != sci:
            names.append(sci)
        out[tax_id] = _TaxonRecord(
            tax_id=tax_id, scientific_name=sci, lineage_names=names
        )
    return out


def _lookup_organism_lineages(
    organisms: list[str],
    *,
    enable_network: bool = True,
) -> dict[str, list[str]]:
    """Map organism label → full NCBI lineage names (Eukaryota…).

    Falls back to empty mapping when offline / lookup fails so the report still
    builds from title-parsed shallow paths.
    """
    result: dict[str, list[str]] = {}
    if not enable_network:
        return result
    unique = []
    seen: set[str] = set()
    for org in organisms:
        key = " ".join(str(org or "").split())
        if not key or key == "Unknown" or key in seen:
            continue
        seen.add(key)
        unique.append(key)
    if not unique:
        return result

    tax_ids: list[str] = []
    org_to_tax: dict[str, str] = {}
    for i, org in enumerate(unique):
        try:
            tax_id = _taxonomy_search_taxid(org)
        except Exception:  # noqa: BLE001
            tax_id = None
        if tax_id:
            org_to_tax[org] = tax_id
            if tax_id not in tax_ids:
                tax_ids.append(tax_id)
        # Be polite to NCBI E-utilities (~3 requests / second).
        if i + 1 < len(unique):
            time.sleep(0.34)
    if not tax_ids:
        return result

    # Batch efetch (comma-separated ids).
    for start in range(0, len(tax_ids), 40):
        batch = tax_ids[start : start + 40]
        ids = ",".join(batch)
        url = (
            f"{_EUTILS}/efetch.fcgi?db=taxonomy&id={ids}"
            f"&retmode=xml&tool={_NCBI_TOOL}&email={_NCBI_EMAIL}"
        )
        try:
            xml_text = _eutils_get(url)
            parsed = _taxonomy_parse_efetch(xml_text)
        except Exception:  # noqa: BLE001
            parsed = {}
        for org, tax_id in org_to_tax.items():
            rec = parsed.get(tax_id)
            if rec is not None:
                result[org] = list(rec.lineage_names)
        if start + 40 < len(tax_ids):
            time.sleep(0.34)
    return result


@dataclass
class _LineageNode:
    name: str
    depth: int
    blast_name: str = ""
    hit_count: int = 0
    max_score: float = 0.0
    leaf_hits: list[BlastHit] = field(default_factory=list)
    is_leaf: bool = False


def _lineage_path_for_hit(
    hit: BlastHit,
    org_lineages: dict[str, list[str]],
) -> list[str]:
    org = _guess_organism(hit.title)
    full = org_lineages.get(org)
    if full:
        return list(full)
    shallow = _organism_lineage_parts(org)
    known_full = [
        list(path)
        for path in org_lineages.values()
        if _is_full_taxonomy_path(path)
    ]
    return _graft_lineage_path(shallow, known_full)


_DOMAIN_ROOTS = frozenset({"Eukaryota", "Bacteria", "Archaea", "Viruses"})


def _is_full_taxonomy_path(path: list[str]) -> bool:
    return bool(path) and path[0] in _DOMAIN_ROOTS


def _graft_lineage_path(path: list[str], known_full: list[list[str]]) -> list[str]:
    """Attach a shallow Genus/species path under a known NCBI taxonomy prefix.

    Without this, a hit that failed Taxonomy lookup creates a second root-level
    ``Beta`` next to the real ``Eukaryota → … → Beta`` branch.
    """
    if not path or _is_full_taxonomy_path(path):
        return path
    leaf = path[-1]
    for full in known_full:
        if leaf in full:
            return full[: full.index(leaf) + 1]
    for full in known_full:
        for i, name in enumerate(path):
            if name in full:
                idx = full.index(name)
                return full[:idx] + path[i:]
    return path


def _build_lineage_rows(
    iteration: BlastIteration,
    *,
    org_lineages: dict[str, list[str]] | None = None,
) -> list[_LineageNode]:
    """Aggregate hits into an NCBI-like indented lineage table."""
    org_lineages = org_lineages or {}
    # Fill gaps: organisms missing from taxonomy still inherit a shared prefix
    # (e.g. Eukaryota…Beta) from siblings that did resolve.
    known_full = [
        list(path)
        for path in org_lineages.values()
        if _is_full_taxonomy_path(path)
    ]
    if known_full:
        for org in {
            _guess_organism(hit.title) for hit in iteration.hits
        }:
            if org in org_lineages and _is_full_taxonomy_path(org_lineages[org]):
                continue
            shallow = org_lineages.get(org) or _organism_lineage_parts(org)
            org_lineages[org] = _graft_lineage_path(list(shallow), known_full)

    nodes: dict[str, _LineageNode] = {}
    order: list[str] = []

    def ensure(path: tuple[str, ...]) -> _LineageNode:
        key = "\0".join(path)
        if key not in nodes:
            blast = _blast_name_along_path(list(path), len(path) - 1)
            nodes[key] = _LineageNode(
                name=path[-1], depth=len(path) - 1, blast_name=blast
            )
            order.append(key)
        return nodes[key]

    for hit in iteration.hits:
        path = tuple(_lineage_path_for_hit(hit, org_lineages))
        if not path:
            continue
        for i in range(len(path)):
            node = ensure(path[: i + 1])
            node.hit_count += 1
            node.max_score = max(node.max_score, hit.max_score)
        leaf = ensure(path)
        leaf.is_leaf = True
        leaf.leaf_hits.append(hit)

    # Depth-first by path so one Beta clade stays contiguous.
    ordered_keys = sorted(nodes.keys(), key=lambda k: (k.split("\0"),))
    return [nodes[key] for key in ordered_keys]


def _pack_hsps_into_lanes(hsps: list[BlastHsp]) -> list[list[BlastHsp]]:
    """Pack HSPs so overlapping intervals go on separate lanes (no paint overlap)."""
    ordered = sorted(hsps, key=lambda h: (int(h.qstart), int(h.qend), -float(h.bit_score)))
    lanes: list[list[BlastHsp]] = []
    lane_ends: list[int] = []
    for hsp in ordered:
        start = int(hsp.qstart)
        end = int(hsp.qend)
        placed = False
        for i, last_end in enumerate(lane_ends):
            if start > last_end:
                lanes[i].append(hsp)
                lane_ends[i] = end
                placed = True
                break
        if not placed:
            lanes.append([hsp])
            lane_ends.append(end)
    return lanes


def _alignment_score_color(bit_score: float) -> str:
    """NCBI Graphic Summary colors keyed by bit score."""
    try:
        score = float(bit_score)
    except (TypeError, ValueError):
        score = 0.0
    if score >= 200:
        return _SCORE_COLOR_GE200
    if score >= 80:
        return _SCORE_COLOR_80_200
    if score >= 50:
        return _SCORE_COLOR_50_80
    if score >= 40:
        return _SCORE_COLOR_40_50
    return _SCORE_COLOR_LT40


def _estimate_desc_height(text: str, *, col_width_px: int, font_px: int = 11) -> int:
    """Rough wrapped-text height so the SVG leaves room for the full description."""
    words = " ".join(str(text or "").split()).split(" ")
    if not words or words == [""]:
        return font_px + 6
    # ~0.55em average glyph width for sans-serif at this size.
    chars_per_line = max(12, int(col_width_px / max(font_px * 0.55, 1)))
    lines = 1
    current = 0
    for word in words:
        wlen = len(word)
        if current == 0:
            current = wlen
            continue
        if current + 1 + wlen <= chars_per_line:
            current += 1 + wlen
        else:
            lines += 1
            current = wlen
            # Very long tokens still consume multiple visual lines.
            while current > chars_per_line:
                lines += 1
                current -= chars_per_line
    line_h = int(round(font_px * 1.3))
    return lines * line_h + 6


def _alignment_score_legend_html() -> str:
    items = [
        (_SCORE_COLOR_LT40, "&lt; 40"),
        (_SCORE_COLOR_40_50, "40 — 50"),
        (_SCORE_COLOR_50_80, "50 — 80"),
        (_SCORE_COLOR_80_200, "80 — 200"),
        (_SCORE_COLOR_GE200, "≥ 200"),
    ]
    swatches = "".join(
        f'<span class="score-swatch"><i style="background:{color}"></i>{label}</span>'
        for color, label in items
    )
    return (
        '<div class="score-legend" role="img" aria-label="Alignment Scores">'
        "<strong>Alignment Scores</strong>"
        f"{swatches}</div>"
    )


def _graphic_svg(iteration: BlastIteration, *, max_subjects: int = 80) -> str:
    query_len = max(1, int(iteration.query_len) or 1)
    hits = iteration.hits[:max_subjects]
    if iteration.query_len <= 0:
        ends = [hsp.qend for hit in iteration.hits for hsp in hit.hsps]
        query_len = max(ends) if ends else 1
        hits = iteration.hits[:max_subjects]

    left = 480  # wide description column; show full titles
    top = 28
    width = 640
    lane_h = 8
    lane_gap = 1
    subject_pad = 5
    query_h = 10
    header_h = query_h + 20
    sep_gap = 4

    packed: list[tuple[BlastHit, list[list[BlastHsp]], int]] = []
    for hit in hits:
        lanes = _pack_hsps_into_lanes(hit.hsps)
        lane_block = max(1, len(lanes)) * (lane_h + lane_gap) - lane_gap
        desc = hit.title or hit.accession
        text_h = _estimate_desc_height(desc, col_width_px=left - 16)
        block_h = max(lane_block, text_h, lane_h)
        packed.append((hit, lanes, block_h))

    height = (
        top
        + header_h
        + sum(block_h + subject_pad + sep_gap for *_, block_h in packed)
        + 20
    )

    def x_of(bp: int) -> float:
        return left + (max(0, bp - 1) / query_len) * width

    parts: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {left + width + 24} {height}" '
        f'width="100%" style="max-width:1400px;background:#fff">'
    ]
    parts.append(
        f'<text x="8" y="{top + 10}" font-size="11" font-family="sans-serif" '
        f'font-weight="600" fill="#333">Description</text>'
    )
    parts.append(
        f'<rect x="{left}" y="{top}" width="{width}" height="{query_h}" '
        f'rx="2" fill="#2a9d8f"/>'
    )
    for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
        bp = int(round(1 + frac * (query_len - 1)))
        x = x_of(bp)
        parts.append(
            f'<line x1="{x:.1f}" y1="{top + query_h}" x2="{x:.1f}" '
            f'y2="{top + query_h + 6}" stroke="#666" stroke-width="1"/>'
        )
        parts.append(
            f'<text x="{x:.1f}" y="{top + query_h + 18}" font-size="10" '
            f'text-anchor="middle" font-family="sans-serif" fill="#444">{bp}</text>'
        )

    y_cursor = top + header_h
    for idx, (hit, lanes, block_h) in enumerate(packed):
        desc = hit.title or hit.accession
        parts.append(
            f'<foreignObject x="4" y="{y_cursor}" width="{left - 12}" height="{block_h}">'
            f'<div xmlns="http://www.w3.org/1999/xhtml" '
            f'style="font:11px/1.3 -apple-system,BlinkMacSystemFont,\'Segoe UI\',sans-serif;'
            f'color:#222;white-space:normal;overflow-wrap:anywhere;word-break:break-word;'
            f'padding:0 6px 0 2px;">'
            f"{html.escape(desc)}</div></foreignObject>"
        )
        for lane_i, lane_hsps in enumerate(lanes):
            y = y_cursor + lane_i * (lane_h + lane_gap)
            sorted_lane = sorted(lane_hsps, key=lambda h: h.qstart)
            if len(sorted_lane) >= 2:
                x1 = x_of(sorted_lane[0].qstart)
                x2 = x_of(sorted_lane[-1].qend)
                mid_y = y + lane_h / 2
                parts.append(
                    f'<line x1="{x1:.1f}" y1="{mid_y:.1f}" x2="{x2:.1f}" '
                    f'y2="{mid_y:.1f}" stroke="#777" stroke-width="1"/>'
                )
            for hsp in sorted_lane:
                x1 = x_of(hsp.qstart)
                x2 = x_of(hsp.qend)
                w = max(2.0, x2 - x1)
                color = _alignment_score_color(hsp.bit_score)
                tip = (
                    f"{hit.accession}: {hsp.qstart}-{hsp.qend} · "
                    f"score={hsp.bit_score:.1f} bits · "
                    f"{hsp.pident:.1f}% id · E={hsp.evalue}"
                )
                parts.append(
                    f'<rect x="{x1:.1f}" y="{y + 1}" width="{w:.1f}" height="{lane_h - 2}" '
                    f'rx="1" fill="{color}">'
                    f"<title>{html.escape(tip)}</title></rect>"
                )
        y_cursor += block_h + subject_pad
        # Thin grey separator between subjects (not after the last).
        if idx < len(packed) - 1:
            y_sep = y_cursor - subject_pad / 2
            parts.append(
                f'<line x1="4" y1="{y_sep:.1f}" x2="{left + width:.1f}" y2="{y_sep:.1f}" '
                f'stroke="#c8c8c8" stroke-width="1"/>'
            )
            y_cursor += sep_gap

    if len(iteration.hits) > max_subjects:
        parts.append(
            f'<text x="{left}" y="{height - 8}" font-size="11" font-family="sans-serif" '
            f'fill="#666">Showing top {max_subjects} of {len(iteration.hits)} subjects</text>'
        )
    parts.append("</svg>")
    return "\n".join(parts)


def _descriptions_table(
    iteration: BlastIteration,
    *,
    org_lineages: dict[str, list[str]] | None = None,
) -> str:
    query_len = max(1, int(iteration.query_len) or 1)
    if iteration.query_len <= 0:
        ends = [hsp.qend for hit in iteration.hits for hsp in hit.hsps]
        if ends:
            query_len = max(ends)
    org_lineages = org_lineages or {}
    # One group per scientific name (not score-order fragments), ordered by
    # best hit score so top organisms stay near the top.
    by_org: dict[str, list[BlastHit]] = {}
    org_best: dict[str, float] = {}
    for hit in iteration.hits:
        org = _guess_organism(hit.title)
        by_org.setdefault(org, []).append(hit)
        org_best[org] = max(org_best.get(org, 0.0), hit.max_score)
    groups = sorted(by_org.items(), key=lambda kv: -org_best.get(kv[0], 0.0))

    blocks: list[str] = []
    for org, hits in groups:
        # Cluster similar descriptions (VulMITE… together), then by score.
        hits_sorted = sorted(
            hits,
            key=lambda h: (_description_cluster_key(h.title), -h.max_score),
        )
        rows: list[str] = []
        for hit in hits_sorted:
            cover = _union_coverage([(h.qstart, h.qend) for h in hit.hsps], query_len)
            lineage = _lineage_path_for_hit(hit, org_lineages)
            lineage_attr = "|".join(lineage)
            rows.append(
                "<tr class='hit-row' "
                f'data-organism="{html.escape(org, quote=True)}" '
                f'data-lineage="{html.escape(lineage_attr, quote=True)}" '
                f'data-accession="{html.escape(hit.accession, quote=True)}" '
                f'data-max="{hit.max_score}" '
                f'data-total="{hit.total_score}" '
                f'data-cover="{cover}" '
                f'data-evalue="{html.escape(hit.best_evalue, quote=True)}" '
                f'data-pident="{hit.best_pident}" '
                f'data-acclen="{hit.length}" '
                f'data-desc="{html.escape(_description_cluster_key(hit.title), quote=True)}">'
                f'<td class="desc" data-col="desc">'
                f"{html.escape(hit.title or hit.accession)}</td>"
                f'<td data-col="organism">{html.escape(org)}</td>'
                f'<td class="num" data-col="max">{_fmt_float(hit.max_score)}</td>'
                f'<td class="num" data-col="total">{_fmt_float(hit.total_score)}</td>'
                f'<td class="num" data-col="cover">{_fmt_pct(cover)}</td>'
                f'<td class="num" data-col="evalue">{html.escape(hit.best_evalue)}</td>'
                f'<td class="num" data-col="pident">{_fmt_pct(hit.best_pident)}</td>'
                f'<td class="num" data-col="acclen">{hit.length or ""}</td>'
                f'<td data-col="acc"><a href="{_ncbi_accession_url(hit.accession)}" '
                f'target="_blank" rel="noopener noreferrer">'
                f"{html.escape(hit.accession)}</a></td>"
                "</tr>"
            )
        blocks.append(
            f'<tbody class="org-group" data-organism="{html.escape(org, quote=True)}">'
            f'<tr class="org-header"><td colspan="9">'
            f'<span class="org-label">Organism</span> '
            f"{html.escape(org)} "
            f'<span class="org-count">({len(hits_sorted)} hit'
            f'{"s" if len(hits_sorted) != 1 else ""})</span>'
            f"</td></tr>"
            f'{"".join(rows)}</tbody>'
        )
    if not blocks:
        return "<p class='empty'>No significant hits.</p>"
    toolbar = (
        '<div class="desc-toolbar">'
        '<label class="desc-option">'
        '<input type="checkbox" id="group-by-organism" checked> '
        "Group by organism"
        "</label>"
        '<button type="button" class="desc-clear-filter" id="clear-org-filter" hidden>'
        "Clear lineage filter"
        "</button>"
        '<span class="sort-hint">Click a column header to sort. '
        "Click again to reverse.</span>"
        "</div>"
    )
    return (
        toolbar
        + '<table class="desc-table" id="desc-table" data-group-by-organism="1">'
        "<thead><tr>"
        '<th class="sortable" data-sort="desc">Description</th>'
        '<th class="sortable" data-sort="organism">Scientific Name</th>'
        '<th class="sortable" data-sort="max">Max Score</th>'
        '<th class="sortable" data-sort="total">Total Score</th>'
        '<th class="sortable" data-sort="cover">Query Cover</th>'
        '<th class="sortable" data-sort="evalue">E value</th>'
        '<th class="sortable" data-sort="pident">Per. Ident</th>'
        '<th class="sortable" data-sort="acclen">Acc. Len</th>'
        '<th class="sortable" data-sort="acc">Accession</th>'
        "</tr></thead>"
        f'{"".join(blocks)}</table>'
    )


def _lineage_html(
    iteration: BlastIteration,
    *,
    org_lineages: dict[str, list[str]] | None = None,
    taxonomy_ok: bool = False,
) -> str:
    rows = _build_lineage_rows(iteration, org_lineages=org_lineages)
    if not rows:
        return "<p class='empty'>No significant hits.</p>"
    n_sel = len(iteration.hits)
    body: list[str] = []
    for node in rows:
        indent = ". " * node.depth
        org_label = f"{indent}{node.name}" if node.depth else node.name
        blast_name = node.blast_name
        score_cell = _fmt_score(node.max_score) if node.is_leaf else ""
        desc = ""
        if node.is_leaf and node.leaf_hits:
            desc = (
                f'<a href="#" class="org-filter" data-organism="'
                f'{html.escape(node.name, quote=True)}">'
                f"{html.escape(node.name)} hits</a>"
            )
        body.append(
            "<tr>"
            f'<td class="org"><a href="#" class="org-filter" data-organism="'
            f'{html.escape(node.name, quote=True)}">'
            f"{html.escape(org_label)}</a></td>"
            f"<td>{html.escape(blast_name)}</td>"
            f'<td class="num">{score_cell}</td>'
            f'<td class="num"><a href="#" class="org-filter" data-organism="'
            f'{html.escape(node.name, quote=True)}">{node.hit_count}</a></td>'
            f"<td>{desc}</td>"
            "</tr>"
        )
    if taxonomy_ok:
        note = (
            "Lineage ranks come from NCBI Taxonomy (Eukaryota → …) looked up for "
            "each scientific name in the BLAST hit titles."
        )
    else:
        note = (
            "NCBI Taxonomy lookup was unavailable, so lineage is inferred only from "
            "scientific names in hit titles (genus / species / subspecies). "
            "That is why high ranks such as <b>Eukaryota</b> may be missing."
        )
    return (
        f'<p class="sel-count">{n_sel} sequence{"s" if n_sel != 1 else ""} selected.</p>'
        f'<p class="note">{note}</p>'
        '<table class="lineage-table">'
        "<thead><tr>"
        "<th>Organism</th><th>Blast Name</th><th>Score</th>"
        "<th>Number of Hits</th><th>Description</th>"
        "</tr></thead>"
        f"<tbody>{''.join(body)}</tbody></table>"
    )


_CSS = """
:root { color-scheme: light; }
body { margin: 0; font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
  font-size: 13px; color: #222; background: #f7f7f7; }
header.page { background: #fff; border-bottom: 1px solid #d0d0d0; padding: 12px 18px; }
header.page h1 { margin: 0 0 4px; font-size: 18px; }
header.page .meta { color: #555; margin: 0; }
.tabs { display: flex; gap: 0; background: #eef2f5; border-bottom: 1px solid #c9d2da;
  padding: 0 12px; }
.tabs button { border: 0; background: transparent; padding: 10px 14px; cursor: pointer;
  font: inherit; color: #234; border-bottom: 2px solid transparent; }
.tabs button.active { border-bottom-color: #1a73a8; color: #0b4f75; font-weight: 600;
  background: #fff; }
.panel { display: none; padding: 14px 18px 28px; background: #fff; min-height: 60vh; }
.panel.active { display: block; }
.desc-table, .lineage-table { width: 100%; border-collapse: collapse; background: #fff; }
.desc-table th, .desc-table td, .lineage-table th, .lineage-table td {
  border-bottom: 1px solid #e3e3e3; padding: 6px 8px; text-align: left; vertical-align: top; }
.desc-table th, .lineage-table th { background: #e8f0f7; }
.desc-table tr:nth-child(even), .lineage-table tr:nth-child(even) { background: #fafbfc; }
.desc-table tr:hover, .lineage-table tr:hover { background: #f5faff; }
.desc-table td.num, .lineage-table td.num { text-align: right; white-space: nowrap;
  font-variant-numeric: tabular-nums; }
.desc-table td.desc {
  max-width: 480px; white-space: normal; overflow-wrap: anywhere; word-break: break-word;
  line-height: 1.35;
}
.desc-table tbody.org-group + tbody.org-group {
  border-top: 2px solid #b0b8c0;
}
.desc-table tr.org-header td {
  background: #eef2f5; color: #333; font-weight: 600;
  border-top: 2px solid #b0b8c0; border-bottom: 1px solid #d0d5da;
  padding-top: 8px; padding-bottom: 8px;
}
.desc-table .org-label {
  display: inline-block; font-size: 11px; text-transform: uppercase;
  letter-spacing: 0.04em; color: #567; margin-right: 6px;
}
.desc-table .org-count { color: #666; font-weight: 500; }
.desc-table th.sortable {
  cursor: pointer; user-select: none; white-space: nowrap;
}
.desc-table th.sortable:hover { background: #dce6ef; }
.desc-table th.sortable.sorted-asc::after { content: " ▲"; font-size: 10px; }
.desc-table th.sortable.sorted-desc::after { content: " ▼"; font-size: 10px; }
.desc-toolbar {
  display: flex; flex-wrap: wrap; align-items: center; gap: 12px 18px;
  margin: 0 0 10px;
}
.desc-option {
  display: inline-flex; align-items: center; gap: 6px;
  font-size: 13px; color: #234; cursor: pointer; user-select: none;
}
.desc-option input { margin: 0; }
.desc-clear-filter {
  border: 1px solid #c9d2da; background: #fff; color: #0b4f75;
  border-radius: 3px; padding: 3px 8px; font: inherit; font-size: 12px;
  cursor: pointer;
}
.desc-clear-filter:hover { background: #f0f6fa; }
.sort-hint { color: #666; font-size: 12px; margin: 0; }
.desc-table[data-group-by-organism="0"] tr.org-header { display: none; }
.desc-table[data-group-by-organism="0"] tbody.org-group + tbody.org-group {
  border-top: 0;
}
.lineage-table td.org { font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
  white-space: pre; }
.sel-count { font-weight: 600; margin: 0 0 8px; }
a { color: #0b57d0; text-decoration: none; }
a:hover { text-decoration: underline; }
.empty { color: #666; font-style: italic; }
.note { color: #555; background: #f8fafc; border: 1px solid #e2e8f0; padding: 8px 10px; }
tr.dimmed { opacity: 0.25; }
.score-legend {
  display: flex; flex-wrap: wrap; align-items: center; gap: 14px 18px;
  background: #e8f5ef; border: 1px solid #c5ddd2; border-radius: 3px;
  padding: 8px 12px; margin: 0 0 12px; font-size: 12px;
}
.score-legend strong { margin-right: 4px; }
.score-swatch { display: inline-flex; align-items: center; gap: 6px; white-space: nowrap; }
.score-swatch i {
  display: inline-block; width: 12px; height: 12px; border: 1px solid #666;
  border-radius: 1px;
}
.graphic-wrap { overflow-x: auto; }
"""

_JS = """
(function () {
  var buttons = document.querySelectorAll(".tabs > button");
  var panels = document.querySelectorAll(".panel");
  function show(id) {
    buttons.forEach(function (b) {
      b.classList.toggle("active", b.getAttribute("data-tab") === id);
    });
    panels.forEach(function (p) {
      p.classList.toggle("active", p.id === id);
    });
  }
  buttons.forEach(function (b) {
    b.addEventListener("click", function () {
      show(b.getAttribute("data-tab"));
    });
  });
  var activeLineageFilter = null;
  var clearFilterBtn = document.getElementById("clear-org-filter");
  function applyLineageFilter(rank) {
    activeLineageFilter = rank || null;
    var rows = document.querySelectorAll("#desc-table tbody tr.hit-row");
    var kept = 0;
    rows.forEach(function (row) {
      var rowOrg = row.getAttribute("data-organism") || "";
      var lineage = (row.getAttribute("data-lineage") || "").split("|");
      var keep = !rank;
      if (rank) {
        keep = lineage.indexOf(rank) >= 0 ||
          rowOrg === rank ||
          rowOrg.indexOf(rank + " ") === 0 ||
          (rank.indexOf(rowOrg) === 0 && rank.length > rowOrg.length);
      }
      if (keep) kept += 1;
      row.classList.toggle("dimmed", !keep);
    });
    if (clearFilterBtn) {
      clearFilterBtn.hidden = !rank;
      if (rank) {
        clearFilterBtn.textContent = kept
          ? ("Clear filter (" + kept + " shown)")
          : "Clear filter (no matches)";
      }
    }
  }
  if (clearFilterBtn) {
    clearFilterBtn.addEventListener("click", function () {
      applyLineageFilter(null);
    });
  }
  document.querySelectorAll(".org-filter").forEach(function (a) {
    a.addEventListener("click", function (ev) {
      ev.preventDefault();
      var org = a.getAttribute("data-organism") || "";
      show("tab-descriptions");
      // Click same taxon again to clear; otherwise filter by lineage rank.
      if (activeLineageFilter === org) applyLineageFilter(null);
      else applyLineageFilter(org);
    });
  });

  var table = document.getElementById("desc-table");
  if (table) {
    var sortState = { key: null, dir: 1 };
    var numericKeys = { max: 1, total: 1, cover: 1, pident: 1, acclen: 1 };
    var groupToggle = document.getElementById("group-by-organism");
    function groupByOrganism() {
      return !groupToggle || groupToggle.checked;
    }
    function evalueKey(text) {
      var v = parseFloat(text);
      return isNaN(v) ? Number.POSITIVE_INFINITY : v;
    }
    function rowValue(row, key) {
      if (key === "evalue") return evalueKey(row.getAttribute("data-evalue") || "");
      if (key === "desc") return (row.getAttribute("data-desc") || "").toLowerCase();
      if (key === "organism") {
        return (row.getAttribute("data-organism") || "").toLowerCase();
      }
      if (key === "acc") {
        return (row.getAttribute("data-accession") || "").toLowerCase();
      }
      if (numericKeys[key]) {
        var n = parseFloat(row.getAttribute("data-" + key) || "0");
        return isNaN(n) ? 0 : n;
      }
      return "";
    }
    function sortRows(rows) {
      if (!sortState.key) return rows;
      var key = sortState.key;
      var dir = sortState.dir;
      rows.sort(function (a, b) {
        var va = rowValue(a, key);
        var vb = rowValue(b, key);
        if (va < vb) return -1 * dir;
        if (va > vb) return 1 * dir;
        return 0;
      });
      return rows;
    }
    function rebuildFlat(rows) {
      while (table.tBodies.length) table.removeChild(table.tBodies[0]);
      var body = document.createElement("tbody");
      body.className = "hit-list";
      rows.forEach(function (row) { body.appendChild(row); });
      table.appendChild(body);
    }
    function rebuildGrouped(rows) {
      var byOrg = {};
      var order = [];
      rows.forEach(function (row) {
        var org = row.getAttribute("data-organism") || "Unknown";
        if (!byOrg[org]) {
          byOrg[org] = [];
          order.push(org);
        }
        byOrg[org].push(row);
      });
      while (table.tBodies.length) table.removeChild(table.tBodies[0]);
      order.forEach(function (org) {
        var body = document.createElement("tbody");
        body.className = "org-group";
        body.setAttribute("data-organism", org);
        var header = document.createElement("tr");
        header.className = "org-header";
        var td = document.createElement("td");
        td.colSpan = 9;
        var safe = org.replace(/&/g, "&amp;").replace(/</g, "&lt;");
        td.innerHTML = '<span class="org-label">Organism</span> ' + safe +
          ' <span class="org-count">(' + byOrg[org].length +
          (byOrg[org].length === 1 ? " hit" : " hits") + ")</span>";
        header.appendChild(td);
        body.appendChild(header);
        byOrg[org].forEach(function (row) { body.appendChild(row); });
        table.appendChild(body);
      });
    }
    function rebuildTable(rows) {
      table.setAttribute(
        "data-group-by-organism",
        groupByOrganism() ? "1" : "0"
      );
      if (groupByOrganism()) rebuildGrouped(rows);
      else rebuildFlat(rows);
    }
    function currentRows() {
      return Array.prototype.slice.call(
        table.querySelectorAll("tr.hit-row")
      );
    }
    if (groupToggle) {
      groupToggle.addEventListener("change", function () {
        var rows = currentRows();
        if (!groupByOrganism() && !sortState.key) {
          sortState.key = "max";
          sortState.dir = -1;
          table.querySelectorAll("th.sortable").forEach(function (h) {
            h.classList.remove("sorted-asc", "sorted-desc");
            if (h.getAttribute("data-sort") === "max") {
              h.classList.add("sorted-desc");
            }
          });
        }
        rebuildTable(sortRows(rows));
      });
    }
    table.querySelectorAll("th.sortable").forEach(function (th) {
      th.addEventListener("click", function () {
        var key = th.getAttribute("data-sort");
        if (!key) return;
        if (sortState.key === key) sortState.dir *= -1;
        else { sortState.key = key; sortState.dir = 1; }
        table.querySelectorAll("th.sortable").forEach(function (h) {
          h.classList.remove("sorted-asc", "sorted-desc");
        });
        th.classList.add(sortState.dir > 0 ? "sorted-asc" : "sorted-desc");
        rebuildTable(sortRows(currentRows()));
      });
    });
  }
})();
"""


def build_online_blast_html(
    result: BlastOutputResult,
    *,
    fetch_taxonomy: bool = True,
) -> str:
    """Render a self-contained NCBI-style multi-tab HTML report.

    When ``fetch_taxonomy`` is True (default), unique organisms are resolved
    through NCBI Taxonomy so Lineage includes ranks such as Eukaryota.
    """
    primary = result.primary
    query_label = (
        (primary.query_def if primary else "")
        or result.query_name
        or "query"
    )
    query_len = primary.query_len if primary else 0
    n_hits = len(primary.hits) if primary else 0
    title = f"Online BLAST — {html.escape(result.query_name or query_label)}"
    meta_bits = [
        f"Program: {html.escape(result.program)}",
        f"Database: {html.escape(result.database or '—')}",
        f"Query length: {query_len or '—'}",
        f"Hits: {n_hits}",
    ]
    if result.rid:
        meta_bits.append(f"RID: {html.escape(result.rid)}")
    if result.version:
        meta_bits.append(html.escape(result.version))

    org_lineages: dict[str, list[str]] = {}
    taxonomy_ok = False
    if primary is not None and fetch_taxonomy:
        organisms = [_guess_organism(hit.title) for hit in primary.hits]
        try:
            org_lineages = _lookup_organism_lineages(organisms, enable_network=True)
            taxonomy_ok = bool(org_lineages)
        except Exception:  # noqa: BLE001
            org_lineages = {}
            taxonomy_ok = False

    if primary is None:
        empty = "<p class='empty'>No BLAST iterations in result.</p>"
        desc = graphic = lineage = empty
        score_legend = ""
    else:
        desc = _descriptions_table(primary, org_lineages=org_lineages)
        graphic = _graphic_svg(primary)
        lineage = _lineage_html(
            primary, org_lineages=org_lineages, taxonomy_ok=taxonomy_ok
        )
        score_legend = _alignment_score_legend_html()

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>{title}</title>
<style>{_CSS}</style>
</head>
<body>
<header class="page">
  <h1>{title}</h1>
  <p class="meta">{" · ".join(meta_bits)}</p>
</header>
<nav class="tabs" role="tablist">
  <button type="button" class="active" data-tab="tab-descriptions">Descriptions</button>
  <button type="button" data-tab="tab-graphic">Graphic Summary</button>
  <button type="button" data-tab="tab-lineage">Lineage</button>
</nav>
<section id="tab-descriptions" class="panel active" role="tabpanel">
  {desc}
</section>
<section id="tab-graphic" class="panel" role="tabpanel">
  {score_legend}
  <div class="graphic-wrap">{graphic}</div>
</section>
<section id="tab-lineage" class="panel" role="tabpanel">
  <h2>Lineage</h2>
  {lineage}
</section>
<script>{_JS}</script>
</body>
</html>
"""


def online_blast_html_path(query_fasta: Path) -> Path:
    query_fasta = Path(query_fasta).expanduser().resolve()
    return query_fasta.parent / f"{query_fasta.name}_ob.html"


def write_online_blast_html(
    result: BlastOutputResult,
    out_path: Path,
    *,
    fetch_taxonomy: bool = True,
) -> Path:
    out_path = Path(out_path).expanduser().resolve()
    out_path.write_text(
        build_online_blast_html(result, fetch_taxonomy=fetch_taxonomy),
        encoding="utf-8",
    )
    return out_path
