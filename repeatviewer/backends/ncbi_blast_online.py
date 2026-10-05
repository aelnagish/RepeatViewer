"""NCBI online BLAST (Common URL API / QBlast).

Submit a FASTA query to NCBI, poll until ready, download XML results, write a
self-contained HTML report (``*_ob.html``), and optionally a tabular hit table
(``*_ob.txt``). Online BLAST no longer paints hits on the Annotation track.

Note: NCBI's web API ``FORMAT_TYPE=Tabular`` currently returns an empty READY
stub, so results are fetched as classic BlastOutput XML and converted.
"""

from __future__ import annotations

import re
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Callable

from repeatviewer.backends.job_cancel import JobCancellation, check_cancelled

NCBI_BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
# NCBI asks clients to identify themselves.
_NCBI_TOOL = "RepeatViewer"
_NCBI_EMAIL = "repeatviewer-noreply@local"

# Common public databases offered in the Online BLAST tab.
# Keys are NCBI ``DATABASE`` API names; labels match the web BLAST wording.
_ONLINE_NT_DATABASES = (
    ("core_nt", "core_nt — Core nucleotide (NCBI default)"),
    ("nt", "nt — Nucleotide collection (nr/nt)"),
    ("refseq_select_rna", "refseq_select_rna — RefSeq Select RNA"),
    ("refseq_rna", "refseq_rna — Reference RNA sequences"),
    ("refseq_reference_genomes", "refseq_reference_genomes — RefSeq reference genomes"),
    ("refseq_genomes", "refseq_genomes — RefSeq Genome Database"),
    ("refseq_genomic", "refseq_genomic — RefSeq genomic sequences"),
    ("wgs", "wgs — Whole-genome shotgun contigs"),
    ("est", "est — Expressed sequence tags"),
    ("tsa", "tsa — Transcriptome Shotgun Assembly"),
    ("htgs", "htgs — High-throughput genomic sequences"),
    ("gss", "gss — Genomic survey sequences"),
    ("sts", "sts — Sequence tagged sites"),
    ("patnt", "patnt — Patent nucleotide sequences"),
    ("pdbnt", "pdbnt — PDB nucleotide structures"),
    ("RefSeq_Gene", "RefSeq_Gene — Human RefSeqGene"),
    ("tls", "tls — Targeted Loci (TLS)"),
    ("sra", "sra — Sequence Read Archive (large / slow)"),
)
_ONLINE_PROTEIN_DATABASES = (
    ("nr", "nr — Non-redundant protein sequences"),
    ("refseq_select_prot", "refseq_select_prot — RefSeq Select proteins"),
    ("refseq_protein", "refseq_protein — NCBI RefSeq proteins"),
    ("swissprot", "swissprot — UniProtKB/Swiss-Prot"),
    ("pdbaa", "pdbaa — PDB protein structures"),
    ("landmark", "landmark — Landmark (model organisms)"),
    ("pataa", "pataa — Patent protein sequences"),
)
ONLINE_BLAST_DATABASES: dict[str, tuple[tuple[str, str], ...]] = {
    "blastn": _ONLINE_NT_DATABASES,
    "blastp": _ONLINE_PROTEIN_DATABASES,
    "blastx": _ONLINE_PROTEIN_DATABASES,
    # Protein query vs translated nucleotide DB.
    "tblastn": _ONLINE_NT_DATABASES,
    # Translated nucleotide query vs translated nucleotide DB.
    "tblastx": _ONLINE_NT_DATABASES,
}

_DEFAULT_ONLINE_DB = {
    "blastn": "core_nt",
    "blastp": "nr",
    "blastx": "nr",
    "tblastn": "core_nt",
    "tblastx": "core_nt",
}

# Web BLAST is not meant for huge submissions.
_MAX_QUERY_CHARS = 100_000
_POLL_SECONDS = 3.0
_MAX_WAIT_SECONDS = 30 * 60

# Match local BLAST ``-outfmt 6`` columns used elsewhere in RepeatViewer,
# plus ``stitle`` (NCBI Hit_def: species / function description).
_TABULAR_COLUMNS = (
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "sstrand",
    "evalue",
    "bitscore",
    "stitle",
)


class NcbiBlastOnlineError(RuntimeError):
    """NCBI online BLAST failed or returned no usable result."""


def default_online_database(program: str) -> str:
    return _DEFAULT_ONLINE_DB.get(str(program).lower(), "nt")


def _user_agent() -> str:
    return f"{_NCBI_TOOL}/1.0 ({_NCBI_EMAIL})"


def _ncbi_get(params: dict[str, str], *, timeout: float = 60.0) -> str:
    query = urllib.parse.urlencode(params)
    req = urllib.request.Request(
        f"{NCBI_BLAST_URL}?{query}",
        headers={"User-Agent": _user_agent()},
        method="GET",
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return resp.read().decode("utf-8", errors="replace")
    except urllib.error.HTTPError as exc:
        body = exc.read().decode("utf-8", errors="replace") if exc.fp else ""
        raise NcbiBlastOnlineError(
            f"NCBI BLAST HTTP {exc.code}: {body[:300] or exc.reason}"
        ) from exc
    except urllib.error.URLError as exc:
        raise NcbiBlastOnlineError(
            f"Could not reach NCBI BLAST (internet required): {exc.reason}"
        ) from exc


def _ncbi_post(params: dict[str, str], *, timeout: float = 120.0) -> str:
    data = urllib.parse.urlencode(params).encode("utf-8")
    req = urllib.request.Request(
        NCBI_BLAST_URL,
        data=data,
        headers={
            "User-Agent": _user_agent(),
            "Content-Type": "application/x-www-form-urlencoded",
        },
        method="POST",
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return resp.read().decode("utf-8", errors="replace")
    except urllib.error.HTTPError as exc:
        body = exc.read().decode("utf-8", errors="replace") if exc.fp else ""
        raise NcbiBlastOnlineError(
            f"NCBI BLAST HTTP {exc.code}: {body[:300] or exc.reason}"
        ) from exc
    except urllib.error.URLError as exc:
        raise NcbiBlastOnlineError(
            f"Could not reach NCBI BLAST (internet required): {exc.reason}"
        ) from exc


def _parse_rid(text: str) -> tuple[str, int]:
    rid_match = re.search(r"^\s*RID\s*=\s*(\S+)", text, re.MULTILINE)
    if not rid_match:
        raise NcbiBlastOnlineError(
            "NCBI BLAST did not return a request ID (RID).\n"
            f"Response starts with:\n{text[:400]}"
        )
    rtoe = 10
    rtoe_match = re.search(r"^\s*RTOE\s*=\s*(\d+)", text, re.MULTILINE)
    if rtoe_match:
        rtoe = max(1, int(rtoe_match.group(1)))
    return rid_match.group(1), rtoe


def _parse_status(text: str) -> str:
    match = re.search(r"^\s*Status\s*=\s*(\S+)", text, re.MULTILINE | re.IGNORECASE)
    if match:
        return match.group(1).upper()
    # Older / HTML-ish replies.
    if "Status=READY" in text or "status=ready" in text.lower():
        return "READY"
    if "Status=WAITING" in text or "status=waiting" in text.lower():
        return "WAITING"
    if "Status=FAILED" in text or "status=failed" in text.lower():
        return "FAILED"
    return "UNKNOWN"


def _read_query_fasta(path: Path) -> tuple[str, int, int]:
    """Return (fasta_text, sequence_count, sequence_characters)."""
    text = path.expanduser().resolve().read_text(encoding="utf-8", errors="replace")
    cleaned = text.strip()
    if not cleaned:
        raise NcbiBlastOnlineError(f"Query FASTA is empty: {path.name}")
    seq_count = sum(1 for line in cleaned.splitlines() if line.startswith(">"))
    seq_chars = sum(
        len(line.strip())
        for line in cleaned.splitlines()
        if line.strip() and not line.startswith(">")
    )
    if seq_chars > _MAX_QUERY_CHARS:
        raise NcbiBlastOnlineError(
            f"Query is too long for NCBI web BLAST ({seq_chars:,} bp; "
            f"limit about {_MAX_QUERY_CHARS:,}). Use local BLAST for large inputs."
        )
    return cleaned, max(1, seq_count), seq_chars


def _xml_text(node: ET.Element | None, tag: str, default: str = "") -> str:
    if node is None:
        return default
    child = node.find(tag)
    if child is None or child.text is None:
        return default
    return str(child.text).strip()


def _short_subject_id(hit: ET.Element) -> str:
    """Prefer accession (XM_… / NR_…) over long ``gi|…|ref|…|`` ids."""
    accession = _xml_text(hit, "Hit_accession")
    if accession:
        return accession.split()[0]
    hit_id = _xml_text(hit, "Hit_id")
    if hit_id:
        # gi|2518588887|ref|XM_010669769.3|  →  XM_010669769.3
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
            # Last non-empty token is usually the accession.
            return parts[-1]
        return hit_id.split()[0]
    hit_def = _xml_text(hit, "Hit_def", "subject")
    return hit_def.split()[0] if hit_def else "subject"


def _hit_stitle(hit: ET.Element) -> str:
    """NCBI description line (species + function); safe for one TSV field."""
    raw = _xml_text(hit, "Hit_def")
    if not raw:
        return ""
    # Collapse whitespace / strip characters that would break tabular rows.
    return " ".join(raw.replace("\t", " ").split())


def _annotation_target_label(accession: str, stitle: str) -> str:
    """Readable Annotation-track label: accession + species/function."""
    acc = (accession or "").strip()
    title = (stitle or "").strip()
    if title and acc:
        if acc in title:
            return title
        return f"{acc} {title}"
    return title or acc or "hit"


def blast_xml_to_tabular(xml_text: str) -> tuple[str, int]:
    """Convert NCBI BlastOutput XML to local BLAST tabular text.

    Returns ``(table_text, hit_row_count)``. Rows use the same columns as
    RepeatViewer's local ``-outfmt 6`` BLAST tables, plus ``stitle``.
    """
    text = (xml_text or "").strip()
    if not text:
        return "", 0
    # Strip a possible HTML wrapper around the XML declaration.
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

    rows: list[str] = []
    for iteration in root.findall(".//Iteration"):
        qseqid = _xml_text(iteration, "Iteration_query-def") or _xml_text(
            iteration, "Iteration_query-ID", "query"
        )
        # Keep only the first token (FASTA id before spaces).
        qseqid = qseqid.split()[0] if qseqid else "query"
        for hit in iteration.findall("Iteration_hits/Hit"):
            sseqid = _short_subject_id(hit)
            stitle = _hit_stitle(hit)
            for hsp in hit.findall("Hit_hsps/Hsp"):
                identity = float(_xml_text(hsp, "Hsp_identity", "0") or 0)
                align_len = int(float(_xml_text(hsp, "Hsp_align-len", "0") or 0))
                gaps = int(float(_xml_text(hsp, "Hsp_gaps", "0") or 0))
                qstart = int(float(_xml_text(hsp, "Hsp_query-from", "0") or 0))
                qend = int(float(_xml_text(hsp, "Hsp_query-to", "0") or 0))
                sstart = int(float(_xml_text(hsp, "Hsp_hit-from", "0") or 0))
                send = int(float(_xml_text(hsp, "Hsp_hit-to", "0") or 0))
                evalue = _xml_text(hsp, "Hsp_evalue", "0") or "0"
                bitscore = _xml_text(hsp, "Hsp_bit-score", "0") or "0"
                pident = (100.0 * identity / align_len) if align_len else 0.0
                mismatch = max(0, align_len - int(identity) - gaps)
                if sstart <= send:
                    sstrand = "plus"
                else:
                    sstrand = "minus"
                    sstart, send = send, sstart
                rows.append(
                    "\t".join(
                        [
                            qseqid,
                            sseqid,
                            f"{pident:.3f}",
                            str(align_len),
                            str(mismatch),
                            str(gaps),
                            str(qstart),
                            str(qend),
                            str(sstart),
                            str(send),
                            sstrand,
                            evalue,
                            bitscore,
                            stitle,
                        ]
                    )
                )
    table = "\n".join(rows) + ("\n" if rows else "")
    return table, len(rows)


def run_ncbi_online_blast(
    query_fasta: Path,
    out_dir: Path,
    *,
    program: str = "blastn",
    database: str = "nt",
    evalue: float = 0.0001,
    hitlist_size: int = 100,
    write_text: bool = True,
    cancel: JobCancellation | None = None,
    progress: Callable[[int, int, str], None] | None = None,
) -> list[Path]:
    """Submit one FASTA to NCBI BLAST; always write ``*_ob.html``.

    Optionally also writes the tabular ``*_ob.txt`` hit table when
    ``write_text`` is True. Returns the list of created output paths
    (HTML first).
    """
    program = str(program or "blastn").lower()
    if program not in ONLINE_BLAST_DATABASES:
        raise NcbiBlastOnlineError(
            f"Unsupported online BLAST program: {program} "
            f"(use blastn, blastp, blastx, tblastn, or tblastx)."
        )
    database = str(database or default_online_database(program)).strip() or (
        default_online_database(program)
    )
    query_path = Path(query_fasta).expanduser().resolve()
    if not query_path.is_file():
        raise NcbiBlastOnlineError(f"Query FASTA not found: {query_path}")
    out_dir = Path(out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    def _emit(done: int, total: int, detail: str) -> None:
        # Wired to the job progress bar and (when enabled) the process window.
        if progress is not None:
            progress(done, total, detail)

    check_cancelled(cancel)
    _emit(0, 1, f"{query_path.name}: connecting to NCBI BLAST…")
    query_text, seq_count, seq_chars = _read_query_fasta(query_path)
    _emit(
        0,
        1,
        f"{query_path.name}: prepared query "
        f"({seq_count} sequence(s), {seq_chars:,} letters); "
        f"program={program}, database={database}, evalue={evalue}",
    )
    _emit(0, 1, f"{query_path.name}: submitting search to NCBI servers…")

    # Do not set FORMAT_TYPE on Put — choose format on Get (XML).
    put_params = {
        "CMD": "Put",
        "PROGRAM": program,
        "DATABASE": database,
        "QUERY": query_text,
        "EXPECT": str(float(evalue)),
        "HITLIST_SIZE": str(max(1, min(500, int(hitlist_size)))),
        "TOOL": _NCBI_TOOL,
        "EMAIL": _NCBI_EMAIL,
    }
    # NCBI may rewrite DATABASE=nt to core_nt unless this flag is set.
    if database == "nt" and program in {"blastn", "tblastn", "tblastx"}:
        put_params["NO_DATABASE_OVERRIDE"] = "true"
    put_body = _ncbi_post(put_params)
    check_cancelled(cancel)
    rid, rtoe = _parse_rid(put_body)
    _emit(
        0,
        1,
        f"{query_path.name}: submitted — NCBI RID {rid} "
        f"(estimated wait ~{rtoe}s)",
    )

    # Honour NCBI's estimated time before the first status poll.
    wait_first = min(float(rtoe), 60.0)
    deadline = time.monotonic() + _MAX_WAIT_SECONDS
    elapsed = 0.0
    while elapsed < wait_first:
        check_cancelled(cancel)
        remaining = int(wait_first - elapsed)
        if remaining == int(wait_first) or remaining % 5 == 0:
            _emit(
                0,
                1,
                f"{query_path.name}: waiting for NCBI queue "
                f"(RID {rid}; ~{max(remaining, 0)}s)…",
            )
        step = min(1.0, wait_first - elapsed)
        time.sleep(step)
        elapsed += step

    status = "WAITING"
    poll = 0
    while time.monotonic() < deadline:
        check_cancelled(cancel)
        poll += 1
        _emit(
            0,
            1,
            f"{query_path.name}: checking NCBI status "
            f"(RID {rid}; poll {poll})…",
        )
        info = _ncbi_get(
            {
                "CMD": "Get",
                "RID": rid,
                "FORMAT_OBJECT": "SearchInfo",
            }
        )
        status = _parse_status(info)
        if status == "READY":
            _emit(0, 1, f"{query_path.name}: NCBI search finished (RID {rid})")
            break
        if status == "FAILED" or (
            status == "UNKNOWN" and "FAILED" in info.upper()
        ):
            raise NcbiBlastOnlineError(
                f"NCBI BLAST search failed (RID {rid}).\n{info[:400]}"
            )
        _emit(
            0,
            1,
            f"{query_path.name}: NCBI still searching "
            f"(RID {rid}; status={status})…",
        )
        time.sleep(_POLL_SECONDS)
    else:
        raise NcbiBlastOnlineError(
            f"NCBI BLAST timed out after {_MAX_WAIT_SECONDS // 60} minutes "
            f"(RID {rid}). Try again later or use local BLAST."
        )

    check_cancelled(cancel)
    _emit(0, 1, f"{query_path.name}: downloading NCBI XML results…")
    raw = _ncbi_get(
        {
            "CMD": "Get",
            "RID": rid,
            "FORMAT_TYPE": "XML",
            "HITLIST_SIZE": str(max(1, min(500, int(hitlist_size)))),
        },
        timeout=180.0,
    )
    if "<BlastOutput" not in raw:
        raise NcbiBlastOnlineError(
            f"NCBI BLAST returned no BlastOutput XML (RID {rid}).\n"
            f"Response starts with:\n{raw[:400]}"
        )

    from repeatviewer.backends.online_blast_report import (
        parse_blast_output_xml,
        write_online_blast_html,
    )

    _emit(0, 1, f"{query_path.name}: building NCBI HTML report…")
    parsed = parse_blast_output_xml(
        raw,
        rid=rid,
        query_name=query_path.name,
        program=program,
        database=database,
    )
    html_path = out_dir / f"{query_path.name}_ob.html"
    write_online_blast_html(parsed, html_path)
    created: list[Path] = [html_path.resolve()]

    table, n_hits = blast_xml_to_tabular(raw)
    if write_text:
        _emit(0, 1, f"{query_path.name}: converting NCBI XML to hit table…")
        out_path = out_dir / f"{query_path.name}_ob.txt"
        header_lines = [
            f"# NCBI online BLAST ({program} vs {database}); RID={rid}",
            f"# query={query_path.name}",
            f"# columns={','.join(_TABULAR_COLUMNS)}",
        ]
        if n_hits == 0:
            header_lines.append(
                "# No significant hits returned "
                "(try a higher -evalue, another database, or local BLAST)."
            )
        out_path.write_text(
            "\n".join(header_lines) + "\n" + table, encoding="utf-8"
        )
        created.append(out_path.resolve())
        _emit(
            0,
            1,
            f"{query_path.name}: wrote {out_path.name} ({n_hits} hit row(s))",
        )
    else:
        # Drop a stale table from a previous run when the user opted out.
        stale = out_dir / f"{query_path.name}_ob.txt"
        stale.unlink(missing_ok=True)

    if n_hits == 0:
        _emit(
            0,
            1,
            f"{query_path.name}: NCBI returned 0 hits "
            f"(RID {rid}; evalue={evalue}, database={database})",
        )
    else:
        _emit(
            0,
            1,
            f"{query_path.name}: NCBI returned {n_hits} hit row(s) (RID {rid})",
        )
    _emit(1, 1, f"{query_path.name}: wrote {html_path.name}")
    return created


# ---------------------------------------------------------------------------
# Annotation cleanup (legacy Online BLAST track features)
# ---------------------------------------------------------------------------

ONLINE_BLAST_SOURCE = "ncbi_online_blast"
ONLINE_BLAST_FEATURE_TYPE = "domain"
ONLINE_BLAST_COLOR = "#1F77B4"  # distinct track colour for NCBI online hits
ONLINE_BLAST_LABEL = "NCBI"


def online_blast_output_path(query_fasta: Path) -> Path:
    """Canonical ``{fasta_basename}_ob.txt`` next to the query."""
    query_fasta = Path(query_fasta).expanduser().resolve()
    return query_fasta.parent / f"{query_fasta.name}_ob.txt"


def online_blast_html_output_path(query_fasta: Path) -> Path:
    """Canonical ``{fasta_basename}_ob.html`` next to the query."""
    from repeatviewer.backends.online_blast_report import online_blast_html_path

    return online_blast_html_path(query_fasta)


def _parse_online_blast_rows(ob_txt: Path) -> list[dict[str, str]]:
    """Parse RepeatViewer online-BLAST tabular rows (skip ``#`` comments)."""
    rows: list[dict[str, str]] = []
    text = Path(ob_txt).expanduser().resolve().read_text(
        encoding="utf-8", errors="replace"
    )
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        parts = stripped.split("\t")
        if len(parts) < 12:
            continue
        rows.append(
            {
                "qseqid": parts[0],
                "sseqid": parts[1],
                "pident": parts[2],
                "length": parts[3],
                "qstart": parts[6],
                "qend": parts[7],
                "sstart": parts[8],
                "send": parts[9],
                "sstrand": parts[10] if len(parts) > 10 else "plus",
                "evalue": parts[11] if len(parts) > 11 else "",
                "bitscore": parts[12] if len(parts) > 12 else "",
                # Description may contain spaces; keep the remainder joined.
                "stitle": (
                    parts[13]
                    if len(parts) == 14
                    else " ".join(parts[13:]) if len(parts) > 14 else ""
                ),
            }
        )
    return rows


def online_blast_features_by_query(
    ob_txt: Path,
    *,
    max_hits: int | None = None,
) -> dict[str, list]:
    """Build ``SequenceFeature`` lists keyed by query id from an ``*_ob.txt``.

    Coordinates are query (consensus) positions — for the Annotation tab, not
    genome maps. Hits are grouped by subject accession; ``max_hits`` caps the
    number of **subjects** (all HSPs for a kept subject are retained).
    """
    from repeatviewer.backends.sequence_features import SequenceFeature
    from repeatviewer.backends.similarity_annotation import similarity_feature_name

    limit = None if max_hits is None else max(1, int(max_hits))
    # (qid, accession) -> {best_bitscore, feats}
    groups: dict[tuple[str, str], dict] = {}
    for hit in _parse_online_blast_rows(ob_txt):
        try:
            qstart = int(float(hit["qstart"]))
            qend = int(float(hit["qend"]))
        except (KeyError, TypeError, ValueError):
            continue
        if qstart <= 0 or qend <= 0:
            continue
        strand = "+"
        if qstart > qend:
            qstart, qend = qend, qstart
            strand = "-"
        sstrand = str(hit.get("sstrand", "plus")).lower()
        if sstrand in {"minus", "-", "rev", "c"}:
            strand = "-"
        try:
            bitscore = float(hit.get("bitscore") or 0.0)
        except (TypeError, ValueError):
            bitscore = 0.0
        try:
            pident = float(hit.get("pident") or 0.0)
        except (TypeError, ValueError):
            pident = 0.0
        qid = str(hit.get("qseqid") or "query").split()[0]
        accession = str(hit.get("sseqid") or "").strip().split()[0] or "subject"
        target = _annotation_target_label(
            accession,
            str(hit.get("stitle") or ""),
        )
        parent = f"{ONLINE_BLAST_LABEL}:{accession}"
        feat = SequenceFeature(
            start=qstart,
            end=qend,
            name=similarity_feature_name(
                ONLINE_BLAST_LABEL,
                target,
                show_database_name=True,
                evalue=hit.get("evalue"),
            ),
            feature_type=ONLINE_BLAST_FEATURE_TYPE,
            color=ONLINE_BLAST_COLOR,
            strand=strand,
            source=ONLINE_BLAST_SOURCE,
            parent=parent,
            pident=pident,
            subject_id=accession,
        )
        key = (qid, accession)
        bucket = groups.get(key)
        if bucket is None:
            groups[key] = {"best_bitscore": bitscore, "feats": [feat]}
        else:
            bucket["feats"].append(feat)
            if bitscore > float(bucket["best_bitscore"]):
                bucket["best_bitscore"] = bitscore

    # Rank subject groups by best HSP bitscore; keep top ``limit`` subjects.
    ranked_groups = sorted(
        groups.items(),
        key=lambda item: float(item[1]["best_bitscore"]),
        reverse=True,
    )
    if limit is not None:
        ranked_groups = ranked_groups[:limit]

    by_query: dict[str, list] = {}
    for (qid, _acc), bucket in ranked_groups:
        by_query.setdefault(qid, []).extend(bucket["feats"])
    return by_query


def enrich_online_blast_features_from_ob(
    features: list,
    ob_txt: Path,
) -> list:
    """Patch legacy online features missing pident/parent from ``*_ob.txt``.

    Matches by query coordinates + subject accession when possible. Returns a
    new list (unchanged entries when no match). Does not write to disk.
    """
    from repeatviewer.backends.sequence_features import SequenceFeature

    online = [
        f
        for f in features
        if str(getattr(f, "source", "") or "") == ONLINE_BLAST_SOURCE
        or str(getattr(f, "name", "") or "").startswith(f"{ONLINE_BLAST_LABEL}:")
    ]
    if not online:
        return list(features)
    if all(float(getattr(f, "pident", 0) or 0) > 0 for f in online):
        return list(features)
    table = Path(ob_txt)
    if not table.is_file():
        return list(features)

    # Index rows by (start, end) and by accession for fuzzy match.
    by_coords: dict[tuple[int, int], list[dict[str, str]]] = {}
    for hit in _parse_online_blast_rows(table):
        try:
            qstart = int(float(hit["qstart"]))
            qend = int(float(hit["qend"]))
        except (KeyError, TypeError, ValueError):
            continue
        if qstart > qend:
            qstart, qend = qend, qstart
        by_coords.setdefault((qstart, qend), []).append(hit)

    updated: list = []
    for feat in features:
        is_online = (
            str(feat.source or "") == ONLINE_BLAST_SOURCE
            or str(feat.name or "").startswith(f"{ONLINE_BLAST_LABEL}:")
        )
        if not is_online or float(getattr(feat, "pident", 0) or 0) > 0:
            updated.append(feat)
            continue
        candidates = by_coords.get((int(feat.start), int(feat.end))) or []
        hit = candidates[0] if candidates else None
        if hit is None:
            updated.append(feat)
            continue
        accession = str(hit.get("sseqid") or "").strip().split()[0]
        try:
            pident = float(hit.get("pident") or 0.0)
        except (TypeError, ValueError):
            pident = 0.0
        updated.append(
            SequenceFeature(
                start=feat.start,
                end=feat.end,
                name=feat.name,
                feature_type=feat.feature_type,
                color=feat.color or ONLINE_BLAST_COLOR,
                strand=feat.strand,
                track=feat.track,
                parent=feat.parent or (f"{ONLINE_BLAST_LABEL}:{accession}" if accession else ""),
                source=ONLINE_BLAST_SOURCE,
                pident=pident,
                subject_id=feat.subject_id or accession,
            )
        )
    return updated


def _match_record_id(qseqid: str, record_ids: list[str]) -> str | None:
    from repeatviewer.backends.sequence_features import (
        fasta_id_token,
        strip_duplicate_uniquifier,
    )

    token = fasta_id_token(qseqid) or qseqid
    if token in record_ids:
        return token
    for rid in record_ids:
        original = strip_duplicate_uniquifier(rid, sibling_ids=record_ids)
        if original == token or fasta_id_token(original) == token:
            return rid
        if rid.split("#", 1)[0] == token or original.split("#", 1)[0] == token:
            return rid
    if len(record_ids) == 1:
        return record_ids[0]
    return None


def merge_online_blast_into_annotations(
    fasta_path: Path,
    ob_txt: Path,
    *,
    max_hits: int | None = None,
) -> list[Path]:
    """Deprecated: Online BLAST no longer paints the Annotation track.

    Kept for compatibility; strips prior NCBI online features and ignores
    ``ob_txt`` / ``max_hits``. Prefer :func:`strip_online_blast_from_annotations`.
    """
    del ob_txt, max_hits
    return strip_online_blast_from_annotations(fasta_path)


def strip_online_blast_from_annotations(fasta_path: Path) -> list[Path]:
    """Remove legacy NCBI Online BLAST features from seqmap JSON / GFF.

    Returns paths of rewritten annotation sidecars (empty if none existed).
    """
    from repeatviewer.backends.sequence_features import (
        SequenceFeature,
        load_features_by_id_from_json,
        load_sequence_map_data,
        save_annotations_multi,
    )

    fasta = Path(fasta_path).expanduser().resolve()
    if not fasta.is_file():
        raise NcbiBlastOnlineError(f"Query FASTA not found: {fasta}")

    data = load_sequence_map_data(fasta)
    record_ids = list(data.record_ids)
    lengths = dict(data.lengths_by_id)
    if not record_ids:
        return []

    json_feats = load_features_by_id_from_json(fasta, record_ids)
    features_by_id: dict[str, list[SequenceFeature]] = {
        rid: list(json_feats.get(rid) or data.features_by_id.get(rid) or [])
        for rid in record_ids
    }
    changed = False
    for rid in record_ids:
        before = features_by_id[rid]
        filtered = [
            feat
            for feat in before
            if str(feat.source or "") != ONLINE_BLAST_SOURCE
            and not str(feat.name or "").startswith(f"{ONLINE_BLAST_LABEL}:")
        ]
        if len(filtered) != len(before):
            changed = True
        features_by_id[rid] = filtered

    if not changed:
        return []

    json_path, gff_path = save_annotations_multi(
        fasta,
        features_by_id=features_by_id,
        lengths_by_id=lengths,
        record_ids=record_ids,
    )
    outputs = [json_path.resolve()]
    if gff_path is not None:
        outputs.append(gff_path.resolve())
    return outputs
