"""Search and download sequences from NCBI (Tools > Download from NCBI).

Three NCBI databases are offered:

* Nucleotide (``nuccore``) and Protein (``protein``): searched with Entrez
  E-utilities (esearch + esummary) and fetched with efetch as FASTA or
  GenBank / GenPept text. Nucleotide FASTA can also carry the NCBI features
  (efetch feature table, converted to GFF3 and imported as annotations).
* Genome Assembly (``assembly``): searched with E-utilities, downloaded as a
  zip from the NCBI Datasets API (genome FASTA, optionally its GFF3).

Pure Python (``urllib``), no Biopython. Everything here runs off the UI
thread; ``progress(done, total, detail)`` and ``cancel()`` let the caller
drive the progress bar and the Stop button.

Developer: Abdullah Yusuf
"""

from __future__ import annotations

import json
import re
import shutil
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterable

EUTILS_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
DATASETS_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/"
NCBI_WEB_URL = "https://www.ncbi.nlm.nih.gov/"
_TOOL = "RepeatViewer"

DB_NUCLEOTIDE = "nuccore"
DB_PROTEIN = "protein"
DB_ASSEMBLY = "assembly"
DATABASES = (DB_NUCLEOTIDE, DB_PROTEIN, DB_ASSEMBLY)
DATABASE_LABELS = {
    DB_NUCLEOTIDE: "Nucleotide",
    DB_PROTEIN: "Protein",
    DB_ASSEMBLY: "Genome Assembly",
}

FORMAT_FASTA = "fasta"
FORMAT_FASTA_FEATURES = "fasta_features"
FORMAT_GENBANK = "genbank"
FORMAT_GENOME = "genome"
FORMAT_GENOME_GFF = "genome_gff"
FORMATS_BY_DB: dict[str, tuple[tuple[str, str], ...]] = {
    DB_NUCLEOTIDE: (
        (FORMAT_FASTA, "FASTA (.fa)"),
        (FORMAT_FASTA_FEATURES, "FASTA with NCBI annotations (.fa)"),
        (FORMAT_GENBANK, "GenBank (.gb)"),
    ),
    DB_PROTEIN: (
        (FORMAT_FASTA, "FASTA (.faa)"),
        (FORMAT_GENBANK, "GenPept (.gp)"),
    ),
    DB_ASSEMBLY: (
        (FORMAT_GENOME, "Genome FASTA (.fa)"),
        (FORMAT_GENOME_GFF, "Genome FASTA + GFF3 annotation (.fa, .gff3)"),
    ),
}

MOLECULE_ANY = ""
MOLECULE_GENOMIC = "genomic"
MOLECULE_MRNA = "mrna"
MOLECULE_RRNA = "rrna"
MOLECULE_LABELS = (
    (MOLECULE_ANY, "Any"),
    (MOLECULE_GENOMIC, "Genomic DNA"),
    (MOLECULE_MRNA, "mRNA"),
    (MOLECULE_RRNA, "rRNA"),
)

DEFAULT_PAGE_SIZE = 100
MAX_PAGE_SIZE = 500
# Records per efetch call when several go into one file.
_EFETCH_BATCH = 20
_CHUNK = 1 << 20
_RETRIES = 3

# Accession shapes (https://www.ncbi.nlm.nih.gov/genbank/acc_prefix/).
_RE_ASSEMBLY = re.compile(r"^GC[AF]_\d{9}(\.\d+)?$", re.I)
_RE_REFSEQ_PROTEIN = re.compile(r"^(AP|NP|XP|YP|WP)_\d+(\.\d+)?$", re.I)
_RE_REFSEQ_NUC = re.compile(r"^(AC|NC|NG|NM|NR|NT|NW|NZ|XM|XR)_\d+(\.\d+)?$", re.I)
_RE_GENBANK_PROTEIN = re.compile(r"^[A-Z]{3}\d{5,7}(\.\d+)?$", re.I)
_RE_GENBANK_NUC = re.compile(r"^[A-Z]{1,2}\d{5,8}(\.\d+)?$", re.I)
_RE_WGS = re.compile(r"^[A-Z]{4,6}\d{8,11}(\.\d+)?$", re.I)


class NcbiError(RuntimeError):
    """NCBI could not be reached or returned something unusable."""


class NcbiCancelled(RuntimeError):
    """The user pressed Stop."""


# --------------------------------------------------------------------- query


def split_accessions(text: str) -> list[str]:
    return [t for t in re.split(r"[\s,;]+", str(text or "").strip()) if t]


def accession_database(token: str) -> str | None:
    """NCBI database an accession belongs to, or ``None`` if it is not one."""
    t = token.strip()
    if _RE_ASSEMBLY.match(t):
        return DB_ASSEMBLY
    if _RE_REFSEQ_PROTEIN.match(t) or _RE_GENBANK_PROTEIN.match(t):
        return DB_PROTEIN
    if _RE_REFSEQ_NUC.match(t) or _RE_GENBANK_NUC.match(t) or _RE_WGS.match(t):
        return DB_NUCLEOTIDE
    return None


def looks_like_accessions(text: str) -> list[str]:
    """Accessions typed or pasted in the search box; [] when it is a free-text search."""
    tokens = split_accessions(text)
    if not tokens:
        return []
    if all(accession_database(t) for t in tokens):
        return tokens
    return []


def guess_database(text: str) -> str | None:
    """Database every accession in ``text`` agrees on (None if mixed or free text)."""
    tokens = looks_like_accessions(text)
    dbs = {accession_database(t) for t in tokens}
    return dbs.pop() if len(dbs) == 1 else None


@dataclass
class NcbiQuery:
    database: str = DB_NUCLEOTIDE
    text: str = ""
    organism: str = ""
    min_length: int = 0
    max_length: int = 0
    molecule: str = MOLECULE_ANY
    refseq_only: bool = False
    latest_only: bool = True


def _quote(term: str) -> str:
    term = term.replace('"', "").strip()
    return f'"{term}"'


def build_term(q: NcbiQuery) -> str:
    """Entrez query string for ``q`` (accessions, or text plus filters)."""
    accs = looks_like_accessions(q.text)
    if accs:
        if q.database == DB_ASSEMBLY:
            return " OR ".join(accs)
        return " OR ".join(f"{a}[accn]" for a in accs)
    parts: list[str] = []
    text = str(q.text or "").strip()
    if text:
        parts.append(f"({text})" if " OR " in text.upper() else text)
    if q.organism.strip():
        parts.append(f"{_quote(q.organism)}[Organism]")
    if q.database in (DB_NUCLEOTIDE, DB_PROTEIN):
        lo, hi = max(0, int(q.min_length)), max(0, int(q.max_length))
        if lo or hi:
            parts.append(f"{lo or 1}:{hi or 999999999999}[SLEN]")
        if q.refseq_only:
            parts.append("refseq[filter]")
    if q.database == DB_NUCLEOTIDE and q.molecule:
        parts.append(f"biomol_{q.molecule}[PROP]")
    if q.database == DB_ASSEMBLY:
        if q.refseq_only:
            parts.append('"latest refseq"[filter]')
        elif q.latest_only:
            parts.append("latest[filter]")
    if not parts:
        raise NcbiError("Type an accession or a search term first.")
    return " AND ".join(parts)


# ------------------------------------------------------------------- records


@dataclass
class NcbiRecord:
    uid: str
    accession: str
    title: str = ""
    organism: str = ""
    length: int = 0
    kind: str = ""
    updated: str = ""
    extra: dict = field(default_factory=dict)

    def web_url(self, database: str) -> str:
        if database == DB_ASSEMBLY:
            return f"{NCBI_WEB_URL}datasets/genome/{self.accession}/"
        return f"{NCBI_WEB_URL}{database}/{self.accession}"


@dataclass
class SearchPage:
    count: int
    records: list[NcbiRecord]
    retstart: int
    term: str
    translated: str = ""


def _to_int(value, default: int = 0) -> int:
    try:
        return int(str(value).strip())
    except (TypeError, ValueError):
        return default


def _date(value) -> str:
    text = str(value or "").strip()
    return text[:10].replace("/", "-") if text else ""


def _assembly_stat(meta: str, category: str) -> int:
    m = re.search(
        rf'<Stat category="{category}" sequence_tag="all">(\d+)</Stat>', str(meta or "")
    )
    return int(m.group(1)) if m else 0


def record_from_summary(database: str, doc: dict) -> NcbiRecord:
    """One esummary JSON document → :class:`NcbiRecord`."""
    uid = str(doc.get("uid", ""))
    if database == DB_ASSEMBLY:
        meta = doc.get("meta", "")
        return NcbiRecord(
            uid=uid,
            accession=str(doc.get("assemblyaccession", "") or uid),
            title=str(doc.get("assemblyname", "")),
            organism=str(doc.get("speciesname") or doc.get("organism") or ""),
            length=_assembly_stat(meta, "total_length"),
            kind=str(doc.get("assemblystatus", "")),
            updated=_date(doc.get("asmreleasedate_genbank") or doc.get("seqreleasedate")),
            extra={
                "assembly_name": str(doc.get("assemblyname", "")),
                "refseq_category": str(doc.get("refseq_category", "")),
                "contig_n50": _to_int(doc.get("contign50")),
                "scaffold_n50": _to_int(doc.get("scaffoldn50")),
                "sequences": _assembly_stat(meta, "scaffold_count")
                or _assembly_stat(meta, "contig_count"),
                "submitter": str(doc.get("submitterorganization", "")),
            },
        )
    biomol = str(doc.get("biomol", "") or "")
    moltype = str(doc.get("moltype", "") or "")
    kind = biomol or moltype
    if database == DB_PROTEIN:
        kind = "protein"
    elif biomol == "genomic" and moltype:
        kind = f"genomic {moltype.upper()}" if moltype in ("dna", "rna") else "genomic"
    topology = str(doc.get("topology", "") or "")
    if topology == "circular":
        kind = f"{kind}, circular" if kind else "circular"
    return NcbiRecord(
        uid=uid,
        accession=str(doc.get("accessionversion") or doc.get("caption") or uid),
        title=str(doc.get("title", "")),
        organism=str(doc.get("organism", "")),
        length=_to_int(doc.get("slen")),
        kind=kind,
        updated=_date(doc.get("updatedate")),
        extra={"taxid": _to_int(doc.get("taxid"))},
    )


# -------------------------------------------------------------------- client


class NcbiClient:
    """Small E-utilities / Datasets client with NCBI's request-rate limits.

    NCBI allows 3 requests per second without an API key and 10 with one.
    One client is shared by every call of a search or download.
    """

    _lock = threading.Lock()
    _last_request = 0.0

    def __init__(self, *, email: str = "", api_key: str = "", timeout: float = 60.0) -> None:
        self.email = str(email or "").strip()
        self.api_key = str(api_key or "").strip()
        self.timeout = float(timeout)

    @property
    def _interval(self) -> float:
        return 0.11 if self.api_key else 0.34

    def _wait_turn(self) -> None:
        with NcbiClient._lock:
            now = time.monotonic()
            delay = NcbiClient._last_request + self._interval - now
            if delay > 0:
                time.sleep(delay)
            NcbiClient._last_request = time.monotonic()

    def _eutils_params(self, params: dict) -> dict:
        out = {k: str(v) for k, v in params.items() if v not in (None, "")}
        out["tool"] = _TOOL
        if self.email:
            out["email"] = self.email
        if self.api_key:
            out["api_key"] = self.api_key
        return out

    def open(self, url: str, *, data: dict | None = None, headers: dict | None = None):
        """urlopen with rate limit and retries on 429 / 5xx / network errors."""
        body = urllib.parse.urlencode(data).encode("utf-8") if data is not None else None
        hdrs = {"User-Agent": f"{_TOOL}/1.0"}
        hdrs.update(headers or {})
        last: Exception | None = None
        for attempt in range(_RETRIES):
            self._wait_turn()
            req = urllib.request.Request(url, data=body, headers=hdrs,
                                         method="POST" if body is not None else "GET")
            try:
                return urllib.request.urlopen(req, timeout=self.timeout)
            except urllib.error.HTTPError as exc:
                last = exc
                if exc.code in (429, 500, 502, 503, 504) and attempt + 1 < _RETRIES:
                    time.sleep(1.5 * (attempt + 1))
                    continue
                detail = ""
                try:
                    detail = exc.read().decode("utf-8", "replace")[:300]
                except Exception:  # noqa: BLE001
                    pass
                raise NcbiError(f"NCBI returned HTTP {exc.code}: {detail or exc.reason}") from exc
            except (urllib.error.URLError, TimeoutError, ConnectionError) as exc:
                last = exc
                if attempt + 1 < _RETRIES:
                    time.sleep(1.5 * (attempt + 1))
                    continue
                reason = getattr(exc, "reason", exc)
                raise NcbiError(f"Could not reach NCBI (internet required): {reason}") from exc
        raise NcbiError(f"NCBI request failed: {last}")

    def eutils(self, endpoint: str, params: dict, *, post: bool = False):
        url = EUTILS_URL + endpoint
        full = self._eutils_params(params)
        if post:
            return self.open(url, data=full)
        return self.open(f"{url}?{urllib.parse.urlencode(full)}")

    def eutils_json(self, endpoint: str, params: dict, *, post: bool = False) -> dict:
        with self.eutils(endpoint, {**params, "retmode": "json"}, post=post) as resp:
            text = resp.read().decode("utf-8", "replace")
        try:
            data = json.loads(text)
        except ValueError as exc:
            raise NcbiError(f"NCBI sent an unreadable answer: {text[:200]}") from exc
        if isinstance(data, dict) and data.get("error"):
            raise NcbiError(f"NCBI: {data['error']}")
        return data

    def datasets(self, path: str, params: dict | None = None):
        url = DATASETS_URL + path
        if params:
            url += "?" + urllib.parse.urlencode(params, doseq=True)
        headers = {"api-key": self.api_key} if self.api_key else {}
        return self.open(url, headers=headers)


# -------------------------------------------------------------------- search


def summaries(client: NcbiClient, database: str, uids: list[str]) -> list[NcbiRecord]:
    if not uids:
        return []
    data = client.eutils_json(
        "esummary.fcgi", {"db": database, "id": ",".join(uids)}, post=len(uids) > 50
    )
    result = data.get("result") or {}
    records = []
    for uid in result.get("uids") or uids:
        doc = result.get(str(uid))
        if isinstance(doc, dict) and not doc.get("error"):
            records.append(record_from_summary(database, doc))
    order = {u: i for i, u in enumerate(uids)}
    records.sort(key=lambda r: order.get(r.uid, len(order)))
    return records


def search(client: NcbiClient, query: NcbiQuery, *, retstart: int = 0,
           retmax: int = DEFAULT_PAGE_SIZE) -> SearchPage:
    """One page of results (esearch + esummary)."""
    term = build_term(query)
    retmax = max(1, min(int(retmax), MAX_PAGE_SIZE))
    data = client.eutils_json(
        "esearch.fcgi",
        {"db": query.database, "term": term, "retstart": int(retstart), "retmax": retmax},
    )
    res = data.get("esearchresult") or {}
    errors = res.get("errorlist") or {}
    phrases = [p for p in (errors.get("phrasesnotfound") or []) if p]
    count = _to_int(res.get("count"))
    ids = [str(i) for i in res.get("idlist") or []]
    if not ids and phrases and count == 0:
        raise NcbiError(f"Nothing found for: {', '.join(phrases)}")
    return SearchPage(
        count=count,
        records=summaries(client, query.database, ids),
        retstart=int(retstart),
        term=term,
        translated=str(res.get("querytranslation", "")),
    )


# ------------------------------------------------------------------ download


def safe_file_stem(text: str, *, limit: int = 60) -> str:
    """File-name-safe stem (no spaces: RepeatViewer tools dislike them)."""
    stem = re.sub(r"[^A-Za-z0-9._-]+", "_", str(text or "").strip()).strip("._-")
    stem = re.sub(r"_+", "_", stem)
    return stem[:limit].rstrip("._-") or "ncbi_download"


def output_suffix(database: str, fmt: str) -> str:
    if database == DB_PROTEIN:
        return ".gp" if fmt == FORMAT_GENBANK else ".faa"
    if fmt == FORMAT_GENBANK:
        return ".gb"
    return ".fa"


def default_file_name(database: str, fmt: str, records: list[NcbiRecord], query_text: str) -> str:
    """Name for a combined download: the accession for one record, else the search."""
    suffix = output_suffix(database, fmt)
    if len(records) == 1:
        return record_file_name(database, fmt, records[0])
    accs = looks_like_accessions(query_text)
    if accs:
        stem = f"{safe_file_stem(accs[0])}_and_{len(records) - 1}_more"
    else:
        stem = safe_file_stem(query_text or "ncbi_download")
    return stem + suffix


def record_file_name(database: str, fmt: str, record: NcbiRecord) -> str:
    if database == DB_ASSEMBLY:
        name = record.extra.get("assembly_name") or record.title
        stem = safe_file_stem(f"{record.accession}_{name}" if name else record.accession)
        return stem + ".fa"
    return safe_file_stem(record.accession) + output_suffix(database, fmt)


def planned_outputs(database: str, fmt: str, records: list[NcbiRecord], *,
                    one_file: bool, file_name: str, out_dir: Path) -> list[Path]:
    """Files a download will write (used for the file-exists check up front)."""
    out_dir = Path(out_dir)
    if database == DB_ASSEMBLY:
        paths = []
        for rec in records:
            fa = out_dir / record_file_name(database, fmt, rec)
            paths.append(fa)
            if fmt == FORMAT_GENOME_GFF:
                paths.append(fa.with_suffix(".gff3"))
        return paths
    if one_file and len(records) > 1:
        return [out_dir / file_name]
    return [out_dir / record_file_name(database, fmt, r) for r in records]


ProgressFn = Callable[[int, int, str], None]
CancelFn = Callable[[], bool]


def _check(cancel: CancelFn | None) -> None:
    if cancel is not None and cancel():
        raise NcbiCancelled("Cancelled")


def _stream(resp, handle, cancel: CancelFn | None, on_bytes: Callable[[int], None] | None = None) -> int:
    total = 0
    while True:
        _check(cancel)
        chunk = resp.read(_CHUNK)
        if not chunk:
            return total
        handle.write(chunk)
        total += len(chunk)
        if on_bytes is not None:
            on_bytes(total)


def _efetch_rettype(database: str, fmt: str) -> str:
    if fmt == FORMAT_GENBANK:
        return "gp" if database == DB_PROTEIN else "gbwithparts"
    return "fasta"


def _efetch_to(client: NcbiClient, database: str, uids: list[str], rettype: str,
               handle, cancel: CancelFn | None) -> int:
    with client.eutils(
        "efetch.fcgi",
        {"db": database, "id": ",".join(uids), "rettype": rettype, "retmode": "text"},
        post=True,
    ) as resp:
        return _stream(resp, handle, cancel)


def _fetch_feature_table(client: NcbiClient, uids: list[str]) -> str:
    with client.eutils(
        "efetch.fcgi",
        {"db": DB_NUCLEOTIDE, "id": ",".join(uids), "rettype": "ft", "retmode": "text"},
        post=True,
    ) as resp:
        return resp.read().decode("utf-8", "replace")


def _fmt_bytes(n: int) -> str:
    size = float(n)
    for unit in ("bytes", "KB", "MB", "GB"):
        if size < 1024 or unit == "GB":
            return f"{size:.0f} {unit}" if unit == "bytes" else f"{size:.1f} {unit}"
        size /= 1024
    return f"{n} bytes"


@dataclass
class DownloadResult:
    files: list[Path]
    fasta_with_features: list[Path] = field(default_factory=list)
    feature_tables: dict = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)


def download_records(
    client: NcbiClient,
    database: str,
    records: list[NcbiRecord],
    out_dir: Path,
    *,
    fmt: str = FORMAT_FASTA,
    one_file: bool = True,
    file_name: str = "",
    targets: dict[str, Path] | None = None,
    progress: ProgressFn | None = None,
    cancel: CancelFn | None = None,
) -> DownloadResult:
    """Download ``records`` into ``out_dir``.

    ``targets`` maps a planned output path (as from :func:`planned_outputs`)
    to the path chosen after the file-exists question (Overwrite / Rename).
    For ``FORMAT_FASTA_FEATURES`` the NCBI feature tables are returned as GFF3
    text in ``feature_tables`` (FASTA path → GFF3); importing them as
    annotations is up to the caller.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    targets = {str(k): Path(v) for k, v in (targets or {}).items()}
    records = list(records)
    if not records:
        raise NcbiError("Tick at least one record to download.")

    def _target(path: Path) -> Path:
        return targets.get(str(path), path)

    def _emit(done: int, detail: str) -> None:
        if progress is not None:
            progress(done, len(records), detail)

    if database == DB_ASSEMBLY:
        return _download_assemblies(client, records, out_dir, fmt=fmt,
                                    target=_target, emit=_emit, cancel=cancel)

    rettype = _efetch_rettype(database, fmt)
    result = DownloadResult(files=[])
    with_features = database == DB_NUCLEOTIDE and fmt == FORMAT_FASTA_FEATURES
    if one_file and len(records) > 1:
        planned = out_dir / (file_name or default_file_name(database, fmt, records, ""))
        final = _target(planned)
        part = final.with_name(final.name + ".part")
        done = 0
        try:
            with part.open("wb") as handle:
                for start in range(0, len(records), _EFETCH_BATCH):
                    batch = records[start:start + _EFETCH_BATCH]
                    _emit(done, f"{batch[0].accession} ({done} done, {len(records) - done} left)")
                    _efetch_to(client, database, [r.uid for r in batch], rettype, handle, cancel)
                    done += len(batch)
            part.replace(final)
        finally:
            part.unlink(missing_ok=True)
        result.files.append(final)
        if with_features:
            result.feature_tables[str(final)] = _features_gff3(client, records, cancel)
            result.fasta_with_features.append(final)
        _emit(len(records), f"{len(records)} done, 0 left")
        return result

    for index, rec in enumerate(records):
        _check(cancel)
        _emit(index, f"{rec.accession} ({index} done, {len(records) - index} left)")
        final = _target(out_dir / record_file_name(database, fmt, rec))
        part = final.with_name(final.name + ".part")
        try:
            with part.open("wb") as handle:
                _efetch_to(client, database, [rec.uid], rettype, handle, cancel)
            part.replace(final)
        finally:
            part.unlink(missing_ok=True)
        result.files.append(final)
        if with_features:
            result.feature_tables[str(final)] = _features_gff3(client, [rec], cancel)
            result.fasta_with_features.append(final)
    _emit(len(records), f"{len(records)} done, 0 left")
    return result


def _features_gff3(client: NcbiClient, records: list[NcbiRecord], cancel: CancelFn | None) -> str:
    chunks = []
    for start in range(0, len(records), _EFETCH_BATCH):
        _check(cancel)
        batch = records[start:start + _EFETCH_BATCH]
        chunks.append(_fetch_feature_table(client, [r.uid for r in batch]))
    return feature_table_to_gff3("\n".join(chunks))


def _download_assemblies(client: NcbiClient, records: list[NcbiRecord], out_dir: Path, *,
                         fmt: str, target, emit, cancel: CancelFn | None) -> DownloadResult:
    result = DownloadResult(files=[])
    want_gff = fmt == FORMAT_GENOME_GFF
    total = len(records)
    for index, rec in enumerate(records):
        _check(cancel)
        left = f"({index} done, {total - index} left)"
        emit(index, f"{rec.accession} {left}")
        include = ["GENOME_FASTA"] + (["GENOME_GFF"] if want_gff else [])
        zip_path = out_dir / f".{safe_file_stem(rec.accession)}.ncbi_download.zip"

        def _bytes(n: int, acc=rec.accession, i=index) -> None:
            emit(i, f"{acc}: {_fmt_bytes(n)} downloaded ({i} done, {total - i} left)")

        try:
            with client.datasets(
                f"genome/accession/{urllib.parse.quote(rec.accession)}/download",
                {"include_annotation_type": include},
            ) as resp, zip_path.open("wb") as handle:
                _stream(resp, handle, cancel, _bytes)
            _check(cancel)
            emit(index, f"{rec.accession}: unpacking {left}")
            fasta_planned = out_dir / record_file_name(DB_ASSEMBLY, fmt, rec)
            fasta_out = target(fasta_planned)
            gff_out = target(fasta_planned.with_suffix(".gff3")) if want_gff else None
            wrote_fa, wrote_gff = extract_assembly_zip(zip_path, fasta_out, gff_out, cancel=cancel)
        finally:
            zip_path.unlink(missing_ok=True)
        if not wrote_fa:
            result.warnings.append(f"{rec.accession}: NCBI sent no genome FASTA.")
            continue
        result.files.append(fasta_out)
        if want_gff:
            if wrote_gff and gff_out is not None:
                result.files.append(gff_out)
            else:
                result.warnings.append(f"{rec.accession}: this assembly has no NCBI annotation (GFF3).")
    emit(total, f"{total} done, 0 left")
    return result


def extract_assembly_zip(zip_path: Path, fasta_out: Path, gff_out: Path | None, *,
                         cancel: CancelFn | None = None) -> tuple[bool, bool]:
    """Pull the genome FASTA (and GFF3) out of an NCBI Datasets zip."""
    wrote_fa = wrote_gff = False
    try:
        archive = zipfile.ZipFile(zip_path)
    except zipfile.BadZipFile as exc:
        raise NcbiError("NCBI Datasets sent a damaged download; try again.") from exc
    with archive:
        names = archive.namelist()
        fastas = sorted(n for n in names if n.startswith("ncbi_dataset/data/")
                        and n.lower().endswith((".fna", ".fa", ".fasta")))
        gffs = sorted(n for n in names if n.startswith("ncbi_dataset/data/")
                      and n.lower().endswith((".gff", ".gff3")))
        if fastas:
            part = fasta_out.with_name(fasta_out.name + ".part")
            try:
                with part.open("wb") as out:
                    for name in fastas:
                        with archive.open(name) as src:
                            while True:
                                _check(cancel)
                                chunk = src.read(_CHUNK)
                                if not chunk:
                                    break
                                out.write(chunk)
                part.replace(fasta_out)
                wrote_fa = True
            finally:
                part.unlink(missing_ok=True)
        if gff_out is not None and gffs:
            with archive.open(gffs[0]) as src, gff_out.open("wb") as out:
                shutil.copyfileobj(src, out, _CHUNK)
            wrote_gff = True
    return wrote_fa, wrote_gff


# --------------------------------------------------- feature table → GFF3

_FT_SKIP_TYPES = {"source"}
_FT_NAME_KEYS = (
    "gene", "product", "locus_tag", "standard_name", "label", "rpt_family",
    "function", "regulatory_class", "ncRNA_class", "note",
)
_GFF_ESCAPE = {";": "%3B", "=": "%3D", "&": "%26", ",": "%2C", "\t": " ", "\n": " "}


def _gff_escape(text: str) -> str:
    return "".join(_GFF_ESCAPE.get(ch, ch) for ch in str(text))


def _ft_seqid(header: str) -> str:
    """``>Feature ref|NC_045512.2|`` → ``NC_045512.2``."""
    token = header[len(">Feature"):].strip().split()[0] if header.strip() != ">Feature" else ""
    parts = [p for p in token.split("|") if p]
    if len(parts) >= 2:
        return parts[1]
    return parts[0] if parts else token


def _ft_pos(text: str) -> int | None:
    text = text.strip().lstrip("<>")
    if "^" in text:
        text = text.split("^")[0]
    return int(text) if text.isdigit() else None


@dataclass
class _FtFeature:
    seqid: str
    ftype: str
    spans: list[tuple[int, int]] = field(default_factory=list)
    strand: str = "+"
    quals: list[tuple[str, str]] = field(default_factory=list)


def parse_feature_table(text: str) -> list[_FtFeature]:
    """Parse NCBI's 5-column feature table (efetch ``rettype=ft``)."""
    features: list[_FtFeature] = []
    seqid = ""
    current: _FtFeature | None = None
    for raw in str(text or "").splitlines():
        if not raw.strip():
            continue
        if raw.startswith(">Feature"):
            seqid = _ft_seqid(raw)
            current = None
            continue
        cols = raw.split("\t")
        if raw.startswith("\t\t\t"):
            if current is not None:
                key = cols[3].strip() if len(cols) > 3 else ""
                value = cols[4].strip() if len(cols) > 4 else ""
                if key:
                    current.quals.append((key, value))
            continue
        if len(cols) < 2:
            continue
        a, b = _ft_pos(cols[0]), _ft_pos(cols[1])
        ftype = cols[2].strip() if len(cols) > 2 else ""
        if a is None or b is None:
            continue
        strand = "-" if a > b else "+"
        span = (min(a, b), max(a, b))
        if ftype:
            current = _FtFeature(seqid=seqid, ftype=ftype, spans=[span], strand=strand)
            features.append(current)
        elif current is not None:
            current.spans.append(span)
    return features


def feature_table_to_gff3(text: str) -> str:
    """NCBI feature table → GFF3 (one row per span; split features share an ID)."""
    lines = ["##gff-version 3"]
    counter = 0
    for feat in parse_feature_table(text):
        if feat.ftype in _FT_SKIP_TYPES or not feat.seqid:
            continue
        counter += 1
        quals = dict(feat.quals)
        name = next((quals[k] for k in _FT_NAME_KEYS if quals.get(k)), feat.ftype)
        attrs = [f"ID=ncbi{counter}", f"Name={_gff_escape(name)}"]
        for key in ("gene", "product", "locus_tag", "protein_id", "note", "rpt_family", "rpt_type"):
            if quals.get(key) and key != "name":
                attrs.append(f"{key}={_gff_escape(quals[key])}")
        attr = ";".join(attrs)
        for start, end in feat.spans:
            lines.append("\t".join([
                feat.seqid, "NCBI", feat.ftype, str(start), str(end), ".", feat.strand, ".", attr,
            ]))
    return "\n".join(lines) + "\n"


def clipboard_accessions(text: str) -> list[str]:
    """Accessions found anywhere in pasted text (for right-click download)."""
    found: list[str] = []
    for token in re.split(r"[^A-Za-z0-9_.]+", str(text or "")):
        token = token.strip(".")
        if token and accession_database(token) and token not in found:
            found.append(token)
    return found


def iter_unique(values: Iterable[str]) -> list[str]:
    seen: list[str] = []
    for v in values:
        if v not in seen:
            seen.append(v)
    return seen
