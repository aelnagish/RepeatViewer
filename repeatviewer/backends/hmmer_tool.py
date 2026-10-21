"""Tools → HMMER: search selected FASTA against a project HMM database (pyhmmer)."""

from __future__ import annotations

import csv
import tempfile
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from repeatviewer.backends.auto_annotation_settings import OrfSettings
from repeatviewer.backends.database_prepare import (
    database_search_readiness,
    pyhmmer_install_hint,
)
from repeatviewer.backends.pyhmmer_search import (
    collect_dna_hmm_domain_hits,
    collect_pyhmmer_domain_hits,
    hmm_uses_nucleotide_query,
    nhmmscan_available,
)
from repeatviewer.backends.sequence_operations import predicted_orf_records
from repeatviewer.job_options import (
    DEFAULT_HMMER_EVALUE,
    HmmerOptions,
    normalize_blast_evalue,
)
from repeatviewer.project import UserDatabase
from repeatviewer.tool_output_paths import kept_input_name


def _hmmer_process_log(message: str) -> None:
    from repeatviewer.external_process import emit_process_log

    emit_process_log(message)


HIT_COLUMNS = (
    "query_id",
    "hmm_name",
    "score",
    "evalue",
    "qstart",
    "qend",
)


def _pyhmmer_available() -> bool:
    try:
        import pyhmmer  # type: ignore[import-untyped]  # noqa: F401
    except ImportError:
        return False
    return True


def _pyhmmer_version_label() -> str:
    try:
        import pyhmmer  # type: ignore[import-untyped]

        return str(getattr(pyhmmer, "__version__", "?"))
    except ImportError:
        return "?"


def resolve_hmm_database(
    databases: list[UserDatabase],
    name: str,
) -> UserDatabase:
    """Return the named HMM user database or raise a clear RuntimeError."""
    wanted = (name or "").strip()
    hmm_dbs = [db for db in databases if db.is_hmm]
    if not wanted:
        if not hmm_dbs:
            raise RuntimeError(
                "No HMM database is configured. Use File → Manage Database "
                "and enable the HMM checkbox for a .hmm / Dfam profile."
            )
        raise RuntimeError(
            "Select an HMM database in Tools → HMMER → Change HMMER parameters."
        )
    for db in hmm_dbs:
        if db.name == wanted:
            problem = database_search_readiness(db)
            if problem:
                raise RuntimeError(f"{db.name}: {problem}")
            return db
    available = ", ".join(db.name for db in hmm_dbs) or "(none)"
    raise RuntimeError(
        f"HMM database {wanted!r} was not found. Available: {available}."
    )


def _fasta_looks_nucleotide(path: Path) -> bool:
    """Heuristic: mostly ACGTN and not short peptide-like sequences."""
    try:
        record = next(SeqIO.parse(str(path), "fasta"))
    except StopIteration:
        return True
    seq = str(record.seq).upper().replace("U", "T")
    if not seq:
        return True
    letters = [c for c in seq if c.isalpha()]
    if not letters:
        return True
    acgtn = sum(1 for c in letters if c in "ACGTN")
    return (acgtn / len(letters)) >= 0.85 and len(letters) >= 60


def _orf_protein_query(
    query_path: Path,
    *,
    max_orfs_per_sequence: int,
    work_dir: Path,
) -> Path:
    """Translate ORFs to a temporary protein FASTA for protein-HMM search."""
    records = predicted_orf_records(
        query_path,
        translate=True,
        settings=OrfSettings(),
    )
    if max_orfs_per_sequence > 0:
        by_source: dict[str, list[SeqRecord]] = defaultdict(list)
        for rec in records:
            source = "query"
            for token in (rec.description or "").split():
                if token.startswith("source="):
                    source = token.split("=", 1)[1]
                    break
            by_source[source].append(rec)
        capped: list[SeqRecord] = []
        for group in by_source.values():
            group.sort(key=lambda r: len(r.seq), reverse=True)
            capped.extend(group[:max_orfs_per_sequence])
        records = capped
    out = work_dir / f"{kept_input_name(query_path)}_hmmer_orfs.faa"
    SeqIO.write(records, str(out), "fasta")
    return out


def _select_hits(
    raw_hits: list[tuple[float, float | None, str, str, int, int]],
    *,
    max_hits: int,
    minimum_bit_score: float,
) -> list[tuple[float, float | None, str, str, int, int]]:
    filtered = [
        hit
        for hit in raw_hits
        if float(hit[0]) >= float(minimum_bit_score)
    ]
    filtered.sort(key=lambda h: (-float(h[0]), h[1] if h[1] is not None else 1e9))
    return filtered[: max(1, int(max_hits))]


def next_hmmer_table_path(out_dir: Path, query_path: Path) -> Path:
    """``{query}_hmmer.txt``, then ``{query}_hmmer1.txt``, ``{query}_hmmer2.txt``, …"""
    out_dir = out_dir.expanduser().resolve()
    stem = f"{kept_input_name(query_path)}_hmmer"
    candidate = out_dir / f"{stem}.txt"
    if not candidate.exists():
        return candidate
    index = 1
    while True:
        candidate = out_dir / f"{stem}{index}.txt"
        if not candidate.exists():
            return candidate
        index += 1


def write_hmmer_hits_table(
    hits: list[tuple[float, float | None, str, str, int, int]],
    output_path: Path,
) -> Path:
    """Write a TSV hit table with a header row."""
    output_path = output_path.expanduser().resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(HIT_COLUMNS)
        for score, evalue, query_id, hmm_name, qstart, qend in hits:
            writer.writerow(
                [
                    query_id,
                    hmm_name,
                    f"{float(score):.3f}",
                    "" if evalue is None else f"{float(evalue):.3g}",
                    int(qstart),
                    int(qend),
                ]
            )
    return output_path


def run_hmmer_search(
    query_path: Path,
    database: UserDatabase,
    *,
    options: HmmerOptions | None = None,
    cancel=None,
) -> list[Path]:
    """Search ``query_path`` against ``database`` and optionally write a hit table."""
    nucleotide_preview = hmm_uses_nucleotide_query(
        database.path.expanduser().resolve()
    )
    if nucleotide_preview and nhmmscan_available():
        pass  # system nhmmscan is enough for DNA HMM DBs such as Dfam
    elif not _pyhmmer_available():
        raise RuntimeError(
            f"pyhmmer is not installed ({pyhmmer_install_hint()})."
        )
    opts = options or HmmerOptions()
    hmm_path = database.path.expanduser().resolve()
    if not hmm_path.is_file():
        raise RuntimeError(f"HMM file not found: {hmm_path}")

    report_e = normalize_blast_evalue(opts.e_value, default=DEFAULT_HMMER_EVALUE)
    domain_e = normalize_blast_evalue(
        opts.domain_e_value, default=DEFAULT_HMMER_EVALUE
    )
    query_path = query_path.expanduser().resolve()
    out_dir = query_path.parent
    created: list[Path] = []

    _hmmer_process_log("—— HMMER (pyhmmer) ——")
    _hmmer_process_log(f"pyhmmer {_pyhmmer_version_label()}")
    _hmmer_process_log(f"Database: {database.name} ({hmm_path})")
    _hmmer_process_log(f"Query: {query_path.name}")
    _hmmer_process_log(
        f"Thresholds: E≤{report_e:g}, domain E≤{domain_e:g}, "
        f"min bit score≥{float(opts.minimum_bit_score):g}, "
        f"max hits={int(opts.max_hits)}"
    )
    _hmmer_process_log(f"Threads: {max(1, int(opts.num_threads))}")

    with tempfile.TemporaryDirectory(prefix="repeatviewer_hmmer_") as tmp:
        tmp_dir = Path(tmp)
        search_query = query_path
        nucleotide_hmm = hmm_uses_nucleotide_query(hmm_path)
        if not nucleotide_hmm and _fasta_looks_nucleotide(query_path):
            _hmmer_process_log(
                "Protein HMM on nucleotide query — predicting ORFs for search…"
            )
            search_query = _orf_protein_query(
                query_path,
                max_orfs_per_sequence=opts.max_orfs_per_sequence,
                work_dir=tmp_dir,
            )
            orf_count = sum(1 for _ in SeqIO.parse(str(search_query), "fasta"))
            cap = int(opts.max_orfs_per_sequence)
            cap_note = f" (max {cap} per sequence)" if cap > 0 else ""
            _hmmer_process_log(
                f"ORF search set: {orf_count} sequence(s){cap_note} → {search_query.name}"
            )
        else:
            mode = "nucleotide" if nucleotide_hmm else "protein"
            _hmmer_process_log(f"Search mode: {mode} query against HMM profile(s)")

        if nucleotide_hmm:
            from repeatviewer.backends.pyhmmer_search import hmm_pressed_indexes_present

            use_nhmmscan = nhmmscan_available() and hmm_pressed_indexes_present(hmm_path)
            engine = "nhmmscan" if use_nhmmscan else "pyhmmer hmmscan"
            _hmmer_process_log(f"Running {engine}…")
            raw_hits = collect_dna_hmm_domain_hits(
                hmm_path,
                search_query,
                report_e=report_e,
                domain_e=domain_e,
                cpus=max(1, int(opts.num_threads)),
                cancel=cancel,
            )
        else:
            _hmmer_process_log("Running hmmscan…")
            raw_hits = collect_pyhmmer_domain_hits(
                hmm_path,
                search_query,
                report_e=report_e,
                domain_e=domain_e,
                cpus=max(1, int(opts.num_threads)),
            )

    _hmmer_process_log(
        f"HMM search finished: {len(raw_hits)} domain hit(s) before filtering"
    )
    selected = _select_hits(
        raw_hits,
        max_hits=opts.max_hits,
        minimum_bit_score=opts.minimum_bit_score,
    )
    _hmmer_process_log(f"After filters: {len(selected)} hit(s) kept")
    if selected:
        best = selected[0]
        _hmmer_process_log(
            f"Top hit: {best[2]} vs {best[3]} "
            f"(score {float(best[0]):.3f}, E={best[1] if best[1] is not None else '—'})"
        )
    if opts.write_text:
        table = next_hmmer_table_path(out_dir, query_path)
        write_hmmer_hits_table(selected, table)
        created.append(table.resolve())
        _hmmer_process_log(f"Wrote: {table.name}")
    elif not selected:
        _hmmer_process_log("No hits — hit table not written (export disabled or empty)")
    _hmmer_process_log("Done.")
    return created
