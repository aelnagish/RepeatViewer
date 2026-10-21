"""pyhmmer / nhmmscan helpers for protein and DNA HMM databases (e.g. Pfam, Dfam)."""

from __future__ import annotations

import shutil
import tempfile
from pathlib import Path


def _parse_hmm_alphabet_from_text(hmm_path: Path) -> str | None:
    """Read the ALPH line from an HMMER3 profile (e.g. Dfam uses DNA)."""
    try:
        with hmm_path.open(encoding="utf-8", errors="replace") as handle:
            for _ in range(80):
                line = handle.readline()
                if not line:
                    break
                stripped = line.strip()
                if stripped.startswith("ALPH"):
                    parts = stripped.split()
                    if len(parts) >= 2:
                        return parts[1].lower()
                if stripped == "//":
                    break
    except OSError:
        return None
    return None


def detect_hmm_alphabet(hmm_path: Path) -> str:
    """Return ``dna``, ``rna``, or ``amino`` for the first profile in the file."""
    hmm_path = hmm_path.expanduser().resolve()
    parsed = _parse_hmm_alphabet_from_text(hmm_path)
    if parsed:
        if parsed in {"dna", "dna/rna"}:
            return "dna"
        if parsed == "rna":
            return "rna"
        if parsed in {"amino", "protein"}:
            return "amino"
        return parsed

    try:
        from pyhmmer.plan7 import HMMFile  # type: ignore[import-untyped]

        with HMMFile(str(hmm_path)) as hmm_file:
            for hmm in hmm_file:
                alphabet = getattr(hmm, "alphabet", None)
                if alphabet is not None:
                    label = str(getattr(alphabet, "name", alphabet)).lower()
                    if "dna" in label:
                        return "dna"
                    if "rna" in label:
                        return "rna"
                    return "amino"
    except Exception:
        pass
    return "amino"


def hmm_uses_nucleotide_query(hmm_path: Path) -> bool:
    """Whether hmmscan should receive nucleotide rather than protein queries."""
    return detect_hmm_alphabet(hmm_path) in {"dna", "rna"}


def _alphabet_from_kind(kind: str):
    """Return a pyhmmer.easel.Alphabet instance for a known kind.

    kind: "dna", "rna", or "amino"
    """
    import pyhmmer  # type: ignore[import-untyped]

    from pyhmmer.easel import Alphabet  # type: ignore[import-untyped]

    kind = (kind or "").lower()
    if kind == "dna":
        candidate = getattr(Alphabet, "dna", None)
    elif kind == "rna":
        candidate = getattr(Alphabet, "rna", None)
    else:
        candidate = getattr(Alphabet, "amino", None) or getattr(Alphabet, "protein", None)

    if candidate is None:
        raise ValueError(f"Unsupported alphabet kind: {kind!r}")

    # pyhmmer exposes some alphabets as callables (e.g. Alphabet.amino()).
    return candidate() if callable(candidate) else candidate


def _decode_name(value: object, *, default: str = "") -> str:
    if value is None:
        return default
    if isinstance(value, bytes):
        return value.decode()
    return str(value)


def _pyhmmer_version_tuple() -> tuple[int, int]:
    import pyhmmer  # type: ignore[import-untyped]

    parts = pyhmmer.__version__.split(".")
    major = int(parts[0]) if parts else 0
    minor = int(parts[1]) if len(parts) > 1 else 0
    return major, minor


def _read_digital_queries(query_path: Path, *, alphabet: object | None = None) -> list[object]:
    from pyhmmer.easel import SequenceFile  # type: ignore[import-untyped]

    kwargs: dict[str, object] = {"digital": True}
    if alphabet is not None:
        kwargs["alphabet"] = alphabet

    with SequenceFile(str(query_path), **kwargs) as seq_file:
        return list(seq_file)


def _scan_kwargs(
    report_e: float,
    *,
    domain_e: float | None = None,
) -> dict[str, float | str]:
    seq_e = float(report_e)
    dom_e = float(domain_e) if domain_e is not None else seq_e
    return {
        "E": seq_e,
        "domE": dom_e,
        "incE": seq_e,
        "incdomE": dom_e,
    }


def _query_name_from_top_hits(top_hits: object) -> str:
    query = getattr(top_hits, "query", None)
    if query is not None:
        return _decode_name(getattr(query, "name", None), default="query")
    return _decode_name(getattr(top_hits, "query_name", None), default="query")


def _domain_reported_evalue(domain: object) -> float | None:
    for attr in ("i_evalue", "evalue"):
        try:
            value = float(getattr(domain, attr))
        except (AttributeError, TypeError, ValueError):
            continue
        if value >= 0:
            return value
    return None


def _append_domain_hits(
    raw_hits: list[tuple[float, float | None, str, str, int, int]],
    *,
    qname: str,
    hit: object,
) -> None:
    hmm_name = _decode_name(getattr(hit, "name", None), default="hit")
    for domain in hit.domains:
        raw_hits.append(
            (
                float(domain.score),
                _domain_reported_evalue(domain),
                qname,
                hmm_name,
                int(domain.env_from) + 1,
                int(domain.env_to),
            )
        )


def _log_scan(message: str) -> None:
    from repeatviewer.external_process import emit_process_log

    emit_process_log(message)


def _collect_with_modern_pyhmmer(
    hmm_path: Path,
    query_path: Path,
    *,
    report_e: float,
    domain_e: float | None = None,
    cpus: int = 1,
) -> list[tuple[float, float | None, str, str, int, int]]:
    import pyhmmer  # type: ignore[import-untyped]

    from pyhmmer.hmmer import hmmscan  # type: ignore[import-untyped]
    from pyhmmer.plan7 import HMMFile  # type: ignore[import-untyped]

    # pyhmmer's digitization can fail when it can't auto-detect query alphabet.
    kind = detect_hmm_alphabet(hmm_path)
    queries: list[object]
    try:
        alphabet = _alphabet_from_kind(kind)
        queries = _read_digital_queries(query_path, alphabet=alphabet)
    except Exception:
        # Fall back to pyhmmer's auto-detection.
        queries = _read_digital_queries(query_path)
    if not queries:
        _log_scan("No sequences in query file.")
        return []

    raw_hits: list[tuple[float, float | None, str, str, int, int]] = []
    thread_count = max(1, int(cpus))
    _log_scan(f"Loaded {len(queries)} query sequence(s) from {query_path.name}")
    with HMMFile(str(hmm_path)) as hmm_file:
        try:
            profiles = hmm_file.optimized_profiles()
        except ValueError:
            # Unpressed .hmm (no .h3*): optimize profiles in memory.
            profiles = list(hmm_file)
        for top_hits in hmmscan(
            queries,
            profiles,
            cpus=thread_count,
            **_scan_kwargs(report_e, domain_e=domain_e),
        ):
            qname = _query_name_from_top_hits(top_hits)
            before = len(raw_hits)
            for hit in top_hits:
                _append_domain_hits(raw_hits, qname=qname, hit=hit)
            added = len(raw_hits) - before
            _log_scan(f"  {qname}: {added} domain hit(s)")
    return raw_hits


def _collect_with_legacy_pyhmmer(
    hmm_path: Path,
    query_path: Path,
    *,
    report_e: float,
    domain_e: float | None = None,
    cpus: int = 1,
) -> list[tuple[float, float | None, str, str, int, int]]:
    import pyhmmer  # type: ignore[import-untyped]

    del cpus  # legacy easy API does not expose a cpus kwarg
    raw_hits: list[tuple[float, float | None, str, str, int, int]] = []
    _log_scan("Using pyhmmer easy API (single-threaded scan)")
    with pyhmmer.easy.HMMFile(str(hmm_path)) as hmm_db:
        with pyhmmer.easy.SeqFile(str(query_path)) as seq_file:
            scan_kwargs = _scan_kwargs(report_e, domain_e=domain_e)
            for top_hits in pyhmmer.hmmscan(hmm_db, seq_file, **scan_kwargs):
                qname = _query_name_from_top_hits(top_hits)
                before = len(raw_hits)
                for hit in top_hits:
                    _append_domain_hits(raw_hits, qname=qname, hit=hit)
                added = len(raw_hits) - before
                _log_scan(f"  {qname}: {added} domain hit(s)")
    return raw_hits


def nhmmscan_available() -> bool:
    """Whether the HMMER ``nhmmscan`` binary is on PATH (best for large DNA DBs like Dfam)."""
    return shutil.which("nhmmscan") is not None


def hmm_pressed_indexes_present(hmm_path: Path) -> bool:
    """True when hmmpress wrote ``.h3m/.h3i/.h3f/.h3p`` beside the profile."""
    hmm_path = hmm_path.expanduser().resolve()
    return all(Path(f"{hmm_path}{suffix}").is_file() for suffix in (".h3m", ".h3i", ".h3f", ".h3p"))


def parse_nhmmscan_tblout(
    text: str,
) -> list[tuple[float, float | None, str, str, int, int]]:
    """Parse ``nhmmscan --tblout`` lines into (score, evalue, query, hmm, env_from, env_to)."""
    hits: list[tuple[float, float | None, str, str, int, int]] = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        # target accession query accession hmmfrom hmmto alifrom alito
        # envfrom envto modlen strand E-value score bias …
        if len(parts) < 15:
            continue
        try:
            hmm_name = parts[0]
            query_id = parts[2]
            env_from = int(parts[8])
            env_to = int(parts[9])
            evalue = float(parts[12])
            score = float(parts[13])
        except (TypeError, ValueError, IndexError):
            continue
        hits.append((score, evalue, query_id, hmm_name, env_from, env_to))
    return hits


def collect_nhmmscan_domain_hits(
    hmm_path: Path,
    query_path: Path,
    *,
    report_e: float,
    cpus: int = 1,
    cancel=None,
) -> list[tuple[float, float | None, str, str, int, int]]:
    """Scan a DNA query against a DNA HMM DB with system ``nhmmscan`` (multi-CPU)."""
    from repeatviewer.external_process import run_logged_process

    hmm_path = hmm_path.expanduser().resolve()
    query_path = query_path.expanduser().resolve()
    cpu_n = max(1, int(cpus))
    with tempfile.TemporaryDirectory(prefix="repeatviewer_nhmmscan_") as tmp:
        tblout = Path(tmp) / "hits.tbl"
        cmd = [
            "nhmmscan",
            "--cpu",
            str(cpu_n),
            "-E",
            f"{float(report_e):g}",
            "--incE",
            f"{float(report_e):g}",
            "--tblout",
            str(tblout),
            "--noali",
            str(hmm_path),
            str(query_path),
        ]
        result = run_logged_process(
            cmd,
            tool="nhmmscan",
            inputs=[hmm_path, query_path],
            check=True,
            cancel=cancel,
        )
        if tblout.is_file():
            return parse_nhmmscan_tblout(tblout.read_text(encoding="utf-8", errors="replace"))
        # Rare: tblout missing but stdout still has a table.
        return parse_nhmmscan_tblout(result.stdout or "")


def collect_pyhmmer_domain_hits(
    hmm_path: Path,
    query_path: Path,
    *,
    report_e: float,
    domain_e: float | None = None,
    cpus: int = 1,
) -> list[tuple[float, float | None, str, str, int, int]]:
    """Return (score, evalue, query_id, hmm_name, seq_from, seq_to) for each domain hit.

    ``report_e`` sets sequence-level ``E`` / ``incE``. ``domain_e`` (default same
    as ``report_e``) sets ``domE`` / ``incdomE``. ``cpus`` is honored on modern
    pyhmmer; the legacy easy API always runs single-threaded.
    """
    import pyhmmer  # type: ignore[import-untyped]  # noqa: F401

    # Mirror BLAST's "Command:" block so Similarity detailed log shows the search.
    cpu_n = max(1, int(cpus))
    _log_scan("Command:")
    _log_scan(
        f"pyhmmer hmmscan -E {report_e:g} --cpu {cpu_n} "
        f"--hmm {hmm_path} --query {query_path}"
    )
    _log_scan("")
    _log_scan("—— program output ——")

    if hasattr(pyhmmer, "easy"):
        return _collect_with_legacy_pyhmmer(
            hmm_path,
            query_path,
            report_e=report_e,
            domain_e=domain_e,
            cpus=cpus,
        )
    return _collect_with_modern_pyhmmer(
        hmm_path,
        query_path,
        report_e=report_e,
        domain_e=domain_e,
        cpus=cpus,
    )


def collect_dna_hmm_domain_hits(
    hmm_path: Path,
    query_path: Path,
    *,
    report_e: float,
    domain_e: float | None = None,
    cpus: int = 1,
    cancel=None,
) -> list[tuple[float, float | None, str, str, int, int]]:
    """DNA HMM search: prefer system ``nhmmscan`` (scales with CPUs on large Dfam DBs).

    Falls back to pyhmmer ``hmmscan`` when ``nhmmscan`` is not installed. pyhmmer's
    DNA path does not usefully parallelize a single query against a pressed DB, so
    Dfam curated searches feel much slower without ``nhmmscan``.
    """
    hmm_path = hmm_path.expanduser().resolve()
    # nhmmscan needs a pressed DB; that is also where multi-CPU helps most (Dfam).
    if nhmmscan_available() and hmm_pressed_indexes_present(hmm_path):
        try:
            return collect_nhmmscan_domain_hits(
                hmm_path,
                query_path,
                report_e=report_e,
                cpus=cpus,
                cancel=cancel,
            )
        except Exception as exc:  # noqa: BLE001 — fall back to pyhmmer
            _log_scan(f"nhmmscan failed ({exc}); falling back to pyhmmer hmmscan")
    return collect_pyhmmer_domain_hits(
        hmm_path,
        query_path,
        report_e=report_e,
        domain_e=domain_e,
        cpus=cpus,
    )


def validate_hmm_search(hmm_path: Path, query_path: Path | None = None) -> None:
    """Run a tiny hmmscan to verify the pressed HMM database is readable."""
    hmm_path = hmm_path.expanduser().resolve()
    if query_path is None:
        with tempfile.TemporaryDirectory(prefix="repeatviewer_hmmcheck_") as tmp:
            tmp_path = Path(tmp)
            if hmm_uses_nucleotide_query(hmm_path):
                query_path = tmp_path / "query.fa"
                query_path.write_text(
                    ">test\nATGCGATCGATCGATCGATCGATCGATCGATCGATCG\n",
                    encoding="utf-8",
                )
            else:
                query_path = tmp_path / "query.faa"
                query_path.write_text(">test\nMKLLVVVFCLGVLA\n", encoding="utf-8")
            hits = collect_pyhmmer_domain_hits(hmm_path, query_path, report_e=10.0)
    else:
        hits = collect_pyhmmer_domain_hits(hmm_path, query_path, report_e=10.0)
    if not hits:
        raise RuntimeError("hmmscan returned no hits for the validation query.")
