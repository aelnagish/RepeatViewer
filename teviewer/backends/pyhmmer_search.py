"""pyhmmer HMM scan helpers (supports pyhmmer 0.10 easy API and 0.12+ plan7/hmmer API)."""

from __future__ import annotations

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


def _scan_kwargs(report_e: float) -> dict[str, float | str]:
    return {
        "E": report_e,
        "domE": report_e,
        "incE": report_e,
        "incdomE": report_e,
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


def _collect_with_modern_pyhmmer(
    hmm_path: Path,
    query_path: Path,
    *,
    report_e: float,
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
        return []

    raw_hits: list[tuple[float, float | None, str, str, int, int]] = []
    with HMMFile(str(hmm_path)) as hmm_file:
        profiles = hmm_file.optimized_profiles()
        for top_hits in hmmscan(queries, profiles, cpus=1, **_scan_kwargs(report_e)):
            qname = _query_name_from_top_hits(top_hits)
            for hit in top_hits:
                _append_domain_hits(raw_hits, qname=qname, hit=hit)
    return raw_hits


def _collect_with_legacy_pyhmmer(
    hmm_path: Path,
    query_path: Path,
    *,
    report_e: float,
) -> list[tuple[float, float | None, str, str, int, int]]:
    import pyhmmer  # type: ignore[import-untyped]

    raw_hits: list[tuple[float, float | None, str, str, int, int]] = []
    with pyhmmer.easy.HMMFile(str(hmm_path)) as hmm_db:
        with pyhmmer.easy.SeqFile(str(query_path)) as seq_file:
            scan_kwargs = _scan_kwargs(report_e)
            for top_hits in pyhmmer.hmmscan(hmm_db, seq_file, **scan_kwargs):
                qname = _query_name_from_top_hits(top_hits)
                for hit in top_hits:
                    _append_domain_hits(raw_hits, qname=qname, hit=hit)
    return raw_hits


def collect_pyhmmer_domain_hits(
    hmm_path: Path,
    query_path: Path,
    *,
    report_e: float,
) -> list[tuple[float, float | None, str, str, int, int]]:
    """Return (score, evalue, query_id, hmm_name, seq_from, seq_to) for each domain hit."""
    import pyhmmer  # type: ignore[import-untyped]  # noqa: F401

    if hasattr(pyhmmer, "easy"):
        return _collect_with_legacy_pyhmmer(hmm_path, query_path, report_e=report_e)
    return _collect_with_modern_pyhmmer(hmm_path, query_path, report_e=report_e)


def validate_hmm_search(hmm_path: Path, query_path: Path | None = None) -> None:
    """Run a tiny hmmscan to verify the pressed HMM database is readable."""
    hmm_path = hmm_path.expanduser().resolve()
    if query_path is None:
        with tempfile.TemporaryDirectory(prefix="teviewer_hmmcheck_") as tmp:
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
