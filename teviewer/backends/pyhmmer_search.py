"""pyhmmer HMM scan helpers (supports pyhmmer 0.10 easy API and 0.12+ plan7/hmmer API)."""

from __future__ import annotations

from pathlib import Path


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


def _read_digital_queries(query_path: Path) -> list[object]:
    from pyhmmer.easel import SequenceFile  # type: ignore[import-untyped]

    with SequenceFile(str(query_path), digital=True) as seq_file:
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
    from pyhmmer.hmmer import hmmscan  # type: ignore[import-untyped]
    from pyhmmer.plan7 import HMMFile  # type: ignore[import-untyped]

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


def validate_hmm_search(hmm_path: Path, query_faa: Path) -> None:
    """Run a tiny hmmscan to verify the pressed HMM database is readable."""
    hits = collect_pyhmmer_domain_hits(hmm_path, query_faa, report_e=10.0)
    if not hits:
        raise RuntimeError("hmmscan returned no hits for the validation query.")
