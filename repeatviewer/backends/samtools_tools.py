"""Sort / index / convert SAM·BAM via pysam (macOS/Linux) or Win-NGS samtools (Windows).

Windows users install the unofficial Win-NGS UCRT64 bundle:
https://github.com/win-ngs/samtools-windows-bundle/releases/tag/v1.23.1-windows-ucrt64
Point Options → samtools executable at ``samtools.exe`` (same folder as its DLLs),
or put that folder on PATH.
"""

from __future__ import annotations

import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path

SAMTOOLS_WIN_BUNDLE_URL = (
    "https://github.com/win-ngs/samtools-windows-bundle/releases/tag/v1.23.1-windows-ucrt64"
)
SAMTOOLS_WIN_DOCS_URL = "https://win-ngs.github.io/"

SAM_SUFFIXES = {".sam"}
BAM_SUFFIXES = {".bam"}
CRAM_SUFFIXES = {".cram"}
ALIGNMENT_SUFFIXES = SAM_SUFFIXES | BAM_SUFFIXES | CRAM_SUFFIXES

OP_SORT = "sort"
OP_INDEX = "index"
OP_CONVERT = "convert"
SAMTOOLS_OPERATIONS = (OP_SORT, OP_INDEX, OP_CONVERT)

# samtools view can rewrite SAM ↔ BAM ↔ CRAM (no sorting).
CONVERT_FORMAT_BAM = "bam"
CONVERT_FORMAT_SAM = "sam"
CONVERT_FORMAT_CRAM = "cram"
CONVERT_FORMATS = (CONVERT_FORMAT_BAM, CONVERT_FORMAT_SAM, CONVERT_FORMAT_CRAM)
DEFAULT_CONVERT_FORMAT = CONVERT_FORMAT_BAM


def normalize_convert_format(raw: object) -> str:
    value = str(raw or DEFAULT_CONVERT_FORMAT).strip().lower()
    if value in CONVERT_FORMATS:
        return value
    return DEFAULT_CONVERT_FORMAT


class SamtoolsError(RuntimeError):
    """User-facing samtools / pysam failure."""


def is_windows() -> bool:
    return sys.platform == "win32" or platform.system() == "Windows"


def pysam_available() -> bool:
    try:
        import pysam  # noqa: F401
    except Exception:  # noqa: BLE001 - missing or broken build
        return False
    return True


def alignment_kind(path: Path) -> str | None:
    suffix = path.suffix.lower()
    if suffix == ".gz" and path.name.lower().endswith(".sam.gz"):
        return "sam"
    if suffix in SAM_SUFFIXES:
        return "sam"
    if suffix in BAM_SUFFIXES:
        return "bam"
    if suffix in CRAM_SUFFIXES:
        return "cram"
    return None


def is_alignment_path(path: Path) -> bool:
    return alignment_kind(path) is not None


# BAM/CRAM indexes stay beside the alignment but are hidden in Files.
_ALIGNMENT_INDEX_SUFFIXES = (".bai", ".crai", ".csi")


def is_alignment_index_file(path: Path) -> bool:
    """True for ``.bai`` / ``.crai`` / ``.csi`` (including ``file.bam.bai``)."""
    return Path(path).suffix.lower() in _ALIGNMENT_INDEX_SUFFIXES


def _alignment_index_candidate_pairs(
    old: Path, new: Path
) -> list[tuple[Path, Path]]:
    """Possible index name pairs when ``old`` is renamed/copied to ``new``."""
    return [
        (old.with_name(old.name + ".bai"), new.with_name(new.name + ".bai")),
        (old.with_suffix(".bai"), new.with_suffix(".bai")),
        (old.with_name(old.name + ".csi"), new.with_name(new.name + ".csi")),
        (old.with_suffix(".csi"), new.with_suffix(".csi")),
        (old.with_name(old.name + ".crai"), new.with_name(new.name + ".crai")),
        (old.with_suffix(".crai"), new.with_suffix(".crai")),
    ]


def alignment_index_companions(path: Path) -> list[Path]:
    """Existing BAM/CRAM index files that belong with this alignment."""
    path = Path(path).expanduser()
    try:
        path = path.resolve()
    except OSError:
        pass
    kind = alignment_kind(path)
    if kind not in ("bam", "cram"):
        return []
    found: list[Path] = []
    seen: set[str] = set()
    for candidate, _dst in _alignment_index_candidate_pairs(path, path):
        try:
            key = str(candidate.resolve())
        except OSError:
            key = str(candidate)
        if key in seen:
            continue
        seen.add(key)
        if candidate.is_file():
            found.append(candidate)
    return found


def alignment_index_renames(old: Path, new: Path) -> list[tuple[Path, Path]]:
    """Keep ``.bai`` / ``.crai`` / ``.csi`` named beside the new BAM/CRAM.

    ``reads.sorted.bam.bai`` follows ``reads.sorted.bam`` → ``other.bam`` as
    ``other.bam.bai``. Same-format no-ops are omitted.
    """
    old, new = Path(old), Path(new)
    if alignment_kind(old) not in ("bam", "cram"):
        return []
    pairs: list[tuple[Path, Path]] = []
    seen_src: set[str] = set()
    for src, dst in _alignment_index_candidate_pairs(old, new):
        if not src.is_file():
            continue
        try:
            src_key = str(src.resolve())
            same = src.resolve() == dst.resolve()
        except OSError:
            src_key = str(src)
            same = src == dst
        if src_key in seen_src or same:
            continue
        seen_src.add(src_key)
        pairs.append((src, dst))
    return pairs


def default_sorted_bam_path(source: Path) -> Path:
    stem = source.name
    lower = stem.lower()
    for suffix in (".sam.gz", ".sam", ".bam", ".cram"):
        if lower.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    return source.with_name(f"{stem}.sorted.bam")


def _alignment_stem(source: Path) -> str:
    stem = source.name
    lower = stem.lower()
    for suffix in (".sam.gz", ".sam", ".bam", ".cram"):
        if lower.endswith(suffix):
            return stem[: -len(suffix)]
    return source.stem


def default_converted_path(
    source: Path, output_format: str = DEFAULT_CONVERT_FORMAT
) -> Path:
    """Beside-input path for Convert (``stem.bam`` / ``.sam`` / ``.cram``).

    Same-format inputs are skipped by :func:`convert_alignment` (no rewrite).
    """
    fmt = normalize_convert_format(output_format)
    stem = _alignment_stem(source)
    return source.with_name(f"{stem}.{fmt}")


def default_converted_bam_path(source: Path) -> Path:
    return default_converted_path(source, CONVERT_FORMAT_BAM)


def resolve_samtools_executable(
    *,
    explicit_path: str | Path | None = None,
    use_system: bool = True,
) -> str | None:
    """Return a samtools executable path, or None if not found.

    On Windows this is the Win-NGS / PATH ``samtools.exe``. On macOS/Linux the
    CLI is optional (pysam is preferred); PATH is still accepted as a fallback.
    """
    raw = str(explicit_path or "").strip()
    if raw:
        path = Path(raw).expanduser()
        if path.is_file():
            return str(path.resolve())
    if not use_system:
        return None
    for name in ("samtools.exe", "samtools") if is_windows() else ("samtools",):
        found = shutil.which(name)
        if found:
            return found
    return None


def backend_label(*, exe: str | None = None) -> str:
    if is_windows():
        return "Win-NGS samtools" if exe else "Win-NGS samtools (not found)"
    if pysam_available():
        return "pysam"
    if exe:
        return "samtools (PATH)"
    return "pysam / samtools (not found)"


def missing_backend_message(*, exe: str | None = None) -> str:
    if is_windows():
        return (
            "samtools was not found.\n\n"
            "On Windows, RepeatViewer uses the Win-NGS samtools bundle:\n"
            f"  {SAMTOOLS_WIN_BUNDLE_URL}\n\n"
            "Download the ZIP, extract it, and either:\n"
            "• set Options → samtools executable to samtools.exe in that folder, or\n"
            "• add that folder to PATH.\n\n"
            f"More Windows NGS tools: {SAMTOOLS_WIN_DOCS_URL}"
        )
    if not pysam_available() and not exe:
        return (
            "Neither pysam nor samtools is available.\n\n"
            "Install pysam (recommended on macOS / Linux):\n"
            "  pip install pysam\n\n"
            "Or install samtools and ensure it is on PATH."
        )
    return "samtools backend is not available."


def _run_cli(exe: str, args: list[str], *, progress=None, label: str = "") -> None:
    cmd = [exe, *args]
    if progress is not None and label:
        progress(0, label)
    # On Windows the Win-NGS DLLs sit next to samtools.exe.
    env = os.environ.copy()
    exe_dir = str(Path(exe).resolve().parent)
    path_key = "PATH"
    env[path_key] = exe_dir + os.pathsep + env.get(path_key, "")
    try:
        completed = subprocess.run(
            cmd,
            check=False,
            capture_output=True,
            text=True,
            env=env,
        )
    except OSError as exc:
        raise SamtoolsError(f"Could not run samtools: {exc}") from exc
    if completed.returncode != 0:
        err = (completed.stderr or completed.stdout or "").strip()
        raise SamtoolsError(err or f"samtools exited with code {completed.returncode}.")


def _require_backend(*, explicit_path: str | None, use_system: bool) -> tuple[str, str | None]:
    """Return (\"pysam\"|\"cli\", exe_or_none)."""
    exe = resolve_samtools_executable(
        explicit_path=explicit_path, use_system=use_system
    )
    if not is_windows() and pysam_available():
        return "pysam", exe
    if exe:
        return "cli", exe
    raise SamtoolsError(missing_backend_message(exe=exe))


def sort_alignment(
    source: Path,
    output: Path | None = None,
    *,
    threads: int = 4,
    index_after: bool = True,
    explicit_path: str | None = None,
    use_system: bool = True,
    progress=None,
) -> list[Path]:
    """Coordinate-sort to BAM (and optionally index)."""
    source = Path(source).expanduser().resolve()
    if not source.is_file():
        raise SamtoolsError(f"Input not found: {source}")
    if alignment_kind(source) is None:
        raise SamtoolsError(f"{source.name} is not a SAM/BAM/CRAM file.")
    out = Path(output).expanduser().resolve() if output else default_sorted_bam_path(source)
    out.parent.mkdir(parents=True, exist_ok=True)
    backend, exe = _require_backend(explicit_path=explicit_path, use_system=use_system)
    threads = max(1, int(threads))
    if progress is not None:
        progress(10, f"Sorting {source.name}…")
    if backend == "pysam":
        import pysam

        pysam.sort("-@", str(threads), "-o", str(out), str(source))
    else:
        assert exe is not None
        _run_cli(
            exe,
            ["sort", "-@", str(threads), "-o", str(out), str(source)],
            progress=progress,
            label=f"Sorting {source.name}…",
        )
    created = [out]
    if index_after:
        if progress is not None:
            progress(70, f"Indexing {out.name}…")
        created.extend(
            index_alignment(
                out,
                threads=threads,
                explicit_path=explicit_path,
                use_system=use_system,
                progress=None,
            )
        )
    if progress is not None:
        progress(100, f"Sorted → {out.name}")
    return created


def index_alignment(
    source: Path,
    *,
    threads: int = 4,
    explicit_path: str | None = None,
    use_system: bool = True,
    progress=None,
) -> list[Path]:
    """Build a ``.bai`` (or CSI) next to a coordinate-sorted BAM/CRAM."""
    source = Path(source).expanduser().resolve()
    if not source.is_file():
        raise SamtoolsError(f"Input not found: {source}")
    kind = alignment_kind(source)
    if kind not in {"bam", "cram"}:
        raise SamtoolsError(
            f"{source.name} must be a sorted BAM or CRAM to index. "
            "Use Sort first for SAM files."
        )
    backend, exe = _require_backend(explicit_path=explicit_path, use_system=use_system)
    threads = max(1, int(threads))
    if progress is not None:
        progress(20, f"Indexing {source.name}…")
    if backend == "pysam":
        import pysam

        # Newer pysam accepts ``-@``; fall back to a single-thread call.
        try:
            pysam.index("-@", str(threads), str(source))
        except TypeError:
            pysam.index(str(source))
        except Exception:
            pysam.index(str(source))
    else:
        assert exe is not None
        _run_cli(
            exe,
            ["index", "-@", str(threads), str(source)],
            progress=progress,
            label=f"Indexing {source.name}…",
        )
    # Prefer .bai; CSI is also valid.
    candidates = [
        source.with_name(source.name + ".bai"),
        source.with_suffix(".bai"),
        source.with_name(source.name + ".csi"),
        source.with_name(source.name + ".crai"),
    ]
    created = [p for p in candidates if p.is_file()]
    if not created:
        raise SamtoolsError(
            f"Index finished but no .bai/.csi was found beside {source.name}."
        )
    if progress is not None:
        progress(100, f"Indexed → {created[0].name}")
    return created


def _view_format_args(output_format: str) -> list[str]:
    """Flags for ``samtools view`` / pysam to select SAM, BAM, or CRAM output."""
    fmt = normalize_convert_format(output_format)
    if fmt == CONVERT_FORMAT_SAM:
        return ["-h"]  # SAM text with header
    if fmt == CONVERT_FORMAT_CRAM:
        return ["-O", "cram"]
    return ["-b"]


def convert_alignment(
    source: Path,
    output: Path | None = None,
    *,
    output_format: str = DEFAULT_CONVERT_FORMAT,
    threads: int = 4,
    explicit_path: str | None = None,
    use_system: bool = True,
    progress=None,
) -> list[Path]:
    """Rewrite alignment format without sorting (``samtools view`` / pysam).

    Input is auto-detected (SAM/BAM/CRAM). Output is BAM (default), SAM, or CRAM.
    Same-format inputs are skipped (returns ``[]``) — no pointless rewrite.
    CRAM writing may require a matching reference in the file header / ``@SQ``.
    """
    source = Path(source).expanduser().resolve()
    if not source.is_file():
        raise SamtoolsError(f"Input not found: {source}")
    kind = alignment_kind(source)
    if kind is None:
        raise SamtoolsError(f"{source.name} is not a SAM/BAM/CRAM file.")
    fmt = normalize_convert_format(output_format)
    if kind == fmt:
        if progress is not None:
            progress(100, f"Skipped {source.name} (already {fmt.upper()})")
        return []
    out = (
        Path(output).expanduser().resolve()
        if output
        else default_converted_path(source, fmt)
    )
    if out.resolve() == source.resolve():
        raise SamtoolsError(
            f"Convert would overwrite the input file ({source.name})."
        )
    out.parent.mkdir(parents=True, exist_ok=True)
    backend, exe = _require_backend(explicit_path=explicit_path, use_system=use_system)
    threads = max(1, int(threads))
    fmt_args = _view_format_args(fmt)
    label_fmt = fmt.upper()
    if progress is not None:
        progress(10, f"Converting {source.name} → {label_fmt}…")
    try:
        if backend == "pysam":
            _pysam_view_to_file(
                source,
                out,
                fmt_args=fmt_args,
                threads=threads,
            )
        else:
            assert exe is not None
            _run_cli(
                exe,
                ["view", "-@", str(threads), *fmt_args, "-o", str(out), str(source)],
                progress=progress,
                label=f"Converting {source.name}…",
            )
    except Exception as exc:  # noqa: BLE001 - surface as SamtoolsError
        if isinstance(exc, SamtoolsError):
            raise
        message = str(exc).strip() or repr(exc)
        if fmt == CONVERT_FORMAT_CRAM and (
            "reference" in message.lower() or "ref" in message.lower()
        ):
            raise SamtoolsError(
                f"CRAM conversion failed for {source.name}. "
                "CRAM usually needs a matching reference sequence in the "
                f"alignment header. Details: {message}"
            ) from exc
        raise SamtoolsError(
            f"Convert {source.name} → {label_fmt} failed: {message}"
        ) from exc
    if not out.is_file() or out.stat().st_size <= 0:
        raise SamtoolsError(f"Convert finished but {out.name} was not created.")
    if progress is not None:
        progress(100, f"Converted → {out.name}")
    return [out]


def _pysam_view_to_file(
    source: Path,
    out: Path,
    *,
    fmt_args: list[str],
    threads: int,
) -> None:
    """Run ``pysam.view … -o out``.

    pysam's default ``catch_stdout=True`` can swallow binary BAM/CRAM and leave
    ``-o`` unwritten (call appears to succeed, no file on disk). Prefer
    ``catch_stdout=False``; fall back to writing returned bytes if needed.
    """
    import pysam

    args = ["-@", str(threads), *fmt_args, "-o", str(out), str(source)]
    try:
        pysam.view(*args, catch_stdout=False)
    except TypeError:
        # Older pysam: no catch_stdout kwarg — may return the file contents.
        data = pysam.view(*args)
        if out.is_file() and out.stat().st_size > 0:
            return
        if isinstance(data, (bytes, bytearray)) and data:
            out.write_bytes(bytes(data))
            return
        if isinstance(data, str) and data:
            # Text SAM from ``-h``; binary BAM must not go through str.
            if "-b" in fmt_args or any(
                a == "cram" or str(a).startswith("cram") for a in fmt_args
            ):
                raise SamtoolsError(
                    "pysam.view returned text instead of writing a binary "
                    f"file ({out.name}). Update pysam or install samtools."
                )
            out.write_text(data, encoding="utf-8")
            return
        raise SamtoolsError(
            f"pysam.view did not write {out.name}."
        ) from None


def convert_to_bam(
    source: Path,
    output: Path | None = None,
    *,
    threads: int = 4,
    explicit_path: str | None = None,
    use_system: bool = True,
    progress=None,
) -> list[Path]:
    """Write BAM without coordinate sorting (``samtools view -b``)."""
    return convert_alignment(
        source,
        output,
        output_format=CONVERT_FORMAT_BAM,
        threads=threads,
        explicit_path=explicit_path,
        use_system=use_system,
        progress=progress,
    )


def run_samtools_operation(
    operation: str,
    sources: list[Path],
    *,
    threads: int = 4,
    index_after_sort: bool = True,
    convert_format: str = DEFAULT_CONVERT_FORMAT,
    explicit_path: str | None = None,
    use_system: bool = True,
    progress=None,
) -> dict:
    """Run sort, index, or convert on one or more alignment files."""
    op = str(operation or "").strip().lower()
    if op not in SAMTOOLS_OPERATIONS:
        raise SamtoolsError(f"Unknown samtools operation: {operation!r}")
    files = [Path(p) for p in sources if Path(p).is_file()]
    if not files:
        raise SamtoolsError("No SAM/BAM/CRAM input files.")
    created: list[Path] = []
    n = len(files)
    for i, path in enumerate(files):
        if progress is not None:
            base = int(100 * i / max(1, n))

            def _p(value: int, label: str, _base=base, _i=i) -> None:
                scaled = _base + int(value * (1 / max(1, n)))
                progress(min(99, scaled), label)

        else:
            _p = None  # type: ignore[assignment]
        if op == OP_SORT:
            created.extend(
                sort_alignment(
                    path,
                    threads=threads,
                    index_after=index_after_sort,
                    explicit_path=explicit_path,
                    use_system=use_system,
                    progress=_p,
                )
            )
        elif op == OP_INDEX:
            created.extend(
                index_alignment(
                    path,
                    threads=threads,
                    explicit_path=explicit_path,
                    use_system=use_system,
                    progress=_p,
                )
            )
        else:
            created.extend(
                convert_alignment(
                    path,
                    output_format=convert_format,
                    threads=threads,
                    explicit_path=explicit_path,
                    use_system=use_system,
                    progress=_p,
                )
            )
    if op == OP_CONVERT and not created:
        fmt = normalize_convert_format(convert_format).upper()
        raise SamtoolsError(
            f"All inputs are already {fmt}; nothing to convert."
        )
    if progress is not None:
        progress(100, f"Done ({len(created)} file(s)).")
    return {
        "operation": op,
        "files": [str(p) for p in created],
        "inputs": [str(p) for p in files],
    }
