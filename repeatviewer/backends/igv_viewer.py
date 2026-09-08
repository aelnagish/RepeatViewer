"""Local igv.js preview for BAM/SAM alignment files."""

from __future__ import annotations

import json
import logging
import os
import re
import secrets
import shutil
import socket
import threading
import urllib.parse
from collections.abc import Callable
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

import pysam

from repeatviewer.backends.gff_features import find_companion_gff

logger = logging.getLogger(__name__)

_ALIGNMENT_SUFFIXES = frozenset({".bam", ".sam"})
_GFF_SUFFIXES = frozenset({".gff", ".gff3", ".gtf"})
_BED_SUFFIXES = frozenset({".bed"})
_TRACK_SUFFIXES = _GFF_SUFFIXES | _BED_SUFFIXES | _ALIGNMENT_SUFFIXES
_MAX_ANNOTATION_TRACKS = 12
_SESSION_ID_RE = re.compile(r"^[a-f0-9]{16}$")

_extra_annotations_provider: Callable[[Path], list[Path]] | None = None

_sessions_root: Path | None = None
_session_dirs: dict[str, Path] = {}
_session_upload_counters: dict[str, int] = {}
_server: ThreadingHTTPServer | None = None
_server_thread: threading.Thread | None = None
_server_port: int | None = None
_server_lock = threading.Lock()


class IgvPreviewError(Exception):
    """User-facing igv.js preview setup failure."""


def is_igv_view_url(url: str | None) -> bool:
    if not url:
        return False
    return parse_igv_session_id(url) is not None


def parse_igv_session_id(url: str) -> str | None:
    try:
        parsed = urllib.parse.urlparse(url)
    except ValueError:
        return None
    if parsed.hostname not in {"127.0.0.1", "localhost"}:
        return None
    parts = [part for part in parsed.path.split("/") if part]
    if len(parts) < 2:
        return None
    session_id, page = parts[0], parts[1]
    if page != "viewer.html" or not _SESSION_ID_RE.fullmatch(session_id):
        return None
    return session_id


def is_alignment_preview_file(path: Path) -> bool:
    return path.suffix.lower() in _ALIGNMENT_SUFFIXES


def set_igv_extra_annotations_provider(
    provider: Callable[[Path], list[Path]] | None,
) -> None:
    """Optional hook: extra GFF/GTF paths to load (e.g. from Files panel selection)."""
    global _extra_annotations_provider
    _extra_annotations_provider = provider


def _query_stem_from_alignment(alignment: Path) -> str:
    stem = alignment.stem
    match = re.match(r"^(.+?)\.minimap2(?:\.|$)", stem)
    if match:
        return match.group(1)
    for suffix in (".sorted", ".aligned", ".mapped"):
        if stem.endswith(suffix):
            return stem[: -len(suffix)]
    return stem


def discover_igv_annotations(alignment: Path) -> list[Path]:
    """Find GFF/GTF files to show alongside an alignment in igv.js."""
    alignment = alignment.expanduser().resolve()
    folder = alignment.parent
    query_stem = _query_stem_from_alignment(alignment)
    found: list[Path] = []
    seen: set[str] = set()

    def add(path: Path | None) -> None:
        if path is None or not path.is_file():
            return
        resolved = path.expanduser().resolve()
        if resolved.suffix.lower() not in _GFF_SUFFIXES:
            return
        key = str(resolved)
        if key in seen:
            return
        seen.add(key)
        found.append(resolved)

    for suffix in (".fa", ".fasta", ".fas", ".fna"):
        add(find_companion_gff(folder / f"{query_stem}{suffix}"))

    prioritized: list[Path] = []
    other: list[Path] = []
    for suffix in (".gff3", ".gff", ".gtf"):
        for path in sorted(folder.glob(f"*{suffix}")):
            if not path.is_file():
                continue
            name_lower = path.name.lower()
            if query_stem.lower() in name_lower or query_stem.lower() in path.stem.lower():
                prioritized.append(path)
            else:
                other.append(path)
    for path in prioritized + other:
        add(path)
        if len(found) >= _MAX_ANNOTATION_TRACKS:
            break

    if _extra_annotations_provider is not None:
        try:
            for path in _extra_annotations_provider(alignment):
                add(Path(path))
                if len(found) >= _MAX_ANNOTATION_TRACKS:
                    break
        except Exception:
            logger.debug("igv extra annotation provider failed", exc_info=True)

    return found[:_MAX_ANNOTATION_TRACKS]


def _igv_annotation_format(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix == ".gff3":
        return "gff3"
    if suffix == ".gtf":
        return "gtf"
    return "gff"


def _register_session(session_id: str, session_dir: Path) -> None:
    _session_dirs[session_id] = session_dir
    _session_upload_counters.setdefault(session_id, 1000)


def _session_dir(session_id: str) -> Path:
    session_dir = _session_dirs.get(session_id)
    if session_dir is None or not session_dir.is_dir():
        raise IgvPreviewError("IGV session expired. Close the preview tab and reopen the alignment.")
    return session_dir


def _next_upload_name(session_id: str, src: Path) -> str:
    counter = _session_upload_counters.get(session_id, 1000) + 1
    _session_upload_counters[session_id] = counter
    safe_stem = re.sub(r"[^\w.\-+#]+", "_", src.stem) or "track"
    return f"upload_{counter}_{safe_stem}{src.suffix.lower()}"


def _session_base_url(session_id: str) -> str:
    port = _ensure_server()
    return f"http://127.0.0.1:{port}/{session_id}"


def _session_file_url(session_id: str, filename: str) -> str:
    return f"{_session_base_url(session_id)}/{filename}"


def _materialize_file(src: Path, dest: Path) -> None:
    """Copy a host file into the igv session (avoids broken symlinks in HTTP serving)."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    src = src.expanduser().resolve()
    if dest.exists() or dest.is_symlink():
        dest.unlink()
    if src == dest.resolve():
        return
    try:
        os.link(src, dest)
        return
    except OSError:
        pass
    shutil.copy2(src, dest)


def _verify_session_file(session_dir: Path, filename: str) -> Path:
    path = session_dir / filename
    if not path.is_file() or path.stat().st_size == 0:
        raise IgvPreviewError(f"Track file is missing or empty: {filename}")
    return path


def _track_config_for_session_file(
    session_dir: Path,
    url_name: str,
    display_name: str,
    *,
    base_url: str,
) -> dict[str, str | int]:
    path = _verify_session_file(session_dir, url_name)
    suffix = path.suffix.lower()
    if suffix in _GFF_SUFFIXES:
        return {
            "name": display_name,
            "type": "annotation",
            "format": _igv_annotation_format(path),
            "url": f"{base_url}/{url_name}",
            "displayMode": "EXPANDED",
            "height": 80,
        }
    if suffix in _BED_SUFFIXES:
        return {
            "name": display_name,
            "type": "annotation",
            "format": "bed",
            "url": f"{base_url}/{url_name}",
            "displayMode": "EXPANDED",
            "height": 50,
        }
    if suffix == ".bam":
        ensure_bam_index(path)
        bai_name = f"{url_name}.bai"
        _verify_session_file(session_dir, bai_name)
        return {
            "name": display_name,
            "type": "alignment",
            "format": "bam",
            "url": f"{base_url}/{url_name}",
            "indexURL": f"{base_url}/{bai_name}",
        }
    if suffix == ".sam":
        sorted_name = f"{url_name}.sorted.bam"
        sorted_bam = session_dir / sorted_name
        sam_to_sorted_bam(path, sorted_bam)
        _verify_session_file(session_dir, sorted_name)
        ensure_bam_index(sorted_bam)
        bai_name = f"{sorted_name}.bai"
        _verify_session_file(session_dir, bai_name)
        return {
            "name": display_name,
            "type": "alignment",
            "format": "bam",
            "url": f"{base_url}/{sorted_name}",
            "indexURL": f"{base_url}/{bai_name}",
        }
    raise IgvPreviewError(
        f"Unsupported track type {suffix or '(none)'}. "
        "Use GFF, GTF, BED, BAM, or SAM."
    )


def add_track_to_igv_session(session_id: str, source: Path) -> dict[str, str | int]:
    """Register a host file in the igv session and return a loadTrack config."""
    session_dir = _session_dir(session_id)
    base_url = _session_base_url(session_id)
    src = source.expanduser().resolve()
    if not src.is_file():
        raise IgvPreviewError(f"File not found:\n{src}")
    if src.suffix.lower() not in _TRACK_SUFFIXES:
        raise IgvPreviewError(
            f"Unsupported file type {src.suffix!r}. "
            "Choose GFF, GTF, BED, BAM, or SAM."
        )
    url_name = _next_upload_name(session_id, src)
    _materialize_file(src, session_dir / url_name)
    return _track_config_for_session_file(
        session_dir,
        url_name,
        src.name,
        base_url=base_url,
    )


def add_uploaded_track_to_session(
    session_id: str,
    filename: str,
    data: bytes,
) -> dict[str, str | int]:
    """Save an uploaded browser file into the session and return loadTrack config."""
    session_dir = _session_dir(session_id)
    name = Path(filename).name
    if not name or name in {".", ".."}:
        raise IgvPreviewError("Invalid upload file name.")
    src = Path(name)
    if src.suffix.lower() not in _TRACK_SUFFIXES:
        raise IgvPreviewError(
            f"Unsupported file type {src.suffix!r}. "
            "Use GFF, GTF, BED, BAM, or SAM."
        )
    url_name = _next_upload_name(session_id, src)
    dest = session_dir / url_name
    dest.write_bytes(data)
    return _track_config_for_session_file(
        session_dir,
        url_name,
        name,
        base_url=_session_base_url(session_id),
    )


def _prepare_annotation_tracks(
    session_dir: Path,
    paths: list[Path],
    *,
    base_url: str,
) -> list[dict[str, str | int]]:
    tracks: list[dict[str, str | int]] = []
    for index, src in enumerate(paths):
        ext = src.suffix.lower()
        session_name = f"annotation_{index}{ext}"
        try:
            _materialize_file(src, session_dir / session_name)
            tracks.append(
                _track_config_for_session_file(
                    session_dir,
                    session_name,
                    src.name,
                    base_url=base_url,
                )
            )
        except (IgvPreviewError, OSError) as exc:
            logger.warning("Skipping annotation track %s: %s", src, exc)
    return tracks


def igv_assets_dir() -> Path:
    return Path(__file__).resolve().parent.parent / "assets" / "igv"


def _pysam_error(action: str, exc: BaseException) -> RuntimeError:
    return RuntimeError(f"pysam {action} failed: {exc}")


def ensure_fasta_index(fasta: Path) -> Path:
    fasta = fasta.expanduser().resolve()
    fai = Path(f"{fasta}.fai")
    if fai.is_file():
        return fai
    try:
        pysam.faidx(str(fasta))
    except Exception as exc:  # noqa: BLE001
        raise _pysam_error("faidx", exc) from exc
    if not fai.is_file():
        raise RuntimeError(f"pysam faidx did not create {fai}")
    return fai


def ensure_bam_index(bam: Path) -> Path:
    bam = bam.expanduser().resolve()
    candidates = (Path(f"{bam}.bai"), bam.with_suffix(bam.suffix + ".bai"))
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    try:
        pysam.index(str(bam))
    except Exception as exc:  # noqa: BLE001
        raise _pysam_error("index", exc) from exc
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise RuntimeError(f"pysam index did not create an index for {bam}")


def sam_to_sorted_bam(sam: Path, sorted_bam: Path) -> None:
    """Convert SAM to coordinate-sorted BAM (required for igv.js)."""
    sam = sam.expanduser().resolve()
    sorted_bam = sorted_bam.expanduser().resolve()
    sorted_bam.parent.mkdir(parents=True, exist_ok=True)
    if sorted_bam.exists():
        sorted_bam.unlink()
    try:
        pysam.sort("-o", str(sorted_bam), str(sam))
    except Exception as exc:  # noqa: BLE001
        raise _pysam_error("sort", exc) from exc
    if not sorted_bam.is_file() or sorted_bam.stat().st_size == 0:
        raise RuntimeError(f"pysam sort did not create {sorted_bam}")


def _sessions_directory() -> Path:
    global _sessions_root
    if _sessions_root is None:
        env_root = os.environ.get("REPEATVIEWER_IGV_SESSIONS", "").strip()
        if env_root:
            root = Path(env_root).expanduser()
        else:
            import tempfile

            root = Path(tempfile.gettempdir()) / "repeatviewer-igv-sessions"
        root.mkdir(parents=True, exist_ok=True)
        _sessions_root = root
    return _sessions_root


def _find_free_port() -> int:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
        sock.bind(("127.0.0.1", 0))
        return int(sock.getsockname()[1])


def _ensure_server() -> int:
    global _server, _server_thread, _server_port
    with _server_lock:
        if _server is not None and _server_port is not None:
            return _server_port

        root = str(_sessions_directory())

        class _Handler(SimpleHTTPRequestHandler):
            def __init__(self, *args, **kwargs) -> None:
                super().__init__(*args, directory=root, **kwargs)

            def log_message(self, format: str, *args) -> None:  # noqa: A003
                logger.debug("igv http: %s", format % args)

            def end_headers(self) -> None:
                self.send_header("Access-Control-Allow-Origin", "*")
                self.send_header("Accept-Ranges", "bytes")
                super().end_headers()

            def do_POST(self) -> None:  # noqa: N802
                parsed = urllib.parse.urlparse(self.path)
                parts = [part for part in parsed.path.split("/") if part]
                if len(parts) != 2 or parts[1] != "upload":
                    self.send_error(404, "Not found")
                    return
                session_id = parts[0]
                if not _SESSION_ID_RE.fullmatch(session_id):
                    self.send_error(400, "Invalid session")
                    return
                query = urllib.parse.parse_qs(parsed.query)
                filename = query.get("name", ["upload.track"])[0]
                try:
                    length = int(self.headers.get("Content-Length", "0"))
                except ValueError:
                    self.send_error(400, "Missing Content-Length")
                    return
                if length <= 0:
                    self.send_error(400, "Empty upload")
                    return
                data = self.rfile.read(length)
                try:
                    track = add_uploaded_track_to_session(session_id, filename, data)
                except IgvPreviewError as exc:
                    payload = str(exc).encode("utf-8")
                    self.send_response(400)
                    self.send_header("Content-Type", "text/plain; charset=utf-8")
                    self.send_header("Content-Length", str(len(payload)))
                    self.end_headers()
                    self.wfile.write(payload)
                    return
                except Exception as exc:  # noqa: BLE001
                    logger.exception("igv upload failed")
                    payload = str(exc).encode("utf-8")
                    self.send_response(500)
                    self.send_header("Content-Type", "text/plain; charset=utf-8")
                    self.send_header("Content-Length", str(len(payload)))
                    self.end_headers()
                    self.wfile.write(payload)
                    return

                body = json.dumps({"track": track}).encode("utf-8")
                self.send_response(200)
                self.send_header("Content-Type", "application/json")
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)

        port = _find_free_port()
        _server = ThreadingHTTPServer(("127.0.0.1", port), _Handler)
        _server.daemon_threads = True
        _server_port = port
        _server_thread = threading.Thread(
            target=_server.serve_forever,
            name="repeatviewer-igv-http",
            daemon=True,
        )
        _server_thread.start()
        logger.info("RepeatViewer igv.js HTTP server on 127.0.0.1:%s", port)
        return port


def _prepare_alignment_track(
    session_dir: Path,
    alignment: Path,
    *,
    base_url: str,
) -> dict[str, str]:
    suffix = alignment.suffix.lower()
    if suffix == ".bam":
        bam_name = "reads.bam"
        _materialize_file(alignment.resolve(), session_dir / bam_name)
        ensure_bam_index(session_dir / bam_name)
        bai_name = f"{bam_name}.bai"
        _verify_session_file(session_dir, bai_name)
        return {
            "name": alignment.name,
            "type": "alignment",
            "format": "bam",
            "url": f"{base_url}/{bam_name}",
            "indexURL": f"{base_url}/{bai_name}",
        }

    sam_name = "reads.sam"
    _materialize_file(alignment.resolve(), session_dir / sam_name)
    sorted_name = "reads.sorted.bam"
    sorted_bam = session_dir / sorted_name
    sam_to_sorted_bam(session_dir / sam_name, sorted_bam)
    _verify_session_file(session_dir, sorted_name)
    ensure_bam_index(sorted_bam)
    bai_name = f"{sorted_name}.bai"
    _verify_session_file(session_dir, bai_name)
    return {
        "name": alignment.name,
        "type": "alignment",
        "format": "bam",
        "url": f"{base_url}/{sorted_name}",
        "indexURL": f"{base_url}/{bai_name}",
    }


def build_igv_view_url(
    alignment: Path,
    reference: Path | None,
    *,
    annotation_paths: list[Path] | None = None,
) -> str:
    """Create a local igv.js session and return its http URL."""
    alignment = alignment.expanduser().resolve()
    if not alignment.is_file():
        raise IgvPreviewError(f"Alignment file not found:\n{alignment}")

    if reference is None or not reference.is_file():
        raise IgvPreviewError(
            "Reference genome FASTA is required.\n"
            "Set it with File → Define genome, then reopen this alignment."
        )
    reference = reference.expanduser().resolve()

    if annotation_paths is None:
        annotations = discover_igv_annotations(alignment)
    else:
        annotations = []
        seen: set[str] = set()
        for raw in annotation_paths:
            path = Path(raw).expanduser().resolve()
            if not path.is_file() or path.suffix.lower() not in _GFF_SUFFIXES:
                continue
            key = str(path)
            if key in seen:
                continue
            seen.add(key)
            annotations.append(path)
        if not annotations:
            annotations = discover_igv_annotations(alignment)

    assets = igv_assets_dir()
    igv_js = assets / "igv.min.js"
    viewer_html = assets / "viewer.html"
    if not igv_js.is_file() or not viewer_html.is_file():
        raise IgvPreviewError(
            "Bundled igv.js assets are missing from the RepeatViewer install "
            f"({assets})."
        )

    session_id = secrets.token_hex(8)
    session_dir = _sessions_directory() / session_id
    session_dir.mkdir(parents=True, exist_ok=True)
    _register_session(session_id, session_dir)

    port = _ensure_server()
    base_url = f"http://127.0.0.1:{port}/{session_id}"

    shutil.copy2(igv_js, session_dir / "igv.min.js")
    shutil.copy2(viewer_html, session_dir / "viewer.html")

    ref_name = "reference.fa"
    fai_name = "reference.fa.fai"
    _materialize_file(reference, session_dir / ref_name)
    ensure_fasta_index(session_dir / ref_name)
    _verify_session_file(session_dir, fai_name)

    try:
        track = _prepare_alignment_track(session_dir, alignment, base_url=base_url)
    except RuntimeError as exc:
        raise IgvPreviewError(str(exc)) from exc
    annotation_tracks = _prepare_annotation_tracks(
        session_dir,
        annotations,
        base_url=base_url,
    )
    config = {
        # Enable IGV's built-in "Save SVG" control in the preview window.
        "showSVGButton": True,
        "reference": {
            "id": reference.stem,
            "name": reference.name,
            "fastaURL": f"{base_url}/{ref_name}",
            "indexURL": f"{base_url}/{fai_name}",
        },
        "tracks": [*annotation_tracks, track],
        "locus": "All",
    }
    (session_dir / "config.json").write_text(
        json.dumps(config, indent=2),
        encoding="utf-8",
    )

    return f"{base_url}/viewer.html"
