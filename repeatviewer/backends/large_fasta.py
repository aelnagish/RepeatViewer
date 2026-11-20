"""Large and compressed FASTA files: index, sequence list, region reads.

Qt-free. A FASTA is indexed once (samtools-style ``.fai``) and then any region
of any sequence is read from disk by offset, so a 4 GB genome never has to be
loaded into memory.

Three storage kinds are handled the same way through :func:`open_seekable`:

- plain FASTA: offsets are file offsets (samtools compatible ``.fai``);
- BGZF (``bgzip``, still a valid ``.gz``): offsets are uncompressed offsets,
  turned into block positions with a ``.gzi`` table (samtools compatible);
- plain gzip: offsets are uncompressed offsets too, but gzip cannot jump, so a
  read decompresses from the start of the file. Fine for libraries and small
  genomes; for large ones, convert to bgzip once (:func:`convert_to_bgzip`).

Positions in the public functions are 1-based and inclusive, like the rest of
the app ("scaffold_5:120,000,001-125,000,000").
"""

from __future__ import annotations

import bisect
import gzip
import hashlib
import io
import os
import struct
import zlib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterator

FASTA_SUFFIXES = (".fa", ".fasta", ".fas", ".fna", ".ffn", ".faa", ".fsa")
COMPRESSED_SUFFIXES = (".gz", ".bgz")

KIND_PLAIN = "plain"
KIND_BGZF = "bgzf"
KIND_GZIP = "gzip"

# Annotation view limit for one sequence (same as sequence_map_viewer).
WHOLE_SEQUENCE_LIMIT = 100_000_000
# Window opened when a longer sequence is double-clicked (decided 16 Nov 2026).
DEFAULT_WINDOW = 5_000_000

_CHUNK = 8 * 1024 * 1024
_BGZF_MAX_BLOCK_DATA = 65280  # what bgzip packs into one block

ProgressFn = Callable[[int, int], None]  # (done bytes, total bytes)
CancelFn = Callable[[], bool]


class IndexCancelled(Exception):
    """Raised when indexing was stopped by the caller."""


class LargeFastaError(ValueError):
    """A FASTA problem the user can act on (shown as is)."""


# --------------------------------------------------------------------------
# Names and kinds
# --------------------------------------------------------------------------


def split_compression(path: Path) -> tuple[str, bool]:
    """Return (name without .gz/.bgz, compressed?)."""
    name = Path(path).name
    low = name.lower()
    for suffix in COMPRESSED_SUFFIXES:
        if low.endswith(suffix):
            return name[: -len(suffix)], True
    return name, False


def inner_suffix(path: Path) -> str:
    """Suffix of the file as if it were not compressed ("x.gff3.gz" -> ".gff3")."""
    base, _ = split_compression(path)
    return Path(base).suffix.lower()


def is_fasta_name(path: Path) -> bool:
    """FASTA by name, compressed or not (.fa, .fasta, .fna, ... and .gz / .bgz)."""
    return inner_suffix(path) in FASTA_SUFFIXES


def storage_kind(path: Path) -> str:
    """KIND_PLAIN, KIND_BGZF or KIND_GZIP, from the first bytes of the file."""
    try:
        with Path(path).open("rb") as fh:
            head = fh.read(18)
    except OSError:
        return KIND_PLAIN
    if len(head) < 2 or head[:2] != b"\x1f\x8b":
        return KIND_PLAIN
    if len(head) >= 18 and head[3] & 4 and head[12:14] == b"BC":
        return KIND_BGZF
    return KIND_GZIP


# --------------------------------------------------------------------------
# Index
# --------------------------------------------------------------------------


@dataclass
class FaiEntry:
    name: str
    length: int
    offset: int  # uncompressed offset of the first base
    line_bases: int = 0  # 0 = lines not regular, read by streaming
    line_bytes: int = 0
    description: str = ""


@dataclass
class FastaIndex:
    path: Path
    kind: str
    entries: list[FaiEntry] = field(default_factory=list)
    blocks: list[tuple[int, int]] = field(default_factory=list)  # BGZF (coffset, uoffset)

    def names(self) -> list[str]:
        return [e.name for e in self.entries]

    def entry(self, name: str) -> FaiEntry:
        for e in self.entries:
            if e.name == name:
                return e
        raise LargeFastaError(f'There is no sequence "{name}" in {self.path.name}.')

    @property
    def total_length(self) -> int:
        return sum(e.length for e in self.entries)

    @property
    def longest(self) -> int:
        return max((e.length for e in self.entries), default=0)


def _cache_dir() -> Path:
    base = os.environ.get("XDG_CACHE_HOME") or str(Path.home() / ".cache")
    return Path(base) / "RepeatViewer" / "fasta_index"


def _cache_stem(path: Path) -> Path:
    key = hashlib.sha1(str(Path(path).resolve()).encode("utf-8", "replace")).hexdigest()[:16]
    return _cache_dir() / f"{key}_{Path(path).name}"


def index_paths(path: Path) -> list[Path]:
    """Places an index may live: next to the FASTA (samtools style), then the cache."""
    path = Path(path)
    return [path.with_name(path.name + ".fai"), _cache_stem(path).with_name(_cache_stem(path).name + ".fai")]


def _fresh(index_file: Path, source: Path) -> bool:
    try:
        return index_file.is_file() and index_file.stat().st_mtime >= source.stat().st_mtime
    except OSError:
        return False


def _read_fai(fai: Path) -> list[FaiEntry]:
    out: list[FaiEntry] = []
    for line in fai.read_text(encoding="utf-8").splitlines():
        cols = line.split("\t")
        if len(cols) >= 5:
            out.append(FaiEntry(cols[0], int(cols[1]), int(cols[2]), int(cols[3]), int(cols[4])))
    return out


def _write_fai(fai: Path, entries: list[FaiEntry]) -> None:
    fai.parent.mkdir(parents=True, exist_ok=True)
    text = "".join(
        f"{e.name}\t{e.length}\t{e.offset}\t{e.line_bases}\t{e.line_bytes}\n" for e in entries
    )
    tmp = fai.with_name(fai.name + ".tmp")
    tmp.write_text(text, encoding="utf-8")
    os.replace(tmp, fai)


def _read_gzi(gzi: Path) -> list[tuple[int, int]]:
    data = gzi.read_bytes()
    (n,) = struct.unpack_from("<Q", data, 0)
    blocks = [(0, 0)]
    for i in range(n):
        c, u = struct.unpack_from("<QQ", data, 8 + 16 * i)
        blocks.append((c, u))
    return blocks


def _write_gzi(gzi: Path, blocks: list[tuple[int, int]]) -> None:
    rest = [b for b in blocks if b != (0, 0)]
    payload = struct.pack("<Q", len(rest)) + b"".join(struct.pack("<QQ", c, u) for c, u in rest)
    tmp = gzi.with_name(gzi.name + ".tmp")
    tmp.write_bytes(payload)
    os.replace(tmp, gzi)


def load_index(path: Path) -> FastaIndex | None:
    """The saved index when it is newer than the FASTA, else None."""
    path = Path(path)
    kind = storage_kind(path)
    for fai in index_paths(path):
        if not _fresh(fai, path):
            continue
        try:
            entries = _read_fai(fai)
        except (OSError, ValueError):
            continue
        if not entries:
            continue
        index = FastaIndex(path, kind, entries)
        if kind == KIND_BGZF:
            gzi = fai.with_name(fai.name[: -len(".fai")] + ".gzi")
            if not _fresh(gzi, path):
                continue
            try:
                index.blocks = _read_gzi(gzi)
            except (OSError, ValueError, struct.error):
                continue
        return index
    return None


def _bgzf_block_table(path: Path, cancel: CancelFn | None = None) -> list[tuple[int, int]]:
    """(compressed offset, uncompressed offset) of every BGZF block, from headers only."""
    blocks: list[tuple[int, int]] = []
    coff = uoff = 0
    size = Path(path).stat().st_size
    with Path(path).open("rb") as fh:
        while coff < size:
            if cancel is not None and cancel():
                raise IndexCancelled()
            fh.seek(coff)
            head = fh.read(18)
            if len(head) < 18 or head[:2] != b"\x1f\x8b":
                raise LargeFastaError(f"{Path(path).name} is not a valid bgzip file.")
            xlen = struct.unpack_from("<H", head, 10)[0]
            fh.seek(coff + 12)
            extra = fh.read(xlen)
            bsize = None
            i = 0
            while i + 4 <= len(extra):
                si1, si2, slen = extra[i], extra[i + 1], struct.unpack_from("<H", extra, i + 2)[0]
                if si1 == 66 and si2 == 67 and slen == 2:
                    bsize = struct.unpack_from("<H", extra, i + 4)[0]
                i += 4 + slen
            if bsize is None:
                raise LargeFastaError(f"{Path(path).name} is not a valid bgzip file.")
            fh.seek(coff + bsize + 1 - 4)
            isize = struct.unpack("<I", fh.read(4))[0]
            blocks.append((coff, uoff))
            coff += bsize + 1
            uoff += isize
    return blocks


def _raw_chunks(path: Path, kind: str, progress: ProgressFn | None,
                cancel: CancelFn | None) -> Iterator[bytes]:
    """Uncompressed bytes of the file in chunks, reporting progress on the file size."""
    total = Path(path).stat().st_size
    with Path(path).open("rb") as raw:
        stream = gzip.GzipFile(fileobj=raw) if kind != KIND_PLAIN else raw
        while True:
            if cancel is not None and cancel():
                raise IndexCancelled()
            chunk = stream.read(_CHUNK)
            if not chunk:
                break
            if progress is not None:
                progress(raw.tell(), total)
            yield chunk


@dataclass
class _Layout:
    """Per-record bookkeeping while scanning, to prove the lines are regular."""

    first: bytearray = field(default_factory=bytearray)
    rel: int = 0  # bytes of the record seen so far
    grid: int = 0  # newlines exactly at line_bytes, 2 * line_bytes, ...
    first_bad: int | None = None  # rel of the first newline off that grid
    bases_at_bad: int = 0


def _nbases(seg: bytes) -> int:
    return len(seg) - seg.count(b"\n") - seg.count(b"\r") - seg.count(b" ")


def _feed(cur: FaiEntry, st: _Layout, seg: bytes) -> None:
    if not seg:
        return
    if cur.line_bytes == 0:  # first line not finished yet
        k = seg.find(b"\n")
        if k < 0:
            st.first += seg
            cur.length += _nbases(seg)
            st.rel += len(seg)
            return
        head = seg[: k + 1]
        st.first += head
        cur.line_bytes = len(st.first)
        cur.line_bases = len(st.first.rstrip(b"\r\n"))
        cur.length += _nbases(head)
        st.rel += len(head)
        st.grid += 1
        seg = seg[k + 1:]
        if not seg:
            return
    lB = cur.line_bytes
    phase = (lB - 1 - st.rel) % lB
    grid = seg[phase::lB].count(b"\n")
    total = seg.count(b"\n")
    if st.first_bad is None and total > grid:
        q = seg.find(b"\n")
        while q >= 0 and (st.rel + q + 1) % lB == 0:
            q = seg.find(b"\n", q + 1)
        st.first_bad = st.rel + q
        st.bases_at_bad = cur.length + _nbases(seg[: q + 1])
    st.grid += grid
    cur.length += _nbases(seg)
    st.rel += len(seg)


def _regular(e: FaiEntry, st: _Layout) -> bool:
    """True when the record is laid out in equal lines, so offsets can be computed."""
    if e.length == 0:
        return True
    lb, lB = e.line_bases, e.line_bytes
    if lb <= 0 or lB < lb:
        return False
    full, rest = divmod(e.length, lb)
    if st.first_bad is not None:
        if st.first_bad < full * lB or st.bases_at_bad != e.length:
            return False
    return st.grid >= full - (0 if rest else 1)


def _scan(chunks: Iterator[bytes]) -> list[FaiEntry]:
    """One pass over the bytes: names, lengths, offsets, line shape, regularity."""
    entries: list[FaiEntry] = []
    layouts: list[_Layout] = []
    pos = 0  # uncompressed offset of chunk start
    in_header = False
    header = bytearray()
    cur: FaiEntry | None = None
    st: _Layout | None = None
    prev_newline = True

    def new_entry(text: str, offset: int) -> tuple[FaiEntry, _Layout]:
        words = text.rstrip("\r").split(None, 1)
        e = FaiEntry(words[0] if words else "", 0, offset, 0, 0,
                     words[1] if len(words) > 1 else "")
        entries.append(e)
        layouts.append(_Layout())
        return e, layouts[-1]

    for chunk in chunks:
        i = 0
        n = len(chunk)
        while i < n:
            if in_header:
                j = chunk.find(b"\n", i)
                if j < 0:
                    header += chunk[i:]
                    break
                header += chunk[i:j]
                cur, st = new_entry(header.decode("utf-8", "replace"), pos + j + 1)
                in_header = False
                header = bytearray()
                i = j + 1
                prev_newline = True
                continue
            if prev_newline and chunk[i] == 62:  # ">" at a line start
                in_header = True
                i += 1
                continue
            h = chunk.find(b"\n>", i)
            seg_end = h + 1 if h >= 0 else n
            seg = chunk[i:seg_end]
            if cur is not None:
                _feed(cur, st, seg)
            elif seg.strip():
                raise LargeFastaError("This file does not start with a FASTA header (\">name\").")
            prev_newline = seg.endswith(b"\n")
            i = seg_end
        pos += n
    if in_header:
        new_entry(header.decode("utf-8", "replace"), pos)
    if cur is not None and st is not None and cur.line_bytes == 0 and st.first:
        cur.line_bytes = len(st.first)
        cur.line_bases = len(st.first.rstrip(b"\r\n"))
        st.grid = 1 if st.first.endswith(b"\n") else 0

    for entry, layout in zip(entries, layouts):
        if not _regular(entry, layout):
            entry.line_bases = entry.line_bytes = 0
    return entries


def build_index(path: Path, progress: ProgressFn | None = None,
                cancel: CancelFn | None = None) -> FastaIndex:
    """Read the whole file once and save the index (next to it, else in the cache)."""
    path = Path(path)
    kind = storage_kind(path)
    entries = _scan(_raw_chunks(path, kind, progress, cancel))
    if not entries:
        raise LargeFastaError(f"No sequences in {path.name}.")
    names = [e.name for e in entries]
    if len(set(names)) != len(names):
        dup = next(n for n in names if names.count(n) > 1)
        raise LargeFastaError(f'{path.name} has more than one sequence named "{dup}".')
    index = FastaIndex(path, kind, entries)
    if kind == KIND_BGZF:
        index.blocks = _bgzf_block_table(path, cancel)
    regular = all(e.line_bases > 0 or e.length == 0 for e in entries)
    # samtools cannot use an index next to a plain-gzip file, and an index of
    # irregular lines is ours only: both stay in the cache.
    beside = regular and kind != KIND_GZIP
    targets = index_paths(path) if beside else index_paths(path)[1:]
    for fai in targets:
        try:
            _write_fai(fai, entries)
            if kind == KIND_BGZF:
                _write_gzi(fai.with_name(fai.name[: -len(".fai")] + ".gzi"), index.blocks)
            break
        except OSError:
            continue
    return index


def get_index(path: Path, progress: ProgressFn | None = None,
              cancel: CancelFn | None = None) -> FastaIndex:
    """The saved index, or a new one. Call from a worker thread for large files."""
    return load_index(path) or build_index(path, progress, cancel)


# --------------------------------------------------------------------------
# Reading regions
# --------------------------------------------------------------------------


class _BgzfReader(io.RawIOBase):
    """Seekable reader of a BGZF file by uncompressed offset (uses the block table)."""

    def __init__(self, path: Path, blocks: list[tuple[int, int]]) -> None:
        super().__init__()
        self._fh = Path(path).open("rb")
        self._blocks = blocks
        self._uoffs = [u for _, u in blocks]
        self._pos = 0
        self._cache_i = -1
        self._cache = b""

    def _block(self, i: int) -> bytes:
        if i != self._cache_i:
            coff = self._blocks[i][0]
            end = self._blocks[i + 1][0] if i + 1 < len(self._blocks) else None
            self._fh.seek(coff)
            raw = self._fh.read() if end is None else self._fh.read(end - coff)
            self._cache = zlib.decompress(raw[: _member_size(raw)], 31)
            self._cache_i = i
        return self._cache

    def seek(self, offset: int, whence: int = 0) -> int:
        self._pos = offset
        return offset

    def tell(self) -> int:
        return self._pos

    def read(self, n: int = -1) -> bytes:
        out = bytearray()
        while n < 0 or len(out) < n:
            i = bisect.bisect_right(self._uoffs, self._pos) - 1
            if i < 0 or i >= len(self._blocks):
                break
            data = self._block(i)
            start = self._pos - self._uoffs[i]
            if start >= len(data):
                if i + 1 >= len(self._blocks):
                    break
                self._pos = self._uoffs[i + 1]
                continue
            take = data[start:] if n < 0 else data[start: start + (n - len(out))]
            out += take
            self._pos += len(take)
        return bytes(out)

    def readable(self) -> bool:
        return True

    def close(self) -> None:
        try:
            self._fh.close()
        finally:
            super().close()


def _member_size(raw: bytes) -> int:
    xlen = struct.unpack_from("<H", raw, 10)[0]
    extra = raw[12: 12 + xlen]
    i = 0
    while i + 4 <= len(extra):
        slen = struct.unpack_from("<H", extra, i + 2)[0]
        if extra[i] == 66 and extra[i + 1] == 67 and slen == 2:
            return struct.unpack_from("<H", extra, i + 4)[0] + 1
        i += 4 + slen
    return len(raw)


def open_seekable(index: FastaIndex):
    """A binary handle with seek(uncompressed offset) / read(n)."""
    if index.kind == KIND_BGZF:
        return _BgzfReader(index.path, index.blocks)
    if index.kind == KIND_GZIP:
        return gzip.open(index.path, "rb")
    return Path(index.path).open("rb")


def fetch(index: FastaIndex, name: str, start: int = 1, end: int | None = None,
          handle=None) -> str:
    """Bases start..end (1-based, inclusive) of one sequence, as stored (case kept)."""
    entry = index.entry(name)
    end = entry.length if end is None else min(int(end), entry.length)
    start = max(1, int(start))
    if end < start:
        return ""
    own = handle is None
    fh = open_seekable(index) if own else handle
    try:
        s0 = start - 1
        if entry.line_bases > 0:
            lb, lB = entry.line_bases, entry.line_bytes
            first = entry.offset + (s0 // lb) * lB + s0 % lb
            last = entry.offset + ((end - 1) // lb) * lB + (end - 1) % lb
            fh.seek(first)
            raw = fh.read(last - first + 1)
            return _clean(raw)
        # Irregular lines: walk from the start of the record.
        fh.seek(entry.offset)
        need = end
        got = bytearray()
        while len(got) < need:
            chunk = fh.read(1 << 20)
            if not chunk:
                break
            k = chunk.find(b"\n>")
            if k >= 0:
                chunk = chunk[: k + 1]
            got += _clean(chunk).encode("ascii")
            if k >= 0:
                break
        return got[s0:end].decode("ascii")
    finally:
        if own:
            fh.close()


def _clean(raw: bytes) -> str:
    return (raw.replace(b"\n", b"").replace(b"\r", b"").replace(b" ", b"")
            .decode("ascii", errors="replace"))


# --------------------------------------------------------------------------
# Regions and output files
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class Region:
    name: str
    start: int  # 1-based, inclusive
    end: int

    @property
    def length(self) -> int:
        return self.end - self.start + 1

    def label(self) -> str:
        return f"{self.name}:{self.start:,}-{self.end:,}"

    def header_id(self) -> str:
        return f"{self.name}:{self.start}-{self.end}"


def parse_region(text: str, index: FastaIndex, current: str | None = None) -> Region:
    """Read "scaffold_5:120,000,000-125,000,000", "scaffold_5", or "1000-2000" (current
    sequence). Clamps to the sequence; raises LargeFastaError with a readable message."""
    from repeatviewer.backends.motif_search import parse_position

    raw = (text or "").strip()
    if not raw:
        raise LargeFastaError("Type a region, e.g. scaffold_5:120,000,000-125,000,000.")
    names = set(index.names())
    if raw in names:
        e = index.entry(raw)
        return Region(raw, 1, max(1, e.length))
    parsed = parse_position(raw)
    if parsed is None:
        raise LargeFastaError(f'"{raw}" is not a region. Example: scaffold_5:1-5,000,000.')
    name, start, end = parsed
    name = name or current
    if not name:
        raise LargeFastaError("This file has several sequences: type the name too, e.g. scaffold_5:1000-2000.")
    e = index.entry(name)
    if start > e.length:
        raise LargeFastaError(f"Position {start:,} is past the end of {name} ({e.length:,} bp).")
    return Region(name, start, min(end, e.length))


def opening_region(index: FastaIndex, name: str, window: int | None = None) -> Region:
    """What a double-click opens: the whole sequence up to the Annotation limit,
    otherwise the first window (decided 16 Nov 2026)."""
    e = index.entry(name)
    if e.length <= WHOLE_SEQUENCE_LIMIT:
        return Region(name, 1, max(1, e.length))
    window = DEFAULT_WINDOW if window is None else window
    return Region(name, 1, min(e.length, max(1, int(window))))


def region_output_path(index: FastaIndex, region: Region) -> Path:
    """`<input without .gz>_<name>.fa` for a whole sequence, else
    `<input without .gz>_<name>_<start>-<end>.fa`, next to the input."""
    base, _ = split_compression(index.path)
    e = index.entry(region.name)
    safe = "".join(c if c.isalnum() or c in "._-" else "_" for c in region.name)
    if region.start == 1 and region.end >= e.length:
        stem = f"{base}_{safe}"
    else:
        stem = f"{base}_{safe}_{region.start}-{region.end}"
    return Path(index.path).with_name(stem + ".fa")


def numbered_path(path: Path) -> Path:
    """`x.fa`, then `x_2.fa`, `x_3.fa` ... (nothing is overwritten)."""
    if not path.exists():
        return path
    n = 2
    while True:
        candidate = path.with_name(f"{path.stem}_{n}{path.suffix}")
        if not candidate.exists():
            return candidate
        n += 1


def write_region(index: FastaIndex, region: Region, out_path: Path, *, line_width: int = 60,
                 progress: ProgressFn | None = None, cancel: CancelFn | None = None) -> Path:
    """Write one region as FASTA. The header keeps genome positions unless the region
    is the whole sequence (then the original name and description)."""
    e = index.entry(region.name)
    whole = region.start == 1 and region.end >= e.length
    header = (f">{e.name} {e.description}".rstrip() if whole
              else f">{region.header_id()} from {Path(index.path).name}")
    out_path = Path(out_path)
    tmp = out_path.with_name(out_path.name + ".part")
    step = max(line_width, (1 << 20) // line_width * line_width)
    with open_seekable(index) as fh, tmp.open("w", encoding="ascii", newline="\n") as out:
        out.write(header + "\n")
        pos = region.start
        while pos <= region.end:
            if cancel is not None and cancel():
                out.close()
                tmp.unlink(missing_ok=True)
                raise IndexCancelled()
            stop = min(region.end, pos + step - 1)
            seq = fetch(index, region.name, pos, stop, handle=fh)
            for i in range(0, len(seq), line_width):
                out.write(seq[i: i + line_width] + "\n")
            if progress is not None:
                progress(stop - region.start + 1, region.length)
            pos = stop + 1
    os.replace(tmp, out_path)
    return out_path


def convert_to_bgzip(path: Path, out_path: Path | None = None,
                     progress: ProgressFn | None = None,
                     cancel: CancelFn | None = None) -> Path:
    """Re-compress a plain-gzip (or plain) FASTA as BGZF, so windows open instantly.
    Default output next to the input: `<stem>_bgzip.fa.gz` (numbered if it exists)."""
    path = Path(path)
    base, _ = split_compression(path)
    if out_path is None:
        stem, suffix = os.path.splitext(base)
        out_path = numbered_path(path.with_name(f"{stem}_bgzip{suffix or '.fa'}.gz"))
    out_path = Path(out_path)
    tmp = out_path.with_name(out_path.name + ".part")
    kind = storage_kind(path)
    pending = bytearray()
    with tmp.open("wb") as out:
        for chunk in _raw_chunks(path, kind, progress, cancel):
            pending += chunk
            while len(pending) >= _BGZF_MAX_BLOCK_DATA:
                out.write(_bgzf_block(bytes(pending[:_BGZF_MAX_BLOCK_DATA])))
                del pending[:_BGZF_MAX_BLOCK_DATA]
        if pending:
            out.write(_bgzf_block(bytes(pending)))
        out.write(_BGZF_EOF)
    os.replace(tmp, out_path)
    return out_path


_BGZF_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")


def _bgzf_block(data: bytes) -> bytes:
    comp = zlib.compressobj(6, zlib.DEFLATED, -15)
    body = comp.compress(data) + comp.flush()
    bsize = len(body) + 25
    head = (b"\x1f\x8b\x08\x04" + b"\x00\x00\x00\x00" + b"\x00\xff" + struct.pack("<H", 6)
            + b"BC" + struct.pack("<H", 2) + struct.pack("<H", bsize))
    tail = struct.pack("<II", zlib.crc32(data) & 0xFFFFFFFF, len(data) & 0xFFFFFFFF)
    return head + body + tail


# --------------------------------------------------------------------------
# Small helpers for the viewers
# --------------------------------------------------------------------------


def open_fasta_text(path) -> io.TextIOBase:
    """Text handle of a FASTA, decompressing .gz / .bgz on the fly (for Biopython)."""
    p = Path(path)
    if split_compression(p)[1] or storage_kind(p) != KIND_PLAIN:
        return gzip.open(p, "rt", encoding="utf-8", errors="replace")
    return p.open("r", encoding="utf-8", errors="replace")


def format_bp(n: int) -> str:
    """12 bp, 4.5 kb, 438.2 Mb, 3.55 Gb."""
    n = int(n)
    if n < 1_000:
        return f"{n} bp"
    if n < 1_000_000:
        return f"{n / 1_000:.1f} kb"
    if n < 1_000_000_000:
        return f"{n / 1_000_000:.1f} Mb"
    return f"{n / 1_000_000_000:.2f} Gb"
