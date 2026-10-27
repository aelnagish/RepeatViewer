"""Primer records and the Primer Database.

The database is an ordinary FASTA file (``primer_database.fa`` in the
RepeatViewer config folder by default, or any FASTA chosen in the Options
tab), so it stays readable, editable and easy to share::

    >CroSat2_F kind=forward note="ladder PCR, 2026-10"
    GACTTGCAAGTTCGAGAAC
    >M13_CroSat2_R kind=reverse
    caggaaacagctatgacGTTGCATCCACTTAGCTTAG

Lower-case 5' letters are a tail that does not bind the template.
"""

from __future__ import annotations

import re
import shlex
from dataclasses import asdict, dataclass, field
from pathlib import Path

from repeatviewer.backends.primers.thermo import PrimerError, clean_oligo, split_tail

KIND_FORWARD = "forward"
KIND_REVERSE = "reverse"
KIND_PROBE = "probe"
KIND_OLIGO = "oligo"
PRIMER_KINDS = (KIND_FORWARD, KIND_REVERSE, KIND_PROBE, KIND_OLIGO)
KIND_LABELS = {
    KIND_FORWARD: "Forward",
    KIND_REVERSE: "Reverse",
    KIND_PROBE: "Probe",
    KIND_OLIGO: "Oligo",
}

_FWD_NAME = re.compile(r"(?:[_\-. ](?:f|fw|fwd|for|forward|left|l)|[_\-.]?F\d*)$", re.I)
_REV_NAME = re.compile(r"(?:[_\-. ](?:r|rv|rev|reverse|right)|[_\-.]?R\d*)$", re.I)
_PROBE_NAME = re.compile(r"(?:[_\-. ](?:p|probe|oligo|int|internal))$", re.I)


def kind_from_name(name: str) -> str:
    """Guess forward / reverse / probe from a name such as ``CroSat2_F``."""
    text = str(name or "").strip()
    if _PROBE_NAME.search(text):
        return KIND_PROBE
    if re.search(r"(?:[_\-. ](?:f|fw|fwd|for|forward|left))$", text, re.I) or re.search(
        r"[_\-.]F\d*$", text
    ):
        return KIND_FORWARD
    if re.search(r"(?:[_\-. ](?:r|rv|rev|reverse|right))$", text, re.I) or re.search(
        r"[_\-.]R\d*$", text
    ):
        return KIND_REVERSE
    return KIND_OLIGO


@dataclass
class Primer:
    name: str
    sequence: str                # full oligo, lower-case tail kept
    kind: str = ""               # empty: guessed from the name (..._F / ..._R / ..._P)
    note: str = ""
    # Where the primer came from (design run, file, typed); informative only.
    source: str = ""
    extra: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.name = str(self.name or "").strip() or "primer"
        self.sequence = clean_oligo(self.sequence)
        if not self.sequence:
            raise PrimerError(f"{self.name}: no sequence.")
        if self.kind not in PRIMER_KINDS:
            self.kind = kind_from_name(self.name)

    @property
    def tail(self) -> str:
        return split_tail(self.sequence)[0]

    @property
    def binding(self) -> str:
        return split_tail(self.sequence)[1]

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict) -> Primer:
        return cls(
            name=data.get("name", ""),
            sequence=data.get("sequence", ""),
            kind=data.get("kind", KIND_OLIGO),
            note=data.get("note", ""),
            source=data.get("source", ""),
            extra=dict(data.get("extra") or {}),
        )

    def fasta_header(self) -> str:
        parts = [self.name.replace(" ", "_"), f"kind={self.kind}"]
        if self.note:
            parts.append(f'note="{self.note.replace(chr(34), chr(39))}"')
        if self.source:
            parts.append(f'source="{self.source.replace(chr(34), chr(39))}"')
        return " ".join(parts)


def _parse_header(header: str) -> tuple[str, dict]:
    text = header.lstrip(">").strip()
    if not text:
        return "", {}
    name, _, rest = text.partition(" ")
    attrs: dict[str, str] = {}
    try:
        tokens = shlex.split(rest)
    except ValueError:
        tokens = rest.split()
    free = []
    for tok in tokens:
        key, eq, value = tok.partition("=")
        if eq and key in ("kind", "note", "source"):
            attrs[key] = value
        else:
            free.append(tok)
    if free and "note" not in attrs:
        attrs["note"] = " ".join(free)
    return name, attrs


def parse_primers_fasta_text(text: str, *, source: str = "") -> list[Primer]:
    primers: list[Primer] = []
    name, attrs, seq = None, {}, []

    def flush() -> None:
        if name is None:
            return
        joined = "".join(seq)
        if not joined.strip():
            return
        kind = attrs.get("kind") or kind_from_name(name)
        primers.append(Primer(name=name or f"primer_{len(primers) + 1}", sequence=joined,
                              kind=kind, note=attrs.get("note", ""),
                              source=attrs.get("source", source)))

    for line in str(text).splitlines():
        line = line.rstrip()
        if line.startswith(">"):
            flush()
            name, attrs = _parse_header(line)
            seq = []
        elif line.strip() and name is not None:
            seq.append(line.strip())
    flush()
    return primers


def parse_primer_text(text: str, *, default_prefix: str = "primer") -> list[Primer]:
    """Primers pasted in any common shape.

    FASTA; ``name<TAB or spaces>sequence`` lines (as copied from a spreadsheet
    or an order sheet); or one bare sequence per line.
    """
    body = str(text or "").strip()
    if not body:
        return []
    if body.lstrip().startswith(">"):
        return parse_primers_fasta_text(body, source="typed")
    out: list[Primer] = []
    for n, raw in enumerate(body.splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = re.split(r"[\t,;]+|\s{2,}|\s+(?=[ACGTUacgtuRYSWKMBDHVN5'-]{8,}$)", line)
        parts = [p for p in parts if p.strip()]
        if len(parts) >= 2:
            name, seq = parts[0].strip(), parts[1].strip()
        else:
            name, seq = f"{default_prefix}_{n}", parts[0].strip()
        out.append(Primer(name=name, sequence=seq, kind=kind_from_name(name), source="typed"))
    return out


def load_primers(path: Path) -> list[Primer]:
    p = Path(path).expanduser()
    if not p.is_file():
        return []
    text = p.read_text(encoding="utf-8", errors="replace")
    if text.lstrip().startswith(">"):
        return parse_primers_fasta_text(text, source=p.name)
    return parse_primer_text(text)


def write_primers(path: Path, primers: list[Primer]) -> Path:
    p = Path(path).expanduser()
    p.parent.mkdir(parents=True, exist_ok=True)
    lines = []
    for pr in primers:
        lines.append(">" + pr.fasta_header())
        lines.append(pr.sequence)
    p.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")
    return p


def unique_name(name: str, taken: set[str]) -> str:
    if name not in taken:
        return name
    n = 2
    while f"{name}_{n}" in taken:
        n += 1
    return f"{name}_{n}"


def default_database_path() -> Path:
    from repeatviewer.app_preferences import config_dir

    return config_dir() / "primers" / "primer_database.fa"


class PrimerDatabase:
    """The saved primers (one FASTA file)."""

    def __init__(self, path: Path | None = None) -> None:
        self.path = Path(path).expanduser() if path else default_database_path()
        self.primers: list[Primer] = load_primers(self.path)

    def reload(self) -> None:
        self.primers = load_primers(self.path)

    def names(self) -> list[str]:
        return [p.name for p in self.primers]

    def get(self, name: str) -> Primer | None:
        return next((p for p in self.primers if p.name == name), None)

    def find_sequence(self, sequence: str) -> Primer | None:
        target = clean_oligo(sequence).upper()
        return next((p for p in self.primers if p.sequence.upper() == target), None)

    def add(self, primers: list[Primer], *, replace_same_name: bool = False) -> list[Primer]:
        """Add primers; identical sequences are skipped, clashing names get _2."""
        added = []
        for pr in primers:
            same = self.find_sequence(pr.sequence)
            if same is not None and (same.name == pr.name or not replace_same_name):
                continue
            existing = self.get(pr.name)
            if existing is not None and replace_same_name:
                self.primers[self.primers.index(existing)] = pr
                added.append(pr)
                continue
            taken = set(self.names())
            if pr.name in taken:
                pr = Primer(name=unique_name(pr.name, taken), sequence=pr.sequence,
                            kind=pr.kind, note=pr.note, source=pr.source, extra=pr.extra)
            self.primers.append(pr)
            added.append(pr)
        return added

    def remove(self, names: list[str]) -> int:
        before = len(self.primers)
        drop = set(names)
        self.primers = [p for p in self.primers if p.name not in drop]
        return before - len(self.primers)

    def save(self) -> Path:
        return write_primers(self.path, self.primers)
