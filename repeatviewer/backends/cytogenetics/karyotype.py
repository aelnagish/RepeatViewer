"""Karyotypes: chromosome set, lengths, centromere positions and bands.

Each karyotype is a JSON file ``<name>.karyotype.json`` in the folder's store.
Lengths can be typed in, imported from a measurement table (KaryoScope TSV or
any table with chromosome / length / centromere index or arm columns) or taken
from the assembly sequence lengths.
"""

from __future__ import annotations

import csv
import io
import json
import os
import re
from dataclasses import asdict, dataclass, field
from pathlib import Path

from ..cyto_images import store_dir

SUFFIX = ".karyotype.json"
ORDERS = ("listed", "length", "ci")
BAND_KINDS = ("Secondary constriction (NOR)", "DAPI+ heterochromatin", "CMA+ band")


def levan_type(ci: float) -> str:
    """Chromosome type from the centromere index (100 * short arm / total), Levan et al. 1964."""
    if ci <= 0:
        return ""
    ci = min(ci, 100.0 - ci) if ci > 50 else ci
    if ci >= 37.5:
        return "m"
    if ci >= 25.0:
        return "sm"
    if ci >= 12.5:
        return "st"
    return "t"


@dataclass
class Band:
    kind: str
    start: float  # fraction of the chromosome from the p telomere (0) to the q telomere (1)
    end: float

    def __post_init__(self) -> None:
        self.start = max(0.0, min(1.0, float(self.start)))
        self.end = max(0.0, min(1.0, float(self.end)))
        if self.start > self.end:
            self.start, self.end = self.end, self.start


@dataclass
class KaryoChromosome:
    name: str
    length: float = 0.0
    centromere_index: float = 50.0  # percent, short arm / total
    copies: int | None = None  # None = the karyotype's ploidy
    assembly_seq: str = ""
    bands: list[Band] = field(default_factory=list)

    @property
    def type(self) -> str:
        return levan_type(self.centromere_index)

    @property
    def centromere_fraction(self) -> float:
        """Centromere position from the p telomere, as a fraction of the length."""
        return max(0.0, min(1.0, self.centromere_index / 100.0))


@dataclass
class Karyotype:
    name: str
    species: str = ""
    basic_number: int = 8
    ploidy: int = 2
    unit: str = "µm"
    chromosomes: list[KaryoChromosome] = field(default_factory=list)
    homologue_letters: str = "abcdefghijkl"

    def copies_of(self, chrom: KaryoChromosome) -> int:
        return self.ploidy if chrom.copies is None else max(0, int(chrom.copies))

    @property
    def total_count(self) -> int:
        return sum(self.copies_of(c) for c in self.chromosomes)

    @property
    def formula(self) -> str:
        base = self.ploidy * self.basic_number
        total = self.total_count if self.chromosomes else base
        head = f"2n = {self.ploidy}x"
        diff = total - base
        if diff:
            head += f" {'+' if diff > 0 else '-'} {abs(diff)}"
        return f"{head} = {total}"

    def homologue_labels(self, chrom: KaryoChromosome) -> list[str]:
        letters = self.homologue_letters
        n = self.copies_of(chrom)
        return [f"{chrom.name}{letters[i] if i < len(letters) else i + 1}" for i in range(n)]

    def all_labels(self) -> list[str]:
        return [label for c in self.chromosomes for label in self.homologue_labels(c)]

    def chromosome(self, name: str) -> KaryoChromosome | None:
        return next((c for c in self.chromosomes if c.name == name), None)

    def set_basic_number(self, n: int) -> None:
        """Add (named by number) or drop chromosomes at the end so there are ``n``."""
        n = max(1, int(n))
        self.basic_number = n
        while len(self.chromosomes) < n:
            used = {c.name for c in self.chromosomes}
            k = len(self.chromosomes) + 1
            while str(k) in used:
                k += 1
            self.chromosomes.append(KaryoChromosome(name=str(k)))
        del self.chromosomes[n:]

    def ordered(self, order: str = "listed") -> list[KaryoChromosome]:
        chroms = list(self.chromosomes)
        if order == "length":
            chroms.sort(key=lambda c: -c.length)
        elif order == "ci":
            chroms.sort(key=lambda c: -c.centromere_index)
        return chroms

    def relative_lengths(self) -> dict[str, float]:
        """Percent of the haploid set length (one copy of each chromosome)."""
        total = sum(c.length for c in self.chromosomes)
        if total <= 0:
            return {c.name: 0.0 for c in self.chromosomes}
        return {c.name: 100.0 * c.length / total for c in self.chromosomes}

    def to_dict(self) -> dict:
        return {"version": 1, **asdict(self)}

    @classmethod
    def from_dict(cls, raw: dict) -> "Karyotype":
        chroms = []
        for item in raw.get("chromosomes") or []:
            if not isinstance(item, dict) or not str(item.get("name", "")).strip():
                continue
            bands = []
            for b in item.get("bands") or []:
                try:
                    bands.append(Band(str(b.get("kind", "")), float(b["start"]), float(b["end"])))
                except (KeyError, TypeError, ValueError, AttributeError):
                    continue
            copies = item.get("copies")
            chroms.append(KaryoChromosome(
                name=str(item["name"]).strip(),
                length=_num(item.get("length")),
                centromere_index=_num(item.get("centromere_index"), 50.0),
                copies=int(copies) if copies not in (None, "") else None,
                assembly_seq=str(item.get("assembly_seq") or ""),
                bands=bands,
            ))
        return cls(
            name=str(raw.get("name") or "Karyotype"),
            species=str(raw.get("species") or ""),
            basic_number=int(_num(raw.get("basic_number"), len(chroms) or 1)),
            ploidy=max(1, int(_num(raw.get("ploidy"), 2))),
            unit=str(raw.get("unit") or "µm"),
            chromosomes=chroms,
            homologue_letters=str(raw.get("homologue_letters") or "abcdefghijkl"),
        )


def _num(value, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def new_karyotype(name: str, basic_number: int, ploidy: int, species: str = "") -> Karyotype:
    k = Karyotype(name=name, species=species, ploidy=max(1, int(ploidy)))
    k.set_basic_number(basic_number)
    return k


# ---------------------------------------------------------------------------
# Store
# ---------------------------------------------------------------------------


def _safe_name(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", name.strip()).strip("._") or "karyotype"


def karyotype_path(folder: Path, name: str) -> Path:
    return store_dir(Path(folder).expanduser().resolve()) / f"{_safe_name(name)}{SUFFIX}"


def list_karyotypes(folder: Path) -> list[Karyotype]:
    directory = store_dir(Path(folder).expanduser().resolve())
    out = []
    for path in sorted(directory.glob(f"*{SUFFIX}")) if directory.is_dir() else []:
        try:
            out.append(Karyotype.from_dict(json.loads(path.read_text(encoding="utf-8"))))
        except (OSError, ValueError, AttributeError):
            continue
    return out


def load_karyotype(folder: Path, name: str) -> Karyotype | None:
    return next((k for k in list_karyotypes(folder) if k.name == name), None)


def save_karyotype(folder: Path, karyotype: Karyotype, *, old_name: str = "") -> Path:
    path = karyotype_path(folder, karyotype.name)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(karyotype.to_dict(), indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    os.replace(tmp, path)
    if old_name and old_name != karyotype.name:
        old = karyotype_path(folder, old_name)
        if old != path:
            try:
                old.unlink()
            except OSError:
                pass
    return path


def delete_karyotype(folder: Path, name: str) -> bool:
    try:
        karyotype_path(folder, name).unlink()
        return True
    except OSError:
        return False


# ---------------------------------------------------------------------------
# Measurements (KaryoScope or any table) and assembly lengths
# ---------------------------------------------------------------------------

_NAME_COLS = ("chromosome", "chr", "chrom", "name", "id", "label")
_LENGTH_COLS = ("length", "total", "total length", "total_length", "tl", "length_um", "length (µm)", "length (um)")
_CI_COLS = ("ci", "centromere index", "centromeric index", "centromere_index", "centromeric_index")
_P_COLS = ("p", "short arm", "short_arm", "s", "p arm", "p_arm", "p length", "p_length")
_Q_COLS = ("q", "long arm", "long_arm", "l", "q arm", "q_arm", "q length", "q_length")


def _find(header: list[str], names: tuple[str, ...]) -> int:
    for i, h in enumerate(header):
        if h in names:
            return i
    for i, h in enumerate(header):
        base = re.sub(r"\s*\(.*\)$", "", h)
        if base in names:
            return i
    return -1


def chromosome_base_name(label: str) -> str:
    """'2b' -> '2', 'Chr5a' -> 'Chr5', 'A' -> 'A' (homologue letter after a digit is dropped)."""
    label = label.strip()
    m = re.match(r"^(.*\d)([a-zA-Z])$", label)
    return m.group(1) if m else label


def parse_measurement_table(text: str) -> list[tuple[str, float, float | None]]:
    """(chromosome, mean length, mean centromere index or None) per chromosome.

    Rows per homologue ("1a", "1b") are averaged into one chromosome ("1").
    """
    lines = [ln for ln in text.splitlines() if ln.strip() and not ln.lstrip().startswith("#")]
    if not lines:
        return []
    sample = "\n".join(lines[:5])
    delimiter = "\t" if "\t" in sample else (";" if sample.count(";") > sample.count(",") else ",")
    rows = list(csv.reader(io.StringIO("\n".join(lines)), delimiter=delimiter))
    header = [h.strip().lower() for h in rows[0]]
    i_name = _find(header, _NAME_COLS)
    i_len, i_ci = _find(header, _LENGTH_COLS), _find(header, _CI_COLS)
    i_p, i_q = _find(header, _P_COLS), _find(header, _Q_COLS)
    if i_name < 0 or (i_len < 0 and (i_p < 0 or i_q < 0)):
        raise ValueError("The table needs a chromosome column and a length column (or p and q arm columns).")
    sums: dict[str, list[float]] = {}
    order: list[str] = []
    for row in rows[1:]:
        if len(row) <= max(i_name, i_len, i_ci, i_p, i_q):
            continue
        name = chromosome_base_name(row[i_name])
        if not name:
            continue
        try:
            p = float(row[i_p]) if i_p >= 0 else None
            q = float(row[i_q]) if i_q >= 0 else None
            length = float(row[i_len]) if i_len >= 0 else (p or 0.0) + (q or 0.0)
            if i_ci >= 0 and row[i_ci].strip():
                ci = float(row[i_ci])
                ci = ci * 100.0 if ci <= 1.0 else ci
            elif p is not None and q is not None and p + q > 0:
                ci = 100.0 * min(p, q) / (p + q)
            else:
                ci = None
        except ValueError:
            continue
        if name not in sums:
            sums[name] = [0.0, 0.0, 0.0, 0.0]  # length sum, n, ci sum, ci n
            order.append(name)
        acc = sums[name]
        acc[0] += length
        acc[1] += 1
        if ci is not None:
            acc[2] += ci
            acc[3] += 1
    return [(n, sums[n][0] / sums[n][1], (sums[n][2] / sums[n][3]) if sums[n][3] else None) for n in order]


def apply_measurements(karyotype: Karyotype, rows: list[tuple[str, float, float | None]], unit: str = "µm") -> list[str]:
    """Set lengths (and centromere indices) by chromosome name. Returns names not in the karyotype."""
    missing = []
    for name, length, ci in rows:
        chrom = karyotype.chromosome(name)
        if chrom is None:
            missing.append(name)
            continue
        chrom.length = round(float(length), 4)
        if ci is not None:
            chrom.centromere_index = round(float(ci), 2)
    karyotype.unit = unit
    return missing


def apply_assembly_lengths(karyotype: Karyotype, lengths_bp: dict[str, int]) -> list[str]:
    """Lengths in Mb from the linked assembly sequences. Returns chromosomes without a linked sequence."""
    missing = []
    for chrom in karyotype.chromosomes:
        bp = lengths_bp.get(chrom.assembly_seq) if chrom.assembly_seq else None
        if bp is None:
            missing.append(chrom.name)
            continue
        chrom.length = round(bp / 1e6, 3)
    karyotype.unit = "Mb"
    return missing


def guess_assembly_links(karyotype: Karyotype, seq_names: list[str]) -> int:
    """Link 'Chr1', 'chr01', 'Csat_chr1' ... to chromosome '1'. Returns how many were set."""
    by_number: dict[str, str] = {}
    for seq in seq_names:
        m = re.search(r"chr(?:omosome)?[_\- ]?0*(\d+)$", seq, re.IGNORECASE)
        if m and m.group(1) not in by_number:
            by_number[m.group(1)] = seq
    n = 0
    for chrom in karyotype.chromosomes:
        if chrom.assembly_seq:
            continue
        key = chrom.name.lstrip("0") or chrom.name
        if key in by_number:
            chrom.assembly_seq = by_number[key]
            n += 1
    return n


# ---------------------------------------------------------------------------
# Bands in arm notation ("q 0.92" = 92 % of the way from the centromere to the q telomere)
# ---------------------------------------------------------------------------


def fraction_to_arm_text(fraction: float, centromere_fraction: float) -> str:
    c = max(0.0, min(1.0, centromere_fraction))
    f = max(0.0, min(1.0, fraction))
    if f <= c:
        return f"p {((c - f) / c if c > 0 else 0.0):.2f}"
    return f"q {((f - c) / (1.0 - c) if c < 1 else 0.0):.2f}"


def arm_text_to_fraction(text: str, centromere_fraction: float) -> float | None:
    """Inverse of ``fraction_to_arm_text``. Returns None when the text cannot be read."""
    m = re.match(r"^\s*([pPqQ])\s*([0-9]*\.?[0-9]+)\s*$", str(text).replace(",", "."))
    if not m:
        return None
    arm, rel = m.group(1).lower(), max(0.0, min(1.0, float(m.group(2))))
    c = max(0.0, min(1.0, centromere_fraction))
    return c - rel * c if arm == "p" else c + rel * (1.0 - c)


def band_summary(chrom: KaryoChromosome) -> str:
    if not chrom.bands:
        return "none"
    names = {"Secondary constriction (NOR)": "NOR", "DAPI+ heterochromatin": "DAPI+", "CMA+ band": "CMA+"}
    out = []
    for b in chrom.bands:
        mid = (b.start + b.end) / 2
        arm = "p" if mid < chrom.centromere_fraction else "q"
        out.append(f"{names.get(b.kind, b.kind)} ({arm})")
    return ", ".join(out)


FASTA_SUFFIXES = (".fa", ".fasta", ".fna", ".fas")


def folder_sequence_lengths(folder: Path) -> dict[str, int]:
    """Record name -> length (bp) for the FASTA files directly in ``folder``."""
    from ..sequence_features import fasta_record_lengths

    lengths: dict[str, int] = {}
    for path in sorted(Path(folder).expanduser().iterdir()):
        if path.is_file() and path.suffix.lower() in FASTA_SUFFIXES:
            try:
                for name, n in fasta_record_lengths(path).items():
                    lengths.setdefault(name, int(n))
            except (OSError, ValueError):
                continue
    return lengths
