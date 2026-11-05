"""Karyotypes: chromosome set, lengths, centromere positions and bands.

Each karyotype is a JSON file ``<name>.karyotype.json`` in the folder's store.
Lengths can be typed in, imported from a measurement table (KaryoScope TSV or
any table with chromosome / length / centromere index or arm columns) or taken
from the assembly sequence lengths.
"""

from __future__ import annotations

import copy
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
SIGNAL_KIND_PREFIX = "FISH "  # a signal placed by hand on the idiogram: kind "FISH <probe>"


def manual_probe(kind: str) -> str:
    """Probe of a signal placed by hand ('FISH CroSat2' -> 'CroSat2'), else ''."""
    return kind[len(SIGNAL_KIND_PREFIX):] if kind.startswith(SIGNAL_KIND_PREFIX) else ""
# How the homologues of one chromosome are named: style -> (menu text, stored suffixes)
HOMOLOGUE_STYLES = {
    "letters": ("a, b, c …", "abcdefghijkl"),
    "capitals": ("A, B, C …", "ABCDEFGHIJKL"),
    "numbers": ("1, 2, 3 …", ",".join(f".{i}" for i in range(1, 13))),
}
GROUP_NAMES = {1: "univalent", 2: "bivalent", 3: "trivalent", 4: "tetravalent", 5: "pentavalent",
               6: "hexavalent", 7: "heptavalent", 8: "octavalent"}


def group_name(copies: int) -> str:
    """'bivalent', 'trivalent' … for the number of homologues of one chromosome."""
    return GROUP_NAMES.get(int(copies), f"{int(copies)} homologues")


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
    color: str = ""  # own colour; empty = the colour of the kind
    satellite: bool = False  # NOR drawn as a secondary constriction with the satellite beyond it

    def __post_init__(self) -> None:
        self.start = max(0.0, min(1.0, float(self.start)))
        self.end = max(0.0, min(1.0, float(self.end)))
        if self.start > self.end:
            self.start, self.end = self.end, self.start


@dataclass
class Homologue:
    """How one homologue differs from its chromosome. None = same as the chromosome."""

    length: float | None = None
    centromere_index: float | None = None
    bands: list[Band] | None = None
    missing: bool = False

    @property
    def separate(self) -> bool:
        return self.length is not None or self.centromere_index is not None or self.bands is not None


@dataclass
class KaryoChromosome:
    name: str
    length: float = 0.0
    centromere_index: float = 50.0  # percent, short arm / total
    copies: int | None = None  # None = the karyotype's ploidy
    assembly_seq: str = ""
    bands: list[Band] = field(default_factory=list)
    homologues: dict[int, Homologue] = field(default_factory=dict)  # index of the copy -> its own values
    b_chromosome: bool = False

    @property
    def type(self) -> str:
        return levan_type(self.centromere_index)

    @property
    def centromere_fraction(self) -> float:
        """Centromere position from the p telomere, as a fraction of the length."""
        return max(0.0, min(1.0, self.centromere_index / 100.0))

    @property
    def arms(self) -> tuple[float, float]:
        """(short arm p, long arm q) in the length unit."""
        p = self.length * self.centromere_fraction
        return p, self.length - p

    def set_arms(self, p: float, q: float) -> None:
        """Set length and centromere index from the two arm lengths (p is always the shorter one)."""
        p, q = sorted((max(0.0, float(p)), max(0.0, float(q))))
        total = p + q
        if total <= 0:
            return
        self.length = round(total, 4)
        self.centromere_index = round(100.0 * p / total, 2)


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
        """Homologue slots of a chromosome (missing ones included)."""
        return self.ploidy if chrom.copies is None else max(0, int(chrom.copies))

    def present_copies(self, chrom: KaryoChromosome) -> int:
        n = self.copies_of(chrom)
        return sum(1 for i in range(n) if not (chrom.homologues.get(i) and chrom.homologues[i].missing))

    def a_chromosomes(self) -> list[KaryoChromosome]:
        return [c for c in self.chromosomes if not c.b_chromosome]

    def b_count(self) -> int:
        return sum(self.present_copies(c) for c in self.chromosomes if c.b_chromosome)

    @property
    def total_count(self) -> int:
        return sum(self.present_copies(c) for c in self.chromosomes)

    @property
    def formula(self) -> str:
        base = self.ploidy * self.basic_number
        a_total = sum(self.present_copies(c) for c in self.a_chromosomes()) if self.chromosomes else base
        head = f"2n = {self.ploidy}x"
        diff = a_total - base
        if diff:
            head += f" {'+' if diff > 0 else '-'} {abs(diff)}"
        nb = self.b_count()
        if nb:
            head += f" + {nb}B"
        return f"{head} = {a_total + nb}"

    def homologue_suffixes(self) -> list[str]:
        """Suffixes added to the chromosome name: a, b, c or .1, .2 or the user's own comma separated list."""
        raw = self.homologue_letters or HOMOLOGUE_STYLES["letters"][1]
        if "," in raw:
            return [s.strip() for s in raw.split(",") if s.strip()]
        return list(raw)

    @property
    def homologue_style(self) -> str:
        for key, (_text, stored) in HOMOLOGUE_STYLES.items():
            if self.homologue_letters == stored:
                return key
        return "custom"

    def homologue_label(self, chrom: KaryoChromosome, i: int) -> str:
        suffixes = self.homologue_suffixes()
        if self.copies_of(chrom) == 1 and chrom.b_chromosome:
            return chrom.name
        return f"{chrom.name}{suffixes[i] if i < len(suffixes) else f'.{i + 1}'}"

    def homologue_labels(self, chrom: KaryoChromosome) -> list[str]:
        """Names of the homologues that are present (a missing one is left out, the others keep their names)."""
        return [self.homologue_label(chrom, i) for i in range(self.copies_of(chrom))
                if not (chrom.homologues.get(i) and chrom.homologues[i].missing)]

    def all_labels(self) -> list[str]:
        return [label for c in self.chromosomes for label in self.homologue_labels(c)]

    def chromosome(self, name: str) -> KaryoChromosome | None:
        return next((c for c in self.chromosomes if c.name == name), None)

    def set_basic_number(self, n: int) -> None:
        """Add (named by number) or drop chromosomes at the end so there are ``n``."""
        n = max(1, int(n))
        self.basic_number = n
        a = [c for c in self.chromosomes if not c.b_chromosome]
        b = [c for c in self.chromosomes if c.b_chromosome]
        while len(a) < n:
            used = {c.name for c in self.chromosomes} | {c.name for c in a}
            k = len(a) + 1
            while str(k) in used:
                k += 1
            a.append(KaryoChromosome(name=str(k)))
        del a[n:]
        self.chromosomes[:] = a + b

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
                band = _band(b)
                if band is not None:
                    bands.append(band)
            copies = item.get("copies")
            homs: dict[int, Homologue] = {}
            for key, h in (item.get("homologues") or {}).items():
                try:
                    hb = h.get("bands")
                    homs[int(key)] = Homologue(
                        length=None if h.get("length") is None else _num(h.get("length")),
                        centromere_index=None if h.get("centromere_index") is None else _num(h.get("centromere_index"), 50.0),
                        bands=None if hb is None else [_band(b) for b in hb if _band(b) is not None],
                        missing=bool(h.get("missing")),
                    )
                except (TypeError, ValueError, AttributeError):
                    continue
            chroms.append(KaryoChromosome(
                name=str(item["name"]).strip(),
                length=_num(item.get("length")),
                centromere_index=_num(item.get("centromere_index"), 50.0),
                copies=int(copies) if copies not in (None, "") else None,
                assembly_seq=str(item.get("assembly_seq") or ""),
                bands=bands,
                homologues=homs,
                b_chromosome=bool(item.get("b_chromosome")),
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


def _band(raw) -> Band | None:
    try:
        return Band(str(raw.get("kind", "")), float(raw["start"]), float(raw["end"]), str(raw.get("color") or ""),
                    bool(raw.get("satellite")))
    except (KeyError, TypeError, ValueError, AttributeError):
        return None


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
    """'2b' -> '2', 'Chr5a' -> 'Chr5', '3.2' -> '3', 'A' -> 'A' (homologue suffix after a digit is dropped)."""
    label = label.strip()
    m = re.match(r"^(.*\d)[._-]\d+$", label)
    if m:
        return m.group(1)
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


# ---------------------------------------------------------------------------
# Chromosome sizes: sizes file, .fai, BED, GFF or FASTA
# ---------------------------------------------------------------------------

SIZES_SUFFIXES = (".sizes", ".genome", ".txt", ".tsv", ".fai", ".bed", ".gff", ".gff3", ".gtf")
SIZES_FILTER = ("Sizes, index, BED, GFF or FASTA (*.sizes *.genome *.txt *.tsv *.fai *.bed *.gff *.gff3 *.gtf "
                "*.fa *.fasta *.fna *.fas *.fa.gz *.fasta.gz *.fna.gz);;All files (*)")
SIZES_OUTPUT_SUFFIX = "_sizes.txt"


def is_fasta(path: Path) -> bool:
    name = Path(path).name.lower()
    if name.endswith(".gz"):
        name = name[:-3]
    return Path(name).suffix in FASTA_SUFFIXES


def sizes_path_for(fasta: Path) -> Path:
    """'genome.fa' -> 'genome.fa_sizes.txt' beside it (full input name kept, then the suffix)."""
    fasta = Path(fasta)
    return fasta.with_name(fasta.name + SIZES_OUTPUT_SUFFIX)


def _newer(a: Path, b: Path) -> bool:
    try:
        return a.is_file() and a.stat().st_mtime >= b.stat().st_mtime
    except OSError:
        return False


def fasta_sizes(path: Path, progress=None, cancel=None) -> dict[str, int]:
    """Record name -> length, read line by line (a whole genome never sits in memory).

    ``progress(percent)`` is called now and then; ``cancel()`` returning True stops with an empty result.
    """
    import gzip

    path = Path(path)
    total = max(1, path.stat().st_size)
    opener = gzip.open if path.name.lower().endswith(".gz") else open
    sizes: dict[str, int] = {}
    name, n, done, next_report = None, 0, 0, total // 50
    with opener(path, "rb") as fh:
        for line in fh:
            done += len(line)
            if line.startswith(b">"):
                if name is not None:
                    sizes[name] = n
                name = line[1:].split(None, 1)[0].decode("utf-8", "replace") if line[1:].strip() else f"seq{len(sizes) + 1}"
                n = 0
            else:
                n += len(line.rstrip(b"\r\n"))
            if done >= next_report:
                next_report += total // 50 or 1
                if cancel is not None and cancel():
                    return {}
                if progress is not None:
                    progress(min(99, int(100 * done / total)))
    if name is not None:
        sizes[name] = n
    return sizes


def _two_columns(text: str) -> dict[str, int]:
    sizes: dict[str, int] = {}
    for line in text.splitlines():
        parts = line.strip().split("\t") if "\t" in line else line.split()
        if len(parts) < 2 or line.startswith(("#", "track", "browser")):
            continue
        try:
            sizes.setdefault(parts[0], int(float(parts[1])))
        except ValueError:
            continue  # header line
    return sizes


def _last_end(text: str, name_col: int, end_col: int) -> dict[str, int]:
    sizes: dict[str, int] = {}
    regions: dict[str, int] = {}
    for line in text.splitlines():
        if line.startswith("##sequence-region"):
            parts = line.split()
            if len(parts) >= 4 and parts[3].isdigit():
                regions[parts[1]] = int(parts[3])
            continue
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        parts = line.split("\t")
        if len(parts) <= max(name_col, end_col):
            continue
        try:
            end = int(float(parts[end_col]))
        except ValueError:
            continue
        sizes[parts[name_col]] = max(end, sizes.get(parts[name_col], 0))
    sizes.update(regions)
    return sizes


def read_sizes(path: Path, progress=None, cancel=None) -> dict[str, int]:
    """Sequence name -> length in bp from a sizes file, .fai, BED (last end), GFF/GTF or FASTA.

    For FASTA a newer ``.fai`` or ``<name>_sizes.txt`` beside it is used when present.
    """
    path = Path(path).expanduser()
    if is_fasta(path):
        for side in (Path(str(path) + ".fai"), sizes_path_for(path)):
            if _newer(side, path):
                return _two_columns(side.read_text(encoding="utf-8", errors="replace"))
        return fasta_sizes(path, progress, cancel)
    text = path.read_text(encoding="utf-8", errors="replace")
    suffix = path.suffix.lower()
    if suffix == ".bed":
        return _last_end(text, 0, 2)
    if suffix in (".gff", ".gff3", ".gtf"):
        return _last_end(text, 0, 4)
    return _two_columns(text)


def write_sizes(path: Path, sizes: dict[str, int]) -> Path:
    path = Path(path)
    path.write_text("".join(f"{name}\t{n}\n" for name, n in sizes.items()), encoding="utf-8")
    return path


def make_sizes_file(fasta: Path, out: Path | None = None, progress=None, cancel=None) -> Path:
    """Write ``<name>_sizes.txt`` (name and length per sequence, tab separated) for a FASTA."""
    fasta = Path(fasta).expanduser().resolve()
    sizes = fasta_sizes(fasta, progress, cancel)
    if not sizes:
        raise ValueError(f"No sequences in {fasta.name}.")
    return write_sizes(out or sizes_path_for(fasta), sizes)


def sizes_sources(folder: Path) -> list[Path]:
    """Files in ``folder`` that can give chromosome lengths: sizes files and indexes first, then BED/GFF, then FASTA."""
    def rank(p: Path) -> int:
        name = p.name.lower()
        if name.endswith((SIZES_OUTPUT_SUFFIX, ".sizes", ".genome", ".fai")):
            return 0
        if is_fasta(p):
            return 2
        return 1

    try:
        files = [p for p in Path(folder).iterdir() if p.is_file() and not p.name.startswith(".")
                 and (is_fasta(p) or p.suffix.lower() in (".sizes", ".genome", ".fai", ".bed", ".gff", ".gff3", ".gtf")
                      or p.name.lower().endswith(SIZES_OUTPUT_SUFFIX))]
    except OSError:
        return []
    return sorted(files, key=lambda p: (rank(p), p.name.lower()))


def match_sequences(karyotype: Karyotype, names: list[str]) -> dict[str, str]:
    """Chromosome name -> sequence name: the current link, the same name, or Chr1 / chr01 / Csat_chr1 for '1'."""
    known = set(names)
    by_number: dict[str, str] = {}
    for seq in names:
        m = re.search(r"(?:chr(?:omosome)?)[_\- ]?0*(\d+)$", seq, re.IGNORECASE)
        if m and m.group(1) not in by_number:
            by_number[m.group(1)] = seq
    out: dict[str, str] = {}
    for chrom in karyotype.chromosomes:
        if chrom.assembly_seq in known:
            out[chrom.name] = chrom.assembly_seq
        elif chrom.name in known:
            out[chrom.name] = chrom.name
        else:
            seq = by_number.get(chrom.name.lstrip("0") or chrom.name)
            if seq:
                out[chrom.name] = seq
    return out


def apply_sizes(karyotype: Karyotype, links: dict[str, str], sizes: dict[str, int]) -> list[str]:
    """Link the sequences and set lengths in Mb. Returns chromosomes left without a length."""
    missing = []
    for chrom in karyotype.chromosomes:
        seq = links.get(chrom.name, "")
        chrom.assembly_seq = seq
        if seq in sizes:
            chrom.length = round(sizes[seq] / 1e6, 3)
        else:
            missing.append(chrom.name)
    karyotype.unit = "Mb"
    return missing


# ---------------------------------------------------------------------------
# Homologues that differ from their chromosome
# ---------------------------------------------------------------------------


def homologue_view(chrom: KaryoChromosome, i: int) -> KaryoChromosome:
    """Chromosome as homologue ``i`` looks (its own length, centromere and bands where it has them)."""
    h = chrom.homologues.get(i)
    if h is None or not h.separate:
        return chrom
    return KaryoChromosome(
        name=chrom.name,
        length=chrom.length if h.length is None else h.length,
        centromere_index=chrom.centromere_index if h.centromere_index is None else h.centromere_index,
        copies=chrom.copies,
        assembly_seq=chrom.assembly_seq,
        bands=chrom.bands if h.bands is None else h.bands,
        b_chromosome=chrom.b_chromosome,
    )


def is_separate(chrom: KaryoChromosome, i: int) -> bool:
    h = chrom.homologues.get(i)
    return bool(h and h.separate)


def is_missing(chrom: KaryoChromosome, i: int) -> bool:
    h = chrom.homologues.get(i)
    return bool(h and h.missing)


def differs(chrom: KaryoChromosome) -> bool:
    """True when any homologue has its own values or is missing."""
    return any(h.separate or h.missing for h in chrom.homologues.values())


def _tidy(chrom: KaryoChromosome, i: int) -> None:
    h = chrom.homologues.get(i)
    if h is not None and not h.separate and not h.missing:
        del chrom.homologues[i]


def set_separate(chrom: KaryoChromosome, i: int, separate: bool) -> None:
    """Give homologue ``i`` its own copy of the chromosome's values, or make it the same as the others again."""
    h = chrom.homologues.setdefault(i, Homologue())
    if separate and not h.separate:
        h.length = chrom.length
        h.centromere_index = chrom.centromere_index
        h.bands = [copy.deepcopy(b) for b in chrom.bands]
    elif not separate:
        h.length = h.centromere_index = None
        h.bands = None
    _tidy(chrom, i)


def set_missing(chrom: KaryoChromosome, i: int, missing: bool) -> None:
    chrom.homologues.setdefault(i, Homologue()).missing = bool(missing)
    _tidy(chrom, i)


def target(chrom: KaryoChromosome, i: int):
    """Object to edit for homologue ``i``: its own Homologue when separate, else the chromosome."""
    return chrom.homologues[i] if i >= 0 and is_separate(chrom, i) else chrom


def drop_extra_homologues(karyotype: Karyotype, chrom: KaryoChromosome) -> None:
    n = karyotype.copies_of(chrom)
    for i in [i for i in chrom.homologues if i >= n]:
        del chrom.homologues[i]


def mirror_bands(chrom: KaryoChromosome, i: int) -> None:
    """Turn the band pattern of homologue ``i`` end to end (an inversion); the homologue becomes separate."""
    if i >= 0:
        set_separate(chrom, i, True)
    owner = target(chrom, i)
    owner.bands = [Band(b.kind, 1.0 - b.end, 1.0 - b.start, b.color, b.satellite) for b in owner.bands]


def copy_bands_to_others(karyotype: Karyotype, chrom: KaryoChromosome, i: int) -> None:
    """Every homologue (and the chromosome) gets the bands of homologue ``i``."""
    bands = homologue_view(chrom, i).bands
    chrom.bands = [copy.deepcopy(b) for b in bands]
    for j, h in list(chrom.homologues.items()):
        if j != i and h.bands is not None:
            h.bands = [copy.deepcopy(b) for b in bands]
    if i in chrom.homologues and chrom.homologues[i].bands is not None:
        chrom.homologues[i].bands = [copy.deepcopy(b) for b in bands]
        h = chrom.homologues[i]
        if (h.length in (None, chrom.length)) and (h.centromere_index in (None, chrom.centromere_index)):
            set_separate(chrom, i, False)


def add_b_chromosome(karyotype: Karyotype, chrom: KaryoChromosome, i: int) -> KaryoChromosome:
    """New B chromosome (one copy) shaped like homologue ``i`` of ``chrom``, added at the end."""
    view = homologue_view(chrom, max(0, i))
    used = {c.name for c in karyotype.chromosomes}
    n = 1
    while f"B{n}" in used:
        n += 1
    new = KaryoChromosome(name=f"B{n}", length=view.length, centromere_index=view.centromere_index, copies=1,
                          bands=[copy.deepcopy(b) for b in view.bands], b_chromosome=True)
    karyotype.chromosomes.append(new)
    return new


# ---------------------------------------------------------------------------
# Band presets per species (band_presets.json in the karyotype store)
# ---------------------------------------------------------------------------


def presets_path(folder: Path) -> Path:
    return karyotype_path(folder, "x").parent / "band_presets.json"


def _load_presets(folder: Path) -> dict:
    try:
        raw = json.loads(presets_path(folder).read_text(encoding="utf-8"))
        return raw if isinstance(raw, dict) else {}
    except (OSError, ValueError):
        return {}


def list_presets(folder: Path, species: str) -> list[str]:
    return sorted((_load_presets(folder).get(species.strip() or "Any species") or {}).keys(), key=str.lower)


def save_preset(folder: Path, species: str, name: str, karyotype: Karyotype) -> Path:
    """Store the bands of every chromosome under ``name`` for ``species``."""
    data = _load_presets(folder)
    bands = {c.name: [asdict(b) for b in c.bands] for c in karyotype.chromosomes if c.bands}
    data.setdefault(species.strip() or "Any species", {})[name.strip()] = bands
    path = presets_path(folder)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    os.replace(tmp, path)
    return path


def apply_preset(folder: Path, species: str, name: str, karyotype: Karyotype) -> list[str]:
    """Replace the bands of the chromosomes named in the preset. Returns preset chromosomes not in the karyotype."""
    preset = (_load_presets(folder).get(species.strip() or "Any species") or {}).get(name) or {}
    missing = []
    for cname, bands in preset.items():
        chrom = karyotype.chromosome(cname)
        if chrom is None:
            missing.append(cname)
            continue
        chrom.bands = [b for b in (_band(x) for x in bands) if b is not None]
    return missing
