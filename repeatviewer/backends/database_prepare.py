"""Validate user databases and build BLAST / CDD / HMM indexes (TE-Aid-style)."""

from __future__ import annotations

import logging
import os
import platform
import subprocess
from repeatviewer.external_process import run_logged_process
import tempfile
import traceback
from dataclasses import dataclass
from pathlib import Path

from repeatviewer.backends.gui_functions import (
    check_cdd_index_files,
    check_cmd_in_path,
    get_original_file_path,
    sanitize_name,
)
from repeatviewer.project import DATABASE_TYPE_DNA, BlastDatabase, UserDatabase

logger = logging.getLogger(__name__)


def _run_logged(cmd, *, inputs=(), **kwargs):
    tool = Path(cmd[0]).name
    return run_logged_process(
        cmd,
        tool=tool,
        inputs=inputs,
        cwd=kwargs.get("cwd"),
        check=bool(kwargs.get("check", False)),
    )

# Full CDD marker (preferred) plus NCBI subset / LE set stems.
CDD_MARKER = "Cdd.pn"
CDD_PROFILE_PREFIX = "cdd_profile"
CDD_SET_STEMS = (
    "Cdd",
    "Cdd_NCBI",
    "Pfam",
    "Smart",
    "Cog",
    "Kog",
    "Prk",
    "Tigr",
)
CDD_PN_NAMES = tuple(f"{stem}.pn" for stem in CDD_SET_STEMS)
CDD_PAL_NAMES = tuple(f"{stem}.pal" for stem in CDD_SET_STEMS)
CDD_INDEX_MINUTES = 20


@dataclass(frozen=True)
class DatabaseTarget:
    label: str
    kind: str  # cdd | protein | dna | hmm
    path: Path

    @property
    def needs_index(self) -> bool:
        return not prepare_index_is_ready(self.kind, self.path)


def _os_type() -> str:
    return platform.system()


def _blast_dir() -> Path:
    return Path(get_original_file_path()) / "blast"


def _bundled_or_system(exe: str, *, use_system_blast: bool) -> str:
    if use_system_blast and check_cmd_in_path(exe):
        return exe
    os_type = _os_type()
    suffix = "_linux" if os_type == "Linux" else "_mac" if os_type == "Darwin" else ".exe"
    bundled = _blast_dir() / f"{exe}{suffix}"
    if bundled.is_file():
        return str(bundled)
    if check_cmd_in_path(exe):
        return exe
    raise FileNotFoundError(
        f"'{exe}' was not found. Install BLAST+ or enable bundled BLAST in Preferences."
    )


def _cdd_child_by_names(cdd_dir: Path, names: tuple[str, ...]) -> Path | None:
    """Return the first existing file among ``names``, preferring order.

    Probes known filenames directly (no ``iterdir``) so large preformatted LE
    folders stay fast when Manage Database switches tabs.
    """
    cdd_dir = cdd_dir.expanduser().resolve()
    if not cdd_dir.is_dir():
        return None
    for name in names:
        candidate = cdd_dir / name
        try:
            if candidate.is_file():
                return candidate
        except OSError:
            continue
    return None


def find_cdd_pn(cdd_dir: Path) -> Path | None:
    """Locate a CDD / subset profile list (``Cdd.pn``, ``Pfam.pn``, …)."""
    return _cdd_child_by_names(cdd_dir, CDD_PN_NAMES)


def find_cdd_preformatted_prefix(cdd_dir: Path) -> Path | None:
    """Locate a preformatted NCBI LE RPS-BLAST alias (``Cdd.pal``, ``Pfam.pal``, …).

    Returns the ``-db`` prefix path (stem without ``.pal``).
    """
    pal = _cdd_child_by_names(cdd_dir, CDD_PAL_NAMES)
    if pal is None:
        return None
    return pal.with_suffix("")


def suggest_cdd_database_name(cdd_dir: Path) -> str:
    """Display name from NCBI set stem (Pfam, Cdd, …) or the folder name."""
    cdd_dir = cdd_dir.expanduser().resolve()
    pn = find_cdd_pn(cdd_dir)
    if pn is not None:
        return pn.stem
    pref = find_cdd_preformatted_prefix(cdd_dir)
    if pref is not None:
        return pref.name
    return cdd_dir.name or ""


def cdd_index_is_ready(cdd_dir: Path) -> bool:
    """True when RepeatViewer-built ``cdd_profile`` or a preformatted LE set is present."""
    cdd_dir = cdd_dir.expanduser().resolve()
    if not cdd_dir.is_dir():
        return False
    if find_cdd_preformatted_prefix(cdd_dir) is not None:
        return True
    return check_cdd_index_files(str(cdd_dir))


def cdd_database_present(cdd_dir: Path) -> bool:
    """True for source (``.pn``) or preformatted LE (``.pal``) CDD / subset folders."""
    cdd_dir = cdd_dir.expanduser().resolve()
    if not cdd_dir.is_dir():
        return False
    return find_cdd_pn(cdd_dir) is not None or find_cdd_preformatted_prefix(cdd_dir) is not None


def blast_index_file(seq_path: Path, *, protein: bool) -> Path:
    suffix = ".pin" if protein else ".nin"
    return Path(f"{seq_path.expanduser().resolve()}{suffix}")


def annotation_blast_prefix_for(seq_path: Path) -> Path:
    """Central index prefix under Documents/RepeatViewer/blast_databases/<name>/<name>.

    Annotation DNA/Protein prep uses this instead of writing .nin/.pin beside the FASTA.
    """
    from repeatviewer.app_preferences import load_blast_db_base_dir
    from repeatviewer.backends.gui_functions import sanitize_name

    safe = sanitize_name(Path(seq_path).stem) or "database"
    return load_blast_db_base_dir().expanduser().resolve() / safe / safe


def blast_search_db_path(seq_path: Path, *, protein: bool) -> Path:
    """Path passed to BLAST ``-db``: prefer central named index, else legacy beside FASTA."""
    seq_path = seq_path.expanduser().resolve()
    named = annotation_blast_prefix_for(seq_path)
    if named_blast_index_exists(named, protein=protein):
        return named
    return seq_path


def _fasta_blast_index_ready(seq_path: Path, *, protein: bool) -> bool:
    seq_path = seq_path.expanduser().resolve()
    if not seq_path.is_file():
        return False
    named = annotation_blast_prefix_for(seq_path)
    if named_blast_index_exists(named, protein=protein):
        return True
    return blast_index_file(seq_path, protein=protein).is_file()


def hmm_index_files(hmm_path: Path) -> tuple[Path, Path, Path, Path]:
    # hmmpress appends .h3m/.h3f/.h3i/.h3p to the full profile path (e.g. Pfam-A.hmm.h3m).
    base = str(hmm_path.expanduser().resolve())
    return (
        Path(f"{base}.h3m"),
        Path(f"{base}.h3f"),
        Path(f"{base}.h3i"),
        Path(f"{base}.h3p"),
    )


def _pyhmmer_importable() -> bool:
    try:
        import pyhmmer  # type: ignore[import-untyped]
        from pyhmmer.plan7 import HMMFile  # type: ignore[import-untyped]  # noqa: F401
    except Exception:
        return False
    return getattr(pyhmmer, "__file__", None) is not None


def pyhmmer_install_hint() -> str:
    if _os_type() == "Windows":
        return 'pip install "pyhmmer>=0.12.0" (64-bit Python required)'
    return "pip install pyhmmer"


def hmm_indexing_available() -> bool:
    return check_cmd_in_path("hmmpress") or _pyhmmer_importable()


def database_search_readiness(db: UserDatabase) -> str | None:
    """Return a user-facing reason when a database cannot be searched yet."""
    path = db.path.expanduser().resolve()
    kind = db.prepare_kind
    problem = describe_database_problem(kind, path)
    if problem:
        return problem.replace("\n", " ")
    if kind == "protein":
        if not _fasta_blast_index_ready(path, protein=True):
            return "Protein database is not indexed. Use File → Manage Database."
    elif kind == "dna":
        if not _fasta_blast_index_ready(path, protein=False):
            return "DNA database is not indexed. Use File → Manage Database."
    elif kind == "cdd":
        if not cdd_index_is_ready(path):
            return "CDD database is not indexed. Use File → Manage Database."
    elif kind == "hmm":
        if not _hmm_prepare_ready(path):
            if not hmm_indexing_available():
                return (
                    "HMM database is not prepared. Install hmmpress or "
                    f"pyhmmer ({pyhmmer_install_hint()}), then use File → Manage Database."
                )
            return (
                "HMM database is not prepared. Use File → Manage Database to index it "
                "(required for large databases such as Dfam)."
            )
    return None


def _hmm_pressed(hmm_path: Path) -> bool:
    return all(p.is_file() for p in hmm_index_files(hmm_path))


def _hmm_prepare_marker(hmm_path: Path) -> Path:
    return Path(f"{hmm_path}.repeatviewer_hmm_ok")


def _mark_hmm_prepared(hmm_path: Path) -> None:
    _hmm_prepare_marker(hmm_path).write_text("ok\n", encoding="utf-8")


def _hmm_prepare_ready(hmm_path: Path) -> bool:
    hmm_path = hmm_path.expanduser().resolve()
    if not hmm_path.is_file():
        return False
    if _hmm_pressed(hmm_path):
        return True
    marker = _hmm_prepare_marker(hmm_path)
    if not marker.is_file():
        return False
    try:
        return marker.stat().st_mtime >= hmm_path.stat().st_mtime
    except OSError:
        return False


def prepare_index_is_ready(kind: str, path: Path) -> bool:
    """Whether database preparation can be skipped (Define database → OK)."""
    path = path.expanduser().resolve()
    if kind == "cdd":
        return path.is_dir() and cdd_index_is_ready(path)
    if kind == "protein":
        return _fasta_blast_index_ready(path, protein=True)
    if kind == "dna":
        return _fasta_blast_index_ready(path, protein=False)
    if kind == "hmm":
        return _hmm_prepare_ready(path)
    return False


def index_is_ready(kind: str, path: Path) -> bool:
    path = path.expanduser().resolve()
    if kind == "cdd":
        return path.is_dir() and cdd_index_is_ready(path)
    if kind == "protein":
        return _fasta_blast_index_ready(path, protein=True)
    if kind == "dna":
        return _fasta_blast_index_ready(path, protein=False)
    if kind == "hmm":
        if not path.is_file():
            return False
        if _hmm_prepare_ready(path):
            return True
        return _pyhmmer_importable()
    return False


def describe_database_problem(kind: str, path: Path) -> str | None:
    path = path.expanduser().resolve()
    if kind == "cdd":
        if not path.is_dir():
            return f"CDD directory not found:\n{path}"
        if not cdd_database_present(path):
            markers = ", ".join(CDD_PN_NAMES[:3]) + ", …"
            pals = ", ".join(CDD_PAL_NAMES[:3]) + ", …"
            return (
                f"CDD database files were not found in:\n{path}\n\n"
                f"Expected a source list ({markers}) or a preformatted "
                f"NCBI LE alias ({pals})."
            )
        return None
    if not path.is_file():
        return f"Database file not found:\n{path}"
    return None


def collect_prepare_targets(
    *,
    cdd_dir: Path | None,
    user_databases: list[UserDatabase],
) -> list[DatabaseTarget]:
    targets: list[DatabaseTarget] = []
    if cdd_dir is not None:
        resolved = cdd_dir.expanduser().resolve()
        if str(resolved):
            targets.append(DatabaseTarget("CDD database", "cdd", resolved))
    for entry in user_databases:
        if not entry.name or not entry.path:
            continue
        kind = entry.prepare_kind
        targets.append(
            DatabaseTarget(
                f"{entry.name} ({entry.display_type_label})",
                kind,
                entry.path.expanduser().resolve(),
            )
        )
    return targets


def _kind_for_db_type(db_type: str, *, is_hmm: bool = False) -> str:
    if is_hmm:
        return "hmm"
    if db_type == DATABASE_TYPE_DNA:
        return "dna"
    return "protein"


def build_cdd_profile_index(
    cdd_dir: Path,
    *,
    use_system_blast: bool,
) -> Path:
    cdd_dir = cdd_dir.expanduser().resolve()
    problem = describe_database_problem("cdd", cdd_dir)
    if problem:
        raise FileNotFoundError(problem)

    # Preformatted NCBI LE packages (Pfam_LE, Cdd_LE, …) are already searchable.
    preformatted = find_cdd_preformatted_prefix(cdd_dir)
    if preformatted is not None:
        return preformatted

    out_prefix = cdd_dir / CDD_PROFILE_PREFIX
    if check_cdd_index_files(str(cdd_dir)):
        return out_prefix

    cdd_pn = find_cdd_pn(cdd_dir)
    if cdd_pn is None:
        raise FileNotFoundError(
            f"CDD profile list (.pn) was not found in:\n{cdd_dir}"
        )
    makeprofiledb = _bundled_or_system("makeprofiledb", use_system_blast=use_system_blast)
    cmd = [
        makeprofiledb,
        "-in",
        str(cdd_pn),
        "-out",
        str(out_prefix),
        "-threshold",
        "9.82",
        "-scale",
        "100.0",
        "-dbtype",
        "rps",
    ]
    logger.info("Building CDD profile index: %s", " ".join(cmd))
    _run_logged(
        cmd,
        inputs=[cdd_pn],
        check=True,
        cwd=str(cdd_dir),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if not check_cdd_index_files(str(cdd_dir)):
        raise RuntimeError("CDD profile index generation finished but index files are missing.")
    return out_prefix


def build_blast_index(
    seq_path: Path,
    *,
    protein: bool,
    use_system_blast: bool,
) -> Path:
    """Build a BLAST index under the central blast_databases folder (not beside the FASTA).

    Returns the ``-db`` prefix path. Legacy in-place indexes beside the FASTA are
    still recognized by :func:`blast_search_db_path` until rebuilt.
    """
    from repeatviewer.app_preferences import load_blast_db_base_dir

    seq_path = seq_path.expanduser().resolve()
    kind = "protein" if protein else "dna"
    problem = describe_database_problem(kind, seq_path)
    if problem:
        raise FileNotFoundError(problem)
    if _fasta_blast_index_ready(seq_path, protein=protein):
        return blast_search_db_path(seq_path, protein=protein)

    molecule_type = "prot" if protein else "nucl"
    return build_named_blast_database(
        seq_path,
        load_blast_db_base_dir(),
        seq_path.stem,
        molecule_type=molecule_type,
        use_system_blast=use_system_blast,
    )


def _named_blast_index_marker(db_prefix: Path, *, protein: bool) -> Path:
    """Index file that proves ``db_prefix`` was built.

    ``.nin`` / ``.pin`` for a single-volume database. makeblastdb splits a
    large genome into volumes (``prefix.00.nin`` ...) and writes an alias
    file ``.nal`` / ``.pal`` instead, so fall back to that.
    """
    suffix = ".pin" if protein else ".nin"
    single = Path(f"{db_prefix}{suffix}")
    if single.is_file():
        return single
    alias = Path(f"{db_prefix}{'.pal' if protein else '.nal'}")
    return alias if alias.is_file() else single


def named_blast_index_exists(db_prefix: Path, *, protein: bool) -> bool:
    return _named_blast_index_marker(db_prefix, protein=protein).is_file()


def _named_blast_index_is_current(db_prefix: Path, source_fasta: Path, *, protein: bool) -> bool:
    marker = _named_blast_index_marker(db_prefix, protein=protein)
    if not marker.is_file():
        return False
    try:
        return marker.stat().st_mtime >= source_fasta.stat().st_mtime
    except OSError:
        return False


def blast_database_is_ready(db: BlastDatabase) -> bool:
    """Whether a saved BLAST database's index files are still present on disk."""
    return named_blast_index_exists(db.db_prefix_path, protein=db.molecule_type == "prot")


def _same_file(a: Path, b: Path) -> bool:
    try:
        return Path(a).expanduser().resolve() == Path(b).expanduser().resolve()
    except OSError:
        return False


def saved_blast_database_for(
    fasta: Path, databases: list[BlastDatabase]
) -> BlastDatabase | None:
    """The registered nucleotide database built from ``fasta``, if still usable.

    Usable means its index files exist and are not older than the FASTA, so an
    edited genome is not searched against a stale index.
    """
    fasta = Path(fasta).expanduser()
    for db in databases:
        if db.molecule_type == "prot":
            continue
        if not _same_file(db.source_fasta, fasta):
            continue
        if _named_blast_index_is_current(db.db_prefix_path, fasta, protein=False):
            return db
    return None


def unique_blast_database_name(
    base_name: str, databases: list[BlastDatabase], base_dir: Path
) -> str:
    """``base_name``, or ``base_name (2)``, ``(3)``... when the name or folder is taken.

    A taken folder is never reused: it may hold an index of a different FASTA
    with the same file name, which the mtime check alone cannot tell apart.
    """
    taken = {db.name.lower() for db in databases}
    root = Path(base_dir).expanduser()
    base_name = base_name.strip() or "database"
    candidate, n = base_name, 1
    while (
        candidate.lower() in taken
        or (root / (sanitize_name(candidate) or "database")).exists()
    ):
        n += 1
        candidate = f"{base_name} ({n})"
    return candidate


def build_named_blast_database(
    source_fasta: Path,
    base_dir: Path,
    name: str,
    *,
    molecule_type: str,
    use_system_blast: bool,
) -> Path:
    """Run makeblastdb for a saved, named BLAST database. Returns the db prefix path.

    Index files are written under ``base_dir / sanitize_name(name)`` — kept
    separate from the source FASTA (which is not copied) so the same folder
    can be safely deleted when the database is removed.
    """
    source_fasta = source_fasta.expanduser().resolve()
    if not source_fasta.is_file():
        raise FileNotFoundError(f"FASTA not found: {source_fasta}")
    protein = molecule_type == "prot"
    safe_name = sanitize_name(name.strip()) or "database"
    dest_dir = base_dir.expanduser().resolve() / safe_name
    db_prefix = dest_dir / safe_name

    if _named_blast_index_is_current(db_prefix, source_fasta, protein=protein):
        return db_prefix

    dest_dir.mkdir(parents=True, exist_ok=True)
    makeblastdb = _bundled_or_system("makeblastdb", use_system_blast=use_system_blast)
    dbtype = "prot" if protein else "nucl"
    cmd = [
        makeblastdb,
        "-in",
        str(source_fasta),
        "-dbtype",
        dbtype,
        "-out",
        str(db_prefix),
    ]
    logger.info("Building named BLAST database: %s", " ".join(cmd))
    try:
        _run_logged(cmd, inputs=[source_fasta], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(format_prepare_error(exc)) from exc
    if not named_blast_index_exists(db_prefix, protein=protein):
        raise RuntimeError(f"BLAST database creation failed for {source_fasta.name}.")
    return db_prefix


def build_hmm_index(hmm_path: Path) -> Path:
    hmm_path = hmm_path.expanduser().resolve()
    problem = describe_database_problem("hmm", hmm_path)
    if problem:
        raise FileNotFoundError(problem)
    if _hmm_pressed(hmm_path):
        return hmm_path

    if not check_cmd_in_path("hmmpress"):
        if _pyhmmer_importable():
            logger.info(
                "hmmpress not found; validating unpressed HMM with pyhmmer: %s", hmm_path
            )
            _run_hmm_test(hmm_path)
            _mark_hmm_prepared(hmm_path)
            return hmm_path
        if _os_type() == "Windows":
            raise FileNotFoundError(
                "hmmpress is not available on Windows. Install pyhmmer 0.12+ "
                '(pip install "pyhmmer>=0.12.0"; 64-bit Python required) to search '
                "unpressed .hmm files, or copy pre-pressed Pfam index files (.h3*) "
                "from a Linux or macOS machine."
            )
        raise FileNotFoundError(
            "hmmpress was not found. Install HMMER (conda install -c bioconda hmmer) "
            "or pyhmmer (pip install pyhmmer) to use HMM databases."
        )
    cmd = ["hmmpress", "-f", str(hmm_path)]
    logger.info("Pressing HMM database: %s", " ".join(cmd))
    try:
        completed = _run_logged(
            cmd, inputs=[hmm_path], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
    except subprocess.CalledProcessError as exc:
        detail = (exc.stderr or exc.stdout or "").strip()
        message = f"hmmpress failed for {hmm_path.name}."
        if detail:
            message = f"{message}\n\n{detail}"
        raise RuntimeError(message) from exc
    if not _hmm_pressed(hmm_path):
        detail = (completed.stderr or completed.stdout or "").strip()
        message = f"HMM press finished but index files are missing for {hmm_path.name}."
        if detail:
            message = f"{message}\n\n{detail}"
        raise RuntimeError(message)
    _mark_hmm_prepared(hmm_path)
    return hmm_path


def _run_blastp_test(db_path: Path, *, use_system_blast: bool) -> None:
    try:
        blastp = _bundled_or_system("blastp", use_system_blast=use_system_blast)
    except FileNotFoundError:
        logger.info("Skipping BLASTP validation (blastp not available in bundled BLAST+).")
        return
    with tempfile.TemporaryDirectory(prefix="repeatviewer_dbcheck_") as tmp:
        query = Path(tmp) / "query.fa"
        query.write_text(">test\nMKLLVVVFCLGVLA\n", encoding="utf-8")
        out = Path(tmp) / "out.txt"
        cmd = [
            blastp,
            "-query",
            str(query),
            "-db",
            str(db_path),
            "-out",
            str(out),
            "-outfmt",
            "6",
            "-max_target_seqs",
            "1",
            "-evalue",
            "10",
        ]
        _run_logged(cmd, inputs=[query, db_path], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)


def _run_blastn_test(db_path: Path, *, use_system_blast: bool) -> None:
    blastn = _bundled_or_system("blastn", use_system_blast=use_system_blast)
    with tempfile.TemporaryDirectory(prefix="repeatviewer_dbcheck_") as tmp:
        query = Path(tmp) / "query.fa"
        query.write_text(">test\nATGCGTACGTAGCTAGCTAG\n", encoding="utf-8")
        out = Path(tmp) / "out.txt"
        cmd = [
            blastn,
            "-query",
            str(query),
            "-db",
            str(db_path),
            "-out",
            str(out),
            "-outfmt",
            "6",
            "-max_target_seqs",
            "1",
            "-evalue",
            "10",
        ]
        _run_logged(cmd, inputs=[query, db_path], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)


def _run_hmm_test(hmm_path: Path) -> None:
    from repeatviewer.backends.pyhmmer_search import validate_hmm_search

    try:
        import pyhmmer  # type: ignore[import-untyped]  # noqa: F401
    except ImportError:
        logger.warning("pyhmmer is not installed; skipped HMM validation search.")
        return

    try:
        validate_hmm_search(hmm_path)
    except Exception as exc:  # noqa: BLE001
        logger.warning("pyhmmer HMM validation search failed: %s", exc)


def validate_database_search(
    target: DatabaseTarget,
    *,
    use_system_blast: bool,
) -> None:
    if target.kind == "protein":
        _run_blastp_test(target.path, use_system_blast=use_system_blast)
    elif target.kind == "dna":
        _run_blastn_test(target.path, use_system_blast=use_system_blast)
    elif target.kind == "hmm":
        _run_hmm_test(target.path)


def prepare_database_target(
    target: DatabaseTarget,
    *,
    use_system_blast: bool,
    validate_search: bool = True,
) -> None:
    problem = describe_database_problem(target.kind, target.path)
    if problem:
        raise FileNotFoundError(problem)

    if target.kind == "cdd":
        build_cdd_profile_index(target.path, use_system_blast=use_system_blast)
    elif target.kind == "protein":
        prefix = build_blast_index(target.path, protein=True, use_system_blast=use_system_blast)
        if validate_search:
            _run_blastp_test(prefix, use_system_blast=use_system_blast)
    elif target.kind == "dna":
        prefix = build_blast_index(target.path, protein=False, use_system_blast=use_system_blast)
        if validate_search:
            _run_blastn_test(prefix, use_system_blast=use_system_blast)
    elif target.kind == "hmm":
        build_hmm_index(target.path)
        if validate_search:
            _run_hmm_test(target.path)


def prepared_cdd_prefix(cdd_dir: Path | None) -> Path | None:
    """Return RPS-BLAST ``-db`` prefix when the CDD / subset folder is ready."""
    if cdd_dir is None:
        return None
    cdd_dir = cdd_dir.expanduser().resolve()
    if not cdd_database_present(cdd_dir):
        return None
    preformatted = find_cdd_preformatted_prefix(cdd_dir)
    if preformatted is not None:
        return preformatted
    if check_cdd_index_files(str(cdd_dir)):
        return cdd_dir / CDD_PROFILE_PREFIX
    return None


def index_build_message(target: DatabaseTarget) -> str:
    if target.kind == "cdd":
        return (
            f"Build CDD profile index for “{target.label}”?\n\n"
            f"Directory:\n{target.path}\n\n"
            f"This usually takes about {CDD_INDEX_MINUTES} minutes."
        )
    if target.kind == "hmm":
        if check_cmd_in_path("hmmpress"):
            detail = (
                "hmmpress will create .hmm.h3m / .hmm.h3f / .hmm.h3i / .hmm.h3p index files."
            )
        elif _pyhmmer_importable():
            detail = (
                "hmmpress is not available; pyhmmer will validate the unpressed .hmm file "
                "(slower to load than a pressed index)."
            )
        else:
            detail = (
                "Install pyhmmer (pip install pyhmmer) or HMMER with hmmpress to prepare "
                "this database."
            )
        return (
            f"Prepare HMM database for “{target.label}”?\n\n"
            f"File:\n{target.path}\n\n"
            f"{detail}"
        )
    blast_mode = "BLASTP (protein)" if target.kind == "protein" else "BLASTN (DNA)"
    return (
        f"Build {blast_mode} index for “{target.label}”?\n\n"
        f"File:\n{target.path}\n\n"
        "makeblastdb will create BLAST index files next to the FASTA."
    )


def format_prepare_error(exc: BaseException) -> str:
    if isinstance(exc, subprocess.CalledProcessError):
        detail = (exc.stderr or exc.stdout or "").strip()
        if detail:
            return f"{exc}\n\n{detail}"
    return str(exc) or traceback.format_exc()
