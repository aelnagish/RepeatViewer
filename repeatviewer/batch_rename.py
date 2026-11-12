"""Edit > Batch Rename: new names and their checks (no Qt, unit tested).

The dialog builds a :class:`RenameOptions` from its fields and shows the
:func:`plan_batch_rename` result as a preview. The main window performs the
plan's moves with :func:`apply_moves`, which renames in two steps (old name to
a temporary name, then to the new name) so swaps and chains such as
a -> b, b -> c work and a failure puts every file back.

Developer: Abdullah Yusuf
"""

from __future__ import annotations

import os
import platform
import re
import uuid
from dataclasses import dataclass, field
from pathlib import Path, PurePath
from typing import Callable, Iterable

# Kept together with the extension before them: x.fa.gz -> ".fa.gz".
_COMPRESSED = {".gz", ".bz2", ".xz", ".zip", ".zst", ".bgz"}
_WINDOWS_BAD = set('<>:"|?*\\')
_MAX_NAME_BYTES = 255


@dataclass(frozen=True)
class RenameOptions:
    find: str = ""
    replace: str = ""
    match_case: bool = False
    regex: bool = False
    prefix: str = ""
    suffix: str = ""
    keep_extension: bool = True

    @property
    def does_something(self) -> bool:
        return bool(self.find or self.prefix or self.suffix)


@dataclass
class RenameRow:
    source: Path
    new_name: str
    problem: str = ""
    companions: list[tuple[Path, Path]] = field(default_factory=list)

    @property
    def target(self) -> Path:
        return self.source.parent / self.new_name

    @property
    def changed(self) -> bool:
        return self.new_name != self.source.name


@dataclass
class RenamePlan:
    rows: list[RenameRow]
    error: str = ""  # e.g. a broken regular expression

    @property
    def changed_rows(self) -> list[RenameRow]:
        return [row for row in self.rows if row.changed]

    @property
    def problem_rows(self) -> list[RenameRow]:
        return [row for row in self.rows if row.problem]

    @property
    def ok(self) -> bool:
        return not self.error and not self.problem_rows and bool(self.changed_rows)

    @property
    def files(self) -> list[tuple[Path, Path]]:
        """(old, new) of the selected files that change name."""
        return [(row.source, row.target) for row in self.changed_rows]

    @property
    def moves(self) -> list[tuple[Path, Path]]:
        """Every rename to perform, paired annotation files included."""
        moves: list[tuple[Path, Path]] = []
        for row in self.changed_rows:
            moves.append((row.source, row.target))
            moves.extend(row.companions)
        return moves


def split_extension(name: str) -> tuple[str, str]:
    """("sample", ".fa.gz") for "sample.fa.gz"; no extension for ".hidden"."""
    if name.startswith(".") and name.count(".") == 1:
        return name, ""
    suffixes = PurePath(name).suffixes
    if not suffixes:
        return name, ""
    ext = suffixes[-1]
    if ext.lower() in _COMPRESSED and len(suffixes) >= 2 and len(suffixes[-2]) <= 7:
        ext = suffixes[-2] + ext
    if " " in ext or len(ext) > 16:
        return name, ""
    stem = name[: -len(ext)]
    if not stem:
        return name, ""
    return stem, ext


def _replace(text: str, options: RenameOptions) -> str:
    if not options.find:
        return text
    flags = 0 if options.match_case else re.IGNORECASE
    if options.regex:
        return re.sub(options.find, options.replace, text, flags=flags)
    pattern = re.compile(re.escape(options.find), flags)
    return pattern.sub(lambda _m: options.replace, text)


def new_name_for(name: str, options: RenameOptions) -> str:
    """The new file name; raises ``re.error`` for a broken regular expression."""
    if options.keep_extension:
        stem, ext = split_extension(name)
    else:
        stem, ext = name, ""
    stem = _replace(stem, options)
    return f"{options.prefix}{stem}{options.suffix}{ext}"


def name_problem(name: str, *, windows: bool | None = None) -> str:
    """Why ``name`` cannot be a file name ("" when it is fine)."""
    if windows is None:
        windows = platform.system() == "Windows"
    if not name.strip():
        return "the new name is empty"
    if name in {".", ".."}:
        return "not a valid file name"
    if "/" in name or "\0" in name:
        return "a file name cannot contain /"
    if windows:
        bad = sorted({c for c in name if c in _WINDOWS_BAD})
        if bad:
            return "a file name cannot contain " + " ".join(bad)
        if name.endswith((" ", ".")):
            return "a file name cannot end with a space or a dot"
    if len(name.encode("utf-8")) > _MAX_NAME_BYTES:
        return "the new name is too long"
    return ""


def _default_case_insensitive() -> bool:
    return platform.system() in {"Darwin", "Windows"}


def plan_batch_rename(
    paths: Iterable[Path],
    options: RenameOptions,
    *,
    companion_files: Callable[[Path], list[Path]] | None = None,
    companion_renames: Callable[[Path, Path], list[tuple[Path, Path]]] | None = None,
    case_insensitive: bool | None = None,
    listdir: Callable[[Path], list[str]] | None = None,
) -> RenamePlan:
    """New names for ``paths`` and whatever stops a file from being renamed.

    ``companion_files(path)`` lists files that always follow ``path`` (the
    annotation sidecar of a FASTA); they are left out of the list and renamed
    by ``companion_renames(old, new)`` together with their file.
    """
    ci = _default_case_insensitive() if case_insensitive is None else case_insensitive

    def key(name: str) -> str:
        return name.casefold() if ci else name

    def list_folder(folder: Path) -> list[str]:
        if listdir is not None:
            return list(listdir(folder))
        try:
            return os.listdir(folder)
        except OSError:
            return []

    unique: list[Path] = []
    seen: set[str] = set()
    for raw in paths:
        path = Path(raw)
        if str(path) in seen:
            continue
        seen.add(str(path))
        unique.append(path)

    followers: set[str] = set()
    if companion_files is not None:
        for path in unique:
            for companion in companion_files(path):
                if str(companion) != str(path):
                    followers.add(str(companion))
    sources = [path for path in unique if str(path) not in followers]

    rows: list[RenameRow] = []
    try:
        for path in sources:
            rows.append(RenameRow(path, new_name_for(path.name, options)))
    except re.error as exc:
        return RenamePlan(
            [RenameRow(path, path.name) for path in sources],
            error=f"Regular expression: {exc}",
        )

    for row in rows:
        if row.changed:
            row.problem = name_problem(row.new_name)
            if not row.problem and companion_renames is not None:
                row.companions = list(companion_renames(row.source, row.target))

    # Names that will exist afterwards, per folder: what is there now, minus
    # what moves away, plus every new name. Any clash is a problem. A row that
    # gets a problem keeps its old name, which can clash in turn, so repeat
    # until nothing changes.
    existing: dict[Path, set[str]] = {}

    def moves_of(row: RenameRow) -> bool:
        return row.changed and not row.problem

    for _ in range(len(rows) + 1):
        before = sum(1 for row in rows if row.problem)
        moving: dict[Path, set[str]] = {}
        for row in rows:
            if moves_of(row):
                for src, _dst in [(row.source, row.target), *row.companions]:
                    moving.setdefault(src.parent, set()).add(key(src.name))
        claimed: dict[tuple[Path, str], RenameRow] = {}
        clashes: list[tuple[RenameRow, RenameRow]] = []
        for row in rows:
            if moves_of(row):
                targets = [row.target, *[dst for _src, dst in row.companions]]
            else:
                targets = [row.source]
            for target in targets:
                slot = (target.parent, key(target.name))
                other = claimed.get(slot)
                if other is not None and other is not row:
                    clashes.append((row, other))
                else:
                    claimed[slot] = row
        for first, second in clashes:
            both_move = moves_of(first) and moves_of(second)
            for row in (first, second):
                if moves_of(row) and not row.problem:
                    row.problem = (
                        "same new name as another file"
                        if both_move
                        else "a file with this name already exists"
                    )
        for row in rows:
            if not moves_of(row):
                continue
            targets = [row.target, *[dst for _src, dst in row.companions]]
            for index, target in enumerate(targets):
                folder = target.parent
                if folder not in existing:
                    existing[folder] = {key(name) for name in list_folder(folder)}
                name_key = key(target.name)
                if name_key in existing[folder] and name_key not in moving.get(folder, set()):
                    row.problem = (
                        "a file with this name already exists"
                        if index == 0
                        else f"its annotation file {target.name} already exists"
                    )
                    break
        if sum(1 for row in rows if row.problem) == before:
            break
    return RenamePlan(rows)


def apply_moves(moves: list[tuple[Path, Path]]) -> str:
    """Rename every (old, new) pair, all or nothing. Returns an error or "".

    Step 1 moves each file to a temporary name in its folder, step 2 to its
    new name, so a -> b together with b -> a (or a case-only change on macOS)
    works. A new name that is still taken stops the rename and every file
    goes back to its old name. Nothing is ever overwritten.
    """
    moves = [(Path(a), Path(b)) for a, b in moves if Path(a) != Path(b)]
    if not moves:
        return ""
    token = uuid.uuid4().hex[:12]
    parked: list[tuple[Path, Path]] = []  # (original, temporary)
    placed: list[tuple[Path, Path]] = []  # (temporary, new)
    try:
        for index, (src, _dst) in enumerate(moves):
            temporary = src.parent / f".rv-rename-{token}-{index}"
            os.rename(src, temporary)
            parked.append((src, temporary))
        for (_src, temporary), (_old, dst) in zip(parked, moves):
            if os.path.lexists(dst):
                raise FileExistsError(f"{dst.name} already exists")
            os.rename(temporary, dst)
            placed.append((temporary, dst))
    except OSError as exc:
        for temporary, dst in reversed(placed):
            try:
                os.rename(dst, temporary)
            except OSError:
                pass
        for src, temporary in reversed(parked):
            try:
                os.rename(temporary, src)
            except OSError:
                pass
        return str(exc) or exc.__class__.__name__
    return ""


def check_moves(moves: list[tuple[Path, Path]], *, case_insensitive: bool | None = None) -> list[str]:
    """Problems that stop an Undo / Redo of a rename (empty when it can go ahead)."""
    ci = _default_case_insensitive() if case_insensitive is None else case_insensitive

    def norm(path: Path) -> str:
        text = str(path)
        return text.casefold() if ci else text

    leaving = {norm(Path(src)) for src, _dst in moves}
    problems: list[str] = []
    for src, dst in moves:
        src, dst = Path(src), Path(dst)
        if not src.exists():
            problems.append(f"{src.name}: no longer there")
        elif dst.exists() and norm(dst) not in leaving:
            problems.append(f"{dst.name}: this name is in use")
    return problems
