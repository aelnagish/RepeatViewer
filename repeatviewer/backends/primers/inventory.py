"""The lab's primer inventory: a Primer Database shared through a synced folder.

Made for Sciebo (or any folder that syncs between computers): nobody ever
writes a file somebody else writes too, so a sync can never produce a
conflicted copy that loses changes.

Folder layout::

    RepeatViewer Primer Database/
        database.json              marker (format, created)
        log/<user>@<computer>.jsonl  one append-only change log per person and computer

Every change (add, edit of some fields, delete) is one JSON line in the
writer's own log. Opening the database reads every log, orders the changes by
time and replays them, so the inventory is the same on every computer once
Sciebo has synced. Two people editing different fields of one primer both
keep their change; on the same field the later edit wins. Every change keeps
who and when, which is the history shown in the Database window.
"""

from __future__ import annotations

import getpass
import json
import os
import re
import socket
import time
import uuid
import zipfile
from dataclasses import dataclass, field
from pathlib import Path
from xml.etree import ElementTree as ET

from repeatviewer.backends.primers.oligos import (
    KIND_FORWARD,
    KIND_OLIGO,
    KIND_PROBE,
    KIND_REVERSE,
    PRIMER_KINDS,
    Primer,
    PrimerDatabase,
    kind_from_name,
    unique_name,
)
from repeatviewer.backends.primers.thermo import PrimerError, clean_oligo

FOLDER_NAME = "RepeatViewer Primer Database"
MARKER = "database.json"
FORMAT = "repeatviewer-primer-inventory-1"

# Record fields (key, label). ``name`` and ``sequence`` are required.
FIELDS: tuple[tuple[str, str], ...] = (
    ("name", "Name"),
    ("sequence", "Sequence"),
    ("kind", "Kind"),
    ("pair", "Pair"),
    ("aliases", "Other names"),
    ("organism", "Organism"),
    ("target", "Target"),
    ("tm_supplier", "Tm (supplier)"),
    ("dissolved", "Dissolved"),
    ("concentration", "Concentration"),
    ("box", "Box"),
    ("freezer", "Freezer"),
    ("room", "Room"),
    ("status", "Status"),
    ("tested", "Tested"),
    ("anneal_temp", "Annealing °C"),
    ("ordered_by", "Ordered by"),
    ("ordered_on", "Ordered on"),
    ("note", "Note"),
)
FIELD_KEYS = tuple(k for k, _ in FIELDS)
FIELD_LABELS = dict(FIELDS)

STATUS_VALUES = ("", "in stock", "running low", "used up", "ordered", "designed")
STATUS_LABELS = {"": "", "in stock": "In stock", "running low": "Running low",
                 "used up": "Used up", "ordered": "Ordered", "designed": "Designed (not ordered)"}
TESTED_VALUES = ("", "worked", "failed", "not tested")
TESTED_LABELS = {"": "", "worked": "Worked", "failed": "Failed", "not tested": "Not tested"}

_COMP = str.maketrans("ACGTRYKMSWBDHVNacgtrykmswbdhvn", "TGCAYRMKSWVHDBNtgcayrmkswvhdbn")


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


# --------------------------------------------------------------------------
# Who is writing
# --------------------------------------------------------------------------


def _safe(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", text).strip("_") or "user"


_LAB_USER = ""


def set_lab_user(name: str) -> None:
    """The name written into the change history (set from the app's settings)."""
    global _LAB_USER
    _LAB_USER = (name or "").strip()


def lab_user() -> str:
    return _LAB_USER or default_user()


def default_user() -> str:
    try:
        return getpass.getuser()
    except Exception:  # noqa: BLE001
        return "user"


_COMPUTER = ""


def computer_name() -> str:
    """A stable name for this computer (the host name can change with the network,
    e.g. on eduroam), kept in the RepeatViewer settings folder."""
    global _COMPUTER
    if _COMPUTER:
        return _COMPUTER
    try:
        from repeatviewer.app_preferences import config_dir

        store = config_dir() / "primers" / "computer_name.txt"
        if store.is_file():
            _COMPUTER = _safe(store.read_text(encoding="utf-8").strip())
            if _COMPUTER:
                return _COMPUTER
        try:
            host = socket.gethostname().split(".")[0]
        except OSError:
            host = ""
        if not host or re.search(r"\d+-\d+-", host):     # network names like 129-250-…
            host = "computer"
        _COMPUTER = _safe(f"{host}-{uuid.uuid4().hex[:4]}")
        store.parent.mkdir(parents=True, exist_ok=True)
        store.write_text(_COMPUTER + "\n", encoding="utf-8")
    except OSError:
        _COMPUTER = _COMPUTER or "computer"
    return _COMPUTER


# --------------------------------------------------------------------------
# Records and changes
# --------------------------------------------------------------------------


@dataclass
class Change:
    eid: str
    ts: float
    user: str
    op: str          # add | update | delete
    rid: str
    fields: dict = field(default_factory=dict)

    def to_json(self) -> str:
        return json.dumps({"eid": self.eid, "ts": self.ts, "user": self.user, "op": self.op,
                           "id": self.rid, "fields": self.fields}, ensure_ascii=False)


@dataclass
class Record:
    rid: str
    data: dict
    created_by: str = ""
    created_ts: float = 0.0
    changed_by: str = ""
    changed_ts: float = 0.0
    history: list[Change] = field(default_factory=list)

    def get(self, key: str) -> str:
        value = self.data.get(key, "")
        return "" if value is None else value

    @property
    def name(self) -> str:
        return str(self.data.get("name", ""))

    @property
    def sequence(self) -> str:
        return str(self.data.get("sequence", ""))

    def to_primer(self) -> Primer:
        extra = {k: self.data[k] for k in FIELD_KEYS
                 if k not in ("name", "sequence", "kind", "note") and self.data.get(k) not in
                 (None, "", False)}
        extra["id"] = self.rid
        return Primer(name=self.name, sequence=self.sequence,
                      kind=kind_of(self.name, str(self.data.get("kind", ""))),
                      note=str(self.data.get("note", "")), source="Primer Database", extra=extra)


def clean_record_fields(fields: dict) -> dict:
    """Normalise user / import input; raises PrimerError on a bad sequence."""
    out = {}
    for key, value in fields.items():
        if key not in FIELD_KEYS:
            continue
        if key == "dissolved":
            out[key] = _to_bool(value)
            continue
        if key == "ordered_by" and isinstance(value, (list, tuple)):
            out[key] = [str(v) for v in value if str(v).strip()]
            continue
        text = "" if value is None else str(value).strip()
        if key == "name":
            text = re.sub(r"\s*([_\-.])\s*", r"\1", text)
            text = re.sub(r"\s+", "_", text)
        elif key == "sequence":
            text = clean_oligo(text)
        elif key == "kind":
            text = text.lower()
            if text not in PRIMER_KINDS:
                text = ""
        elif key == "status":
            text = _pick(text, STATUS_VALUES, STATUS_LABELS)
        elif key == "tested":
            text = _pick(text, TESTED_VALUES, TESTED_LABELS)
        elif key == "tm_supplier":
            text = text.replace(",", ".")
        out[key] = text
    return out


def _pick(text: str, values: tuple, labels: dict) -> str:
    low = text.strip().lower()
    if low in values:
        return low
    for k, lab in labels.items():
        if lab.lower() == low:
            return k
    return low if low else ""


def _to_bool(value) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in ("1", "true", "yes", "y", "x", "✓", "ja", "wahr")


# --------------------------------------------------------------------------
# Lists (organisms, targets, places)
# --------------------------------------------------------------------------

LIST_PREFIX = "list:"
LIST_FIELDS: tuple[tuple[str, str, str], ...] = (
    # field, plural label, singular label
    ("organism", "Organisms", "organism"),
    ("target", "Targets", "target"),
    ("box", "Boxes", "box"),
    ("freezer", "Freezers", "freezer"),
    ("room", "Rooms", "room"),
)
LIST_LABELS = {k: (pl, sg) for k, pl, sg in LIST_FIELDS}


def list_entry_id(field_key: str, value: str) -> str:
    slug = re.sub(r"\s+", " ", (value or "").strip().lower())
    return f"{LIST_PREFIX}{field_key}:{slug}"


# --------------------------------------------------------------------------
# People
# --------------------------------------------------------------------------

PERSON_PREFIX = "person:"
PERSON_FIELDS: tuple[tuple[str, str], ...] = (
    ("name", "Name"), ("initials", "Initials"), ("role", "Role"), ("email", "Email"),
    ("active", "Active"),
)
PERSON_ROLES = ("", "Group leader", "Postdoc", "PhD student", "Master student",
                "Bachelor student", "Technician", "Guest")
_NAME_SPLIT = re.compile(r"\s*(?:&|,|;|/|\+|\band\b|\bund\b)\s*", re.I)


@dataclass
class Person:
    pid: str
    name: str
    initials: str = ""
    role: str = ""
    email: str = ""
    active: bool = True
    created_by: str = ""
    history: list = field(default_factory=list)

    @property
    def label(self) -> str:
        return self.name + ("" if self.active else " (former)")


def person_id(name: str) -> str:
    """Stable id from the name, so the same person added on two computers
    before Sciebo has synced is still one person."""
    slug = re.sub(r"[^a-z0-9]+", "-", (name or "").strip().lower()).strip("-")
    return PERSON_PREFIX + (slug or "unnamed")


def split_names(text: str) -> list[str]:
    """``"Anna & Nusrat"`` -> ``["Anna", "Nusrat"]``."""
    out = []
    for part in _NAME_SPLIT.split(str(text or "")):
        part = re.sub(r"\s+", " ", part).strip(" .")
        if part and part.lower() not in [x.lower() for x in out]:
            out.append(part)
    return out


def initials_of(name: str) -> str:
    words = [w for w in re.split(r"[\s\-]+", name) if w]
    if len(words) >= 2:
        return (words[0][0] + words[-1][0]).upper()
    return name[:2].upper()


# --------------------------------------------------------------------------
# The shared database
# --------------------------------------------------------------------------


def is_inventory_folder(path: Path | str | None) -> bool:
    if not path:
        return False
    p = Path(path).expanduser()
    return p.is_dir() and ((p / MARKER).is_file() or (p / "log").is_dir())


def find_inventory(folder: Path | str | None) -> Path | None:
    """The database in ``folder``: the folder itself, or one inside it (one level
    down), so choosing the lab's shared folder that holds it also works."""
    if not folder:
        return None
    p = Path(folder).expanduser()
    if is_inventory_folder(p):
        return p
    if not p.is_dir():
        return None
    named = p / FOLDER_NAME
    if is_inventory_folder(named):
        return named
    try:
        found = [c for c in sorted(p.iterdir()) if c.is_dir() and is_inventory_folder(c)]
    except OSError:
        return None
    return found[0] if len(found) == 1 else None


def create_inventory(parent: Path) -> Path:
    """Make (or reuse) ``<parent>/RepeatViewer Primer Database``."""
    parent = Path(parent).expanduser()
    folder = parent if parent.name == FOLDER_NAME or is_inventory_folder(parent) \
        else parent / FOLDER_NAME
    (folder / "log").mkdir(parents=True, exist_ok=True)
    marker = folder / MARKER
    if not marker.is_file():
        marker.write_text(json.dumps({"format": FORMAT,
                                      "created": time.strftime("%Y-%m-%d %H:%M"),
                                      "created_by": default_user()}, indent=1) + "\n",
                          encoding="utf-8")
    return folder


class PrimerInventory(PrimerDatabase):
    """A Primer Database folder shared by the lab (see the module docstring).

    It behaves like :class:`PrimerDatabase` (``primers``, ``get``, ``add``,
    ``remove``, ``save``) so Map, the pickers and Design use it unchanged;
    the Database window uses the record API (``records``, ``update`` …).
    """

    def __init__(self, folder: Path, *, user: str | None = None) -> None:
        self.path = Path(folder).expanduser()
        if not is_inventory_folder(self.path):
            raise PrimerError(f"Not a RepeatViewer primer database folder: {self.path}")
        self.set_user(user or "")
        self._signature: tuple = ()
        self.records: dict[str, Record] = {}
        self.persons: dict[str, Person] = {}
        self.lists: dict[str, dict[str, str]] = {}
        self.primers: list[Primer] = []
        self.skipped_lines = 0
        self.reload()

    def set_user(self, user: str) -> None:
        self.user = (user or "").strip() or lab_user()
        self._log = self.path / "log" / f"{_safe(self.user)}@{computer_name()}.jsonl"

    # ------------------------------------------------------------ reading

    def log_files(self) -> list[Path]:
        folder = self.path / "log"
        # Sciebo "conflicted copy" files are logs too: their lines are read and
        # de-duplicated, so nothing is lost even then.
        return sorted(p for p in folder.glob("*.jsonl*") if p.is_file())

    def signature(self) -> tuple:
        sig = []
        for p in self.log_files():
            try:
                st = p.stat()
            except OSError:
                continue
            sig.append((p.name, st.st_size, int(st.st_mtime_ns)))
        return tuple(sig)

    def changed_on_disk(self) -> bool:
        return self.signature() != self._signature

    def reload(self) -> None:
        changes: dict[str, Change] = {}
        skipped = 0
        for p in self.log_files():
            try:
                text = p.read_text(encoding="utf-8", errors="replace")
            except OSError:
                continue
            for line in text.splitlines():
                line = line.strip()
                if not line:
                    continue
                try:
                    d = json.loads(line)
                    ch = Change(str(d["eid"]), float(d["ts"]), str(d.get("user", "")),
                                str(d["op"]), str(d["id"]), dict(d.get("fields") or {}))
                except (ValueError, KeyError, TypeError):
                    skipped += 1   # half-synced last line; it completes on the next sync
                    continue
                changes[ch.eid] = ch
        records: dict[str, Record] = {}
        for ch in sorted(changes.values(), key=lambda c: (c.ts, c.eid)):
            rec = records.get(ch.rid)
            if ch.op == "add":
                if rec is None:
                    rec = Record(ch.rid, {}, created_by=ch.user, created_ts=ch.ts)
                    records[ch.rid] = rec
                rec.data.update(ch.fields)
            elif ch.op == "update" and rec is not None:
                rec.data.update(ch.fields)
            elif ch.op == "delete" and rec is not None:
                rec.data["_deleted"] = True
            else:
                continue
            rec.changed_by, rec.changed_ts = ch.user, ch.ts
            rec.history.append(ch)
        self.persons = {}
        for k, r in records.items():
            if not k.startswith(PERSON_PREFIX) or r.data.get("_deleted"):
                continue
            d = r.data
            self.persons[k] = Person(k, str(d.get("name") or k[len(PERSON_PREFIX):]),
                                     str(d.get("initials", "")), str(d.get("role", "")),
                                     str(d.get("email", "")),
                                     bool(d.get("active", True)), r.created_by,
                                     list(r.history))
        self.lists = {k: {} for k, _p, _s in LIST_FIELDS}
        for k, r in records.items():
            if k.startswith(LIST_PREFIX) and not r.data.get("_deleted"):
                fk, value = str(r.data.get("field", "")), str(r.data.get("value", "")).strip()
                if fk in self.lists and value:
                    self.lists[fk][k] = value
        self.records = {k: r for k, r in records.items() if not r.data.get("_deleted")
                        and not k.startswith((PERSON_PREFIX, LIST_PREFIX)) and r.sequence}
        self.skipped_lines = skipped
        self._signature = self.signature()
        self._rebuild_primers()

    def _rebuild_primers(self) -> None:
        out = []
        for rec in sorted(self.records.values(), key=lambda r: r.name.lower()):
            try:
                out.append(rec.to_primer())
            except PrimerError:
                continue
        self.primers = out

    # ------------------------------------------------------------ writing

    def _write(self, changes: list[Change]) -> None:
        if not changes:
            return
        self._log.parent.mkdir(parents=True, exist_ok=True)
        with open(self._log, "a", encoding="utf-8") as fh:
            for ch in changes:
                fh.write(ch.to_json() + "\n")
            fh.flush()
            try:
                os.fsync(fh.fileno())
            except OSError:
                pass
        self.reload()

    def _change(self, op: str, rid: str, fields: dict) -> Change:
        # Strictly increasing per writer, so two changes in one call keep their order.
        ts = time.time()
        last = getattr(self, "_last_ts", 0.0)
        if ts <= last:
            ts = last + 1e-6
        self._last_ts = ts
        return Change(uuid.uuid4().hex, ts, self.user, op, rid, fields)

    def add_records(self, rows: list[dict]) -> list[str]:
        """Add records (cleaned dicts with at least name and sequence). Returns ids."""
        taken = {r.name for r in self.records.values()}
        changes, ids = [], []
        rows = [clean_record_fields(row) for row in rows]
        self._resolve_people_in(rows)
        for data in rows:
            if not data.get("name") or not data.get("sequence"):
                raise PrimerError("Every primer needs a name and a sequence.")
            data["name"] = unique_name(data["name"], taken)
            taken.add(data["name"])
            rid = uuid.uuid4().hex[:12]
            changes.append(self._change("add", rid, data))
            ids.append(rid)
        self._write(changes)
        return ids

    def update(self, rid: str, fields: dict) -> bool:
        """Change some fields of one record; only real changes are written."""
        rec = self.records.get(rid)
        if rec is None:
            raise PrimerError("This primer was deleted meanwhile (by someone else?).")
        data = clean_record_fields(fields)
        self._resolve_people_in([data])
        if "name" in data:
            if not data["name"]:
                raise PrimerError("The name cannot be empty.")
            clash = next((r for r in self.records.values()
                          if r.name == data["name"] and r.rid != rid), None)
            if clash is not None:
                raise PrimerError(f"Another primer is already called {data['name']}.")
        if "sequence" in data and not data["sequence"]:
            raise PrimerError("The sequence cannot be empty.")
        diff = {k: v for k, v in data.items() if rec.data.get(k, "" if k != "dissolved"
                                                                 else False) != v}
        if not diff:
            return False
        self._write([self._change("update", rid, diff)])
        return True

    def delete(self, rids: list[str]) -> int:
        changes = [self._change("delete", rid, {}) for rid in rids if rid in self.records]
        self._write(changes)
        return len(changes)

    def record_by_name(self, name: str) -> Record | None:
        return next((r for r in self.records.values() if r.name == name), None)

    # ------------------------------------------------- PrimerDatabase API

    def reload_if_changed(self) -> bool:
        if self.changed_on_disk():
            self.reload()
            return True
        return False

    def add(self, primers: list[Primer], *, replace_same_name: bool = False) -> list[Primer]:
        rows, added = [], []
        seqs = {r.sequence.upper() for r in self.records.values()}
        for pr in primers:
            if pr.sequence.upper() in seqs:
                continue
            seqs.add(pr.sequence.upper())
            row = {"name": pr.name, "sequence": pr.sequence, "kind": pr.kind, "note": pr.note}
            for k in FIELD_KEYS:
                if k in pr.extra and k not in row:
                    row[k] = pr.extra[k]
            rows.append(row)
            added.append(pr)
        if rows:
            ids = self.add_records(rows)
            added = [self.records[i].to_primer() for i in ids if i in self.records]
        return added

    def remove(self, names: list[str]) -> int:
        rids = [r.rid for r in self.records.values() if r.name in set(names)]
        return self.delete(rids)

    def save(self) -> Path:
        # Every change is already written when it is made.
        return self.path

    # ------------------------------------------------------------ helpers

    def writers(self) -> list[str]:
        """Everyone who has written to this database (one log per person+computer)."""
        seen = []
        for p in self.log_files():
            who = p.name.split("@", 1)[0]
            if who not in seen:
                seen.append(who)
        return seen

    # ------------------------------------------------------------ lists

    def values(self, field_key: str) -> list[tuple[str, int]]:
        """Every value of a list field with how many primers use it: the saved
        list plus whatever primers already use. Sorted by name."""
        counts: dict[str, int] = {}
        canon: dict[str, str] = {}
        for value in self.lists.get(field_key, {}).values():
            canon.setdefault(value.lower(), value)
            counts.setdefault(value.lower(), 0)
        for r in self.records.values():
            v = str(r.data.get(field_key, "") or "").strip()
            if v:
                canon.setdefault(v.lower(), v)
                counts[v.lower()] = counts.get(v.lower(), 0) + 1
        return sorted(((canon[k], n) for k, n in counts.items()), key=lambda x: x[0].lower())

    def records_with(self, field_key: str, value: str) -> list[Record]:
        low = (value or "").strip().lower()
        return [r for r in self.records.values()
                if str(r.data.get(field_key, "") or "").strip().lower() == low]

    def add_value(self, field_key: str, value: str) -> bool:
        value = re.sub(r"\s+", " ", str(value or "")).strip()
        if field_key not in self.lists or not value:
            raise PrimerError("Give a name.")
        if any(v.lower() == value.lower() for v, _n in self.values(field_key)):
            return False
        self._write([self._change("add", list_entry_id(field_key, value),
                                  {"field": field_key, "value": value})])
        return True

    def rename_value(self, field_key: str, old: str, new: str) -> int:
        """Rename (or merge into an existing value) on every primer. Returns primers changed."""
        new = re.sub(r"\s+", " ", str(new or "")).strip()
        if not new:
            raise PrimerError("Give the new name.")
        changes = [self._change("update", r.rid, {field_key: new})
                   for r in self.records_with(field_key, old)]
        old_id, new_id = list_entry_id(field_key, old), list_entry_id(field_key, new)
        if old_id in self.lists.get(field_key, {}) and old_id != new_id:
            changes.append(self._change("delete", old_id, {}))
        if new_id not in self.lists.get(field_key, {}):
            changes.append(self._change("add", new_id, {"field": field_key, "value": new}))
        elif self.lists[field_key][new_id] != new:
            changes.append(self._change("update", new_id, {"value": new}))
        self._write(changes)
        return len([c for c in changes if c.rid in self.records])

    def remove_value(self, field_key: str, value: str, *, clear: bool = False) -> int:
        """Remove from the list. Used by primers: only with ``clear`` (the field is
        emptied on them). Returns primers changed, or -1 when refused."""
        used = self.records_with(field_key, value)
        if used and not clear:
            return -1
        changes = [self._change("update", r.rid, {field_key: ""}) for r in used]
        entry = list_entry_id(field_key, value)
        if entry in self.lists.get(field_key, {}):
            changes.append(self._change("delete", entry, {}))
        self._write(changes)
        return len(used)

    # ------------------------------------------------------------ people

    def people(self, *, include_former: bool = True) -> list[Person]:
        out = [p for p in self.persons.values() if include_former or p.active]
        return sorted(out, key=lambda p: (not p.active, p.name.lower()))

    def person_by_name(self, name: str) -> Person | None:
        return self.persons.get(person_id(name)) or next(
            (p for p in self.persons.values() if p.name.lower() == (name or "").strip().lower()),
            None)

    def ensure_people(self, names: list[str]) -> list[str]:
        """Person ids for these names; people not there yet are added."""
        ids, changes = [], []
        for name in names:
            name = re.sub(r"\s+", " ", str(name)).strip()
            if not name:
                continue
            p = self.person_by_name(name)
            pid = p.pid if p is not None else person_id(name)
            if p is None and pid not in [c.rid for c in changes]:
                changes.append(self._change("add", pid, {"name": name,
                                                         "initials": initials_of(name),
                                                         "active": True}))
            if pid not in ids:
                ids.append(pid)
        self._write(changes)
        return ids

    def add_person(self, name: str, **fields) -> str:
        name = re.sub(r"\s+", " ", str(name or "")).strip()
        if not name:
            raise PrimerError("A person needs a name.")
        if self.person_by_name(name) is not None:
            raise PrimerError(f"{name} is already in the list of people.")
        data = {"name": name, "initials": fields.get("initials") or initials_of(name),
                "role": fields.get("role", ""), "email": fields.get("email", ""),
                "active": bool(fields.get("active", True))}
        pid = person_id(name)
        self._write([self._change("add", pid, data)])
        return pid

    def update_person(self, pid: str, fields: dict) -> bool:
        p = self.persons.get(pid)
        if p is None:
            raise PrimerError("This person was removed meanwhile.")
        data = {k: v for k, v in fields.items() if k in dict(PERSON_FIELDS)}
        if "name" in data:
            data["name"] = re.sub(r"\s+", " ", str(data["name"])).strip()
            if not data["name"]:
                raise PrimerError("A person needs a name.")
            other = self.person_by_name(data["name"])
            if other is not None and other.pid != pid:
                raise PrimerError(f"{data['name']} is already in the list of people.")
        if "active" in data:
            data["active"] = bool(data["active"])
        current = {"name": p.name, "initials": p.initials, "role": p.role, "email": p.email,
                   "active": p.active}
        diff = {k: v for k, v in data.items() if current.get(k) != v}
        if not diff:
            return False
        self._write([self._change("update", pid, diff)])
        return True

    def primers_of(self, pid: str) -> list[Record]:
        return [r for r in self.records.values() if pid in self.ordered_by_ids(r)]

    def remove_person(self, pid: str) -> str:
        """Delete a person nobody's primers refer to; otherwise mark them former.
        Returns "deleted" or "former"."""
        if pid not in self.persons:
            return ""
        if self.primers_of(pid):
            self.update_person(pid, {"active": False})
            return "former"
        self._write([self._change("delete", pid, {})])
        return "deleted"

    def merge_people(self, keep: str, drop: str) -> int:
        """Two entries for one person: primers of ``drop`` go to ``keep``."""
        if keep == drop or keep not in self.persons or drop not in self.persons:
            return 0
        changes = []
        for rec in self.primers_of(drop):
            ids = [keep if x == drop else x for x in self.ordered_by_ids(rec)]
            ids = list(dict.fromkeys(ids))
            changes.append(self._change("update", rec.rid, {"ordered_by": ids}))
        changes.append(self._change("delete", drop, {}))
        self._write(changes)
        return len(changes) - 1

    def ordered_by_ids(self, rec: Record) -> list[str]:
        value = rec.data.get("ordered_by")
        if isinstance(value, list):
            return [str(v) for v in value]
        if value:   # older text value not turned into people yet
            return [p.pid if (p := self.person_by_name(n)) else person_id(n)
                    for n in split_names(str(value))]
        return []

    def names_of(self, ids) -> list[str]:
        out = []
        for pid in ids or []:
            p = self.persons.get(pid)
            out.append(p.name if p else pid[len(PERSON_PREFIX):] if
                       str(pid).startswith(PERSON_PREFIX) else str(pid))
        return out

    def ordered_by_text(self, rec: Record) -> str:
        value = rec.data.get("ordered_by")
        if isinstance(value, str):      # not turned into people yet
            return ", ".join(p.name if (p := self.person_by_name(n)) else n
                             for n in split_names(value))
        return ", ".join(self.names_of(self.ordered_by_ids(rec)))

    def _resolve_people_in(self, rows: list[dict]) -> None:
        """``ordered_by`` given as names (text) -> person ids, adding new people."""
        names: list[str] = []
        for row in rows:
            v = row.get("ordered_by")
            if isinstance(v, str):
                names += split_names(v)
        if names:
            self.ensure_people(names)
        for row in rows:
            v = row.get("ordered_by")
            if isinstance(v, str):
                row["ordered_by"] = [self.person_by_name(n).pid for n in split_names(v)
                                     if self.person_by_name(n) is not None]
            elif isinstance(v, list):
                row["ordered_by"] = [x for x in v if x in self.persons]

    def migrate_people(self) -> int:
        """Older databases: "Anna & Nusrat" text -> two people. Returns primers changed.

        Also: "Former member" used to be a role; it is now only the status, so
        such people get an empty role and are marked as former members."""
        for p in list(self.persons.values()):
            if p.role.strip().lower() == "former member":
                self.update_person(p.pid, {"role": "", "active": False})
        todo = [r for r in self.records.values()
                if isinstance(r.data.get("ordered_by"), str) and r.data.get("ordered_by")]
        if not todo:
            return 0
        names = [n for r in todo for n in split_names(r.data["ordered_by"])]
        self.ensure_people(names)
        changes = []
        for r in todo:
            ids = [self.person_by_name(n).pid for n in split_names(r.data["ordered_by"])
                   if self.person_by_name(n) is not None]
            changes.append(self._change("update", r.rid, {"ordered_by": ids}))
        self._write(changes)
        return len(changes)


def default_inventory_folder() -> Path:
    from repeatviewer.app_preferences import config_dir

    return config_dir() / "primers" / FOLDER_NAME


def ensure_default_inventory() -> Path:
    """The database on this computer; the old ``primer_database.fa`` is moved in once."""
    from repeatviewer.backends.primers.oligos import default_database_path, load_primers

    folder = default_inventory_folder()
    fresh = not is_inventory_folder(folder)
    create_inventory(folder)
    old = default_database_path()
    if fresh and old.is_file():
        try:
            db = PrimerInventory(folder)
            db.add(load_primers(old))
            old.rename(old.with_name(old.name + ".moved_to_database_folder"))
        except (OSError, PrimerError):
            pass
    return folder


def open_database(path: Path | str | None = None, *, user: str | None = None) -> PrimerDatabase:
    """The Primer Database at ``path``: a database folder (shared or on this
    computer), or an older single FASTA file. No path: the one on this computer."""
    from repeatviewer.backends.primers.oligos import default_database_path

    if path and Path(path).expanduser() == default_database_path():
        path = None   # the old default FASTA now lives in the database folder
    if not path:
        return PrimerInventory(ensure_default_inventory(), user=user)
    if is_inventory_folder(path):
        return PrimerInventory(Path(path), user=user)
    return PrimerDatabase(Path(path))


def copy_into(src: PrimerDatabase, dst: PrimerInventory) -> int:
    """Copy every primer (with all its details) into another database folder."""
    have = {r.sequence.upper() for r in dst.records.values()}
    rows = []
    if isinstance(src, PrimerInventory):
        for r in src.records.values():
            if r.sequence.upper() not in have:
                rows.append({k: v for k, v in r.data.items() if k in FIELD_KEYS})
    else:
        for p in src.primers:
            if p.sequence.upper() not in have:
                rows.append({"name": p.name, "sequence": p.sequence, "kind": p.kind,
                             "note": p.note})
    if rows:
        dst.add_records(rows)
    return len(rows)


# --------------------------------------------------------------------------
# Pairs, duplicates
# --------------------------------------------------------------------------

_PAIR_SUFFIX = re.compile(
    r"(?P<sep>[_\-. ]?)(?P<tag>fwd|for|forward|fw|left|rev|reverse|rv|right|f|r)(?P<num>\d*)$",
    re.I)


def pair_stem(name: str) -> tuple[str, str]:
    """``fern_5S_rRNA_rev2`` -> (``fern_5S_rRNA``, ``reverse``); no direction -> (name, "").

    The direction needs a separator before it (``_F``, ``-rev``, ``.fw``), or a
    digit before an upper-case F / R (``CroSat2F``), so names like ``Primer`` or
    ``ITS4`` are not taken for a direction.
    """
    text = (name or "").strip()
    m = _PAIR_SUFFIX.search(text)
    if not m or m.start() == 0:
        return text, ""
    tag = m["tag"]
    if not m["sep"]:
        before = text[m.start("tag") - 1]
        if not (before.isdigit() and tag in ("F", "R", "FW", "RV", "Fw", "Rv")):
            return text, ""
    low = tag.lower()
    kind = KIND_FORWARD if low in ("f", "fw", "fwd", "for", "forward", "left") else KIND_REVERSE
    return text[:m.start()].rstrip("_-. "), kind


def kind_of(name: str, kind: str = "") -> str:
    if kind in PRIMER_KINDS and kind != KIND_OLIGO:
        return kind
    stem, k = pair_stem(name)
    return k or direction_key(name)[2] or kind_from_name(name)


_DIR_TOKEN = re.compile(
    r"(?<=[_\-.\d])(?P<tag>fwd|for|forward|fw|left|rev|reverse|rv|right|f|r)(?=$|[_\-.\d(])",
    re.I)


def direction_key(name: str) -> tuple[str, str, str]:
    """(exact key, loose key, kind) from the last direction word in the name.

    ``Pisum_PGene_f3_exon`` and ``Pisum_PGene_r3_exon`` share the exact key
    ``Pisum_PGene_<>3_exon``; ``fern_5S_fwd`` and ``fern_5S_rev2`` only share
    the loose key (digits right after the direction word dropped).
    """
    text = (name or "").strip()
    last = None
    for m in _DIR_TOKEN.finditer(text):
        last = m
    if last is None:
        stem, kind = pair_stem(text)
        if not kind:
            return "", "", ""
        return stem.lower() + "|<>", stem.lower() + "|<>", kind
    low = last["tag"].lower()
    kind = KIND_FORWARD if low in ("f", "fw", "fwd", "for", "forward", "left") else KIND_REVERSE
    before, after = text[:last.start()], text[last.end():]
    exact = (before + "<>" + after).lower()
    loose = (before + "<>" + re.sub(r"^\d+", "", after)).lower()
    return exact, loose, kind


def partners(primers: list[Primer]) -> dict[str, list[str]]:
    """Name -> partner names: an explicit ``pair`` field, else from the names
    (``_for``/``_rev``, ``_f1``/``_r1``, ``insert1f``/``insert1r`` …)."""
    exact: dict[str, dict[str, list[str]]] = {}
    loose: dict[str, dict[str, list[str]]] = {}
    keys: dict[str, tuple[str, str, str]] = {}
    for p in primers:
        ek, lk, k = direction_key(p.name)
        if p.kind in (KIND_FORWARD, KIND_REVERSE) and k and p.kind != k:
            k = p.kind
        if not k:
            continue
        keys[p.name] = (ek, lk, k)
        exact.setdefault(ek, {}).setdefault(k, []).append(p.name)
        loose.setdefault(lk, {}).setdefault(k, []).append(p.name)
    out: dict[str, list[str]] = {}
    for p in primers:
        explicit = [x.strip() for x in str(p.extra.get("pair", "")).split(",") if x.strip()]
        if explicit:
            out[p.name] = explicit
            continue
        if p.name not in keys:
            continue
        ek, lk, k = keys[p.name]
        other = KIND_REVERSE if k == KIND_FORWARD else KIND_FORWARD
        names = exact.get(ek, {}).get(other) or loose.get(lk, {}).get(other) or []
        if names:
            out[p.name] = list(names)
    # A primer someone else names as its partner pairs back.
    for name, names in list(out.items()):
        for n in names:
            if n not in out:
                out[n] = [name]
    return out


def is_saved_pair(a: str, b: str, pairs: dict[str, list[str]]) -> bool:
    return b in pairs.get(a, []) or a in pairs.get(b, [])


@dataclass
class Duplicate:
    kind: str        # same | reverse complement | contained
    a: str
    b: str
    detail: str = ""


def find_duplicates(primers: list[Primer], *, min_contained: int = 18) -> list[Duplicate]:
    """Same sequence, reverse complement, or one binding part inside another."""
    out: list[Duplicate] = []
    seqs = [(p.name, p.binding.upper()) for p in primers if p.binding]
    by_seq: dict[str, list[str]] = {}
    for name, s in seqs:
        by_seq.setdefault(s, []).append(name)
    for s, names in by_seq.items():
        for other in names[1:]:
            out.append(Duplicate("same sequence", names[0], other))
        rc = revcomp(s)
        if rc != s and rc in by_seq and names[0] < by_seq[rc][0]:
            out.append(Duplicate("reverse complement", names[0], by_seq[rc][0]))
    uniq = sorted(by_seq.items(), key=lambda kv: len(kv[0]))
    for i, (s, names) in enumerate(uniq):
        if len(s) < min_contained:
            continue
        for t, tnames in uniq[i + 1:]:
            if len(t) > len(s) and s in t:
                out.append(Duplicate("contained", names[0], tnames[0],
                                     f"{names[0]} lies inside {tnames[0]}"))
    return out


# --------------------------------------------------------------------------
# Importing the lab's sheet (CSV, TSV, XLSX: e.g. Google Sheets downloads)
# --------------------------------------------------------------------------

_SYNONYMS: dict[str, tuple[str, ...]] = {
    "name": ("primer name", "oligo name", "name", "primer", "id"),
    "sequence": ("primer sequence", "oligo sequence", "sequence", "seq", "5'-3'", "5′-3′",
                 "sequence 5'-3'"),
    "organism": ("organism", "species", "plant", "taxon"),
    "target": ("amplified sequence", "target", "amplicon", "gene", "region", "locus"),
    "tm_supplier": ("tm", "tm [°c]", "tm (°c)", "melting temperature"),
    "dissolved": ("h2o added", "dissolved", "resuspended", "water added"),
    "concentration": ("concentration", "conc", "stock"),
    "box": ("located in box", "box", "position", "location"),
    "freezer": ("freezer", "fridge", "storage", "shelf"),
    "room": ("room", "lab"),
    "ordered_by": ("ordered by", "orderer", "who", "owner"),
    "ordered_on": ("ordered on", "order date", "date"),
    "status": ("status", "stock status"),
    "tested": ("tested", "validated", "works"),
    "anneal_temp": ("annealing", "annealing temperature", "ta"),
    "pair": ("pair", "partner"),
    "kind": ("kind", "direction", "orientation"),
    "note": ("note", "notes", "comment", "comments", "remarks"),
}


def _norm_header(text: str) -> str:
    return re.sub(r"\s+", " ", str(text or "").strip().lower())


def map_columns(header: list[str]) -> dict[int, str]:
    """Column index -> field key (best synonym match; each field once)."""
    out: dict[int, str] = {}
    used: set[str] = set()
    norm = [_norm_header(h) for h in header]
    # Exact synonym first, then "starts with" (e.g. "h2o added for 100 pmol/µl …").
    for exact in (True, False):
        for i, h in enumerate(norm):
            if i in out or not h:
                continue
            for key, syns in _SYNONYMS.items():
                if key in used:
                    continue
                if any(h == s for s in syns) if exact else any(
                        h.startswith(s) or s in h.split(" ") for s in syns if len(s) > 2):
                    out[i] = key
                    used.add(key)
                    break
    return out


def read_table(path: Path) -> list[list[str]]:
    """Rows of a CSV / TSV / XLSX file (first sheet)."""
    path = Path(path).expanduser()
    suffix = path.suffix.lower()
    if suffix in (".xlsx", ".xlsm"):
        return _read_xlsx(path)
    import csv

    text = path.read_text(encoding="utf-8-sig", errors="replace")
    sample = text[:5000]
    delim = "\t" if suffix == ".tsv" or sample.count("\t") > sample.count(",") else ","
    if delim == "," and sample.count(";") > sample.count(","):
        delim = ";"
    return [list(r) for r in csv.reader(text.splitlines(), delimiter=delim)]


_NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}


def _read_xlsx(path: Path) -> list[list[str]]:
    try:
        zf = zipfile.ZipFile(path)
    except zipfile.BadZipFile as exc:
        raise PrimerError(f"{path.name} is not a readable .xlsx file.") from exc
    with zf:
        shared: list[str] = []
        if "xl/sharedStrings.xml" in zf.namelist():
            root = ET.fromstring(zf.read("xl/sharedStrings.xml"))
            for si in root.findall("m:si", _NS):
                shared.append("".join(t.text or "" for t in si.iter(
                    "{%s}t" % _NS["m"])))
        sheets = sorted(n for n in zf.namelist()
                        if re.fullmatch(r"xl/worksheets/sheet\d+\.xml", n))
        if not sheets:
            raise PrimerError(f"{path.name} has no worksheet.")
        root = ET.fromstring(zf.read(sheets[0]))
    rows: list[list[str]] = []
    for row in root.iter("{%s}row" % _NS["m"]):
        cells: dict[int, str] = {}
        for c in row.findall("m:c", _NS):
            ref = c.get("r", "")
            col = _col_index(re.sub(r"\d", "", ref)) if ref else len(cells)
            t = c.get("t", "")
            v = c.find("m:v", _NS)
            if t == "s" and v is not None:
                val = shared[int(v.text)] if v.text and int(v.text) < len(shared) else ""
            elif t == "inlineStr":
                val = "".join(x.text or "" for x in c.iter("{%s}t" % _NS["m"]))
            elif t == "b" and v is not None:
                val = "TRUE" if v.text == "1" else "FALSE"
            else:
                val = v.text if v is not None and v.text is not None else ""
            cells[col] = val
        width = max(cells) + 1 if cells else 0
        rows.append([cells.get(i, "") for i in range(width)])
    return rows


def _col_index(letters: str) -> int:
    n = 0
    for ch in letters.upper():
        n = n * 26 + (ord(ch) - 64)
    return n - 1


@dataclass
class ImportPlan:
    columns: dict[int, str]
    header_row: int
    new: list[dict]
    updates: list[tuple[str, dict]]      # (record id, fields to fill)
    same: list[str]                      # already there, nothing new
    renamed: list[tuple[str, str]]       # (sheet name, name in the database)
    problems: list[str]
    notes: list[str] = field(default_factory=list)   # added, but worth a look

    def summary(self) -> str:
        parts = [f"{len(self.new)} new"]
        if self.updates:
            parts.append(f"{len(self.updates)} already in, with new details")
        if self.same:
            parts.append(f"{len(self.same)} already in")
        if self.problems:
            parts.append(f"{len(self.problems)} rows skipped")
        if self.notes:
            parts.append(f"{len(self.notes)} to check")
        return ", ".join(parts)


def plan_import(rows: list[list[str]], db: PrimerDatabase) -> ImportPlan:
    """Match sheet rows to the database: new primers, details to fill, problems."""
    header_row, columns = -1, {}
    for i, row in enumerate(rows[:30]):
        # A header can span two rows ("H2O added for 100 pmol/µl" above
        # "according to synthesis report"): join each cell with the one above.
        above = rows[i - 1] if i > 0 else []
        joined = [" ".join(x for x in ((above[j] if j < len(above) else ""), cell)
                           if str(x).strip()) for j, cell in enumerate(row)]
        cols = map_columns(row)
        for j, key in map_columns(joined).items():
            if j not in cols and key not in cols.values():
                cols[j] = key
        if "name" in cols.values() and "sequence" in cols.values():
            header_row, columns = i, cols
            break
    if header_row < 0:
        raise PrimerError("No header row with a name and a sequence column was found "
                          "(for example “Primer name” and “Primer sequence”).")
    existing: dict[str, Record] = {}
    if isinstance(db, PrimerInventory):
        for rec in db.records.values():
            existing[rec.sequence.upper()] = rec
    else:
        for p in db.primers:
            existing[p.sequence.upper()] = Record("", {"name": p.name, "sequence": p.sequence})
    plan = ImportPlan(columns, header_row, [], [], [], [], [])
    seen: dict[str, str] = {}
    names_used = {r.name for r in existing.values()}
    for n, row in enumerate(rows[header_row + 1:], start=header_row + 2):
        raw = {key: (row[i] if i < len(row) else "") for i, key in columns.items()}
        if not any(str(v).strip() for v in raw.values()):
            continue
        if not str(raw.get("sequence", "")).strip():
            if str(raw.get("name", "")).strip():
                plan.problems.append(f"Row {n}: {raw.get('name')} has no sequence.")
            continue
        if not str(raw.get("name", "")).strip():
            plan.problems.append(f"Row {n}: a sequence without a name.")
            continue
        try:
            data = clean_record_fields(raw)
        except PrimerError as exc:
            plan.problems.append(f"Row {n}: {raw.get('name')}: {exc}")
            continue
        if "ordered_on" in data:
            data["ordered_on"] = _excel_date(data["ordered_on"])
        data = {k: v for k, v in data.items() if v not in ("", None)}
        if not data.get("status") and data.get("box"):
            data["status"] = "in stock"      # it has a place in a box
        key = data["sequence"].upper()
        if key in seen:
            plan.problems.append(f"Row {n}: {data['name']} has the same sequence as "
                                 f"{seen[key]} (row above).")
            continue
        seen[key] = data["name"]
        rec = existing.get(key)
        if rec is None:
            if data["name"] in names_used:
                new_name = unique_name(data["name"], names_used)
                plan.notes.append(
                    f"Row {n}: {data['name']} is there already with another sequence; this "
                    f"one is added as {new_name}. Check which sequence is right.")
                data["name"] = new_name
            names_used.add(data["name"])
            plan.new.append(data)
            continue
        if rec.name != data["name"]:
            plan.renamed.append((data["name"], rec.name))
            aliases = [a.strip() for a in str(rec.get("aliases")).split(",") if a.strip()]
            if data["name"] not in aliases:
                aliases.append(data["name"])
            data["aliases"] = ", ".join(aliases)
        fill = {k: v for k, v in data.items() if k not in ("name", "sequence")
                and (k == "aliases" or rec.get(k) in ("", None, False, []))}
        if fill and rec.rid:
            plan.updates.append((rec.rid, fill))
        else:
            plan.same.append(rec.name)
    return plan


def apply_import(plan: ImportPlan, db: PrimerDatabase) -> int:
    if isinstance(db, PrimerInventory):
        db.add_records(plan.new)
        for rid, fill in plan.updates:
            try:
                db.update(rid, fill)
            except PrimerError:
                continue
        return len(plan.new)
    added = db.add([Primer(d["name"], d["sequence"], kind=d.get("kind", ""),
                           note=d.get("note", ""), source="import",
                           extra={k: v for k, v in d.items()
                                  if k not in ("name", "sequence", "kind", "note")})
                    for d in plan.new])
    db.save()
    return len(added)


def _excel_date(text: str) -> str:
    """Excel serial day numbers (e.g. 46093) -> 2026-03-12; anything else unchanged."""
    if re.fullmatch(r"\d{5}(\.0+)?", text or ""):
        import datetime as _dt

        day = _dt.date(1899, 12, 30) + _dt.timedelta(days=int(float(text)))
        return day.isoformat()
    return text


def export_rows(db: PrimerDatabase) -> list[list[str]]:
    """Header + one row per primer, every field (for TSV / sharing)."""
    rows = [[lab for _k, lab in FIELDS] + ["Added by", "Last changed"]]
    if isinstance(db, PrimerInventory):
        recs = sorted(db.records.values(), key=lambda r: r.name.lower())
        for r in recs:
            vals = []
            for k in FIELD_KEYS:
                v = r.get(k)
                if k == "dissolved":
                    v = "yes" if v else ""
                elif k == "status":
                    v = STATUS_LABELS.get(v, v)
                elif k == "tested":
                    v = TESTED_LABELS.get(v, v)
                elif k == "ordered_by":
                    v = db.ordered_by_text(r)
                vals.append(str(v))
            vals += [r.created_by, time.strftime("%Y-%m-%d %H:%M",
                                                 time.localtime(r.changed_ts)) + f" {r.changed_by}"]
            rows.append(vals)
    else:
        for p in db.primers:
            d = {"name": p.name, "sequence": p.sequence, "kind": p.kind, "note": p.note, **p.extra}
            rows.append([str(d.get(k, "")) for k in FIELD_KEYS] + ["", ""])
    return rows


def describe_change(ch: Change, db: "PrimerInventory | None" = None) -> str:
    when = time.strftime("%Y-%m-%d %H:%M", time.localtime(ch.ts))
    if ch.op == "add":
        return f"{when}  {ch.user}: added"
    if ch.op == "delete":
        return f"{when}  {ch.user}: deleted"
    parts = []
    for k, v in ch.fields.items():
        label = FIELD_LABELS.get(k, k)
        if k == "dissolved":
            v = "yes" if v else "no"
        elif k == "ordered_by" and isinstance(v, list):
            v = ", ".join(db.names_of(v)) if db is not None else ", ".join(v)
        elif k == "active":
            v = "yes" if v else "no (former member)"
        parts.append(f"{label} = {v}" if v != "" else f"{label} cleared")
    return f"{when}  {ch.user}: " + "; ".join(parts)


__all__ = [
    "FIELDS", "FOLDER_NAME", "KIND_PROBE", "PrimerInventory", "Record", "apply_import",
    "create_inventory", "describe_change", "export_rows", "find_duplicates",
    "is_inventory_folder", "is_saved_pair", "open_database", "pair_stem", "partners",
    "plan_import", "read_table",
]
