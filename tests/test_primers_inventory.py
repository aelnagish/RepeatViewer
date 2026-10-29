"""Primer Database as a lab inventory: shared folder (Sciebo-safe logs), sheet
import, pairs, duplicates, checking a new sequence against the stock, and the
Database window."""

from __future__ import annotations

import json
import time
import random
import shutil
from pathlib import Path

import pytest

pytest.importorskip("primer3")

from repeatviewer.backends.primers import inventory as inv  # noqa: E402
from repeatviewer.backends.primers import jobs as pj  # noqa: E402
from repeatviewer.backends.primers.oligos import Primer, PrimerDatabase  # noqa: E402
from repeatviewer.backends.primers.thermo import PrimerError  # noqa: E402

SHEET = Path(__file__).parent / "data" / "lab_primer_inventory.csv"


@pytest.fixture()
def shared(tmp_path):
    """A Sciebo-like folder: the same database opened by two people."""
    folder = inv.create_inventory(tmp_path / "Sciebo")
    return folder


def _imported(folder: Path, user: str = "Abdullah") -> inv.PrimerInventory:
    db = inv.PrimerInventory(folder, user=user)
    inv.apply_import(inv.plan_import(inv.read_table(SHEET), db), db)
    return db


# --------------------------------------------------------------------------
# Shared folder
# --------------------------------------------------------------------------


def test_create_and_basic_changes(shared):
    assert shared.name == inv.FOLDER_NAME and inv.is_inventory_folder(shared)
    assert json.loads((shared / inv.MARKER).read_text())["format"] == inv.FORMAT
    db = inv.PrimerInventory(shared, user="Abdullah")
    rid, = db.add_records([{"name": " CroSat2 _F ", "sequence": "gacttGCAAGTTCGAGAACTG",
                            "box": "A1", "status": "In stock", "dissolved": "TRUE"}])
    rec = db.records[rid]
    assert rec.name == "CroSat2_F"
    assert rec.get("status") == "in stock" and rec.get("dissolved") is True
    assert db.primers[0].tail == "gactt" and db.primers[0].kind == "forward"
    assert db.update(rid, {"box": "A2"}) and not db.update(rid, {"box": "A2"})
    assert [c.op for c in db.records[rid].history] == ["add", "update"]
    with pytest.raises(PrimerError, match="already called"):
        db.add_records([{"name": "x", "sequence": "ACGTACGTACGTACGTAA"}])
        x = db.record_by_name("x")
        db.update(x.rid, {"name": "CroSat2_F"})
    assert db.delete([rid]) == 1 and rid not in db.records
    # Only the writer's own log file exists.
    assert [p.name for p in db.log_files()] == [f"Abdullah@{inv.computer_name()}.jsonl"]


def test_two_people_never_lose_changes(shared):
    a = _imported(shared, "Abdullah")
    b = inv.PrimerInventory(shared, user="Anna")
    rid = a.record_by_name("pBV1_for").rid
    # Different fields of one primer: both kept.
    a.update(rid, {"box": "BI-9"})
    b.update(rid, {"tested": "worked"})
    for db in (a, b):
        db.reload()
        assert db.records[rid].get("box") == "BI-9"
        assert db.records[rid].get("tested") == "worked"
    # Same field: the later edit wins, on both computers.
    a.update(rid, {"status": "running low"})
    b.update(rid, {"status": "used up"})
    a.reload()
    b.reload()
    assert a.records[rid].get("status") == b.records[rid].get("status") == "used up"
    assert {c.user for c in a.records[rid].history} == {"Abdullah", "Anna"}
    assert set(a.writers()) == {"Abdullah", "Anna"}


def test_half_synced_line_and_conflicted_copy(shared):
    db = _imported(shared)
    n = len(db.records)
    log = db.log_files()[0]
    # A Sciebo "conflicted copy" of the same log: every change is read once.
    shutil.copy(log, log.with_name(log.stem + " (conflicted copy 2026-10-29).jsonl"))
    # A line still being synced (cut off).
    with open(log, "a", encoding="utf-8") as fh:
        fh.write('{"eid": "abc", "ts": 1')
    other = inv.PrimerInventory(shared, user="Nici")
    assert len(other.records) == n and other.skipped_lines == 1
    assert other.changed_on_disk() is False
    db.update(db.record_by_name("ITS4").rid, {"box": "D1-3"})
    assert other.changed_on_disk() and other.reload_if_changed()


def test_open_database_and_default_migration(tmp_path, monkeypatch):
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "cfg"))
    old = tmp_path / "cfg" / "RepeatViewer" / "primers" / "primer_database.fa"
    old.parent.mkdir(parents=True)
    old.write_text(">Old_F\nACGTGGTCAAGCTTGACC\n>Old_R\nGGTACCTTGACAAGCTTG\n")
    db = inv.open_database(None)
    assert isinstance(db, inv.PrimerInventory)
    assert db.path == inv.default_inventory_folder()
    assert sorted(db.names()) == ["Old_F", "Old_R"]
    assert not old.exists()  # moved in once
    assert sorted(inv.open_database(None).names()) == ["Old_F", "Old_R"]
    fasta = tmp_path / "x.fa"
    fasta.write_text(">A_F\nACGTGGTCAAGCTTGACC\n")
    assert type(inv.open_database(fasta)) is PrimerDatabase


def test_primer_database_api_on_inventory(shared):
    """Map, pickers and Design use the plain PrimerDatabase API."""
    db = inv.PrimerInventory(shared, user="Abdullah")
    added = db.add([Primer("Gene_F", "ACGTGGTCAAGCTTGACC"), Primer("Gene_R", "GGTACCTTGACAAGCTTG"),
                    Primer("Same", "ACGTGGTCAAGCTTGACC")])
    assert [p.name for p in added] == ["Gene_F", "Gene_R"]
    assert db.get("Gene_F").extra["id"] and db.find_sequence("ggtaccttgacaagcttg")
    assert db.remove(["Gene_R"]) == 1 and db.names() == ["Gene_F"]
    assert db.save() == shared


# --------------------------------------------------------------------------
# The lab sheet
# --------------------------------------------------------------------------


def test_import_lab_sheet_csv(shared):
    db = inv.PrimerInventory(shared, user="Abdullah")
    plan = inv.plan_import(inv.read_table(SHEET), db)
    assert set(plan.columns.values()) >= {"name", "sequence", "organism", "target",
                                          "tm_supplier", "dissolved", "ordered_by",
                                          "ordered_on", "box", "freezer", "room"}
    assert len(plan.new) == 28 and not plan.problems
    inv.apply_import(plan, db)
    r = db.record_by_name("SAIL_877_F01_F")     # "SAIL_877_F01 _F" in the sheet
    assert r is not None and r.get("box") == "E1-4" and r.get("organism") == "A.thaliana"
    p = db.record_by_name("pBV1_for")
    assert p.get("tm_supplier") == "53.7" and p.get("dissolved") is True
    assert p.get("room") == "Sterilbankraum" and p.get("freezer").startswith("-20")
    # Importing again adds nothing; a renamed copy is recognised by its sequence.
    rows = inv.read_table(SHEET)
    rows[3][0] = "BetaSat_F"
    again = inv.plan_import(rows, db)
    assert not again.new and again.renamed == [("BetaSat_F", "pBV1_for")]
    inv.apply_import(again, db)
    assert "BetaSat_F" in db.record_by_name("pBV1_for").get("aliases")


def test_import_xlsx(shared, tmp_path):
    openpyxl = pytest.importorskip("openpyxl")
    wb = openpyxl.Workbook()
    ws = wb.active
    rows = inv.read_table(SHEET)
    for r in rows:
        ws.append([True if v == "TRUE" else v for v in r])
    ws.cell(row=4, column=8).value = 46093          # a real date cell
    path = tmp_path / "Molbot Primer Inventory.xlsx"
    wb.save(path)
    db = inv.PrimerInventory(shared, user="Abdullah")
    plan = inv.plan_import(inv.read_table(path), db)
    assert len(plan.new) == 28
    first = plan.new[0]
    assert first["name"] == "pBV1_for" and first["dissolved"] is True
    assert first["ordered_on"] == "2026-03-12"


def test_import_problems_are_listed(shared, tmp_path):
    bad = tmp_path / "bad.tsv"
    bad.write_text("Name\tSequence\tBox\nA_F\tACGTACGTACGTACGTAC\tB1\nB_F\t\tB2\n"
                   "\tACGTTTTTGGGGCCCAAA\tB3\nC_F\tACGT ACGT XXXX\tB4\n"
                   "D_F\tACGTACGTACGTACGTAC\tB5\n")
    plan = inv.plan_import(inv.read_table(bad), inv.PrimerInventory(shared))
    assert [d["name"] for d in plan.new] == ["A_F"]
    assert len(plan.problems) == 4
    with pytest.raises(PrimerError, match="header row"):
        none = tmp_path / "none.csv"
        none.write_text("a,b\n1,2\n")
        inv.plan_import(inv.read_table(none), inv.PrimerInventory(shared))


def test_names_with_direction_inside():
    prs = [Primer(n, "ACGTGGTCAAGCTTGACC" + str(i) * 0 + "ACGT"[i % 4] * 3)
           for i, n in enumerate(["Pisum_PGene_f3_exon", "Pisum_PGene_r3_exon",
                                  "IGene_insert1f", "IGene_insert1r", "IGene_insert2r",
                                  "RGene_f1", "RGene_r1", "RGene_f2", "RGene_r2", "rbcLarev"])]
    pr = inv.partners(prs)
    assert pr["Pisum_PGene_f3_exon"] == ["Pisum_PGene_r3_exon"]
    assert pr["IGene_insert1f"] == ["IGene_insert1r"] and "IGene_insert2r" not in pr
    assert pr["RGene_f1"] == ["RGene_r1"] and pr["RGene_r2"] == ["RGene_f2"]
    assert "rbcLarev" not in pr
    assert inv.kind_of("Pisum_PGene_r3_exon") == "reverse"


def test_same_name_other_sequence_is_flagged(shared, tmp_path):
    sheet = tmp_path / "s.tsv"
    sheet.write_text("Primer name\tPrimer sequence\nITS4\tTCCTCCGCTTATTGATATGC\n"
                     "ITS4\tTCCTCCCGCTTATTGATATGC\n")
    plan = inv.plan_import(inv.read_table(sheet), inv.PrimerInventory(shared))
    assert [d["name"] for d in plan.new] == ["ITS4", "ITS4_2"]
    assert plan.notes and "Check which sequence is right" in plan.notes[0]
    assert "1 to check" in plan.summary() and not plan.problems


def test_pairs_and_duplicates(shared):
    db = _imported(shared)
    pairs = inv.partners(db.primers)
    assert pairs["pBV1_for"] == ["pBV1_rev"]
    assert sorted(pairs["fern_5S_rRNA_fwd"]) == ["fern_5S_rRNA_rev1", "fern_5S_rRNA_rev2"]
    assert pairs["Family04_FW"] == ["Family04_RV"]
    assert "ITS4" not in pairs  # no direction in the name
    # An explicit pair wins.
    db.update(db.record_by_name("ITS1-F").rid, {"pair": "ITS4"})
    assert inv.partners(db.primers)["ITS1-F"] == ["ITS4"]
    assert inv.find_duplicates(db.primers) == []
    db.add_records([{"name": "copy", "sequence": db.record_by_name("ITS4").sequence},
                    {"name": "rc", "sequence": inv.revcomp(db.record_by_name("pTS5_for").sequence)},
                    {"name": "longer", "sequence": "AAA" + db.record_by_name("pBV1_rev").sequence}])
    kinds = {(d.kind, frozenset((d.a, d.b))) for d in inv.find_duplicates(db.primers)}
    assert ("same sequence", frozenset(("ITS4", "copy"))) in kinds
    assert ("reverse complement", frozenset(("pTS5_for", "rc"))) in kinds
    assert ("contained", frozenset(("pBV1_rev", "longer"))) in kinds


def test_export_rows(shared):
    db = _imported(shared)
    rows = inv.export_rows(db)
    assert rows[0][:2] == ["Name", "Sequence"] and len(rows) == 29
    row = next(r for r in rows if r[0] == "pBV1_for")
    assert "BI-8" in row and "yes" in row


# --------------------------------------------------------------------------
# A new sequence against the stock
# --------------------------------------------------------------------------


def _new_plant(db: inv.PrimerInventory, path: Path) -> Path:
    rng = random.Random(4)

    def rnd(n):
        return "".join(rng.choice("ACGT") for _ in range(n))

    g = db.record_by_name
    f5 = g("fern_5S_rRNA_fwd").sequence
    f5_mm = f5[:5] + ("A" if f5[5] != "A" else "C") + f5[6:]     # 1 mismatch, far from 3'
    seq = (rnd(100) + g("pBV1_for").sequence + rnd(400) + inv.revcomp(g("pBV1_rev").sequence)
           + rnd(3000) + f5_mm + rnd(250) + inv.revcomp(g("fern_5S_rRNA_rev1").sequence)
           + rnd(3000) + g("ITS1-F").sequence + rnd(300))
    path.write_text(">new_plant_contig1\n" + seq + "\n")
    return path


def test_check_a_sequence_against_the_stock(shared, tmp_path):
    db = _imported(shared)
    fasta = _new_plant(db, tmp_path / "new_plant.fa")
    ms = pj.MapSettings(mode=pj.MAP_DATABASE, max_mismatches=2, three_prime_exact=3,
                        product_max=3000, search_probe=False)
    main, _files, payload = pj.run_map([fasta], ms, database_path=shared)
    st = payload["stock"]
    first, second = st["ready"][0], st["ready"][1]
    assert (first["forward"], first["reverse"], first["saved_pair"]) == (
        "pBV1_for", "pBV1_rev", True)
    assert first["mm_forward"] == first["mm_reverse"] == 0 and first["sizes"] == [436]
    assert first["forward_where"].startswith("box BI-8")
    assert (second["forward"], second["saved_pair"], second["mm_forward"]) == (
        "fern_5S_rRNA_fwd", True, 1)
    assert [x["name"] for x in st["single"]] == ["ITS1-F"]
    assert st["single"][0]["where"].startswith("box D1-2")
    html = main.read_text()
    assert "From your Primer Database" in html and "saved pair" in html
    assert "Single primers that bind" in html and "box BI-8" in html
    assert "—" not in html


def test_check_finds_nothing(shared, tmp_path):
    _imported(shared)
    fasta = tmp_path / "other.fa"
    rng = random.Random(9)
    fasta.write_text(">x\n" + "".join(rng.choice("ACGT") for _ in range(3000)) + "\n")
    main, _f, payload = pj.run_map([fasta], pj.MapSettings(mode=pj.MAP_DATABASE,
                                                           max_mismatches=1),
                                   database_path=shared)
    assert payload["stock"]["ready"] == [] and payload["stock"]["single"] == []
    assert "Design new primers" in main.read_text()


def test_designed_primers_go_in_as_designed(shared, tmp_path):
    from primer_testdata import make_data

    from repeatviewer.backends.primers import design as dz

    data = make_data(tmp_path / "d", genome_flank=2000)
    _m, _f, payload = pj.run_design([data.gene], dz.DesignSettings(), annotate=False)
    added = pj.save_pairs_to_database(payload, shared)
    assert len(added) == 10
    db = inv.PrimerInventory(shared)
    rec = db.record_by_name(added[0].name)
    assert rec.get("status") == "designed" and rec.get("target") == "gene1"
    assert inv.partners(db.primers)[added[0].name]


# --------------------------------------------------------------------------
# The Database window
# --------------------------------------------------------------------------


@pytest.fixture()
def settings(qapp, tmp_path):
    from PySide6.QtCore import QCoreApplication, QSettings

    org, app = QCoreApplication.organizationName(), QCoreApplication.applicationName()
    QCoreApplication.setOrganizationName("RVInventoryTests")
    QCoreApplication.setApplicationName(f"t{id(tmp_path)}")
    QSettings().clear()
    QSettings().setValue("primers/lab_user", "Abdullah")
    yield
    QSettings().clear()
    inv.set_lab_user("")
    QCoreApplication.setOrganizationName(org)
    QCoreApplication.setApplicationName(app)


def _row_of(d, name):
    for r in range(d._table.rowCount()):
        if d._table.item(r, 0).text() == name:
            return r
    raise AssertionError(name)


def test_window_import_edit_filter_history(qapp, settings, shared):
    from repeatviewer.views.primer_database_dialog import PrimerDatabaseDialog

    d = PrimerDatabaseDialog(shared)
    assert d._table.rowCount() == 0
    assert d.import_file(SHEET, confirm=False) == 28
    assert d._table.rowCount() == 28
    assert "User: <b>Abdullah</b>" in d._where.text()
    # Filters: organism and free text over every field.
    d._organism.setCurrentIndex(d._organism.findData("Nicotiana sp."))
    assert d._table.rowCount() == 4
    d._organism.setCurrentIndex(0)
    d._filter.setText("BI-8")
    assert d._table.rowCount() == 8
    d._filter.clear()
    # Pair column from the names.
    assert d._table.item(_row_of(d, "pBV1_for"), 1).text() == "pBV1_rev"
    # Edit one primer in the form.
    d._table.selectRow(_row_of(d, "pBV1_for"))
    qapp.processEvents()
    assert d._fields["box"].currentText() == "BI-8" and not d._save_btn.isEnabled()
    d._fields["tested"].setCurrentIndex(d._fields["tested"].findData("worked"))
    d._fields["anneal_temp"].setText("55")
    assert d._save_btn.isEnabled()
    d._save_form()
    rec = d.database().record_by_name("pBV1_for")
    assert rec.get("tested") == "worked" and rec.get("anneal_temp") == "55"
    assert rec.get("box") == "BI-8"           # untouched fields stay
    text = d._details.toPlainText()
    assert "History" in text and "Abdullah: Tested = worked" in text
    assert "supplier 53.7 °C" in text and "With pBV1_rev" in text
    assert "Hairpin" in text and "Self-dimer" in text and "3′ end" in text
    # Several at once: only the changed field is written to each.
    d._table.clearSelection()
    for name in ("ITS1-F", "ITS4"):
        d._table.selectRow(_row_of(d, name))
    from PySide6.QtWidgets import QAbstractItemView

    d._table.setSelectionMode(QAbstractItemView.SelectionMode.MultiSelection)
    d._table.clearSelection()
    d._table.selectRow(_row_of(d, "ITS1-F"))
    d._table.selectRow(_row_of(d, "ITS4"))
    qapp.processEvents()
    assert len(d._form_ids) == 2 and not d._fields["name"].isEnabled()
    d._fields["status"].setCurrentIndex(d._fields["status"].findData("running low"))
    d._save_form()
    db = d.database()
    assert db.record_by_name("ITS1-F").get("status") == "running low"
    assert db.record_by_name("ITS4").get("status") == "running low"
    assert db.record_by_name("ITS4").get("box") == "D1-2"
    d.close()


def test_window_sees_changes_from_others(qapp, settings, shared):
    from repeatviewer.views.primer_database_dialog import PrimerDatabaseDialog

    _imported(shared, "Anna")
    d = PrimerDatabaseDialog(shared)
    d.show()
    assert d._table.rowCount() == 28
    other = inv.PrimerInventory(shared, user="Nici")
    other.add_records([{"name": "NewSat_F", "sequence": "ACGTGGTCAAGCTTGACCAA", "box": "Z1"}])
    d._poll_disk()
    assert d._table.rowCount() == 29
    assert "2 people have made changes" in d._where.toolTip()   # Anna and Nici
    d.close()


def test_window_find_duplicates_and_shared_move(qapp, settings, shared, tmp_path, monkeypatch):
    from PySide6.QtWidgets import QFileDialog, QMessageBox

    from repeatviewer.views.primer_database_dialog import PrimerDatabaseDialog

    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "cfg"))
    d = PrimerDatabaseDialog(None)
    assert isinstance(d.database(), inv.PrimerInventory)
    assert "Database is on <b>this computer</b>" in d._where.text() and "On this computer only" in d._where.toolTip()
    d.add_primers([Primer("Loc_F", "ACGTGGTCAAGCTTGACCAA"), Primer("Loc_R", "TTGGTCAAGCTTGACCACGT")])
    moved = []
    d.databasePathChanged.connect(moved.append)
    monkeypatch.setattr(QFileDialog, "getExistingDirectory", lambda *a, **k: str(tmp_path / "Sciebo2"))
    (tmp_path / "Sciebo2").mkdir()
    monkeypatch.setattr(QMessageBox, "question", lambda *a, **k: QMessageBox.StandardButton.Yes)
    monkeypatch.setattr(QMessageBox, "information", lambda *a, **k: None)
    d._use_shared_folder()
    target = tmp_path / "Sciebo2" / inv.FOLDER_NAME
    assert moved == [str(target)] and sorted(d.database().names()) == ["Loc_F", "Loc_R"]
    assert "<b>a shared folder</b>" in d._where.text() and "Sciebo" in d._where.toolTip()
    # Loc_R is the reverse complement of Loc_F.
    monkeypatch.setattr(QMessageBox, "exec", lambda self: 0)
    d._find_duplicates()
    assert len(d._selected_ids()) == 2
    d.close()


# --------------------------------------------------------------------------
# People
# --------------------------------------------------------------------------


def test_people_from_the_sheet(shared):
    db = _imported(shared)
    names = [p.name for p in db.people()]
    assert names == ["Anna", "Loreen", "Nici", "Sara"]
    fern = db.record_by_name("fern_5S_rRNA_fwd")
    assert fern.get("ordered_by") == [inv.person_id("Loreen"), inv.person_id("Anna")]
    assert db.ordered_by_text(fern) == "Loreen, Anna"
    assert len(db.primers_of(inv.person_id("Anna"))) == 20   # Anna, Loreen & Anna, Sara & Anna
    assert inv.split_names("Anna & Nusrat") == ["Anna", "Nusrat"]
    assert inv.split_names("Anna und Daniel, Sara/Nici") == ["Anna", "Daniel", "Sara", "Nici"]


def test_people_add_edit_merge_remove(shared):
    db = _imported(shared)
    pid = db.add_person("Abdullah Yusuf", role="PhD student", email="a@rwth-aachen.de")
    assert db.persons[pid].initials == "AY"
    with pytest.raises(PrimerError, match="already"):
        db.add_person("abdullah  yusuf")
    db.update_person(pid, {"role": "Postdoc"})
    assert db.persons[pid].role == "Postdoc"
    rid = db.record_by_name("ITS4").rid
    db.update(rid, {"ordered_by": [inv.person_id("Anna"), pid]})
    assert db.ordered_by_text(db.records[rid]) == "Anna, Abdullah Yusuf"
    # Text still works (names are matched, new names become people).
    db.update(rid, {"ordered_by": "Anna & Nusrat"})
    assert db.ordered_by_text(db.records[rid]) == "Anna, Nusrat"
    assert db.person_by_name("Nusrat") is not None
    # Renaming a person shows the new name on every primer.
    db.update_person(inv.person_id("Nici"), {"name": "Nicola"})
    assert db.ordered_by_text(db.record_by_name("pBV1_for")) == "Nicola"
    # Merge: two entries for one person.
    dup = db.add_person("Anna K.")
    db.update(db.record_by_name("ITS1-F").rid, {"ordered_by": [dup]})
    assert db.merge_people(inv.person_id("Anna"), dup) == 1
    assert dup not in db.persons
    assert db.ordered_by_text(db.record_by_name("ITS1-F")) == "Anna"
    # Remove: with primers -> former member; without -> gone.
    assert db.remove_person(inv.person_id("Sara")) == "former"
    assert db.persons[inv.person_id("Sara")].active is False
    assert db.remove_person(pid) == "deleted" and pid not in db.persons
    # History speaks names, not ids.
    hist = [inv.describe_change(c, db) for c in db.records[rid].history]
    assert any("Ordered by = Anna, Nusrat" in h for h in hist)


def test_same_person_added_on_two_computers_is_one(shared):
    a = inv.PrimerInventory(shared, user="Abdullah")
    b = inv.PrimerInventory(shared, user="Anna")
    a.add_person("Daniel")
    b.add_person("Daniel")       # before Sciebo synced: b did not see a's entry
    a.reload()
    assert [p.name for p in a.people()] == ["Daniel"]


def test_old_text_is_migrated(shared):
    db = inv.PrimerInventory(shared, user="Abdullah")
    # An older database: ordered_by written as text.
    db._write([db._change("add", "r1", {"name": "X_F", "sequence": "ACGTACGTACGTACGTAA",
                                        "ordered_by": "Anna & Daniel"})])
    assert db.ordered_by_text(db.records["r1"]) == "Anna, Daniel"   # readable before
    assert db.migrate_people() == 1 and db.migrate_people() == 0
    assert db.records["r1"].get("ordered_by") == [inv.person_id("Anna"), inv.person_id("Daniel")]


def test_window_people(qapp, settings, shared, monkeypatch):
    from PySide6.QtWidgets import QInputDialog, QMessageBox

    from repeatviewer.views.primer_database_dialog import (
        PeopleDialog,
        PeoplePicker,
        PrimerDatabaseDialog,
    )

    _imported(shared)
    d = PrimerDatabaseDialog(shared)
    # Person filter.
    d._person.setCurrentIndex(d._person.findData(inv.person_id("Nici")))
    assert d._table.rowCount() == 8
    assert d._person.currentText() == "Nici (8)"
    d._person.setCurrentIndex(0)
    # The form shows names; the picker sets people.
    d._table.selectRow(_row_of(d, "ITS4"))
    qapp.processEvents()
    field = d._fields["ordered_by"]
    assert field._text.text() == "Anna"
    db = d.database()
    picker = PeoplePicker(db, field.value(), d)
    monkeypatch.setattr(QInputDialog, "getText", lambda *a, **k: ("Daniel", True))
    picker._add()
    assert picker.chosen() == [inv.person_id("Anna"), inv.person_id("Daniel")]
    field.set_value(picker.chosen())
    d._form_touched()
    assert d._save_btn.isEnabled()
    d._save_form()
    assert db.ordered_by_text(db.record_by_name("ITS4")) == "Anna, Daniel"
    from repeatviewer.views.primer_database_dialog import _COL

    assert d._table.item(_row_of(d, "ITS4"), _COL["ordered_by"]).text() == "Anna, Daniel"
    # Calculated Tm always with one decimal, next to the supplier's Tm.
    r = _row_of(d, "pBV1_for")
    assert d._table.item(r, _COL["tm"]).text().count(".") == 1
    assert d._table.item(r, _COL["tm_supplier"]).text() == "53.7"
    assert d.windowTitle() == "Primer Database: 28 primers"
    # People window: edit and mark as former.
    pd = PeopleDialog(db, d)
    assert pd._table.rowCount() == 5
    r = next(i for i in range(5) if pd._table.item(i, 0).text() == "Daniel")
    pd._table.selectRow(r)
    pd._role.setEditText("Master student")
    pd._save_person()
    assert db.persons[inv.person_id("Daniel")].role == "Master student"
    monkeypatch.setattr(QMessageBox, "question", lambda *a, **k: QMessageBox.StandardButton.Yes)
    pd._table.selectRow(r)
    pd._remove_selected()
    assert db.persons[inv.person_id("Daniel")].active is False   # has a primer: kept
    d.close()


def test_add_primers_window(qapp, settings, shared):
    from repeatviewer.backends.primers.thermo import ThermoSettings
    from repeatviewer.views.primer_database_dialog import AddPrimersDialog, PrimerDatabaseDialog

    db = _imported(shared)
    db.ensure_people(["Abdullah"])
    inv.set_lab_user("Abdullah")
    dlg = AddPrimersDialog(db, ThermoSettings(), None)
    assert dlg._table.rowCount() == 2 and not dlg._save.isEnabled()
    # The person adding is "Ordered by" by default; status Ordered; today.
    assert dlg._ordered_by.value() == [inv.person_id("Abdullah")]
    assert dlg._status.currentData() == "ordered"
    n = dlg.paste_text("CroSat2_F\tGACTTGCAAGTTCGAGAAC\nCroSat2_R\tGTTGCATCCACTTAGCTTAG\n"
                       f"Copy\t{db.record_by_name('ITS4').sequence}\nBad\tACGTXXZ\n")
    assert n == 4
    rows = dlg._table.rowCount()
    texts = {dlg._table.item(r, 0).text(): dlg._table.item(r, 5).text() for r in range(rows)}
    assert texts["CroSat2_F"] and "already" not in texts["CroSat2_F"]
    assert texts["Bad"] and texts["Bad"] != "OK"
    assert "already saved as ITS4" in texts["Copy"]
    assert dlg._table.item(0, 3).text().count(".") == 1
    dlg._organism.setEditText("Crocus sativus")
    dlg._target.setEditText("satDNA")
    dlg._box.setEditText("S1-1")
    assert dlg._save.text() == "Add 2"
    primers = dlg.primers()
    assert [p.name for p in primers] == ["CroSat2_F", "CroSat2_R"]
    d = PrimerDatabaseDialog(shared)
    assert d.add_primers(primers) == 2
    rec = d.database().record_by_name("CroSat2_F")
    assert rec.get("organism") == "Crocus sativus" and rec.get("box") == "S1-1"
    assert rec.get("status") == "ordered" and d.database().ordered_by_text(rec) == "Abdullah"
    assert rec.get("ordered_on") == time.strftime("%Y-%m-%d")
    d.close()


def test_window_remembers_its_layout(qapp, settings, shared):
    from repeatviewer.views.primer_database_dialog import _COL, PrimerDatabaseDialog

    _imported(shared)
    d = PrimerDatabaseDialog(shared)
    d.show()
    d.resize(1100, 650)
    d._table.setColumnWidth(_COL["name"], 222)
    d._split_bottom.setSizes([300, 500])
    qapp.processEvents()
    bottom = d._split_bottom.sizes()
    d.close()
    e = PrimerDatabaseDialog(shared)
    e.show()
    qapp.processEvents()
    assert e.height() == 650   # (the test screen is only 800 px wide)
    assert e._table.columnWidth(_COL["name"]) == 222
    e._fill()   # a refill keeps the user's column widths
    assert e._table.columnWidth(_COL["name"]) == 222
    got = e._split_bottom.sizes()
    assert abs(got[0] / sum(got) - bottom[0] / sum(bottom)) < 0.05
    e.close()


# --------------------------------------------------------------------------
# Lab lists and plain-language checks
# --------------------------------------------------------------------------


def test_lists_add_rename_merge_remove(shared):
    db = _imported(shared)
    assert ("Nicotiana sp.", 4) in db.values("organism")
    assert db.add_value("organism", "Crocus sativus") is True
    assert db.add_value("organism", "crocus  SATIVUS") is False     # already there
    assert ("Crocus sativus", 0) in db.values("organism")
    # Rename changes every primer; the list follows.
    assert db.rename_value("organism", "sect. Beta sp.", "Beta vulgaris") == 2
    assert db.record_by_name("pBV1_for").get("organism") == "Beta vulgaris"
    assert "sect. Beta sp." not in [v for v, _n in db.values("organism")]
    # Merge = rename into an existing value.
    db.rename_value("box", "AI-9", "BI-8")
    assert ("BI-8", 17) in db.values("box")
    # Remove: refused while used, unless cleared on the primers.
    assert db.remove_value("organism", "Crocus sativus") == 0
    assert db.remove_value("room", "Sterilbankraum") == -1
    assert db.remove_value("room", "Sterilbankraum", clear=True) == 28
    assert db.values("room") == []
    # Other computers see the list after syncing.
    other = inv.PrimerInventory(shared, user="Anna")
    assert ("Beta vulgaris", 2) in other.values("organism")


def test_plain_language_checks():
    from repeatviewer.backends.primers.thermo import ThermoSettings, analyze_oligo, heterodimer
    from repeatviewer.views.primer_database_dialog import BAD, OK, pair_verdicts, primer_verdicts

    t = ThermoSettings()
    good = analyze_oligo("AGCACCACCAGAGAGGAAGTACAGT", t)
    rows = {label: (status, text) for status, label, text in primer_verdicts(good)}
    assert rows["Hairpin"][0] == OK and rows["Self-dimer"][0] == OK and rows["3′ end"][0] == OK
    # A palindromic primer dimerizes strongly.
    pal = analyze_oligo("GCGCGAATTCGCGCGAATTCGCGC", t)
    rows = {label: status for status, label, _t in primer_verdicts(pal)}
    assert rows["Self-dimer"] == BAD
    a = analyze_oligo("AGCACCACCAGAGAGGAAGTACAGT", t)
    b = analyze_oligo("ATAATTATAAATATTAAATTTATAT", t)
    rows = {label: status for status, label, _t in pair_verdicts(
        a, b, heterodimer(a.sequence, b.sequence, t))}
    assert rows["Tm difference"] == BAD


def test_lab_window(qapp, settings, shared, monkeypatch):
    from PySide6.QtWidgets import QInputDialog

    from repeatviewer.views.primer_database_dialog import ListsDialog, PrimerDatabaseDialog

    _imported(shared)
    d = PrimerDatabaseDialog(shared)
    lab = ListsDialog(d.database(), d, start="places")
    assert [lab.tabs.tabText(i) for i in range(lab.tabs.count())] == [
        "People", "Places", "Organisms", "Targets"]
    assert lab.tabs.currentIndex() == 1
    boxes = lab.pages["box"]
    monkeypatch.setattr(QInputDialog, "getText", lambda *a, **k: ("Z-1", True))
    boxes._add()
    assert ("Z-1", 0) in d.database().values("box")
    lab.accept()
    d._refill_choices()
    assert d._fields["box"].findText("Z-1") >= 0     # offered in the form now
    d.close()


def test_check_only_the_selected_primers(shared, tmp_path, qapp, settings):
    db = _imported(shared)
    fasta = _new_plant(db, tmp_path / "new_plant.fa")
    ms = pj.MapSettings(mode=pj.MAP_DATABASE, max_mismatches=2, product_max=3000,
                        primer_names=["pBV1_for", "pBV1_rev"])
    _main, _f, payload = pj.run_map([fasta], ms, database_path=shared)
    st = payload["stock"]
    assert st["primers"] == 2
    assert [(r["forward"], r["reverse"]) for r in st["ready"]] == [("pBV1_for", "pBV1_rev")]
    # The Map window shows what it will test, and can go back to all primers.
    from repeatviewer.job_options import PrimerOptions
    from repeatviewer.views.primer_dialogs import MapPrimersDialog

    m = MapPrimersDialog(initial=PrimerOptions(map=ms.to_dict()), files=[])
    m.show()
    qapp.processEvents()
    assert m._form_main.form.isRowVisible(m._only_field)
    assert "2 selected" in m._only_label.text()
    assert m.map_settings().primer_names == ["pBV1_for", "pBV1_rev"]
    m._clear_only()
    assert m.map_settings().primer_names == []
    assert not m._form_main.form.isRowVisible(m._only_field)
    m.close()
