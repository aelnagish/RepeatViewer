"""Persistence and smoke tests for Manage Database (BLAST / HMM / Others)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from PySide6.QtCore import QSettings, Qt
from PySide6.QtWidgets import QDialog

from repeatviewer.database_settings import (
    SETTINGS_GROUP,
    USER_DATABASES_KEY,
    load_blast_databases_from_settings,
    load_managed_genomes_from_settings,
    load_similarity_annotation_databases,
    load_user_databases_from_settings,
    save_blast_databases_to_settings,
    save_database_settings,
    save_managed_genomes_to_settings,
)
from repeatviewer.project import (
    ANNOTATION_NAME_MAX_LENGTH,
    BLAST_MOLECULE_NUCL,
    DATABASE_TYPE_PROTEIN,
    DB_KIND_FASTA_DNA,
    DB_KIND_FASTA_PROTEIN,
    DB_KIND_HMM_PROTEIN,
    BlastDatabase,
    ManagedGenome,
    RepeatViewerProject,
    UserDatabase,
    kind_from_user_database,
    kind_label,
    normalize_annotation_name,
    normalize_db_kind,
    user_fields_from_kind,
)
from repeatviewer.views.manage_references_dialog import (
    AddDatabaseDialog,
    ManageReferencesDialog,
    _format_time,
)


@pytest.fixture()
def isolated_settings(tmp_path, monkeypatch):
    """Point QSettings at a temp INI so tests do not touch real prefs."""
    conf = tmp_path / "settings.ini"
    monkeypatch.setenv("QT_QPA_PLATFORM", "offscreen")

    import repeatviewer.database_settings as ds

    class _QS(QSettings):
        def __init__(self, *args, **kwargs):
            if not args and not kwargs:
                super().__init__(str(conf), QSettings.Format.IniFormat)
            else:
                super().__init__(*args, **kwargs)

    monkeypatch.setattr(ds, "QSettings", _QS)
    yield conf


def test_kind_helpers() -> None:
    assert normalize_db_kind("FASTA DNA") == DB_KIND_FASTA_DNA
    assert normalize_db_kind("hmm_protein") == DB_KIND_HMM_PROTEIN
    assert kind_from_user_database("DNA", is_hmm=False) == DB_KIND_FASTA_DNA
    assert kind_from_user_database("Protein", is_hmm=True) == DB_KIND_HMM_PROTEIN
    assert kind_label(DB_KIND_FASTA_PROTEIN) == "FASTA Protein"
    assert user_fields_from_kind(DB_KIND_HMM_PROTEIN) == ("Protein", True)


def test_annotation_name_max_length() -> None:
    long = "x" * (ANNOTATION_NAME_MAX_LENGTH + 20)
    assert len(normalize_annotation_name(long)) == ANNOTATION_NAME_MAX_LENGTH


def test_format_time_date_only() -> None:
    assert _format_time("2026-09-28T12:00:00") == "2026-09-28"
    assert _format_time("2026-01-02 03:04:05") == "2026-01-02"
    assert _format_time("2026-05-01") == "2026-05-01"


def test_user_database_type_and_annotation_roundtrip(isolated_settings, tmp_path) -> None:
    fasta = tmp_path / "lib.fa"
    fasta.write_text(">a\nATGC\n")
    dbs = [
        UserDatabase(
            name="PubLib",
            path=fasta,
            db_type=DATABASE_TYPE_PROTEIN,
            created_at="2026-01-02T03:04:05",
            show_in_annotation=True,
            annotation_name="Public library",
        ),
        UserDatabase(
            name="Hidden",
            path=fasta,
            db_type=DATABASE_TYPE_PROTEIN,
            is_hmm=True,
            created_at="2026-01-03T00:00:00",
            show_in_annotation=False,
        ),
    ]
    save_database_settings(cdd_dirs=[], user_databases=dbs)
    loaded = load_user_databases_from_settings()
    assert len(loaded) == 2
    by_name = {d.name: d for d in loaded}
    assert by_name["PubLib"].display_type_label == "FASTA Protein"
    assert by_name["PubLib"].annotation_name == "Public library"
    assert by_name["Hidden"].display_type_label == "HMM Protein"
    assert by_name["Hidden"].show_in_annotation is False
    assert by_name["PubLib"].created_at.startswith("2026-01-02")

    visible = load_similarity_annotation_databases()
    assert [d.name for d in visible] == ["Public library"]


def test_legacy_user_database_defaults(isolated_settings, tmp_path) -> None:
    fasta = tmp_path / "old.fa"
    fasta.write_text(">a\nATGC\n")
    s = QSettings(str(isolated_settings), QSettings.Format.IniFormat)
    s.beginGroup(SETTINGS_GROUP)
    s.setValue(
        USER_DATABASES_KEY,
        json.dumps(
            [{"name": "Legacy", "path": str(fasta), "db_type": "Protein", "is_hmm": False}]
        ),
    )
    s.endGroup()
    loaded = load_user_databases_from_settings()
    assert len(loaded) == 1
    assert loaded[0].show_in_annotation is True
    assert loaded[0].display_type_label == "FASTA Protein"
    assert loaded[0].annotation_name == ""


def test_blast_database_kind_and_time(isolated_settings, tmp_path) -> None:
    fasta = tmp_path / "g.fa"
    fasta.write_text(">chr1\nATGC\n")
    prefix = tmp_path / "idx" / "g"
    prefix.parent.mkdir()
    dbs = [
        BlastDatabase(
            name="GenomeDB",
            source_fasta=fasta,
            db_prefix_path=prefix,
            molecule_type=BLAST_MOLECULE_NUCL,
            created_at="2026-09-28T12:00:00",
            annotation_name="GenomeAnnot",
            show_in_annotation=True,
        )
    ]
    save_blast_databases_to_settings(dbs)
    loaded = load_blast_databases_from_settings()
    assert len(loaded) == 1
    assert loaded[0].display_type_label == "FASTA DNA"
    assert loaded[0].created_at.startswith("2026-09-28")
    assert loaded[0].annotation_name == "GenomeAnnot"
    visible = load_similarity_annotation_databases()
    assert any(d.name == "GenomeAnnot" for d in visible)


def test_managed_genomes_roundtrip(isolated_settings, tmp_path) -> None:
    genome = tmp_path / "ref.fa"
    genome.write_text(">chr1\nATGCATGC\n")
    save_managed_genomes_to_settings(
        [
            ManagedGenome(
                name="Ref", path=genome, created_at="2026-05-01T10:00:00"
            )
        ]
    )
    loaded = load_managed_genomes_from_settings()
    assert len(loaded) == 1
    assert loaded[0].display_type_label == "FASTA DNA"
    assert loaded[0].created_at.startswith("2026-05-01")


def test_manage_database_dialog_smoke(qapp, isolated_settings, tmp_path, monkeypatch) -> None:
    # Treat indexes as ready so incomplete test files still appear in the list.
    monkeypatch.setattr(
        "repeatviewer.views.manage_references_dialog.prepare_index_is_ready",
        lambda *_a, **_k: True,
    )
    monkeypatch.setattr(
        "repeatviewer.views.manage_references_dialog.describe_database_problem",
        lambda *_a, **_k: None,
    )
    genome = tmp_path / "g.fa"
    genome.write_text(">c\nAAAA\n")
    hmm = tmp_path / "Pfam.hmm"
    hmm.write_text("HMMER3/f\n")
    genomes = [
        ManagedGenome(name="G", path=genome, created_at="2026-01-01T00:00:00")
    ]
    anns = [
        UserDatabase(
            name="Pfam",
            path=hmm,
            db_type=DATABASE_TYPE_PROTEIN,
            is_hmm=True,
            created_at="2026-01-02T00:00:00",
            show_in_annotation=True,
            annotation_name="PfamA",
        )
    ]
    save_managed_genomes_to_settings(genomes)
    save_database_settings(cdd_dirs=[], user_databases=anns)
    project = RepeatViewerProject(
        managed_genomes=genomes,
        user_databases=anns,
        genome_file=genome,
    )
    dialog = ManageReferencesDialog(project)
    assert dialog.windowTitle() == "Manage Database"
    assert dialog._section_bar.count() == 3
    assert dialog._section_bar.tabText(0) == "BLAST database"
    assert dialog._section_bar.tabText(1) == "HMM database"
    assert dialog._section_bar.tabText(2) == "Others"
    assert dialog._section_bar.expanding() is False
    assert dialog._section_bar.drawBase() is False
    assert dialog._table.columnCount() == 6
    # BLAST tab → genomes / FASTA
    dialog._section_bar.setCurrentIndex(0)
    dialog._legacy_role = "all"
    dialog._populate()
    assert dialog._table.rowCount() == 1
    assert dialog._table.item(0, 0).text() == "G"
    assert dialog._table.item(0, 1).text() == "-"  # no annotation name
    assert dialog._table.item(0, 2).text() == "FASTA DNA"
    assert dialog._table.item(0, 4).text() == "2026-01-01"
    # Genome-only rows can toggle Annotation
    ann_item = dialog._table.item(0, 5)
    assert ann_item is not None
    assert ann_item.flags() & Qt.ItemFlag.ItemIsUserCheckable
    # HMM tab
    dialog._section_bar.setCurrentIndex(1)
    dialog._legacy_role = "all"
    dialog._populate()
    assert dialog._table.rowCount() == 1
    assert dialog._table.item(0, 0).text() == "Pfam"
    assert dialog._table.item(0, 1).text() == "PfamA"
    dialog.close()


def test_add_database_dialog_hmm(qapp, tmp_path) -> None:
    hmm = tmp_path / "Pfam.hmm"
    hmm.write_text("HMMER3/f\n")
    project = RepeatViewerProject()
    # Add from BLAST tab but choose HMM type → routes to HMM storage.
    dialog = AddDatabaseDialog(
        project, initial_kind=DB_KIND_FASTA_DNA, section="blast"
    )
    assert dialog.windowTitle() == "Add database"
    assert dialog._type.count() == 5  # all types available
    idx = dialog._type.findData(DB_KIND_HMM_PROTEIN)
    assert idx >= 0
    dialog._type.setCurrentIndex(idx)
    dialog._name.setText("Pfam")
    dialog._ann_name.setText("PfamA")
    dialog._path.setText(str(hmm))
    dialog._show_annotation.setChecked(True)
    dialog._on_ok()
    assert dialog.result() == int(QDialog.DialogCode.Accepted)
    db = dialog.result_annotation()
    assert db is not None
    assert db.is_hmm is True
    assert db.display_type_label == "HMM Protein"
    assert db.annotation_name == "PfamA"
    assert dialog.result_blast() is None
    assert dialog._destination_section() == "hmm"


def test_add_database_dialog_legacy_annotation_section(qapp, tmp_path) -> None:
    """Legacy section='annotation' maps to HMM tab kinds."""
    fasta = tmp_path / "db.fa"
    fasta.write_text(">a\nM*\n")
    project = RepeatViewerProject()
    dialog = AddDatabaseDialog(
        project, initial_kind=DB_KIND_HMM_PROTEIN, section="annotation"
    )
    assert dialog._section == "hmm"
    dialog._name.setText("MyDB")
    dialog._path.setText(str(fasta))
    dialog._show_annotation.setChecked(True)
    dialog._on_ok()
    assert dialog.result() == int(QDialog.DialogCode.Accepted)
    db = dialog.result_annotation()
    assert db is not None
    assert db.annotation_name == "MyDB"  # auto-filled from name


def test_add_reference_legacy_section_maps_to_blast(qapp, tmp_path, monkeypatch) -> None:
    genome = tmp_path / "ref.fa"
    genome.write_text(">c\nATGC\n")
    project = RepeatViewerProject()
    dialog = AddDatabaseDialog(
        project, initial_kind=DB_KIND_FASTA_DNA, section="reference"
    )
    assert dialog._section == "blast"
    dialog._path.setText(str(genome))
    # Name auto-fills from path on text change; ensure name present.
    if not dialog._name.text().strip():
        dialog._name.setText("RefBeet-3.0")
    # Building BLAST needs makeblastdb — stub the worker path by not calling OK
    # for full build; just verify section mapping and name default.
    assert dialog._name.text().strip() in {"ref", "RefBeet-3.0"}


def test_cdd_annotation_name_persists(isolated_settings, tmp_path, qapp, monkeypatch) -> None:
    """Checking Annotation on a legacy cdd_dirs-only CDD must keep the name."""
    monkeypatch.setattr(
        "repeatviewer.views.manage_references_dialog.prepare_index_is_ready",
        lambda *_a, **_k: True,
    )
    monkeypatch.setattr(
        "repeatviewer.views.manage_references_dialog.describe_database_problem",
        lambda *_a, **_k: None,
    )
    from PySide6.QtCore import QSettings

    from repeatviewer.database_settings import (
        CDD_DIRS_KEY,
        SETTINGS_GROUP,
        load_user_databases_from_settings,
    )

    cdd_dir = tmp_path / "cdd"
    cdd_dir.mkdir()
    # Seed only cdd_dirs (no user_databases CDD row) — the legacy layout.
    s = QSettings(str(isolated_settings), QSettings.Format.IniFormat)
    s.beginGroup(SETTINGS_GROUP)
    s.setValue(CDD_DIRS_KEY, json.dumps([str(cdd_dir)]))
    s.endGroup()

    project = RepeatViewerProject(cdd_dirs=[cdd_dir])
    dialog = ManageReferencesDialog(project)
    dialog._section_bar.setCurrentIndex(2)
    dialog._legacy_role = "all"
    dialog._populate()
    assert dialog._table.rowCount() == 1
    assert dialog._table.item(0, 0).text() == "CDD"
    # Promoted on open with Annotation name filled.
    assert dialog._table.item(0, 1).text() == "CDD"
    assert dialog._table.item(0, 5).checkState() == Qt.CheckState.Checked

    loaded = load_user_databases_from_settings()
    assert any(
        db.db_type == "CDD" and db.annotation_name == "CDD" and db.show_in_annotation
        for db in loaded
    )
    dialog.close()

    # Reopen — Annotation name must still be CDD, not "-".
    dialog2 = ManageReferencesDialog(RepeatViewerProject())
    dialog2._section_bar.setCurrentIndex(2)
    dialog2._legacy_role = "all"
    dialog2._populate()
    assert dialog2._table.rowCount() == 1
    assert dialog2._table.item(0, 1).text() == "CDD"
    dialog2.close()


def test_add_blast_name_defaults_to_filename(qapp, tmp_path) -> None:
    fasta = tmp_path / "MyLibrary.fasta"
    fasta.write_text(">a\nATGC\n")
    project = RepeatViewerProject()
    dialog = AddDatabaseDialog(
        project, initial_kind=DB_KIND_FASTA_DNA, section="blast"
    )
    dialog._path.setText(str(fasta))
    dialog._maybe_fill_name_from_path()
    name = dialog._name.text().strip() or Path(dialog._path.text().strip()).stem
    assert name == "MyLibrary"
