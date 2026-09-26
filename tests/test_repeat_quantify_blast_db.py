"""Repeat Quantify reuses saved BLAST databases instead of indexing again."""

import os
import time
from pathlib import Path

import pytest

from repeatviewer.backends import database_prepare as dp
from repeatviewer.backends import repeat_quantify as rq
from repeatviewer.project import BlastDatabase


def _fasta(path: Path) -> Path:
    path.write_text(">s1\nACGTACGTACGT\n")
    return path


def _index(prefix: Path, suffix: str = ".nin", *, age: float = 0.0) -> Path:
    """Fake makeblastdb output; ``age`` > 0 makes it older than now."""
    prefix.parent.mkdir(parents=True, exist_ok=True)
    marker = Path(f"{prefix}{suffix}")
    marker.write_text("")
    if age:
        old = time.time() - age
        os.utime(marker, (old, old))
    return marker


def _db(name, fasta, prefix, molecule_type="nucl"):
    return BlastDatabase(
        name=name, source_fasta=fasta, db_prefix_path=prefix,
        molecule_type=molecule_type,
    )


# --- saved_blast_database_for ----------------------------------------------

def test_saved_database_for_the_same_fasta_is_found(tmp_path):
    genome = _fasta(tmp_path / "scaffold_5.fasta")
    prefix = tmp_path / "dbs" / "scaffold_5.fasta" / "scaffold_5.fasta"
    _index(prefix)
    db = _db("scaffold_5.fasta", genome, prefix)
    assert dp.saved_blast_database_for(genome, [db]) is db


def test_another_fasta_with_the_same_name_is_not_matched(tmp_path):
    (tmp_path / "a").mkdir()
    genome = _fasta(tmp_path / "a" / "scaffold_5.fasta")
    other = tmp_path / "b"
    other.mkdir()
    other_fasta = _fasta(other / "scaffold_5.fasta")
    prefix = tmp_path / "dbs" / "x" / "x"
    _index(prefix)
    assert dp.saved_blast_database_for(genome, [_db("x", other_fasta, prefix)]) is None


def test_index_older_than_the_fasta_is_not_reused(tmp_path):
    genome = _fasta(tmp_path / "g.fa")
    prefix = tmp_path / "dbs" / "g" / "g"
    _index(prefix, age=3600)
    assert dp.saved_blast_database_for(genome, [_db("g", genome, prefix)]) is None


def test_missing_index_files_are_not_reused(tmp_path):
    genome = _fasta(tmp_path / "g.fa")
    prefix = tmp_path / "dbs" / "g" / "g"
    assert dp.saved_blast_database_for(genome, [_db("g", genome, prefix)]) is None


def test_protein_databases_are_ignored(tmp_path):
    genome = _fasta(tmp_path / "g.fa")
    prefix = tmp_path / "dbs" / "g" / "g"
    _index(prefix, ".pin")
    db = _db("g", genome, prefix, molecule_type="prot")
    assert dp.saved_blast_database_for(genome, [db]) is None


def test_multi_volume_genome_index_counts_as_built(tmp_path):
    # Large genomes are split into volumes; makeblastdb then writes .nal only.
    genome = _fasta(tmp_path / "big.fa")
    prefix = tmp_path / "dbs" / "big" / "big"
    _index(prefix, ".nal")
    db = _db("big", genome, prefix)
    assert dp.blast_database_is_ready(db)
    assert dp.saved_blast_database_for(genome, [db]) is db


# --- unique_blast_database_name --------------------------------------------

def test_unique_name_keeps_a_free_name(tmp_path):
    assert dp.unique_blast_database_name("g.fa", [], tmp_path) == "g.fa"


def test_unique_name_skips_taken_names_and_folders(tmp_path):
    taken = [_db("G.fa", tmp_path / "x.fa", tmp_path / "x")]
    (tmp_path / "g.fa_(2)").mkdir()  # folder left by an older database
    assert dp.unique_blast_database_name("g.fa", taken, tmp_path) == "g.fa (3)"


# --- resolve_blast_database ------------------------------------------------

def _no_build(*args, **kwargs):
    raise AssertionError("makeblastdb must not run")


def test_resolve_reuses_saved_database_without_building(tmp_path, monkeypatch):
    genome = _fasta(tmp_path / "g.fa")
    prefix = tmp_path / "dbs" / "g" / "g"
    _index(prefix)
    monkeypatch.setattr(dp, "build_named_blast_database", _no_build)
    got, created = rq.resolve_blast_database(
        genome, [_db("g", genome, prefix)], save_index=True, base_dir=tmp_path / "dbs"
    )
    assert got == prefix
    assert created is None


def test_resolve_builds_and_returns_a_new_database(tmp_path, monkeypatch):
    genome = _fasta(tmp_path / "g.fa")
    calls = []

    def _build(source, base_dir, name, *, molecule_type, use_system_blast):
        calls.append((source, name, molecule_type, use_system_blast))
        prefix = Path(base_dir) / name / name
        _index(prefix)
        return prefix

    monkeypatch.setattr(dp, "build_named_blast_database", _build)
    got, created = rq.resolve_blast_database(
        genome, [], save_index=True, base_dir=tmp_path / "dbs", use_system_blast=True
    )
    assert calls == [(genome.resolve(), "g.fa", "nucl", True)]
    assert created.name == "g.fa"
    assert created.source_fasta == genome.resolve()
    assert created.db_prefix_path == got
    assert created.created_at  # shown in the Created column


def test_resolve_without_saving_uses_a_temporary_index(tmp_path, monkeypatch):
    genome = _fasta(tmp_path / "g.fa")
    monkeypatch.setattr(dp, "build_named_blast_database", _no_build)
    assert rq.resolve_blast_database(
        genome, [], save_index=False, base_dir=tmp_path / "dbs"
    ) == (None, None)


# --- quantify_library_blastn with a saved database -------------------------

def _capture_blastn(monkeypatch):
    seen = {}

    def _run(cmd, cancel=None, check=True):
        seen["cmd"] = cmd
        Path(cmd[cmd.index("-out") + 1]).write_text("")

    monkeypatch.setattr(rq, "run_subprocess_with_cancel", _run)
    monkeypatch.setattr(rq, "check_database", _no_build)
    return seen


def test_blastn_searches_the_saved_database(tmp_path, monkeypatch):
    lib = _fasta(tmp_path / "lib.fa")
    genome = _fasta(tmp_path / "g.fa")
    prefix = tmp_path / "dbs" / "g" / "g"
    _index(prefix)
    seen = _capture_blastn(monkeypatch)
    rq.quantify_library_blastn(lib, genome, tmp_path / "out", db_prefix=prefix)
    cmd = seen["cmd"]
    assert cmd[cmd.index("-db") + 1] == str(prefix)


def test_blastn_with_a_deleted_saved_database_says_what_to_do(tmp_path, monkeypatch):
    lib = _fasta(tmp_path / "lib.fa")
    genome = _fasta(tmp_path / "g.fa")
    _capture_blastn(monkeypatch)
    with pytest.raises(RuntimeError, match="Add/Manage databases"):
        rq.quantify_library_blastn(
            lib, genome, tmp_path / "out", db_prefix=tmp_path / "gone" / "g"
        )


def test_quantify_library_passes_the_database_and_system_blast(tmp_path, monkeypatch):
    got = {}

    def _fake(library, target, out_dir, **kwargs):
        got.update(kwargs)
        return rq.QuantifyResult(engine="blastn", target_length=1)

    monkeypatch.setattr(rq, "quantify_library_blastn", _fake)
    rq.quantify_library(
        tmp_path / "l.fa", tmp_path / "g.fa", tmp_path, engine="blastn",
        blast_db_prefix="/dbs/g/g", use_system_blast=True,
    )
    assert got["db_prefix"] == "/dbs/g/g"
    assert got["use_system_blast"] is True
