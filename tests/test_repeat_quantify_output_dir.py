"""Repeat quantify output folders: timestamped per run, never overwritten,
with a run_parameters.json record."""

import json
from datetime import datetime
from pathlib import Path

from repeatviewer.backends.repeat_quantify import (
    RUN_PARAMETERS_FILENAME,
    create_repeat_quant_output_dir,
    engine_params,
    repeat_quant_output_dir,
    run_stamp,
    write_run_parameters,
)

STAMP = "2026-09-17_1432"
PARAMS = {
    "blast_word_size": 11,
    "blast_evalue": 0.001,
    "blast_min_identity": 0.0,
    "blast_min_coverage": 0.0,
    "rm_div": 40,
    "rm_sensitivity": "sensitive",
}


def test_run_stamp_format():
    assert run_stamp(datetime(2026, 9, 7, 9, 5)) == "2026-09-07_0905"


def test_default_folder_is_next_to_the_genome(tmp_path):
    genome = tmp_path / "data" / "Csat_v1.fa"
    assert repeat_quant_output_dir(genome, STAMP) == (
        tmp_path / "data" / f"Csat_v1_repeat_quant_{STAMP}"
    )


def test_chosen_folder_gets_the_run_folder_inside(tmp_path):
    genome = tmp_path / "data" / "Csat_v1.fa"
    base = tmp_path / "results"
    assert repeat_quant_output_dir(genome, STAMP, base) == (
        base / f"Csat_v1_repeat_quant_{STAMP}"
    )


def test_create_makes_missing_parents_and_never_reuses_a_folder(tmp_path):
    genome = tmp_path / "Csat_v1.fa"
    base = tmp_path / "new" / "results"
    first = create_repeat_quant_output_dir(genome, STAMP, base)
    second = create_repeat_quant_output_dir(genome, STAMP, base)
    third = create_repeat_quant_output_dir(genome, STAMP, base)
    assert first.name == f"Csat_v1_repeat_quant_{STAMP}"
    assert second.name == f"Csat_v1_repeat_quant_{STAMP}_2"
    assert third.name == f"Csat_v1_repeat_quant_{STAMP}_3"
    assert all(p.is_dir() for p in (first, second, third))


def test_same_genome_name_from_two_folders_gets_two_run_folders(tmp_path):
    base = tmp_path / "results"
    a = create_repeat_quant_output_dir(tmp_path / "a" / "genome.fa", STAMP, base)
    b = create_repeat_quant_output_dir(tmp_path / "b" / "genome.fa", STAMP, base)
    assert a != b


def test_engine_params_keep_only_what_the_engine_uses():
    assert set(engine_params("blastn", PARAMS)) == {
        "blast_word_size", "blast_evalue", "blast_min_identity", "blast_min_coverage",
    }
    assert set(engine_params("repeatmasker", PARAMS)) == {"rm_div", "rm_sensitivity"}
    assert engine_params("both", PARAMS) == PARAMS


def test_run_parameters_json_is_written_and_updated(tmp_path):
    lib = tmp_path / "CroSat.fa"
    genome = tmp_path / "Csat_v1.fa"
    out = create_repeat_quant_output_dir(genome, STAMP)
    started = datetime(2026, 9, 17, 14, 32, 5)

    path = write_run_parameters(
        out, library_fasta=lib, target_fasta=genome, engine="blastn",
        threads=8, params=PARAMS, started=started,
    )
    assert path == out / RUN_PARAMETERS_FILENAME
    record = json.loads(path.read_text())
    assert record["status"] == "running"
    assert record["finished"] is None
    assert record["library"] == str(lib.resolve())
    assert record["genome"] == str(genome.resolve())
    assert record["engine"] == "blastn"
    assert record["threads"] == 8
    assert "rm_div" not in record["parameters"]
    assert record["started"] == "2026-09-17T14:32:05"
    assert record["repeatviewer_version"]

    report = out / "Csat_v1_repeat_quant.html"
    report.write_text("<html></html>")
    write_run_parameters(
        out, library_fasta=lib, target_fasta=genome, engine="blastn",
        threads=8, params=PARAMS, started=started, status="completed",
        finished=datetime(2026, 9, 17, 14, 40), outputs=[report],
    )
    record = json.loads(path.read_text())
    assert record["status"] == "completed"
    assert record["outputs"] == ["Csat_v1_repeat_quant.html"]
    assert "error" not in record
    assert not list(out.glob("*.tmp"))
