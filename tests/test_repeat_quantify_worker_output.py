"""RepeatQuantifyWorker writes every run to its own folder with a JSON record."""

import json
from pathlib import Path

import pytest

from repeatviewer.backends import repeat_quantify as rq
from repeatviewer.backends.job_cancel import JobCancelledError

mw = pytest.importorskip("repeatviewer.main_window")

STAMP = "2026-09-17_1432"


@pytest.fixture
def fastas(tmp_path):
    lib = tmp_path / "lib.fa"
    genome = tmp_path / "genome.fa"
    for f in (lib, genome):
        f.write_text(">s\nACGT\n")
    return lib, genome


def _fake_quantify(library, target, out_dir, **kwargs):
    (Path(out_dir) / f"{Path(target).name}_repeat_quant.tsv").write_text("x\n")
    return []


def _fake_summary(results, path, **kwargs):
    Path(path).write_text("<html></html>")
    return Path(path)


def _run(worker):
    got = {}
    worker.finished.connect(lambda out: got.setdefault("finished", out))
    worker.failed.connect(lambda msg: got.setdefault("failed", msg))
    worker.run()
    return got


def test_repeat_runs_get_separate_folders(qapp, fastas, monkeypatch):
    lib, genome = fastas
    monkeypatch.setattr(rq, "quantify_library", _fake_quantify)
    monkeypatch.setattr(rq, "write_summary_html", _fake_summary)
    for _ in range(2):
        got = _run(mw.RepeatQuantifyWorker(
            lib, [genome], engine="blastn", threads=2,
            params={"blast_word_size": 11, "rm_div": 40}, stamp=STAMP,
        ))
        assert "failed" not in got
    first = genome.parent / f"genome.fa_repeat_quant_{STAMP}"
    second = genome.parent / f"genome.fa_repeat_quant_{STAMP}_2"
    for folder in (first, second):
        record = json.loads((folder / "run_parameters.json").read_text())
        assert record["status"] == "completed"
        assert record["parameters"] == {"blast_word_size": 11}
        assert sorted(record["outputs"]) == [
            "genome.fa_repeat_quant.html", "genome.fa_repeat_quant.tsv",
        ]
    assert got["finished"] == [str(second / "genome.fa_repeat_quant.html")]


def test_chosen_output_base_is_used(qapp, fastas, tmp_path, monkeypatch):
    lib, genome = fastas
    monkeypatch.setattr(rq, "quantify_library", _fake_quantify)
    monkeypatch.setattr(rq, "write_summary_html", _fake_summary)
    base = tmp_path / "results"
    _run(mw.RepeatQuantifyWorker(
        lib, [genome], engine="blastn", threads=1, stamp=STAMP, output_base=base,
    ))
    assert (base / f"genome.fa_repeat_quant_{STAMP}" / "run_parameters.json").is_file()


@pytest.mark.parametrize(
    "error, status", [(RuntimeError("boom"), "failed"), (JobCancelledError(), "cancelled")]
)
def test_failed_or_cancelled_runs_are_recorded(qapp, fastas, monkeypatch, error, status):
    lib, genome = fastas

    def _raise(*a, **k):
        raise error

    monkeypatch.setattr(rq, "quantify_library", _raise)
    _run(mw.RepeatQuantifyWorker(lib, [genome], engine="blastn", threads=1, stamp=STAMP))
    record = json.loads(
        (genome.parent / f"genome.fa_repeat_quant_{STAMP}" / "run_parameters.json").read_text()
    )
    assert record["status"] == status
    assert record["finished"]
    if status == "failed":
        assert "boom" in record["error"]
