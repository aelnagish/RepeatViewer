"""Tools → HMMER: options, dialog layout, and search runner."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import (
    QApplication,
    QGroupBox,
    QLabel,
    QPushButton,
    QRadioButton,
)

from repeatviewer.backends.hmmer_tool import (
    next_hmmer_table_path,
    run_hmmer_search,
    write_hmmer_hits_table,
)
from repeatviewer.job_options import (
    ALIGN_INPUT_DEFINE,
    ALIGN_INPUT_SELECTED,
    DEFAULT_HMMER_EVALUE,
    HMMER_RIGHT_CLICK_DIALOG,
    HmmerOptions,
    load_hmmer_options,
    save_hmmer_options,
)
from repeatviewer.project import DATABASE_TYPE_PROTEIN, UserDatabase
from repeatviewer.views.tool_dialogs import HmmerParamsDialog


def test_hmmer_options_round_trip(monkeypatch) -> None:
    store: dict[str, dict] = {}

    def _save(group: str, values: dict) -> None:
        store[group] = dict(values)

    def _load(group: str) -> dict:
        return dict(store.get(group, {}))

    monkeypatch.setattr("repeatviewer.job_options._save_settings_group", _save)
    monkeypatch.setattr("repeatviewer.job_options._load_settings_group", _load)
    opts = HmmerOptions(
        database_name="Pfam",
        e_value=1e-4,
        domain_e_value=1e-4,
        num_threads=2,
        max_hits=50,
        write_text=False,
        minimum_bit_score=10.0,
        max_orfs_per_sequence=20,
        show_process_window=True,
        close_process_window_when_finished=True,
        right_click_action=HMMER_RIGHT_CLICK_DIALOG,
    )
    save_hmmer_options(opts)
    got = load_hmmer_options()
    assert got.database_name == "Pfam"
    assert got.e_value == 1e-4
    assert got.domain_e_value == 1e-4
    assert got.num_threads == 2
    assert got.max_hits == 50
    assert got.write_text is False
    assert got.minimum_bit_score == 10.0
    assert got.max_orfs_per_sequence == 20
    assert got.show_process_window is True
    assert got.close_process_window_when_finished is True
    assert got.right_click_action == HMMER_RIGHT_CLICK_DIALOG
    assert got.input_source == ALIGN_INPUT_SELECTED
    assert got.input_paths == []


def test_next_hmmer_table_path_avoids_existing(tmp_path: Path) -> None:
    query = tmp_path / "seq.fa"
    query.write_text(">q\nA\n", encoding="utf-8")
    first = next_hmmer_table_path(tmp_path, query)
    assert first.name == "seq.fa_hmmer.txt"
    first.write_text("existing", encoding="utf-8")
    second = next_hmmer_table_path(tmp_path, query)
    assert second.name == "seq.fa_hmmer1.txt"
    second.write_text("existing", encoding="utf-8")
    third = next_hmmer_table_path(tmp_path, query)
    assert third.name == "seq.fa_hmmer2.txt"


def test_write_hmmer_hits_table(tmp_path: Path) -> None:
    out = tmp_path / "q_hmmer.txt"
    write_hmmer_hits_table(
        [(42.5, 1e-10, "orf1", "PF00001", 10, 80)],
        out,
    )
    text = out.read_text(encoding="utf-8")
    assert text.startswith("query_id\thmm_name\tscore\tevalue\tqstart\tqend\n")
    assert "orf1\tPF00001\t42.500\t1e-10\t10\t80\n" in text


def test_run_hmmer_emits_process_log(tmp_path: Path, monkeypatch) -> None:
    query = tmp_path / "query.fa"
    query.write_text(">q\nATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n", encoding="utf-8")
    hmm = tmp_path / "dna.hmm"
    hmm.write_text(
        "HMMER3/f [3.3 | November 2019]\nNAME  test\nALPH  DNA\nLENG  10\n",
        encoding="utf-8",
    )
    db = UserDatabase(name="dna_hmm", path=hmm, db_type="DNA", is_hmm=True)
    monkeypatch.setattr(
        "repeatviewer.backends.hmmer_tool.database_search_readiness",
        lambda _db: None,
    )
    monkeypatch.setattr(
        "repeatviewer.backends.hmmer_tool.collect_pyhmmer_domain_hits",
        lambda *a, **k: [(1.0, 0.1, "q", "test", 1, 5)],
    )
    monkeypatch.setattr(
        "repeatviewer.backends.hmmer_tool.hmm_uses_nucleotide_query",
        lambda _p: True,
    )
    monkeypatch.setattr(
        "repeatviewer.backends.hmmer_tool._pyhmmer_available",
        lambda: True,
    )
    lines: list[str] = []
    from repeatviewer.external_process import process_log_sink

    with process_log_sink(lines.append):
        run_hmmer_search(query, db, options=HmmerOptions(write_text=False))
    joined = "\n".join(lines)
    assert "Running hmmscan" in joined
    assert "domain hit(s) before filtering" in joined
    assert "Done." in joined


def test_run_hmmer_dna_path_writes_table(tmp_path: Path, monkeypatch) -> None:
    query = tmp_path / "query.fa"
    query.write_text(">q\nATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n", encoding="utf-8")
    hmm = tmp_path / "dna.hmm"
    hmm.write_text(
        "HMMER3/f [3.3 | November 2019]\nNAME  test\nALPH  DNA\nLENG  10\n",
        encoding="utf-8",
    )
    db = UserDatabase(
        name="dna_hmm",
        path=hmm,
        db_type="DNA",
        is_hmm=True,
    )
    monkeypatch.setattr(
        "repeatviewer.backends.hmmer_tool.database_search_readiness",
        lambda _db: None,
    )
    hits = [(12.0, 0.001, "q", "test", 1, 20)]

    def _fake_collect(hmm_path, query_path, *, report_e, domain_e=None, cpus=1):
        assert Path(hmm_path) == hmm.resolve()
        assert Path(query_path) == query.resolve()
        assert cpus == 3
        assert report_e == DEFAULT_HMMER_EVALUE
        return hits

    monkeypatch.setattr(
        "repeatviewer.backends.hmmer_tool.collect_pyhmmer_domain_hits",
        _fake_collect,
    )
    monkeypatch.setattr(
        "repeatviewer.backends.hmmer_tool.hmm_uses_nucleotide_query",
        lambda _p: True,
    )
    monkeypatch.setattr(
        "repeatviewer.backends.hmmer_tool._pyhmmer_available",
        lambda: True,
    )
    outs = run_hmmer_search(
        query,
        db,
        options=HmmerOptions(num_threads=3, max_hits=10, write_text=True),
    )
    assert len(outs) == 1
    assert outs[0].name == "query.fa_hmmer.txt"
    assert "test" in outs[0].read_text(encoding="utf-8")


def test_run_hmmer_numbered_table_when_default_exists(
    tmp_path: Path, monkeypatch
) -> None:
    query = tmp_path / "query.fa"
    query.write_text(">q\nATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n", encoding="utf-8")
    (tmp_path / "query.fa_hmmer.txt").write_text("old\n", encoding="utf-8")
    hmm = tmp_path / "dna.hmm"
    hmm.write_text(
        "HMMER3/f [3.3 | November 2019]\nNAME  test\nALPH  DNA\nLENG  10\n",
        encoding="utf-8",
    )
    db = UserDatabase(name="dna_hmm", path=hmm, db_type="DNA", is_hmm=True)
    monkeypatch.setattr(
        "repeatviewer.backends.hmmer_tool.database_search_readiness",
        lambda _db: None,
    )
    monkeypatch.setattr(
        "repeatviewer.backends.hmmer_tool.collect_pyhmmer_domain_hits",
        lambda *a, **k: [(1.0, 0.1, "q", "test", 1, 5)],
    )
    monkeypatch.setattr(
        "repeatviewer.backends.hmmer_tool.hmm_uses_nucleotide_query",
        lambda _p: True,
    )
    monkeypatch.setattr(
        "repeatviewer.backends.hmmer_tool._pyhmmer_available",
        lambda: True,
    )
    outs = run_hmmer_search(query, db, options=HmmerOptions(write_text=True))
    assert len(outs) == 1
    assert outs[0].name == "query.fa_hmmer1.txt"
    assert (tmp_path / "query.fa_hmmer.txt").read_text(encoding="utf-8") == "old\n"


def test_run_hmmer_protein_path_uses_orfs(tmp_path: Path, monkeypatch) -> None:
    query = tmp_path / "te.fa"
    query.write_text(
        ">te\n" + ("ATGCGATCGATC" * 20) + "\n",
        encoding="utf-8",
    )
    hmm = tmp_path / "prot.hmm"
    hmm.write_text(
        "HMMER3/f [3.3 | November 2019]\nNAME  RT\nALPH  amino\nLENG  50\n",
        encoding="utf-8",
    )
    db = UserDatabase(
        name="prot_hmm",
        path=hmm,
        db_type=DATABASE_TYPE_PROTEIN,
        is_hmm=True,
    )
    orf_faa = tmp_path / "orfs.faa"
    orf_faa.write_text(">te_ORF1\nMKLLVVV\n", encoding="utf-8")
    called: dict[str, Path] = {}

    def _fake_orfs(query_path, *, max_orfs_per_sequence, work_dir):
        called["query"] = Path(query_path)
        return orf_faa

    def _fake_collect(hmm_path, query_path, *, report_e, domain_e=None, cpus=1):
        called["search_query"] = Path(query_path)
        return [(8.0, 0.01, "te_ORF1", "RT", 1, 7)]

    monkeypatch.setattr(
        "repeatviewer.backends.hmmer_tool._orf_protein_query", _fake_orfs
    )
    monkeypatch.setattr(
        "repeatviewer.backends.hmmer_tool.collect_pyhmmer_domain_hits",
        _fake_collect,
    )
    monkeypatch.setattr(
        "repeatviewer.backends.hmmer_tool.hmm_uses_nucleotide_query",
        lambda _p: False,
    )
    monkeypatch.setattr(
        "repeatviewer.backends.hmmer_tool._fasta_looks_nucleotide",
        lambda _p: True,
    )
    monkeypatch.setattr(
        "repeatviewer.backends.hmmer_tool._pyhmmer_available",
        lambda: True,
    )
    outs = run_hmmer_search(query, db, options=HmmerOptions(write_text=True))
    assert called["query"] == query.resolve()
    assert called["search_query"] == orf_faa
    assert outs and outs[0].is_file()


def test_hmmer_dialog_sections_and_advanced(qapp: QApplication) -> None:
    dlg = HmmerParamsDialog(initial=HmmerOptions())
    dlg.show()
    qapp.processEvents()
    labels = [lab.text() for lab in dlg.findChildren(QLabel)]
    assert "<b>Inputs</b>" in labels
    assert "<b>Settings</b>" in labels
    radios = [r.text() for r in dlg.findChildren(QRadioButton)]
    assert "Use selected sequences" in radios
    assert "Define inputs" in radios
    texts = [b.text() for b in dlg.findChildren(QPushButton)]
    assert any(t.startswith("Advanced Settings") for t in texts)
    assert "Manage databases…" in texts
    assert dlg.input_source() == ALIGN_INPUT_SELECTED
    assert dlg._hmmer_advanced_panel.isHidden()
    group_titles = [g.title() for g in dlg.findChildren(QGroupBox)]
    assert "Right-click action" in group_titles
    assert "Executables" not in group_titles
    assert "Run HMMER with saved parameters" in radios
    compact_h = dlg.height()
    dlg._advanced_settings_btn.click()
    qapp.processEvents()
    assert dlg._hmmer_advanced_panel.isVisible()
    assert dlg.height() >= compact_h
    dlg._input_define.setChecked(True)
    qapp.processEvents()
    assert dlg.input_source() == ALIGN_INPUT_DEFINE
    assert not dlg._define_inputs.isHidden()
    opts = dlg.options()
    assert opts.write_text is True
    assert opts.e_value == DEFAULT_HMMER_EVALUE
    dlg.close()
