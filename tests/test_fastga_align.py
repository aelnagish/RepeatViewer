"""FastGA in the Align window: modes (self / all pairs / against a reference),
reference index built once, version-aware soft-mask flag, option migration,
and the FastGA tab of the Alignment window."""

from __future__ import annotations

import stat
import sys
from pathlib import Path

import pytest
from PySide6.QtWidgets import QPushButton, QTabWidget

from repeatviewer.backends import fastga_dotplot as fg
from repeatviewer import job_options
from repeatviewer.job_options import (
    ALIGNER_FASTGA,
    ALIGNER_MAFFT,
    ALIGNER_MUSCLE,
    AlignOptions,
    FASTGA_MODE_PAIRS,
    FASTGA_MODE_REFERENCE,
    FASTGA_MODE_SELF,
    FastGaOptions,
)

USAGE_15 = """Usage: FastGA [-vkMS] [-L:<log:path>] [-T<int(8)>] [-P<dir($TMPDIR)>] [<format(-paf)>]
              [-f<int(10)>] [-c<int(85)> [-s<int(1000)>] [-l<int(100)>] [-i<float(.7)]
                <source1:path>[<precursor>] (#[<mask>[.1bed]])*
      -M: Use GIX soft masks (if no # arguments).
"""
USAGE_13 = """Usage: FastGA [-vkNS] [-L:<log:path>] [-T<int(8)>] [-P<dir($TMPDIR)>] [<format(-paf)>]
                <source1:path>[<precursor>] [ <source2:path>[<precursor>] ]
      -N: Do not soft mask.
"""


# ------------------------------------------------------------ soft masking


def test_mask_style_is_read_from_the_usage_text():
    assert fg.mask_style_from_usage(USAGE_15) == fg.MASK_STYLE_ON_FLAG
    assert fg.mask_style_from_usage(USAGE_13) == fg.MASK_STYLE_OFF_FLAG
    assert fg.mask_style_from_usage("Usage: FastGA [-vk] <source1:path>") is None


def test_expand_fastga_input_paths_files_and_folders(tmp_path):
    folder = tmp_path / "genomes"
    folder.mkdir()
    a = folder / "a.fa"
    b = folder / "b.fasta"
    a.write_text(">a\nACGT\n")
    b.write_text(">b\nACGT\n")
    (folder / "notes.txt").write_text("skip")
    lone = tmp_path / "lone.fa"
    lone.write_text(">l\nAAAA\n")
    got = fg.expand_fastga_input_paths([lone, folder])
    assert got == [lone.resolve(), a.resolve(), b.resolve()]


def test_soft_mask_flags_per_version():
    on, off = fg.MASK_STYLE_ON_FLAG, fg.MASK_STYLE_OFF_FLAG
    # FastGA 1.5+: -N would be "an illegal option" -- never sent.
    assert fg.soft_mask_flags(True, on) == ["-M"]
    assert fg.soft_mask_flags(False, on) == []
    # FastGA <= 1.4: masks by default, -N switches it off.
    assert fg.soft_mask_flags(True, off) == []
    assert fg.soft_mask_flags(False, off) == ["-N"]
    # None = leave FastGA's own default; unknown version = send nothing.
    assert fg.soft_mask_flags(None, on) == []
    assert fg.soft_mask_flags(False, None) == []


# ------------------------------------------------------- names and summary


def test_reference_output_names_keep_the_full_input_name(tmp_path):
    q = tmp_path / "Csat_chr1.fa"
    ref = tmp_path / "genome.fasta"
    assert fg.reference_output_path(q, ref, fg.FORMAT_PAF).name == (
        "Csat_chr1.fa_FastGA_vs_genome.fasta.paf"
    )
    assert fg.reference_output_path(q, ref, fg.FORMAT_ALN).name.endswith(".1aln")
    assert fg.reference_combined_output_path(q, ref).name == (
        "Csat_chr1.fa_multi_FastGA_vs_genome.fasta.paf"
    )


def test_reference_is_never_its_own_query(tmp_path):
    a, b, ref = (tmp_path / n for n in ("a.fa", "b.fa", "ref.fa"))
    for p in (a, b, ref):
        p.write_text(">x\nACGT\n")
    assert fg.reference_queries([a, ref, b, a], ref) == [a, b]
    assert fg.reference_queries([ref], ref) == []


def test_summary_line(tmp_path):
    a, b, ref = (tmp_path / n for n in ("a.fa", "b.fa", "ref.fa"))
    s = fg.fastga_mode_summary
    assert s([a, b], FASTGA_MODE_SELF) == "2 genomes selected → 2 self-comparisons."
    assert "3 comparisons" in s([a, b, ref], FASTGA_MODE_PAIRS)
    assert "needs two" in s([a], FASTGA_MODE_PAIRS)
    assert "Choose the reference" in s([a, b], FASTGA_MODE_REFERENCE, None)
    assert s([a, b, ref], FASTGA_MODE_REFERENCE, ref) == (
        "2 query genomes → 2 runs against ref.fa (one result file per query)."
    )
    assert "into one file" in s([a, b], FASTGA_MODE_REFERENCE, ref, combined=True)
    assert "Only the reference" in s([ref], FASTGA_MODE_REFERENCE, ref)


def test_step_count_includes_the_index_build(tmp_path):
    a, b, ref = (tmp_path / n for n in ("a.fa", "b.fa", "ref.fa"))
    o = FastGaOptions(mode=FASTGA_MODE_REFERENCE, reference_fasta=str(ref))
    assert fg.fastga_step_count([a, b, ref], o) == 3
    assert fg.fastga_step_count([a, b, ref], FastGaOptions(mode=FASTGA_MODE_PAIRS)) == 3
    assert fg.fastga_step_count([a, b], FastGaOptions()) == 2


# ----------------------------------------------------------- saved options


def test_old_pairwise_setting_becomes_the_pairs_mode(monkeypatch):
    monkeypatch.setattr(job_options, "_load_settings_group", lambda g: {"pairwise": "true"})
    assert job_options.load_fastga_options().mode == FASTGA_MODE_PAIRS
    monkeypatch.setattr(job_options, "_load_settings_group", lambda g: {})
    assert job_options.load_fastga_options().mode == FASTGA_MODE_SELF
    monkeypatch.setattr(
        job_options, "_load_settings_group",
        lambda g: {"mode": "reference", "pairwise": "true", "reference_fasta": "/r.fa"},
    )
    loaded = job_options.load_fastga_options()
    assert (loaded.mode, loaded.reference_fasta) == (FASTGA_MODE_REFERENCE, "/r.fa")


def test_saved_options_keep_temp_dir_and_write_pairwise_for_old_builds(monkeypatch):
    saved = {}
    monkeypatch.setattr(job_options, "_save_settings_group", lambda g, v: saved.update(v))
    job_options.save_fastga_options(
        FastGaOptions(mode=FASTGA_MODE_PAIRS, temp_dir="/scratch/fastga")
    )
    assert saved["temp_dir"] == "/scratch/fastga"
    assert saved["pairwise"] is True and saved["mode"] == FASTGA_MODE_PAIRS


def test_fastga_can_be_the_right_click_aligner_but_never_the_msa_one(monkeypatch):
    monkeypatch.setattr(
        job_options, "_load_settings_group",
        lambda g: {"right_click_aligner": "fastga", "combined_aligner": "fastga"},
    )
    loaded = job_options.load_align_options()
    assert loaded.right_click_aligner == ALIGNER_FASTGA
    assert loaded.combined_aligner == ALIGNER_MAFFT


# ------------------------------------------ batch run with a stand-in FastGA

_FAKE_FASTGA = r"""#!/bin/sh
[ $# -eq 0 ] && { echo 'Usage: FastGA [-vkMS] [-T<int(8)>] <source1:path>' >&2; exit 1; }
echo "FastGA $*" >> "$FAKE_LOG"
src=""; out1=""
for a in "$@"; do
  case "$a" in
    -1:*) out1="${a#-1:}" ;;
    -*) ;;
    *) src="$src $(basename "$a")" ;;
  esac
done
if [ -n "$out1" ]; then echo "aln$src" > "$out1"; else echo "paf$src"; fi
"""
_FAKE_GIXMAKE = r"""#!/bin/sh
echo "GIXmake $*" >> "$FAKE_LOG"
eval core=\${$#}
touch "$core.gix" "$core.1gdb" "$(dirname "$core")/.$(basename "$core").ktab.1"
"""


@pytest.fixture
def fake_tools(tmp_path, monkeypatch):
    if sys.platform == "win32":
        pytest.skip("shell stand-ins")
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    for name, body in (("FastGA", _FAKE_FASTGA), ("GIXmake", _FAKE_GIXMAKE)):
        exe = bin_dir / name
        exe.write_text(body)
        exe.chmod(exe.stat().st_mode | stat.S_IEXEC)
    log = tmp_path / "calls.log"
    monkeypatch.setenv("FAKE_LOG", str(log))
    monkeypatch.setattr(fg, "fastga_executable", lambda: str(bin_dir / "FastGA"))
    monkeypatch.setattr(fg, "gixmake_executable", lambda: str(bin_dir / "GIXmake"))
    fg._mask_style_cache.clear()
    data = tmp_path / "data"
    data.mkdir()
    files = []
    for name in ("q1.fa", "q2.fa", "ref.fa"):
        (data / name).write_text(">x\nACGTACGT\n")
        files.append(data / name)
    return files, log


def _calls(log: Path, tool: str) -> list[str]:
    return [l for l in log.read_text().splitlines() if l.startswith(tool + " ")]


def test_reference_mode_indexes_the_reference_once(fake_tools):
    (q1, q2, ref), log = fake_tools
    opts = FastGaOptions(mode=FASTGA_MODE_REFERENCE, reference_fasta=str(ref),
                         threads=2, keep_index=True)
    steps = []
    out = fg.run_fastga_batch([q1, q2, ref], opts, on_step=steps.append)
    assert [p.name for p in out] == [
        "q1.fa_FastGA_vs_ref.fa.paf", "q2.fa_FastGA_vs_ref.fa.paf",
    ]
    assert out[0].read_text().strip() == "paf q1.fa ref_fa.gix"
    assert steps == [1, 2, 3]
    assert len(_calls(log, "GIXmake")) == 1
    runs = _calls(log, "FastGA")
    assert len(runs) == 2 and all("-k" not in r.split() for r in runs)
    # Kept index sits in a hidden folder (not listed in Files) and is reused.
    assert (ref.parent / ".repeatviewer_fastga" / "ref_fa.gix").is_file()
    fg.run_fastga_batch([q1], opts)
    assert len(_calls(log, "GIXmake")) == 1


def test_reference_index_is_removed_when_not_kept(fake_tools):
    (q1, _q2, ref), log = fake_tools
    opts = FastGaOptions(mode=FASTGA_MODE_REFERENCE, reference_fasta=str(ref),
                         threads=2, keep_index=False)
    fg.run_fastga_batch([q1], opts)
    assert len(_calls(log, "GIXmake")) == 1
    assert sorted(p.name for p in ref.parent.iterdir()) == [
        "q1.fa", "q1.fa_FastGA_vs_ref.fa.paf", "q2.fa", "ref.fa",
    ]


def test_reference_mode_into_one_file_and_1aln_falls_back_to_paf(fake_tools):
    (q1, q2, ref), _log = fake_tools
    opts = FastGaOptions(mode=FASTGA_MODE_REFERENCE, reference_fasta=str(ref),
                         reference_combined=True, output_format=fg.FORMAT_ALN)
    (out,) = fg.run_fastga_batch([q1, q2], opts)
    assert out.name == "q1.fa_multi_FastGA_vs_ref.fa.paf"
    assert out.read_text().splitlines() == ["paf q1.fa ref_fa.gix", "paf q2.fa ref_fa.gix"]
    assert not list(ref.parent.glob(".repeatviewer_fastga_part_*"))


def test_soft_mask_skips_the_prebuilt_index_and_sends_M(fake_tools):
    (q1, _q2, ref), log = fake_tools
    opts = FastGaOptions(mode=FASTGA_MODE_REFERENCE, reference_fasta=str(ref),
                         soft_mask=True)
    fg.run_fastga_batch([q1], opts)
    assert _calls(log, "GIXmake") == []
    (run,) = _calls(log, "FastGA")
    assert "-M" in run.split() and "-N" not in run.split()
    assert run.endswith(str(ref))


def test_self_and_pairs_modes_are_unchanged(fake_tools):
    (q1, q2, _ref), _log = fake_tools
    out = fg.run_fastga_batch([q1, q2], FastGaOptions(mode=FASTGA_MODE_SELF))
    assert [p.name for p in out] == ["q1.fa_FastGA.paf", "q2.fa_FastGA.paf"]
    (pairs,) = fg.run_fastga_batch([q1, q2], FastGaOptions(mode=FASTGA_MODE_PAIRS))
    assert pairs.name == "q1.fa_multi_FastGA.paf"
    assert pairs.read_text().strip() == "paf q1.fa q2.fa"


def test_reference_mode_errors(fake_tools):
    (q1, _q2, ref), _log = fake_tools
    with pytest.raises(fg.FastgaError, match="Choose a reference"):
        fg.run_fastga_batch([q1], FastGaOptions(mode=FASTGA_MODE_REFERENCE))
    with pytest.raises(fg.FastgaError, match="reference itself is skipped"):
        fg.run_fastga_batch(
            [ref], FastGaOptions(mode=FASTGA_MODE_REFERENCE, reference_fasta=str(ref))
        )


# ------------------------------------------------------ the Alignment window


def _align_dialog(qapp, tmp_path, **kw):
    from repeatviewer.views.tool_dialogs import AlignParamsDialog

    files = []
    for name in ("q1.fa", "q2.fa", "ref.fa"):
        (tmp_path / name).write_text(">x\nACGT\n")
        files.append(tmp_path / name)
    kw.setdefault("fastga_files", files)
    return AlignParamsDialog(**kw), files


@pytest.fixture(autouse=True)
def _silence_fastga_reminder(monkeypatch):
    """FastGA tab open schedules an install reminder; keep tests non-interactive."""
    monkeypatch.setattr(
        "repeatviewer.views.tool_dialogs.AlignParamsDialog._remind_fastga_if_needed",
        lambda self: None,
    )

def test_align_window_has_a_fastga_tab_and_follows_the_rules(qapp, tmp_path):
    dlg, _files = _align_dialog(qapp, tmp_path)
    tabs = dlg.findChild(QTabWidget)
    assert [tabs.tabText(i) for i in range(tabs.count())] == [
        "MAFFT", "MUSCLE", "FastGA", "Options",
    ]
    texts = [b.text() for b in dlg.findChildren(QPushButton)]
    assert [t for t in texts if t in ("Cancel", "Save", "Run")] == ["Cancel", "Save", "Run"]
    assert "OK" not in texts
    # Run follows the tab shown -- by widget, not by position.
    for index, expected in ((0, ALIGNER_MAFFT), (1, ALIGNER_MUSCLE), (2, ALIGNER_FASTGA)):
        tabs.setCurrentIndex(index)
        assert dlg.selected_aligner() == expected
    tabs.setCurrentIndex(3)
    assert dlg.selected_aligner() == ALIGNER_MAFFT  # right-click choice


def test_mafft_and_muscle_define_inputs_rows_keep_height_after_add(qapp, tmp_path):
    """Add another input must grow the dialog — not squeeze path/Browse rows."""
    from PySide6.QtCore import QEventLoop, QTimer

    dlg, _files = _align_dialog(qapp, tmp_path)
    dlg.show()
    qapp.processEvents()
    tabs = dlg.findChild(QTabWidget)
    fastga = dlg.fastga_panel()

    for tab_index, define_radio, widget, resize in (
        (
            0,
            dlg._mafft_input_define,
            dlg._mafft_define_inputs,
            dlg._resize_for_mafft_advanced,
        ),
        (
            1,
            dlg._muscle_input_define,
            dlg._muscle_define_inputs,
            dlg._resize_for_muscle_advanced,
        ),
        (
            2,
            fastga._input_define,
            fastga._define_inputs,
            dlg._resize_for_fastga_advanced,
        ),
    ):
        tabs.setCurrentIndex(tab_index)
        define_radio.setChecked(True)
        qapp.processEvents()
        resize(False)
        qapp.processEvents()
        assert len(widget._rows) == 1
        before = widget._rows[0]._path.height()
        assert before >= 20
        height_before = dlg.height()
        widget.add_row(focus=True)
        # Drain deferred refit (paths_changed → QTimer.singleShot(0, …)).
        loop = QEventLoop()
        QTimer.singleShot(0, loop.quit)
        loop.exec()
        qapp.processEvents()
        assert len(widget._rows) == 2
        h0 = widget._rows[0]._path.height()
        h1 = widget._rows[1]._path.height()
        assert h0 >= before - 2, (tab_index, h0, before)
        assert h1 >= before - 2, (tab_index, h1, before)
        assert abs(h0 - h1) <= 2
        assert dlg.height() >= height_before, (tab_index, dlg.height(), height_before)


def test_align_dialog_reopens_last_tab(qapp, tmp_path):
    from repeatviewer.job_options import AlignOptions, DIALOG_TAB_OPTIONS

    dlg, _files = _align_dialog(
        qapp, tmp_path, initial=AlignOptions(last_dialog_tab=ALIGNER_MUSCLE)
    )
    assert dlg.current_dialog_tab() == ALIGNER_MUSCLE
    assert dlg.selected_aligner() == ALIGNER_MUSCLE
    assert dlg.options().last_dialog_tab == ALIGNER_MUSCLE

    dlg2, _f2 = _align_dialog(
        qapp, tmp_path, initial=AlignOptions(last_dialog_tab=DIALOG_TAB_OPTIONS)
    )
    assert dlg2.current_dialog_tab() == ALIGNER_MAFFT


def test_mafft_add_to_msa_mode_shows_two_paths_and_prefills_msa(qapp, tmp_path):
    from repeatviewer.job_options import (
        ALIGN_RUN_MODE_ADD_TO_MSA,
        MAFFT_ADD_METHOD_ADDFRAGMENTS,
        AlignOptions,
    )
    from repeatviewer.views.tool_dialogs import TOOL_PARAMS_DIALOG_WIDTH_WIDE

    msa = tmp_path / "msa.fa"
    msa.write_text(">a\nAC--GT\n>b\nACGTGT\n")
    seqs = tmp_path / "extra.fa"
    seqs.write_text(">c\nACGT\n")
    dlg, _files = _align_dialog(
        qapp,
        tmp_path,
        initial=AlignOptions(
            mafft_run_mode=ALIGN_RUN_MODE_ADD_TO_MSA,
            mafft_add_method=MAFFT_ADD_METHOD_ADDFRAGMENTS,
            mafft_add_keeplength=True,
        ),
        default_msa=msa,
    )
    dlg.show()
    qapp.processEvents()
    assert dlg.mafft_run_mode() == ALIGN_RUN_MODE_ADD_TO_MSA
    assert dlg._mafft_add_msa.path() == msa.resolve()
    assert dlg._mafft_add_sequences.isVisible()
    assert dlg._mafft_add_sequences._extra_btn is not None
    assert dlg._mafft_add_sequences._extra_btn.text() == "Paste…"
    assert dlg._mafft_define_inputs.isHidden()
    assert dlg._align_needs_wide_inputs() is True
    dlg._resize_for_mafft_advanced(False)
    qapp.processEvents()
    assert dlg.width() == TOOL_PARAMS_DIALOG_WIDTH_WIDE

    assert dlg.apply_path_from_panel(seqs) is True
    got_seq, got_msa = dlg.mafft_add_paths()
    assert got_seq == seqs.resolve()
    assert got_msa == msa.resolve()

    opts = dlg.options()
    assert opts.mafft_run_mode == ALIGN_RUN_MODE_ADD_TO_MSA
    assert opts.mafft_combined is False
    assert opts.mafft_add_method == MAFFT_ADD_METHOD_ADDFRAGMENTS
    assert opts.mafft_add_keeplength is True


def test_paste_sequence_to_add_writes_fasta_beside_msa(qapp, tmp_path, monkeypatch):
    from repeatviewer.job_options import ALIGN_RUN_MODE_ADD_TO_MSA, AlignOptions
    from repeatviewer.views.tool_dialogs import NewSequenceDialog

    msa = tmp_path / "msa.fa"
    msa.write_text(">a\nAC--GT\n>b\nACGTGT\n")
    dlg, _files = _align_dialog(
        qapp,
        tmp_path,
        initial=AlignOptions(mafft_run_mode=ALIGN_RUN_MODE_ADD_TO_MSA),
        default_msa=msa,
    )

    class _FakePaste(NewSequenceDialog):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            assert kwargs.get("show_metadata") is False
            assert not self._form_host.isVisibleTo(self)
            assert not self._summary.isVisibleTo(self)

        def exec(self):  # noqa: A003
            return NewSequenceDialog.DialogCode.Accepted

        def values(self):
            return {
                "name": "pasted_seq",
                "sequence": "ACGTACGT",
                "description": "",
                "organism": "",
                "molecule": "DNA",
            }

    monkeypatch.setattr(
        "repeatviewer.views.tool_dialogs.NewSequenceDialog",
        _FakePaste,
    )
    dlg._paste_sequence_to_add(dlg._mafft_add_sequences)
    got_seq, got_msa = dlg.mafft_add_paths()
    assert got_msa == msa.resolve()
    assert got_seq is not None
    assert got_seq.parent == msa.parent
    assert got_seq.read_text().startswith(">pasted_seq\n")
    assert "ACGTACGT" in got_seq.read_text().replace("\n", "")


def test_sequences_to_add_dialog_hides_metadata_fields(qapp):
    from repeatviewer.views.tool_dialogs import NewSequenceDialog

    dlg = NewSequenceDialog(title="Sequences to add", show_metadata=False)
    assert dlg._form_host.isHidden()
    assert dlg._summary.isHidden()
    assert dlg.width() == NewSequenceDialog._PASTE_DEFAULT_WIDTH
    assert dlg.height() == NewSequenceDialog._PASTE_DEFAULT_HEIGHT
    dlg._text.setPlainText(">my_seq desc\nACGT\n")
    qapp.processEvents()
    assert dlg._buttons.button(dlg._buttons.StandardButton.Ok).isEnabled()
    values = dlg.values()
    assert values["name"] == "my_seq"
    assert values["description"] == ""
    assert values["organism"] == ""
    dlg.close()


def test_muscle_add_to_msa_mode_mirrors_mafft_inputs(qapp, tmp_path):
    from repeatviewer.job_options import ALIGN_RUN_MODE_ADD_TO_MSA, AlignOptions

    msa = tmp_path / "aln.fa"
    msa.write_text(">a\nAAAA\n>b\nAA-A\n")
    dlg, _files = _align_dialog(
        qapp,
        tmp_path,
        initial=AlignOptions(muscle_run_mode=ALIGN_RUN_MODE_ADD_TO_MSA),
        initial_tab=ALIGNER_MUSCLE,
        default_msa=msa,
    )
    assert dlg.muscle_run_mode() == ALIGN_RUN_MODE_ADD_TO_MSA
    assert dlg._muscle_add_msa.path() == msa.resolve()
    assert dlg._muscle_input_selected.isHidden()
    assert dlg.options().muscle_run_mode == ALIGN_RUN_MODE_ADD_TO_MSA


def test_fastga_reference_path_matches_define_inputs_width(qapp, tmp_path):
    """Reference genome path box / Browse must match Define-inputs width."""
    from repeatviewer.job_options import FASTGA_INPUT_DEFINE
    from repeatviewer.views.tool_dialogs import TOOL_PARAMS_DIALOG_WIDTH_WIDE

    ref = tmp_path / "ref.fa"
    ref.write_text(">x\nACGT\n")
    dlg, _files = _align_dialog(
        qapp,
        tmp_path,
        initial_tab=ALIGNER_FASTGA,
        fastga_initial=FastGaOptions(
            mode=FASTGA_MODE_REFERENCE, reference_fasta=str(ref)
        ),
    )
    dlg.show()
    qapp.processEvents()
    panel = dlg.fastga_panel()
    panel._mode_buttons[FASTGA_MODE_REFERENCE].setChecked(True)
    panel._input_define.setChecked(True)
    panel._sync()
    panel._sync_field_column()
    qapp.processEvents()
    assert panel.input_source() == FASTGA_INPUT_DEFINE
    define_row = panel._define_inputs._rows[0]._row
    ref_path = panel._reference.line_edit()
    define_path = define_row.line_edit()
    assert abs(ref_path.width() - define_path.width()) <= 2
    ref_browse = panel._reference._browse_btn
    define_browse = define_row._browse_btn
    assert abs(
        ref_browse.mapTo(dlg, ref_browse.rect().topLeft()).x()
        - define_browse.mapTo(dlg, define_browse.rect().topLeft()).x()
    ) <= 2


def test_fastga_reference_mode_widens_dialog_like_define_inputs(qapp, tmp_path):
    """Compare against a reference uses the same 720 width as Define inputs."""
    from repeatviewer.job_options import FASTGA_MODE_SELF
    from repeatviewer.views.tool_dialogs import (
        TOOL_PARAMS_DIALOG_WIDTH,
        TOOL_PARAMS_DIALOG_WIDTH_WIDE,
    )

    dlg, _files = _align_dialog(qapp, tmp_path, initial_tab=ALIGNER_FASTGA)
    dlg.show()
    qapp.processEvents()
    panel = dlg.fastga_panel()
    panel._mode_buttons[FASTGA_MODE_SELF].setChecked(True)
    panel._input_selected.setChecked(True)
    panel._sync()
    qapp.processEvents()
    dlg._resize_for_fastga_advanced(False)
    qapp.processEvents()
    assert dlg.width() == TOOL_PARAMS_DIALOG_WIDTH

    panel._mode_buttons[FASTGA_MODE_REFERENCE].setChecked(True)
    panel._sync()
    qapp.processEvents()
    dlg._resize_for_fastga_advanced(False)
    qapp.processEvents()
    assert dlg.width() == TOOL_PARAMS_DIALOG_WIDTH_WIDE
    assert dlg._align_needs_wide_inputs() is True


def test_files_panel_click_fills_the_reference_only_in_reference_mode(qapp, tmp_path):
    dlg, (q1, _q2, ref) = _align_dialog(
        qapp, tmp_path, initial_tab=ALIGNER_FASTGA,
        fastga_initial=FastGaOptions(mode=FASTGA_MODE_SELF),
    )
    panel = dlg.fastga_panel()
    assert dlg.apply_path_from_panel(ref) is False
    panel._mode_buttons[FASTGA_MODE_REFERENCE].setChecked(True)
    assert dlg.apply_path_from_panel(tmp_path / "notes.txt") is False
    assert dlg.apply_path_from_panel(ref) is True
    assert panel.reference() == ref
    # On another tab a panel click is not taken.
    dlg.findChild(QTabWidget).setCurrentIndex(0)
    assert dlg.apply_path_from_panel(q1) is False


def test_define_inputs_accepts_panel_files_and_folders(qapp, tmp_path):
    from repeatviewer.job_options import FASTGA_INPUT_DEFINE

    dlg, (q1, q2, ref) = _align_dialog(
        qapp,
        tmp_path,
        initial_tab=ALIGNER_FASTGA,
        fastga_initial=FastGaOptions(),
    )
    panel = dlg.fastga_panel()
    # Define inputs is for this session only; the dialog always opens on Use selected.
    assert panel.input_source() != FASTGA_INPUT_DEFINE
    panel._input_define.setChecked(True)
    assert panel.input_source() == FASTGA_INPUT_DEFINE
    assert not panel._define_inputs.isHidden()
    assert dlg.apply_path_from_panel(q1) is True
    # Next folder click should fill a new row, not overwrite the focused file row.
    panel._define_inputs.add_row(focus=True)
    assert dlg.apply_folder_from_panel(tmp_path) is True
    opts = dlg.fastga_options()
    assert opts.input_source == FASTGA_INPUT_DEFINE
    assert Path(opts.input_paths[0]) == q1.resolve()
    assert Path(opts.input_paths[1]) == tmp_path.resolve()
    assert q2.resolve() in dlg.fastga_run_files()
    assert ref.resolve() in dlg.fastga_run_files()


def test_fastga_options_round_trip_and_fastga_never_becomes_the_msa_aligner(qapp, tmp_path):
    ref = tmp_path / "ref.fa"
    dlg, _files = _align_dialog(
        qapp, tmp_path,
        initial=AlignOptions(combined_aligner=ALIGNER_MUSCLE),
        fastga_initial=FastGaOptions(
            mode=FASTGA_MODE_REFERENCE, reference_fasta=str(ref),
            reference_combined=True, soft_mask=True,
        ),
    )
    dlg._right_click_buttons[ALIGNER_FASTGA].setChecked(True)
    align = dlg.options()
    assert align.right_click_aligner == ALIGNER_FASTGA
    assert align.combined_aligner == ALIGNER_MUSCLE
    got = dlg.fastga_options()
    assert (got.mode, Path(got.reference_fasta), got.reference_combined, got.soft_mask) == (
        FASTGA_MODE_REFERENCE, ref, True, True,
    )
    assert got.temp_dir == ""
    dlg._restore_defaults()
    d = dlg.fastga_options()
    assert (d.mode, d.soft_mask, d.temp_dir) == (FASTGA_MODE_SELF, False, "")
    assert d.symmetric is False and d.seed_freq == 10
    assert Path(d.reference_fasta) == ref  # kept when no setting genome is supplied


def test_fastga_reference_defaults_to_setting_genome(qapp, tmp_path):
    genome = tmp_path / "setting_genome.fa"
    genome.write_text(">g\nACGTACGT\n")
    other = tmp_path / "other.fa"
    other.write_text(">o\nTTTT\n")
    dlg, _files = _align_dialog(
        qapp,
        tmp_path,
        fastga_initial=FastGaOptions(mode=FASTGA_MODE_REFERENCE),
        fastga_default_reference=genome,
    )
    panel = dlg.fastga_panel()
    assert panel.reference() == genome.resolve()
    assert Path(dlg.fastga_options().reference_fasta) == genome.resolve()

    panel._reference.set_path(other)
    dlg._restore_defaults()
    assert panel.reference() == genome.resolve()


def test_fastga_saved_reference_is_not_overridden_by_setting_genome(qapp, tmp_path):
    genome = tmp_path / "setting_genome.fa"
    genome.write_text(">g\nACGTACGT\n")
    saved = tmp_path / "saved_ref.fa"
    saved.write_text(">s\nAAAA\n")
    dlg, _files = _align_dialog(
        qapp,
        tmp_path,
        fastga_initial=FastGaOptions(
            mode=FASTGA_MODE_REFERENCE, reference_fasta=str(saved)
        ),
        fastga_default_reference=genome,
    )
    assert dlg.fastga_panel().reference() == saved.resolve()


def test_every_fastga_field_has_a_help_button(qapp, tmp_path):
    from PySide6.QtWidgets import QToolButton

    dlg, _files = _align_dialog(qapp, tmp_path)
    helps = [
        b for b in dlg.fastga_panel().findChildren(QToolButton) if b.text() == "?"
    ]
    assert len(helps) >= 6


def test_old_seed_freq_default_migrates_to_fastga_default(monkeypatch):
    monkeypatch.setattr(
        job_options, "_load_settings_group", lambda g: {"seed_freq": "200"}
    )
    assert job_options.load_fastga_options().seed_freq == 10


def test_fastga_advanced_settings_button_toggles_options(qapp, tmp_path):
    dlg, _files = _align_dialog(qapp, tmp_path, initial_tab=ALIGNER_FASTGA)
    dlg.show()
    qapp.processEvents()
    panel = dlg.fastga_panel()
    btn = dlg._advanced_settings_btn
    assert not btn.isHidden()
    assert btn.text() == "Advanced Settings ▸"
    assert panel.advanced_visible() is False
    dlg._shrink_to_fastga_content()
    compact_h = dlg.height()
    btn.click()
    qapp.processEvents()
    assert btn.text() == "Advanced Settings ▾"
    assert panel.advanced_visible() is True
    assert panel._advanced_panel.isVisible()
    assert panel._identity.isVisible()
    assert dlg.height() > compact_h
    btn.click()
    qapp.processEvents()
    assert panel.advanced_visible() is False
    assert dlg.height() <= compact_h + 20
    dlg.close()


def test_fastga_reference_mode_refits_dialog_height(qapp, tmp_path):
    """Compare-against-reference must grow for the reference rows, without a blank gap."""
    dlg, _files = _align_dialog(
        qapp,
        tmp_path,
        initial_tab=ALIGNER_FASTGA,
        fastga_initial=FastGaOptions(mode=FASTGA_MODE_SELF),
    )
    dlg.show()
    qapp.processEvents()
    dlg._fit_current_align_tab()
    qapp.processEvents()
    panel = dlg.fastga_panel()
    h_self = dlg.height()
    panel._mode_buttons[FASTGA_MODE_REFERENCE].click()
    qapp.processEvents()
    qapp.processEvents()
    h_ref = dlg.height()
    assert h_ref > h_self
    assert panel._reference_field.isVisible()
    assert not panel.advanced_visible()
    # No leftover Advanced-sized blank under Settings.
    gap = dlg.height() - dlg._dialog_chrome_height() - panel._basic_panel.height()
    assert gap < 40
    dlg.close()