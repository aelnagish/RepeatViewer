"""Tools > Primers windows: rules (button order, Options tab, Files-panel fill),
settings round trips, saved options, and the job worker."""

from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("primer3")

from PySide6.QtCore import QCoreApplication, QSettings  # noqa: E402
from PySide6.QtWidgets import QCheckBox, QLabel, QPushButton  # noqa: E402

from primer_testdata import make_data  # noqa: E402

from repeatviewer.backends.primers import design as dz  # noqa: E402
from repeatviewer.backends.primers import jobs as pj  # noqa: E402
from repeatviewer.job_options import (  # noqa: E402
    PRIMER_DIALOG_TAB_REGIONS,
    PRIMER_RIGHT_CLICK_DIALOG,
    PrimerOptions,
    load_primer_options,
    save_primer_options,
)


@pytest.fixture()
def isolated_settings(qapp, tmp_path):
    org, app = QCoreApplication.organizationName(), QCoreApplication.applicationName()
    QCoreApplication.setOrganizationName("RVPrimerTests")
    QCoreApplication.setApplicationName(f"t{id(tmp_path)}")
    QSettings().clear()
    yield
    QSettings().clear()
    QCoreApplication.setOrganizationName(org)
    QCoreApplication.setApplicationName(app)


@pytest.fixture()
def data(tmp_path):
    return make_data(tmp_path / "d", genome_flank=5_000)


def _buttons_right(dialog) -> list[str]:
    bar = dialog._save_button.parentWidget()
    return [b.text() for b in bar.findChildren(QPushButton) if b.isVisibleTo(dialog)]


def test_design_dialog_follows_tool_rules(qapp, isolated_settings, data):
    from repeatviewer.views.primer_dialogs import PrimerDesignDialog

    d = PrimerDesignDialog(initial=PrimerOptions(), files=[data.gene])
    texts = _buttons_right(d)
    # Cancel, Save, Run in this order; never OK.
    assert texts[-3:] == ["Cancel", "Save", "Run"] and "OK" not in texts
    assert [d._tab_bar.tabText(i) for i in range(d._tab_bar.count())][-1] == "Options"
    assert [p.name for p in d.template_files()] == ["gene.fa"]
    # A Files-panel click changes the template.
    assert d.apply_path_from_panel(data.monomers)
    assert [p.name for p in d.template_files()] == ["monomers.fa"]
    # A non-FASTA click does nothing.
    assert not d.apply_path_from_panel(data.root / "missing.fa")


def test_design_dialog_mode_rows_and_presets(qapp, isolated_settings, data):
    from repeatviewer.views.primer_dialogs import PrimerDesignDialog

    d = PrimerDesignDialog(initial=PrimerOptions(), files=[data.gene])
    d.show()
    qapp.processEvents()
    form = d._primers_form.form
    for mode in dz.DESIGN_MODES:
        d._mode.setCurrentIndex(d._mode.findData(mode))
        qapp.processEvents()
        visible = {id(w) for w in d._mode_rows[mode]}
        for rows in d._mode_rows.values():
            for w in rows:
                assert form.isRowVisible(w) == (id(w) in visible)
        s = d.design_settings()
        preset = dz.mode_preset(mode)
        assert (s.size_min, s.tm_opt, s.product_min) == (
            preset["size_min"], preset["tm_opt"], preset.get("product_min", s.product_min))
    d.close()


def test_design_dialog_settings_round_trip(qapp, isolated_settings, data):
    from repeatviewer.views.primer_dialogs import PrimerDesignDialog

    start = dz.DesignSettings(
        mode=dz.MODE_SPECIFIC, target_pattern="TEfamA", min_coverage=0.8, max_off_target=0.05,
        forward_tail="gtaaaacgacggccagt", mask_tandem=False, mask_file=str(data.primers),
        check_genome=True, genome_fasta=str(data.genome), genome_max_products=3,
        target_start=10, target_end=50, num_return=7,
    )
    opts = PrimerOptions(design=start.to_dict(), right_click_action=PRIMER_RIGHT_CLICK_DIALOG,
                         last_dialog_tab=PRIMER_DIALOG_TAB_REGIONS)
    d = PrimerDesignDialog(initial=opts, files=[data.te_aln])
    assert d._tab_bar.currentIndex() == 1  # reopens on Regions
    out = d.options()
    s = out.design_settings()
    for key in ("mode", "target_pattern", "min_coverage", "max_off_target", "forward_tail",
                "mask_tandem", "mask_file", "check_genome", "genome_fasta", "genome_max_products",
                "target_start", "target_end", "num_return"):
        assert getattr(s, key) == pytest.approx(getattr(start, key)) if isinstance(
            getattr(start, key), float) else getattr(s, key) == getattr(start, key), key
    assert out.right_click_action == PRIMER_RIGHT_CLICK_DIALOG


def test_restore_defaults_uses_mode_preset(qapp, isolated_settings, data):
    from repeatviewer.views.primer_dialogs import PrimerDesignDialog

    opts = PrimerOptions(design=dz.DesignSettings(mode=dz.MODE_QPCR, tm_opt=70.0).to_dict())
    d = PrimerDesignDialog(initial=opts, files=[data.gene])
    assert d.design_settings().tm_opt == 70.0
    d._restore_defaults()
    assert d.design_settings().tm_opt == dz.mode_preset(dz.MODE_QPCR)["tm_opt"]


def test_map_dialog_round_trip_and_files(qapp, isolated_settings, data, tmp_path):
    from repeatviewer.views.primer_dialogs import MapPrimersDialog

    ms = pj.MapSettings(use_forward=True, forward="Gene_F", use_reverse=True,
                        reverse="ACGTACGTACGTACGT", max_mismatches=1, three_prime_exact=4,
                        test_on=pj.TEST_ON_CONSENSUS, consensus_threshold=0.75,
                        use_product_range=True, product_min=50, product_max=900)
    opts = PrimerOptions(map=ms.to_dict(), database_path=str(tmp_path / "db.fa"))
    m = MapPrimersDialog(initial=opts, files=[data.gene])
    back = m.map_settings()
    for key in ("forward", "reverse", "max_mismatches", "three_prime_exact", "test_on",
                "consensus_threshold", "product_min", "product_max"):
        assert getattr(back, key) == getattr(ms, key), key
    assert _buttons_right(m)[-3:] == ["Cancel", "Save", "Run"]
    assert m.apply_path_from_panel(data.array)
    assert [p.name for p in m.run_files()] == ["array.fa"]


def test_checkbox_text_is_right_of_the_box(qapp, isolated_settings, data):
    """Rule 16: checkbox text on the right (default Qt layout, never right-to-left)."""
    from repeatviewer.views.primer_dialogs import MapPrimersDialog, PrimerDesignDialog

    for dlg in (PrimerDesignDialog(initial=PrimerOptions(), files=[data.gene]),
                MapPrimersDialog(initial=PrimerOptions(), files=[data.gene])):
        for box in dlg.findChildren(QCheckBox):
            assert not box.layoutDirection() == 1, box.text()  # Qt.RightToLeft
        for label in dlg.findChildren(QLabel):
            assert "—" not in label.text()


def test_database_dialog_lists_and_analyses(qapp, isolated_settings, data):
    from repeatviewer.views.primer_dialogs import PrimerDatabaseDialog

    db_path = data.primers
    d = PrimerDatabaseDialog(db_path)
    assert d._table.rowCount() == 4
    d._table.selectRow(0)
    qapp.processEvents()
    text = d._details.toPlainText()
    assert "Tm" in text and "hairpin" in text
    d._filter.setText("Sat")
    assert d._table.rowCount() == 2


def test_primer_options_persist(qapp, isolated_settings):
    opts = PrimerOptions(design=dz.DesignSettings(mode=dz.MODE_FISH, fish_label="Cy5").to_dict(),
                         database_path="/x/db.fa", save_design_to_database=True,
                         right_click_action=PRIMER_RIGHT_CLICK_DIALOG)
    save_primer_options(opts)
    back = load_primer_options()
    assert back.design_settings().mode == dz.MODE_FISH
    assert back.design_settings().fish_label == "Cy5"
    assert back.database_path == "/x/db.fa" and back.save_design_to_database
    assert back.right_click_action == PRIMER_RIGHT_CLICK_DIALOG
    # Unknown / old keys never break loading.
    QSettings().setValue("primers/design_json", "{not json")
    assert load_primer_options().design_settings().mode == dz.MODE_PCR


def test_worker_runs_and_reports(qapp, isolated_settings, data):
    from repeatviewer.primer_actions import PRIMER_JOB_DESIGN, PRIMER_JOB_MAP, PrimerJobWorker

    got, failed, steps = [], [], []
    w = PrimerJobWorker(PRIMER_JOB_DESIGN, [data.gene], dz.DesignSettings(), report_path=None,
                        save_to_database=True, database_path=str(data.root / "db.fa"))
    w.finished.connect(got.append)
    w.failed.connect(failed.append)
    w.progress.connect(lambda v, t: steps.append(v))
    w.run()
    assert not failed and got[0]["pairs"] == 5 and got[0]["saved_to_database"] == 10
    assert Path(got[0]["report"]).name == "gene.fa_primers.html"
    assert max(steps) <= 3 and steps[0] == 0
    w2 = PrimerJobWorker(PRIMER_JOB_MAP, [data.gene], pj.MapSettings(forward=""), report_path=None)
    w2.failed.connect(failed.append)
    w2.run()
    assert failed and "choose a saved primer" in failed[0]


def test_worker_cancel(qapp, isolated_settings, data):
    from repeatviewer.backends.job_cancel import CANCELLED_MESSAGE, JobCancellation
    from repeatviewer.primer_actions import PRIMER_JOB_DESIGN, PrimerJobWorker

    cancel = JobCancellation()
    cancel.request()
    failed = []
    w = PrimerJobWorker(PRIMER_JOB_DESIGN, [data.gene], dz.DesignSettings(), report_path=None,
                        cancel=cancel)
    w.failed.connect(failed.append)
    w.run()
    assert failed == [CANCELLED_MESSAGE]


def test_menu_texts_and_file_family(qapp):
    from repeatviewer.browser_file_groups import file_family_key
    from repeatviewer.primer_actions import PRIMER_MODE_MENU

    # Rule 11: nouns start upper case in menu items.
    for label, _mode in PRIMER_MODE_MENU:
        first = label.split()[0]
        assert first[0].isupper() or first == "qPCR", label
    assert file_family_key("gene.fa_primers.html") == "gene.fa"
    assert file_family_key("gene.fa_primer_map_2.html") == "gene.fa"
