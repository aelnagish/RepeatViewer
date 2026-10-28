"""Tools > Primers windows: rules (button order, Options tab, Files-panel fill),
settings round trips, saved options, and the job worker."""

from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("primer3")

from PySide6.QtCore import QCoreApplication, QSettings, Qt  # noqa: E402
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
    from PySide6.QtWidgets import QMenu

    from repeatviewer.browser_file_groups import file_family_key
    from repeatviewer.primer_actions import PRIMER_MANUAL, PRIMER_MENU_ITEMS, PrimerActionsMixin

    class _Host(PrimerActionsMixin):
        pass

    menu = QMenu()
    _Host()._populate_primers_menu(menu)
    texts = [a.text() for a in menu.actions() if not a.isSeparator()]
    assert texts == list(PRIMER_MENU_ITEMS) == ["Design", "Map", "Database"]
    assert all("." not in t and "…" not in t for t in texts)
    assert PRIMER_MANUAL.is_file()
    html = PRIMER_MANUAL.read_text(encoding="utf-8")
    assert "—" not in html
    for label in dz.MODE_LABELS.values():
        assert label in html, label
    assert file_family_key("gene.fa_primers.html") == "gene.fa"
    assert file_family_key("gene.fa_primer_map_2.html") == "gene.fa"


def test_mode_list_is_grouped_and_headers_not_selectable(qapp, isolated_settings, data):
    from repeatviewer.views.primer_dialogs import PrimerDesignDialog

    grouped = [m for _f, groups in dz.MODE_FAMILIES for _g, modes in groups for m in modes]
    assert sorted(grouped) == sorted(dz.ALL_MODES) and len(grouped) == len(set(grouped))
    d = PrimerDesignDialog(initial=PrimerOptions(), files=[data.gene])
    model = d._mode.model()
    rows = [(model.item(i).text(), model.item(i).data(Qt.ItemDataRole.UserRole),
             bool(model.item(i).flags() & Qt.ItemFlag.ItemIsSelectable))
            for i in range(model.rowCount())]
    assert [r[0] for r in rows[:3]] == ["PCR primers", "General", "Standard PCR"]
    for text, mode, selectable in rows:
        assert selectable == bool(mode), text
    assert [m for _t, m, _s in rows if m] == grouped
    assert "FISH probes" in [r[0] for r in rows]
    assert d._mode.currentData() == dz.MODE_PCR


def test_inputs_define_template(qapp, isolated_settings, data):
    from repeatviewer.views.primer_dialogs import MapPrimersDialog, PrimerDesignDialog

    d = PrimerDesignDialog(initial=PrimerOptions(), files=[data.gene])
    d.show()
    qapp.processEvents()
    form = d._primers_form.form
    assert d._input_selected.isChecked()
    assert form.isRowVisible(d._template_label_field) and not form.isRowVisible(d._template_field)
    d._input_define.setChecked(True)
    qapp.processEvents()
    assert form.isRowVisible(d._template_field) and not form.isRowVisible(d._template_label_field)
    assert d.template_files() == []
    assert d.apply_path_from_panel(data.locus)
    assert [p.name for p in d.template_files()] == ["locus.fa"]
    d._input_selected.setChecked(True)
    assert [p.name for p in d.template_files()] == ["gene.fa"]
    d.close()
    # No selection: the window opens on Define inputs.
    e = PrimerDesignDialog(initial=PrimerOptions(), files=[])
    assert e._input_define.isChecked()
    m = MapPrimersDialog(initial=PrimerOptions(), files=[])
    assert m._input_define.isChecked() and m.run_files() == []
    assert m.apply_path_from_panel(data.array)
    assert [p.name for p in m.run_files()] == ["array.fa"]


def test_methylation_fields(qapp, isolated_settings, data):
    from repeatviewer.views.primer_dialogs import PrimerDesignDialog

    d = PrimerDesignDialog(initial=PrimerOptions(), files=[data.promoter])
    d.show()
    form = d._primers_form.form
    d._mode.setCurrentIndex(d._mode.findData(dz.MODE_MSP))
    qapp.processEvents()
    assert form.isRowVisible(d._meth_system.parentWidget()) or d._meth_system.isVisible()
    assert d.design_settings().methylation_system == dz.METH_PLANT
    _set = d._msp_context.findData("CHG")
    d._msp_context.setCurrentIndex(_set)
    assert d.design_settings().msp_context == "CHG"
    # Animal: CHG is not offered, the site falls back to CG.
    d._meth_system.setCurrentIndex(d._meth_system.findData(dz.METH_ANIMAL))
    s = d.design_settings()
    assert s.methylation_system == dz.METH_ANIMAL and s.msp_context == "CG"
    d._mode.setCurrentIndex(d._mode.findData(dz.MODE_PCR))
    qapp.processEvents()
    assert not d._meth_system.isVisible()
    d.close()


def test_painting_mode_swaps_tabs_and_round_trips(qapp, isolated_settings, data, tmp_path):
    from repeatviewer.backends.primers import painting as pt
    from repeatviewer.views.primer_dialogs import PrimerDesignDialog

    d = PrimerDesignDialog(initial=PrimerOptions(), files=[data.gene])
    d.show()
    qapp.processEvents()
    tabs = lambda: [d._tab_bar.tabText(i) for i in range(d._tab_bar.count())  # noqa: E731
                    if d._tab_bar.isTabVisible(i)]
    assert tabs() == ["Primers", "Regions", "Specificity", "Options"]
    d._mode.setCurrentIndex(d._mode.findData(dz.MODE_PAINTING))
    qapp.processEvents()
    assert tabs() == ["Chorus2", "Selection", "Options"]
    assert d._tab_bar.currentIndex() == d._TAB_CHORUS2 and d.is_painting()
    assert d._paint_mode.currentData() == dz.MODE_PAINTING
    # Import: a BED click fills the BED field, a FASTA click the genome.
    d._pt_import.setChecked(True)
    bed = tmp_path / "chorus.bed"
    bed.write_text("chr1\t1\t45\t" + "A" * 45 + "\n")
    assert d.apply_path_from_panel(bed)
    assert d.apply_path_from_panel(data.genome)
    d._pt_target.setCurrentIndex(d._pt_target.findData(pt.TARGET_POOL))
    d._pt_pool.setValue(1234)
    d._pt_regions.setText("chr1, chr2:1-5000")
    d._pt_p5.setText("TAATACGACTCACTATAGGG")
    opts = d.options()
    ps = opts.painting_settings()
    assert (ps.source, ps.import_bed, ps.genome_fasta) == (pt.SOURCE_IMPORT, str(bed),
                                                           str(data.genome))
    assert (ps.target, ps.pool_size, ps.regions, ps.primer5) == (
        pt.TARGET_POOL, 1234, "chr1, chr2:1-5000", "TAATACGACTCACTATAGGG")
    assert opts.design_settings().mode == dz.MODE_PAINTING
    assert opts.last_dialog_tab == "chorus2"
    d.close()
    # Reopen: comes back on painting, Chorus2 tab, same values.
    e = PrimerDesignDialog(initial=opts, files=[data.gene])
    assert e.is_painting() and e._tab_bar.currentIndex() == e._TAB_CHORUS2
    assert e.painting_settings() == ps
    # Back to a primer mode from the Chorus2 tab's Mode box.
    e._paint_mode.setCurrentIndex(e._paint_mode.findData(dz.MODE_PCR))
    assert not e.is_painting() and e._tab_bar.currentIndex() == e._TAB_PRIMERS
    assert e._mode.currentData() == dz.MODE_PCR
    for label in e.findChildren(QLabel):
        assert "—" not in label.text()
