"""Chromosome painting as a centre tab: probes, picture, settings panel."""

import json
import random

import pytest
from PySide6.QtCore import QCoreApplication, QSettings

from repeatviewer.backends.primers import painting as pt


@pytest.fixture()
def isolated_settings(qapp, tmp_path):
    org, app = QCoreApplication.organizationName(), QCoreApplication.applicationName()
    QCoreApplication.setOrganizationName("RVPaintingPage")
    QCoreApplication.setApplicationName(f"t{id(tmp_path)}")
    QSettings().clear()
    yield
    QSettings().clear()
    QCoreApplication.setOrganizationName(org)
    QCoreApplication.setApplicationName(app)


def _bed(tmp_path, counts=(("Chr1", 120_000_000, 60), ("Chr2", 95_000_000, 40))):
    r = random.Random(2)
    lines = []
    for chrom, length, n in counts:
        for _ in range(n):
            start = r.randrange(1, length - 50)
            lines.append(f"{chrom}\t{start}\t{start + 45}\t"
                         + "".join(r.choice("ACGT") for _ in range(45)))
    path = tmp_path / "probes.bed"
    path.write_text("\n".join(lines) + "\n")
    return path


def test_importing_a_bed_fills_the_table_and_the_picture(qapp, isolated_settings, tmp_path):
    from repeatviewer.views.painting_page import PaintingPage

    page = PaintingPage()
    assert page.tab_name() == "Painting Probes"
    assert page.chorus_section.isVisibleTo(page)          # running Chorus2 by default

    sources = [page.source.itemData(i) for i in range(page.source.count())]
    page.source.setCurrentIndex(sources.index(pt.SOURCE_IMPORT))
    assert not page.chorus_section.isVisibleTo(page)      # nothing to run here
    assert page.run_button.text() == "Select"

    page.set_probes(pt.parse_probe_bed(_bed(tmp_path)))
    assert page.table.rowCount() == 2
    assert page.table.item(0, 0).text() == "Chr1" and page.table.item(0, 1).text() == "60"
    assert "100 probes on 2 chromosomes" in page.verdict.text()
    assert [row[0] for row in page.picture._rows] == ["Chr1", "Chr2"]
    page.picture.resize(700, 100)
    page.picture.grab()
    page.deleteLater()


def test_amount_rows_follow_the_choice_and_settings_save(qapp, isolated_settings, tmp_path):
    from repeatviewer.views.painting_page import PaintingPage

    page = PaintingPage()
    targets = [page.target.itemData(i) for i in range(page.target.count())]
    assert not page.density_row.isVisibleTo(page) and not page.pool_row.isVisibleTo(page)
    page.target.setCurrentIndex(targets.index(pt.TARGET_DENSITY))
    assert page.density_row.isVisibleTo(page) and not page.pool_row.isVisibleTo(page)
    page.target.setCurrentIndex(targets.index(pt.TARGET_POOL))
    assert page.pool_row.isVisibleTo(page) and not page.density_row.isVisibleTo(page)

    page.length.setValue(50)
    page.homology.setValue(80.0)
    page.spacing.setValue(100)
    page.pool.setValue(12000)
    page.pool_name.setText("CroSat pool")
    page.primer5.setText("ACGTACGT")
    s = page.settings()
    assert (s.length, s.homology, s.spacing, s.pool_size) == (50, 80.0, 100, 12000)
    assert s.pool_name == "CroSat pool" and s.primer5 == "ACGTACGT"
    assert s.target == pt.TARGET_POOL

    page.save_button.click()
    saved = json.loads(QSettings().value("primers/painting_page"))
    assert saved["length"] == 50 and saved["pool_name"] == "CroSat pool"
    page.defaults_button.click()
    assert page.length.value() == pt.PaintingSettings().length
    page.deleteLater()


def test_run_hands_the_settings_to_the_job(qapp, isolated_settings, tmp_path):
    from repeatviewer.views.painting_page import PaintingPage

    started: list[dict] = []
    page = PaintingPage(start_job=started.append)
    page.length.setValue(45)
    page.threads.setValue(8)
    page.run_button.click()
    assert len(started) == 1
    assert started[0]["threads"] == 8 and started[0]["source"] == pt.SOURCE_RUN
    page.deleteLater()
