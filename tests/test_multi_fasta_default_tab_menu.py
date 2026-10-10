"""Default preview tab Actions menu for multi-sequence FASTA."""

from __future__ import annotations

from PySide6.QtWidgets import QApplication, QMenu

from repeatviewer.widgets.preview_controls import (
    AnnotationZoomBar,
    PreviewZoomBar,
    _add_default_fasta_tab_menu,
)
from repeatviewer.widgets.preview_zoom import (
    DefaultFastaPreviewTab,
    save_default_fasta_preview_tab,
)


def test_multi_fasta_default_tab_menu_text_and_overview(qapp: QApplication) -> None:
    menu = QMenu()
    action = _add_default_fasta_tab_menu(menu)
    action.set_multi_fasta_choices(True)
    submenu = action.menu()
    assert submenu is not None
    labels = [a.text() for a in submenu.actions() if a.isCheckable()]
    assert labels == ["Text", "Overview"]
    assert "Dotplot (Self)" not in labels
    assert "Annotation" not in labels


def test_single_fasta_default_tab_menu_keeps_overview_dotplot(
    qapp: QApplication,
) -> None:
    menu = QMenu()
    action = _add_default_fasta_tab_menu(menu)
    action.set_multi_fasta_choices(False)
    submenu = action.menu()
    assert submenu is not None
    labels = [a.text() for a in submenu.actions() if a.isCheckable()]
    assert labels == ["Text", "Overview", "Dotplot (Self)"]


def test_annotation_zoom_bar_hides_actions_in_overview_mode(
    qapp: QApplication,
) -> None:
    """Overview mounts −/+/Select/Fit on the map; Actions stays on Preview title."""
    bar = AnnotationZoomBar()
    bar.set_overview_mode(True)
    assert not bar._overview_host.isHidden()
    assert bar._actions_btn.isHidden()
    bar.deleteLater()


def test_preview_zoom_bar_multi_fasta_default_tabs(qapp: QApplication) -> None:
    save_default_fasta_preview_tab(DefaultFastaPreviewTab.ANNOTATION)
    bar = PreviewZoomBar()
    bar.set_multi_fasta_default_tabs(True)
    bar.set_default_preview_tab_visible(True)
    submenu = bar._default_tab_action.menu()
    assert submenu is not None
    labels = [a.text() for a in submenu.actions() if a.isCheckable()]
    assert labels == ["Text", "Overview"]
    bar.deleteLater()


def test_overview_actions_chrome_is_default_tab_only(qapp: QApplication) -> None:
    """Overview title-bar Actions matches Dotplot: no font/%, only Default preview tab."""
    bar = PreviewZoomBar()
    bar.set_chrome_mode("dotplot")
    bar.set_default_preview_tab_visible(True)
    bar.set_multi_fasta_default_tabs(True)
    assert not bar._actions_btn.isHidden()
    assert bar._spin_cluster.isHidden()
    assert bar._font_cluster.isHidden()
    assert not bar._reset_scale_action.isVisible()
    assert not bar._edit_text_action.isVisible()
    assert not bar._save_text_action.isVisible()
    assert not bar._default_tab_sep.isVisible()
    assert bar._default_tab_action.isVisible()
    bar.deleteLater()

