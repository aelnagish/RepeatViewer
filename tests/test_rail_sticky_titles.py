"""Page titles stay pinned above the scroll (Annotation Annotations pattern)."""

from __future__ import annotations

from PySide6.QtWidgets import QLabel, QScrollArea, QWidget


def _assert_title_outside_scroll(shell: QWidget, title: str) -> None:
    labels = [lab for lab in shell.findChildren(QLabel) if lab.text() == title]
    assert labels, f"missing page title {title!r}"
    label = labels[0]
    ancestor: QWidget | None = label
    while ancestor is not None and ancestor is not shell:
        assert not isinstance(ancestor, QScrollArea), (
            f"{title!r} must stay fixed above the scroll area"
        )
        ancestor = ancestor.parentWidget()
    assert shell.findChildren(QScrollArea), f"{title!r} page should still scroll"


def test_annotation_rail_sticky_page_titles(qapp) -> None:
    from repeatviewer.widgets.sequence_annotation_rail import SequenceAnnotationRail

    rail = SequenceAnnotationRail()
    cases = (
        (rail._show_page, "Display"),
        (rail._annotations_page, "Annotations"),
        (rail._find_page, "Find"),
        (rail._others_page, "Others"),
        (rail._primers_page, "Primers"),
    )
    for page, title in cases:
        _assert_title_outside_scroll(page, title)
    assert getattr(rail, "_tracks_divider", None) is not None
    rail.deleteLater()


def test_dotplot_rail_sticky_display_title(qapp) -> None:
    from repeatviewer.widgets.dotplot_rail import DotplotRail

    rail = DotplotRail()
    _assert_title_outside_scroll(rail._display_page, "Display")
    # Body controls still live under the scroll, not under the title host.
    scroll = rail._display_page.findChild(QScrollArea)
    assert scroll is not None and scroll.widget() is not None
    assert rail.line_width.isAncestorOf(scroll.widget()) or scroll.widget().isAncestorOf(
        rail.line_width
    )
    rail.deleteLater()


def test_tree_rail_sticky_display_and_find_titles(qapp) -> None:
    from repeatviewer.widgets.tree_display_rail import TreeDisplayRail

    content = QWidget()
    find = QWidget()
    rail = TreeDisplayRail(content, find_page=find)
    display_page = rail._stack.widget(0)
    find_page = rail._stack.widget(1)
    assert display_page is not None and find_page is not None
    _assert_title_outside_scroll(display_page, "Display")
    _assert_title_outside_scroll(find_page, "Find")
    rail.deleteLater()
