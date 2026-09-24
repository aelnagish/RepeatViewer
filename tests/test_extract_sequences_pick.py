"""Extract Sequences dialog — pick-list multi-select and search."""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QApplication

from repeatviewer.views.tool_dialogs import ExtractSequencesDialog


def _dialog(qapp: QApplication) -> ExtractSequencesDialog:
    del qapp
    records = [
        ("seq_a", "alpha", 10),
        ("seq_b", "beta", 20),
        ("seq_c", "gamma", 30),
        ("other", "delta", 40),
    ]
    dialog = ExtractSequencesDialog(
        records=records,
        source_name="demo.fasta",
        parent=None,
    )
    return dialog


def test_pick_filter_hides_non_matching(qapp: QApplication) -> None:
    dialog = _dialog(qapp)
    dialog._apply_pick_filter("seq_")
    visible = [item.text() for item in dialog._visible_items()]
    assert len(visible) == 3
    assert all("seq_" in text for text in visible)
    dialog._apply_pick_filter("")
    assert len(dialog._visible_items()) == 4
    dialog.close()


def test_select_all_visible_respects_filter(qapp: QApplication) -> None:
    dialog = _dialog(qapp)
    dialog._apply_pick_filter("seq_a")
    dialog._select_all_visible()
    assert dialog.selected_ids() == {"seq_a"}
    dialog._apply_pick_filter("")
    # Hidden-then-shown again: still only seq_a checked.
    assert dialog.selected_ids() == {"seq_a"}
    dialog.close()


def test_extended_selection_syncs_checkboxes(qapp: QApplication) -> None:
    dialog = _dialog(qapp)
    first = dialog._list.item(0)
    third = dialog._list.item(2)
    assert first is not None and third is not None
    dialog._pick_syncing = True
    try:
        dialog._list.clearSelection()
        first.setSelected(True)
        third.setSelected(True)
    finally:
        dialog._pick_syncing = False
    dialog._on_pick_selection_changed()
    assert first.checkState() == Qt.CheckState.Checked
    assert third.checkState() == Qt.CheckState.Checked
    assert dialog.selected_ids() == {"seq_a", "seq_c"}
    dialog.close()


def test_checkbox_clicks_accumulate_without_clearing(qapp: QApplication) -> None:
    """Ticking boxes one-by-one must not wipe earlier ticks (regression)."""
    from PySide6.QtCore import QEvent, QPoint
    from PySide6.QtGui import QMouseEvent

    dialog = _dialog(qapp)
    list_widget = dialog._list
    ids: list[str] = []
    for row in range(4):
        item = list_widget.item(row)
        assert item is not None
        rect = list_widget.visualItemRect(item)
        # Click near the left edge where the checkbox lives.
        pos = QPoint(rect.left() + 8, rect.center().y())
        press = QMouseEvent(
            QEvent.Type.MouseButtonPress,
            pos.toPointF(),
            Qt.MouseButton.LeftButton,
            Qt.MouseButton.LeftButton,
            Qt.KeyboardModifier.NoModifier,
        )
        list_widget.mousePressEvent(press)
        list_widget.mouseReleaseEvent(
            QMouseEvent(
                QEvent.Type.MouseButtonRelease,
                pos.toPointF(),
                Qt.MouseButton.LeftButton,
                Qt.MouseButton.NoButton,
                Qt.KeyboardModifier.NoModifier,
            )
        )
        ids.append(str(item.data(Qt.ItemDataRole.UserRole)))
        assert item.checkState() == Qt.CheckState.Checked
        assert set(ids) <= dialog.selected_ids()
    assert dialog.selected_ids() == {"seq_a", "seq_b", "seq_c", "other"}
    dialog.close()


def test_checkbox_click_does_not_fill_gap(qapp: QApplication) -> None:
    """Tick rows 1–3 then 5: the skipped row 4 must stay unticked."""
    records = [
        ("s1", "", 1),
        ("s2", "", 2),
        ("s3", "", 3),
        ("s4", "", 4),
        ("s5", "", 5),
    ]
    dialog = ExtractSequencesDialog(
        records=records, source_name="gap.fasta", parent=None
    )
    for row in (0, 1, 2, 4):
        item = dialog._list.item(row)
        assert item is not None
        dialog._toggle_pick_checkbox(item)
    assert dialog.selected_ids() == {"s1", "s2", "s3", "s5"}
    gap = dialog._list.item(3)
    assert gap is not None
    assert gap.checkState() == Qt.CheckState.Unchecked
    assert not gap.isSelected()
    dialog.close()


def test_select_none_clears(qapp: QApplication) -> None:
    dialog = _dialog(qapp)
    dialog._select_all_visible()
    assert len(dialog.selected_ids()) == 4
    dialog._select_none()
    assert dialog.selected_ids() == set()
    dialog.close()
