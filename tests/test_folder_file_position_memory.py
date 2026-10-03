"""Returning to a folder restores the last selected Files-list position."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import QApplication, QAbstractItemView

from repeatviewer.widgets.file_list_table import FileListTableWidget


class _FolderFileMemoryHarness:
    """Minimal stand-in for MainWindow folder file-position memory helpers."""

    def __init__(self) -> None:
        self._current_folder: Path | None = None
        self._file_list = FileListTableWidget()
        self._last_file_by_folder: dict[str, str] = {}
        self._last_scroll_by_folder: dict[str, int] = {}
        self.opened_preview: Path | None = None

    def _folder_memory_key(self, folder: Path | None) -> str | None:
        if folder is None:
            return None
        try:
            return str(folder.expanduser().resolve())
        except OSError:
            return None

    def _preview_path_from_item(self, item) -> Path | None:
        if item is None:
            return None
        from PySide6.QtCore import Qt

        path_raw = item.data(Qt.ItemDataRole.UserRole)
        if not path_raw:
            return None
        path = Path(str(path_raw))
        return path if path.is_file() else None

    def _remember_file_list_position(self, folder: Path | None = None) -> None:
        target = folder if folder is not None else self._current_folder
        key = self._folder_memory_key(target)
        if key is None:
            return
        current = self._file_list.currentItem()
        path = self._preview_path_from_item(current)
        if path is not None:
            self._last_file_by_folder[key] = str(path.expanduser().resolve())
        if self._file_list.rowCount() > 0:
            self._last_scroll_by_folder[key] = int(
                self._file_list.verticalScrollBar().value()
            )

    def _reveal_file_list_item(self, item, *, scroll_hint) -> None:
        item.setSelected(True)
        self._file_list.setCurrentItem(item)
        self._file_list.scrollToItem(item, scroll_hint)

    def _open_file_in_preview(self, path: Path) -> None:
        self.opened_preview = path

    def _restore_file_list_position(
        self, folder: Path, *, open_preview: bool = True
    ) -> bool:
        key = self._folder_memory_key(folder)
        if key is None:
            return False
        remembered = self._last_file_by_folder.get(key)
        match = None
        if remembered:
            row = self._file_list._path_row_index().get(remembered)
            if row is not None:
                match = self._file_list.item(row)
        if match is not None:
            self._reveal_file_list_item(
                match,
                scroll_hint=QAbstractItemView.ScrollHint.EnsureVisible,
            )
            if open_preview:
                path = self._preview_path_from_item(match)
                if path is not None:
                    self._open_file_in_preview(path)
        scroll = self._last_scroll_by_folder.get(key)
        if scroll is not None:
            vsb = self._file_list.verticalScrollBar()
            vsb.setValue(min(int(scroll), vsb.maximum()))
            if match is not None:
                self._file_list.scrollToItem(
                    match, QAbstractItemView.ScrollHint.EnsureVisible
                )
        return match is not None


def test_returning_to_folder_restores_selected_file(tmp_path: Path) -> None:
    app = QApplication.instance() or QApplication([])
    del app

    folder_a = tmp_path / "A"
    folder_b = tmp_path / "B"
    folder_a.mkdir()
    folder_b.mkdir()
    files_a = []
    for name in ("one.fa", "two.fa", "three.fa"):
        path = folder_a / name
        path.write_text(f">{name}\nACGT\n")
        files_a.append(path)
    (folder_b / "other.fa").write_text(">o\nACGT\n")

    harness = _FolderFileMemoryHarness()
    harness._current_folder = folder_a
    for path in files_a:
        harness._file_list.add_file_row(path)

    # Select the middle file in folder A.
    mid = harness._file_list.item(1)
    harness._file_list.setCurrentItem(mid)
    harness._remember_file_list_position(folder_a)
    assert harness._last_file_by_folder[str(folder_a.resolve())].endswith("two.fa")

    # Leave to folder B (clear list like a folder switch).
    harness._file_list.clear_rows()
    harness._current_folder = folder_b
    harness._file_list.add_file_row(folder_b / "other.fa")
    harness._file_list.setCurrentItem(harness._file_list.item(0))

    # Return to folder A and rebuild the list.
    harness._file_list.clear_rows()
    harness._current_folder = folder_a
    for path in files_a:
        harness._file_list.add_file_row(path)

    assert harness._restore_file_list_position(folder_a, open_preview=True)
    current = harness._preview_path_from_item(harness._file_list.currentItem())
    assert current is not None
    assert current.name == "two.fa"
    assert harness.opened_preview is not None
    assert harness.opened_preview.name == "two.fa"
