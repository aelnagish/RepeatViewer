"""Shared overwrite / cancel / rename prompt when a target path already exists."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import (
    QDialog,
    QInputDialog,
    QMessageBox,
    QWidget,
)


def prompt_rename_text(
    parent: QWidget | None,
    *,
    title: str,
    label: str,
    initial: str,
) -> tuple[str, bool]:
    """Wide rename prompt so long filenames fit without horizontal scrolling."""
    dialog = QInputDialog(parent)
    dialog.setWindowTitle(title)
    dialog.setLabelText(label)
    dialog.setTextValue(initial)
    dialog.setInputMode(QInputDialog.InputMode.TextInput)
    dialog.resize(560, dialog.sizeHint().height())
    dialog.setMinimumWidth(520)
    ok = dialog.exec() == QDialog.DialogCode.Accepted
    return dialog.textValue(), ok


def confirm_existing_file_path(
    parent: QWidget | None,
    proposed: Path,
    *,
    title: str,
    informative: str | None = None,
    force_html_suffix: bool = False,
) -> Path | None:
    """If ``proposed`` already exists, ask to overwrite, rename, or cancel.

    Returns the path to use, or ``None`` if the user cancels.
    When ``force_html_suffix`` is True (MSA / genome plots), renamed paths
    without ``.html`` / ``.htm`` get ``.html`` appended.
    """
    proposed = proposed.expanduser().resolve()
    if not proposed.exists():
        return proposed

    box = QMessageBox(parent)
    box.setWindowTitle(title)
    box.setIcon(QMessageBox.Icon.Warning)
    box.setText(f"“{proposed.name}” already exists in this folder.")
    box.setInformativeText(
        informative
        or "Overwrite the existing file, or save under a different name?"
    )
    overwrite_btn = box.addButton("Overwrite", QMessageBox.ButtonRole.DestructiveRole)
    rename_btn = box.addButton("Rename…", QMessageBox.ButtonRole.AcceptRole)
    box.addButton(QMessageBox.StandardButton.Cancel)
    box.setDefaultButton(rename_btn)
    box.exec()
    clicked = box.clickedButton()
    if clicked == overwrite_btn:
        return proposed
    if clicked != rename_btn:
        return None

    new_name, ok = prompt_rename_text(
        parent,
        title=f"{title} — rename",
        label="New file name:",
        initial=proposed.name,
    )
    if not ok or not new_name.strip():
        return None
    name = new_name.strip()
    if Path(name).name != name or "/" in name or "\\" in name:
        QMessageBox.warning(
            parent,
            title,
            "Enter a file name only (no folder path).",
        )
        return confirm_existing_file_path(
            parent,
            proposed,
            title=title,
            informative=informative,
            force_html_suffix=force_html_suffix,
        )
    path = (proposed.parent / name).expanduser().resolve()
    if force_html_suffix and path.suffix.lower() not in {".html", ".htm"}:
        path = path.with_suffix(".html")
    if path == proposed:
        return confirm_existing_file_path(
            parent,
            proposed,
            title=title,
            informative=informative,
            force_html_suffix=force_html_suffix,
        )
    if path.exists():
        return confirm_existing_file_path(
            parent,
            path,
            title=title,
            informative=informative,
            force_html_suffix=force_html_suffix,
        )
    return path
