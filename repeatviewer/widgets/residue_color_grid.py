"""Residue → color button grid (shared by MSA plot and Annotation colors)."""

from __future__ import annotations

from PySide6.QtCore import QSize, Qt
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (
    QColorDialog,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QVBoxLayout,
    QWidget,
)


class ResidueColorGrid(QWidget):
    """Residue → color buttons, optionally laid out by category."""

    # Wide enough for "Negative charge" / "Positive charge"; keep buttons aligned.
    LABEL_WIDTH = 140

    def __init__(
        self,
        parent: QWidget,
        *,
        residues: list[str] | None = None,
        colors: dict[str, str],
        categories: list[tuple[str, list[str] | tuple[str, ...]]] | None = None,
        button_labels: dict[str, str] | None = None,
        linked_residues: dict[str, str] | None = None,
    ) -> None:
        super().__init__(parent)
        self._buttons: dict[str, QPushButton] = {}
        self._button_labels = dict(button_labels or {})
        # Bidirectional links (e.g. T ↔ U).
        self._linked: dict[str, str] = {}
        for left, right in (linked_residues or {}).items():
            self._linked[str(left).upper()] = str(right).upper()
            self._linked[str(right).upper()] = str(left).upper()
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(6)

        if categories:
            grid = QGridLayout()
            grid.setContentsMargins(0, 0, 0, 0)
            grid.setHorizontalSpacing(6)
            grid.setVerticalSpacing(6)
            grid.setColumnMinimumWidth(0, self.LABEL_WIDTH)
            grid.setColumnStretch(1, 1)
            for row_index, (category_name, group) in enumerate(categories):
                label = QLabel(category_name, self)
                label.setFixedWidth(self.LABEL_WIDTH)
                label.setAlignment(
                    Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
                )
                # Use the theme text color so labels stay readable in dark mode.
                label.setStyleSheet("color: palette(window-text);")
                buttons = QWidget(self)
                row = QHBoxLayout(buttons)
                row.setContentsMargins(0, 0, 0, 0)
                row.setSpacing(6)
                for residue in group:
                    row.addWidget(self._make_residue_button(residue, colors))
                row.addStretch(1)
                grid.addWidget(label, row_index, 0, Qt.AlignmentFlag.AlignRight)
                grid.addWidget(buttons, row_index, 1)
            outer.addLayout(grid)
        else:
            grid = QGridLayout()
            grid.setContentsMargins(0, 0, 0, 0)
            grid.setHorizontalSpacing(6)
            grid.setVerticalSpacing(4)
            items = list(residues or [])
            columns = 8 if len(items) > 8 else max(1, len(items))
            for index, residue in enumerate(items):
                grid.addWidget(
                    self._make_residue_button(residue, colors),
                    index // columns,
                    index % columns,
                )
            outer.addLayout(grid)

    def _make_residue_button(self, residue: str, colors: dict[str, str]) -> QPushButton:
        color = colors.get(residue, "#c7d1d0")
        caption = self._button_labels.get(residue, residue)
        button = QPushButton(caption, self)
        button.setFixedSize(QSize(34, 28))
        button.setCursor(Qt.CursorShape.PointingHandCursor)
        if not caption or caption == residue:
            tip_name = residue
        else:
            tip_name = f"{residue} ({caption})"
        button.setToolTip(f"Choose color for {tip_name}")
        self._apply_button_color(button, color)
        button.clicked.connect(lambda checked=False, r=residue: self._pick_color(r))
        self._buttons[residue] = button
        return button

    @staticmethod
    def contrast_text_color(background: QColor) -> str:
        """Black or white text for readable letters on light yellows etc."""

        def _channel(value: int) -> float:
            x = value / 255.0
            return x / 12.92 if x <= 0.03928 else ((x + 0.055) / 1.055) ** 2.4

        luminance = (
            0.2126 * _channel(background.red())
            + 0.7152 * _channel(background.green())
            + 0.0722 * _channel(background.blue())
        )
        return "#111111" if luminance > 0.45 else "#ffffff"

    @classmethod
    def _apply_button_color(cls, button: QPushButton, color: str) -> None:
        qcolor = QColor(color)
        if not qcolor.isValid():
            qcolor = QColor("#c7d1d0")
        text = cls.contrast_text_color(qcolor)
        button.setStyleSheet(
            f"QPushButton {{ background-color: {qcolor.name()}; color: {text}; "
            f"border: 1px solid #777; border-radius: 3px; font-weight: 600; }}"
        )
        button.setProperty("hex_color", qcolor.name())

    def _set_residue_color(self, residue: str, color: str) -> None:
        button = self._buttons.get(residue)
        if button is None:
            return
        self._apply_button_color(button, color)
        linked = self._linked.get(residue.upper())
        if linked and linked in self._buttons:
            self._apply_button_color(self._buttons[linked], color)

    def _pick_color(self, residue: str) -> None:
        button = self._buttons[residue]
        current = QColor(str(button.property("hex_color") or "#c7d1d0"))
        tip = self._button_labels.get(residue) or residue
        picked = QColorDialog.getColor(current, self.window(), f"Color for {tip}")
        if picked.isValid():
            self._set_residue_color(residue, picked.name())

    def colors(self) -> dict[str, str]:
        out = {
            residue: str(button.property("hex_color") or "#c7d1d0")
            for residue, button in self._buttons.items()
        }
        # Keep linked pairs identical in the saved map.
        for left, right in self._linked.items():
            if left in out and right in out:
                out[right] = out[left]
        return out

    def set_colors(self, colors: dict[str, str]) -> None:
        for residue, button in self._buttons.items():
            if residue in colors:
                self._apply_button_color(button, colors[residue])
        # Enforce links after bulk apply.
        for left, right in list(self._linked.items()):
            if left in self._buttons and right in self._buttons:
                color = str(self._buttons[left].property("hex_color") or "")
                if color:
                    self._apply_button_color(self._buttons[right], color)
