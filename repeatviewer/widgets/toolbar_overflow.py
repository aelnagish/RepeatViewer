"""Hide whole toolbar controls from right to left when width is too small.

Never leaves a partially clipped button/label visible.
"""

from __future__ import annotations

from PySide6.QtCore import QEvent, QObject
from PySide6.QtWidgets import QLayout, QWidget


class ToolbarOverflowHider(QObject):
    """Show as many left-side units as fit; hide the rest as whole units.

    ``units`` is left-to-right: the rightmost unit is hidden first when the
    host shrinks. Units that are ``set_unit_allowed(False)`` stay hidden
    (e.g. Clean/Cons outside MSA mode) and are skipped by overflow logic.
    """

    # Leave a few pixels so the last visible control is not flush-clipped.
    _SAFETY_PAD = 8

    def __init__(
        self,
        host: QWidget,
        units: list[list[QWidget]],
        *,
        spacing: int = 6,
        layout: QLayout | None = None,
    ) -> None:
        super().__init__(host)
        self._host = host
        self._units = [list(unit) for unit in units if unit]
        self._spacing = max(0, int(spacing))
        self._layout = layout if layout is not None else host.layout()
        self._allowed = [True] * len(self._units)
        self._width_sources: list[QWidget] = []
        self._syncing = False
        host.installEventFilter(self)
        if self._layout is not None:
            self._layout.setSizeConstraint(QLayout.SizeConstraint.SetNoConstraint)
        host.setMinimumWidth(0)

    def add_width_source(self, widget: QWidget | None) -> None:
        """Also watch ``widget`` for resize (e.g. map toolbar slot / pane)."""
        if widget is None or widget is self._host:
            return
        if widget in self._width_sources:
            return
        self._width_sources.append(widget)
        widget.installEventFilter(self)

    def set_unit_allowed(self, index: int, allowed: bool) -> None:
        if not (0 <= index < len(self._units)):
            return
        self._allowed[index] = bool(allowed)
        self.sync()

    def set_unit_allowed_by_widget(self, widget: QWidget, allowed: bool) -> None:
        for index, unit in enumerate(self._units):
            if widget in unit:
                self.set_unit_allowed(index, allowed)
                return

    def eventFilter(self, watched, event) -> bool:  # noqa: N802
        if watched is self._host or watched in self._width_sources:
            if event.type() in (
                QEvent.Type.Resize,
                QEvent.Type.Show,
                QEvent.Type.LayoutRequest,
            ):
                self.sync()
        return False

    def sync(self) -> None:
        if self._syncing or self._host is None:
            return
        self._syncing = True
        try:
            self._apply()
        finally:
            self._syncing = False

    def _widget_width(self, widget: QWidget) -> int:
        # Prefer an explicit fixed width (stable after stylesheet polish).
        max_w = int(widget.maximumWidth())
        min_w = int(widget.minimumWidth())
        if 0 < min_w == max_w < QWIDGETSIZE_MAX:
            return min_w
        hint = widget.sizeHint()
        width = max(0, int(hint.width()))
        if width <= 0:
            width = max(0, int(widget.minimumSizeHint().width()))
        if width <= 0:
            width = max(0, int(widget.width()))
        if min_w > 0:
            width = max(width, min_w)
        return width

    def _unit_width(self, unit: list[QWidget]) -> int:
        widths = [self._widget_width(w) for w in unit]
        if not widths:
            return 0
        return sum(widths) + self._spacing * (len(widths) - 1)

    def _margins_horizontal(self) -> tuple[int, int]:
        if self._layout is None:
            return 0, 0
        left, _top, right, _bottom = self._layout.getContentsMargins()
        return int(left), int(right)

    def _available_width(self) -> int:
        """Width the host can actually use (never larger than a narrow parent)."""
        left_m, right_m = self._margins_horizontal()
        margins = left_m + right_m
        candidates: list[int] = []
        host_w = max(0, int(self._host.width()) - margins)
        if host_w > 0:
            candidates.append(host_w)
        for src in self._width_sources:
            sw = max(0, int(src.width()) - margins)
            if sw > 0:
                candidates.append(sw)
        # Walk parents: in split panes the host can still report a wide
        # sizeHint while a narrower ancestor clips it.
        parent = self._host.parentWidget()
        depth = 0
        while parent is not None and depth < 6:
            pw = max(0, int(parent.width()) - margins)
            if pw > 0:
                candidates.append(pw)
            parent = parent.parentWidget()
            depth += 1
        if not candidates:
            return 0
        # The tightest constraint wins — never pretend we have the full
        # preferred toolbar width when the map column is narrower.
        return max(0, min(candidates) - self._SAFETY_PAD)

    def _apply(self) -> None:
        available = self._available_width()
        allowed_indexes = [i for i, ok in enumerate(self._allowed) if ok]
        # Before the host has a real width (not shown / not laid out yet),
        # keep every allowed unit visible — do not hide everything at width 0.
        if available <= 0 or (
            not self._host.isVisible() and available < 80
        ):
            hide_from = len(allowed_indexes)
        else:
            hide_from = len(allowed_indexes)
            while hide_from > 0:
                visible = allowed_indexes[:hide_from]
                needed = 0
                for pos, index in enumerate(visible):
                    needed += self._unit_width(self._units[index])
                    if pos:
                        needed += self._spacing
                if needed <= available:
                    break
                hide_from -= 1

        visible_set = set(allowed_indexes[:hide_from])
        for index, unit in enumerate(self._units):
            show = index in visible_set
            for widget in unit:
                if widget.isHidden() == show:
                    widget.setVisible(show)


# Match Qt's QWIDGETSIZE_MAX without importing private headers.
QWIDGETSIZE_MAX = 16777215
