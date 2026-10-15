"""Hide whole toolbar controls from right to left when width is too small.

Never leaves a partially clipped button/label visible.
"""

from __future__ import annotations

from PySide6.QtCore import QEvent, QObject
from PySide6.QtWidgets import QLayout, QWidget


def _qt_alive(obj: object | None) -> bool:
    """False when a QObject's C++ side is already gone (shiboken)."""
    if obj is None:
        return False
    try:
        from shiboken6 import isValid

        if not bool(isValid(obj)):
            return False
    except Exception:
        pass
    # Even when isValid is unavailable / optimistic, probe a cheap call —
    # deleted Overview toolbar slots raise RuntimeError here.
    try:
        obj.objectName()  # type: ignore[attr-defined]
        return True
    except RuntimeError:
        return False


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
        if widget is None or widget is self._host or not _qt_alive(widget):
            return
        self._prune_width_sources()
        if widget in self._width_sources:
            return
        self._width_sources.append(widget)
        widget.installEventFilter(self)

    def remove_width_source(self, widget: QWidget | None) -> None:
        """Stop watching ``widget`` (call before destroying an Overview panel)."""
        if widget is None:
            return
        if widget in self._width_sources:
            self._width_sources.remove(widget)
        if _qt_alive(widget):
            try:
                widget.removeEventFilter(self)
            except RuntimeError:
                pass

    def clear_width_sources(self) -> None:
        """Drop every external width source (Overview unbind / tab replace)."""
        for src in list(self._width_sources):
            self.remove_width_source(src)
        self._width_sources.clear()

    def _prune_width_sources(self) -> None:
        alive: list[QWidget] = []
        for src in self._width_sources:
            try:
                if _qt_alive(src):
                    alive.append(src)
            except RuntimeError:
                continue
        self._width_sources = alive

    def set_unit_allowed(self, index: int, allowed: bool) -> None:
        if not (0 <= index < len(self._units)):
            return
        self._allowed[index] = bool(allowed)
        self._prune_width_sources()
        self.sync()

    def set_unit_allowed_by_widget(self, widget: QWidget, allowed: bool) -> None:
        # Drop destroyed Overview map hosts before measuring widths.
        self._prune_width_sources()
        for index, unit in enumerate(self._units):
            if widget in unit:
                self.set_unit_allowed(index, allowed)
                return

    def eventFilter(self, watched, event) -> bool:  # noqa: N802
        if not _qt_alive(self) or not _qt_alive(self._host):
            return False
        self._prune_width_sources()
        if watched is self._host or watched in self._width_sources:
            if event.type() in (
                QEvent.Type.Resize,
                QEvent.Type.Show,
                QEvent.Type.LayoutRequest,
            ):
                self.sync()
        return False

    def sync(self) -> None:
        if self._syncing or not _qt_alive(self) or not _qt_alive(self._host):
            return
        self._syncing = True
        try:
            self._apply()
        except RuntimeError:
            # Host / width source deleted mid-sync (tab replace).
            self._prune_width_sources()
        finally:
            self._syncing = False

    def _widget_width(self, widget: QWidget) -> int:
        if not _qt_alive(widget):
            return 0
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
        widths = [self._widget_width(w) for w in unit if _qt_alive(w)]
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
        if not _qt_alive(self._host):
            return 0
        left_m, right_m = self._margins_horizontal()
        margins = left_m + right_m
        candidates: list[int] = []
        host_w = max(0, int(self._host.width()) - margins)
        if host_w > 0:
            candidates.append(host_w)
        self._prune_width_sources()
        for src in self._width_sources:
            try:
                sw = max(0, int(src.width()) - margins)
            except RuntimeError:
                continue
            if sw > 0:
                candidates.append(sw)
        # Walk parents: in split panes the host can still report a wide
        # sizeHint while a narrower ancestor clips it.
        parent = self._host.parentWidget()
        depth = 0
        while parent is not None and depth < 6:
            if not _qt_alive(parent):
                break
            try:
                pw = max(0, int(parent.width()) - margins)
            except RuntimeError:
                break
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
        if not _qt_alive(self._host):
            return
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
                if not _qt_alive(widget):
                    continue
                if widget.isHidden() == show:
                    widget.setVisible(show)


# Match Qt's QWIDGETSIZE_MAX without importing private headers.
QWIDGETSIZE_MAX = 16777215
