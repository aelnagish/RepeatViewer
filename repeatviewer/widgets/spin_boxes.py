"""Spin boxes that ignore mouse-wheel changes unless the user uses the arrows."""

from __future__ import annotations

from PySide6.QtCore import QPoint, Qt
from PySide6.QtGui import QMouseEvent, QWheelEvent
from PySide6.QtWidgets import (
    QDoubleSpinBox,
    QSlider,
    QSpinBox,
    QStyle,
    QStyleOptionSlider,
    QStyleOptionSpinBox,
)


def _hit_spin_stepper(widget: QSpinBox | QDoubleSpinBox, pos: QPoint) -> bool:
    """True when ``pos`` is on the up/down stepper buttons."""
    option = QStyleOptionSpinBox()
    widget.initStyleOption(option)
    style = widget.style()
    up = style.subControlRect(
        QStyle.ComplexControl.CC_SpinBox,
        option,
        QStyle.SubControl.SC_SpinBoxUp,
        widget,
    )
    down = style.subControlRect(
        QStyle.ComplexControl.CC_SpinBox,
        option,
        QStyle.SubControl.SC_SpinBoxDown,
        widget,
    )
    return up.contains(pos) or down.contains(pos)


class NoWheelSpinBox(QSpinBox):
    """Do not change value when the scroll wheel passes over the widget.

    Wheel events are ignored so a parent (e.g. Preview zoom bar) can still
    handle Shift+scroll zoom after the spin was clicked.

    Keyboard ↑/↓ keep auto-repeating while held. Focus is cleared only after
    clicking the stepper buttons so Shift+scroll is not stuck on the spin.
    """

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self._clear_focus_on_stepper_release = False

    def wheelEvent(self, event: QWheelEvent) -> None:  # noqa: N802 (Qt API)
        event.ignore()

    def mousePressEvent(self, event: QMouseEvent) -> None:  # noqa: N802 (Qt API)
        self._clear_focus_on_stepper_release = (
            event.button() == Qt.MouseButton.LeftButton
            and _hit_spin_stepper(self, event.position().toPoint())
        )
        super().mousePressEvent(event)

    def mouseReleaseEvent(self, event: QMouseEvent) -> None:  # noqa: N802 (Qt API)
        super().mouseReleaseEvent(event)
        if self._clear_focus_on_stepper_release:
            self._clear_focus_on_stepper_release = False
            self.clearFocus()


class NoWheelDoubleSpinBox(QDoubleSpinBox):
    """Do not change value when the scroll wheel passes over the widget."""

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self._clear_focus_on_stepper_release = False

    def wheelEvent(self, event: QWheelEvent) -> None:  # noqa: N802 (Qt API)
        event.ignore()

    def mousePressEvent(self, event: QMouseEvent) -> None:  # noqa: N802 (Qt API)
        self._clear_focus_on_stepper_release = (
            event.button() == Qt.MouseButton.LeftButton
            and _hit_spin_stepper(self, event.position().toPoint())
        )
        super().mousePressEvent(event)

    def mouseReleaseEvent(self, event: QMouseEvent) -> None:  # noqa: N802 (Qt API)
        super().mouseReleaseEvent(event)
        if self._clear_focus_on_stepper_release:
            self._clear_focus_on_stepper_release = False
            self.clearFocus()


class NoWheelSlider(QSlider):
    """Ignore wheel events so scrolling the panel does not change the value.

    Groove clicks jump the handle to the clicked position (AbsoluteSet), so the
    right end of a zoom slider is maximum zoom — not a page-step nudge.
    """

    def wheelEvent(self, event: QWheelEvent) -> None:  # noqa: N802 (Qt API)
        event.ignore()

    def mousePressEvent(self, event: QMouseEvent) -> None:  # noqa: N802 (Qt API)
        if event.button() == Qt.MouseButton.LeftButton and not self.isSliderDown():
            opt = QStyleOptionSlider()
            self.initStyleOption(opt)
            style = self.style()
            handle = style.subControlRect(
                QStyle.ComplexControl.CC_Slider,
                opt,
                QStyle.SubControl.SC_SliderHandle,
                self,
            )
            pos = event.position().toPoint()
            if not handle.contains(pos):
                groove = style.subControlRect(
                    QStyle.ComplexControl.CC_Slider,
                    opt,
                    QStyle.SubControl.SC_SliderGroove,
                    self,
                )
                if self.orientation() == Qt.Orientation.Horizontal:
                    span = max(1, groove.width() - handle.width())
                    x = pos.x() - groove.x() - handle.width() // 2
                    value = QStyle.sliderValueFromPosition(
                        self.minimum(),
                        self.maximum(),
                        max(0, min(span, x)),
                        span,
                        opt.upsideDown,
                    )
                else:
                    span = max(1, groove.height() - handle.height())
                    y = pos.y() - groove.y() - handle.height() // 2
                    value = QStyle.sliderValueFromPosition(
                        self.minimum(),
                        self.maximum(),
                        max(0, min(span, y)),
                        span,
                        opt.upsideDown,
                    )
                self.setSliderPosition(value)
                self.setValue(value)
                event.accept()
                # Continue so dragging the handle from this point still works.
        super().mousePressEvent(event)
