"""Trackpad pinch and Shift+scroll zoom for the multi-sequence / MSA overview."""

from __future__ import annotations

from pathlib import Path

import pytest
from PySide6.QtCore import QPoint, QPointF, Qt
from PySide6.QtGui import QNativeGestureEvent, QPointingDevice, QWheelEvent
from PySide6.QtWidgets import QApplication

from repeatviewer.widgets.multi_sequence_overview import (
    MultiSequenceOverviewPanel,
    _WHEEL_ZOOM_PER_NOTCH,
    _pinch_zoom_factor,
    _shift_wheel_zoom_factor,
)


def _pinch(widget, value: float, kind=Qt.NativeGestureType.ZoomNativeGesture):
    local = QPointF(widget.width() / 2.0, widget.height() / 2.0)
    glob = QPointF(widget.mapToGlobal(local.toPoint()))
    ev = QNativeGestureEvent(
        kind,
        QPointingDevice.primaryPointingDevice(),
        2,
        local,
        local,
        glob,
        value,
        QPointF(0, 0),
    )
    QApplication.sendEvent(widget, ev)
    return ev


def _wheel(angle_y: int) -> QWheelEvent:
    return QWheelEvent(
        QPointF(10, 10),
        QPointF(10, 10),
        QPoint(0, 0),
        QPoint(0, angle_y),
        Qt.MouseButton.NoButton,
        Qt.KeyboardModifier.ShiftModifier,
        Qt.ScrollPhase.NoScrollPhase,
        False,
    )


@pytest.fixture(params=[True, False], ids=["msa", "multi"])
def panel(qapp: QApplication, tmp_path: Path, request):
    path = tmp_path / "aln.fa"
    seq = "ACGT" * 500
    if request.param:
        path.write_text(f">a\n{seq}\n>b\n{seq}\n")
    else:
        path.write_text(f">a\n{seq}\n>b\n{seq[:1200]}\n")
    p = MultiSequenceOverviewPanel(path)
    p.resize(900, 500)
    p.load_path(path, msa_mode=request.param)
    p.show()
    qapp.processEvents()
    p._bars.fit_to_bars_width(p._bars_viewport_width())
    yield p
    p.deleteLater()


def test_pinch_out_zooms_in_and_pinch_in_zooms_out(panel) -> None:
    fit = panel._bars.pixels_per_base()
    for target in (panel._bars, panel._names):
        before = panel._bars.pixels_per_base()
        ev = _pinch(target, 0.1)
        assert ev.isAccepted()
        assert panel._bars.pixels_per_base() == pytest.approx(
            before * _pinch_zoom_factor(0.1)
        )
    _pinch(panel._bars, -5.0)
    _pinch(panel._bars, -5.0)
    assert panel._bars.pixels_per_base() == pytest.approx(fit)


def test_smart_zoom_toggles_fit(panel) -> None:
    fit = panel._bars.pixels_per_base()
    _pinch(panel._bars, 0.0, Qt.NativeGestureType.SmartZoomNativeGesture)
    assert panel._bars.pixels_per_base() > fit * 2
    _pinch(panel._bars, 0.0, Qt.NativeGestureType.SmartZoomNativeGesture)
    assert panel._bars.pixels_per_base() == pytest.approx(fit)


def test_shift_wheel_zoom_scales_with_delta() -> None:
    assert _shift_wheel_zoom_factor(_wheel(120)) == pytest.approx(
        _WHEEL_ZOOM_PER_NOTCH
    )
    assert _shift_wheel_zoom_factor(_wheel(-120)) == pytest.approx(
        1 / _WHEEL_ZOOM_PER_NOTCH
    )
    small = _shift_wheel_zoom_factor(_wheel(12))
    assert 1.0 < small < 1.05
    assert _shift_wheel_zoom_factor(_wheel(0)) is None
