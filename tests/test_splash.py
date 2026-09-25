"""Start-up banner: tips, "Next tip", closing, and the Settings switch."""

from __future__ import annotations

from PySide6.QtCore import QEvent, QPointF, Qt
from PySide6.QtGui import QMouseEvent
from PySide6.QtWidgets import QWidget

from repeatviewer.widgets import splash as sp


def _click(widget, pos: QPointF) -> None:
    for kind in (QEvent.Type.MouseButtonPress,):
        widget.mousePressEvent(QMouseEvent(kind, pos, widget.mapToGlobal(pos.toPoint()),
                                           Qt.MouseButton.LeftButton, Qt.MouseButton.LeftButton,
                                           Qt.KeyboardModifier.NoModifier))


def test_bundled_tips_and_banner_exist():
    tips = sp.load_tips()
    assert len(tips) >= 10 and not any(t.startswith("#") for t in tips)
    assert sp.BANNER.is_file()


def test_tips_file_comments_and_fallback(tmp_path):
    f = tmp_path / "tips.txt"
    f.write_text("# comment\n\nFirst tip.\n  Second tip.  \n")
    assert sp.load_tips(f) == ["First tip.", "Second tip."]
    assert sp.load_tips(tmp_path / "missing.txt")


def test_next_tip_click_and_close(qapp):
    s = sp.RepeatViewerSplash("1.0.0", ["a", "b", "c"], seed=1)
    s.show()
    s.set_status("Loading tools…")
    s.grab()                                   # paints version, tip and status
    first = s.tip
    _click(s, s._next_rect.center())
    assert s.tip != first and s.isVisible()
    _click(s, QPointF(sp.WIDTH / 2, 60))       # anywhere else closes it
    assert not s.isVisible()


def test_finish_waits_for_minimum_time(qapp):
    window = QWidget()
    window.show()
    s = sp.RepeatViewerSplash("1.0.0", ["a"])
    s.show()
    s.finish_after_minimum(window, min_ms=10_000)
    assert s.isVisible()                       # still up: timer pending
    s.finish_after_minimum(window, min_ms=0)
    qapp.processEvents()
    assert not s.isVisible()


def test_splash_setting_round_trip(qapp):
    from repeatviewer import app_preferences as ap

    ap.save_show_splash(False)
    assert ap.load_show_splash() is False
    ap.save_show_splash(True)
    assert ap.load_show_splash() is True


def test_banner_shows_major_minor_version():
    assert sp.short_version("1.0.0") == "1.0"
    assert sp.short_version("2.3.1") == "2.3"
    assert sp.short_version("2") == "2"
