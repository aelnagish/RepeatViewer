"""Unit tests for frozen-friendly Flexidot invocation."""

from __future__ import annotations

import subprocess

import pytest

from repeatviewer.backends import gui_functions


def test_run_flexidot_inprocess_help(monkeypatch: pytest.MonkeyPatch) -> None:
    def fake_main():
        print("usage: flexidot [-h] -f {png,pdf,svg,html}")
        raise SystemExit(0)

    monkeypatch.setattr(
        "flexidot.app.main",
        fake_main,
        raising=False,
    )
    # Import path used inside _run_flexidot_inprocess
    import flexidot.app as flexidot_app

    monkeypatch.setattr(flexidot_app, "main", fake_main)

    result = gui_functions._run_flexidot_inprocess(["-h"])
    assert result.returncode == 0
    assert "html" in result.stdout.lower()


def test_run_flexidot_argv_uses_inprocess_when_frozen(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(gui_functions.sys, "frozen", True, raising=False)
    calls: list[list[str]] = []

    def fake_inprocess(argv: list[str]) -> subprocess.CompletedProcess[str]:
        calls.append(list(argv))
        return subprocess.CompletedProcess(["flexidot", *argv], 0, "", "")

    monkeypatch.setattr(gui_functions, "_run_flexidot_inprocess", fake_inprocess)
    monkeypatch.setattr(
        gui_functions,
        "run_subprocess_with_cancel",
        lambda *a, **k: (_ for _ in ()).throw(AssertionError("subprocess used")),
    )

    completed = gui_functions._run_flexidot_argv(["-i", "x.fa", "-m", "0", "-k", "5"])
    assert completed.returncode == 0
    assert calls == [["-i", "x.fa", "-m", "0", "-k", "5"]]


def test_flexidot_supports_html_true_when_frozen(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(gui_functions.sys, "frozen", True, raising=False)
    assert gui_functions.flexidot_supports_html() is True
    assert gui_functions.flexidot_supports_repeatviewer_integration() is True
