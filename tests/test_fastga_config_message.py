"""FastGA configuration reminder messages."""

from __future__ import annotations

import sys

from repeatviewer.backends.fastga_dotplot import (
    fastga_configuration_message,
    install_hint,
)


def test_fastga_configuration_message_highlights_supported_version() -> None:
    rich = fastga_configuration_message(rich=True)
    plain = install_hint()
    assert "<b>Supported version: FastGA 1.5 or newer</b>" in rich
    assert "Supported version: FastGA 1.5 or newer" in plain
    assert "Open a Terminal window" in rich
    assert "Install FastGA with conda into RepeatViewer" in rich
    assert "<b>If</b>" in rich
    assert "<b>else:</b>" in rich
    assert "<b>Then</b>, in RepeatViewer:" in rich
    assert "rv_env" in rich
    assert "fastga=1.5.1" in rich
    assert "conda install -n rv_env -c bioconda fastga=1.5.1" in rich
    assert "conda create -n rv_env -c bioconda fastga=1.5.1" in rich
    assert "repeatviewer_env" not in rich
    assert "Homebrew" not in rich
    assert rich.count("<b>") == 4  # version + If + else + Then
    assert "External tools" in rich


def test_fastga_configuration_message_is_platform_specific(monkeypatch) -> None:
    monkeypatch.setattr(sys, "platform", "win32")
    win = fastga_configuration_message(rich=True)
    assert "WSL2" in win
    assert "rv_env" in win
    assert "conda install -n rv_env -c bioconda fastga=1.5.1" in win

    monkeypatch.setattr(sys, "platform", "darwin")
    mac = fastga_configuration_message(rich=True)
    assert "Install FastGA with conda into RepeatViewer" in mac
    assert "<b>Then</b>, in RepeatViewer:" in mac
    assert "WSL2" not in mac

    monkeypatch.setattr(sys, "platform", "linux")
    linux = fastga_configuration_message(rich=True)
    assert "conda create -n rv_env -c bioconda fastga=1.5.1" in linux
    assert "<b>Then</b>, in RepeatViewer:" in linux
    assert "WSL2" not in linux
