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
    assert "FastGA 1.5" in rich
    assert "<b>Supported version: FastGA 1.5 or newer</b>" in rich
    assert "FastGA 1.5" in plain


def test_fastga_configuration_message_is_platform_specific(monkeypatch) -> None:
    monkeypatch.setattr(sys, "platform", "win32")
    win = fastga_configuration_message(rich=True)
    assert "WSL2" in win
    assert "bioconda" in win

    monkeypatch.setattr(sys, "platform", "darwin")
    mac = fastga_configuration_message(rich=True)
    assert "brew install fastga" in mac
    assert "conda install -c bioconda fastga" in mac
    assert "WSL2" not in mac

    monkeypatch.setattr(sys, "platform", "linux")
    linux = fastga_configuration_message(rich=True)
    assert "conda install -c bioconda fastga" in linux
    assert "thegenemyers/FASTGA" in linux
    assert "WSL2" not in linux
