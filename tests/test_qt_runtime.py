"""Qt / WebEngine runtime flags and stderr noise filtering."""

from __future__ import annotations

import os

import pytest

from repeatviewer import qt_runtime


def test_webengine_flags_disable_gpu_on_darwin(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(qt_runtime.sys, "platform", "darwin")
    monkeypatch.delenv("QTWEBENGINE_CHROMIUM_FLAGS", raising=False)
    # Avoid creating the dictionaries directory under site-packages during tests
    # when PySide6 is present — still exercise flag assembly.
    monkeypatch.setenv("QTWEBENGINE_DICTIONARIES_PATH", "/tmp/rv-test-dicts")
    qt_runtime._configure_webengine()
    flags = os.environ["QTWEBENGINE_CHROMIUM_FLAGS"]
    for required in (
        "--disable-gpu",
        "--disable-gpu-compositing",
        "--in-process-gpu",
        "--use-gl=disabled",
        "--log-level=3",
    ):
        assert required in flags, f"missing {required} in {flags!r}"


def test_stderr_filter_drops_gpuinfo_noise() -> None:
    class _Sink:
        def __init__(self) -> None:
            self.chunks: list[str] = []

        def write(self, text: str) -> int:
            self.chunks.append(text)
            return len(text)

        def flush(self) -> None:
            return None

        def fileno(self) -> int:
            return -1

        def isatty(self) -> bool:
            return False

    sink = _Sink()
    filtered = qt_runtime._QtNoiseStderrFilter(sink)
    filtered.write("GPUInfo not initialized on GpuInfoUpdate\n")
    filtered.write("real error line\n")
    filtered.flush()
    assert sink.chunks == ["real error line\n"]


def test_macos_fd_skip_includes_gpuinfo() -> None:
    blobs = b" ".join(qt_runtime._MACOS_STDERR_SKIP)
    assert b"GPUInfo not initialized" in blobs
    assert b"GpuInfoUpdate" in blobs
