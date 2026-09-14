"""Windows / WSL path mapping for conda-backed CLI tools."""

from __future__ import annotations

from repeatviewer.conda_tools import (
    _conda_lib_dir,
    _map_windows_path_to_linux,
    _rewrite_arg_for_wsl,
    _unc_wsl_path_to_linux,
    _windows_drive_to_wsl,
)


def test_unc_wsl_path_to_linux() -> None:
    assert (
        _unc_wsl_path_to_linux(
            r"\\wsl$\Ubuntu\home\qian\miniconda3\envs\TETrimmer"
        )
        == "/home/qian/miniconda3/envs/TETrimmer"
    )
    assert (
        _unc_wsl_path_to_linux(
            r"\\wsl.localhost\Ubuntu\home\qian\miniconda3\envs\TETrimmer"
        )
        == "/home/qian/miniconda3/envs/TETrimmer"
    )
    assert _unc_wsl_path_to_linux(r"C:\Users\QJIAN") is None


def test_windows_drive_to_wsl() -> None:
    assert (
        _windows_drive_to_wsl(r"C:\Users\QJIAN\file.fa")
        == "/mnt/c/Users/QJIAN/file.fa"
    )
    assert (
        _windows_drive_to_wsl("C:/Users/QJIAN/file.fa")
        == "/mnt/c/Users/QJIAN/file.fa"
    )
    assert _windows_drive_to_wsl("C:\\") == "/mnt/c"
    assert _windows_drive_to_wsl("/home/qian") is None


def test_rewrite_fastga_glued_flags() -> None:
    assert (
        _rewrite_arg_for_wsl(r"-PC:\Users\QJIAN\tmp")
        == "-P/mnt/c/Users/QJIAN/tmp"
    )
    assert (
        _rewrite_arg_for_wsl(r"-1:C:\Users\QJIAN\out.1aln")
        == "-1:/mnt/c/Users/QJIAN/out.1aln"
    )
    assert _rewrite_arg_for_wsl("-i0.700") == "-i0.700"
    assert (
        _rewrite_arg_for_wsl(r"C:\Users\QJIAN\seq.fa")
        == "/mnt/c/Users/QJIAN/seq.fa"
    )


def test_map_prefers_unc_over_drive() -> None:
    assert (
        _map_windows_path_to_linux(
            r"\\wsl$\Ubuntu\home\qian\miniconda3\envs\TETrimmer"
        )
        == "/home/qian/miniconda3/envs/TETrimmer"
    )


def test_conda_lib_dir() -> None:
    assert (
        _conda_lib_dir("/home/qian/miniconda3/envs/TETrimmer/bin")
        == "/home/qian/miniconda3/envs/TETrimmer/lib"
    )
