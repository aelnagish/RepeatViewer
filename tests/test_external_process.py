from __future__ import annotations

import logging
import sys
import subprocess
import pytest
from pathlib import Path

from repeatviewer.external_process import log_command, log_process_output, run_logged_process


def test_command_and_terminal_output_are_logged(caplog, tmp_path: Path) -> None:
    fasta = tmp_path / "sequence with spaces.fa"
    command = ["external-tool", "--query", str(fasta), "--value", "a b"]

    with caplog.at_level(logging.INFO, logger="repeatviewer.commands"):
        log_command("ExternalTool", command, inputs=[fasta], cwd=tmp_path)
        log_process_output(
            "ExternalTool",
            stdout="progress 1\nprogress 2\n",
            stderr="warning text\n",
            output_path=tmp_path / "result.txt",
        )

    text = caplog.text
    assert "RUN tool=ExternalTool" in text
    assert "sequence with spaces.fa'" in text
    assert "ExternalTool stdout | progress 1" in text
    assert "ExternalTool stdout | progress 2" in text
    assert "ExternalTool stderr | warning text" in text
    assert "ExternalTool output |" in text


def test_result_stdout_can_be_omitted_from_log(caplog) -> None:
    with caplog.at_level(logging.INFO, logger="repeatviewer.commands"):
        log_process_output(
            "MAFFT",
            stdout=">sequence\nACGT\n",
            stderr="alignment progress\n",
            stdout_is_result=True,
        )

    assert ">sequence" not in caplog.text
    assert "MAFFT stderr | alignment progress" in caplog.text


def test_run_logs_process_lifecycle_and_streams(caplog) -> None:
    command = [
        sys.executable,
        "-c",
        "import sys; print('step one', flush=True); "
        "print('warning one', file=sys.stderr, flush=True)",
    ]
    with caplog.at_level(logging.INFO, logger="repeatviewer.commands"):
        result = run_logged_process(command, tool="demo", check=True)

    assert result.returncode == 0
    assert "START tool=demo pid=" in caplog.text
    assert "demo stdout | step one" in caplog.text
    assert "demo stderr | warning one" in caplog.text
    assert "FINISH tool=demo pid=" in caplog.text
    assert "exit=0" in caplog.text
    assert "elapsed=" in caplog.text


def test_failed_run_logs_error_output_and_exit_code(caplog) -> None:
    command = [
        sys.executable,
        "-c",
        "import sys; print('fatal detail', file=sys.stderr, flush=True); sys.exit(7)",
    ]
    with caplog.at_level(logging.INFO, logger="repeatviewer.commands"):
        with pytest.raises(subprocess.CalledProcessError) as caught:
            run_logged_process(command, tool="broken", check=True)

    assert caught.value.returncode == 7
    assert "broken stderr | fatal detail" in caplog.text
    assert "FINISH tool=broken" in caplog.text
    assert "exit=7" in caplog.text
