import os
import shutil
import subprocess
import sys

sys.path.insert(0, r"C:\Users\QJIAN\Personal_folder\RepeatViewer\RepeatViewer\RepeatViewer")

print("HOME=", os.environ.get("HOME"))
print("USERPROFILE=", os.environ.get("USERPROFILE"))

cmd = [
    "wsl",
    "-e",
    "bash",
    "-lc",
    "export PATH=/home/qian/miniconda3/envs/TETrimmer/bin:$PATH; "
    "command -v FastGA; "
    "file /home/qian/miniconda3/envs/TETrimmer/bin/FastGA; "
    "FastGA 2>&1 | head -8",
]
print("running", cmd)
r = subprocess.run(
    cmd,
    capture_output=True,
    text=True,
    timeout=30,
    encoding="utf-8",
    errors="replace",
)
print("rc", r.returncode)
print("stdout", r.stdout)
print("stderr", r.stderr[:800] if r.stderr else None)

print("--- HOME inside wsl -lc ---")
r = subprocess.run(
    ["wsl", "-e", "bash", "-lc", "echo HOME=$HOME; echo USER=$USER; id -un; getent passwd $(id -un)"],
    capture_output=True,
    text=True,
    timeout=20,
    encoding="utf-8",
    errors="replace",
)
print("stdout", repr(r.stdout))
print("stderr", r.stderr[:300] if r.stderr else "")

print("--- conda_tools ---")
from repeatviewer.conda_tools import (
    describe_conda_tools_status,
    resolve_wsl_conda_bin,
    uses_wsl_conda,
    which_tool,
    wsl_available,
)

spec = "/home/qian/miniconda3/envs/TETrimmer"
print("wsl_available", wsl_available())
print("resolve_wsl_conda_bin", resolve_wsl_conda_bin(spec))
print("uses_wsl_conda", uses_wsl_conda(spec))
print("which FastGA", which_tool("FastGA", spec=spec))
print("which fastga", which_tool("fastga", spec=spec))
print("status", describe_conda_tools_status(spec))
print("shutil.which FastGA", shutil.which("FastGA"))
print("env name resolve", resolve_wsl_conda_bin("TETrimmer"))
print("status env name", describe_conda_tools_status("TETrimmer"))
