import subprocess
import sys

sys.path.insert(0, r"C:\Users\QJIAN\Personal_folder\RepeatViewer\RepeatViewer\RepeatViewer")

from repeatviewer.conda_tools import _wsl_run, resolve_wsl_conda_bin

print("=== env name TETrimmer via existing helper ===")
print(resolve_wsl_conda_bin("TETrimmer"))

print("=== dump _wsl_run for env-name script ===")
name = "TETrimmer"
script = f"""
set -e
name={name}
for root in "$HOME/miniconda3" "$HOME/anaconda3" "$HOME/mambaforge" "$HOME/miniforge3" \\
            "$HOME/micromamba" /opt/conda /opt/miniconda3 /opt/anaconda3; do
  if [ -d "$root/envs/$name/bin" ]; then echo "$root/envs/$name/bin"; exit 0; fi
done
if command -v conda >/dev/null 2>&1; then
  base=$(conda info --base 2>/dev/null || true)
  if [ -n "$base" ] && [ -d "$base/envs/$name/bin" ]; then
    echo "$base/envs/$name/bin"
    exit 0
  fi
fi
exit 1
"""
print("script repr first 80:", repr(script[:80]))
print("contains CR", "\\r" in repr(script))
r = _wsl_run(script, timeout=25)
print("rc", r.returncode)
print("stdout", repr(r.stdout))
print("stderr", repr(r.stderr[:500] if r.stderr else r.stderr))

print("=== one-liner env name ===")
one = (
    'for root in "$HOME/miniconda3" "$HOME/anaconda3"; do '
    'if [ -d "$root/envs/TETrimmer/bin" ]; then echo "$root/envs/TETrimmer/bin"; exit 0; fi; '
    "done; echo HOME=$HOME; ls -ld $HOME/miniconda3/envs/TETrimmer/bin; exit 1"
)
r = _wsl_run(one, timeout=20)
print("rc", r.returncode)
print("stdout", repr(r.stdout))
print("stderr", repr((r.stderr or "")[:400]))

print("=== wslpath ===")
r = subprocess.run(
    ["wsl", "wslpath", "-a", "-u", r"C:\Users\QJIAN"],
    capture_output=True,
    timeout=8,
    check=False,
)
print("bytes", r.stdout)
print("rc", r.returncode)
