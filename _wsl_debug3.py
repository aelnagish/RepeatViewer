import subprocess
from pathlib import Path

print("=== wslpath variants ===")
win = r"C:\Users\QJIAN"
for cmd in (
    ["wsl", "wslpath", "-a", "-u", win],
    ["wsl", "wslpath", "-a", "-u", "C:/Users/QJIAN"],
    ["wsl", "-e", "wslpath", "-a", "-u", "C:/Users/QJIAN"],
    ["wsl", "-e", "bash", "-c", "wslpath -a -u 'C:/Users/QJIAN'"],
    ["wsl", "-e", "bash", "-c", r"wslpath -a -u 'C:\Users\QJIAN'"],
):
    r = subprocess.run(cmd, capture_output=True, timeout=10, check=False)
    print(cmd)
    print(" rc", r.returncode, "stdout", r.stdout, "stderr", r.stderr[:200] if r.stderr else None)

print("=== multiline script cat -A ===")
script = r"""
set -e
name=TETrimmer
echo START
for root in "$HOME/miniconda3" "$HOME/anaconda3"; do
  echo "try $root/envs/$name/bin"
  if [ -d "$root/envs/$name/bin" ]; then echo FOUND; echo "$root/envs/$name/bin"; exit 0; fi
done
echo AFTER_LOOP
exit 1
"""
r = subprocess.run(
    ["wsl", "-e", "bash", "-lc", script],
    capture_output=True,
    timeout=20,
    check=False,
    text=True,
    encoding="utf-8",
    errors="replace",
)
print("rc", r.returncode)
print("stdout", repr(r.stdout))
print("stderr", repr(r.stderr))

print("=== bash -c not login ===")
r = subprocess.run(
    ["wsl", "-e", "bash", "-c", script],
    capture_output=True,
    timeout=20,
    check=False,
    text=True,
    encoding="utf-8",
    errors="replace",
)
print("rc", r.returncode)
print("stdout", repr(r.stdout))
