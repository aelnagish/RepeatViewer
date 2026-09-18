# Linux packaging (folder + tar.gz)

Portable **RepeatViewer** binary folder, wrapped in a `.tar.gz`.

## Build

Must run **on Linux**, matching the CPU you will distribute for (`x86_64` or
`aarch64`). PyInstaller does not cross-compile from macOS/Windows.

From the repo root, in an environment where RepeatViewer already runs:

```bash
# 1. Environment with app deps (Python ≥ 3.10)
conda create -n repeatviewer "python>=3.10" -y
conda activate repeatviewer
# or: python3 -m venv venv && source venv/bin/activate

# 2. Install RepeatViewer (pulls orfipy / pyhmmer / pyahocorasick / Flexidot, …)
pip install -e .

# 3. Pack
bash packaging/linux/build_linux.sh
```

Outputs:

- `dist/RepeatViewer/RepeatViewer` — run this binary
- `dist/RepeatViewer-<version>-linux-<arch>.tar.gz` — share this archive

Quick dependency check after packing:

```bash
./dist/RepeatViewer/RepeatViewer --repeatviewer-import-check
```

Expect `OK` lines for `orfipy_core`, `pyhmmer`, and `pytrf`.

## Install (end users)

```bash
tar -xzf RepeatViewer-*-linux-*.tar.gz
cd RepeatViewer
./RepeatViewer
```

Optional desktop shortcut: edit `RepeatViewer.desktop` (`Exec=` / `Icon=` to
**absolute** paths under the unpacked folder), then copy it to
`~/.local/share/applications/`.

## What’s included

- PySide6 UI + Qt WebEngine (HTML previews)
- Bundled **Linux-only** tools for the **build CPU**:
  - BLAST+ and MAFFT: packed on **x86_64** only (no aarch64 builds in-repo)
  - MUSCLE: only `linux/<arch>/` for the host (`x86_64` or `aarch64`)
- Flexidot (in-process when frozen — no second window)
- AliView runtime (`aliview.jar` + launcher); installer scripts are omitted
- Annotation deps **orfipy**, **pyahocorasick**, and **pyhmmer** (full package
  trees, including the host’s `pyhmmer.platform.*` extension)

macOS/Windows tool binaries are never packed.

## Optional conda tools

Set **File → Setting → External tools (conda)** to an env name (e.g.
`repeatviewer`) or the env folder path. RepeatViewer prepends that env’s `bin`
to `PATH` so RepeatMasker / FastGA / bedtools work from the packed binary.

## Notes / known limits

- **Build on the same OS/CPU** you ship for (glibc-linked; very old distros may
  fail to start Qt/WebEngine).
- **aarch64**: no bundled BLAST/MAFFT — install system or conda copies if needed.
- If the GUI fails to start with Qt platform plugin errors, install common X11 /
  Wayland packages for your distro (e.g. `libxcb-*`, `libxkbcommon`, `libegl`).
