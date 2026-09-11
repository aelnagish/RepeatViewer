# Linux packaging (folder + tar.gz)

Portable **RepeatViewer** binary folder, wrapped in a `.tar.gz`.

## Build

Must run **on Linux**, matching the CPU you will distribute for (`x86_64` or `aarch64`).

From the repo root, in an environment where `pip install -e .` already works:

```bash
bash packaging/linux/build_linux.sh
```

Outputs:

- `dist/RepeatViewer/RepeatViewer`
- `dist/RepeatViewer-<version>-linux-<arch>.tar.gz`

## Install (end users)

```bash
tar -xzf RepeatViewer-*-linux-*.tar.gz
cd RepeatViewer
./RepeatViewer
```

Optional: edit `RepeatViewer.desktop` (`Exec=` / `Icon=` absolute paths) and copy it to `~/.local/share/applications/`.

## What’s included

- PySide6 UI + Qt WebEngine (HTML previews)
- Bundled Linux BLAST+, MAFFT, MUSCLE
- Flexidot (via the frozen app re-entry flag)

macOS/Windows tool binaries are omitted to save space.

**Note:** Bundled MAFFT is x86_64-oriented; MUSCLE ships both `x86_64` and `aarch64`. Build on the arch you intend to ship.

## Optional conda tools

Set **File → Setting → External tools (conda)** to an env name (e.g. `tetrimmer`) or the env folder path. RepeatViewer prepends that env’s `bin` to `PATH` so RepeatMasker / FastGA / bedtools work from the packed binary.
