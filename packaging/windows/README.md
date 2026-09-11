# Windows packaging (folder + zip)

Double-clickable **RepeatViewer.exe** in a portable folder, wrapped in a `.zip`.

## Build

Must run **on Windows** (PyInstaller does not cross-compile a Windows `.exe` from macOS/Linux).

From the repo root, in an environment where `pip install -e .` already works:

```powershell
powershell -ExecutionPolicy Bypass -File packaging/windows/build_windows.ps1
```

Outputs:

- `dist/RepeatViewer/RepeatViewer.exe` (+ supporting files)
- `dist/RepeatViewer-<version>-windows-x64.zip`

## Install (end users)

1. Unzip the archive.
2. Double-click **RepeatViewer.exe**.

SmartScreen may warn on unsigned builds: **More info → Run anyway**.

## What’s included

- PySide6 UI + Qt WebEngine (HTML previews)
- Bundled **Windows-only** BLAST+, MAFFT, MUSCLE (mac/linux tool trees are not packed)
- Flexidot (via the frozen app re-entry flag)
- AliView runtime (`aliview.jar`); installer scripts are omitted

macOS/Linux tool binaries are omitted from the Windows zip to save space.

## Optional tools via WSL2 conda

RepeatMasker, FastGA, bedtools, and similar bioconda tools are **Linux** binaries. On Windows, RepeatViewer keeps the same **External tools (conda)** preference, but resolves it against **conda inside WSL2** and launches those tools with `wsl`.

1. Install [WSL2](https://learn.microsoft.com/windows/wsl/install) (e.g. Ubuntu).
2. Inside WSL, install conda/mamba and your tools env, e.g.:

   ```bash
   conda create -n tetrimmer -c bioconda repeatmasker fastga bedtools -y
   ```

3. In RepeatViewer: **File → Setting → External tools (conda / WSL2)** → set `tetrimmer` (or `~/miniconda3/envs/tetrimmer`).
4. Click **Test** — status should show `Using WSL2 …/bin (like …)`.

Bundled BLAST/MAFFT/MUSCLE stay on Windows and are not routed through WSL.
