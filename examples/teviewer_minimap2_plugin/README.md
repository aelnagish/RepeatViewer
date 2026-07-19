# Minimap2 plugin for TEviewer

Map selected FASTA files to a reference genome with [minimap2](https://github.com/lh3/minimap2).

## Platforms

| OS | minimap2 source |
|----|-----------------|
| **Linux** | Bundled binary in `src/teviewer_minimap2/binaries/linux/` |
| **Windows** | Bundled binary in `src/teviewer_minimap2/binaries/win/` |
| **macOS** | Homebrew: `brew install minimap2` |

## Install

Bundled with TEviewer when you install from the repository root:

```bash
pip install -e .
```

Or install only the plugin package:

```bash
cd examples/teviewer_minimap2_plugin
pip install -e .
```

Restart TEviewer. **Minimap2** appears on the function toolbar after the plugin separator.

## Usage

1. Set **Project genome** (File → Project settings), or choose a reference FASTA in the Minimap2 dialog.
2. Select one or more **FASTA** files in the Files panel.
3. Click **Minimap2** on the toolbar.
4. Choose preset (`map-ont`, `map-hifi`, …), threads, and output format (PAF or SAM).
5. Output is written next to each input: `{name}.minimap2.paf` (or `.sam`).

## Replace bundled binaries

Copy new builds into:

- `examples/teviewer_minimap2_plugin/src/teviewer_minimap2/binaries/linux/minimap2`
- `examples/teviewer_minimap2_plugin/src/teviewer_minimap2/binaries/win/minimap2.exe`

Then reinstall: `pip install -e .`
