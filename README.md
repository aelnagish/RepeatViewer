# RepeatViewer (TEviewer)

**RepeatViewer** is a desktop application for working with transposable element (TE) consensus sequences: browse folders of FASTA and tool outputs, run analyses from a toolbar, and inspect results in an integrated preview pane.

Built with **PySide6**, it targets the manual curation workflow—proofreading consensus sequences, checking terminal repeats, comparing dotplots, running BLAST against a reference genome, and cleaning alignments before publication or database submission.

The app name shown in the window is **RepeatViewer**. The installable Python package and command-line entry point are **`teviewer`**.

## What you can do

- Open any folder tree of TE FASTA files, BLAST tables, MSAs, GFF annotations, and HTML reports
- Define a single reference genome or a **pan-genome** for BLAST and chromosome-scale plots
- Run curation and annotation tools without leaving the app; outputs land **beside the source FASTA** in the same folder
- Preview sequences (text, annotation map, self dotplot), Plotly HTML, PDFs, and images inline
- Detach Folders, Files, or Preview to a second monitor; **View → Reset layout** restores the default docking

RepeatViewer works well on TEtrimmer output folders, curated consensus libraries, or any organized set of TE sequences you already have on disk.

## Browser workspace

Three docked panels form the main interface:

| Panel | Purpose |
|-------|---------|
| **Folders** | Workspace roots and subfolders; back / forward / up; filter by name |
| **Files** | Sortable table (name, length, GC%, sequence count); filter by name; optional columns |
| **Preview** | FASTA tabs (Text · Annotation · Dotplot), HTML reports, alignments, PDFs, images |

Set the reference under **File → Define genome** before BLAST, genome maps, or genome-aware TEAid plots.

## Typical workflow

1. **Add folders** — drag a workspace root into the Folders panel or use the folder context menu.
2. **Select a consensus FASTA** — the Files panel lists folder contents; click a file to preview it.
3. **Define genome** — **File → Define genome** (single FASTA or pan-genome list).
4. **Run a tool** — use the top toolbar or **Sequence** menu (BLAST, Align, TEAid, Flexidot, …).
5. **Review outputs** — new files appear in the same folder; recent results are highlighted in green or blue.
6. **Iterate** — extend ends, clean the MSA, rebuild consensus, re-BLAST, compare dotplots.

Background jobs show progress in the status bar; **Stop** cancels the active task.

## Analysis tools

Available from the toolbar and **Sequence** menu:

| Category | Capabilities |
|----------|--------------|
| **BLAST** | Bundled BLASTN and BLASTP; nucleotide/protein databases; hit tables; BLASTN GFF and interactive genome maps |
| **Align** | Bundled MAFFT and MUSCLE; configurable per-aligner settings, default engine selection, and combined multi-sequence alignment |
| **MSA clean** | Remove gappy columns; trim divergent rows; combined clean |
| **Consensus** | Build consensus from an alignment |
| **TEAid** | Interactive HTML TE annotation report |
| **Flexidot** | Self, paired, and multi-sequence dotplots; in-app HTML dotplots with **Go-to → Annotation** |
| **Extend** | Extend alignment ends (preset or custom distances) |
| **TE visualization** | Genome map and TE divergence plot |
| **Tree** | Phylogenetic tree from alignment |
| **FastGA** | FastGA alignment workflow |
| **TSD search** | Target site duplication search |
| **Minimap2** | Minimap2 alignment when available |
| **Benchmarking** | Confusion-matrix benchmark against reference classification tables |
| **Plugins** | Install custom command plugins (**Plugins** menu) |

Bundled **BLAST+**, **MAFFT**, and **MUSCLE** binaries are included for macOS, Linux, and Windows. The corresponding **Use system …** options in **File → Setting** select an installation from `PATH` instead.

## Install

**Requirements:** Python 3.10+, conda or venv recommended.

```bash
git clone https://github.com/qjiangzhao/TEviewer.git
cd TEviewer
conda create -n teviewer python=3.10 -y
conda activate teviewer
pip install -e .
```

### Optional external programs

| Program | Role |
|---------|------|
| Reference **genome FASTA** | BLAST, extend, genome maps, TEAid genome tracks |
| **AliView** | Open selected FASTA externally (**Sequence → Open in AliView**) |

Flexidot (with HTML output and RepeatViewer integration) is installed as a dependency of `pip install -e .`.

## Run

```bash
teviewer
```

Or from the repository root:

```bash
python TEviewer.py
```

Open a folder and genome immediately:

```bash
teviewer -i /path/to/consensus/folder -g /path/to/genome.fa
```

| Flag | Meaning |
|------|---------|
| `-i`, `--input-dir` | Workspace folder opened on launch |
| `-g`, `--genome` | Reference genome FASTA |

Workspace roots, window layout, genomes, and preferences are saved between sessions.

## Shortcuts

| Key | Action |
|-----|--------|
| **Ctrl+F** / **⌘F** | Focus file or folder filter (depends on active panel) |
| **F5** | Refresh current folder |
| **Delete** | Delete selected files (with confirmation) |

## Project structure

```
TEviewer/
  TEviewer.py              # entry script
  teviewer/
    main_window.py         # browser UI, menus, job runner
    project.py             # genome, pan-genome, database settings
    backends/              # analysis engines (BLAST, TEAid, Flexidot, …)
    views/                 # tool and settings dialogs
    widgets/               # folder tree, file table, preview tabs
    plugins/               # plugin API
  examples/                # sample plugin
```

## Companion to TEtrimmer

[TEtrimmer](https://github.com/qjiangzhao/TEtrimmer) automates TE discovery and consensus building from genome assemblies. **RepeatViewer** is the interactive companion for what comes next: inspecting candidates, fixing consensus sequences, running targeted BLAST and annotation, and preparing libraries for release.

The two tools are independent applications; RepeatViewer does not require a TEtrimmer installation.

## Citation

If RepeatViewer contributes to your work, cite **TEtrimmer** (*Nature Communications*, 2025) when applicable and mention **RepeatViewer (TEviewer)** for manual curation and visualization.
