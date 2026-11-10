# Installation

RepeatViewer runs on macOS, Linux and Windows and needs **Python 3.10 or
newer**. Installing from the source code is the recommended way for now;
ready-made packages for each system are in preparation.

## macOS and Linux

Open **Terminal** and run:

```bash
git clone https://github.com/qjiangzhao/RepeatViewer.git
cd RepeatViewer
conda create -n repeatviewer "python>=3.10" -y
conda activate repeatviewer
pip install -e .
```

Start RepeatViewer with:

```bash
conda activate repeatviewer
repeatviewer
```

From the repository folder, `python RepeatViewer.py` works too.

`pip install -e .` also installs the RepeatViewer build of Flexidot, which makes
the interactive HTML dotplots. The Flexidot on PyPI only writes images.

## Windows

Install with the normal Windows Python (not inside WSL), in **PowerShell** (not
Command Prompt, not Ubuntu/WSL). Put RepeatViewer in a folder you can write to,
such as `Documents`, not under `C:\Program Files`.

1. Install [Python 3.10 or later](https://www.python.org/downloads/) from
   **python.org**. Keep **py launcher** ticked and tick **Add python.exe to
   PATH**. Other Pythons (from Rtools, MSYS2, Git for Windows, Cygwin) do not
   work: PySide6 cannot be installed in their environments.
2. Install the Microsoft C++ Build Tools, needed by one dependency (orfipy).
3. Download the repository from GitHub (**Code → Download ZIP**) and unzip it,
   for example to `Documents\RepeatViewer-main`.
4. In PowerShell, go to that folder, create an environment called `.venv`, and
   install.

The commands for steps 2 to 4:

```powershell
winget install -e --id Microsoft.VisualStudio.BuildTools --override "--passive --wait --add Microsoft.VisualStudio.Workload.VCTools --includeRecommended"
```

```powershell
cd $HOME\Documents\RepeatViewer-main
py -0p
py -3 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -U pip
.\.venv\Scripts\python.exe -m pip install -e .
```

`py -0p` lists the python.org Pythons the launcher can use; if it is empty,
install Python first. Start RepeatViewer with:

```powershell
.\.venv\Scripts\python.exe RepeatViewer.py
```

> **Do not name the environment `repeatviewer`.** The project already has a
> folder called `repeatviewer` with the program code; an environment with that
> name is written into it, and deleting the environment deletes RepeatViewer.

If `.venv\Scripts\python.exe` does not exist after step 4, the environment was
made by a non-Windows Python. Remove the `.venv` folder, install Python from
python.org, and try again. The README has the full Windows notes.

## What is bundled

These programs ship inside RepeatViewer for macOS, Linux and Windows. You do not
need to install them:

| Program | Used by |
|---------|---------|
| BLAST+ (blastn, blastp, blastx, makeblastdb, rpstblastn) | BLAST, TE-Aid, Repeat Quantify, annotation |
| MAFFT | Align, Tree, boundary prediction |
| MUSCLE | Align |
| minimap2 2.31 | Mapping |
| IQ-TREE 3.1.4 | Tree |

Python libraries installed with RepeatViewer cover the rest: Primer3
(primer design), pyhmmer (HMMER searches), pytrf (tandem repeats), orfipy
(ORFs), Biopython.

## Optional external programs

Some tools call programs that are too large to bundle or are only needed by a
few people. Install the ones you need into a conda environment:

| Program | Needed for | Install |
|---------|------------|---------|
| RepeatMasker | Repeat Quantify with the RepeatMasker engine | `conda install -c bioconda -c conda-forge repeatmasker` |
| FastGA 1.5 or newer | Align → FastGA | `conda install -c bioconda fastga=1.5.1` |
| HISAT2 2.2 or newer | Mapping short reads | `conda install -c conda-forge -c bioconda hisat2=2.2.3` |
| FastME, FastTree | Tree with these methods | `conda install -c bioconda fastme fasttree` |
| pysam or samtools | **Tools → Samtools**, sorting and indexing mapping results | `pip install pysam` (macOS, Linux); Windows: the [Win-NGS samtools bundle](https://github.com/win-ngs/samtools-windows-bundle/releases) |
| bedtools | Some extraction steps | `conda install -c bioconda bedtools` |
| AliView | **Tools → Open in AliView** | [ormbunkar.se/aliview](https://ormbunkar.se/aliview/) |
| Chorus2 | Chromosome painting probes | see [Primers and PCR](12-primers.md) |

You can install them into the same `repeatviewer` environment, or keep them in
a separate one. Either way, tell RepeatViewer where they are:
**File → Setting → General → External tools**, enter the environment name
(for example `repeatviewer`) or browse to the environment folder, and click
**Test** to see which programs are found.

- **macOS and Linux**: RepeatViewer adds that environment's `bin` folder to the
  search path of every tool it runs. This matters most when RepeatViewer is not
  started from a terminal, because then your shell's conda setup is not loaded.
- **Windows**: the setting points at a conda environment **inside WSL2**. Linux
  programs (RepeatMasker, FastGA) run through `wsl`; the bundled BLAST, MAFFT,
  MUSCLE, minimap2 and IQ-TREE run natively.

Most tool windows also have a **Use … from the system PATH or defined Conda
environment** box on their **Options** tab, if you prefer your own copy of a
bundled program.

## Updating

From the repository folder:

```bash
git pull
pip install -e .
```

Run `pip install -e .` again after every update, because new versions sometimes
add Python packages.
