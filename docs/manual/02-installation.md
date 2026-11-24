# Installation

RepeatViewer runs on macOS, Linux and Windows. There are two ways to install
it:

- the **ready-made app**, which needs no Python and is the quickest way to
  start;
- the **source code**, if you want to change RepeatViewer or prefer a conda or
  pip installation.

## Ready-made app

Download the file for your system from the
[latest release](https://github.com/qjiangzhao/RepeatViewer/releases/latest):

| System | File name ends in | To install |
|--------|-------------------|------------|
| macOS, Apple Silicon (M1 and later) | `macos-arm64.dmg` | Open the `.dmg` and drag RepeatViewer into Applications |
| macOS, Intel | `macos-x86_64.dmg` | Open the `.dmg` and drag RepeatViewer into Applications |
| Windows, 64-bit | `windows-x64.zip` | Unzip, then start `RepeatViewer.exe` |
| Linux, x86_64 | `linux-x86_64.tar.gz` | Unpack, then run `./RepeatViewer` in the unpacked folder |

The apps are not code-signed, so your system asks for confirmation the first
time you start RepeatViewer:

- **macOS 15 or newer**: open the app once (macOS refuses), then go to
  **System Settings → Privacy & Security**, scroll down to **Security** and
  click **Open Anyway**.
- **macOS 14 or older**: right-click the app in Applications and choose
  **Open**.
- **Windows**: click **More info**, then **Run anyway**.

The app already contains Python and the programs listed under **What is
bundled** below. The optional external programs further down still need a
conda environment.

## From source

Installing from source needs **Python 3.10 or newer**. First get the source
code, in one of two ways:

- **From a release.** On the
  [releases page](https://github.com/qjiangzhao/RepeatViewer/releases), open
  the latest release and download **Source code (zip)** or **Source code
  (tar.gz)** under **Assets**. Unpack it. The folder carries the version in
  its name, for example `RepeatViewer-1.2.0`.
- **With git.** `git clone https://github.com/qjiangzhao/RepeatViewer.git`
  makes a folder called `RepeatViewer`, which you can update later with
  `git pull`.

### macOS and Linux

Open **Terminal**, go to the source folder (here the one made by git) and run:

```bash
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

From the source folder, `python RepeatViewer.py` works too.

`pip install -e .` also installs the RepeatViewer build of Flexidot, which makes
the interactive HTML dotplots. The Flexidot on PyPI only writes images.

### Windows

Install with the normal Windows Python (not inside WSL), in **PowerShell** (not
Command Prompt, not Ubuntu/WSL). Put RepeatViewer in a folder you can write to,
such as `Documents`, not under `C:\Program Files`.

1. Install [Python 3.10 or later](https://www.python.org/downloads/) from
   **python.org**. Keep **py launcher** ticked and tick **Add python.exe to
   PATH**. Other Pythons (from Rtools, MSYS2, Git for Windows, Cygwin) do not
   work: PySide6 cannot be installed in their environments.
2. Install the Microsoft C++ Build Tools, needed by one dependency (orfipy).
3. Download **Source code (zip)** from the latest release (see above) and unzip
   it, for example to `Documents\RepeatViewer-1.2.0`. The green
   **Code → Download ZIP** button on the repository page works too; its folder
   is called `RepeatViewer-main`.
4. In PowerShell, go to that folder, create an environment called `.venv`, and
   install.

The commands for steps 2 to 4 (in the `cd` line, use the name of the folder you
unzipped):

```powershell
winget install -e --id Microsoft.VisualStudio.BuildTools --override "--passive --wait --add Microsoft.VisualStudio.Workload.VCTools --includeRecommended"
```

```powershell
cd $HOME\Documents\RepeatViewer-1.2.0
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
| pysam | Large BAM / CRAM genome browser and **Tools → Samtools** | Included with `pip install -e .` on macOS / Linux. Windows: the [Win-NGS samtools bundle](https://github.com/win-ngs/samtools-windows-bundle/releases) (no pysam wheels) |
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

**Ready-made app.** Download the new release and replace the old copy: on
macOS, drag the new app into Applications and choose **Replace**; on Windows
and Linux, unpack the new archive and use that folder from now on.

**Source from git.** In the source folder:

```bash
git pull
pip install -e .
```

**Source from a release.** Download and unpack the new source, then run
`pip install -e .` in the new folder with the same environment active. On
Windows the `.venv` environment sits inside the source folder, so repeat steps
3 and 4 in the new folder.

Run `pip install -e .` after every source update, because new versions
sometimes add Python packages.
