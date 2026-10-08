# RepeatViewer User Manual

**Version 1.1.0**

RepeatViewer is a desktop application for the manual curation of transposable
element (TE) consensus sequences. It gives you a file browser, a preview pane,
and a toolbar of analysis tools in one window, so that a curation pass —
inspect, align, clean, extend, rebuild the consensus, check again — happens
without leaving the application or assembling command lines by hand.

This manual assumes you are familiar with TE biology and with tools such as
BLAST, MAFFT and RepeatMasker. It does not assume you have used RepeatViewer
before.

---

## Contents

1. [About RepeatViewer](#1-about-repeatviewer)
2. [Installation](#2-installation)
3. [Quick start: your first curation pass](#3-quick-start-your-first-curation-pass)
4. [The workspace](#4-the-workspace)
5. [Setting up a project](#5-setting-up-a-project)
6. [How tools behave](#6-how-tools-behave)
7. [Sequence operations](#7-sequence-operations)
8. [Analysis tools](#8-analysis-tools)
9. [Visualization](#9-visualization)
10. [Benchmarking](#10-benchmarking)
11. [Workflows](#11-workflows)
12. [Plugins](#12-plugins)
13. [Preferences](#13-preferences)
14. [Command line](#14-command-line)
15. [Output file reference](#15-output-file-reference)
16. [Troubleshooting](#16-troubleshooting)
17. [Citation](#17-citation)

---

## 1. About RepeatViewer

### What it is for

Automated TE discovery pipelines produce candidate consensus sequences that are
rarely publication-ready. Boundaries are wrong, terminal repeats are truncated,
alignments contain junk rows, and classifications need checking. That last mile
is manual, iterative, and involves many small tool invocations on the same
file.

RepeatViewer is built for that loop. You point it at a folder of candidates,
click a file to see it, run a tool from the toolbar, and the output appears in
the same folder — highlighted, and ready to be the input to the next step.

### What it is not

RepeatViewer does not discover TEs and does not build consensus sequences from
a raw genome. That is the job of [TEtrimmer](https://github.com/qjiangzhao/TEtrimmer),
of which RepeatViewer is the interactive companion. The two are independent
applications; RepeatViewer works equally well on any organised folder of TE
sequences, and does not require a TEtrimmer installation.

### Naming

The window title and the application are called **RepeatViewer**. The Python
package and the command you type in a terminal are **`repeatviewer`**. Older
material may refer to the project as TEviewer.

---

## 2. Installation

Choose the packaged application if you just want to use RepeatViewer, or the
source installation if you want to modify it or write plugins.

### macOS

1. Download the `.dmg`.
2. Open it and drag **RepeatViewer** to Applications.
3. On first launch, if macOS refuses to open an unsigned build, **right-click
   the app → Open**, then confirm. This is only needed once.

### Windows

1. Download and unzip the release archive.
2. Run **RepeatViewer.exe**.

Bundled BLAST, MAFFT and MUSCLE run natively on Windows. Linux-only tools
(RepeatMasker, FastGA) require WSL2 — see **External tools** below.

### Linux

1. Unpack the `.tar.gz`.
2. Run `./RepeatViewer`.

### From source

Requires Python 3.10 or newer. A conda environment or venv is strongly
recommended.

```bash
git clone https://github.com/qjiangzhao/RepeatViewer.git
cd RepeatViewer
conda create -n repeatviewer python=3.10 -y
conda activate repeatviewer
pip install -e .
repeatviewer
```

`pip install -e .` also installs the RepeatViewer build of Flexidot, which
produces the interactive HTML dotplots used by the preview pane. The stock
PyPI Flexidot only writes static images and will not give you
**Go-to → Annotation**.

### What is bundled, and what is not

BLAST+, MAFFT and MUSCLE binaries ship inside the application for macOS, Linux
and Windows. You do not need to install them.

These are optional and are detected on `PATH` at runtime:

| Program | Needed for |
|---------|------------|
| Reference **genome FASTA** | BLAST, MSA Extend, genome maps, genome-aware TE-Aid |
| **RepeatMasker** | Repeat Quantify (RepeatMasker engine), TE library benchmarking |
| **FastGA** | The FastGA alignment tool |
| **bedtools** | Some extension and extraction paths |
| **IQ-TREE / FastTree** | Maximum-likelihood trees. Without them, RepeatViewer falls back to Biopython distance trees (neighbour-joining / UPGMA), which need no installation |
| **AliView** | **Tools → Open in AliView** |

### External tools

If you keep the optional programs in a conda environment, tell RepeatViewer
where it is: **File → Setting → External tools**, then enter the environment
name (for example `tetrimmer`) or browse to the environment folder.

- **macOS and Linux** — RepeatViewer prepends that environment's `bin` to
  `PATH`. This matters most for the packaged `.app`, which does not otherwise
  inherit your shell's conda setup.
- **Windows** — the same setting points at a conda environment **inside WSL2**.
  Linux CLIs are launched through `wsl`; the bundled BLAST, MAFFT and MUSCLE
  stay native.

Use the **check** button next to the setting to confirm the environment
resolves and to see which tools are visible.

---

## 3. Quick start: your first curation pass

This walks through one complete pass on a single candidate. Fifteen minutes.

### Step 1 — Add your folder

**File → Add Folder…** (⇧⌘O / Ctrl+Shift+O), or the **+** button on the
Folders panel header, or drag a folder from your file manager onto the Folders
panel.

Pick the folder holding your candidate FASTA files — a TEtrimmer output
directory, or any folder you have organised yourself. It becomes a *workspace
root* and is remembered between sessions.

### Step 2 — Manage Database

**File → Manage Database…** (⇧⌘G / Ctrl+Shift+G).

Three tabs:

- **BLAST database** — FASTA DNA / FASTA Protein (indexes under
  `Documents/RepeatViewer/blast_databases`). FASTA DNA can be set as the main
  reference genome; a `.chrom.sizes` sidecar is written for every FASTA DNA
  BLAST database.
- **HMM database** — HMM DNA / HMM Protein for pyhmmer / HMMER
- **Others** — CDD (RPS-BLAST)

Each entry has a **Name** (defaults to the file name) and an optional
**Annotation name** (max 15 characters; shown as `-` when empty). Check
**Annotation** to show it on the Annotation tab. Click column headers to sort.
**Add…** from any tab lists all types; the chosen Type decides which tab stores
the entry.

### Step 3 — Look at a candidate

Click a FASTA file in the Files panel. The Preview pane opens it with three
sub-tabs:

- **Text** — the sequence itself.
- **Annotation** — a feature map: ORFs, terminal repeats, and any annotations
  you have added.
- **Dotplot** — a self dotplot, which is the fastest way to see internal
  tandem structure or an unresolved terminal repeat.

### Step 4 — Find the copies

With the file selected, click **BLAST** on the toolbar. RepeatViewer searches
your consensus against the reference genome and writes two files beside the
source: `<yourfile>_b.txt`, the hit table you can read, and
`<yourfile>_b.txt.fa`, the hit sequences pulled out of the genome. The
second one is what you align in the next step.

New outputs are highlighted **green** in the Files panel; the inputs that
produced them are highlighted **blue**. You do not have to hunt for what just
appeared.

### Step 5 — Build and clean an alignment

Click **Align** to align the BLAST hits into a multiple sequence alignment
(`hits.fa_mafft.fa` or `hits.fa_muscle.fa`; `_1`, `_2`, … before `.fa` if that
name already exists).

Real alignments of genomic copies contain rows that should not be there —
fragments, unrelated insertions, mostly-gap sequences. Click **MSA Clean** to
remove them. The combined clean runs three steps: remove gappy columns, trim
rows with excessive divergence, and trim rows that are mostly gaps.

### Step 6 — Fix the boundaries

If the alignment shows your consensus is cut short — the alignment block runs
to the edge and keeps going — click **MSA Extend** to pull in flanking genomic
sequence on either side, then clean and re-inspect.

This is the step that most often needs repeating. Extend, clean, look, extend
again.

### Step 7 — Rebuild the consensus

Click **Consensus** to collapse the cleaned alignment back into a single
sequence, using a majority threshold you control.

### Step 8 — Check your work

Run **TE-Aid** for an interactive HTML report of the new consensus — coverage,
terminal repeats, ORFs, genome hits — and **Flexidot** for a dotplot comparing
the old and new versions.

If the boundaries are now clean and the terminal repeats are complete, that
candidate is done. If not, go back to step 6.

---

## 4. The workspace

RepeatViewer's main window is three docked panels, a toolbar, a menu bar and a
status bar.

### Folders panel

Your workspace roots and their subfolders.

| Control | What it does |
|---------|--------------|
| **+** | Add a folder as a new workspace root |
| **◀ ▶** | Back and forward through folders you have visited |
| **↑** | Go to the parent folder |
| **Filter by name…** | Show only folders matching the text |

Each folder shows a count in the form `(7/76)` — files matching the current
filter, out of files in the folder.

Right-click a folder for: **Add Folder…**, **Create Folder…**, **Rename
folder**, **Paste**, colour assignment, and removing a root from the list
(which does not delete it from disk).

### Files panel

The contents of the selected folder, as a sortable table. Available columns:
name, length, number of sequences, GC%, type, description, modified date and
size. Not all are shown by default — use **Columns** to choose.

Clicking a file previews it. Checkboxes and ⌘/Ctrl-click select several files
for tools that accept more than one input.

Right-click a file for preview options, rename, duplicate, cut/copy/paste,
**New Folder with Selection**, **Reveal in Finder**, delete, and — on FASTA
files — sequence-specific actions such as separating a multi-FASTA into
individual files or adding the file as a BLAST database.

### Preview panel

Shows whatever you have selected:

- **FASTA** — Text, Annotation and Dotplot sub-tabs
- **HTML reports** — TE-Aid, Flexidot, genome maps, quantification reports,
  rendered live and interactive
- **Alignments**, **PDFs**, **images**, **plain text**

Controls on the header: zoom percentage, **Actions**, and **Open in browser**
for HTML reports you want in a full-size window.

Preview opens in the active tab by default. **Open in Additional Tab** from
the right-click menu keeps the current one and adds another, which is how you
compare a before and after.

### Status bar

Bottom of the window. Shows either the name of the file you last previewed —
as a link that reveals it in your file manager — or the result of the last
job. While a job runs it shows progress and a **Stop** button.

### Toolbar

The tools you use most, as buttons. The default set is BLAST, Align, Flexidot,
TE-Aid, MSA Clean, MSA Extend, Consensus and Visualization.

Change which buttons appear and in what order under **File → Setting →
Toolbar**. Every built-in tool can be shown, hidden or reordered, and
separators can be inserted.

Several toolbar buttons are **two-button**: left-click runs with the current
settings, right-click opens the parameter dialog. For MSA Clean and TE
annotation you can configure what right-click does.

### Layout

Folders, Files and Preview can be resized, rearranged, or detached into their
own windows — useful on a second monitor. **View** toggles each panel;
**View → Reset layout** restores the default arrangement.

---

## 5. Setting up a project

### Manage Database

**File → Manage Database…** (⇧⌘G / Ctrl+Shift+G)

One window with three tabs:

- **BLAST database** — FASTA DNA / FASTA Protein indexes under
  `Documents/RepeatViewer/blast_databases` (Name defaults to the FASTA file
  name). FASTA DNA entries can be **Set as main reference**; a chromosome
  sizes file (`.chrom.sizes`) is generated for every FASTA DNA BLAST database.
- **HMM database** — HMM DNA / HMM Protein profiles
- **Others** — CDD directories for RPS-BLAST

**Columns:** Name, Annotation name (max 15; `-` when empty), Type,
Path / source (middle-elided like the Files panel), Time (YYYY-MM-DD),
Annotation checkbox. Entries that are not ready are hidden (no Status column).
Click a column header to sort. The current main reference is shown below the
table (not as a column).

**Types**

| Type | Tab | Used by |
|------|-----|---------|
| FASTA DNA | BLAST | blastn; Set as main reference; optional Annotation |
| FASTA Protein | BLAST | blastp, blastx |
| HMM DNA / HMM Protein | HMM | HMMER / pyhmmer |
| CDD | Others | RPS-BLAST (Annotation tab) |

- **Add…** — all types are available from every tab; Type decides where the
  entry is stored. Name and Annotation name have `?` help. Path field matches
  other fields in width. For FASTA DNA optionally Set as main reference.
- **Remove** — delete selected entries.
- **Set as main reference** — marks a FASTA DNA BLAST (or legacy genome) as
  the current genome.
- **Annotation** column — toggle which databases appear on the Annotation tab
  (fills Annotation name from Name if empty, capped at 15 characters).
  FASTA DNA reference genomes can be enabled for annotation the same way.

The same window opens from BLAST / HMMER / Preferences **Manage** buttons.

---

## 6. How tools behave

Five rules hold across every tool. Learning them once saves reading each
section closely.

**Outputs land beside their input.** A tool run on
`/data/run2/CroSat2.fasta` writes into `/data/run2/`. RepeatViewer never
collects results into a separate output directory, because the next tool in
your curation loop takes the previous output as input and you want them
together.

**New files are highlighted.** Outputs turn green in the Files panel, and the
inputs that produced them turn blue. The highlight clears when you navigate
away.

**Jobs run in the background.** The interface stays usable — you can browse
folders and preview files while BLAST runs. Progress and a **Stop** button
appear in the status bar. Only one job runs at a time; starting a second
before the first finishes is refused.

**Settings persist.** Parameters you set in a tool dialog are remembered for
the next run and between sessions. Each dialog's **Options** tab has a reset.

**Selection drives everything.** Tools act on what is selected in the Files
panel. Most accept several files; some require exactly one. If nothing
suitable is selected the tool says so instead of guessing.

---

## 7. Sequence operations

From the **Sequence** menu, and from the right-click menu in the Files panel.

### New Sequence

Create a FASTA from pasted text, with name, description, organism and type
(nucleotide or protein).

### Combine Sequences to One File

Select several FASTA files and merge them into one multi-FASTA — the usual way
to assemble a library from individual curated consensus files.

### Extract Sequences from File

Pick individual records out of a multi-FASTA. Choose **Write all selected into
one file** to get a single output, or leave it unchecked for one file per
record.

### Separate sequences

On a multi-FASTA, right-click offers **Separate sequence into N files** or
**into N folders** — the reverse of Combine, for working through a library
record by record.

### Reverse Complement

Writes `<name>_rc.fa`. Use it when a candidate was built on the wrong strand,
which is common and worth checking before you curate boundaries.

### ORF Prediction

Finds open reading frames and writes `<name>_orfs.fa`.

| Parameter | Meaning |
|-----------|---------|
| Min ORF length | Shortest ORF reported |
| Max ORF length | Longest ORF reported |
| Add other start codons | Include non-ATG starts |

Predicted ORFs also appear in the **Annotation** preview tab, which is usually
where you want to look at them.

### Translate

Writes the protein translation as `<name>_translated.faa`.

---

## 8. Analysis tools

### BLAST

**Toolbar: BLAST · Tools → Blast**

Searches your sequences against the reference genome or a saved database.
Three programs, each on its own tab: **BLASTn** (nucleotide), **BLASTp** and
**BLASTx** (protein). Each program tab uses **Inputs** / **Settings** /
collapsible **Advanced Settings** (same pattern as Align). There is no Mode —
each selected query file is blasted separately. The **Options** tab sets the
right-click program, system BLAST+, and the process window (shown while the
tool runs and closed automatically on success by default).

| Parameter | Meaning |
|-----------|---------|
| Inputs | Use selected Files-panel FASTA, or Define inputs for this run only (paths not remembered) |
| Database | Reference genome, or a saved BLAST database |
| `-evalue` | Expectation threshold (NCBI BLAST+) |
| Threads | `-num_threads`; defaults to your CPU count minus the UI reserve |
| `-word_size` | Seed length. Smaller is more sensitive and much slower on repetitive genomes |
| Export text / GFF | Keep `*_b.txt` and/or write `*_b.gff` |
| Max sequence number | Cap on hit sequences kept for the extracted FASTA / MSA (default 100) |
| Top longest sequences | When capping: keep this many longest hits first, then fill at random up to Max sequence number (default 100) |
| `-gapopen`, `-gapextend` | Advanced: gap penalties |
| `-penalty` | Advanced: mismatch penalty (BLASTn) |
| `-perc_identity` | Advanced: minimum percent identity (BLASTn) |
| `-qcov_hsp_perc` | Advanced: minimum query coverage per HSP |

**Outputs:**

| File | Contents |
|------|----------|
| `<name>_b.txt` | Tabular hit table |
| `<name>_b.txt.fa` | Hit sequences extracted from the genome — the input to Align |
| `<name>_b.gff` | Hit positions, when GFF output is enabled |

The hit table is ordinary tabular BLAST output — readable in the preview pane,
and the input for the TE density and divergence plots. **Max sequence number**
and **Top longest sequences** cap how many hits are carried into the `.fa`, which
matters because a high-copy element can otherwise produce an alignment too
large to work with. BLAST+ option meanings follow the
[NCBI BLAST Command Line Applications User Manual](https://www.ncbi.nlm.nih.gov/books/NBK1762/).

### Align

**Toolbar: Align · Tools → Align** (MAFFT, MUSCLE, FastGA, Align all
selected files in one MSA, Change alignment parameters)

Builds a multiple sequence alignment with **MAFFT** or **MUSCLE**, both
bundled, or a whole-genome alignment with **FastGA** (see [FastGA](#fastga)).
The window has the tabs MAFFT, MUSCLE, FastGA and Options; **Run** runs the
tab you are on. The **Options** tab sets what right-clicking the Align button
runs (MAFFT, MUSCLE or FastGA) and what right-clicking the FastGA button does.

MAFFT parameters:

| Parameter | Meaning |
|-----------|---------|
| Mode | Combine selected sequences and align, or align each selected file separately |
| Inputs | Use selected Files-panel FASTA, or Define inputs for this run only (paths not remembered) |
| Strategy | Algorithm preset (FFT-NS-1/2, L/G/E-INS-i, Auto, Custom) |
| Sequence type | Auto-detect, nucleotide, or protein |
| Threads | Parallelism (`--thread`) |
| `--retree` / `--maxiterate` | Advanced: guide-tree rebuilds and refinement cycles (Custom) |
| `--op` / `--ep` | Advanced: gap-opening penalty and offset |
| `--reorder` / `--adjustdirectionaccurately` | Advanced: output order and strand correction |

Strategy meanings follow the [MAFFT algorithms](https://mafft.cbrc.jp/alignment/software/algorithms/algorithms.html) page: progressive FFT-NS methods are faster; L/G/E-INS-i are iterative and more accurate.

MUSCLE parameters:

| Parameter | Meaning |
|-----------|---------|
| Mode | Combine selected sequences and align, or align each selected file separately |
| Inputs | Use selected Files-panel FASTA, or Define inputs for this run only (paths not remembered) |
| Algorithm | `-align` (high accuracy) or `-super5` (large sets) |
| Threads | Parallelism (`-threads`) |
| Guide tree permutation | Advanced: `-perm` (`none`, `abc`, `acb`, `bca`) |
| Perturbation seed | Advanced: `-perturb` (`0` = off) |

Ensemble options (`-stratified`, `-diversified`, `-resample`, …) are not
exposed here because Align writes a single MSA for the viewer.

**Output:** one file at a time → `<name>_mafft.fa` / `<name>_muscle.fa`;
combined sequences → `<name>_cmafft.fa` / `<name>_cmuscle.fa` (then
`_1`, `_2`, … before `.fa` if that name already exists).

### MSA Clean

**Toolbar: MSA Clean · Tools → MSA clean**

Removes the rows and columns that make a genomic-copy alignment unusable.
Three steps, individually or combined:

| Step | Removes |
|------|---------|
| **Remove gappy columns** | Alignment columns that are mostly gaps |
| **Clean rows with divergence** | Sequences too divergent from the alignment consensus |
| **Clean rows with gaps** | Sequences that are mostly gaps |

| Parameter | Applies to |
|-----------|-----------|
| Divergence threshold | Row divergence cleaning |
| Window size (bp) | Sliding window for divergence |
| Gap proportion threshold | Row gap cleaning |
| Gap threshold | Column removal |
| Similarity check gap threshold | Similarity comparison |
| Similarity threshold | Similarity comparison |
| Min nucleotides per column | Column retention floor |

**Combined MSA clean** runs the checked steps in the order you set — drag to
reorder. **Right-click action** sets what a right-click on the toolbar button
does, so your most common variant is one click away.

**Outputs:** `<name>_gap.fa` (columns), `<name>_rowdiv.fa` (divergent rows),
`<name>_rowgap.fa` (gappy rows).

Each step writes its own file, so you can see what each removed and back out
one step without redoing the others.

### MSA Extend

**Toolbar: MSA Extend · Tools → MSA Extend**

Extends an alignment outward by pulling flanking sequence from the genome.
This is how you recover truncated boundaries and complete terminal repeats.

| Parameter | Meaning |
|-----------|---------|
| Genome / Pan-genome | Source of flanking sequence |
| Extend left of MSA | Bases to add on the 5′ side |
| Extend right of MSA | Bases to add on the 3′ side |

Requires a defined genome. The usual rhythm is extend, clean, inspect the new
edges, and extend again until the alignment block ends in genomic sequence
that no longer aligns.

### Consensus

**Toolbar: Consensus · Tools → Consensus**

Collapses an alignment to a single sequence.

| Parameter | Meaning |
|-----------|---------|
| Sequence type | Auto-detect, nucleotide, or protein (Biopython alignment consensus) |
| Majority threshold | Fraction of rows that must agree for a residue to be called |
| Ambiguous base / residue | **N** (DNA/RNA) or **X** (protein) where no residue reaches the threshold |

Consensus needs either a **multiple sequence alignment** — every row the same
length — or a single-sequence FASTA. An unaligned multi-FASTA is rejected with
a message telling you to align it first.

**Output:** a consensus FASTA whose name keeps the full input filename and
records the threshold, e.g. `CroSat2.fasta_0.8_co.fa`. Several thresholds can
therefore coexist in one folder and be compared.

### TE-Aid

**Toolbar: TE-Aid · Tools → TE-Aid**

An interactive HTML report for one consensus: self-alignment structure,
terminal repeats, ORFs, and genome hit coverage in one page. This is the
clearest single view of whether a consensus is finished. The parameters
dialog follows the BLAST / Flexidot layout: **Inputs**, **Settings**,
collapsible **Advanced Settings**, and **Options**.

| Section | Parameter | Meaning |
|---------|-----------|---------|
| Inputs | Use selected / Define inputs | Files-panel selection, or explicit FASTA paths for this run only |
| Settings | Genome | Project genome or a saved nucleotide BLAST database |
| Settings | Genome BLAST E-value | Threshold for genome BLAST (finds copies of the input) |
| Settings | Self BLAST word size | Seed length for self-BLAST (finds repetitive regions) |
| Settings | Threads | Threads for genome BLAST and self-BLAST |
| Advanced | Line transparency | Opacity of genome BLAST hit lines |
| Advanced | Full-/mid-length hits | Color bands by hit length vs consensus |
| Options | Process output | Live genome BLAST and self-BLAST progress window |

**Output:** `<name>_TEAid.html`, rendered in the preview pane.

### Flexidot

**Toolbar: Flexidot · Tools → Flexidot**

Dotplots: self, pairwise, or all-against-all for a multi-FASTA. The parameters
dialog follows the same layout as Align: **Mode**, **Inputs**, **Outputs**,
**Settings**, and collapsible **Advanced Settings**.

| Section | Parameter | Meaning |
|---------|-----------|---------|
| Mode | Combine / separate | One dotplot for all selected files, or one run per file |
| Inputs | Use selected / Define inputs | Files-panel selection, or explicit FASTA paths for this run only |
| Outputs | Output format (`-f`) | HTML for interactive, or a static image |
| Settings | k-mer (`-k`) | Match length. Lower finds more, and more noise |
| Settings | Sequence type (`-t`) | DNA (nuc) or protein (aa) |
| Settings | Show GFF features | Shade with companion annotation GFF (`-g`) |
| Advanced | Allowed substitutions (`-S`) | Mismatches tolerated within a word |
| Advanced | Maximum N content (`--max_n`) | Skip windows above this N fraction when wobble is on |
| Advanced | Match line thickness (`-A`) | Line width in the plot |

Choose HTML output to get **Go-to → Annotation**: clicking a region in the
dotplot jumps the Annotation preview tab to that coordinate. This is the
feature that makes a dotplot actionable rather than decorative, and it needs
the RepeatViewer build of Flexidot.

**Outputs:** `<name>_flexidot.html`, `<name>_flexidot.gff`,
`<name>_flexidot_legend.png`.

### TSD search

**Tools → TSD search**

Finds target site duplications flanking an element in an alignment — direct
repeats immediately outside the element boundary, and strong evidence that
your boundary is correct.

| Parameter | Meaning |
|-----------|---------|
| Left boundary (MSA column) | Where the element starts |
| Right boundary (MSA column) | Where the element ends |
| Search radius | How far outside the boundary to look |
| Max TSD similarity | Similarity ceiling for candidate pairs |
| Poly-N threshold | Reject candidates above this N content |
| Cluster distance (eps) | Clustering tolerance when grouping candidates |

Set the boundaries from what you see in the Annotation or alignment view, then
widen the radius if nothing is found.

### TE annotation

**Toolbar: TE annotation · Tools → TE annotation**

Classifies elements using the databases configured under **File → Define
database**. Two annotation sources:

- **Auto annotation** — structural and domain-based, including CDD protein
  domains
- **Similarity annotation** — comparison against a reference library

**Element colors** sets how annotation classes are drawn in the Annotation
preview tab. **Right-click action** chooses what right-clicking the toolbar
button runs.

If the dialog reports no databases, set them up under **File → Define
database** first.

### Repeat Quantify

**Toolbar: Repeat Quantify · Tools → Repeat Quantify**

Measures how much of a genome a library of monomers or satellites accounts
for. Aimed at tandem repeats and satellites, where "how many copies and how
much of the genome" is the question.

**Inputs:**

- **Library** — the monomer or satellite FASTA
- **Genome** — the assembly to quantify against. **Add another genome** runs
  the same library against several genomes, producing one report each

The dialog is non-modal: the Folders and Files panels stay live while it is
open, and clicking a FASTA there fills whichever input field last had focus.

**Engines**, chosen by tab:

- **BLASTn** — fast, uses the bundled BLAST
- **RepeatMasker** — requires RepeatMasker on `PATH`; also produces a
  divergence landscape
- **Both** — runs each and reports them side by side

| BLASTn parameter | Meaning |
|------------------|---------|
| Word size | Seed length; 11 is a good balance, 7 can hang on large genomes |
| E-value | Expectation threshold |
| Min identity | `-perc_identity`; 0 means no minimum |
| Min coverage | `-qcov_hsp_perc` per hit; 0 means no minimum |

| RepeatMasker parameter | Meaning |
|------------------------|---------|
| Divergence cutoff (`-div`) | Maximum divergence to count |
| Sensitivity | RepeatMasker speed/sensitivity preset |

**Output:** a folder `<genome>_repeat_quant/` beside the genome, containing an
HTML report plus TSV and GFF tables. The report gives, per family:

| Column | Meaning |
|--------|---------|
| family | Library family matched |
| monomer_len | Length of the library monomer (period size), bp |
| copy_number | Estimated copies = covered_bp ÷ monomer length; shown as a fraction when a family covers less than one full monomer |
| hits | Raw alignment fragment count |
| arrays | Hits merged into contiguous blocks (≈ distinct array loci) |
| covered_bp | Target bp covered, overlaps counted once |
| %target | covered_bp as a percentage of target length |
| mean_array_bp | Average tandem array length = covered_bp ÷ arrays |
| max_array_bp | Largest single contiguous array, bp |
| mean_identity | Length-weighted mean % identity to the consensus (high = young or conserved) |
| mean_divergence | RepeatMasker Kimura substitution % (high = old or diverged) |
| scaffolds | Number of distinct target sequences hit |

If a family name appears more than once in your library, the report warns you:
hits are merged into one row and copy number uses a single monomer length.
Give each entry a unique name to keep them separate.

The genome BLAST index is cached between runs, so quantifying the same genome
repeatedly does not re-index it each time. **File → Setting → General →
Diagnostics → Clear cache** removes them, and shows how much space they
currently use.

### FastGA

**Toolbar: FastGA · Tools → Align → FastGA · Tools → FastGA**

Whole-genome alignment for large sequences, where BLAST is too slow. FastGA
lives in the **Align** window, on its own **FastGA** tab; the FastGA button and
**Tools → FastGA → Change FastGA parameters** open that tab. Requires FastGA
1.5 or newer from conda into a shared env, for example:

```bash
# if rv_env exists
conda install -n rv_env -c bioconda fastga=1.5.1
# else
conda create -n rv_env -c bioconda fastga=1.5.1
```

Then point **File → Setting → External tools** at `rv_env`. The Options
tab shows where it was found. On Windows it runs in WSL.

**Mode**

| Mode | Runs | Output (next to the input) |
|------|------|----------------------------|
| Each file against itself | one run per selected FASTA | `<file>_FastGA.paf` |
| Selected files against each other | every pair (A, B, C → A–B, A–C, B–C) | `<first>_multi_FastGA.paf` |
| Selected files against a reference | one run per selected FASTA | `<file>_FastGA_vs_<reference>.paf`, or with **Write all into one file** `<first>_multi_FastGA_vs_<reference>.paf` |

Against a reference: select the files to compare, open the window, then click
the reference in the Files panel (or Browse). The files compared are the ones
selected when the window opened; the reference is skipped if it is among them.
The reference is indexed once (GIXmake) and reused for every run. With **Keep
genome index** the index is kept in a hidden `.repeatviewer_fastga` folder
next to the reference and reused while the reference is unchanged; an index
takes about 14 GB per gigabase. The grey line under the reference says what
Run will do ("3 FASTA selected → 3 runs against genome.fa").

| Parameter | Meaning |
|-----------|---------|
| Output format | PAF (`-paf`, `-pafx`, `-pafm`, `-pafs`, `-pafS`), PSL, or binary `.1aln`. PAF and `.1aln` open in the dotplot view |
| `-T` threads | Parallelism |
| `-i` min identity | Identity floor |
| `-l` min length | Shortest alignment reported (≥ 100) |
| `-f` seed cutoff | Seeds more frequent than this are dropped as repetitive; RepeatViewer uses 200 (FastGA: 10) so repeats stay visible |
| `-c` / `-s` | Chain coverage / chain gap |
| Symmetric seeding (`-S`) | Result independent of argument order; for repeat structure |
| Use soft mask | Seed only in uppercase bases. Sent as `-M` (FastGA 1.5+) or by leaving out `-N` (older FastGA), whichever the installed version takes |
| Keep genome index (`-k`) | Keep `.1gdb` / `.gix` for later runs |
| Temp dir (`-P`) | Folder for FastGA's large temporary files; empty = hidden folder next to the input |

### Tree

**Toolbar: Tree (add it in File → Setting → Toolbar) · Tools → Tree**

Builds a phylogenetic tree from the selected FASTA file(s). Several selected
files are combined into one tree. RepeatViewer prepares the input, checks it
and writes the report; established programs do the computing:

- **MAFFT** (bundled) aligns unaligned input.
- **IQ-TREE** (maximum likelihood), **FastME** (minimum evolution) or
  **FastTree** (approximate maximum likelihood, for more than ~1,000
  sequences) builds the tree.

The tree programs are not bundled. Install them into the conda environment
set in **File → Setting → External tools**:

```bash
conda install -c bioconda iqtree fastme fasttree
```

If the chosen program is missing, the run stops with this hint. It never
falls back silently to a weaker method.

Left-click opens the dialog; right-click builds with saved parameters (or
opens the dialog, chosen in the Options tab). Without a FASTA selection the
button only asks you to select one; **Tools → Tree → Change Tree parameters**
opens the dialog anyway, to adjust settings in advance. **Tools → Tree → Tree
with default parameters** runs with the saved settings, like TE-Aid. The round
**?** buttons next to the MAFFT strategy, the method and the substitution model
explain the choices; options that do not apply to the current choices are
greyed out.

| Tab | Setting | Default | Why |
|-----|---------|---------|-----|
| Alignment | Alignment | Align if not already aligned | Equal-length input is used as an alignment |
| | MAFFT strategy | L-INS-i | Accurate up to ~200 sequences; FFT-NS-2 for larger sets |
| | Adjust direction | On | TE/satellite copies from the other strand are reverse-complemented instead of looking extremely divergent |
| | Minimum length | 70 % of the median length | Fragments (5'-truncated copies) are excluded, listed in the report and saved to `.excluded.fa`. The grey line below says how many. Use *of the longest sequence* when many sequences may be fragments, or all should be about the same length (monomers) |
| | More options: `--op` / `--ep` | 1.53 / 0 | MAFFT defaults; change only with a reason |
| Tree | Method | IQ-TREE | |
| | Substitution model | ModelFinder, best by BIC | The chosen model goes into the report |
| | Independent searches | 5 | Separate IQ-TREE searches (seeds seed, seed+1, …); the best is kept and the report says how many found the same topology. ModelFinder runs once, in the first search |
| | Random seed | 12345 | Recorded, so the run can be repeated exactly |
| | Threads | Preferences value | Maximum; IQ-TREE chooses the best number up to it (`-T AUTO`) |
| Support and rooting | SH-aLRT / UFBoot2 | 1000 / 1000, `-bnni` | Shown as SH-aLRT/UFBoot; commonly trusted at ≥ 80 / ≥ 95 |
| | Bootstrap (FastME) | 100 | Trusted at ≥ 70 |
| | Root | None | Outgroup (tick one or more sequences) or midpoint |
| Options | Right-click action, process window, MAFFT / tree programs from `PATH` | | |

**Output** — like every tool, the results are written next to the input and
keep its full name. When the run finishes, the input is highlighted light blue,
the new files green, and the report opens (deep blue):

| File | Content |
|------|---------|
| `<input>_tree_report.html` | Tree, run report, warnings, methods text, citations (opened automatically) |
| `<input>_tree.treefile` | Newick with support values, for FigTree / iTOL / ggtree |
| `<input>_tree.svg` / `.pdf` / `.png` | Tree figure |
| `<input>_tree.aln.fa` | The alignment the tree was built from |

Supporting files go into `<input>_tree_details/`: the IQ-TREE report and log
(`.iqtree` / `.log`, plus `.modelfinder.*` when the best search was not the
first), the methods paragraph with references (`_methods.txt`), program
versions and exact command lines (`_commands.txt`), tree names ↔ original
FASTA headers (`.names.tsv`) and the excluded fragments (`.excluded.fa`).
A second run on the same input writes `<input>_tree_2…`; nothing is
overwritten.

Rooting keeps each support value on its branch (a common mistake when
rerooting Newick trees). If several outgroup sequences do not form one clade,
the report warns that the root position is uncertain.

**If a run fails**, no half-finished results are left. The message gives the
reason in plain words, and the full log (commands and program output) is saved
next to the input as `<input>_tree_failed.log`, ready to send with a bug report.

### Mapping

**Toolbar: Mapping (add it in File → Setting → Toolbar) · Tools → Mapping**

Places the selected sequences or reads on a reference and shows what differs,
like Geneious' *Map to Reference*. The dialog follows Align: **Minimap2**,
**HISAT2**, and **Options** tabs. Select the files to map in the Files panel
(FASTA or FASTQ, gzipped too), or use **Define inputs**; click the reference in
the Files panel (or Browse).

Two mappers; **Auto** (Options → mapper for right-click) picks one from the input:

- **minimap2** — long reads (ONT, HiFi) or assembled contigs against a genome.
- **HISAT2** — short DNA / RNA-seq (Illumina), including spliced RNA.

minimap2 **2.31** is bundled with RepeatViewer where platform builds exist.
HISAT2 is **not** bundled (binaries are large); install **2.2 or newer**
(recommended **2.2.3**) via conda or PATH:

```
conda install -c conda-forge -c bioconda hisat2=2.2.3
```

To use your own minimap2 instead of the bundle, tick
**Use Minimap2 from the system PATH or defined Conda environment** in
Mapping → Options.

A mapper you chose that is missing stops the run with this hint. If Auto's
choice is missing, the other read mapper is used and the report says so.
HISAT2 builds its index once and keeps it next to the reference in the hidden
folder `.repeatviewer_hisat2`.

| Tab | Setting | Default | Why |
|-----|---------|---------|-----|
| Minimap2 | Inputs / Reference | Use selected; project genome if set | Or Define inputs / Browse |
| | Map mode (`-x`) | map-ont | Real minimap2 preset (no Auto); includes RNA/cDNA splice presets |
| | Secondaries | Primary only | `--secondary=no` / `-N 50` |
| | Threads | CPU count | minimap2 `-t` / HISAT2 `-p` |
| | Output Format | PAF | Or SAM |
| | Advanced `-k` / `-w` | From preset | e.g. map-ont 15/10, splice 15/5, sr 21/11 |
| | Advanced CIGAR / CS | On for PAF | `-c` and `--cs=long` (contig view needs CS) |
| HISAT2 | Library | RNA-seq (spliced) | Or DNA-seq (`--no-spliced-alignment`) |
| | Reads | Single-end | Or Paired-end Illumina R1/R2 (`-1`/`-2` from filenames) |
| | Strandness | Unstranded | `--rna-strandness` F / R (SE) or FR / RF (PE); not file pairing |
| | Report | Best place only | Or `-k 20` |
| | Threads | CPU count | |
| | Advanced `--dta` | Off | Longer anchors for StringTie (RNA) |
| | Advanced intron lengths | 20 / 500000 | `--min-intronlen` / `--max-intronlen` |
| | Advanced `--no-softclip` | Off | Disallow soft-clipping |
| Options | Right-click mapper | Auto | Auto, minimap2, or HISAT2 |
| | Process window, executables | | Use Minimap2 from the system PATH or defined Conda environment |

**Map mode (`-x`).** Choose a real minimap2 preset (there is no Auto in minimap2).
Default is `map-ont`. Changing Map mode fills `-k` / `-w`.

- **Long genomic:** `map-ont` (ONT), `map-pb` (PacBio CLR), `map-hifi` (HiFi/CCS),
  `lr:hq` (Nanopore Q20), `map-iclr` (Illumina Complete Long Reads)
- **RNA / cDNA:** `splice` (long RNA/cDNA, strand unknown), `splice:hq`
  (PacBio Iso-seq / accurate long RNA), `splice:sr` (short RNA-seq)
- **Assemblies / short DNA:** `asm5` / `asm10` / `asm20` (~0.1 / 1 / 5% divergence),
  `sr` (short genomic reads)

**The contig view.** The result, `<input>_map.rvmap`, is written next to the
input and opens in the contig view: the reference annotations, the consensus
(differences from the reference in red), the coverage and the reference, then
every placed sequence with the bases that differ from the reference coloured.
Deletions are drawn as a line, insertions as a purple mark. Ctrl/Cmd + wheel
zooms, **Bases** zooms until the bases are written out, **Fit** (or a
double-click) shows the whole reference. Four tables below:

- **Sequences** — position, strand, identity, mismatches, deleted and inserted
  bases of every placement.
- **Features** — every reference feature in every placed sequence: *intact*,
  *n diff* (changed), *partial* (only part covered) or *lost* (half or more
  deleted). Hover a cell for its identity and its position on the sequence;
  click it to jump there.
- **Variants** — every difference with its count, depth and frequency, the
  feature it falls in and its cytosine context (CG / CHG), where
  methylation-driven C→T changes are expected.
- **Not used** — the sequences that were not placed, with the reason.

**Results** (top right) opens the report and the other files. Everything
except the `.rvmap` is in `<input>_map_details/`:

| File | Content |
|------|---------|
| `<input>_map_report.html` | Run report, coverage plots, feature and variant summaries, methods text, citations, commands |
| `<input>_map_consensus.fa` | Consensus, with the reference annotations moved onto it (`.RV.gff3`) |
| `<input>_map_mapped.fa` | The placed sequences, turned to the reference strand, with the lifted annotations |
| `<input>_map_unused.fa` | The sequences not placed, with the reason in the header |
| `<input>_map_variants.tsv` / `_features.tsv` / `_reads.tsv` | The three tables |
| `<input>_map.aln.fa` / `.paf` / `.sam` | The mapper's own output |

No pysam or samtools is needed: minimap2 and HISAT2 output is read directly.
A second run on the same input writes `<input>_map_2…`; nothing is
overwritten. If a run fails, nothing half-finished is left and the full log is
saved as `<input>_map_failed.log`.

### Open in AliView

**Tools → Open in AliView**

Opens the selected FASTA in AliView for hands-on alignment editing, when you
need to fix something by eye that no automatic clean will get right. Requires
AliView installed; set the path under **File → Setting → Use system sequence
viewer**.

---

## 9. Visualization

**Toolbar: Visualization · Tools → Visualization**

Three genome-scale plots. Like Repeat Quantify, these dialogs are non-modal —
leave one open and click files in the Files panel to fill its fields.

### TE density plot (genome map)

Where TEs sit along the chromosomes.

| Input | Accepts |
|-------|---------|
| TE masking file | RepeatMasker `.out` or a BLAST hit table |
| TE consensus file | TE consensus FASTA |
| Other annotation file | GFF / GTF / BED / VCF |
| DNA methylation file | Methylation coverage data |

**Add another …** under each field layers several sources on one plot; track
order follows field order, first row at the top.

| Option | Meaning |
|--------|---------|
| Y axis content | What the height shows, e.g. total TE length per window |
| Color by | Grouping, e.g. TE superfamily |
| Separate forward / reverse strands | Split strands into two tracks |
| Collapse overlapped TEs | Merge overlapping annotations |
| Show TE difference track | Compare two TE sources |
| Display content | Which annotation features to draw (gene, mRNA, CDS, exon, transcript …) |
| Methylation coverage min | Coverage floor for methylation points |
| Window size | Binning resolution |
| Chromosome order | Sort order for chromosomes |
| Chromosomes to show | Subset to plot |

**Outputs:** `<name>_{blast|RM}_den_plot.html` (or
`<name>_multi_dens_plot.html` when multiple TE sources are used), and the
same names with `_region` before `.html` for a zoomed region.

### TE divergence plot

The age distribution of a TE family — a histogram of divergence, where a
sharp low-divergence peak means a recent burst and a broad distribution means
sustained older activity.

| Parameter | Meaning |
|-----------|---------|
| TE input files | RepeatMasker `.out` or BLAST hit tables |
| Bin size | Histogram bin width (divergence = 1 − identity, range 0–1) |

**Output:** `<name>_{blast|RM}_div_plot.html` (or
`<name>_multi_div_plot.html` for multiple sources)

### TE abundance plot

Rank-abundance of families in a RepeatMasker run — which families dominate the
genome.

| Parameter | Meaning |
|-----------|---------|
| RepeatMasker .out files | One or more; each becomes a separate curve |
| Y axis | Percentage of genome occupied, or occupied nucleotides (bp) |
| Top families | How many families to show, ranked by occupied bp |

**Output:** `<name>_RM_abund_plot.html` (or
`<name>_multi_abund_plot.html` for multiple sources)

### MSA plot

**Tools → Visualization → MSA plot** (or the MSA plot entry where exposed)

Overview image of an aligned FASTA. The dialog follows Align: **Inputs**,
**Settings**, and collapsible **Advanced Settings** for residue colors.

| Section | Parameter | Meaning |
|---------|-----------|---------|
| Inputs | Use selected / Define inputs | Files-panel selection, or an explicit aligned FASTA |
| Settings | Sequence type | Auto-detect, Nucleotide, or Protein |
| Settings | Output format / Quality | PNG or PDF; DPI for PNG |
| Settings | Sequence height / distance | Row height and gap between rows (cm) |
| Advanced | Nucleotide / amino-acid colors | Presets and per-residue colour pickers (scrollable) |

**Output:** `<name>_msa_plot.png` or `<name>_msa_plot.pdf`

---

## 10. Benchmarking

**Toolbar: Benchmarking · Tools → Benchmarking**

For comparing a curated library against a reference — measuring whether your
curation improved anything.

### Confusion matrix benchmarker

Compares a TEtrimmer output folder against a reference RepeatMasker `.out`,
scoring classification agreement.

| Input | Meaning |
|-------|---------|
| TEtrimmer folder | Folder of curated output to score |
| Reference RM .out | Ground-truth RepeatMasker annotation |

**Outputs:** `<base>_confusion_matrix.txt`,
`<base>_heatmap_proportion_plot.pdf`, `<base>_spider_plot.pdf`,
`<base>_sankey.html`.

### TE library benchmarker

Runs RepeatMasker with your library against a standard library and reports how
they compare.

| Input | Meaning |
|-------|---------|
| Query library FASTA | The library you built |
| Standard library FASTA | The reference to compare against |
| Existing RM .out | Optional: reuse a RepeatMasker run instead of repeating it |
| RepeatMasker threads | Parallelism |

**Output:** `<name>_library_benchmark.html`

Supplying an existing `.out` is worth doing — RepeatMasker on a large genome
is slow, and re-running it to change only the report is wasted time.

---

## 11. Workflows

**Workflows → Design and manage your workflows**

A workflow chains tools so a multi-step curation pass runs as one action. If
you always follow BLAST with Align, MSA Clean and Consensus, make it a
workflow and run it on a whole folder.

### How chaining works

Each tool declares what it accepts and what it produces, in three logical
formats: **fasta**, **html** and **plot**. RepeatViewer only lets you connect
steps whose formats match, so an invalid chain cannot be built.

Tools producing `html` or `plot` are *terminal* — TE-Aid and Flexidot end a
workflow, because nothing consumes their output.

| Tool | Accepts | Produces | Needs genome |
|------|---------|----------|--------------|
| BLAST | fasta | fasta | yes |
| MSA Extend | fasta | fasta | yes |
| Align | fasta | fasta | no |
| MSA Clean (combined, and each individual step) | fasta | fasta | no |
| Consensus | fasta | fasta | no |
| ORF Prediction | fasta | fasta | no |
| Reverse Complement | fasta | fasta | no |
| TE-Aid | fasta | html | yes |
| Flexidot | fasta | plot | yes |

Steps using the genome will not run unless one is defined.

Each step runs with that tool's saved settings — set the parameters in the
tool's own dialog first, then build the workflow.

---

## 12. Plugins

**Plugins menu → Manage plugins / Design your plugin**

Plugins add toolbar tools that use the same file selection and job machinery as
the built-in ones. A plugin receives the selected files and the project
(genome, databases, thread count), and returns the paths it created.

The **Design your plugin** dialog builds a simple command-wrapping plugin
without writing code — enough for wrapping a command-line tool you use often.

For plugins written as Python packages, see **[PLUGINS.md](PLUGINS.md)** for
the API, and **[PLUGINS_TE_BENCHMARK.md](PLUGINS_TE_BENCHMARK.md)** for a
worked example. Sample plugins are in `examples/`.

---

## 13. Preferences

**File → Setting** (⌘, on macOS)

### General

**Appearance** — Theme (Auto follows the operating system), accent colour.

**Sequence tools** — path to an external sequence viewer (AliView);
colour-blind-safe nucleotide colours.

**Files & preview** — font size for the Folders and Files panels; default
preview zoom and fit mode (fit window, fit width, original size); which FASTA
sub-tab opens on a single click (Text, Annotation or Dotplot); whether to show
length and GC% in the preview.

> GC% requires reading every base. On whole-scaffold or whole-genome FASTA
> files that costs a few seconds per file — turn it off if you browse large
> sequences.

**Performance** — threads for tools, and threads kept free for the interface.
Tool thread counts derive from your CPU count minus this reserve, so the
window stays responsive while jobs run.

**Startup** — reopen the last folder and genome; restore window size and
position.

**External tools** — the conda (or WSL2 conda) environment providing
RepeatMasker, FastGA and friends, with a button to check what resolves.

**Diagnostics** — enable DEBUG logging, open the log file, open the config
folder, **Copy diagnostics**, which puts version, OS and folder locations on the
clipboard for a bug report, and **Clear cache**, which shows the current cache
size and deletes it — including the BLAST indexes Repeat Quantify keeps between
runs.

### Toolbar

Choose which tools appear on the toolbar, in what order, with separators.
**Reset to default** restores BLAST, Align, Flexidot, TE-Aid, MSA Clean, MSA
Extend, Consensus and Visualization.

---

## 14. Command line

```bash
repeatviewer [-i FOLDER] [-g GENOME.fa]
```

| Flag | Meaning |
|------|---------|
| `-i`, `--input-dir` | Workspace folder to open on launch |
| `-g`, `--genome` | Reference genome FASTA |
| `--version` | Print the version and exit |

From a source checkout you can also run `python RepeatViewer.py`.

Workspace roots, window layout, genomes and preferences are saved between
sessions, so the flags are for scripted or repeated launches rather than
everyday use.

### Keyboard shortcuts

| Key | Action |
|-----|--------|
| ⇧⌘O / Ctrl+Shift+O | Add Folder |
| ⇧⌘N / Ctrl+Shift+N | New Folder |
| ⇧⌘G / Ctrl+Shift+G | Manage Database |
| ⌘F / Ctrl+F | Focus the file or folder filter, depending on the active panel |
| F5 | Refresh the current folder |
| Delete | Delete selected files, with confirmation |
| ⌘Q / Ctrl+Q | Quit |

---

## 15. Output file reference

All outputs are written beside their input file unless noted.
`<name>` is the **full input filename including its extension** (same pattern
as BLAST), e.g. `reads.fa` → `reads.fa_mafft.fa`, `reads.fa_flexidot.html`.

| Pattern | Produced by | Contents |
|---------|-------------|----------|
| `<name>_b.txt` | BLAST | Tabular hit table |
| `<name>_b.txt.fa` | BLAST | Hit sequences extracted from the genome |
| `<name>_b.gff` | BLAST | Hit positions (optional) |
| `<name>_mafft.fa` / `<name>_muscle.fa` | Align (each file) | MSA of one FASTA (`_1`, `_2`, … on conflict) |
| `<name>_cmafft.fa` / `<name>_cmuscle.fa` | Align (combined) | MSA after merging selected sequences (`_1`, `_2`, … on conflict) |
| `<name>_gap.fa` | MSA Clean | Gappy columns removed |
| `<name>_rowdiv.fa` | MSA Clean | Divergent rows removed |
| `<name>_rowgap.fa` | MSA Clean | Gappy rows removed |
| `<file.fasta>_<threshold>_co.fa` | Consensus | Consensus sequence |
| `<name>_TEAid.html` | TE-Aid | Interactive TE report |
| `<name>_flexidot.html` | Flexidot | Interactive dotplot |
| `<name>_flexidot.gff` | Flexidot | Match coordinates |
| `<name>_flexidot_legend.png` | Flexidot | Plot legend |
| `<name>_rc.fa` | Reverse Complement | Reverse-complemented sequence |
| `<name>_orfs.fa` | ORF Prediction | Predicted ORFs |
| `<name>_translated.faa` | Translate | Protein translation |
| `<genome>_repeat_quant/` | Repeat Quantify | Folder of reports and tables |
| `<genome>_repeat_quant.html` | Repeat Quantify | Per-family report |
| `<genome>_repeat_quant.tsv` | Repeat Quantify | BLASTn summary table |
| `<genome>_repeat_quant_rm.tsv` | Repeat Quantify | RepeatMasker summary table |
| `<genome>_repeat_quant.gff` | Repeat Quantify | Array coordinates |
| `<genome>_repeat_landscape.html` | Repeat Quantify | Divergence landscape |
| `<genome>_repeat_landscape.tsv` | Repeat Quantify | Landscape data |
| `<name>_{blast|RM}_den_plot.html` | Visualization | Chromosome-scale TE density map |
| `<name>_multi_dens_plot.html` | Visualization | Density map from multiple TE sources |
| `<name>_{blast|RM}_den_plot_region.html` | Visualization | Zoomed region density map |
| `<name>_{blast|RM}_div_plot.html` | TE divergence | Divergence histogram |
| `<name>_{blast|RM|multi}_abund_plot.html` | TE abundance | Family abundance curves |
| `<name>_library_benchmark.html` | Benchmarking | Library comparison |
| `<base>_confusion_matrix.txt` | Benchmarking | Classification matrix |
| `<base>_heatmap_proportion_plot.pdf` | Benchmarking | Proportion heatmap |
| `<base>_spider_plot.pdf` | Benchmarking | Radar comparison |
| `<base>_sankey.html` | Benchmarking | Classification flow |

---

## 16. Troubleshooting

**A tool says the genome is not set.** **File → Manage Database…**. BLAST, MSA
Extend, genome maps and genome-aware TE-Aid all need it.

**RepeatMasker or FastGA is not found.** They are not bundled. Clicking the
FastGA tab shows platform-specific setup steps. Use a shared conda env
`rv_env` with **FastGA 1.5 or newer**
(`conda install -n rv_env -c bioconda fastga=1.5.1` if the env exists, else
`conda create -n rv_env -c bioconda fastga=1.5.1`) and point
**File → Setting → External tools** at it, then use Test to confirm. On Windows
FastGA must live in a WSL2 environment.

**BLAST or another tool is very slow.** Check the word size — small values are
far more sensitive and dramatically slower on repetitive genomes. For Repeat
Quantify, 11 is a reasonable BLASTn default and 7 can effectively hang. Also
check **File → Setting → Performance**: too high a UI thread reserve leaves
few threads for the job.

**A job cannot be stopped.** **Stop** in the status bar cancels the running
job. Long-running external processes may take a moment to terminate.

**Preview is blank for an HTML report.** Use **Open in browser** on the
preview header. Some environments have limited web rendering support.

**Recent outputs are not visible.** Press **F5** to refresh, or check you are
looking at the folder holding the input file — outputs go there, not into a
central results folder.

**"appears 2× in the library" in a Repeat Quantify report.** Two records in
your library share a name. Their hits are merged into one row and copy number
uses a single monomer length. Rename them to keep them separate.

**macOS refuses to open the app.** For an unsigned build, right-click the app
→ **Open**, then confirm. Only needed the first time.

**Reporting a bug.** **File → Setting → General → Copy diagnostics** puts
version, OS and folder locations on the clipboard. Enable DEBUG logging and
attach the log file from **Open log file**.

---

## 17. Citation

If RepeatViewer contributes to your work, please cite **TEtrimmer**
(*Nature Communications*, 2025) where applicable, and mention **RepeatViewer**
for manual curation and visualization.
