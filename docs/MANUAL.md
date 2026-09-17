# RepeatViewer User Manual

**Version 1.0.0**

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
| Reference **genome FASTA** | BLAST, MSA Extend, genome maps, genome-aware TEAid |
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

### Step 2 — Define the genome

**File → Define genome** (⇧⌘G / Ctrl+Shift+G).

Select the reference genome FASTA these TEs came from. Several tools refuse to
run without it — BLAST against the genome, MSA Extend, genome maps, and the
genome tracks in TEAid — because they need somewhere to look the copies up.

You can define a **pan-genome** instead: a list of several genome FASTAs,
searched together.

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
source: `<yourfile>_blast.txt`, the hit table you can read, and
`<yourfile>_blast.txt.fa`, the hit sequences pulled out of the genome. The
second one is what you align in the next step.

New outputs are highlighted **green** in the Files panel; the inputs that
produced them are highlighted **blue**. You do not have to hunt for what just
appeared.

### Step 5 — Build and clean an alignment

Click **Align** to align the BLAST hits into a multiple sequence alignment
(`<yourfile>_msa.fa`).

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

Run **TEAid** for an interactive HTML report of the new consensus — coverage,
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
- **HTML reports** — TEAid, Flexidot, genome maps, quantification reports,
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
TEAid, MSA Clean, MSA Extend, Consensus and Visualization.

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

### Define genome

**File → Define genome** (⇧⌘G / Ctrl+Shift+G)

Sets the reference genome. Two modes:

- **Single genome** — one FASTA.
- **Pan-genome** — several FASTAs treated as one search space.

RepeatViewer prepares an index and reads chromosome names and sizes the first
time; this takes a moment on a large assembly and is then cached.

Tools that require a genome will tell you to set one rather than failing
silently.

### Define database

**File → Define database**

Configures the annotation databases used by TE annotation:

- **CDD** — conserved domain database, for protein domain annotation
- **Similarity databases** — reference libraries for classification

Saved BLAST databases are also managed here. Any FASTA in the Files panel can
be turned into one with **Add as BLAST database…** from its right-click menu,
after which it appears in the **Database** dropdown in the BLAST dialog.

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
**BLASTx** (protein).

| Parameter | Meaning |
|-----------|---------|
| Database | Reference genome, or a saved BLAST database |
| `-evalue` | Expectation threshold |
| `-word_size` | Seed length. Smaller is more sensitive and much slower on repetitive genomes |
| `-gapopen`, `-gapextend` | Gap penalties |
| `-penalty` | Mismatch penalty (BLASTn) |
| `-perc_identity` | Minimum percent identity |
| `-qcov_hsp_perc` | Minimum query coverage per HSP |
| Threads | Defaults to your CPU count minus the UI reserve |
| Max MSA lines | Cap on hits carried into an alignment |
| Top longest lines | Prefer the longest hits when capping |
| Minimum bit score (GFF) | Score floor for hits written to GFF |

**Outputs:**

| File | Contents |
|------|----------|
| `<name>_blast.txt` | Tabular hit table |
| `<name>_blast.txt.fa` | Hit sequences extracted from the genome — the input to Align |
| `<name>_blast.txt.gff` | Hit positions, when GFF output is enabled |

The hit table is ordinary tabular BLAST output — readable in the preview pane,
and the input for the TE density and divergence plots. **Max MSA lines** and
**Top longest lines** cap how many hits are carried into the `.fa`, which
matters because a high-copy element can otherwise produce an alignment too
large to work with.

### Align

**Toolbar: Align · Tools → Align**

Builds a multiple sequence alignment with **MAFFT** or **MUSCLE**, both
bundled. The **Options** tab selects the default engine.

MAFFT parameters:

| Parameter | Meaning |
|-----------|---------|
| MAFFT strategy | Algorithm preset |
| Sequence type | Nucleotide or protein |
| `--retree` | Guide-tree rebuilds |
| `--maxiterate` | Refinement cycles |
| `--op` | Gap-opening penalty |
| `--ep` | Offset value |
| Threads | Parallelism |

MUSCLE exposes its algorithm choice.

**Output:** `<name>_msa.fa`

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
| Majority threshold | Fraction of rows that must agree for a base to be called |
| Ambiguous base | Character written where no base reaches the threshold |

Consensus needs either a **multiple sequence alignment** — every row the same
length — or a single-sequence FASTA. An unaligned multi-FASTA is rejected with
a message telling you to align it first.

**Output:** a consensus FASTA whose name keeps the full input filename and
records the threshold, e.g. `CroSat2.fasta_0.8_co.fa`. Several thresholds can
therefore coexist in one folder and be compared.

### TEAid

**Toolbar: TEAid · Tools → TEAid**

An interactive HTML report for one consensus: self-alignment structure,
terminal repeats, ORFs, and genome hit coverage in one page. This is the
clearest single view of whether a consensus is finished.

| Parameter | Meaning |
|-----------|---------|
| Genome BLAST E-value | Threshold for the genome hit track |
| Self BLAST word size | Seed length for the self-comparison |

**Output:** `<name>_TEAid.html`, rendered in the preview pane.

### Flexidot

**Toolbar: Flexidot · Tools → Flexidot**

Dotplots: self, pairwise, or all-against-all for a multi-FASTA.

| Parameter | Meaning |
|-----------|---------|
| Word size (k-mer, `-k`) | Match length. Lower finds more, and more noise |
| Output format (`-f`) | HTML for interactive, or a static image |
| Allowed substitutions (`-S`) | Mismatches tolerated within a word |
| Maximum N content (`--max_n`) | Skip windows above this N fraction |
| Match line thickness (`-A`) | Line width in the plot |

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

**Toolbar: FastGA · Tools → FastGA**

Fast genome-scale alignment, for comparing large sequences where BLAST is too
slow. Requires FastGA on `PATH`.

| Parameter | Meaning |
|-----------|---------|
| Output format | Alignment output type |
| Min identity | Identity floor |
| Min alignment length | Length floor |
| Repeat cutoff (`-f`) | Repeat masking threshold |
| Threads | Parallelism |

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

Tools producing `html` or `plot` are *terminal* — TEAid and Flexidot end a
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
| TEAid | fasta | html | yes |
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
**Reset to default** restores BLAST, Align, Flexidot, TEAid, MSA Clean, MSA
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
| ⇧⌘G / Ctrl+Shift+G | Define genome |
| ⌘F / Ctrl+F | Focus the file or folder filter, depending on the active panel |
| F5 | Refresh the current folder |
| Delete | Delete selected files, with confirmation |
| ⌘Q / Ctrl+Q | Quit |

---

## 15. Output file reference

All outputs are written beside their input file unless noted.

| Pattern | Produced by | Contents |
|---------|-------------|----------|
| `<name>_blast.txt` | BLAST | Tabular hit table |
| `<name>_blast.txt.fa` | BLAST | Hit sequences extracted from the genome |
| `<name>_blast.txt.gff` | BLAST | Hit positions (optional) |
| `<name>_msa.fa` | Align | Multiple sequence alignment |
| `<name>_gap.fa` | MSA Clean | Gappy columns removed |
| `<name>_rowdiv.fa` | MSA Clean | Divergent rows removed |
| `<name>_rowgap.fa` | MSA Clean | Gappy rows removed |
| `<file.fasta>_<threshold>_co.fa` | Consensus | Consensus sequence |
| `<name>_TEAid.html` | TEAid | Interactive TE report |
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

**A tool says the genome is not set.** **File → Define genome**. BLAST, MSA
Extend, genome maps and genome-aware TEAid all need it.

**RepeatMasker or FastGA is not found.** They are not bundled. Install them
and either put them on `PATH` or point **File → Setting → External tools** at
the conda environment holding them, then use the check button to confirm. On
Windows they must live in a WSL2 environment.

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
