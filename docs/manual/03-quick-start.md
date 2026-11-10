# Quick start: one curation pass

This chapter walks through one complete pass on one TE candidate: find its
copies in the genome, align them, clean the alignment, fix the ends, rebuild
the consensus and check it. It takes about fifteen minutes and touches most of
the tools you will use every day.

## Step 1: add your folder

Choose **File → Add Folder…** (**⇧⌘O** / **Ctrl+Shift+O**), click **+** in the
Folders panel title bar, or drag a folder from Finder or Explorer onto the
Folders panel.

Pick the folder that holds your candidate FASTA files: a TEtrimmer output
folder, or any folder you organised yourself. It is added as a top-level folder
in the Folders panel and is remembered between sessions.

## Step 2: register your genome

Most curation tools compare your candidate with the genome it came from.
Choose **File → Manage Database…** (**⇧⌘G** / **Ctrl+Shift+G**), click
**Add…**, choose the genome FASTA with type **FASTA DNA**, and click **Set as
reference genome**. RepeatViewer builds the BLAST index once, in the
background. See [Manage Database](07-search.md#manage-database) for details.

## Step 3: look at the candidate

Click a FASTA file in the Files panel. The Preview panel shows it. For a single
sequence there are three views:

- **Text**: the sequence itself.
- **Overview**: a map of the sequence with its annotations (LTRs, TIRs, ORFs,
  domains). The side panel on the right detects, adds and edits features.
- **Dotplot (Self)**: the sequence against itself. Terminal repeats show as
  short diagonals in the corners; tandem repeats as parallel lines.

Run **Feature detection** in the side panel (**Annotation** tab) to find LTRs,
TIRs, tandem repeats and ORFs in a few seconds. See
[Annotations](06-annotations.md).

## Step 4: find the copies

With the file selected, click **BLAST** on the toolbar. RepeatViewer searches
the candidate against the reference genome and writes two files next to it:

- `<name>_bn.txt`: the hit table.
- `<name>_bn.txt.fa`: the hit sequences cut out of the genome. This is the
  input for the next step.

New files are highlighted **green** in the Files panel, and the file they were
made from **light blue**, so you can see at once what just appeared.

## Step 5: align and clean

Select the `_bn.txt.fa` file and click **Align**. MAFFT (bundled) aligns the
copies and writes `<name>_mafft.fa`, which opens in the alignment viewer.

Genomic copies always include rows that do not belong: fragments, unrelated
insertions, rows that are mostly gaps. Click **MSA Clean** to remove them.
Each cleaning step writes its own file, so you can see what each step removed.

## Step 6: fix the ends

If the alignment runs to its edge and the copies keep agreeing, the candidate is
too short. Click **MSA Extend** to add genomic sequence on both sides of every
copy, then clean again and look at the new ends. Repeat until the alignment ends
where the copies stop agreeing.

The **Boundary prediction** section in the alignment viewer's side panel
suggests where the element starts and ends (from LTRs, TIRs, alignment
identity, genome coverage or a poly-A tail) and draws these positions on the
alignment.

## Step 7: rebuild the consensus

Click **Consensus**. The alignment is collapsed into one sequence at the
majority threshold you set (default 0.8), for example
`CroSat2.fa_mafft.fa_0.8_co.fa`. Several thresholds can sit next to each other
in the folder for comparison.

## Step 8: check the result

- **TE-Aid** makes an HTML report of the new consensus: genome copies, their
  coverage along the consensus, self-dotplot, terminal repeats and ORFs.
- **Dotplot** compares the old and the new consensus.
- **Feature detection** on the new consensus should now find complete terminal
  repeats.

If the ends are clean and the terminal repeats complete, the candidate is
finished. If not, go back to step 6.

> **Tip.** Once the steps are routine, build them into a workflow
> (**Workflows → Design and manage your workflows**) and run them on a whole
> folder at once. See [Automation](13-automation.md).
