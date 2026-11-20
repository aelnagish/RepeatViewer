# Primers and PCR

RepeatViewer designs PCR primers and FISH probes for repeats, shows primers on
the sequence, runs in-silico PCR, draws the expected gel, and keeps the lab's
primer inventory. Everything runs on your computer; nothing is sent online.

| Where | What |
|-------|------|
| **Tools → Primers and Probes → Design** | Design primers (**Primers**), FISH probes for a repeat (**Repetitive Probes**) or painting probes (**Chromosome Painting Probes**) |
| **Tools → Primers and Probes → Map** | Test primers or probes against sequences, like Geneious' *Map Primers* |
| **Tools → Primers and Probes → Database** | The lab's primer inventory |
| **Tools → Primers and Probes → In-silico PCR**, **Restriction Digest**, **Gel** | The products of a PCR, a digest of the sequence, the expected gel |
| **Tools → Primers and Probes → Sanger Check**, **Cloning** | Check a clone against its Sanger reads; plan putting a product into a vector |

Every item opens the same tab, **Primers and Probes**, on the matching section.

The toolbar can also show a **Primers and Probes** button (**File → Setting →
Toolbar → Add buttons…**). It offers the same items; right-click designs with
the saved settings.

## Quick start

1. Click the FASTA you want primers for in the Files panel.
2. Choose **Tools → Primers and Probes → Design → Primers**. The **Primers and
   Probes** tab opens on the **Design** section, with your sequence.
3. Drag on the map over the part to amplify, keep **Mode** on *Standard PCR*
   and click **Design**.
4. The candidate pairs appear in the list under the map and as arrows above the
   sequence. Tick the pairs you want and click **Keep Ticked**: their primers
   are now on the sequence.
5. Click the **PCR** section to see the products of the ticked primers and the
   expected gel.

## The workspace

**Primers and Probes** is one tab. The row on top holds its sections: **Design**,
**Painting**, **Map**, **Database**, **PCR**, **Digest**, **Gel**, **Sanger
Check** and **Cloning**. A section keeps what you did in it while you look at
another one.

**Design**, **PCR** and **Digest** work on one sequence:

- **Sequence** shows the file. Click a FASTA's name in the Files panel to work
  on it, or press **Browse**. Ticking a file's box does not change it.
  **Record** chooses the sequence in a file that holds several.
- The **map** shows the sequence. What a section adds is drawn in rows of its
  own above the bases, so no letter is covered: enzyme names, one row for each
  candidate pair, and the primer you are editing. A region you chose is a light
  tint behind them.
- The **list under the map** holds what Design found (*candidate*) and the
  primers already on the sequence (*on the sequence*); in **PCR** it holds the
  products. Drag the border between the map and the list to give one of them
  more room.
- The **panel on the right** holds the settings of the section.

The other sections are whole pages of their own. **Map** and **Database** open
without a sequence.

> **Tip.** A path field in the panel (for example *Avoid file*) is filled by a
> click in the Files panel: put the cursor in the field first, then click the
> file. To open a file in a tab of its own while the workspace is in front,
> right-click it and choose **Open in Additional Tab**.

## Design

The settings are in the panel on the right of the **Design** section.

**Run on**: **Open sequence** designs on the open sequence. The candidates appear in
the list under the map and nothing is written until you press **Keep Ticked**
or **Open Report**. **Selected files** runs the design as a job on every FASTA
ticked in the Files panel and writes the report files beside them, with the
same settings.

**Mode** is the kind of design. Changing it loads sensible starting values
(size, Tm, product size); **Restore Defaults** brings them back.

| Mode | Use it for | Template | You get |
|------|------------|----------|---------|
| Standard PCR | An ordinary pair | Any sequence | Ranked pairs and a virtual gel |
| qPCR | Quantitative PCR | The target; optionally a single-copy reference gene FASTA | Short amplicons (70 to 180 bp), an optional TaqMan probe, and a matched pair on the reference gene for ΔCt |
| Ladder PCR | Showing that a satellite sits in tandem arrays | Monomers (one or many) | Outward-facing primers inside one monomer. On genomic DNA they give a ladder: 1×, 2×, 3× … the monomer |
| Satellite repeat oligos | FISH probes for a satellite or repeat family | Monomers or copies, or a TAREAN dimer consensus | 35 to 60 nt oligos present in at least *Min coverage* of the copies, without hairpins, with a Tm for your hybridisation conditions, ready to order with a 5′ label (FAM, Cy3, Cy5 …) |
| Chromosome painting (Chorus2) | Oligo pools that paint whole chromosomes or regions | The genome FASTA | Thousands of single-copy 45 nt oligos along the chromosomes, an order file and a density plot (see below) |

For ladder PCR and FISH oligos, use many monomers (for example extracted from an
array) rather than one: the tool then keeps the parts shared by most copies. A
dimer consensus (TAREAN's `consensus_dimer.fasta`, the monomer written twice) is
recognised and cut back to one monomer.

**Satellite repeat oligos** have their own settings:

| Field | Meaning |
|---|---|
| Avoid the other selected files | Select the target family first, then the other families (one FASTA each). Probes are designed on the first file and every probe that matches the others is ranked last; the report names the sequences it matches |
| Leave out low-complexity and inverted repeats | Low-complexity stretches (DUST, as `dustmasker -level 30`), simple repeats and inverted repeats (arms of 12 nt or more, as EMBOSS `einverted`) are found first and no probe may touch them; the report lists them |
| Hybridisation | Formamide, SSC and temperature of the hybridisation (default 50 % formamide, 2× SSC, 37 °C) |
| Stringent wash | The same for the stringent wash (default 20 % formamide, 0.1× SSC, 42 °C); probes that would melt off in it are not offered |
| Probe | Probe concentration in the mix (default 100 nM) |

For these probes **Tm** is the Tm under the hybridisation conditions, not in a
PCR buffer: SantaLucia (1998) nearest-neighbour values with the Owczarzy et al.
(2004) Na⁺ correction (1× SSC = 195 mM Na⁺), minus 0.60 °C per % formamide
(Hutton 1977). The default window (48 to 58 °C, about 52 °C for 50 to 60 nt)
follows the TAREAN probe design protocol. The report gives each probe's Tm in
the wash and how far it stays above the wash temperature.

> **Tip.** Never use a PCR product amplified from genomic DNA as a hybridisation
> probe, even with one band of the expected size: clone and sequence it first,
> and keep short simple repeats in the clone out of the probe.

**Primer Conditions**:

| Field | Meaning |
|-------|---------|
| Primer size | Minimum, optimum and maximum length (nt) |
| Tm (°C) | Minimum, optimum and maximum melting temperature of the binding part |
| GC | Allowed GC content (%) |
| Product size | Allowed amplicon length (bp) |
| Max Tm difference | Largest Tm gap between the two primers of a pair |
| Pairs to return | How many pairs the report lists |

**Advanced Settings ▸** adds poly-X, GC clamp, 3′ stability, probe size and Tm,
5′ tails (M13, T7, restriction sites, or your own) and the Tm calculation (salt,
Mg²⁺, dNTP and oligo concentration).

**Options**: **Ask which pairs to add here** (for a run on selected files);
**Add them to the Database**. **Restore Defaults** brings the starting values
of the mode back, **Save** keeps the settings for the next time.

### Design results

Next to the input:

```
gene.fa_primers.html              the report (opens in Preview)
gene.fa_primers_details/
   gene.fa_primers_primers.fa     every primer, with order-ready names
   gene.fa_primers_pairs.tsv      table of pairs
   gene.fa_primers_gel.svg        virtual gel
   gene.fa_primers.json           all results
```

The report lists the pairs (Tm, GC, hairpin and dimer ΔG, product size, notes),
the virtual gel and a **Methods text** paragraph for a paper. When nothing is
found it explains why (for example how many primers failed on Tm or GC). Common
fixes: widen the Tm range or the product size; for satellite oligos lower *Min
coverage* or allow more mismatches.

## Design on the sequence

1. **Region**: drag on the map or the bases; the range appears in *From* and
   *to* (or type it). Choose what the region is for:
   - **Amplify**: the primers sit outside the region, one on each side, and the
     product covers all of it.
   - **Inside**: both primers and the product stay inside the region.
   - **Avoid**: no primer may sit in the region (a SNP, a methylated site, a
     repeat).

   For probes only **Inside** and **Avoid** are offered: a probe sits in the
   region or keeps out of it.
2. **Settings**: the mode, primer size, Tm, GC and product size; **Advanced**
   has the finer rules and the salt.
3. **Design**: the pairs appear in the list under the map, best first, and
   each pair gets a row of arrows above the sequence (the selected pair is
   filled, the others are outlined). Click a pair to see it; tick the ones you
   want and press **Keep Ticked**. They are added to this sequence's primers.
   Nothing is written before that. **Open Report** writes the usual primer
   report.

When the region and the settings cannot work together (for example a 1,500 bp
region to amplify with products of at most 1,000 bp), an amber note says why
before you run, and **Fit** changes the product range. Nothing is blocked: the
settings guide the automatic design, your choice wins.

## Primers on the sequence

The list under the map shows the primers kept with this FASTA as *on the
sequence*. Each is drawn as an arrow where it binds (blue: points right;
orange: points left).

- **Add Primer** types a new one, **From Database** picks saved ones.
- **Click** a primer to see its details in the panel on the right and jump to
  it: Tm, hairpin and self-dimer ΔG, where and how well it binds. While a
  primer is selected the design settings are folded to one line, **Design
  Settings**; **Show** opens them again. Clicking a candidate also brings them
  back.
- The selected primer is drawn as a bar in its own row right above its bases.
  **Drag** the bar to move the primer, or drag one of its ends to make it
  longer or shorter. A base that does not match the template is red on the
  primer.
- **Snap to Tm** shortens or lengthens the primer along the template to reach a
  Tm; **Fix Mismatches to Match Template** makes it match this sequence;
  **Edit Primer** opens the primer with its tail and checks; **Sites in Genome**
  shows where it binds in a genome.
- **Right-click** a primer in the list to edit it, copy its sequence, add it to
  the Primer Database or remove it.
- The **tick** marks the primers that the **PCR** section uses.

A **tail** is the 5′ part that does not bind (lower-case letters). It is copied
into the product from cycle 3 on, so two Tm values are shown: **annealed** (the
binding part, mismatches included) and **full** (the whole oligo).

**Edit Primer** shows the primer on its best site, with mismatches in red and
matching 3′ bases in green; **Site** steps through other sites; **− +** make
either end shorter or longer along the template; **5′ tail** adds a restriction
site (with padding bases so the enzyme cuts), a sequencing or promoter tail, or
your own; and the dimer and hairpin checks give ΔG at 37 °C (below −9 kcal/mol
for dimers or −3 for hairpins is a warning).

Primers are saved in a hidden `name.primers_RV.json` next to the FASTA and
follow it when it is renamed or moved.

## In-silico PCR

**Tools → Primers and Probes → In-silico PCR**, or the **PCR** section of the
workspace. The products are listed under the map and drawn in the expected gel
in the panel.

![The PCR window: templates and primers on the left, products in the middle, a quick gel on the right.](img/pcr.png)

- **Templates**: click FASTA files in the Files panel to add them (a whole
  genome works). Each template is one lane.
- **Primers**: from a sequence's Primers tab, **From Database…**, or
  **Paste…**. All combinations are tested, as in a real tube, including one
  primer at both ends.
- **Kind of PCR**:
  - **Standard**: every product.
  - **Satellite ladder**: products that step by a common size are labelled 1×,
    2×, 3× (the monomer of a tandem repeat).
- **Mismatches per primer** and **3′ end exact** decide which sites count;
  **Longest product**, **Polymerase** and **Annealing** set the reaction.

The products table lists each product with its size, where it is and the
mismatches of each primer; **Expected yield** is relative (fewer mismatches,
shorter products and repeated sites give more DNA). The strip above the table
shows the product you click: the forward primer (blue), the reverse primer
(orange) and the stretch between. Small arrows are other places the primers
bind. **Zoom to Product** and **Show All Products** frame the strip; **Copy
Product Sequence** copies it.

Buttons at the bottom:

| Button | Opens |
|--------|-------|
| **Gel** | The full gel document (below) |
| **Products File** | The products FASTA, with primers and the template's annotations |
| **Table** | The products table |
| **Bench Card** | Mastermix volumes and the cycling program for the bench |
| **Sites in Genome** | Every binding site of each primer in a genome |

**Outputs**, next to the first template: `<template>_PCR_<primers>.fa` (the
products) and `…_details/` (`_products.tsv`, `_history.json`,
`_result.json`). The products file opens like any multi-FASTA, with the primer
sites and the template's features on each product.

![The products of a satellite ladder PCR: 1×, 2×, 3× … monomers, with the primers and the monomer annotations carried over.](img/pcr_products.png)

**Bench Card**: volumes per reaction and for all reactions plus extra (template
added to each tube); annealing from the lower primer Tm (Taq 3 °C below,
Phusion or Q5 3 °C above); extension from the longest product wanted (Taq
1 min per kb, Phusion or Q5 20 s per kb); freezer places from the Primer
Database. Change the numbers on the right, then **Save PDF** or **Print**.

**Sites in Genome**: each genome sequence drawn on one scale, forward sites as
ticks above the bar, reverse sites below (darker = fewer mismatches), and for
a qPCR pair the number of products per haploid genome, the value to expect in a
copy-number assay.

## The gel

A `.rvgel` file opens as an agarose gel, in the **Gel** section of the workspace
(every **Open in a Gel** button, **Tools → Primers and Probes → Gel**, or a
click on a gel in the Files panel while the workspace is in front). Every lane is a ladder, the products of
a PCR run, or products cut with a restriction enzyme. The gel is saved as it
changes.

![The in-silico gel of a satellite ladder PCR: GeneRuler 1 kb Plus next to the 1×, 2×, 3× … products.](img/gel.png)

| Control | What it does |
|---------|--------------|
| **Add Lane**, **Duplicate**, **Remove**, **◀ ▶** | Edit and order the lanes |
| **Ladder** | The size marker (GeneRuler, NEB and others) |
| **Agarose** | Which sizes separate: 0.7 % for 1 to 12 kb, 3 % for 50 to 1,000 bp |
| **Run** | How far the front got; small bands can run off a long run |
| **Look** | Dark (SYBR, GelRed) or light (stained) gel |
| **Brightness by amount of DNA** | A band's brightness follows its mass (amount × length), as with SYBR, so short products look fainter |
| **Show primer-dimer** | Draw the primer-dimer smear |
| **Cut with** (per lane) | Digest the lane's products with a restriction enzyme |

Hover a band to see what it is, click to select it, double-click to open the
products file.

**Compare with Gel Photo…** puts a picture of your real gel (PNG, JPG or TIFF)
next to the in-silico one. Lanes and bands are found for you. Choose which
photo lane holds the ladder (its bands give the size scale) and which gel lane
each photo lane should match. Every band is then named: **matched** (within
7 %), **off-target** (a predicted off-target product), **unexplained** (try the
PCR with more mismatches), and expected bands **not seen**.

## Map: test primers on sequences

**Tools → Primers and Probes → Map → Primers** (or **Repetitive Probes**). It
opens as the **Map** section of the workspace, with its settings on the right.

![The Map Primers window.](img/primer_map.png)

- **Inputs** as in Design. **Test primers on**: every sequence, or the
  consensus of an alignment.
- **Primers**: **Test specific primers** (type a sequence or pick a saved one),
  or **Search for saved primers that match** (from the Primer Database or a
  primers FASTA; forward, reverse, probe; pairs only).
- **Products**: optional included and target regions, product size range,
  optimal size.
- **Mismatches**: maximum mismatches, and how many 3′ bases must match exactly.
- **Advanced Settings**: Tm conditions, and counting sites and products in a
  genome (both strands; every site within the allowed mismatches is found, also
  on a whole genome).

**Output:** `<name>_primer_map.html` with binding sites, predicted products and
a gel.

## Database: the lab's primer inventory

The Database is the record of every primer the lab has: sequence and Tm, what it
amplifies and in which organism, where the tube is, whether it is in stock and
whether it worked, who ordered it and when, and every change since.

![The Primer Database window.](img/primer_database.png)

It opens as the **Database** section of the workspace.

### The page

- **Filter** searches every field; the boxes next to it show one organism, one
  person or one status.
- Select a primer to edit it below the table, then **Save Changes**. Select
  several to set a field on all of them at once (for example a new box, or
  *Used up*).
- **History** lists every change to the selected primer: who, when, what.
- **Pair** is found from the names (`_for`/`_rev`, `_fwd`/`_rev1`,
  `_FW`/`_RV`, `_F`/`_R`); write it in *Pairs with* when the names do not show
  it.
- **Two Tm columns**: *Tm Primer3* is calculated by RepeatViewer; *Tm supplier*
  is the value from the synthesis report. They differ because each uses its own
  method and concentrations.
- **Lab** holds the lists the fields choose from: **People** (with roles;
  rename, merge two entries for one person, or remove; people who ordered
  primers stay as former members), **Places** (rooms, freezers, boxes),
  **Organisms** and **Targets**.
- **Find Duplicates**: the same sequence under two names, reverse complements,
  one primer inside another.
- **Export**: the whole table (TSV), a FASTA, or an order sheet.
- **Import**: a CSV, TSV or Excel sheet. The header row is found and known
  columns are read (name, sequence, organism, target, Tm, dissolved volume,
  ordered by, ordered on, box, freezer, room, status, tested, annealing,
  concentration, pair, note). Before importing it shows what is new, what is
  already there (same sequence under another name is kept as *Other names*)
  and which rows are skipped. Importing again later only adds what is new.

### Sharing with the lab

1. The button at the bottom left says where the database is:
   **Database: This Computer ▾**, or the name of a shared folder.
2. One person: **Database ▾ → Shared Database…** and choose the lab's folder
   in a synced cloud folder (Sciebo, Dropbox). If there is no database there
   yet, RepeatViewer offers to make one and to copy your primers into it.
3. Everyone else: **Database ▾ → Shared Database…** and choose the same lab
   folder in their own synced folder.
4. If the shared database cannot be found (not synced yet, moved), a bar says
   so and offers **Try Again** or **Use This Computer**; nothing switches
   silently.
5. The first time you change something, you are asked for your name; it
   appears in the history (**Lab → Change Your Name…**).

Each person and computer writes only its own change log, so the cloud service
never has to choose between two versions of a file and no change is lost. If
two people change the same field of one primer, the later change wins.

### Check a new sequence against the stock

Before designing and ordering for a new plant, see what the lab already has:
select the primers to test (none = all) and click **Check a Sequence**. Give
the **Sequence** (click a FASTA in the Files panel) and how strict to be
(product size; mismatches, by default up to 2 but with the last 3 bases at the
3′ end matching). The report starts with **From your Primer Database**:

- **Ready pairs**: two saved primers that give a product, with size,
  mismatches and where the tubes are. Saved pairs, exact matches and pairs that
  worked before come first.
- **Single primers that bind**: no saved partner gives a product; one new
  partner is enough.

## Chromosome painting (Chorus2)

Choose **Tools → Primers and Probes → Design → Chromosome Painting Probes**. It
opens as the **Painting** section of the workspace.
[Chorus2](https://github.com/zhangtaolab/Chorus2) (Zhang et al. 2021) finds
oligos that occur once in the genome along every chromosome. RepeatViewer runs
it, or reads what it produced on a server, and then picks the pool you order.

**Chorus2 tab**:

- **Run Chorus2 here**: give the **Genome**; optionally paint only one
  chromosome or region. A genome of a few Gb takes hours and much memory: a
  server is better.
- **Import results from the server**: choose the Chorus2 BED.
- Chorus2 settings: probe length (45 nt), maximum homology (75 %), minimum dTm
  (10 °C), k-mer step (5 bp), ploidy, threads.
- **NGS Filter** (optional): genomic reads give each probe a k-mer score;
  probes in the lowest and highest score quantiles (default 10 % and 90 %) are
  dropped.
- **Copy Server Command**: the exact command lines for the server.

Install Chorus2 once, in its own conda environment (it needs an older pandas):

```bash
conda create -n chorus2 -c conda-forge -c bioconda python=3.9 "pandas<2" numpy cython pyfasta "primer3-py<2" pybigwig pybedtools kmer-jellyfish bwa
conda run -n chorus2 pip install --no-deps git+https://github.com/zhangtaolab/Chorus2.git
```

On a server:

```bash
conda activate chorus2
Chorus2 -g genome.fa -i genome.fa -s chorus2_out -t 16 -l 45 --homology 75 -d 10 --step 5 --ploidy 2 --docker True
```

Copy `chorus2_out/genome.fa.bed` back (or `probes_kmer_score.bed` after the
NGS filter) and import it. On Windows, run Chorus2 on a Linux server and import.

**Selection tab**:

| Field | Meaning |
|-------|---------|
| Min gap between probes | Neighbouring probes at least this far apart (25 bp) |
| Strand | **Both** (every second probe reverse-complemented, the Chorus2 default) or **Plus** |
| Amount | **All probes**, a **Density** (probes per kb; painting usually needs 0.1 to 1 per kb), or a **Pool size** (total number, spread evenly) |
| Regions | Keep probes only here, for example `chr1, chr2:1-5000000` |
| Pool name, 5′ / 3′ adapter | Used in oligo names; adapters (for example T7) go into the order file only |

**Results**:

```
genome.fa_painting.html             report: counts, density plot, per chromosome, methods
genome.fa_painting_details/
   genome.fa_painting_probes.bed      probes (BED, for example for IGV)
   genome.fa_painting_probes.fa       probe sequences
   genome.fa_painting_order.tsv       order file (adapter + probe + adapter)
   genome.fa_painting.json            all results
   chorus2/                           Chorus2's own files (when run here)
```

Painting pools are not added to the Primer Database (they hold thousands of
oligos); use the order file. Too few probes for the pool size: lower the gap,
widen the regions or the NGS score window, or rerun Chorus2 with a smaller
k-mer step or a higher maximum homology.

## How it works

Primer picking, Tm, and the hairpin and dimer checks use
[Primer3](https://github.com/primer3-org/primer3) (primer3-py). The
binding-site search (Map, coverage on copies, genome counts) is RepeatViewer's
own, on both strands: small inputs are scored at every position; on a genome,
every version of each primer's 3′ end within the allowed mismatches is looked
up, so the same sites are found. Ladder PCR and FISH oligos build on these.
Ladder primers are picked on the dimer with Primer3's template mispriming check
switched off, because every primer in a tandem repeat binds the template many
times.
Chromosome painting uses Chorus2 (MIT licence) in its own conda environment.

References: Untergasser et al. (2012) *Nucleic Acids Research* 40:e115
(Primer3). Owczarzy R et al. (2004) *Biochemistry* 43:3537–3554 (salt
correction). Hutton JR (1977) *Nucleic Acids Research* 4:3537–3555 (formamide).
Morgulis A et al. (2006) *Journal of Computational Biology* 13:1028–1040 (DUST).
Zhang T, Liu G, Zhao H, Braz GT, Jiang J (2021) *Plant
Biotechnology Journal* 19:1967–1978 (Chorus2).
