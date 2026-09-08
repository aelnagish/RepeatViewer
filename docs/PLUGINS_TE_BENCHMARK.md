# TE benchmarking plugin (TEtrimmer)

This wraps **`tetrimmer/TE_anno_benchmarker.py`** from the [TEtrimmer](https://github.com/qjiangzhao/TEtrimmer) repository (`calculate_overlap` CLI).

It compares **query** RepeatMasker `.out` files to a **reference** `.out` (curated TE library) and writes:

| Output | Description |
|--------|-------------|
| `{query}_sankey.html` | Interactive Sankey (open in RepeatViewer Preview) |
| `{query}_heatmap_proportion_plot.pdf` | Heatmap |
| `{query}_confusion_matrix.txt` | Confusion metrics |
| `{query}_confusion_matrix_TP_TN_FP_FN.txt` | TP/TN/FP/FN |
| `{query}_spider_plot.pdf` | Spider plot |
| `{query}_TE_length.txt` | Length table |

Outputs are saved **next to each query `.out` file** (same folder as in the Files panel).

## Prerequisites

1. **TEtrimmer repository** cloned locally (contains `tetrimmer/TE_anno_benchmarker.py`).
2. **bedtools** on `PATH`.
3. **Genome FASTA** set in RepeatViewer **File → Project settings** (same as Blast/TEAid).
4. **Reference** RepeatMasker `.out` from your gold-standard TE library (set once in the plugin dialog).

## Install the plugin

```bash
cd /path/to/RepeatViewer/examples/repeatviewer_benchmark_plugin
pip install -e .
```

Optional — point to TEtrimmer without using the dialog every time:

```bash
export TETRIMMER_HOME=/path/to/TEtrimmer
```

Restart RepeatViewer. A **TE benchmark** button appears on the toolbar (after built-in tools).

## How to run in RepeatViewer

1. In **Files**, turn off **FASTA only** if needed (or keep it on and only select `.out` files).
2. Select one or more **query** `.out` files (RepeatMasker output to evaluate).
3. Click **TE benchmark** (or **Sequence → Plugins → TE benchmark**).
4. In the dialog:
   - **TEtrimmer folder** — root of the TEtrimmer git repo (e.g. `~/Documents/GitHub/TEtrimmer`)
   - **Reference RM .out** — curated library RepeatMasker output
   - Optional: simpler TE classes, debug mode
5. Confirm. Each selected query runs sequentially; the Sankey HTML opens in Preview when done.

## Standalone CLI (same backend)

From the TEtrimmer repo:

```bash
python tetrimmer/TE_anno_benchmarker.py \
  -q query.out \
  -r reference.out \
  -g /path/to/genome.fa \
  -o /path/to/output_dir
```

## Package layout

```
examples/repeatviewer_benchmark_plugin/
  pyproject.toml          # entry point: repeatviewer.plugins → te_benchmark
  src/repeatviewer_benchmark/
    plugin.py             # RepeatViewerPlugin subclass
    dialog.py             # TEtrimmer path + reference .out
    tetrimmer_runner.py   # import TE_anno_benchmarker by file path
```

## Customization

- Change `file_suffixes` in `plugin.py` if you use another naming convention.
- Set `supports_parallel = True` only if you verify thread-safety of matplotlib/plotly in your environment.
- For batch query + single reference, the default `run()` per file is the intended workflow.

## Troubleshooting

| Problem | Fix |
|---------|-----|
| Plugin not listed | `pip install -e .` in the plugin folder; restart RepeatViewer |
| TEtrimmer not found | Set **TEtrimmer folder** in dialog or `TETRIMMER_HOME` |
| bedtools missing | `conda install -c bioconda bedtools` or system package |
| No genome | Project settings → Genome FASTA |
| Wrong file type | Select `.out` files, not FASTA |
