# RepeatViewer plugins

Plugins add **new toolbar tools** that use the same file list and job flow as Blast, TE-Aid, Flexidot, etc.

> **Note:** Plugin loading is enabled in the current build. Install additional plugins with their own `pip install -e .` packages (see examples below).

## What plugins receive

- **Selected files** from the Files panel (by default only FASTA; can opt into all file types).
- **`RepeatViewerProject`** — genome path, input dir, BLAST settings, `batch_workers`, etc.
- Each file’s folder is the natural place for outputs (same rule as built-in tools).

## What plugins must provide

Subclass `RepeatViewerPlugin` in `repeatviewer.plugins.api` and implement:

| Member | Purpose |
|--------|---------|
| `id` | Unique name, e.g. `my_tool` |
| `label` | Text on the toolbar button |
| `validate(files, project)` | Return `None` or an error string |
| `run(file_path, project)` | Process one file; return list of **new file paths** (strings) |

Optional:

- `tooltip` — hover help
- `fasta_only = False` — allow PNG/HTML/etc. in the selection
- `supports_parallel = False` — force one file at a time
- Override `run_batch(...)` for custom multi-file logic (e.g. one combined run)

## Register a plugin (pip package)

In your package `pyproject.toml`:

```toml
[project.entry-points."repeatviewer.plugins"]
my_tool = "my_package.repeatviewer_plugin:MyPlugin"
```

Install the package in the same environment as RepeatViewer:

```bash
pip install -e ./my-repeatviewer-plugin
```

Restart RepeatViewer. The plugin appears on the function toolbar after a separator.

Open **RepeatViewer → Setting → Toolbar** to choose visible tools and button order.

## Examples

| Example | Purpose |
|---------|---------|
| `examples/repeatviewer_sequence_stats_plugin/` | Minimal plugin (sequence stats text file) |
| `examples/repeatviewer_repeatmasker_plugin/` | **Run RepeatMasker** (`-lib` or `-species`; mask selection or defined genome) |
| `examples/repeatviewer_benchmark_plugin/` | **TEtrimmer TE benchmarking** (RepeatMasker `.out` vs reference) — see **[PLUGINS_TE_BENCHMARK.md](PLUGINS_TE_BENCHMARK.md)** |

```bash
cd examples/repeatviewer_sequence_stats_plugin
pip install -e .
# or: cd examples/repeatviewer_repeatmasker_plugin && pip install -e .
repeatviewer
```

## Result format

`run` should return created output paths. The UI will:

- Highlight new outputs (green) and inputs (blue)
- Refresh the folder tree
- Open HTML/PNG/PDF previews when appropriate
- Show failures in the status bar and a warning dialog if any file failed

## Threading

`run` / `run_batch` execute on a **background thread**. Do not touch Qt widgets from a plugin. Use files + `progress_cb` only.

## Genome / external tools

Check paths in `validate()`:

```python
def validate(self, files, project):
    if not project.genome_file:
        return "Set genome FASTA in Project settings."
    return None
```

## Versioning

Pin a compatible RepeatViewer range in your plugin’s dependencies, e.g. `repeatviewer>=0.2.0,<0.3.0`.

## Future ideas (not implemented yet)

- Plugin settings dialog hook
- Right-click menu entries
- Preview widgets for custom file types

Contributions to extend the API are welcome.
