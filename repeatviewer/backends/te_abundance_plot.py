"""Cumulative TE-family abundance plot from a RepeatMasker .out file."""

from __future__ import annotations

import html
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from repeatviewer.backends.blast_genome_plot import (
    _figure_to_preview_html,
    is_repeatmasker_out_file,
    read_repeatmasker_out,
    te_hit_table_stem,
    te_plot_output_suffix,
)


def te_abundance_plot_html_path_for_sources(sources: list[Path] | Path) -> Path:
    """Expected abundance-plot HTML path for one or more RepeatMasker .out files."""
    paths = (
        [sources.expanduser().resolve()]
        if isinstance(sources, Path)
        else [p.expanduser().resolve() for p in sources if p]
    )
    if not paths:
        raise ValueError("Select at least one RepeatMasker .out file.")
    stem = te_hit_table_stem(paths[0])
    suffix = te_plot_output_suffix(paths, kind="abund")
    return paths[0].parent / f"{stem}{suffix}"

TE_ABUNDANCE_PERCENT = "percent"
TE_ABUNDANCE_BP = "bp"
TE_ABUNDANCE_COLORS = (
    "#3b82f6",
    "#ef4444",
    "#10b981",
    "#f59e0b",
    "#8b5cf6",
    "#06b6d4",
    "#ec4899",
    "#84cc16",
    "#f97316",
    "#6366f1",
)


def _merged_interval_bp(intervals: list[tuple[int, int]]) -> int:
    """Inclusive union length for intervals from one family/chromosome/strand."""
    if not intervals:
        return 0
    ordered = sorted((min(a, b), max(a, b)) for a, b in intervals)
    total = 0
    current_start, current_end = ordered[0]
    for start, end in ordered[1:]:
        if start <= current_end:
            current_end = max(current_end, end)
        else:
            total += current_end - current_start + 1
            current_start, current_end = start, end
    return total + current_end - current_start + 1


def family_occupied_bp(hits: pd.DataFrame) -> pd.Series:
    """Occupied bp per RepeatMasker column-10 name; opposite strands stay separate."""
    if hits.empty:
        return pd.Series(dtype=np.int64)
    totals: dict[str, int] = {}
    # Use RepeatMasker's 10th column (matching repeat), never column 11
    # (class/family). ``family`` is retained as a compatibility fallback.
    name_column = "repeat_name" if "repeat_name" in hits.columns else "family"
    grouped = hits.groupby([name_column, "chrom", "sstrand"], dropna=False, sort=False)
    for (family, _chrom, _strand), group in grouped:
        intervals = list(zip(group["start"].astype(int), group["end"].astype(int)))
        name = str(family).strip() or "Unclassified"
        totals[name] = totals.get(name, 0) + _merged_interval_bp(intervals)
    return pd.Series(totals, dtype=np.int64).sort_values(ascending=False, kind="stable")


def repeatmasker_genome_size(hits: pd.DataFrame) -> int:
    """Infer total query length using RepeatMasker's end + bases-left fields."""
    if hits.empty or "query_size" not in hits.columns:
        return 0
    sizes = pd.to_numeric(hits["query_size"], errors="coerce").fillna(0)
    per_query = pd.DataFrame({"chrom": hits["chrom"].astype(str), "size": sizes}).groupby(
        "chrom", sort=False
    )["size"].max()
    return int(per_query.clip(lower=0).sum())


def _inject_page_controls(
    page: str,
    *,
    datasets: list[dict[str, Any]],
) -> str:
    sections: list[str] = []
    for index, dataset in enumerate(datasets):
        rows = "".join(
            "<tr>"
            f"<td>{rank}</td>"
            f'<td class="te-name" title="{html.escape(str(family), quote=True)}">'
            f"{html.escape(str(family))}</td>"
            f"<td>{int(bp):,}</td>"
            f"<td>{float(percent):.2f}</td>"
            "</tr>"
            for rank, (family, bp, percent) in enumerate(
                zip(dataset["families"], dataset["occupied"], dataset["family_percent"]),
                start=1,
            )
        )
        label = html.escape(str(dataset["label"]))
        color = html.escape(str(dataset["color"]), quote=True)
        sections.append(
            f'<section class="te-dataset" data-trace-index="{index}" '
            f'style="--dataset-color:{color}">'
            f'<h3><input type="color" class="te-color-picker" value="{color}" '
            f'title="Change line color for {label}" aria-label="Color for {label}">'
            f"{label}</h3>"
            '<table><colgroup><col style="width:9%"><col style="width:46%">'
            '<col style="width:25%"><col style="width:20%"></colgroup>'
            '<thead><tr><th>Rank</th><th>TE family</th><th>Occupied bp</th>'
            f'<th>Genome %</th></tr></thead><tbody>{rows}</tbody></table></section>'
        )
    controls = """
<style>
.te-abundance-export{position:fixed;right:14px;top:14px;z-index:100000;display:flex;gap:6px;
background:#fff;border:1px solid #d0d7de;border-radius:8px;padding:8px;box-shadow:0 4px 16px rgba(0,0,0,.15)}
.te-abundance-export button{font-size:12px;padding:3px 8px;cursor:pointer}
.te-abundance-list{{position:fixed;right:14px;top:70px;z-index:99999;width:440px;min-width:300px;min-height:100px;
max-width:calc(100vw - 28px);max-height:calc(100vh - 88px);resize:both;
overflow:auto;background:#fff;border:1px solid #d0d7de;border-radius:8px;padding:8px;
font:9px Arial,sans-serif;box-shadow:0 4px 16px rgba(0,0,0,.12)}}
.te-abundance-list .te-list-title{{display:flex;align-items:center;gap:8px;cursor:move;user-select:none;margin-bottom:8px}}
.te-abundance-list h2{{font-size:11px;margin:0;flex:1}}
.te-abundance-list .te-list-toggle{{font-size:9px;padding:2px 6px;cursor:pointer}}
.te-abundance-list.te-collapsed{{min-width:220px;min-height:0;height:auto!important;resize:none;overflow:hidden}}
.te-abundance-list.te-collapsed .te-list-body{{display:none}}
.te-abundance-list h3{{font-size:10px;margin:7px 0 4px;display:flex;align-items:center;gap:5px}}
.te-abundance-list .te-dataset{{border-left:3px solid var(--dataset-color);padding-left:6px;margin-bottom:12px}}
.te-abundance-list .te-color-picker{{width:22px;height:18px;padding:0;border:1px solid #ccc;border-radius:3px;cursor:pointer;flex:0 0 auto;background:none}}
.te-abundance-list table{{width:100%;border-collapse:collapse;table-layout:fixed}}
.te-abundance-list th{{position:sticky;top:0;background:#f6f8fa;text-align:right;padding:3px;border-bottom:1px solid #d0d7de;font-size:9px}}
.te-abundance-list th:nth-child(2){{text-align:left}}
.te-abundance-list td{{padding:2px 3px;border-bottom:1px solid #eef0f2;text-align:right;white-space:nowrap;vertical-align:top;font-size:9px}}
.te-abundance-list td:first-child{{width:34px}}
.te-abundance-list .te-name{{text-align:left;white-space:normal;overflow:visible;overflow-wrap:anywhere;word-break:break-word;line-height:1.25}}
@media print{{.te-abundance-export,.modebar{{display:none!important}}html,body{{margin:0!important;padding:0!important}}
.te-abundance-list{{position:absolute;box-shadow:none;max-height:none;overflow:visible}}}}
</style>
<div class="te-abundance-export">
<button type="button" id="te-sf-save-svg">Save SVG</button>
<button type="button" id="te-sf-save-pdf">Save PDF</button>
</div>
<aside class="te-abundance-list" id="te-abundance-list">
<div class="te-list-title" id="te-list-drag-handle">
<h2>Ranked TE families</h2><button type="button" class="te-list-toggle" id="te-list-toggle">Hide</button>
</div>
<div class="te-list-body">__TE_ABUNDANCE_SECTIONS__</div>
</aside>
<script>
(() => {
  const gd = document.querySelector('.plotly-graph-div');
  if (!gd) return;
  document.getElementById('te-sf-save-svg').addEventListener('click', () =>
    Plotly.downloadImage(gd, {format:'svg', filename:'te-abundance-plot'}));
  document.getElementById('te-sf-save-pdf').addEventListener('click', () => window.print());
  const list = document.getElementById('te-abundance-list');
  const handle = document.getElementById('te-list-drag-handle');
  const toggle = document.getElementById('te-list-toggle');
  toggle.addEventListener('click', (event) => {
    event.stopPropagation();
    list.classList.toggle('te-collapsed');
    toggle.textContent = list.classList.contains('te-collapsed') ? 'Show' : 'Hide';
  });
  let dragging = false;
  let dragX = 0;
  let dragY = 0;
  handle.addEventListener('mousedown', (event) => {
    if (event.target === toggle) return;
    dragging = true;
    const rect = list.getBoundingClientRect();
    list.style.left = rect.left + 'px';
    list.style.top = rect.top + 'px';
    list.style.right = 'auto';
    dragX = event.clientX - rect.left;
    dragY = event.clientY - rect.top;
    event.preventDefault();
  });
  window.addEventListener('mousemove', (event) => {
    if (!dragging) return;
    list.style.left = Math.max(4, Math.min(window.innerWidth - 80, event.clientX - dragX)) + 'px';
    list.style.top = Math.max(4, Math.min(window.innerHeight - 40, event.clientY - dragY)) + 'px';
  });
  window.addEventListener('mouseup', () => { dragging = false; });
  list.querySelectorAll('.te-color-picker').forEach((input) => {
    input.addEventListener('input', () => {
      const section = input.closest('.te-dataset');
      if (!section) return;
      const idx = Number(section.getAttribute('data-trace-index'));
      if (!Number.isInteger(idx) || idx < 0) return;
      const color = input.value;
      section.style.setProperty('--dataset-color', color);
      Plotly.restyle(gd, {
        'line.color': color,
        'marker.color': color,
      }, [idx]);
    });
  });
})();
</script>
"""
    controls = (
        controls.replace("{{", "{")
        .replace("}}", "}")
        .replace("__TE_ABUNDANCE_SECTIONS__", "".join(sections))
    )
    return page.replace("</body>", controls + "\n</body>") if "</body>" in page else page + controls


def write_te_abundance_plot_html(
    repeatmasker_out: Path | list[Path],
    *,
    output_html: Path | None = None,
    metric: str = TE_ABUNDANCE_PERCENT,
    top_n: int = 100,
) -> Path:
    """Write one or more independently ranked cumulative family-abundance curves."""
    sources = [repeatmasker_out] if isinstance(repeatmasker_out, Path) else repeatmasker_out
    sources = [source.expanduser().resolve() for source in sources]
    if not sources:
        raise ValueError("Select at least one RepeatMasker .out file.")
    datasets: list[dict[str, Any]] = []
    for index, source in enumerate(sources):
        if not is_repeatmasker_out_file(source):
            raise ValueError(f"TE abundance plot requires RepeatMasker .out files: {source}")
        hits = read_repeatmasker_out(source)
        abundance = family_occupied_bp(hits).iloc[: max(1, int(top_n))]
        if abundance.empty:
            raise ValueError(f"No RepeatMasker TE intervals were found in {source.name}.")
        families = abundance.index.astype(str).tolist()
        occupied = abundance.to_numpy(dtype=np.int64)
        cumulative_bp = np.cumsum(occupied)
        genome_size = repeatmasker_genome_size(hits)
        if metric == TE_ABUNDANCE_PERCENT and genome_size <= 0:
            raise ValueError(f"Could not infer genome length from {source.name}.")
        family_percent = (
            occupied.astype(float) / float(genome_size) * 100.0
            if genome_size > 0
            else np.zeros(len(occupied), dtype=float)
        )
        y = (
            cumulative_bp.astype(float) / float(genome_size) * 100.0
            if metric == TE_ABUNDANCE_PERCENT
            else cumulative_bp
        )
        datasets.append(
            {
                "source": source,
                "label": source.stem,
                "color": TE_ABUNDANCE_COLORS[index % len(TE_ABUNDANCE_COLORS)],
                "families": families,
                "occupied": occupied,
                "family_percent": family_percent,
                "y": y,
                "genome_size": genome_size,
            }
        )

    y_title = (
        "Cumulative occupied genome (%)"
        if metric == TE_ABUNDANCE_PERCENT
        else "Cumulative occupied nucleotides (bp)"
    )
    value_suffix = "%" if metric == TE_ABUNDANCE_PERCENT else " bp"
    fig = go.Figure()
    for dataset in datasets:
        ranks = np.arange(1, len(dataset["families"]) + 1)
        customdata = np.column_stack(
            [
                np.asarray(dataset["families"], dtype=object),
                dataset["occupied"],
                dataset["family_percent"],
            ]
        )
        fig.add_trace(
            go.Scatter(
                x=ranks,
                y=dataset["y"],
                mode="lines+markers",
                name=dataset["label"],
                line={"color": dataset["color"], "width": 2},
                marker={"color": dataset["color"], "size": 7},
                customdata=customdata,
                hovertemplate=(
                    f"<b>{html.escape(str(dataset['label']))}</b><br>"
                    "TE family: %{customdata[0]}<br>Rank: %{x}"
                    "<br>Family occupied: %{customdata[1]:,} bp"
                    "<br>Family genome: %{customdata[2]:.2f}%"
                    f"<br>Cumulative: %{{y:,.4g}}{value_suffix}<extra></extra>"
                ),
                showlegend=len(datasets) > 1,
            )
        )
    fig.update_layout(
        title=(
            f"TE family abundance — {sources[0].name}"
            if len(sources) == 1
            else f"TE family abundance comparison — {len(sources)} RepeatMasker outputs"
        ),
        xaxis={
            "title": "Accumulated TE family number (abundance rank)",
            "showline": True,
            "linecolor": "#000000",
            "linewidth": 1,
            "ticks": "outside",
            "tickcolor": "#000000",
            "showgrid": False,
            "zeroline": False,
        },
        yaxis={
            "title": y_title,
            "showline": True,
            "linecolor": "#000000",
            "linewidth": 1,
            "ticks": "outside",
            "tickcolor": "#000000",
            "showgrid": False,
            "zeroline": False,
        },
        template="plotly_white",
        plot_bgcolor="#ffffff",
        paper_bgcolor="#ffffff",
        hovermode="closest",
        width=1200,
        height=700,
        margin={"l": 90, "r": 500, "t": 85, "b": 75},
        showlegend=len(datasets) > 1,
        legend={"title": "RepeatMasker outputs", "orientation": "h", "y": 1.08},
        meta={
            "sources": [str(source) for source in sources],
            "genomeSizes": [dataset["genome_size"] for dataset in datasets],
            "topN": top_n,
            "metric": metric,
            "families": json.dumps([dataset["families"] for dataset in datasets]),
        },
    )
    destination = (
        output_html.expanduser().resolve()
        if output_html is not None
        else te_abundance_plot_html_path_for_sources(sources)
    )
    page = _inject_page_controls(
        _figure_to_preview_html(fig, destination),
        datasets=datasets,
    )
    destination.write_text(page, encoding="utf-8")
    return destination
