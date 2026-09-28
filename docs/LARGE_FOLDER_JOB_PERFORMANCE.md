# UI performance: large folders, Annotation, and Overview

## 1. Large-folder jobs (file list)

### Why this exists

BLAST (or another tool) on **one** sequence can finish quickly, but the
application still has to update the file list and open the result. In a folder
with hundreds or thousands of files (for example ~1700), that UI work used to
touch the **whole** list — rebuild rows, restyle highlights, wait on metadata —
so **starting** and **finishing** felt slow even when the tool itself was fast.

The delay was the **file-list UI**, not BLAST.

### Rule (short)

When a folder is large, a small job must not pay costs proportional to the
folder size. Work should scale with **how many outputs** the job produced, not
with **how many files are already in the folder**.

**Quick check:** one-file BLAST in a ~1700-file folder should start quickly,
run normally, and show the result without a long freeze.

### What went wrong before

- Clearing or applying job highlights restyled **every** row (very expensive
  at ~1700 files; path-index rebuilds made it worse).
- Stopping background metadata enrichment could **block the GUI** for up to
  about a second when a job started (`thread.wait(1000)`).
- After finish, the app often **rescanned / rebuilt** the full folder and
  opened preview more than once.
- Optimizations that only kicked in for **many outputs** missed the common
  case: **one** output inside a **huge** folder.

### Rules to follow

**Thresholds** — see `_LARGE_FILE_LIST` and `_LARGE_JOB_OUTPUT_BATCH` in
`repeatviewer/main_window.py`.

- **Folder size** drives list rebuild, metadata, and highlight cost.
- **Output batch size** drives preview picking and input tinting.
- A **single** output in a huge folder still needs the **folder-size** fast
  path. Batch-size gates alone are not enough.

**Never do on the GUI thread for large folders**

- Full file-list rebuild / `scandir` when only a few known outputs were added.
- Restyle every row to clear or apply job highlights.
- Rebuild the path index once per row (O(n²)) — tint by **key set** instead
  (`reapply_file_colors_for_keys`, `_apply_job_highlights_to_keys`).
- `N × path.resolve()` or `N × is_file()` over the whole folder or all outputs.
- Sync FASTA metadata/stat while starting or finishing a job.
- `_stop_metadata_enrichment()` with a blocking wait on job start — use
  `wait_ms=0`.
- Open preview more than once (e.g. populate with `open_preview=True` **and**
  `_auto_open_outputs`).
- Connect worker `finished` / `failed` to bare lambdas or closures (can run
  off the GUI thread and crash) — use `@Slot` methods and `QueuedConnection`.

**Required patterns**

- **Job start:** clear only previous highlight keys; stop enrichment with
  `wait_ms=0`.
- **Job finish refresh:** never auto-open inside the refresh — the caller
  opens once.
- **Large same-folder + few new outputs:** incremental add; skip
  scandir/rebuild.
- **Progress:** emit the **completed** count when a work item finishes;
  reserve **+1** step until results are shown on screen.
- **One-file jobs:** status like `Running BLAST…` with **no** filename.

### Checklist (jobs / file list)

1. One-file BLAST in a ~1700-file folder: start feels instant; finish does not
   rebuild the whole list.
2. Only new or affected rows are tinted; highlights are cleared by key, not by
   walking every row.
3. Exactly one preview open after the job.
4. No new O(n) resolve / stat / `is_file` over the folder on the GUI thread for
   small jobs.

### Main code locations

- `repeatviewer/main_window.py` — job start/finish, populate/refresh, progress
- `repeatviewer/widgets/file_list_table.py` — keyed color reapply, bulk add
- `repeatviewer/backends/curation.py` — worker progress step counts

---

## 2. Annotation tab and Overview tab (always keep smooth)

### Why this exists

Annotation (single-sequence map) and Overview (multi-FASTA bar list) can show
very long sequences or many records. If they parse, layout, or paint **more
than the visible region** on every interaction, the tabs freeze. Opening Text
must not force a full Annotation parse either.

### Rule (short)

Annotation and Overview must stay interactive: **lazy load**, **viewport-only
work**, and **no full-sequence paint or layout on every drag/zoom**. Prefer
deferred updates over blocking the GUI thread.

**Quick check:** open a large FASTA on Text first (Annotation should not parse
yet); open Annotation/Overview and pan/zoom/select without stuttering; switch
to another file and Annotation/Dotplot must not keep stale data from the
previous sequence.

### Rules to follow

**Lazy load**

- Do **not** load Annotation/Overview data when the file is only opened on the
  Text tab. Load once on first Annotation open (`_ensure_annotation_data_loaded`
  in `preview_tabs.py`).
- Create panels with `autoload=False` where applicable; call `load_path` when
  the tab is actually shown.
- When switching files, clear stale Annotation loaded data and reset Dotplot
  path so the previous long sequence cannot be re-synced into a small file.

**Viewport-only rendering (Annotation)**

- Draw nucleotide/amino-acid **letters only for the visible viewport**, not the
  whole sequence (see comments around `_MAX_ANNOTATION_DISPLAY_LEN` /
  `shows_letters()` in `sequence_map_viewer.py`).
- Batch bar/feature painting for the visible range; avoid per-base widget
  creation or thousands of individual draw calls outside the clip.
- Respect the hard display limit (`_MAX_ANNOTATION_DISPLAY_LEN`, currently
  100 Mbp). Do not silently remove that gate without a measured plan.

**Keep interaction smooth**

- During drag-select / pan / resize: **paint immediately**, **defer** spin box,
  status, and heavy chrome updates (timers / queued emits).
- Prefer damage/update of changed regions over full-widget `repaint()` of the
  entire map.
- Fit-to-window and layout sync should be deferred (`QTimer.singleShot` /
  `_schedule_fit_window`) so they do not fight Go-to or user scroll.

**Overview (multi-FASTA)**

- Overview is length/feature bars for many records — keep zoom and scroll tied
  to the **track viewport**, not a full pixel width of every base of every
  sequence.
- Prefer streaming / header-based length scans (`fasta_record_lengths`) over
  loading every full sequence into memory for the overview list.
- Filter, zoom, and extract should not rebuild the entire overview from scratch
  on every keystroke without debounce where appropriate.

**Never do**

- Parse GFF + full sequence map on every file-list click “just in case”.
- Paint letters for the entire multi-Mb sequence at letter zoom in one pass
  outside the viewport.
- Sync Dotplot’s active sequence from Annotation data that belongs to a
  **different** file path.
- Block the GUI with a synchronous full-library read when opening Overview if
  a streaming length pass is enough.

### Checklist (Annotation / Overview)

1. Opening a large FASTA on **Text** does not load Annotation/Overview yet.
2. First open of Annotation/Overview shows progress if the load is heavy, then
   remains usable for pan/zoom/select.
3. Dragging a selection does not stall the spin boxes / status every pixel.
4. Switching from a long sequence to a short one clears stale Annotation data;
   Dotplot does not keep the “too long” / old-sequence state.
5. Overview of a multi-FASTA library scrolls and zooms without allocating a
   full-base-resolution canvas for every record.

### Main code locations

- `repeatviewer/widgets/preview_tabs.py` — lazy Annotation/Overview load, path
  switch cleanup, Dotplot sync guards
- `repeatviewer/widgets/sequence_map_viewer.py` — Annotation map, viewport
  letters, display limit, deferred selection/fit
- `repeatviewer/widgets/multi_sequence_overview.py` — Overview bars / zoom /
  viewport
- `repeatviewer/backends/sequence_features.py` — `fasta_record_lengths` and
  related loaders

---

## 3. Annotation Save vs Export (JSON first)

### Why

Save used to write both `{stem}.seqmap_RV.json` and `{stem}.RV.gff3`. GFF3 is
visible in the file list and is only needed for interchange. The canonical
edit store is JSON.

### Rules

- **Save** writes **JSON only** (`.seqmap_RV.json`).
- **Export annotation…** (FASTA right-click) writes **GFF3 / GTF / BED** on
  demand.
- JSON may hold annotations for **multiple sequences**, including **duplicate
  FASTA headers**, via `format_version: 2` blocks with `record_index` (0-based
  order in the file). In memory, duplicate names become `name`, `name#2`, …
- On export, if original headers collide, warn and offer suggested unique
  seqids (`name`, `name_2`, …) the user can edit before writing.
- Rename / move / delete a FASTA also applies to its JSON (and GFF if present).
- Combine / extract / separate offer a tick/question to also process
  annotations.
