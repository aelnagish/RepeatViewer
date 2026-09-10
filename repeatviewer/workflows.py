"""User-defined tool workflows (chain tools by compatible file formats)."""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path

from PySide6.QtCore import QSettings

# Logical format tags used for chain compatibility (not file suffixes).
FMT_FASTA = "fasta"
FMT_HTML = "html"
FMT_PLOT = "plot"

FASTA_SUFFIXES = {".fa", ".fasta", ".fas", ".faa", ".fna"}

_GROUP = "workflows"
_KEY = "workflows_json"
_KEY_DEFAULTS_SEEDED = "defaults_seeded"


@dataclass(frozen=True)
class WorkflowToolInfo:
    """Metadata for one chainable tool step."""

    tool_id: str
    display_name: str
    accepts: frozenset[str]
    produces: frozenset[str]
    needs_genome: bool = False


# v1 chainable tools. Terminal producers (html/plot only) may only be last step.
WORKFLOW_TOOLS: dict[str, WorkflowToolInfo] = {
    "blast": WorkflowToolInfo(
        "blast",
        "BLAST",
        frozenset({FMT_FASTA}),
        frozenset({FMT_FASTA}),
        needs_genome=True,
    ),
    "extend": WorkflowToolInfo(
        "extend",
        "MSA Extend",
        frozenset({FMT_FASTA}),
        frozenset({FMT_FASTA}),
        needs_genome=True,
    ),
    "align": WorkflowToolInfo(
        "align",
        "Align (MAFFT)",
        frozenset({FMT_FASTA}),
        frozenset({FMT_FASTA}),
    ),
    "msa_clean": WorkflowToolInfo(
        "msa_clean",
        "Combined MSA clean",
        frozenset({FMT_FASTA}),
        frozenset({FMT_FASTA}),
    ),
    "msa_clean_gappy": WorkflowToolInfo(
        "msa_clean_gappy",
        "Remove gappy columns",
        frozenset({FMT_FASTA}),
        frozenset({FMT_FASTA}),
    ),
    "msa_clean_div": WorkflowToolInfo(
        "msa_clean_div",
        "Clean rows with divergence",
        frozenset({FMT_FASTA}),
        frozenset({FMT_FASTA}),
    ),
    "msa_clean_gap": WorkflowToolInfo(
        "msa_clean_gap",
        "Clean rows with gaps",
        frozenset({FMT_FASTA}),
        frozenset({FMT_FASTA}),
    ),
    "consensus_generation": WorkflowToolInfo(
        "consensus_generation",
        "Consensus",
        frozenset({FMT_FASTA}),
        frozenset({FMT_FASTA}),
    ),
    "teaid": WorkflowToolInfo(
        "teaid",
        "TEAid",
        frozenset({FMT_FASTA}),
        frozenset({FMT_HTML}),
        needs_genome=True,
    ),
    "flexidot": WorkflowToolInfo(
        "flexidot",
        "Flexidot",
        frozenset({FMT_FASTA}),
        frozenset({FMT_PLOT}),
        needs_genome=True,
    ),
    "orf_prediction": WorkflowToolInfo(
        "orf_prediction",
        "ORF Prediction",
        frozenset({FMT_FASTA}),
        frozenset({FMT_FASTA}),
    ),
    "reverse_complement": WorkflowToolInfo(
        "reverse_complement",
        "Reverse complement",
        frozenset({FMT_FASTA}),
        frozenset({FMT_FASTA}),
    ),
}

WORKFLOW_TOOL_ORDER: tuple[str, ...] = tuple(WORKFLOW_TOOLS.keys())

# Menu categories for the designer (matches Sequence / Tools / Plugins).
WORKFLOW_CATEGORY_SEQUENCES: tuple[str, ...] = (
    "reverse_complement",
    "orf_prediction",
)
WORKFLOW_CATEGORY_TOOLS: tuple[str, ...] = (
    "blast",
    "align",
    "extend",
    "msa_clean",
    "msa_clean_gappy",
    "msa_clean_div",
    "msa_clean_gap",
    "consensus_generation",
    "teaid",
    "flexidot",
)


def is_plugin_tool_id(tool_id: str) -> bool:
    return str(tool_id).startswith("plugin:")


def plugin_id_from_tool_id(tool_id: str) -> str:
    return str(tool_id)[len("plugin:") :]


def resolve_tool_info(tool_id: str) -> WorkflowToolInfo | None:
    """Built-in tool metadata, or a FASTA→FASTA stub for ``plugin:{id}``."""
    if tool_id in WORKFLOW_TOOLS:
        return WORKFLOW_TOOLS[tool_id]
    if is_plugin_tool_id(tool_id):
        from repeatviewer.plugins.registry import get_registry

        plugin = get_registry().get(plugin_id_from_tool_id(tool_id))
        if plugin is None:
            return None
        return WorkflowToolInfo(
            tool_id,
            getattr(plugin, "label", None) or plugin_id_from_tool_id(tool_id),
            frozenset({FMT_FASTA}),
            frozenset({FMT_FASTA}),
            needs_genome=False,
        )
    return None


def list_plugin_workflow_tools() -> list[WorkflowToolInfo]:
    from repeatviewer.plugins.registry import get_registry

    tools: list[WorkflowToolInfo] = []
    for plugin in get_registry().all_plugins():
        tool_id = f"plugin:{plugin.id}"
        info = resolve_tool_info(tool_id)
        if info is not None:
            tools.append(info)
    return tools


def workflow_id_from_name(name: str) -> str:
    text = re.sub(r"[^a-zA-Z0-9]+", "_", (name or "").strip().lower())
    text = text.strip("_")
    return text or "workflow"


def tools_compatible(prev_id: str, next_id: str) -> bool:
    """True when previous tool output formats intersect next tool accepted formats."""
    prev = resolve_tool_info(prev_id)
    nxt = resolve_tool_info(next_id)
    if prev is None or nxt is None:
        return False
    return bool(prev.produces & nxt.accepts)


def validate_workflow_steps(tool_ids: list[str]) -> str | None:
    """Return an error message if the chain is invalid, else None."""
    if not tool_ids:
        return "Add at least one step."
    for tool_id in tool_ids:
        if resolve_tool_info(tool_id) is None:
            return f"Unknown tool: {tool_id}"
    for i in range(len(tool_ids) - 1):
        if not tools_compatible(tool_ids[i], tool_ids[i + 1]):
            a = resolve_tool_info(tool_ids[i]).display_name  # type: ignore[union-attr]
            b = resolve_tool_info(tool_ids[i + 1]).display_name  # type: ignore[union-attr]
            return f"Cannot chain {a} → {b}: output format does not match the next input."
    return None


def allowed_next_tools(after_tool_id: str | None) -> list[WorkflowToolInfo]:
    """Tools that may be appended after ``after_tool_id`` (or any if None)."""
    candidates: list[WorkflowToolInfo] = [
        WORKFLOW_TOOLS[tid] for tid in WORKFLOW_TOOL_ORDER
    ]
    try:
        candidates.extend(list_plugin_workflow_tools())
    except Exception:  # noqa: BLE001 — registry may be unavailable in tests
        pass
    if after_tool_id is None:
        return candidates
    return [info for info in candidates if tools_compatible(after_tool_id, info.tool_id)]


def allowed_next_in_category(
    after_tool_id: str | None,
    category_tool_ids: tuple[str, ...] | list[str],
) -> list[WorkflowToolInfo]:
    out: list[WorkflowToolInfo] = []
    for tool_id in category_tool_ids:
        info = WORKFLOW_TOOLS.get(tool_id)
        if info is None:
            continue
        if after_tool_id is None or tools_compatible(after_tool_id, tool_id):
            out.append(info)
    return out


def filter_chainable_fasta(
    paths: list[str | Path],
    *,
    exclude: set[Path] | None = None,
) -> list[Path]:
    """Keep existing FASTA outputs suitable as the next step's inputs.

    Skips genome reference links (``*.lnk.fa``) and any paths in ``exclude``
    (typically configured genome FASTA files). Those are BLAST sidecars, not
    hit FASTAs, and break MSA Extend / Align if chained.
    """
    excluded: set[str] = set()
    for raw in exclude or ():
        try:
            excluded.add(str(Path(raw).expanduser().resolve()))
        except OSError:
            continue

    out: list[Path] = []
    seen: set[str] = set()
    for raw in paths:
        path = Path(str(raw)).expanduser()
        name_lower = path.name.lower()
        if name_lower.endswith(".lnk.fa") or name_lower.endswith(".lnk.fasta"):
            continue
        if path.suffix.lower() not in FASTA_SUFFIXES:
            continue
        if not path.is_file() and not path.is_symlink():
            continue
        try:
            resolved = path.resolve()
        except OSError:
            continue
        key = str(resolved)
        if key in seen or key in excluded:
            continue
        seen.add(key)
        # Prefer the on-disk path (symlink) when present; else resolved target.
        out.append(path if path.exists() else resolved)
    return out


@dataclass
class WorkflowStep:
    tool_id: str

    def to_json(self) -> dict:
        return {"tool_id": self.tool_id}

    @staticmethod
    def from_json(raw: dict) -> "WorkflowStep":
        return WorkflowStep(tool_id=str(raw.get("tool_id", "")).strip())


@dataclass
class WorkflowSpec:
    id: str
    name: str
    steps: list[WorkflowStep] = field(default_factory=list)

    def tool_ids(self) -> list[str]:
        return [step.tool_id for step in self.steps]

    def to_json(self) -> dict:
        return {
            "id": self.id,
            "name": self.name,
            "steps": [step.to_json() for step in self.steps],
        }

    @staticmethod
    def from_json(raw: dict) -> "WorkflowSpec":
        name = str(raw.get("name", "")).strip() or "Workflow"
        workflow_id = str(raw.get("id", "")).strip() or workflow_id_from_name(name)
        steps: list[WorkflowStep] = []
        for item in raw.get("steps", []) or []:
            if not isinstance(item, dict):
                continue
            step = WorkflowStep.from_json(item)
            if step.tool_id in WORKFLOW_TOOLS or is_plugin_tool_id(step.tool_id):
                if resolve_tool_info(step.tool_id) is not None:
                    steps.append(step)
        return WorkflowSpec(id=workflow_id, name=name, steps=steps)

    def needs_genome(self) -> bool:
        return any(
            (info := resolve_tool_info(s.tool_id)) is not None and info.needs_genome
            for s in self.steps
        )


def default_bee_workflow() -> WorkflowSpec:
    """Built-in BEE pipeline: BLAST → MSA Extend → Align (MAFFT)."""
    return WorkflowSpec(
        id="bee",
        name="BEE",
        steps=[
            WorkflowStep(tool_id="blast"),
            WorkflowStep(tool_id="extend"),
            WorkflowStep(tool_id="align"),
        ],
    )


def ensure_default_workflows() -> list[WorkflowSpec]:
    """Seed the default BEE workflow once when none are stored yet."""
    specs = load_workflow_specs()
    if specs:
        return specs
    s = QSettings()
    s.beginGroup(_GROUP)
    already_seeded = bool(s.value(_KEY_DEFAULTS_SEEDED, False))
    if already_seeded:
        s.endGroup()
        return []
    s.setValue(_KEY_DEFAULTS_SEEDED, True)
    s.endGroup()
    bee = default_bee_workflow()
    _write_all([bee])
    return [bee]


def load_workflow_specs() -> list[WorkflowSpec]:
    s = QSettings()
    s.beginGroup(_GROUP)
    raw = s.value(_KEY, "")
    s.endGroup()
    if not raw:
        return []
    try:
        data = json.loads(str(raw))
    except (TypeError, ValueError, json.JSONDecodeError):
        return []
    if not isinstance(data, list):
        return []
    specs: list[WorkflowSpec] = []
    for item in data:
        if not isinstance(item, dict):
            continue
        try:
            specs.append(WorkflowSpec.from_json(item))
        except Exception:  # noqa: BLE001
            continue
    return specs


def get_workflow_spec(workflow_id: str) -> WorkflowSpec | None:
    for spec in load_workflow_specs():
        if spec.id == workflow_id:
            return spec
    return None


def _write_all(specs: list[WorkflowSpec]) -> None:
    s = QSettings()
    s.beginGroup(_GROUP)
    s.setValue(_KEY, json.dumps([item.to_json() for item in specs]))
    s.endGroup()


def save_workflow_order(workflow_ids: list[str]) -> None:
    """Persist the designer list order (menu order follows this)."""
    by_id = {spec.id: spec for spec in load_workflow_specs()}
    ordered: list[WorkflowSpec] = []
    seen: set[str] = set()
    for workflow_id in workflow_ids:
        spec = by_id.get(workflow_id)
        if spec is None or workflow_id in seen:
            continue
        ordered.append(spec)
        seen.add(workflow_id)
    for spec in load_workflow_specs():
        if spec.id not in seen:
            ordered.append(spec)
    _write_all(ordered)


def save_workflow_spec(spec: WorkflowSpec) -> None:
    error = validate_workflow_steps(spec.tool_ids())
    if error:
        raise ValueError(error)
    if not spec.name.strip():
        raise ValueError("Workflow name is required.")
    specs = load_workflow_specs()
    updated = False
    for index, existing in enumerate(specs):
        if existing.id == spec.id:
            specs[index] = spec
            updated = True
            break
    if not updated:
        # Avoid colliding ids when renaming display name only.
        existing_ids = {item.id for item in specs}
        if spec.id in existing_ids:
            base = spec.id
            n = 2
            while f"{base}_{n}" in existing_ids:
                n += 1
            spec = WorkflowSpec(id=f"{base}_{n}", name=spec.name, steps=list(spec.steps))
        specs.append(spec)
    _write_all(specs)


def delete_workflow_spec(workflow_id: str) -> bool:
    specs = load_workflow_specs()
    kept = [item for item in specs if item.id != workflow_id]
    if len(kept) == len(specs):
        return False
    _write_all(kept)
    return True


def duplicate_workflow_spec(workflow_id: str) -> WorkflowSpec | None:
    source = get_workflow_spec(workflow_id)
    if source is None:
        return None
    existing_ids = {item.id for item in load_workflow_specs()}
    base_name = f"{source.name} (copy)"
    new_id = workflow_id_from_name(base_name)
    n = 2
    candidate = new_id
    while candidate in existing_ids:
        candidate = f"{new_id}_{n}"
        n += 1
    copy = WorkflowSpec(
        id=candidate,
        name=base_name,
        steps=[WorkflowStep(tool_id=step.tool_id) for step in source.steps],
    )
    save_workflow_spec(copy)
    return copy


def run_workflow(
    spec: WorkflowSpec,
    seq_files: list[Path],
    project,
    progress_cb,
    *,
    cancel=None,
    max_workers: int = 4,
    blast_options=None,
    extend=None,
    teaid_options=None,
    msa_clean_options=None,
    consensus_options=None,
    flexidot_options=None,
    orf_prediction_settings=None,
    align_options=None,
    use_system_mafft: bool = False,
    use_system_muscle: bool = False,
    align_engine: str = "mafft",
) -> dict:
    """Run workflow steps sequentially; feed FASTA outputs into the next step."""
    from repeatviewer.backends.curation import run_batch
    from repeatviewer.backends.job_cancel import check_cancelled
    from repeatviewer.backends.sequence_aligner import run_alignment
    from repeatviewer.job_options import ALIGNER_MAFFT

    error = validate_workflow_steps(spec.tool_ids())
    if error:
        raise ValueError(error)
    if not seq_files:
        raise ValueError("No FASTA sequences selected.")

    n_steps = len(spec.steps)
    inputs = [Path(p).expanduser().resolve() for p in seq_files]
    created_outputs: list[str] = []
    source_files: list[str] = [str(p) for p in inputs]
    failures: list[str] = []
    genome_exclude: set[Path] = set()
    try:
        for _label, genome_path in project.genome_targets():
            genome_exclude.add(Path(genome_path).expanduser().resolve())
    except Exception:  # noqa: BLE001
        pass
    try:
        for entry in getattr(project, "pan_genomes", []) or []:
            genome_exclude.add(Path(entry.path).expanduser().resolve())
    except Exception:  # noqa: BLE001
        pass
    if getattr(project, "genome_file", None):
        try:
            genome_exclude.add(Path(project.genome_file).expanduser().resolve())
        except Exception:  # noqa: BLE001
            pass

    def _step_progress(completed: int, label: str, *, step_index: int, tool_name: str) -> None:
        # Report overall step index; detail includes tool + file.
        progress_cb(step_index, f"{spec.name}: {step_index + 1}/{n_steps} {tool_name} — {label}")

    for step_index, step in enumerate(spec.steps):
        check_cancelled(cancel)
        info = resolve_tool_info(step.tool_id)
        if info is None:
            raise RuntimeError(f"Unknown tool: {step.tool_id}")
        tool_name = info.display_name
        progress_cb(step_index, f"{spec.name}: {step_index + 1}/{n_steps} {tool_name}…")

        step_outs: list[str] = []
        if is_plugin_tool_id(step.tool_id):
            from repeatviewer.plugins.registry import get_registry

            plugin = get_registry().get(plugin_id_from_tool_id(step.tool_id))
            if plugin is None:
                raise RuntimeError(f"Plugin not found: {step.tool_id}")
            runtime_values = {}
            spec_obj = getattr(plugin, "spec", None)
            if spec_obj is not None:
                runtime_values = {
                    d.name: d.value for d in getattr(spec_obj, "var_defs", [])
                }
            for path in inputs:
                check_cancelled(cancel)
                _step_progress(0, path.name, step_index=step_index, tool_name=tool_name)
                try:
                    outs = plugin.run(
                        path, project, runtime_values=runtime_values
                    )
                    step_outs.extend(str(p) for p in (outs or []))
                except Exception as exc:  # noqa: BLE001
                    failures.append(f"{tool_name}/{path.name}: {exc}")
        elif step.tool_id == "align":
            engine = align_engine or ALIGNER_MAFFT
            options = align_options
            if options is None:
                from repeatviewer.job_options import AlignOptions

                options = AlignOptions()
            align_kwargs = {
                "num_threads": options.num_threads,
                "strategy": options.strategy,
                "sequence_type": options.sequence_type,
                "retree": options.retree,
                "max_iterate": options.max_iterate,
                "gap_open_penalty": options.gap_open_penalty,
                "offset": options.offset,
                "reorder": options.reorder,
                "adjust_direction": options.adjust_direction,
                "muscle_mode": options.muscle_mode,
                "muscle_num_threads": options.muscle_num_threads,
            }
            for path in inputs:
                check_cancelled(cancel)
                _step_progress(0, path.name, step_index=step_index, tool_name=tool_name)
                try:
                    output = run_alignment(
                        path,
                        **align_kwargs,
                        use_system_mafft=use_system_mafft,
                        use_system_muscle=use_system_muscle,
                        engine=engine,
                    )
                    step_outs.append(str(output))
                except Exception as exc:  # noqa: BLE001
                    failures.append(f"{tool_name}/{path.name}: {exc}")
        else:
            batch_kwargs: dict = {
                "max_workers": max_workers,
                "cancel": cancel,
            }
            if blast_options is not None:
                batch_kwargs["blast_options"] = blast_options
            if extend is not None:
                batch_kwargs["extend"] = extend
            if teaid_options is not None:
                batch_kwargs["teaid_options"] = teaid_options
            if msa_clean_options is not None:
                batch_kwargs["msa_clean_options"] = msa_clean_options
            if consensus_options is not None:
                batch_kwargs["consensus_options"] = consensus_options
            if flexidot_options is not None:
                batch_kwargs["flexidot_options"] = flexidot_options
            if orf_prediction_settings is not None:
                batch_kwargs["orf_prediction_settings"] = orf_prediction_settings

            def _nested(completed: int, label: str, *, idx=step_index, name=tool_name) -> None:
                _step_progress(completed, label, step_index=idx, tool_name=name)

            result = run_batch(
                step.tool_id,
                inputs,
                project,
                _nested,
                **batch_kwargs,
            )
            step_outs.extend(result.get("created_outputs") or [])
            for fail in result.get("failures") or []:
                failures.append(f"{tool_name}: {fail}")

        created_outputs.extend(step_outs)

        if step_index < n_steps - 1:
            next_inputs = filter_chainable_fasta(step_outs, exclude=genome_exclude)
            if not next_inputs:
                raise RuntimeError(
                    f"{tool_name} produced no FASTA output to feed the next step."
                )
            inputs = next_inputs

    progress_cb(n_steps, f"{spec.name}: finished")
    return {
        "func_name": f"workflow:{spec.id}",
        "total": n_steps,
        "failures": failures,
        "created_outputs": created_outputs,
        "source_files": source_files,
        "workflow_name": spec.name,
    }
