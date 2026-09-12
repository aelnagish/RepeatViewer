"""User-defined command plugins configured inside RepeatViewer."""

from __future__ import annotations

import json
import re
import shlex
import subprocess
from repeatviewer.external_process import run_logged_process
from dataclasses import dataclass
from pathlib import Path

from PySide6.QtCore import QEvent, QEventLoop, QSettings, Qt, QTimer, Signal
from PySide6.QtWidgets import (
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.widgets.spin_boxes import (
    NoWheelDoubleSpinBox as QDoubleSpinBox,
    NoWheelSpinBox as QSpinBox,
)


from repeatviewer.plugins.api import RepeatViewerPlugin
from repeatviewer.project import RepeatViewerProject

_GROUP = "plugins"
_KEY = "custom_commands_json"
# ${name} (user vars) and {RV_*} / ${RV_*} (built-in specials).
_VAR_PATTERN = re.compile(
    r"\$\{([A-Za-z_][A-Za-z0-9_]*)\}|\{(RV_[A-Za-z0-9_]+)\}"
)
# Replacement may keep a fixed suffix glued to the variable, e.g. ${RV_STEM}_b.txt
_VAR_REPLACE_PATTERN = re.compile(
    r"\$\{([A-Za-z_][A-Za-z0-9_]*)\}([A-Za-z0-9._+-]*)|"
    r"\{(RV_[A-Za-z0-9_]+)\}([A-Za-z0-9._+-]*)"
)
# Stored name for `${RV_INPUT}_b.text` as its own parameter.
_AFFIX_SEP = "|"
_AFFIX_DISPLAY = re.compile(
    r"\$\{([A-Za-z_][A-Za-z0-9_]*)\}([A-Za-z0-9._+-]+)"
)

# Built-in placeholders filled by RepeatViewer at run time (no Define-variables row).
SPECIAL_VARIABLES: dict[str, str] = {
    "RV_INPUT": (
        "Absolute path of the first selected file in the Files panel.\n"
        "Example: if you select contig.fa → /path/to/contig.fa"
    ),
    "RV_INPUTS": (
        "Absolute paths of all selected files, separated by spaces.\n"
        "Useful when a command accepts multiple inputs."
    ),
    "RV_GENOME": (
        "Absolute path of the genome FASTA from File → Define genome.\n"
        "Empty if no genome is defined."
    ),
    "RV_FOLDER": (
        "Absolute path of the folder currently shown in the Files panel."
    ),
    "RV_OUTDIR": (
        "Working directory for the run: parent folder of RV_INPUT if a file "
        "is selected, otherwise the current Files-panel folder."
    ),
    "RV_STEM": (
        "Name of RV_INPUT without the final file extension.\n"
        "Example: /data/TE123.fasta → TE123\n"
        "Use this to build output names like ${RV_OUTDIR}/${RV_STEM}_B.txt\n"
        "(Prefer ${RV_STEM}_B.txt over ${RV_INPUT}_B.txt, which keeps the .fasta suffix.)"
    ),
    "RV_NAME": (
        "File name of RV_INPUT including the extension.\n"
        "Example: /data/TE123.fasta → TE123.fasta"
    ),
}


def encode_affixed(base: str, suffix: str) -> str:
    return f"{base}{_AFFIX_SEP}{suffix}"


def split_affixed(name: str) -> tuple[str, str] | None:
    """Return (base, suffix) for affixed params like RV_INPUT|_b.text or ${RV_INPUT}_b.text."""
    text = str(name or "").strip()
    if not text:
        return None
    if _AFFIX_SEP in text:
        base, suffix = text.split(_AFFIX_SEP, 1)
        base = base.strip()
        suffix = suffix.strip()
        if base and suffix:
            return base, suffix
        return None
    m = _AFFIX_DISPLAY.fullmatch(text)
    if m:
        return m.group(1), m.group(2)
    return None


def is_affixed_variable(name: str) -> bool:
    return split_affixed(name) is not None


def is_special_variable(name: str) -> bool:
    """True only for bare specials (RV_INPUT), not affixed (${RV_INPUT}_b.text)."""
    if split_affixed(name) is not None:
        return False
    return normalize_var_name(name) in SPECIAL_VARIABLES


def normalize_var_name(raw: str) -> str:
    """Normalize UI / stored names. Affixed forms stay as base|suffix."""
    text = str(raw or "").strip()
    if not text:
        return ""
    affixed = split_affixed(text)
    if affixed is not None:
        return encode_affixed(affixed[0], affixed[1])
    m = re.fullmatch(r"\$\{([A-Za-z_][A-Za-z0-9_]*)\}", text)
    if m:
        return m.group(1)
    m = re.fullmatch(r"\{(RV_[A-Za-z0-9_]+)\}", text)
    if m:
        return m.group(1)
    if re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", text):
        return text
    return text


def display_var_name(name: str) -> str:
    """Show `${name}` or `${RV_INPUT}_b.text` in the designer UI."""
    key = normalize_var_name(name)
    if not key:
        return ""
    affixed = split_affixed(key)
    if affixed is not None:
        return f"${{{affixed[0]}}}{affixed[1]}"
    return f"${{{key}}}"


def special_variable_mode(name: str) -> str:
    """UI mode for a special or affixed-from-special variable."""
    affixed = split_affixed(name)
    key = affixed[0] if affixed is not None else normalize_var_name(name)
    if key in {"RV_INPUT", "RV_GENOME", "RV_INPUTS"}:
        return "file"
    if key in {"RV_FOLDER", "RV_OUTDIR"}:
        return "folder"
    return "text"


def special_variables_help_text() -> str:
    lines = [
        "Special variables are filled automatically when you run a wrapper.",
        "Write them as ${RV_INPUT} (or {RV_INPUT}) in the command template.",
        "",
        "A glued suffix is its own parameter you can edit at run time:",
        "  ${RV_INPUT}_b.text  →  separate field, default = RV_INPUT path + _b.text",
        "  ${RV_OUTDIR}/${RV_STEM}_b.txt  →  use RV_STEM for name-without-extension",
        "",
    ]
    for key, desc in SPECIAL_VARIABLES.items():
        lines.append(f"${{{key}}}")
        for part in desc.split("\n"):
            lines.append(f"  {part}")
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def normalize_var_mode(mode: str) -> str:
    """Map stored / legacy modes to file, folder, text, integer, or float."""
    key = (mode or "text").strip().lower()
    if key in {"path", "input", "output"}:
        return "file"
    if key in {"number", "float", "double"}:
        return "float"
    if key in {"int", "integer"}:
        return "integer"
    if key in {"file", "folder", "text", "integer", "float"}:
        return key
    return "text"


def is_numeric_var_mode(mode: str) -> bool:
    return normalize_var_mode(mode) in {"integer", "float"}


def is_path_var_mode(mode: str) -> bool:
    return normalize_var_mode(mode) in {"file", "folder"}


def is_file_var_mode(mode: str) -> bool:
    return normalize_var_mode(mode) == "file"


def is_folder_var_mode(mode: str) -> bool:
    return normalize_var_mode(mode) == "folder"


def plugin_id_from_label(label: str) -> str:
    """Derive a stable internal plugin id from the user-visible name."""
    text = re.sub(r"[^a-zA-Z0-9]+", "_", (label or "").strip().lower())
    text = text.strip("_")
    return text or "plugin"


@dataclass
class CommandVarDef:
    name: str
    mode: str  # file|folder|text|integer|float
    value: str = ""


@dataclass
class UserCommandPluginSpec:
    id: str
    label: str
    command_template: str
    output_vars: list[str]
    var_defs: list[CommandVarDef]
    fasta_only: bool = True
    tooltip: str = ""
    show_on_toolbar: bool = False
    input_vars: list[str] | None = None

    def to_json(self) -> dict:
        return {
            "id": self.id,
            "label": self.label,
            "command_template": self.command_template,
            "input_vars": list(self.input_vars or []),
            "output_vars": list(self.output_vars),
            "var_defs": [
                {
                    "name": v.name,
                    "mode": normalize_var_mode(v.mode),
                    "value": v.value,
                }
                for v in self.var_defs
            ],
            "fasta_only": bool(self.fasta_only),
            "tooltip": self.tooltip,
            "show_on_toolbar": bool(self.show_on_toolbar),
        }

    @staticmethod
    def from_json(raw: dict) -> "UserCommandPluginSpec":
        defs = []
        for item in raw.get("var_defs", []):
            name = normalize_var_name(str(item.get("name", "")).strip())
            if not name:
                continue
            mode = normalize_var_mode(str(item.get("mode", "text")))
            value = str(item.get("value", ""))
            defs.append(
                CommandVarDef(
                    name=name,
                    mode=mode,
                    value=value,
                )
            )
        label = str(raw.get("label", "")).strip()
        plugin_id = str(raw.get("id", "")).strip() or plugin_id_from_label(label)
        return UserCommandPluginSpec(
            id=plugin_id,
            label=label,
            command_template=str(raw.get("command_template", "")),
            input_vars=[
                normalize_var_name(str(v).strip())
                for v in raw.get("input_vars", [])
                if normalize_var_name(str(v).strip())
            ],
            output_vars=[
                normalize_var_name(str(v).strip())
                for v in raw.get("output_vars", [])
                if normalize_var_name(str(v).strip())
            ],
            var_defs=defs,
            fasta_only=bool(raw.get("fasta_only", True)),
            tooltip=str(raw.get("tooltip", "")).strip(),
            show_on_toolbar=bool(raw.get("show_on_toolbar", False)),
        )

    def is_workflow_ready(self) -> bool:
        """True when this wrapper can be chained as a workflow step."""
        has_output = bool(self.output_vars)
        has_input = bool(self.input_vars)
        if not has_input:
            names = set(template_variables(self.command_template))
            has_input = bool(names & {"RV_INPUT", "RV_INPUTS"})
        return has_input and has_output


def template_variables(template: str) -> list[str]:
    names: list[str] = []
    seen: set[str] = set()
    for match in _VAR_PATTERN.finditer(template or ""):
        name = match.group(1) or match.group(2)
        if not name or name in seen:
            continue
        names.append(name)
        seen.add(name)
    return names


def template_affixed_variables(template: str) -> list[str]:
    """`${RV_INPUT}_b.text`-style tokens as their own parameters (base|suffix)."""
    names: list[str] = []
    seen: set[str] = set()
    for match in _VAR_REPLACE_PATTERN.finditer(template or ""):
        if match.group(1) is not None:
            base, suffix = match.group(1), match.group(2) or ""
        else:
            base, suffix = match.group(3), match.group(4) or ""
        if not base or not suffix:
            continue
        key = encode_affixed(base, suffix)
        if key in seen:
            continue
        seen.add(key)
        names.append(key)
    return names


def user_template_variables(template: str) -> list[str]:
    """Non-special variables that need Define-variables rows (includes affixed)."""
    names: list[str] = []
    seen: set[str] = set()
    for name in template_variables(template):
        if is_special_variable(name) or name in seen:
            continue
        names.append(name)
        seen.add(name)
    for name in template_affixed_variables(template):
        if name in seen:
            continue
        names.append(name)
        seen.add(name)
    return names


def expand_affixed_default(base_value: str, suffix: str) -> str:
    base = str(base_value or "").strip()
    if not base:
        return ""
    return f"{base}{suffix}"


def resolve_special_variables(
    *,
    selected_files: list[Path] | None = None,
    genome_file: Path | None = None,
    current_folder: Path | None = None,
    work_dir: Path | None = None,
) -> dict[str, str]:
    """Build values for RV_* placeholders from the current RepeatViewer context."""
    selected = [Path(p).expanduser().resolve() for p in (selected_files or [])]
    first = selected[0] if selected else None
    folder = (
        Path(current_folder).expanduser().resolve()
        if current_folder is not None
        else None
    )
    outdir: Path | None = work_dir
    if outdir is None:
        if first is not None and first.is_file():
            outdir = first.parent
        elif first is not None and first.is_dir():
            outdir = first
        elif folder is not None:
            outdir = folder
    outdir_text = (
        str(Path(outdir).expanduser().resolve()) if outdir is not None else ""
    )

    genome = ""
    if genome_file is not None:
        gpath = Path(genome_file).expanduser()
        if gpath.is_file():
            genome = str(gpath.resolve())

    return {
        "RV_INPUT": str(first) if first is not None else "",
        "RV_INPUTS": " ".join(str(p) for p in selected),
        "RV_GENOME": genome,
        "RV_FOLDER": str(folder) if folder is not None else "",
        "RV_OUTDIR": outdir_text,
        "RV_STEM": first.stem if first is not None else "",
        "RV_NAME": first.name if first is not None else "",
    }


def parse_var_defs(text: str) -> list[CommandVarDef]:
    defs: list[CommandVarDef] = []
    for raw_line in (text or "").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        name, sep, rhs = line.partition("=")
        if not sep:
            continue
        name = name.strip()
        mode, sep2, value = rhs.strip().partition(":")
        mode = normalize_var_mode(mode.strip() or "text")
        value = value if sep2 else ""
        if name:
            defs.append(CommandVarDef(name=name, mode=mode, value=value))
    return defs


def format_var_defs(defs: list[CommandVarDef]) -> str:
    return "\n".join(f"{v.name}={v.mode}:{v.value}" for v in defs)


def load_user_command_plugin_specs() -> list[UserCommandPluginSpec]:
    s = QSettings()
    s.beginGroup(_GROUP)
    raw = s.value(_KEY, "[]")
    s.endGroup()
    try:
        data = json.loads(str(raw))
    except (TypeError, json.JSONDecodeError):
        return []
    if not isinstance(data, list):
        return []
    out: list[UserCommandPluginSpec] = []
    for item in data:
        if isinstance(item, dict):
            spec = UserCommandPluginSpec.from_json(item)
            if spec.id and spec.label and spec.command_template:
                out.append(spec)
    return out


def save_user_command_plugin_spec(spec: UserCommandPluginSpec) -> None:
    specs = load_user_command_plugin_specs()
    replaced = False
    for idx, existing in enumerate(specs):
        if existing.id == spec.id:
            specs[idx] = spec
            replaced = True
            break
    if not replaced:
        specs.append(spec)
    _write_user_command_plugin_specs(specs)


def delete_user_command_plugin_spec(plugin_id: str) -> bool:
    specs = load_user_command_plugin_specs()
    new_specs = [item for item in specs if item.id != plugin_id]
    if len(new_specs) == len(specs):
        return False
    _write_user_command_plugin_specs(new_specs)
    return True


def get_user_command_plugin_spec(plugin_id: str) -> UserCommandPluginSpec | None:
    for spec in load_user_command_plugin_specs():
        if spec.id == plugin_id:
            return spec
    return None


def user_command_plugin_ids() -> set[str]:
    return {spec.id for spec in load_user_command_plugin_specs()}


def _write_user_command_plugin_specs(specs: list[UserCommandPluginSpec]) -> None:
    s = QSettings()
    s.beginGroup(_GROUP)
    s.setValue(_KEY, json.dumps([item.to_json() for item in specs]))
    s.endGroup()


class UserCommandPlugin(RepeatViewerPlugin):
    needs_selected_files = False
    fasta_only = False
    supports_parallel = False

    def __init__(self, spec: UserCommandPluginSpec) -> None:
        self.spec = spec
        self.id = spec.id
        self.label = spec.label
        self.tooltip = spec.tooltip
        self.fasta_only = False
        self._defs = {
            normalize_var_name(d.name): CommandVarDef(
                name=normalize_var_name(d.name),
                mode=d.mode,
                value=d.value,
            )
            for d in spec.var_defs
            if normalize_var_name(d.name)
        }
        for name in template_affixed_variables(spec.command_template):
            if name not in self._defs:
                self._defs[name] = CommandVarDef(
                    name=name,
                    mode=special_variable_mode(name),
                    value="",
                )
        self._active_run_dialog: UserCommandRunDialog | None = None

    def validate(self, files: list[Path], project: RepeatViewerProject) -> str | None:
        del files, project
        vars_in_cmd = set(user_template_variables(self.spec.command_template))
        missing = [v for v in vars_in_cmd if v not in self._defs and not is_affixed_variable(v)]
        if missing:
            return f"Missing variable definitions: {', '.join(sorted(display_var_name(v) for v in missing))}"
        opts = getattr(self, "_job_options", {}) or {}
        runtime = {
            str(k): str(v) for k, v in (opts.get("runtime_values") or {}).items()
        }
        if runtime:
            unset = []
            for name in template_variables(self.spec.command_template):
                if not str(runtime.get(name, "")).strip():
                    unset.append(display_var_name(name))
            for name in template_affixed_variables(self.spec.command_template):
                if not str(runtime.get(name, "")).strip():
                    unset.append(display_var_name(name))
            if unset:
                return f"Please set values for: {', '.join(unset)}"
        return None

    def prepare_job(self, parent, files: list[Path], project: RepeatViewerProject) -> bool:
        selected_files: list[str] = []
        current_folder = ""
        if parent is not None:
            getter = getattr(parent, "_selected_files", None)
            if callable(getter):
                try:
                    selected_files = [
                        str(Path(p).expanduser().resolve()) for p in getter()
                    ]
                except Exception:  # noqa: BLE001
                    selected_files = []
            folder = getattr(parent, "_current_folder", None)
            if folder is not None:
                try:
                    current_folder = str(Path(folder).expanduser().resolve())
                except Exception:  # noqa: BLE001
                    current_folder = str(folder)

        genome_file = ""
        if project.genome_file is not None and project.genome_file.is_file():
            genome_file = str(project.genome_file.expanduser().resolve())

        dialog = UserCommandRunDialog(
            self.spec,
            self._defs,
            parent=parent,
            selected_files=[Path(p) for p in selected_files],
            genome_file=Path(genome_file) if genome_file else None,
            current_folder=Path(current_folder) if current_folder else None,
        )
        self._active_run_dialog = dialog
        if parent is not None:
            parent._plugin_run_dialog = dialog
        dialog.setModal(False)
        dialog.setWindowModality(Qt.WindowModality.NonModal)
        # Keep the runtime dialog above the main window while the user
        # clicks files/paths in the panel.
        dialog.setWindowFlag(Qt.WindowType.WindowStaysOnTopHint, True)
        loop = QEventLoop()
        dialog.finished.connect(loop.quit)
        dialog.show()
        dialog.raise_()
        dialog.activateWindow()
        try:
            loop.exec()
            if dialog.result() != QDialog.DialogCode.Accepted:
                return False
            self._job_options = {
                "runtime_values": dialog.runtime_values(),
                "selected_files": selected_files,
                "genome_file": genome_file,
                "current_folder": current_folder,
            }
            return True
        finally:
            self._active_run_dialog = None
            if parent is not None:
                parent._plugin_run_dialog = None

    def run(self, file_path: Path, project: RepeatViewerProject, **job_options) -> list[str]:
        anchor = file_path.expanduser().resolve()
        if anchor.is_dir():
            work_dir = anchor
        elif anchor.is_file():
            work_dir = anchor.parent if anchor.parent.exists() else Path.cwd()
        else:
            work_dir = Path.cwd()
        runtime_values = {
            str(k): str(v) for k, v in (job_options.get("runtime_values") or {}).items()
        }
        selected_files = [
            Path(p) for p in (job_options.get("selected_files") or []) if str(p).strip()
        ]
        if not selected_files and anchor.is_file():
            selected_files = [anchor]
        # Workflow / batch: map the step FASTA onto declared input variables.
        if anchor.is_file():
            for name in self.spec.input_vars or []:
                if not str(runtime_values.get(name, "")).strip():
                    runtime_values[name] = str(anchor)
            if not str(runtime_values.get("RV_INPUT", "")).strip():
                runtime_values["RV_INPUT"] = str(anchor)
        genome_raw = job_options.get("genome_file") or (
            str(project.genome_file) if project.genome_file else ""
        )
        genome_path = Path(genome_raw).expanduser() if genome_raw else None
        folder_raw = job_options.get("current_folder") or ""
        folder_path = Path(folder_raw).expanduser() if folder_raw else None

        specials = resolve_special_variables(
            selected_files=selected_files,
            genome_file=genome_path if genome_path and genome_path.is_file() else None,
            current_folder=folder_path if folder_path and folder_path.is_dir() else None,
            work_dir=work_dir,
        )
        # Dialog-edited values win; resolved specials are only fallbacks.
        for key, value in specials.items():
            if key not in runtime_values or not str(runtime_values.get(key, "")).strip():
                if value:
                    runtime_values[key] = value
        outdir_raw = str(runtime_values.get("RV_OUTDIR", "")).strip()
        if outdir_raw:
            work_dir = Path(outdir_raw).expanduser().resolve()
        else:
            rv_input = str(runtime_values.get("RV_INPUT", "")).strip()
            if rv_input:
                input_path = Path(rv_input).expanduser()
                if input_path.is_file():
                    work_dir = input_path.parent.resolve()
                elif input_path.is_dir():
                    work_dir = input_path.resolve()

        context = build_command_context(
            self._defs,
            runtime_values=runtime_values,
            output_vars=set(self.spec.output_vars),
            work_dir=work_dir,
        )
        context.update({k: v for k, v in runtime_values.items() if str(v).strip()})
        # Ensure STEM/NAME/OUTDIR exist for suffixes like ${RV_STEM}_b.txt
        context = enrich_specials_from_input(context)
        commands = [
            _replace_vars(line, context).strip()
            for line in command_template_lines(self.spec.command_template)
        ]
        commands = [cmd for cmd in commands if cmd]
        if not commands:
            raise RuntimeError("No commands to run.")

        before_names: dict[str, set[str]] = {}
        for name in self.spec.output_vars:
            value = context.get(name, "")
            if not value:
                continue
            expected = Path(value).expanduser().resolve()
            defn = self._defs.get(name)
            if defn is not None and is_folder_var_mode(defn.mode):
                expected.mkdir(parents=True, exist_ok=True)
            parent = expected.parent if expected.parent != Path(".") else work_dir
            before_names[name] = _snapshot_directory_names(parent)

        for index, command in enumerate(commands, start=1):
            proc = run_logged_process(
                command,
                tool=f"plugin:{self.spec.label}",
                cwd=work_dir,
                shell=True,
            )
            if proc.returncode != 0:
                stderr = (proc.stderr or proc.stdout or "").strip()
                prefix = f"Command line {index} failed"
                raise RuntimeError(stderr or f"{prefix} with exit code {proc.returncode}")

        outputs: list[str] = []
        missing_outputs: list[str] = []
        for name in self.spec.output_vars:
            value = context.get(name, "")
            if not value:
                missing_outputs.append(name)
                continue
            expected = Path(value).expanduser().resolve()
            defn = self._defs.get(name)
            found = find_created_output_file(
                expected,
                before_names=before_names.get(name),
                is_folder=is_folder_var_mode(defn.mode) if defn else False,
            )
            if found is not None:
                outputs.append(str(found))
                context[name] = str(found)
            else:
                missing_outputs.append(f"{name} ({expected})")
        if missing_outputs:
            raise RuntimeError(
                "Expected output was not created: " + ", ".join(missing_outputs)
            )
        return outputs


def enrich_specials_from_input(values: dict[str, str]) -> dict[str, str]:
    """Fill RV_STEM / RV_NAME / RV_OUTDIR from RV_INPUT when missing."""
    out = {str(k): str(v) for k, v in values.items()}
    raw = str(out.get("RV_INPUT", "")).strip()
    if not raw:
        return out
    path = Path(raw.split()[0]).expanduser()
    if not str(out.get("RV_STEM", "")).strip():
        out["RV_STEM"] = path.stem
    if not str(out.get("RV_NAME", "")).strip():
        out["RV_NAME"] = path.name
    if not str(out.get("RV_OUTDIR", "")).strip():
        parent = path.parent
        if str(parent) not in {"", "."}:
            out["RV_OUTDIR"] = str(parent)
    return out


def _quote_for_shell(value: str) -> str:
    """Quote a substituted value when the shell would otherwise split it."""
    if not value:
        return value
    if re.search(r"[\s\"'\\$`;&|<>()\[\]{}?*~]", value):
        return shlex.quote(value)
    return value


def _replace_vars(template: str, context: dict[str, str]) -> str:
    def _repl(match: re.Match[str]) -> str:
        if match.group(1) is not None:
            key = match.group(1)
            suffix = match.group(2) or ""
        else:
            key = match.group(3)
            suffix = match.group(4) or ""
        if key is None:
            return match.group(0)
        if suffix:
            affixed_key = encode_affixed(key, suffix)
            affixed_val = str(context.get(affixed_key, "")).strip()
            if affixed_val:
                return _quote_for_shell(affixed_val)
        if key not in context:
            return match.group(0)
        value = str(context[key])
        if not value.strip():
            return match.group(0)
        return _quote_for_shell(f"{value}{suffix}")

    return _VAR_REPLACE_PATTERN.sub(_repl, template or "")


def command_template_lines(template: str) -> list[str]:
    lines: list[str] = []
    for raw in (template or "").splitlines():
        line = raw.strip()
        if line and not line.startswith("#"):
            lines.append(line)
    return lines


def build_command_context(
    defs: dict[str, CommandVarDef],
    *,
    runtime_values: dict[str, str],
    output_vars: set[str],
    work_dir: Path,
) -> dict[str, str]:
    context: dict[str, str] = {}
    pending = dict(defs)
    for _ in range(8):
        changed = False
        for name, defn in pending.items():
            if name in context:
                continue
            raw_value = runtime_values.get(name, "")
            if not str(raw_value).strip():
                raw_value = defn.value
            if not str(raw_value).strip():
                continue
            raw_value = str(raw_value)
            if is_numeric_var_mode(defn.mode):
                context[name] = raw_value
                changed = True
                continue
            expanded = _replace_vars(raw_value, context)
            if name in output_vars or is_path_var_mode(defn.mode):
                out_path = Path(expanded).expanduser()
                if not out_path.is_absolute():
                    out_path = work_dir / out_path
                context[name] = str(out_path.resolve())
            else:
                context[name] = expanded
            changed = True
        if not changed:
            break
    return context


def _snapshot_directory_names(directory: Path) -> set[str]:
    try:
        return {entry.name for entry in directory.iterdir()}
    except OSError:
        return set()


_COMMON_OUTPUT_SUFFIXES = (
    ".1aln",
    ".paf",
    ".psl",
    ".out",
    ".txt",
    ".gff",
    ".gff3",
    ".bed",
)


def candidate_output_paths(expected: Path) -> list[Path]:
    """Possible on-disk paths for a declared output (exact path + common extensions)."""
    path = expected.expanduser()
    parent = path.parent if path.parent != Path(".") else Path.cwd()
    name = path.name
    stem = path.stem if path.suffix else name

    candidates: list[Path] = []
    seen: set[str] = set()

    def add(candidate: Path) -> None:
        resolved = candidate.expanduser()
        key = str(resolved)
        if key in seen:
            return
        seen.add(key)
        candidates.append(resolved)

    add(path if path.is_absolute() else parent / name)

    extension_bases = [stem] if not path.suffix else [name, stem]
    for base in extension_bases:
        for ext in _COMMON_OUTPUT_SUFFIXES:
            if base.endswith(ext):
                continue
            add(parent / f"{base}{ext}")

    return candidates


def find_created_output_file(
    expected: Path,
    *,
    before_names: set[str] | None = None,
    is_folder: bool = False,
) -> Path | None:
    """Resolve an output path after a command finishes."""
    path = expected.expanduser()
    if path.is_file():
        return path.resolve()
    if is_folder and path.is_dir():
        return path.resolve()

    for candidate in candidate_output_paths(path):
        if candidate.is_file():
            return candidate.resolve()

    parent = path.parent if path.parent != Path(".") else Path.cwd()
    if not parent.is_dir():
        return None

    stem = path.stem if path.suffix else path.name
    expected_name = path.name

    if before_names is not None:
        for entry in sorted(parent.iterdir(), key=lambda p: p.name.lower()):
            if not entry.is_file():
                continue
            if entry.name in before_names:
                continue
            if entry.name == expected_name:
                return entry.resolve()
            if entry.stem == stem and entry.name != expected_name:
                return entry.resolve()
            if not path.suffix and entry.name.startswith(f"{stem}."):
                return entry.resolve()

    return None


def existing_output_files(
    spec: UserCommandPluginSpec,
    defs: dict[str, CommandVarDef],
    *,
    runtime_values: dict[str, str],
    work_dir: Path,
) -> list[Path]:
    """Output files that already exist on disk before running a plugin."""
    context = build_command_context(
        defs,
        runtime_values=runtime_values,
        output_vars=set(spec.output_vars),
        work_dir=work_dir.expanduser().resolve(),
    )
    existing: list[Path] = []
    seen: set[str] = set()
    for name in spec.output_vars:
        value = context.get(name, "")
        if not value:
            continue
        expected = Path(value).expanduser().resolve()
        defn = defs.get(name)
        is_folder = is_folder_var_mode(defn.mode) if defn is not None else False
        if is_folder:
            if expected.is_dir():
                key = str(expected)
                if key not in seen:
                    seen.add(key)
                    existing.append(expected)
            continue
        for candidate in candidate_output_paths(expected):
            if candidate.is_file():
                resolved = candidate.resolve()
                key = str(resolved)
                if key not in seen:
                    seen.add(key)
                    existing.append(resolved)
                break
    return existing


def resolve_plugin_output_paths(
    spec: UserCommandPluginSpec,
    defs: dict[str, CommandVarDef],
    *,
    runtime_values: dict[str, str],
    work_dir: Path,
) -> list[Path]:
    """Resolve marked output variables to absolute paths for overwrite checks."""
    context = build_command_context(
        defs,
        runtime_values=runtime_values,
        output_vars=set(spec.output_vars),
        work_dir=work_dir.expanduser().resolve(),
    )
    paths: list[Path] = []
    for name in spec.output_vars:
        value = context.get(name, "")
        if not value:
            continue
        expected = Path(value).expanduser().resolve()
        defn = defs.get(name)
        found = find_created_output_file(
            expected,
            is_folder=is_folder_var_mode(defn.mode) if defn is not None else False,
        )
        paths.append(found if found is not None else expected)
    return paths


class _FileInputRow(QWidget):
    activated = Signal(object)
    path_changed = Signal(str)

    def __init__(
        self,
        *,
        initial: str = "",
        parent: QWidget | None = None,
        allow_multiple: bool = False,
    ) -> None:
        super().__init__(parent)
        self._allow_multiple = allow_multiple
        self._edit = QLineEdit(initial, self)
        self._browse = QPushButton("Browse…", self)
        self._browse.clicked.connect(self._pick_path)
        self._edit.installEventFilter(self)
        self._edit.textChanged.connect(self.path_changed.emit)
        row = QHBoxLayout(self)
        row.setContentsMargins(0, 0, 0, 0)
        row.addWidget(self._edit, 1)
        row.addWidget(self._browse, 0)

    def eventFilter(self, obj, event) -> bool:  # noqa: N802
        if event.type() == QEvent.Type.FocusIn and obj is self._edit:
            self.activated.emit(self)
        return super().eventFilter(obj, event)

    def _pick_path(self) -> None:
        start = self._edit.text().strip().split()[0] if self._edit.text().strip() else str(Path.home())
        if self._allow_multiple:
            picked, _ = QFileDialog.getOpenFileNames(
                self, "Select files", start, "All files (*)"
            )
            if picked:
                self._edit.setText(" ".join(picked))
            return
        picked, _ = QFileDialog.getOpenFileName(self, "Select file", start, "All files (*)")
        if picked:
            self._edit.setText(picked)

    def set_path(self, path: Path) -> bool:
        resolved = path.expanduser().resolve()
        if not resolved.is_file():
            return False
        if self._allow_multiple:
            existing = self._edit.text().strip()
            text = str(resolved)
            if existing:
                parts = existing.split()
                if text not in parts:
                    parts.append(text)
                self._edit.setText(" ".join(parts))
            else:
                self._edit.setText(text)
            return True
        self._edit.setText(str(resolved))
        return True

    def focus_field(self) -> None:
        self._edit.setFocus(Qt.FocusReason.OtherFocusReason)

    def set_value(self, text: str) -> None:
        self._edit.setText(text)

    def value(self) -> str:
        return self._edit.text().strip()


class _FolderInputRow(QWidget):
    activated = Signal(object)
    path_changed = Signal(str)

    def __init__(self, *, initial: str = "", parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._edit = QLineEdit(initial, self)
        self._browse = QPushButton("Browse…", self)
        self._browse.clicked.connect(self._pick_path)
        self._edit.installEventFilter(self)
        self._edit.textChanged.connect(self.path_changed.emit)
        row = QHBoxLayout(self)
        row.setContentsMargins(0, 0, 0, 0)
        row.addWidget(self._edit, 1)
        row.addWidget(self._browse, 0)

    def eventFilter(self, obj, event) -> bool:  # noqa: N802
        if event.type() == QEvent.Type.FocusIn and obj is self._edit:
            self.activated.emit(self)
        return super().eventFilter(obj, event)

    def _pick_path(self) -> None:
        start = self._edit.text().strip() or str(Path.home())
        picked = QFileDialog.getExistingDirectory(self, "Select folder", start)
        if picked:
            self._edit.setText(picked)

    def set_path(self, path: Path) -> bool:
        resolved = path.expanduser().resolve()
        if not resolved.is_dir():
            return False
        self._edit.setText(str(resolved))
        return True

    def focus_field(self) -> None:
        self._edit.setFocus(Qt.FocusReason.OtherFocusReason)

    def set_value(self, text: str) -> None:
        self._edit.setText(text)

    def value(self) -> str:
        return self._edit.text().strip()


class _TextInputRow(QWidget):
    activated = Signal(object)
    text_changed = Signal(str)

    def __init__(self, *, initial: str = "", parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._edit = QLineEdit(initial, self)
        self._edit.installEventFilter(self)
        self._edit.textChanged.connect(self.text_changed.emit)
        row = QHBoxLayout(self)
        row.setContentsMargins(0, 0, 0, 0)
        row.addWidget(self._edit, 1)

    def eventFilter(self, obj, event) -> bool:  # noqa: N802
        if event.type() == QEvent.Type.FocusIn and obj is self._edit:
            self.activated.emit(self)
        return super().eventFilter(obj, event)

    def set_text(self, text: str) -> None:
        self._edit.setText(text)

    def set_value(self, text: str) -> None:
        self._edit.setText(text)

    def focus_field(self) -> None:
        self._edit.setFocus(Qt.FocusReason.OtherFocusReason)

    def value(self) -> str:
        return self._edit.text().strip()


class _IntegerInputRow(QWidget):
    activated = Signal(object)

    def __init__(self, *, initial: str = "", parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._spin = QSpinBox(self)
        self._spin.setRange(-2_147_483_648, 2_147_483_647)
        self._spin.setSingleStep(1)
        if initial:
            try:
                self._spin.setValue(int(float(initial)))
            except ValueError:
                pass
        self._spin.installEventFilter(self)
        row = QHBoxLayout(self)
        row.setContentsMargins(0, 0, 0, 0)
        row.addWidget(self._spin, 1)

    def eventFilter(self, obj, event) -> bool:  # noqa: N802
        if event.type() == QEvent.Type.FocusIn and obj is self._spin:
            self.activated.emit(self)
        return super().eventFilter(obj, event)

    def focus_field(self) -> None:
        self._spin.setFocus(Qt.FocusReason.OtherFocusReason)

    def value(self) -> str:
        return str(self._spin.value())


class _FloatInputRow(QWidget):
    activated = Signal(object)

    def __init__(self, *, initial: str = "", parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._spin = QDoubleSpinBox(self)
        self._spin.setRange(-1e12, 1e12)
        self._spin.setDecimals(6)
        self._spin.setSingleStep(1.0)
        if initial:
            try:
                self._spin.setValue(float(initial))
            except ValueError:
                pass
        self._spin.installEventFilter(self)
        row = QHBoxLayout(self)
        row.setContentsMargins(0, 0, 0, 0)
        row.addWidget(self._spin, 1)

    def eventFilter(self, obj, event) -> bool:  # noqa: N802
        if event.type() == QEvent.Type.FocusIn and obj is self._spin:
            self.activated.emit(self)
        return super().eventFilter(obj, event)

    def focus_field(self) -> None:
        self._spin.setFocus(Qt.FocusReason.OtherFocusReason)

    def value(self) -> str:
        text = f"{self._spin.value():.6f}".rstrip("0").rstrip(".")
        return text or "0"


class UserCommandRunDialog(QDialog):
    """Runtime parameter dialog for a user-designed command wrapper."""

    def __init__(
        self,
        spec: UserCommandPluginSpec,
        defs: dict[str, CommandVarDef],
        *,
        parent=None,
        selected_files: list[Path] | None = None,
        genome_file: Path | None = None,
        current_folder: Path | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle(spec.label)
        self.setMinimumWidth(560)
        self._active_path_row: _FileInputRow | _FolderInputRow | None = None
        self._rows: dict[str, QWidget] = {}
        self._modes: dict[str, str] = {}
        self._outputs = set(spec.output_vars)
        self._var_order: list[str] = []
        self._syncing_derived = False
        self._default_templates: dict[str, str] = {}
        self._auto_defaults: dict[str, str] = {}
        self._specials = resolve_special_variables(
            selected_files=selected_files,
            genome_file=genome_file,
            current_folder=current_folder,
        )

        title = QLabel(spec.label)
        title_font = title.font()
        title_font.setBold(True)
        title_font.setPointSize(title_font.pointSize() + 1)
        title.setFont(title_font)

        self._command_template = spec.command_template.strip()
        self._command_box = QTextEdit(self)
        self._command_box.setReadOnly(True)
        self._command_box.setPlainText(self._command_template)
        self._command_box.setFixedHeight(
            max(72, 24 * max(1, len(command_template_lines(self._command_template))))
        )
        self._command_box.setLineWrapMode(QTextEdit.LineWrapMode.WidgetWidth)

        form = QFormLayout()
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)
        vars_in_cmd = template_variables(spec.command_template)
        affixed_in_cmd = template_affixed_variables(spec.command_template)
        if not vars_in_cmd and not affixed_in_cmd:
            form.addRow("(No variables)", QLabel("This command has no ${variables}."))
        for name in vars_in_cmd:
            self._var_order.append(name)
            if is_special_variable(name):
                mode = special_variable_mode(name)
                initial = self._specials.get(name, "")
                self._modes[name] = mode
                if mode == "folder":
                    row: QWidget = _FolderInputRow(initial=initial, parent=self)
                    row.activated.connect(self._on_path_row_activated)
                elif mode == "file":
                    row = _FileInputRow(
                        initial=initial,
                        parent=self,
                        allow_multiple=(name == "RV_INPUTS"),
                    )
                    row.activated.connect(self._on_path_row_activated)
                else:
                    row = _TextInputRow(initial=initial, parent=self)
                self._rows[name] = row
                form.addRow(f"${{{name}}} ({mode}):", row)
                continue
            defn = defs.get(name)
            if defn is None:
                form.addRow(
                    f"${{{name}}}:",
                    QLabel("(undefined — edit wrapper to add this variable)"),
                )
                continue
            mode = normalize_var_mode(defn.mode)
            self._modes[name] = mode
            template = str(defn.value or "")
            self._default_templates[name] = template
            initial = self._expand_default(template, mode)
            self._auto_defaults[name] = initial
            if mode == "integer":
                row = _IntegerInputRow(initial=initial or defn.value, parent=self)
            elif mode == "float":
                row = _FloatInputRow(initial=initial or defn.value, parent=self)
            elif mode == "folder":
                row = _FolderInputRow(initial=initial, parent=self)
                row.activated.connect(self._on_path_row_activated)
            elif mode == "file":
                row = _FileInputRow(initial=initial, parent=self)
                row.activated.connect(self._on_path_row_activated)
            else:
                row = _TextInputRow(initial=initial, parent=self)
            self._rows[name] = row
            if name in self._outputs:
                kind = f"{mode} (output)"
            else:
                kind = mode
            form.addRow(f"${{{name}}} ({kind}):", row)

        for affixed_name in affixed_in_cmd:
            self._add_affixed_row(affixed_name, defs, form)

        input_row = self._rows.get("RV_INPUT")
        if isinstance(input_row, _FileInputRow):
            input_row.path_changed.connect(self._on_rv_input_changed)
        # Also refresh affixed defaults when any special base changes.
        for base in ("RV_GENOME", "RV_OUTDIR", "RV_STEM", "RV_FOLDER"):
            base_row = self._rows.get(base)
            if hasattr(base_row, "path_changed"):
                base_row.path_changed.connect(lambda _t="": self._refresh_affixed_defaults())  # type: ignore[union-attr]
            elif hasattr(base_row, "text_changed"):
                base_row.text_changed.connect(lambda _t="": self._refresh_affixed_defaults())  # type: ignore[union-attr]

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout = QVBoxLayout(self)
        layout.addWidget(title)
        if spec.tooltip:
            desc = QLabel(spec.tooltip)
            desc.setWordWrap(True)
            layout.addWidget(desc)
        layout.addWidget(QLabel("Command:", self))
        layout.addWidget(self._command_box)
        layout.addLayout(form)
        layout.addWidget(buttons)

    def _add_affixed_row(
        self,
        affixed_name: str,
        defs: dict[str, CommandVarDef],
        form: QFormLayout,
    ) -> None:
        parts = split_affixed(affixed_name)
        if parts is None:
            return
        base, suffix = parts
        defn = defs.get(affixed_name)
        mode = normalize_var_mode(defn.mode) if defn is not None else special_variable_mode(base)
        self._modes[affixed_name] = mode
        self._var_order.append(affixed_name)
        template = str(defn.value or "") if defn is not None else ""
        self._default_templates[affixed_name] = template
        if template:
            initial = self._expand_default(template, mode)
        else:
            initial = expand_affixed_default(self._specials.get(base, ""), suffix)
        self._auto_defaults[affixed_name] = initial
        if mode == "folder":
            row: QWidget = _FolderInputRow(initial=initial, parent=self)
            row.activated.connect(self._on_path_row_activated)
        elif mode == "file":
            row = _FileInputRow(initial=initial, parent=self)
            row.activated.connect(self._on_path_row_activated)
        else:
            row = _TextInputRow(initial=initial, parent=self)
        self._rows[affixed_name] = row
        label = display_var_name(affixed_name)
        kind = f"{mode} (output)" if affixed_name in self._outputs else mode
        form.addRow(f"{label} ({kind}):", row)

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        QTimer.singleShot(0, self._focus_first_field)

    def _focus_first_field(self) -> None:
        for name in self._var_order:
            row = self._rows.get(name)
            if row is None or not hasattr(row, "focus_field"):
                continue
            row.focus_field()  # type: ignore[union-attr]
            if isinstance(row, (_FileInputRow, _FolderInputRow)):
                self._active_path_row = row
            return

    def _on_path_row_activated(self, row: QWidget) -> None:
        if isinstance(row, (_FileInputRow, _FolderInputRow)):
            self._active_path_row = row

    def _work_dir_hint(self) -> Path | None:
        outdir = str(self._specials.get("RV_OUTDIR", "")).strip()
        if outdir:
            path = Path(outdir).expanduser()
            if path.is_dir():
                return path.resolve()
        rv_input = str(self._specials.get("RV_INPUT", "")).strip()
        if rv_input:
            path = Path(rv_input.split()[0]).expanduser()
            if path.is_file():
                return path.parent.resolve()
            if path.is_dir():
                return path.resolve()
        return None

    def _expand_default(self, template: str, mode: str) -> str:
        raw = str(template or "").strip()
        if not raw:
            return ""
        if "${" not in raw and "{RV_" not in raw:
            expanded = raw
        else:
            needed = set(template_variables(raw))
            if any(
                is_special_variable(name) and not str(self._specials.get(name, "")).strip()
                for name in needed
            ):
                return ""
            expanded = _replace_vars(raw, self._specials).strip()
            if not expanded or "${" in expanded or "{RV_" in expanded:
                return ""
        if is_path_var_mode(mode):
            path = Path(expanded).expanduser()
            if not path.is_absolute():
                base = self._work_dir_hint()
                if base is not None:
                    path = base / path
            try:
                return str(path.resolve())
            except OSError:
                return str(path)
        return expanded

    def _refresh_specials_from_rows(self) -> None:
        for key in SPECIAL_VARIABLES:
            row = self._rows.get(key)
            if row is not None and hasattr(row, "value"):
                self._specials[key] = row.value()  # type: ignore[union-attr]

    def _refresh_affixed_defaults(self) -> None:
        self._refresh_specials_from_rows()
        for name, row in list(self._rows.items()):
            parts = split_affixed(name)
            if parts is None or not hasattr(row, "value"):
                continue
            current = row.value()  # type: ignore[union-attr]
            previous = self._auto_defaults.get(name, "")
            if current and current != previous:
                continue
            base, suffix = parts
            template = self._default_templates.get(name, "")
            mode = self._modes.get(name, "file")
            if template:
                new_value = self._expand_default(template, mode)
            else:
                base_val = self._specials.get(base, "")
                if not base_val and hasattr(self._rows.get(base), "value"):
                    base_val = self._rows[base].value()  # type: ignore[union-attr]
                new_value = expand_affixed_default(base_val, suffix)
            self._auto_defaults[name] = new_value
            if hasattr(row, "set_value"):
                row.set_value(new_value)  # type: ignore[union-attr]

    def _refresh_auto_defaults(self) -> None:
        self._refresh_specials_from_rows()
        for name, template in self._default_templates.items():
            if split_affixed(name) is not None:
                continue
            if not template or ("${" not in template and "{RV_" not in template):
                continue
            row = self._rows.get(name)
            if row is None or not hasattr(row, "value"):
                continue
            current = row.value()  # type: ignore[union-attr]
            previous = self._auto_defaults.get(name, "")
            if current and current != previous:
                continue
            mode = self._modes.get(name, "text")
            new_value = self._expand_default(template, mode)
            self._auto_defaults[name] = new_value
            if hasattr(row, "set_value"):
                row.set_value(new_value)  # type: ignore[union-attr]
        self._refresh_affixed_defaults()

    def _on_rv_input_changed(self, text: str) -> None:
        if self._syncing_derived:
            return
        raw = str(text or "").strip()
        if not raw:
            stem = ""
            name = ""
            outdir = ""
        else:
            path = Path(raw.split()[0]).expanduser()
            stem = path.stem
            name = path.name
            outdir = str(path.parent) if path.parent != Path(".") else ""
        self._syncing_derived = True
        try:
            self._specials["RV_INPUT"] = raw
            self._specials["RV_STEM"] = stem
            self._specials["RV_NAME"] = name
            if outdir and not str(self._specials.get("RV_OUTDIR", "")).strip():
                self._specials["RV_OUTDIR"] = outdir
            stem_row = self._rows.get("RV_STEM")
            if isinstance(stem_row, _TextInputRow):
                stem_row.set_text(stem)
            name_row = self._rows.get("RV_NAME")
            if isinstance(name_row, _TextInputRow):
                name_row.set_text(name)
            out_row = self._rows.get("RV_OUTDIR")
            if isinstance(out_row, _FolderInputRow) and outdir:
                if not out_row.value():
                    out_row.set_value(outdir)
                    self._specials["RV_OUTDIR"] = outdir
            self._refresh_auto_defaults()
        finally:
            self._syncing_derived = False

    def apply_path_from_panel(self, path: Path) -> bool:
        if self._active_path_row is None:
            return False
        if isinstance(self._active_path_row, _FileInputRow) and not path.is_file():
            return False
        if isinstance(self._active_path_row, _FolderInputRow) and not path.is_dir():
            return False
        if not self._active_path_row.set_path(path):
            return False
        self.raise_()
        return True

    def runtime_values(self) -> dict[str, str]:
        values: dict[str, str] = {}
        for name, row in self._rows.items():
            if hasattr(row, "value"):
                values[name] = row.value()  # type: ignore[union-attr]
        return values

    def accept(self) -> None:
        missing = []
        for name, row in self._rows.items():
            value = row.value() if hasattr(row, "value") else ""  # type: ignore[union-attr]
            if not str(value).strip():
                missing.append(display_var_name(name) or name)
        if missing:
            QMessageBox.warning(
                self,
                self.windowTitle(),
                f"Please set values for: {', '.join(missing)}",
            )
            return
        super().accept()
