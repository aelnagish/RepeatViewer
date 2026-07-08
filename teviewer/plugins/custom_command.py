"""User-defined command plugins configured inside TEviewer."""

from __future__ import annotations

import json
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path

from PySide6.QtCore import QEvent, QEventLoop, QSettings, Qt, Signal
from PySide6.QtWidgets import (
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
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

from teviewer.plugins.api import TEviewerPlugin
from teviewer.project import TEviewerProject

_GROUP = "plugins"
_KEY = "custom_commands_json"
_VAR_PATTERN = re.compile(r"\$\{([A-Za-z_][A-Za-z0-9_]*)\}")


def normalize_var_mode(mode: str) -> str:
    """Map stored / legacy modes to path, text, or number."""
    key = (mode or "text").strip().lower()
    if key in {"file", "input", "output"}:
        return "path"
    if key in {"path", "text", "number"}:
        return key
    return "text"


@dataclass
class CommandVarDef:
    name: str
    mode: str  # path|text|number
    value: str = ""
    output_suffix: str = ""  # e.g. .1aln when tools append an extension


@dataclass
class UserCommandPluginSpec:
    id: str
    label: str
    command_template: str
    output_vars: list[str]
    var_defs: list[CommandVarDef]
    fasta_only: bool = True
    tooltip: str = ""

    def to_json(self) -> dict:
        return {
            "id": self.id,
            "label": self.label,
            "command_template": self.command_template,
            "output_vars": list(self.output_vars),
            "var_defs": [
                {
                    "name": v.name,
                    "mode": normalize_var_mode(v.mode),
                    "value": v.value,
                    "output_suffix": v.output_suffix,
                }
                for v in self.var_defs
            ],
            "fasta_only": bool(self.fasta_only),
            "tooltip": self.tooltip,
        }

    @staticmethod
    def from_json(raw: dict) -> "UserCommandPluginSpec":
        defs = []
        for item in raw.get("var_defs", []):
            name = str(item.get("name", "")).strip()
            if not name:
                continue
            mode = normalize_var_mode(str(item.get("mode", "text")))
            value = str(item.get("value", ""))
            output_suffix = str(item.get("output_suffix", ""))
            defs.append(
                CommandVarDef(
                    name=name,
                    mode=mode,
                    value=value,
                    output_suffix=output_suffix,
                )
            )
        return UserCommandPluginSpec(
            id=str(raw.get("id", "")).strip(),
            label=str(raw.get("label", "")).strip(),
            command_template=str(raw.get("command_template", "")),
            output_vars=[str(v).strip() for v in raw.get("output_vars", []) if str(v).strip()],
            var_defs=defs,
            fasta_only=bool(raw.get("fasta_only", True)),
            tooltip=str(raw.get("tooltip", "")).strip(),
        )


def template_variables(template: str) -> list[str]:
    names: list[str] = []
    seen: set[str] = set()
    for name in _VAR_PATTERN.findall(template or ""):
        if name not in seen:
            names.append(name)
            seen.add(name)
    return names


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


class UserCommandPlugin(TEviewerPlugin):
    needs_selected_files = False
    fasta_only = False
    supports_parallel = False

    def __init__(self, spec: UserCommandPluginSpec) -> None:
        self.spec = spec
        self.id = spec.id
        self.label = spec.label
        self.tooltip = spec.tooltip
        self.fasta_only = False
        self._defs = {d.name: d for d in spec.var_defs}
        self._active_run_dialog: UserCommandRunDialog | None = None

    def validate(self, files: list[Path], project: TEviewerProject) -> str | None:
        vars_in_cmd = set(template_variables(self.spec.command_template))
        missing = [v for v in vars_in_cmd if v not in self._defs]
        if missing:
            return f"Missing variable definitions: {', '.join(sorted(missing))}"
        return None

    def prepare_job(self, parent, files: list[Path], project: TEviewerProject) -> bool:
        dialog = UserCommandRunDialog(self.spec, self._defs, parent=parent)
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
            self._job_options = {"runtime_values": dialog.runtime_values()}
            return True
        finally:
            self._active_run_dialog = None
            if parent is not None:
                parent._plugin_run_dialog = None

    def run(self, file_path: Path, project: TEviewerProject, **job_options) -> list[str]:
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
        context = build_command_context(
            self._defs,
            runtime_values=runtime_values,
            output_vars=set(self.spec.output_vars),
            work_dir=work_dir,
        )
        commands = [
            _replace_vars(line, context).strip()
            for line in command_template_lines(self.spec.command_template)
        ]
        commands = [cmd for cmd in commands if cmd]
        if not commands:
            raise RuntimeError("No commands to run.")

        for index, command in enumerate(commands, start=1):
            proc = subprocess.run(
                command,
                shell=True,
                cwd=str(work_dir),
                capture_output=True,
                text=True,
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
            suffix = defn.output_suffix if defn is not None else ""
            found = find_created_output_file(expected, suffix=suffix)
            if found is not None:
                outputs.append(str(found))
                context[name] = str(found)
            else:
                missing_outputs.append(f"{name} ({expected})")
        if missing_outputs:
            hint = (
                " If the tool adds a file extension automatically, set “Added suffix” "
                "for that output variable in the plugin design."
            )
            raise RuntimeError(
                "Expected output file(s) were not created: "
                + ", ".join(missing_outputs)
                + hint
            )
        return outputs


def _replace_vars(template: str, context: dict[str, str]) -> str:
    def _repl(match: re.Match[str]) -> str:
        key = match.group(1)
        return context.get(key, match.group(0))

    return _VAR_PATTERN.sub(_repl, template or "")


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
            if normalize_var_mode(defn.mode) == "number":
                context[name] = raw_value
                changed = True
                continue
            expanded = _replace_vars(raw_value, context)
            if name in output_vars or normalize_var_mode(defn.mode) == "path":
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


_COMMON_OUTPUT_SUFFIXES = (
    ".1aln",
    ".paf",
    ".psl",
    ".sam",
    ".bam",
    ".out",
    ".txt",
    ".gff",
    ".gff3",
    ".bed",
)


def normalize_output_suffix(suffix: str) -> str:
    text = suffix.strip()
    if not text:
        return ""
    return text if text.startswith(".") else f".{text}"


def candidate_output_paths(expected: Path, *, suffix: str = "") -> list[Path]:
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

    user_suffix = normalize_output_suffix(suffix)
    if user_suffix and not name.endswith(user_suffix):
        add(parent / f"{stem}{user_suffix}")

    extension_bases = [stem] if not path.suffix else [name, stem]
    for base in extension_bases:
        for ext in _COMMON_OUTPUT_SUFFIXES:
            if base.endswith(ext):
                continue
            add(parent / f"{base}{ext}")

    return candidates


def find_created_output_file(expected: Path, *, suffix: str = "") -> Path | None:
    for candidate in candidate_output_paths(expected, suffix=suffix):
        if candidate.is_file():
            return candidate.resolve()
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
        suffix = defn.output_suffix if defn is not None else ""
        for candidate in candidate_output_paths(expected, suffix=suffix):
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
        suffix = defn.output_suffix if defn is not None else ""
        found = find_created_output_file(expected, suffix=suffix)
        paths.append(found if found is not None else expected)
    return paths


class _PathInputRow(QWidget):
    activated = Signal(object)

    def __init__(self, *, initial: str = "", parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._edit = QLineEdit(initial, self)
        self._browse = QPushButton("Browse…", self)
        self._browse.clicked.connect(self._pick_path)
        self._edit.installEventFilter(self)
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
        picked, _ = QFileDialog.getOpenFileName(self, "Select file or folder", start, "All files (*)")
        if picked:
            self._edit.setText(picked)

    def set_path(self, path: Path) -> None:
        self._edit.setText(str(path.expanduser().resolve()))

    def value(self) -> str:
        return self._edit.text().strip()


class _TextInputRow(QWidget):
    activated = Signal(object)

    def __init__(self, *, initial: str = "", parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._edit = QLineEdit(initial, self)
        self._edit.installEventFilter(self)
        row = QHBoxLayout(self)
        row.setContentsMargins(0, 0, 0, 0)
        row.addWidget(self._edit, 1)

    def eventFilter(self, obj, event) -> bool:  # noqa: N802
        if event.type() == QEvent.Type.FocusIn and obj is self._edit:
            self.activated.emit(self)
        return super().eventFilter(obj, event)

    def value(self) -> str:
        return self._edit.text().strip()


class _NumberInputRow(QWidget):
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

    def value(self) -> str:
        return str(self._spin.value())


class UserCommandRunDialog(QDialog):
    """Runtime parameter dialog for a user-designed command plugin."""

    def __init__(
        self,
        spec: UserCommandPluginSpec,
        defs: dict[str, CommandVarDef],
        *,
        parent=None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle(spec.label)
        self.setMinimumWidth(560)
        self._active_path_row: _PathInputRow | None = None
        self._rows: dict[str, QWidget] = {}
        self._modes: dict[str, str] = {}
        self._outputs = set(spec.output_vars)

        title = QLabel(spec.label)
        title_font = title.font()
        title_font.setBold(True)
        title_font.setPointSize(title_font.pointSize() + 1)
        title.setFont(title_font)

        hint = QLabel(
            "Set a value for each variable (defaults are pre-filled). For path fields, click a "
            "file or folder in the panel to fill the active path."
        )
        hint.setWordWrap(True)

        command_box = QTextEdit(self)
        command_box.setReadOnly(True)
        command_box.setPlainText(spec.command_template.strip())
        command_box.setFixedHeight(max(72, 24 * max(1, len(command_template_lines(spec.command_template)))))
        command_box.setLineWrapMode(QTextEdit.LineWrapMode.WidgetWidth)

        form = QFormLayout()
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)
        vars_in_cmd = template_variables(spec.command_template)
        if not vars_in_cmd:
            form.addRow("(No variables)", QLabel("This command has no ${variables}."))
        for name in vars_in_cmd:
            defn = defs.get(name)
            if defn is None:
                form.addRow(
                    f"{name}:",
                    QLabel("(undefined — edit plugin to add this variable)"),
                )
                continue
            mode = normalize_var_mode(defn.mode)
            self._modes[name] = mode
            if mode == "number":
                row: QWidget = _NumberInputRow(initial=defn.value, parent=self)
            elif mode == "path":
                row = _PathInputRow(initial=defn.value, parent=self)
                row.activated.connect(self._on_path_row_activated)
            else:
                row = _TextInputRow(initial=defn.value, parent=self)
            self._rows[name] = row
            kind = "path (output)" if name in self._outputs else mode
            form.addRow(f"{name} ({kind}):", row)

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
        layout.addWidget(command_box)
        layout.addWidget(hint)
        layout.addLayout(form)
        layout.addWidget(buttons)

        path_rows = [row for row in self._rows.values() if isinstance(row, _PathInputRow)]
        if path_rows:
            self._active_path_row = path_rows[0]

    def _on_path_row_activated(self, row: _PathInputRow) -> None:
        self._active_path_row = row

    def apply_path_from_panel(self, path: Path) -> bool:
        if self._active_path_row is None:
            path_rows = [row for row in self._rows.values() if isinstance(row, _PathInputRow)]
            if not path_rows:
                return False
            self._active_path_row = path_rows[0]
        self._active_path_row.set_path(path)
        # Reorder z-stack after panel click; keeps dialog visible.
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
                missing.append(name)
        if missing:
            QMessageBox.warning(
                self,
                self.windowTitle(),
                f"Please set values for: {', '.join(missing)}",
            )
            return
        super().accept()
