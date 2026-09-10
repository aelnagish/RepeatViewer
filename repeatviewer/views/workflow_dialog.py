"""Dialog to design and manage tool workflows."""

from __future__ import annotations

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QAbstractItemView,
    QDialog,
    QHBoxLayout,
    QInputDialog,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QMenu,
    QMessageBox,
    QPushButton,
    QSizePolicy,
    QVBoxLayout,
)

from repeatviewer.workflows import (
    WORKFLOW_CATEGORY_SEQUENCES,
    WORKFLOW_CATEGORY_TOOLS,
    WorkflowSpec,
    WorkflowStep,
    allowed_next_in_category,
    allowed_next_tools,
    delete_workflow_spec,
    duplicate_workflow_spec,
    list_plugin_workflow_tools,
    load_workflow_specs,
    resolve_tool_info,
    save_workflow_order,
    save_workflow_spec,
    validate_workflow_steps,
    workflow_id_from_name,
    ensure_default_workflows,
)

_ROLE_WORKFLOW_ID = Qt.ItemDataRole.UserRole
_ROLE_WORKFLOW_NAME = Qt.ItemDataRole.UserRole + 1
_ROLE_STEP_TOOL_ID = Qt.ItemDataRole.UserRole
_ROLE_STEP_KEEP = Qt.ItemDataRole.UserRole + 1


class WorkflowDesignerDialog(QDialog):
    """Design and manage named tool pipelines."""

    workflows_changed = Signal()

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Design and manage your workflows")
        self.setMinimumSize(760, 440)
        self._dirty = False
        self._current_id: str | None = None
        self._updating_steps = False
        self._updating_workflows = False
        self._suppress_workflow_select = False

        # --- Shared list row (aligned bottoms) ---
        self._workflow_list = QListWidget(self)
        self._workflow_list.setToolTip("Drag to change workflow menu order.")
        self._workflow_list.setDragDropMode(QAbstractItemView.DragDropMode.InternalMove)
        self._workflow_list.setDefaultDropAction(Qt.DropAction.MoveAction)
        self._workflow_list.currentItemChanged.connect(self._on_workflow_selected)
        self._workflow_list.model().rowsMoved.connect(self._on_workflows_reordered)

        self._steps_list = QListWidget(self)
        self._steps_list.setToolTip(
            "Ordered pipeline. Drag to reorder. "
            "Tick Keep to retain that step’s output files after the run "
            "(the last step is always kept)."
        )
        self._steps_list.setDragDropMode(QAbstractItemView.DragDropMode.InternalMove)
        self._steps_list.setDefaultDropAction(Qt.DropAction.MoveAction)
        self._steps_list.currentRowChanged.connect(lambda _r: self._update_controls())
        self._steps_list.model().rowsMoved.connect(self._on_steps_reordered)
        self._steps_list.itemChanged.connect(self._on_step_item_changed)

        # One column each so the Workflows list and Delete share the same width.
        self._new_btn = QPushButton("New", self)
        self._rename_btn = QPushButton("Rename", self)
        self._duplicate_btn = QPushButton("Duplicate", self)
        self._delete_btn = QPushButton("Delete", self)
        self._new_btn.clicked.connect(self._new_workflow)
        self._rename_btn.clicked.connect(self._rename_workflow)
        self._duplicate_btn.clicked.connect(self._duplicate_workflow)
        self._delete_btn.clicked.connect(self._delete_workflow)

        left_btns = QHBoxLayout()
        left_btns.setContentsMargins(0, 0, 0, 0)
        left_btns.setSpacing(6)
        for btn in (
            self._new_btn,
            self._rename_btn,
            self._duplicate_btn,
            self._delete_btn,
        ):
            btn.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
            left_btns.addWidget(btn)

        # Same QPushButton style/font as New/Rename (QToolButton draws smaller).
        self._seq_btn = QPushButton("Sequences", self)
        self._tools_btn = QPushButton("Tools", self)
        self._plugins_btn = QPushButton("Plugins", self)
        for btn in (self._seq_btn, self._tools_btn, self._plugins_btn):
            btn.setSizePolicy(QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Fixed)
        self._remove_btn = QPushButton("Remove", self)
        self._remove_btn.clicked.connect(self._remove_step)

        right_btns = QHBoxLayout()
        right_btns.setContentsMargins(0, 0, 0, 0)
        right_btns.setSpacing(6)
        right_btns.addWidget(self._seq_btn)
        right_btns.addWidget(self._tools_btn)
        right_btns.addWidget(self._plugins_btn)
        right_btns.addStretch(1)
        right_btns.addWidget(self._remove_btn)

        left_col = QVBoxLayout()
        left_col.setSpacing(6)
        left_col.addWidget(QLabel("Workflows", self))
        left_col.addWidget(self._workflow_list, 1)
        left_col.addLayout(left_btns)

        right_col = QVBoxLayout()
        right_col.setSpacing(6)
        right_col.addWidget(QLabel("Steps (run in order)", self))
        right_col.addWidget(
            QLabel("Tick a step to keep its output files after the workflow finishes.", self)
        )
        right_col.addWidget(self._steps_list, 1)
        right_col.addLayout(right_btns)

        columns = QHBoxLayout()
        columns.setSpacing(12)
        columns.addLayout(left_col, 1)
        columns.addLayout(right_col, 2)

        self._status = QLabel("", self)
        self._status.setWordWrap(True)
        self._status.setStyleSheet("color: #c44;")

        hint = QLabel(
            "Steps use each tool’s currently saved parameters. "
            "Only tools whose input format matches the previous output can be added. "
            "Unticked intermediate outputs are deleted after a successful run; "
            "the final step is always kept.",
            self,
        )
        hint.setWordWrap(True)

        self._close_btn = QPushButton("Close", self)
        self._save_btn = QPushButton("Save", self)
        self._close_btn.clicked.connect(self.reject)
        self._save_btn.clicked.connect(self._save_current)
        buttons = QHBoxLayout()
        buttons.addStretch(1)
        buttons.addWidget(self._close_btn)
        buttons.addWidget(self._save_btn)

        layout = QVBoxLayout(self)
        layout.setSpacing(10)
        layout.addLayout(columns, 1)
        layout.addWidget(self._status)
        layout.addWidget(hint)
        layout.addLayout(buttons)

        self._reload_workflows()
        self._update_controls()

    def _last_step_tool_id(self) -> str | None:
        if self._steps_list.count() == 0:
            return None
        return str(
            self._steps_list.item(self._steps_list.count() - 1).data(_ROLE_STEP_TOOL_ID)
            or ""
        ) or None

    def _make_step_item(
        self,
        tool_id: str,
        *,
        keep_output: bool = True,
        label: str | None = None,
    ) -> QListWidgetItem:
        info = resolve_tool_info(tool_id)
        text = label or (info.display_name if info else tool_id)
        item = QListWidgetItem(text)
        item.setFlags(
            item.flags()
            | Qt.ItemFlag.ItemIsUserCheckable
            | Qt.ItemFlag.ItemIsEnabled
            | Qt.ItemFlag.ItemIsSelectable
            | Qt.ItemFlag.ItemIsDragEnabled
        )
        item.setData(_ROLE_STEP_TOOL_ID, tool_id)
        item.setData(_ROLE_STEP_KEEP, bool(keep_output))
        item.setCheckState(
            Qt.CheckState.Checked if keep_output else Qt.CheckState.Unchecked
        )
        item.setToolTip(
            "Keep this step’s output files after the workflow finishes. "
            "The last step is always kept."
        )
        return item

    def _on_step_item_changed(self, item: QListWidgetItem) -> None:
        if self._updating_steps:
            return
        row = self._steps_list.row(item)
        last = self._steps_list.count() - 1
        keep = item.checkState() == Qt.CheckState.Checked
        if row == last and not keep:
            self._updating_steps = True
            try:
                item.setCheckState(Qt.CheckState.Checked)
                item.setData(_ROLE_STEP_KEEP, True)
            finally:
                self._updating_steps = False
            self._renumber_steps()
            return
        item.setData(_ROLE_STEP_KEEP, keep)
        self._renumber_steps()
        self._mark_dirty()

    def _rebuild_category_menus(self) -> None:
        after = self._last_step_tool_id()
        self._fill_menu_button(
            self._seq_btn,
            allowed_next_in_category(after, WORKFLOW_CATEGORY_SEQUENCES),
        )
        self._fill_menu_button(
            self._tools_btn,
            allowed_next_in_category(after, WORKFLOW_CATEGORY_TOOLS),
        )
        plugins: list = []
        try:
            plugin_allowed = {
                info.tool_id for info in allowed_next_tools(after)
            }
            plugins = [
                info
                for info in list_plugin_workflow_tools()
                if info.tool_id in plugin_allowed
            ]
        except Exception:  # noqa: BLE001
            plugins = []
        self._fill_menu_button(self._plugins_btn, plugins)

    def _fill_menu_button(self, button: QPushButton, tools) -> None:
        menu = QMenu(button)
        if not tools:
            empty = menu.addAction("(none available)")
            empty.setEnabled(False)
        else:
            for info in tools:
                action = menu.addAction(info.display_name)
                action.triggered.connect(
                    lambda checked=False, tid=info.tool_id: self._add_tool_id(tid)
                )
        button.setMenu(menu)
        button.setEnabled(self._current_id is not None and bool(tools))

    def _renumber_workflows(self) -> None:
        self._updating_workflows = True
        try:
            for row in range(self._workflow_list.count()):
                item = self._workflow_list.item(row)
                name = str(item.data(_ROLE_WORKFLOW_NAME) or item.text())
                item.setText(f"{row + 1}. {name}")
        finally:
            self._updating_workflows = False

    def _reload_workflows(self, select_id: str | None = None) -> None:
        self._workflow_list.blockSignals(True)
        self._workflow_list.clear()
        specs = ensure_default_workflows()
        for spec in specs:
            item = QListWidgetItem(spec.name)
            item.setData(_ROLE_WORKFLOW_ID, spec.id)
            item.setData(_ROLE_WORKFLOW_NAME, spec.name)
            self._workflow_list.addItem(item)
        self._renumber_workflows()
        target = select_id or self._current_id
        if target:
            for row in range(self._workflow_list.count()):
                item = self._workflow_list.item(row)
                if item.data(_ROLE_WORKFLOW_ID) == target:
                    self._workflow_list.setCurrentRow(row)
                    break
            else:
                if self._workflow_list.count():
                    self._workflow_list.setCurrentRow(0)
        elif self._workflow_list.count():
            self._workflow_list.setCurrentRow(0)
        self._workflow_list.blockSignals(False)

        item = self._workflow_list.currentItem()
        if item is None:
            self._current_id = None
            self._steps_list.clear()
            self._dirty = False
            self._update_controls()
            return
        workflow_id = str(item.data(_ROLE_WORKFLOW_ID) or "")
        self._current_id = workflow_id
        by_id = {s.id: s for s in specs}
        spec = by_id.get(workflow_id)
        if spec is None:
            name = str(item.data(_ROLE_WORKFLOW_NAME) or item.text())
            spec = WorkflowSpec(id=workflow_id, name=name, steps=[])
        self._dirty = False
        self._load_steps_for(spec)
        self._update_controls()

    def _on_workflows_reordered(self, *_args) -> None:
        if self._updating_workflows:
            return
        self._renumber_workflows()
        ids = [
            str(self._workflow_list.item(i).data(_ROLE_WORKFLOW_ID) or "")
            for i in range(self._workflow_list.count())
        ]
        ids = [wid for wid in ids if wid]
        saved = {spec.id for spec in load_workflow_specs()}
        if any(wid in saved for wid in ids):
            save_workflow_order([wid for wid in ids if wid in saved])
            self.workflows_changed.emit()

    def _current_spec_from_ui(self) -> WorkflowSpec | None:
        if self._current_id is None:
            return None
        name = "Workflow"
        for row in range(self._workflow_list.count()):
            item = self._workflow_list.item(row)
            if item.data(_ROLE_WORKFLOW_ID) == self._current_id:
                name = str(item.data(_ROLE_WORKFLOW_NAME) or "Workflow")
                break
        steps = []
        for row in range(self._steps_list.count()):
            item = self._steps_list.item(row)
            tool_id = str(item.data(_ROLE_STEP_TOOL_ID) or "")
            if not tool_id:
                continue
            keep = item.checkState() == Qt.CheckState.Checked
            # Final step is always retained on disk.
            if row == self._steps_list.count() - 1:
                keep = True
            steps.append(WorkflowStep(tool_id=tool_id, keep_output=keep))
        return WorkflowSpec(id=self._current_id, name=name, steps=steps)

    def _renumber_steps(self) -> None:
        self._updating_steps = True
        try:
            count = self._steps_list.count()
            for row in range(count):
                item = self._steps_list.item(row)
                tool_id = str(item.data(_ROLE_STEP_TOOL_ID) or "")
                info = resolve_tool_info(tool_id)
                label = info.display_name if info else tool_id
                keep = bool(item.data(_ROLE_STEP_KEEP))
                if row == count - 1:
                    keep = True
                    item.setData(_ROLE_STEP_KEEP, True)
                    item.setCheckState(Qt.CheckState.Checked)
                    suffix = " — keep (final)"
                elif keep:
                    item.setCheckState(Qt.CheckState.Checked)
                    suffix = " — keep"
                else:
                    item.setCheckState(Qt.CheckState.Unchecked)
                    suffix = ""
                item.setText(f"{row + 1}. {label}{suffix}")
                item.setFlags(
                    Qt.ItemFlag.ItemIsUserCheckable
                    | Qt.ItemFlag.ItemIsEnabled
                    | Qt.ItemFlag.ItemIsSelectable
                    | Qt.ItemFlag.ItemIsDragEnabled
                )
        finally:
            self._updating_steps = False

    def _load_steps_for(self, spec: WorkflowSpec) -> None:
        self._steps_list.clear()
        self._updating_steps = True
        try:
            for step in spec.steps:
                item = self._make_step_item(
                    step.tool_id,
                    keep_output=step.keep_output,
                )
                self._steps_list.addItem(item)
        finally:
            self._updating_steps = False
        self._renumber_steps()
        self._validate_inline()

    def _on_steps_reordered(self, *_args) -> None:
        if self._updating_steps:
            return
        self._renumber_steps()
        self._mark_dirty()

    def _confirm_leave_current(self, *, message: str | None = None) -> bool:
        """Ask to save/discard when leaving a dirty workflow. True = continue."""
        if not self._dirty or self._current_id is None:
            return True
        reply = QMessageBox.question(
            self,
            "Unsaved changes",
            message
            or "Save changes to the current workflow before continuing?",
            QMessageBox.StandardButton.Save
            | QMessageBox.StandardButton.Discard
            | QMessageBox.StandardButton.Cancel,
        )
        if reply == QMessageBox.StandardButton.Cancel:
            return False
        if reply == QMessageBox.StandardButton.Save:
            return self._save_current()
        self._dirty = False
        return True

    def _on_workflow_selected(
        self,
        current: QListWidgetItem | None,
        previous: QListWidgetItem | None,
    ) -> None:
        if getattr(self, "_suppress_workflow_select", False):
            return

        leaving_id = self._current_id
        if previous is not None:
            leaving_id = str(previous.data(_ROLE_WORKFLOW_ID) or "") or leaving_id

        if self._dirty and leaving_id:
            # Keep editing the workflow being left until the user chooses.
            self._current_id = leaving_id
            if not self._confirm_leave_current(
                message="Save changes to the current workflow before switching?"
            ):
                self._workflow_list.blockSignals(True)
                for row in range(self._workflow_list.count()):
                    item = self._workflow_list.item(row)
                    if item.data(_ROLE_WORKFLOW_ID) == leaving_id:
                        self._workflow_list.setCurrentRow(row)
                        break
                self._workflow_list.blockSignals(False)
                return

        if current is None:
            self._current_id = None
            self._steps_list.clear()
            self._update_controls()
            return
        workflow_id = str(current.data(_ROLE_WORKFLOW_ID) or "")
        self._current_id = workflow_id
        specs = {s.id: s for s in load_workflow_specs()}
        spec = specs.get(workflow_id)
        if spec is None:
            name = str(current.data(_ROLE_WORKFLOW_NAME) or current.text())
            spec = WorkflowSpec(id=workflow_id, name=name, steps=[])
        self._dirty = False
        self._load_steps_for(spec)
        self._update_controls()

    def _validate_inline(self) -> None:
        tool_ids = [
            str(self._steps_list.item(i).data(_ROLE_STEP_TOOL_ID) or "")
            for i in range(self._steps_list.count())
        ]
        error = validate_workflow_steps(tool_ids) if tool_ids else "Add at least one step."
        if tool_ids and error is None:
            self._status.clear()
        else:
            self._status.setText(error or "")

    def _mark_dirty(self) -> None:
        self._dirty = True
        self._rebuild_category_menus()
        self._validate_inline()
        self._update_controls()

    def _update_controls(self) -> None:
        has = self._current_id is not None
        self._rename_btn.setEnabled(has)
        self._duplicate_btn.setEnabled(has)
        self._delete_btn.setEnabled(has)
        self._remove_btn.setEnabled(has and self._steps_list.currentRow() >= 0)
        self._save_btn.setEnabled(has and self._dirty)
        self._rebuild_category_menus()

    def _new_workflow(self) -> None:
        if not self._confirm_leave_current():
            return
        existing = {s.id for s in load_workflow_specs()}
        for row in range(self._workflow_list.count()):
            existing.add(str(self._workflow_list.item(row).data(_ROLE_WORKFLOW_ID)))
        base_name = "Untitled"
        name = base_name
        workflow_id = workflow_id_from_name(name)
        n = 2
        while workflow_id in existing:
            name = f"{base_name} {n}"
            workflow_id = workflow_id_from_name(name)
            n += 1
        spec = WorkflowSpec(id=workflow_id, name=name, steps=[])
        item = QListWidgetItem(name)
        item.setData(_ROLE_WORKFLOW_ID, workflow_id)
        item.setData(_ROLE_WORKFLOW_NAME, name)
        self._suppress_workflow_select = True
        self._workflow_list.blockSignals(True)
        try:
            self._workflow_list.addItem(item)
            self._renumber_workflows()
            self._workflow_list.setCurrentItem(item)
        finally:
            self._workflow_list.blockSignals(False)
            self._suppress_workflow_select = False
        self._current_id = workflow_id
        self._dirty = True
        self._load_steps_for(spec)
        self._update_controls()

    def _rename_workflow(self) -> None:
        if self._current_id is None:
            return
        item = self._workflow_list.currentItem()
        if item is None:
            return
        current_name = str(item.data(_ROLE_WORKFLOW_NAME) or "")
        name, ok = QInputDialog.getText(
            self, "Rename workflow", "Name:", text=current_name
        )
        if not ok:
            return
        name = name.strip()
        if not name or name == current_name:
            return
        item.setData(_ROLE_WORKFLOW_NAME, name)
        self._renumber_workflows()
        # Persist rename immediately when the workflow already exists so menus
        # and toolbar labels update without waiting for Save.
        from repeatviewer.workflows import get_workflow_spec

        if get_workflow_spec(self._current_id) is not None:
            spec = self._current_spec_from_ui()
            if spec is not None:
                error = validate_workflow_steps(spec.tool_ids())
                if error is None:
                    try:
                        save_workflow_spec(spec)
                    except Exception as exc:  # noqa: BLE001
                        QMessageBox.warning(self, "Rename workflow", str(exc))
                        item.setData(_ROLE_WORKFLOW_NAME, current_name)
                        self._renumber_workflows()
                        return
                    self._dirty = False
                    self._update_controls()
                    self.workflows_changed.emit()
                    return
        self._mark_dirty()

    def _duplicate_workflow(self) -> None:
        if self._current_id is None:
            return
        if not self._confirm_leave_current():
            return
        copy = duplicate_workflow_spec(self._current_id)
        if copy is None:
            QMessageBox.warning(self, "Duplicate", "Save the workflow before duplicating.")
            return
        self._reload_workflows(select_id=copy.id)
        self.workflows_changed.emit()

    def _delete_workflow(self) -> None:
        if self._current_id is None:
            return
        reply = QMessageBox.question(
            self,
            "Delete workflow",
            "Delete this workflow?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
        )
        if reply != QMessageBox.StandardButton.Yes:
            return
        delete_workflow_spec(self._current_id)
        self._dirty = False
        self._current_id = None
        self._reload_workflows()
        self.workflows_changed.emit()

    def _add_tool_id(self, tool_id: str) -> None:
        if self._current_id is None:
            return
        info = resolve_tool_info(tool_id)
        if info is None:
            return
        after = self._last_step_tool_id()
        if after is not None:
            from repeatviewer.workflows import tools_compatible

            if not tools_compatible(after, tool_id):
                QMessageBox.warning(
                    self,
                    "Add step",
                    f"Cannot add {info.display_name} after the previous step.",
                )
                return
        item = self._make_step_item(info.tool_id, keep_output=True)
        self._steps_list.addItem(item)
        self._steps_list.setCurrentItem(item)
        self._renumber_steps()
        self._mark_dirty()

    def _remove_step(self) -> None:
        row = self._steps_list.currentRow()
        if row < 0:
            return
        self._steps_list.takeItem(row)
        self._renumber_steps()
        self._mark_dirty()

    def _save_current(self) -> bool:
        spec = self._current_spec_from_ui()
        if spec is None:
            return False
        error = validate_workflow_steps(spec.tool_ids())
        if error:
            QMessageBox.warning(self, "Save workflow", error)
            return False
        try:
            save_workflow_spec(spec)
        except Exception as exc:  # noqa: BLE001
            QMessageBox.warning(self, "Save workflow", str(exc))
            return False
        ids = [
            str(self._workflow_list.item(i).data(_ROLE_WORKFLOW_ID) or "")
            for i in range(self._workflow_list.count())
        ]
        save_workflow_order([wid for wid in ids if wid])
        self._dirty = False
        self._reload_workflows(select_id=spec.id)
        self.workflows_changed.emit()
        return True

    def reject(self) -> None:
        if self._dirty:
            reply = QMessageBox.question(
                self,
                "Unsaved changes",
                "Discard unsaved changes?",
                QMessageBox.StandardButton.Discard | QMessageBox.StandardButton.Cancel,
            )
            if reply != QMessageBox.StandardButton.Discard:
                return
        super().reject()
