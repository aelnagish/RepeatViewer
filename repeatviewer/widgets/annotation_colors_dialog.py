"""Dialog for editing annotation element colors."""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (
    QAbstractItemView,
    QDialog,
    QDialogButtonBox,
    QFormLayout,
    QFrame,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QVBoxLayout,
    QWidget,
    QColorDialog,
)

from repeatviewer.backends.annotation_colors import (
    DEFAULT_ADD_FEATURE_BUTTONS,
    annotation_feature_type_labels_for_colors,
    default_feature_type_colors,
    feature_type_color,
    load_add_feature_buttons,
    load_feature_type_colors,
    save_add_feature_buttons,
    save_feature_type_colors,
    save_similarity_database_colors,
    similarity_database_color,
    snapshot_similarity_database_colors,
    slug_feature_type,
)
from repeatviewer.backends.msa_plot import (
    AMINO_ACID_CATEGORIES,
    AMINO_ACID_COLOR_PRESETS,
    DEFAULT_AMINO_ACID_COLORS,
    DEFAULT_NUCLEOTIDE_COLORS,
    NUCLEOTIDE_CATEGORIES,
    NUCLEOTIDE_COLOR_PRESETS,
)
from repeatviewer.backends.similarity_annotation import SIMILARITY_DB_COLORS
from repeatviewer.database_settings import load_similarity_annotation_databases
from repeatviewer.widgets.preview_zoom import (
    default_amino_acid_backgrounds,
    default_nucleotide_backgrounds,
    load_custom_amino_acid_backgrounds,
    load_custom_nucleotide_backgrounds,
    save_color_blind_nucleotide_scheme,
    save_custom_amino_acid_backgrounds,
    save_custom_nucleotide_backgrounds,
)
from repeatviewer.widgets.residue_color_grid import ResidueColorGrid


class _ColorPickerRow(QWidget):
    def __init__(
        self,
        label: str,
        color: str,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self._color = QColor(color)
        self._button = QPushButton(self)
        self._button.setFixedSize(56, 22)
        self._button.setToolTip(f"Choose color for {label}")
        self._button.clicked.connect(self._pick_color)
        self._update_button()

        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(QLabel(label, self))
        layout.addStretch(1)
        layout.addWidget(self._button)

    def _update_button(self) -> None:
        hex_color = self._color.name(QColor.NameFormat.HexRgb)
        self._button.setStyleSheet(f"background-color: {hex_color};")

    def _pick_color(self) -> None:
        chosen = QColorDialog.getColor(self._color, self, "Choose color")
        if not chosen.isValid():
            return
        self._color = chosen
        self._update_button()

    def set_color_hex(self, color: str) -> None:
        self._color = QColor(color)
        self._update_button()

    def color_hex(self) -> str:
        return self._color.name(QColor.NameFormat.HexRgb)


def _tile_letter(label: str) -> str:
    text = str(label or "").strip()
    return text[0].upper() if text else "?"


class AnnotationElementColorsPanel(QWidget):
    """Embeddable colors for nucleotides, amino acids, features, and databases.

    Pass ``show_residue_colors=False`` for Sequence Annotation Advanced Settings
    (feature types + similarity databases only; no DNA / protein palettes).
    """

    def __init__(
        self,
        parent: QWidget | None = None,
        *,
        show_residue_colors: bool = True,
    ) -> None:
        super().__init__(parent)
        self._show_residue_colors = bool(show_residue_colors)
        databases = load_similarity_annotation_databases()

        body_layout = QVBoxLayout(self)
        body_layout.setContentsMargins(0, 0, 0, 0)
        body_layout.setSpacing(12)

        self._nuc_grid: ResidueColorGrid | None = None
        self._aa_grid: ResidueColorGrid | None = None
        if self._show_residue_colors:
            # --- Nucleotides (MSA figure-2 style) ---
            body_layout.addWidget(self._section_label("Nucleotide colors"))
            nuc_colors = dict(DEFAULT_NUCLEOTIDE_COLORS)
            nuc_colors.update(load_custom_nucleotide_backgrounds())
            if "U" in nuc_colors:
                nuc_colors["T"] = nuc_colors["U"]
            self._nuc_grid = ResidueColorGrid(
                self,
                colors=nuc_colors,
                categories=list(NUCLEOTIDE_CATEGORIES),
                linked_residues={"T": "U"},
            )
            nuc_presets = QHBoxLayout()
            nuc_presets.setContentsMargins(0, 0, 0, 0)
            nuc_presets.setSpacing(6)
            nuc_preset_label = QLabel("Preset:", self)
            nuc_preset_label.setFixedWidth(ResidueColorGrid.LABEL_WIDTH)
            nuc_presets.addWidget(nuc_preset_label)
            nuc_tips = {
                "Default": "Standard A/C/G/T colours",
                "Color-blind": (
                    "Wong / Paul Tol palette (same as Annotation → nucleotide "
                    "Color-blind)"
                ),
            }
            for label, palette in NUCLEOTIDE_COLOR_PRESETS:
                btn = QPushButton(label, self)
                btn.setToolTip(
                    nuc_tips.get(label, f"Apply {label} nucleotide colours")
                )
                btn.clicked.connect(
                    lambda _checked=False, name=label, colors=palette: (
                        self._apply_nuc_preset(name, colors)
                    )
                )
                nuc_presets.addWidget(btn)
            nuc_presets.addStretch(1)
            body_layout.addLayout(nuc_presets)
            body_layout.addWidget(self._nuc_grid)

            # --- Amino acids (MSA style) ---
            body_layout.addWidget(self._section_label("Amino-acid colors"))
            aa_colors = dict(DEFAULT_AMINO_ACID_COLORS)
            aa_colors.update(load_custom_amino_acid_backgrounds())
            self._aa_grid = ResidueColorGrid(
                self,
                colors=aa_colors,
                categories=list(AMINO_ACID_CATEGORIES),
            )
            aa_presets = QHBoxLayout()
            aa_presets.setContentsMargins(0, 0, 0, 0)
            aa_presets.setSpacing(6)
            aa_preset_label = QLabel("Preset:", self)
            aa_preset_label.setFixedWidth(ResidueColorGrid.LABEL_WIDTH)
            aa_presets.addWidget(aa_preset_label)
            aa_tips = {
                "Clustal": "Clustal X default amino-acid colours",
                "Zappo": "Zappo physicochemical palette",
                "Taylor": "Taylor / Lesk palette",
            }
            for label, palette in AMINO_ACID_COLOR_PRESETS:
                btn = QPushButton(label, self)
                btn.setToolTip(aa_tips.get(label, f"Apply {label} amino-acid colours"))
                btn.clicked.connect(
                    lambda _checked=False, colors=palette: self._aa_grid.set_colors(
                        colors
                    )
                )
                aa_presets.addWidget(btn)
            aa_presets.addStretch(1)
            body_layout.addLayout(aa_presets)
            body_layout.addWidget(self._aa_grid)

        # --- Feature types (label + colour swatch, no letter) ---
        body_layout.addWidget(self._section_label("Feature types"))
        feature_labels = list(annotation_feature_type_labels_for_colors())
        feature_categories = [
            (display, (feature_type,)) for display, feature_type in feature_labels
        ]
        feature_colors = {
            feature_type: feature_type_color(feature_type)
            for _display, feature_type in feature_labels
        }
        self._feature_grid = ResidueColorGrid(
            self,
            colors=feature_colors,
            categories=feature_categories,
            button_labels={feature_type: "" for _d, feature_type in feature_labels},
        )
        body_layout.addWidget(self._feature_grid)

        # --- Similarity databases (label + colour swatch, no letter) ---
        body_layout.addWidget(self._section_label("Similarity databases"))
        if databases:
            db_categories = [(entry.name, (entry.name,)) for entry in databases]
            db_colors = {
                entry.name: similarity_database_color(entry.name, index)
                for index, entry in enumerate(databases)
            }
            self._database_grid: ResidueColorGrid | None = ResidueColorGrid(
                self,
                colors=db_colors,
                categories=db_categories,
                button_labels={entry.name: "" for entry in databases},
            )
            body_layout.addWidget(self._database_grid)
        else:
            self._database_grid = None
            body_layout.addWidget(
                QLabel("No databases defined. Add them in File → Manage Database.", self)
            )

        self._old_feature_colors = load_feature_type_colors()
        self._old_database_colors = snapshot_similarity_database_colors(
            [entry.name for entry in databases]
        )

    def _apply_nuc_preset(self, name: str, colors: dict[str, str]) -> None:
        if self._nuc_grid is None:
            return
        self._nuc_grid.set_colors(colors)
        # Keep the Annotation sequence-view scheme flag in sync with this preset.
        save_color_blind_nucleotide_scheme(name == "Color-blind")

    def _section_label(self, text: str) -> QLabel:
        label = QLabel(f"<b>{text}</b>")
        label.setTextFormat(Qt.TextFormat.RichText)
        return label

    def restore_defaults(self) -> None:
        if self._nuc_grid is not None:
            self._nuc_grid.set_colors(default_nucleotide_backgrounds())
            save_color_blind_nucleotide_scheme(False)
        if self._aa_grid is not None:
            self._aa_grid.set_colors(default_amino_acid_backgrounds())
        defaults = default_feature_type_colors()
        self._feature_grid.set_colors(
            {
                key: defaults.get(key, "#7f8c8d")
                for key in self._feature_grid.colors()
            }
        )
        if self._database_grid is not None:
            db_defaults = {
                name: SIMILARITY_DB_COLORS[index % len(SIMILARITY_DB_COLORS)]
                for index, name in enumerate(self._database_grid.colors())
            }
            self._database_grid.set_colors(db_defaults)

    def new_feature_type_colors(self) -> dict[str, str]:
        return self._feature_grid.colors()

    def new_database_colors(self) -> dict[str, str]:
        if self._database_grid is None:
            return {}
        return self._database_grid.colors()

    def new_nucleotide_backgrounds(self) -> dict[str, str]:
        if self._nuc_grid is None:
            return {}
        return self._nuc_grid.colors()

    def new_amino_acid_backgrounds(self) -> dict[str, str]:
        if self._aa_grid is None:
            return {}
        return self._aa_grid.colors()

    def old_feature_type_colors(self) -> dict[str, str]:
        return dict(self._old_feature_colors)

    def old_database_colors(self) -> dict[str, str]:
        return dict(self._old_database_colors)

    def save(self) -> None:
        save_feature_type_colors(self.new_feature_type_colors())
        save_similarity_database_colors(self.new_database_colors())
        if self._show_residue_colors:
            save_custom_nucleotide_backgrounds(self.new_nucleotide_backgrounds())
            save_custom_amino_acid_backgrounds(self.new_amino_acid_backgrounds())


class AnnotationElementColorsDialog(QDialog):
    """Element colors for nucleotides, amino acids, features, and databases.

    Font size and “split tracks by source” live on the Annotation Display rail.
    """

    def __init__(
        self,
        parent: QWidget | None = None,
        *,
        show_apply_to_sequence: bool = True,
        font_point_size: int = 13,
        split_tracks_by_source: bool = True,
        min_font_pt: int = 6,
        max_font_pt: int = 36,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Annotation settings")
        self.resize(620, 700)

        # Kept for compatibility with older call sites; Display options moved to
        # the Annotation rail. Values are echoed back unchanged.
        _ = show_apply_to_sequence
        _ = min_font_pt
        _ = max_font_pt
        self._font_point_size = int(font_point_size)
        self._split_tracks_by_source = bool(split_tracks_by_source)

        color_box = QGroupBox("Element colors", self)
        color_box_layout = QVBoxLayout(color_box)
        color_box_layout.setContentsMargins(12, 12, 12, 12)
        self._panel = AnnotationElementColorsPanel(color_box)
        color_box_layout.addWidget(self._panel)

        body = QWidget(self)
        body_layout = QVBoxLayout(body)
        body_layout.setContentsMargins(12, 10, 12, 10)
        body_layout.setSpacing(12)
        body_layout.addWidget(color_box)
        body_layout.addStretch(1)

        scroll = QScrollArea(self)
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setWidget(body)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok
            | QDialogButtonBox.StandardButton.Cancel
            | QDialogButtonBox.StandardButton.RestoreDefaults,
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        buttons.button(QDialogButtonBox.StandardButton.RestoreDefaults).clicked.connect(
            self._restore_defaults
        )

        root = QVBoxLayout(self)
        root.setContentsMargins(12, 10, 12, 10)
        root.setSpacing(8)
        root.addWidget(scroll, stretch=1)
        root.addWidget(buttons)

    def _restore_defaults(self) -> None:
        self._panel.restore_defaults()

    def font_point_size(self) -> int:
        return int(self._font_point_size)

    def split_tracks_by_source(self) -> bool:
        return bool(self._split_tracks_by_source)

    def apply_to_current_sequence(self) -> bool:
        return True

    def new_feature_type_colors(self) -> dict[str, str]:
        return self._panel.new_feature_type_colors()

    def new_database_colors(self) -> dict[str, str]:
        return self._panel.new_database_colors()

    def new_nucleotide_backgrounds(self) -> dict[str, str]:
        return self._panel.new_nucleotide_backgrounds()

    def new_amino_acid_backgrounds(self) -> dict[str, str]:
        return self._panel.new_amino_acid_backgrounds()

    def old_feature_type_colors(self) -> dict[str, str]:
        return self._panel.old_feature_type_colors()

    def old_database_colors(self) -> dict[str, str]:
        return self._panel.old_database_colors()

    def accept(self) -> None:
        self._panel.save()
        super().accept()


class SimilarityDatabaseColorsDialog(QDialog):
    """Edit colors for similarity annotation databases only."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Database colors")
        self.resize(360, 320)

        self._database_rows: dict[str, _ColorPickerRow] = {}
        databases = load_similarity_annotation_databases()

        root = QVBoxLayout(self)
        scroll = QScrollArea(self)
        scroll.setWidgetResizable(True)
        body = QWidget(scroll)
        body_layout = QVBoxLayout(body)
        body_layout.setSpacing(8)

        database_form = QFormLayout()
        database_form.setSpacing(8)
        if databases:
            for index, entry in enumerate(databases):
                row = _ColorPickerRow(
                    entry.name,
                    similarity_database_color(entry.name, index),
                    body,
                )
                self._database_rows[entry.name] = row
                database_form.addRow(row)
        else:
            database_form.addRow(
                QLabel("No databases defined. Add them in File → Manage Database.", body)
            )
        body_layout.addLayout(database_form)
        body_layout.addStretch(1)
        scroll.setWidget(body)
        root.addWidget(scroll, stretch=1)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok
            | QDialogButtonBox.StandardButton.Cancel
            | QDialogButtonBox.StandardButton.RestoreDefaults,
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        buttons.button(QDialogButtonBox.StandardButton.RestoreDefaults).clicked.connect(
            self._restore_defaults
        )
        root.addWidget(buttons)

    def _restore_defaults(self) -> None:
        for index, name in enumerate(self._database_rows):
            self._database_rows[name].set_color_hex(
                SIMILARITY_DB_COLORS[index % len(SIMILARITY_DB_COLORS)]
            )

    def new_database_colors(self) -> dict[str, str]:
        return {name: row.color_hex() for name, row in self._database_rows.items()}

    def accept(self) -> None:
        if self._database_rows:
            save_similarity_database_colors(self.new_database_colors())
        super().accept()


class AddFeatureButtonsDialog(QDialog):
    """Create, delete, and reorder buttons on the Add feature row."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Manage feature buttons")
        self.resize(420, 380)

        self._list = QListWidget(self)
        self._list.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self._list.setDragDropMode(QAbstractItemView.DragDropMode.InternalMove)
        self._list.setDefaultDropAction(Qt.DropAction.MoveAction)
        self._list.model().rowsMoved.connect(self._renumber_items)

        add_btn = QPushButton("Add…", self)
        add_btn.clicked.connect(self._add_button)
        delete_btn = QPushButton("Delete", self)
        delete_btn.clicked.connect(self._delete_selected)
        up_btn = QPushButton("Move up", self)
        up_btn.clicked.connect(lambda: self._move_selected(-1))
        down_btn = QPushButton("Move down", self)
        down_btn.clicked.connect(lambda: self._move_selected(1))

        side = QVBoxLayout()
        side.setSpacing(6)
        for btn in (add_btn, delete_btn, up_btn, down_btn):
            side.addWidget(btn)
        side.addStretch(1)

        body = QHBoxLayout()
        body.addWidget(self._list, stretch=1)
        body.addLayout(side)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok
            | QDialogButtonBox.StandardButton.Cancel
            | QDialogButtonBox.StandardButton.RestoreDefaults,
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        restore = buttons.button(QDialogButtonBox.StandardButton.RestoreDefaults)
        restore.setText("Restore")
        restore.clicked.connect(self._restore_defaults)

        root = QVBoxLayout(self)
        root.addLayout(body, stretch=1)
        root.addWidget(buttons)

        self._set_buttons(load_add_feature_buttons())

    def _item_text(self, index: int, label: str, feature_type: str) -> str:
        return f"{index}. {label}  ({feature_type})"

    def _renumber_items(self, *_args) -> None:
        for index in range(self._list.count()):
            item = self._list.item(index)
            data = item.data(Qt.ItemDataRole.UserRole)
            if not isinstance(data, (list, tuple)) or len(data) < 2:
                continue
            label, feature_type = str(data[0]), str(data[1])
            item.setText(self._item_text(index + 1, label, feature_type))

    def _set_buttons(self, buttons: list[tuple[str, str]]) -> None:
        self._list.clear()
        for index, (label, feature_type) in enumerate(buttons, start=1):
            item = QListWidgetItem(self._item_text(index, label, feature_type))
            item.setData(Qt.ItemDataRole.UserRole, (label, feature_type))
            self._list.addItem(item)

    def _current_buttons(self) -> list[tuple[str, str]]:
        buttons: list[tuple[str, str]] = []
        for index in range(self._list.count()):
            item = self._list.item(index)
            data = item.data(Qt.ItemDataRole.UserRole)
            if isinstance(data, (list, tuple)) and len(data) >= 2:
                buttons.append((str(data[0]), str(data[1])))
        return buttons

    def _add_button(self) -> None:
        dialog = QDialog(self)
        dialog.setWindowTitle("Add feature button")
        label_edit = QLineEdit(dialog)
        name_edit = QLineEdit(dialog)

        def _sync_name_from_label(text: str) -> None:
            if name_edit.isModified():
                return
            name_edit.blockSignals(True)
            name_edit.setText(slug_feature_type(text))
            name_edit.blockSignals(False)

        label_edit.textChanged.connect(_sync_name_from_label)
        form = QFormLayout()
        form.addRow("Button label:", label_edit)
        form.addRow("Names stored in GFF file:", name_edit)
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel,
            parent=dialog,
        )
        buttons.accepted.connect(dialog.accept)
        buttons.rejected.connect(dialog.reject)
        root = QVBoxLayout(dialog)
        root.addLayout(form)
        root.addWidget(buttons)
        label_edit.setFocus()
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        label = label_edit.text().strip()
        if not label:
            QMessageBox.warning(self, "Add feature button", "Enter a button label.")
            return
        feature_type = name_edit.text().strip() or slug_feature_type(label)
        existing = {ftype.lower() for _label, ftype in self._current_buttons()}
        if feature_type.lower() in existing:
            QMessageBox.warning(
                self,
                "Add feature button",
                f"A button for “{feature_type}” already exists.",
            )
            return
        index = self._list.count() + 1
        item = QListWidgetItem(self._item_text(index, label, feature_type))
        item.setData(Qt.ItemDataRole.UserRole, (label, feature_type))
        self._list.addItem(item)
        self._list.setCurrentItem(item)

    def _delete_selected(self) -> None:
        row = self._list.currentRow()
        if row < 0:
            return
        if self._list.count() <= 1:
            QMessageBox.warning(
                self,
                "Delete feature button",
                "Keep at least one Add feature button.",
            )
            return
        self._list.takeItem(row)
        self._renumber_items()

    def _move_selected(self, delta: int) -> None:
        row = self._list.currentRow()
        if row < 0:
            return
        target = row + delta
        if target < 0 or target >= self._list.count():
            return
        item = self._list.takeItem(row)
        self._list.insertItem(target, item)
        self._list.setCurrentRow(target)
        self._renumber_items()

    def _restore_defaults(self) -> None:
        self._set_buttons(list(DEFAULT_ADD_FEATURE_BUTTONS))

    def buttons(self) -> list[tuple[str, str]]:
        return self._current_buttons()

    def accept(self) -> None:
        buttons = self._current_buttons()
        if not buttons:
            QMessageBox.warning(
                self,
                "Manage feature buttons",
                "Keep at least one Add feature button.",
            )
            return
        save_add_feature_buttons(buttons)
        super().accept()
