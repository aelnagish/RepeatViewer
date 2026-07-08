"""Dialog for editing annotation element colors."""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (
    QDialog,
    QDialogButtonBox,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QScrollArea,
    QVBoxLayout,
    QWidget,
    QColorDialog,
)

from teviewer.backends.annotation_colors import (
    ANNOTATION_FEATURE_TYPE_LABELS,
    default_feature_type_colors,
    feature_type_color,
    load_feature_type_colors,
    save_feature_type_colors,
    save_similarity_database_colors,
    similarity_database_color,
    snapshot_similarity_database_colors,
)
from teviewer.backends.similarity_annotation import SIMILARITY_DB_COLORS
from teviewer.database_settings import load_similarity_annotation_databases
from teviewer.widgets.preview_zoom import (
    default_nucleotide_backgrounds,
    load_custom_nucleotide_backgrounds,
    save_custom_nucleotide_backgrounds,
)


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
        self._button.setStyleSheet(
            f"QPushButton {{ background-color: {hex_color}; border: 1px solid palette(mid); }}"
        )

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


class AnnotationElementColorsDialog(QDialog):
    """Edit colors for manual/auto feature types and similarity databases."""

    def __init__(
        self,
        parent: QWidget | None = None,
        *,
        show_apply_to_sequence: bool = True,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Element colors")
        self.resize(420, 520)

        self._feature_rows: dict[str, _ColorPickerRow] = {}
        self._nucleotide_rows: dict[str, _ColorPickerRow] = {}
        self._database_rows: dict[str, _ColorPickerRow] = {}
        databases = load_similarity_annotation_databases()

        root = QVBoxLayout(self)

        scroll = QScrollArea(self)
        scroll.setWidgetResizable(True)
        body = QWidget(scroll)
        body_layout = QVBoxLayout(body)
        body_layout.setSpacing(10)

        body_layout.addWidget(self._section_label("Feature types"))
        feature_form = QFormLayout()
        feature_form.setSpacing(8)
        for display_label, feature_type in ANNOTATION_FEATURE_TYPE_LABELS:
            row = _ColorPickerRow(display_label, feature_type_color(feature_type), body)
            self._feature_rows[feature_type] = row
            feature_form.addRow(row)
        body_layout.addLayout(feature_form)

        body_layout.addWidget(self._section_label("Nucleotides"))
        nucleotide_form = QFormLayout()
        nucleotide_form.setSpacing(8)
        custom_nucleotide = load_custom_nucleotide_backgrounds()
        scheme_nucleotide = default_nucleotide_backgrounds()
        for base in ("A", "C", "G", "T"):
            color = custom_nucleotide.get(base, scheme_nucleotide[base])
            row = _ColorPickerRow(base, color, body)
            self._nucleotide_rows[base] = row
            nucleotide_form.addRow(row)
        body_layout.addLayout(nucleotide_form)

        body_layout.addWidget(self._section_label("Similarity databases"))
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
                QLabel("No databases defined. Add them in File → Define database.", body)
            )
        body_layout.addLayout(database_form)
        body_layout.addStretch(1)
        scroll.setWidget(body)
        root.addWidget(scroll, stretch=1)

        # Kept for compatibility with older call sites; this dialog no longer
        # exposes a per-sequence apply toggle and always applies to current sequence.
        _ = show_apply_to_sequence

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

        self._old_feature_colors = load_feature_type_colors()
        self._old_database_colors = snapshot_similarity_database_colors(
            [entry.name for entry in databases]
        )

    def _section_label(self, text: str) -> QLabel:
        label = QLabel(f"<b>{text}</b>")
        label.setTextFormat(Qt.TextFormat.RichText)
        return label

    def _restore_defaults(self) -> None:
        defaults = default_feature_type_colors()
        for feature_type, row in self._feature_rows.items():
            row.set_color_hex(defaults.get(feature_type, "#7f8c8d"))
        for index, name in enumerate(self._database_rows):
            row = self._database_rows[name]
            row.set_color_hex(SIMILARITY_DB_COLORS[index % len(SIMILARITY_DB_COLORS)])
        for base, row in self._nucleotide_rows.items():
            row.set_color_hex(default_nucleotide_backgrounds()[base])

    def apply_to_current_sequence(self) -> bool:
        return True

    def new_feature_type_colors(self) -> dict[str, str]:
        return {key: row.color_hex() for key, row in self._feature_rows.items()}

    def new_database_colors(self) -> dict[str, str]:
        return {name: row.color_hex() for name, row in self._database_rows.items()}

    def new_nucleotide_backgrounds(self) -> dict[str, str]:
        return {base: row.color_hex() for base, row in self._nucleotide_rows.items()}

    def old_feature_type_colors(self) -> dict[str, str]:
        return dict(self._old_feature_colors)

    def old_database_colors(self) -> dict[str, str]:
        return dict(self._old_database_colors)

    def accept(self) -> None:
        save_feature_type_colors(self.new_feature_type_colors())
        save_similarity_database_colors(self.new_database_colors())
        save_custom_nucleotide_backgrounds(self.new_nucleotide_backgrounds())
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
                QLabel("No databases defined. Add them in File → Define database.", body)
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
