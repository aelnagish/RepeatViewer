"""Export Annotation JSON to GFF3 / GTF / BED."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QScrollArea,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends.sequence_features import (
    duplicate_original_sequence_ids,
    export_features_bed_multi,
    export_features_gff_multi,
    export_features_gtf_multi,
    features_gff_path,
    load_features_by_id_from_json,
    load_sequence_map_data,
    resolve_sidecar_path,
    suggest_unique_export_ids,
)
from repeatviewer.widgets.file_conflict import confirm_existing_file_path

_FORMAT_SUFFIX = {
    "gff3": ".gff3",
    "gtf": ".gtf",
    "bed": ".bed",
}


def fasta_has_exportable_annotations(fasta_path: Path) -> bool:
    """True when a JSON sidecar exists with at least one feature."""
    try:
        resolved = fasta_path.expanduser().resolve()
    except OSError:
        return False
    if resolve_sidecar_path(resolved) is None:
        return False
    try:
        data = load_sequence_map_data(resolved)
    except Exception:  # noqa: BLE001
        return False
    record_ids = list(data.record_ids)
    json_feats = load_features_by_id_from_json(resolved, record_ids)
    feats = json_feats if json_feats and any(json_feats.values()) else data.features_by_id
    return bool(feats and any(feats.values()))


class _DuplicateSeqidDialog(QDialog):
    """Edit export seqids when FASTA headers are not unique."""

    def __init__(
        self,
        record_ids: list[str],
        suggestions: list[str],
        duplicates: list[str],
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Modify sequence names for export")
        self.setMinimumWidth(420)
        layout = QVBoxLayout(self)
        layout.addWidget(
            QLabel(
                "Some sequence names appear more than once. "
                "GFF3 / GTF / BED need unique names — edit the suggestions below "
                f"(duplicates: {', '.join(duplicates)})."
            )
        )
        scroll = QScrollArea(self)
        scroll.setWidgetResizable(True)
        host = QWidget(scroll)
        host_layout = QVBoxLayout(host)
        self._edits: list[QLineEdit] = []
        for idx, (internal, suggested) in enumerate(
            zip(record_ids, suggestions, strict=True)
        ):
            row = QHBoxLayout()
            row.addWidget(QLabel(f"#{idx + 1} {internal}:"))
            edit = QLineEdit(suggested, host)
            self._edits.append(edit)
            row.addWidget(edit, stretch=1)
            wrap = QWidget(host)
            wrap.setLayout(row)
            host_layout.addWidget(wrap)
        host_layout.addStretch()
        scroll.setWidget(host)
        scroll.setMaximumHeight(280)
        layout.addWidget(scroll)
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel,
            parent=self,
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)
        self._record_ids = list(record_ids)

    def seqid_map(self) -> dict[str, str]:
        mapping: dict[str, str] = {}
        used: set[str] = set()
        for rid, edit in zip(self._record_ids, self._edits, strict=True):
            name = edit.text().strip()
            if not name:
                raise ValueError("Export names cannot be empty.")
            if name in used:
                raise ValueError(f"Export name “{name}” is still duplicated.")
            used.add(name)
            mapping[rid] = name
        return mapping


class ExportAnnotationDialog(QDialog):
    """Format + editable filename; conflict uses overwrite/rename prompt."""

    def __init__(self, fasta_path: Path, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._fasta = fasta_path.expanduser().resolve()
        self.setWindowTitle("Export annotation")
        self.setMinimumWidth(360)

        layout = QVBoxLayout(self)
        form = QFormLayout()
        self._format = QComboBox(self)
        self._format.addItem("GFF3", "gff3")
        self._format.addItem("GTF", "gtf")
        self._format.addItem("BED", "bed")
        form.addRow("Format:", self._format)
        self._filename = QLineEdit(self)
        form.addRow("File name:", self._filename)
        layout.addLayout(form)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel,
            parent=self,
        )
        buttons.accepted.connect(self._accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

        self._record_ids: list[str] = []
        self._features_by_id: dict = {}
        self._lengths_by_id: dict = {}
        self._written: Path | None = None
        self._syncing_name = False
        self._load()
        self._sync_filename_from_format()
        self._format.currentIndexChanged.connect(self._on_format_changed)

    def _load(self) -> None:
        data = load_sequence_map_data(self._fasta)
        self._record_ids = list(data.record_ids)
        self._features_by_id = dict(data.features_by_id)
        self._lengths_by_id = dict(data.lengths_by_id)
        json_feats = load_features_by_id_from_json(self._fasta, self._record_ids)
        if json_feats and any(json_feats.values()):
            self._features_by_id = json_feats

    def _default_out(self, fmt: str) -> Path:
        if fmt == "gff3":
            return features_gff_path(self._fasta)
        if fmt == "gtf":
            return self._fasta.parent / f"{self._fasta.stem}.RV.gtf"
        return self._fasta.parent / f"{self._fasta.stem}.RV.bed"

    def _current_fmt(self) -> str:
        return str(self._format.currentData() or "gff3")

    def _sync_filename_from_format(self) -> None:
        self._syncing_name = True
        try:
            self._filename.setText(self._default_out(self._current_fmt()).name)
        finally:
            self._syncing_name = False

    def _on_format_changed(self, _index: int = 0) -> None:
        if self._syncing_name:
            return
        current = self._filename.text().strip()
        suffix = _FORMAT_SUFFIX[self._current_fmt()]
        if current:
            name = Path(current).with_suffix(suffix).name
            self._syncing_name = True
            try:
                self._filename.setText(name)
            finally:
                self._syncing_name = False
        else:
            self._sync_filename_from_format()

    def _proposed_path(self) -> Path | None:
        name = self._filename.text().strip()
        if not name:
            QMessageBox.warning(self, "Export annotation", "Enter a file name.")
            return None
        if Path(name).name != name or "/" in name or "\\" in name:
            QMessageBox.warning(
                self,
                "Export annotation",
                "Enter a file name only (no folder path).",
            )
            return None
        fmt = self._current_fmt()
        suffix = _FORMAT_SUFFIX[fmt]
        path = (self._fasta.parent / name).expanduser().resolve()
        if path.suffix.lower() != suffix:
            path = path.with_suffix(suffix)
        return path

    def _resolve_seqid_map(self) -> dict[str, str] | None:
        originals = [rid.split("#", 1)[0] for rid in self._record_ids]
        dups = duplicate_original_sequence_ids(originals)
        if not dups:
            return None
        suggestions = suggest_unique_export_ids(originals)
        dialog = _DuplicateSeqidDialog(
            self._record_ids, suggestions, dups, parent=self
        )
        if dialog.exec() != QDialog.DialogCode.Accepted:
            raise RuntimeError("cancelled")
        return dialog.seqid_map()

    def _accept(self) -> None:
        if not self._features_by_id or not any(self._features_by_id.values()):
            QMessageBox.information(
                self, "Export annotation", "No features to export."
            )
            return
        out = self._proposed_path()
        if out is None:
            return
        fmt = self._current_fmt()
        suffix = _FORMAT_SUFFIX[fmt]
        confirmed = confirm_existing_file_path(
            self,
            out,
            title="Export annotation",
            informative=(
                "Overwrite the existing file, or save under a different name?"
            ),
            required_suffix=suffix,
        )
        if confirmed is None:
            return
        try:
            seqid_map = self._resolve_seqid_map()
        except RuntimeError:
            return
        except ValueError as exc:
            QMessageBox.warning(self, "Export annotation", str(exc))
            return
        try:
            if fmt == "gff3":
                export_features_gff_multi(
                    self._fasta,
                    features_by_id=self._features_by_id,
                    lengths_by_id=self._lengths_by_id,
                    output_path=confirmed,
                    seqid_map=seqid_map,
                )
            elif fmt == "gtf":
                export_features_gtf_multi(
                    self._fasta,
                    features_by_id=self._features_by_id,
                    lengths_by_id=self._lengths_by_id,
                    output_path=confirmed,
                    seqid_map=seqid_map,
                )
            else:
                export_features_bed_multi(
                    self._fasta,
                    features_by_id=self._features_by_id,
                    lengths_by_id=self._lengths_by_id,
                    output_path=confirmed,
                    seqid_map=seqid_map,
                )
        except OSError as exc:
            QMessageBox.critical(self, "Export annotation", str(exc))
            return
        self._written = confirmed
        self.accept()

    def written_path(self) -> Path | None:
        return self._written
