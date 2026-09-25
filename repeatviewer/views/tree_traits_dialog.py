"""Tree tab → Trait table → Create or update: settings for pre-filling the table.

Follows the RepeatViewer dialog rules: Cancel / Save / Run, the file-path field
in the shared path style (click a file in the Files panel to fill it), and the
round "?" help.
"""

from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import (
    QDialog,
    QFormLayout,
    QLabel,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends import tree_traits as tt
from repeatviewer.views.tool_dialogs import QSpinBox, dialog_should_run, make_cancel_save_run_bar
from repeatviewer.widgets.help_button import field_with_help
from repeatviewer.widgets.path_browse import PathBrowseRow

INTRO = (
    "Fills a table next to the tree with what RepeatViewer knows about each sequence: "
    "position, strand and tandem array (from the FASTA headers), methylation (from the "
    "file below, optional) and, for trees of several FASTA files, the source file and "
    "whether FISH images are attached to it. Columns you added yourself are kept."
)
METHYLATION_HELP = (
    "<b>Methylation file</b>: a Bismark coverage file (<b>.cov</b> or <b>.cov.gz</b>) or a "
    "modkit <b>bedMethyl</b> file. For each sequence, the % of methylated reads over its "
    "position is written to the <b>methylation</b> column.<br>"
    "Needs positions in the headers (chrom:start-end), as RepeatViewer writes them when it "
    "extracts copies from a genome. The chromosome names must match the file.<br>"
    "Leave empty to skip."
)
ARRAY_HELP = (
    "<b>Array gap</b>: copies on the same chromosome closer than this belong to one tandem "
    "array. Each array is named by its span, e.g. Chr1:12.40-12.91 Mb."
)
_METH_SUFFIXES = (".cov", ".cov.gz", ".bed", ".bed.gz", ".bedmethyl", ".bedmethyl.gz",
                  ".txt", ".txt.gz", ".tsv", ".tsv.gz")


class TraitTableDialog(QDialog):
    def __init__(self, *, settings: dict | None = None, tree_path: Path | None = None,
                 parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Trait table")
        settings = settings or {}
        intro = QLabel(INTRO, self)
        intro.setWordWrap(True)
        self._meth = PathBrowseRow(
            self,
            initial=settings.get("methylation_file") or None,
            placeholder="Optional: click a .cov or bedMethyl file in Files, or Browse",
            browse_title="Choose a methylation file",
            file_filter="Methylation (*.cov *.cov.gz *.bed *.bed.gz *.bedmethyl *.bedmethyl.gz);;All files (*)",
            start_dir=str(tree_path.parent) if tree_path else "",
        )
        self._meth.set_default_tooltip("Click a file in the Files panel to fill this path, or use Browse.")
        self._gap = QSpinBox(self)
        self._gap.setRange(1, 10_000)
        self._gap.setSuffix(" kb")
        self._gap.setValue(max(1, int(settings.get("array_gap", tt.DEFAULT_ARRAY_GAP)) // 1000))
        self._gap.setFixedWidth(110)

        form = QFormLayout()
        form.setSpacing(10)
        form.addRow("Methylation file:", field_with_help(self._meth, "Methylation file",
                                                          METHYLATION_HELP, parent=self))
        form.addRow("Array gap:", field_with_help(self._gap, "Array gap", ARRAY_HELP,
                                                   parent=self, stretch=True))
        target = tt.traits_path(tree_path).name if tree_path else "…_traits.tsv"
        where = QLabel(f"Writes: {target}", self)
        where.setStyleSheet("color: palette(mid);")

        footer = make_cancel_save_run_bar(self)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(18, 16, 18, 14)
        layout.setSpacing(12)
        layout.addWidget(intro)
        layout.addLayout(form)
        layout.addWidget(where)
        layout.addStretch(1)
        layout.addWidget(footer)
        self.resize(560, 300)

    def settings(self) -> dict:
        path = self._meth.path()
        return {
            "methylation_file": str(path) if path else "",
            "array_gap": self._gap.value() * 1000,
        }

    def should_run(self) -> bool:
        return dialog_should_run(self)

    def apply_path_from_panel(self, path: Path) -> bool:
        """Fill the methylation file from a Files-panel click."""
        resolved = Path(path).expanduser()
        if not resolved.is_file() or not resolved.name.lower().endswith(_METH_SUFFIXES):
            return False
        self._meth.set_path(resolved)
        return True
