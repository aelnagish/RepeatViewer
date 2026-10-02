"""Edit > Batch Rename window: live preview and the Rename button."""

from __future__ import annotations

from pathlib import Path


def _files(folder: Path, *names: str) -> list[Path]:
    paths = []
    for name in names:
        path = folder / name
        path.write_text(">x\nACGT\n")
        paths.append(path.resolve())
    return paths


def test_preview_and_rename_button(qapp, tmp_path):
    from repeatviewer.views.batch_rename_dialog import BatchRenameDialog

    paths = _files(tmp_path, "contig_1.fa", "contig_2.fa", "Chr_2.fa")
    (tmp_path / "contig_1.seqmap_RV.json").write_text("{}")
    dialog = BatchRenameDialog(paths)
    assert not dialog._rename_btn.isEnabled()          # nothing to rename yet

    dialog._find.setText("contig")
    dialog._replace.setText("Chr")
    plan = dialog.plan()                               # flushes the typing pause
    rows = {row.source.name: row for row in plan.rows}
    assert rows["contig_2.fa"].problem == "a file with this name already exists"
    assert not dialog._rename_btn.isEnabled()
    assert "already exists" in dialog._status.text()

    dialog._prefix.setText("Cs_")
    dialog._find.setText("contig_")
    plan = dialog.plan()
    assert plan.ok and dialog._rename_btn.isEnabled()
    assert sorted(new.name for _old, new in plan.files) == [
        "Cs_Chr1.fa",
        "Cs_Chr2.fa",
        "Cs_Chr_2.fa",
    ]
    moves = dict((old.name, new.name) for old, new in plan.moves)
    assert moves["contig_1.seqmap_RV.json"] == "Cs_Chr1.seqmap_RV.json"


def test_broken_regular_expression_is_shown(qapp, tmp_path):
    from repeatviewer.views.batch_rename_dialog import BatchRenameDialog

    dialog = BatchRenameDialog(_files(tmp_path, "a.fa"))
    dialog._regex.setChecked(True)
    dialog._find.setText("(")
    plan = dialog.plan()
    assert plan.error and not dialog._rename_btn.isEnabled()
    assert dialog._status.text().startswith("Regular expression")
