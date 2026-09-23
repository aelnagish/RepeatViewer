"""Default BEEA workflow naming and step normalization."""

from __future__ import annotations

from repeatviewer.workflows import (
    WorkflowStep,
    _normalize_beea_steps,
    default_beea_workflow,
    ensure_default_workflows,
)


def test_default_beea_name_and_steps() -> None:
    spec = default_beea_workflow()
    assert spec.id == "bee"
    assert spec.name == "Blast-Extract-Extend-Align (BEEA)"
    assert [step.tool_id for step in spec.steps] == ["blast", "extend", "align"]


def test_normalize_beea_appends_align() -> None:
    steps = _normalize_beea_steps(
        [WorkflowStep(tool_id="blast"), WorkflowStep(tool_id="extend")]
    )
    assert [step.tool_id for step in steps] == ["blast", "extend", "align"]


def test_ensure_default_workflows_renames_legacy_bee(monkeypatch) -> None:
    stored: list[dict] = [
        {
            "id": "bee",
            "name": "BEE",
            "steps": [{"tool_id": "blast"}, {"tool_id": "extend"}],
        },
        {"id": "other", "name": "Other", "steps": [{"tool_id": "align"}]},
    ]
    written: list = []

    monkeypatch.setattr(
        "repeatviewer.workflows.load_workflow_specs",
        lambda: [
            __import__("repeatviewer.workflows", fromlist=["WorkflowSpec"]).WorkflowSpec.from_json(
                item
            )
            for item in stored
        ],
    )
    monkeypatch.setattr(
        "repeatviewer.workflows._write_all",
        lambda specs: written.append([s.to_json() for s in specs]),
    )

    specs = ensure_default_workflows()
    assert specs[0].name == "Blast-Extract-Extend-Align (BEEA)"
    assert [step.tool_id for step in specs[0].steps] == ["blast", "extend", "align"]
    assert specs[0].id == "bee"
    assert written, "expected migration to persist"
