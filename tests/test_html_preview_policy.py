import os
from pathlib import Path

import pytest

from repeatviewer.backends import html_preview_policy as policy
from repeatviewer.backends import webengine_probe


def _fake_webengine_is_safe(*, force_probe: bool = False) -> bool:
    """Stand-in for the real probe: honors the env override, skips the subprocess."""
    del force_probe
    flag = os.environ.get("REPEATVIEWER_EMBED_HTML", "").strip().lower()
    if flag in {"0", "false", "no", "off"}:
        return False
    return True


@pytest.fixture(autouse=True)
def _clear_html_env(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("REPEATVIEWER_EMBED_HTML", raising=False)
    monkeypatch.delenv("REPEATVIEWER_FLEXIDOT_WEBENGINE", raising=False)
    # Policy defers to the real crash probe when no env override is set; pin it
    # to "safe" (unless the test sets REPEATVIEWER_EMBED_HTML=0 itself) so these
    # tests exercise the policy logic, not this machine's actual WebEngine
    # health or a slow subprocess probe.
    monkeypatch.setattr(webengine_probe, "webengine_is_safe", _fake_webengine_is_safe)


def test_flexidot_html_uses_webengine_by_default() -> None:
    path = Path("/tmp/gene_flexidot.html")
    assert policy.should_use_webengine_for_html(path) is True


def test_flexidot_html_can_opt_out_to_browser(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("REPEATVIEWER_FLEXIDOT_WEBENGINE", "0")
    path = Path("/tmp/gene_flexidot.html")
    assert policy.should_use_webengine_for_html(path) is False


def test_tea_id_html_still_uses_webengine() -> None:
    path = Path("/tmp/teaid_report.html")
    assert policy.should_use_webengine_for_html(path) is True


def test_embed_html_disabled_blocks_all(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("REPEATVIEWER_EMBED_HTML", "0")
    assert policy.should_use_webengine_for_html(Path("/tmp/teaid_report.html")) is False
    assert policy.should_use_webengine_for_html(Path("/tmp/gene_flexidot.html")) is False


def test_unsafe_webengine_probe_does_not_block_embedded_preview(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Background probe failure must not force browser-only mode by default."""
    monkeypatch.setattr(webengine_probe, "webengine_is_safe", lambda **_: False)
    assert policy.should_use_webengine_for_html(Path("/tmp/teaid_report.html")) is True
    assert policy.may_attempt_embedded_html_preview() is True


def test_embed_html_env_opt_out_blocks_preview(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("REPEATVIEWER_EMBED_HTML", "0")
    assert policy.may_attempt_embedded_html_preview() is False
    assert policy.should_use_webengine_for_html(Path("/tmp/teaid_report.html")) is False
    assert "REPEATVIEWER_EMBED_HTML=0" in policy.embedded_html_disabled_reason()
