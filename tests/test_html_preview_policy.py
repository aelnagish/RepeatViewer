from pathlib import Path

import pytest

from teviewer.backends import html_preview_policy as policy


@pytest.fixture(autouse=True)
def _clear_html_env(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("TEVIEWER_EMBED_HTML", raising=False)
    monkeypatch.delenv("TEVIEWER_FLEXIDOT_WEBENGINE", raising=False)


def test_flexidot_html_uses_webengine_by_default() -> None:
    path = Path("/tmp/gene_flexidot.html")
    assert policy.should_use_webengine_for_html(path) is True


def test_flexidot_html_can_opt_out_to_browser(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("TEVIEWER_FLEXIDOT_WEBENGINE", "0")
    path = Path("/tmp/gene_flexidot.html")
    assert policy.should_use_webengine_for_html(path) is False


def test_tea_id_html_still_uses_webengine() -> None:
    path = Path("/tmp/teaid_report.html")
    assert policy.should_use_webengine_for_html(path) is True


def test_embed_html_disabled_blocks_all(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("TEVIEWER_EMBED_HTML", "0")
    assert policy.should_use_webengine_for_html(Path("/tmp/teaid_report.html")) is False
    assert policy.should_use_webengine_for_html(Path("/tmp/gene_flexidot.html")) is False
