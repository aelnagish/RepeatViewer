from pathlib import Path

from PySide6.QtCore import QUrl

from teviewer.widgets.html_preview import _local_file_url, _webengine_safe_html_path


def test_local_file_url_encodes_hash_in_filename() -> None:
    path = Path("/tmp/AmaS-Au_Beta_nana#Class_I__SINE_flexidot.html")
    url = _local_file_url(path)
    assert isinstance(url, QUrl)
    assert url.isLocalFile()
    assert url.fragment() == ""
    assert "%23" in url.toString()
    assert url.toLocalFile().endswith("AmaS-Au_Beta_nana#Class_I__SINE_flexidot.html")


def test_webengine_safe_html_path_strips_hash(tmp_path: Path) -> None:
    src = tmp_path / "Ama5-lic_Beta_lomatogona#Class_I__SINE.fa_TEAid.html"
    src.write_text("<html><body>teaid</body></html>", encoding="utf-8")
    safe = _webengine_safe_html_path(src)
    assert safe != src
    assert "#" not in safe.name
    assert safe.is_file()
    assert "%23" not in _local_file_url(safe).toString()
