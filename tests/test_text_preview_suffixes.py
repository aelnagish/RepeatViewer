"""Text-preview recognition for common bioinformatics sidecar files."""

from pathlib import Path

from repeatviewer.widgets.preview_tabs import (
    _is_text_like_preview,
    is_previewable_file,
)


def test_chrom_sizes_file_is_previewed_as_text() -> None:
    path = Path("/tmp/KWS23200NT_v1.0.fasta.chrom.sizes")

    assert is_previewable_file(path)
    assert _is_text_like_preview(path)
