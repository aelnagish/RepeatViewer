"""Tool Options → Citation strings (Zotero RepeatViewer → Tools)."""

from __future__ import annotations

from repeatviewer.tool_citations import (
    CITATION_ALIGN,
    CITATION_BLAST,
    CITATION_DOTPLOT,
    CITATION_MAPPING,
    CITATION_SEQUENCE_ANNOTATION,
    CITATION_TEAID,
    CITATION_TREE,
)


def test_blast_citation_is_blast_plus_only() -> None:
    assert "BLAST+" in CITATION_BLAST
    assert "10.1186/1471-2105-10-421" in CITATION_BLAST
    assert "Ye J" not in CITATION_BLAST
    assert "10.1093/nar/gkl361" not in CITATION_BLAST


def test_align_citations_mafft_muscle_fastga() -> None:
    assert "MAFFT" in CITATION_ALIGN and "10.1093/molbev/mst010" in CITATION_ALIGN
    assert "MUSCLE" in CITATION_ALIGN and "10.1093/nar/gkh340" in CITATION_ALIGN
    assert "FastGA" in CITATION_ALIGN and "10.1093/bioadv/vbaf238" in CITATION_ALIGN


def test_dotplot_flexidot() -> None:
    assert "FlexiDot" in CITATION_DOTPLOT
    assert "10.1093/bioinformatics/bty395" in CITATION_DOTPLOT


def test_teaid_goubert_guide() -> None:
    assert "beginner" in CITATION_TEAID.lower() or "Beginner" in CITATION_TEAID
    assert "10.1186/s13100-022-00267-1" in CITATION_TEAID


def test_tree_iqtree_fastme_fasttree() -> None:
    assert "IQ-TREE 2" in CITATION_TREE and "10.1093/molbev/msaa015" in CITATION_TREE
    assert "FastME" in CITATION_TREE and "10.1093/molbev/msv150" in CITATION_TREE
    assert "FastTree" in CITATION_TREE and "10.1093/molbev/msp077" in CITATION_TREE


def test_mapping_minimap2_hisat() -> None:
    assert "Minimap2" in CITATION_MAPPING and "10.1093/bioinformatics/bty191" in CITATION_MAPPING
    assert "HISAT" in CITATION_MAPPING and "10.1038/nmeth.3317" in CITATION_MAPPING


def test_sequence_annotation_orf_tandem_tetrimmer() -> None:
    assert "orfipy" in CITATION_SEQUENCE_ANNOTATION
    assert "10.1093/bioinformatics/btab090" in CITATION_SEQUENCE_ANNOTATION
    assert "Pytrf" in CITATION_SEQUENCE_ANNOTATION
    assert "10.1186/s12859-025-06168-3" in CITATION_SEQUENCE_ANNOTATION
    assert "TEtrimmer" in CITATION_SEQUENCE_ANNOTATION
    assert "10.1038/s41467-025-63889-y" in CITATION_SEQUENCE_ANNOTATION


def test_make_citation_box_smoke(qapp) -> None:
    from PySide6.QtWidgets import QWidget

    from repeatviewer.views.tool_dialogs import make_citation_box

    parent = QWidget()
    box = make_citation_box(parent, CITATION_BLAST)
    assert box.title() == "Citation"
    label = box.findChild(type(box), "ToolCitationBody")
    # Object name is on the QLabel inside the box.
    from PySide6.QtWidgets import QLabel

    labels = box.findChildren(QLabel)
    bodies = [lb for lb in labels if lb.objectName() == "ToolCitationBody"]
    assert bodies and "BLAST+" in bodies[0].text()
