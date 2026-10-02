"""Cytogenetic image store: attach, query, rename, remove, caption."""

from __future__ import annotations

from pathlib import Path

from repeatviewer.backends.cyto_images import (
    STORE_DIR_NAME,
    CytoImage,
    FeatureRef,
    ImageTarget,
    add_image,
    build_caption,
    feature_keys_with_images,
    image_count_for_sequence,
    image_file_path,
    images_for_feature,
    images_for_sequence,
    load_images,
    remove_image,
    rename_sequence,
    rename_sequences,
    unlink_target,
)


def _setup(tmp_path: Path) -> tuple[Path, Path, Path]:
    fasta = tmp_path / "CroSat1.fasta"
    fasta.write_text(">CroSat1\nACGT\n", encoding="utf-8")
    other = tmp_path / "CroSat3.fasta"
    other.write_text(">CroSat3\nACGT\n", encoding="utf-8")
    image = tmp_path / "incoming" / "fish.png"
    image.parent.mkdir()
    image.write_bytes(b"\x89PNG fake")
    return fasta, other, image


def test_add_image_copies_into_hidden_store(tmp_path: Path) -> None:
    fasta, _other, source = _setup(tmp_path)
    img = add_image(tmp_path, source, targets=[ImageTarget(fasta.name)])
    stored = image_file_path(tmp_path, img)
    assert stored.parent.name == STORE_DIR_NAME
    assert stored.read_bytes() == source.read_bytes()
    assert img.original_name == "fish.png"
    assert image_count_for_sequence(fasta) == 1
    # The source stays where it was.
    assert source.is_file()


def test_empty_metadata_is_not_stored(tmp_path: Path) -> None:
    fasta, _other, source = _setup(tmp_path)
    add_image(tmp_path, source, targets=[ImageTarget(fasta.name)], metadata={"species": "", "notes": " "})
    assert load_images(tmp_path)[0].metadata == {}


def test_multi_probe_image_is_found_from_both_sequences(tmp_path: Path) -> None:
    fasta, other, source = _setup(tmp_path)
    add_image(
        tmp_path,
        source,
        targets=[
            ImageTarget(fasta.name, dye="Cy3", color="#ff3b30"),
            ImageTarget(other.name, dye="FITC", color="#34c759"),
        ],
    )
    assert len(images_for_sequence(fasta)) == 1
    assert len(images_for_sequence(other)) == 1
    assert images_for_sequence(fasta)[0].id == images_for_sequence(other)[0].id


def test_feature_lookup_exact_then_by_name(tmp_path: Path) -> None:
    fasta, _other, source = _setup(tmp_path)
    add_image(
        tmp_path,
        source,
        targets=[ImageTarget(fasta.name, feature=FeatureRef("probeA", 10, 50))],
    )
    assert len(images_for_feature(fasta, "probeA", 10, 50)) == 1
    # Annotation resized later: still found by its name.
    assert len(images_for_feature(fasta, "probeA", 12, 60)) == 1
    assert images_for_feature(fasta, "probeB", 10, 50) == []
    exact, names = feature_keys_with_images(fasta)
    assert ("probeA", 10, 50) in exact
    assert names == {"probeA"}


def test_rename_sequence_follows_file_rename(tmp_path: Path) -> None:
    fasta, _other, source = _setup(tmp_path)
    add_image(tmp_path, source, targets=[ImageTarget(fasta.name)])
    new = tmp_path / "CroSat1_v2.fasta"
    fasta.rename(new)
    assert rename_sequence(tmp_path, fasta.name, new.name) == 1
    assert image_count_for_sequence(new) == 1
    assert image_count_for_sequence(fasta) == 0


def test_rename_sequences_swaps_links_at_once(tmp_path: Path) -> None:
    fasta, other, source = _setup(tmp_path)
    add_image(tmp_path, source, targets=[ImageTarget(fasta.name)])
    second = tmp_path / "incoming" / "fish2.png"
    second.write_bytes(b"\x89PNG fake 2")
    add_image(tmp_path, second, targets=[ImageTarget(other.name)])
    # Batch Rename swapped the two names: each link must follow its own file.
    assert rename_sequences(tmp_path, {fasta.name: other.name, other.name: fasta.name}) == 2
    files = {img.targets[0].sequence: img.original_name for img in load_images(tmp_path)}
    assert files == {other.name: "fish.png", fasta.name: "fish2.png"}
    assert rename_sequences(tmp_path, {}) == 0


def test_remove_and_unlink(tmp_path: Path) -> None:
    fasta, other, source = _setup(tmp_path)
    img = add_image(tmp_path, source, targets=[ImageTarget(fasta.name), ImageTarget(other.name)])
    stored = image_file_path(tmp_path, img)
    unlink_target(tmp_path, img.id, fasta.name)
    assert image_count_for_sequence(fasta) == 0
    assert image_count_for_sequence(other) == 1
    assert stored.is_file()
    unlink_target(tmp_path, img.id, other.name)  # last link -> image goes
    assert load_images(tmp_path) == []
    assert not stored.exists()
    assert remove_image(tmp_path, "nope") is False


def test_caption_multicolour() -> None:
    image = CytoImage(
        id="x",
        file="x.png",
        targets=[
            ImageTarget("CroSat1.fasta", dye="Cy3", color="#ff3b30"),
            ImageTarget("CroSat3.fasta", dye="FITC", color="#34c759"),
        ],
        metadata={
            "technique": "FISH",
            "species": "Crocus sativus",
            "stage": "Mitotic metaphase",
            "counterstain": "DAPI",
        },
    )
    assert build_caption(image) == (
        "FISH of CroSat1 (Cy3, red) and CroSat3 (FITC, green) on Crocus sativus "
        "mitotic metaphase chromosomes, DAPI counterstain."
    )


def test_caption_without_metadata_uses_probe_name() -> None:
    image = CytoImage(id="x", file="x.png", targets=[ImageTarget("CroSat1.fasta")])
    assert build_caption(image) == "CroSat1."
