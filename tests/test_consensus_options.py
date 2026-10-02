from repeatviewer.job_options import (
    ConsensusOptions,
    MAFFT_SEQUENCE_AUTO,
    MAFFT_SEQUENCE_NUCLEOTIDE,
    MAFFT_SEQUENCE_PROTEIN,
    load_consensus_options,
    normalize_consensus_sequence_type,
    save_consensus_options,
)


def test_normalize_consensus_sequence_type() -> None:
    assert normalize_consensus_sequence_type("auto") == MAFFT_SEQUENCE_AUTO
    assert normalize_consensus_sequence_type("nucleotide") == MAFFT_SEQUENCE_NUCLEOTIDE
    assert normalize_consensus_sequence_type("dna") == MAFFT_SEQUENCE_NUCLEOTIDE
    assert normalize_consensus_sequence_type("protein") == MAFFT_SEQUENCE_PROTEIN
    assert normalize_consensus_sequence_type("aa") == MAFFT_SEQUENCE_PROTEIN
    assert normalize_consensus_sequence_type("unknown") == MAFFT_SEQUENCE_AUTO


def test_save_and_load_consensus_sequence_type(monkeypatch) -> None:
    saved: dict = {}

    def _save(group, values):
        saved[group] = values

    def _load(group):
        return saved.get(group, {})

    monkeypatch.setattr("repeatviewer.job_options._save_settings_group", _save)
    monkeypatch.setattr("repeatviewer.job_options._load_settings_group", _load)

    save_consensus_options(
        ConsensusOptions(
            threshold=0.75,
            ambiguous="X",
            sequence_type=MAFFT_SEQUENCE_PROTEIN,
        )
    )
    loaded = load_consensus_options()
    assert loaded.threshold == 0.75
    assert loaded.ambiguous == "X"
    assert loaded.sequence_type == MAFFT_SEQUENCE_PROTEIN
