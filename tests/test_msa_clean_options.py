from teviewer.job_options import MsaCleanOptions


def test_msa_clean_defaults_to_gappy_columns() -> None:
    options = MsaCleanOptions()

    assert options.use_divergence is False
    assert options.use_gap_row is False
    assert options.use_gappy_columns is True


def test_msa_clean_accepts_combined_methods() -> None:
    options = MsaCleanOptions(
        use_divergence=True,
        use_gap_row=True,
        use_gappy_columns=True,
    )

    assert options.use_divergence is True
    assert options.use_gap_row is True
    assert options.use_gappy_columns is True


def test_msa_clean_keeps_at_least_one_method_enabled() -> None:
    options = MsaCleanOptions(
        use_divergence=False,
        use_gap_row=False,
        use_gappy_columns=False,
    )

    assert options.use_gappy_columns is True
