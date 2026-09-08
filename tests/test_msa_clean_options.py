from repeatviewer.job_options import DEFAULT_MSA_CLEAN_ORDER, MsaCleanOptions


def test_msa_clean_defaults_to_gappy_columns_then_divergence() -> None:
    options = MsaCleanOptions()

    assert options.use_divergence is True
    assert options.use_gap_row is False
    assert options.use_gappy_columns is True
    assert options.method_order == DEFAULT_MSA_CLEAN_ORDER


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


def test_msa_clean_preserves_custom_method_order() -> None:
    order = ("msa_clean_gappy", "msa_clean_gap", "msa_clean_div")

    options = MsaCleanOptions(method_order=order)

    assert options.method_order == order
