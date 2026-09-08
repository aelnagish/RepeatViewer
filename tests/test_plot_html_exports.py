from repeatviewer.backends.blast_genome_plot import _plotly_interactive_controls_html


def test_plot_controls_offer_svg_and_pdf_exports() -> None:
    html = _plotly_interactive_controls_html(
        scale_label="Y axis zoom",
        panel_mode="divergence",
    )

    assert 'id="te-sf-save-svg"' in html
    assert 'id="te-sf-save-pdf"' in html
    assert "Plotly.downloadImage(gd, {format: 'svg'" in html
    assert "window.print();" in html
    assert "@media print" in html


def test_te_density_controls_offer_svg_and_pdf_exports() -> None:
    """TE density plots use the genome-map controls rather than divergence controls."""
    html = _plotly_interactive_controls_html(
        scale_label="Genome track bar scale",
        panel_mode="genome",
        prefs_path="/tmp/density.prefs.json",
    )

    assert 'id="te-sf-save-svg"' in html
    assert 'id="te-sf-save-pdf"' in html
    assert "Plotly.downloadImage(gd, {format: 'svg'" in html
    assert "window.print();" in html


def test_te_density_controls_accept_plotly3_binary_arrays() -> None:
    """Plotly.js v3 stores large y arrays as {{dtype,bdata}}; panel must still init."""
    html = _plotly_interactive_controls_html(
        scale_label="Genome track bar scale",
        panel_mode="genome",
    )

    assert "function hasTraceValues" in html
    assert "function traceNumericArray" in html
    assert "!hasTraceValues(tr.y)" in html
    assert "Array.isArray(tr.y) || tr.y.length === 0" not in html
