from repeatviewer.backends.blast_genome_plot import (
    PLOTLY_HTML_CONFIG,
    _plotly_interactive_controls_html,
)


def test_density_plot_wheel_does_not_scroll_zoom() -> None:
    """Plain wheel scrolls the page; Shift+scroll zooms X via custom handler."""
    assert PLOTLY_HTML_CONFIG.get("scrollZoom") is False
    html = _plotly_interactive_controls_html(
        scale_label="Genome track bar scale",
        panel_mode="genome",
    )
    assert "applyShiftWheelXZoom" in html
    assert "teWheelWantsXZoom" in html
    assert "teShiftDown" in html
    assert "gd.addEventListener('wheel'" in html
    # macOS Shift+scroll often arrives as deltaX — handler must use either axis.
    assert "Number(e.deltaX)" in html
    assert "capture: true" in html
    # Zoom keeps the genome position under the cursor fixed.
    assert "anchorX - frac * newSpan" in html


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
