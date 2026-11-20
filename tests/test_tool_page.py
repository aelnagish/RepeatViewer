"""A tool window shown as a tab: the Primer Database."""

from PySide6.QtWidgets import QDialog, QPushButton, QTabWidget, QVBoxLayout


def test_tool_page_holds_the_window_and_hides_its_close(qapp):
    from PySide6.QtCore import Qt

    from repeatviewer.views.tool_page import ToolDialogPage

    dialog = QDialog()
    lay = QVBoxLayout(dialog)
    close = QPushButton("Close", dialog)
    cancel = QPushButton("Cancel", dialog)
    other = QPushButton("Add…", dialog)
    for b in (close, cancel, other):
        lay.addWidget(b)

    page = ToolDialogPage(dialog, "Primer Database")
    assert page.tab_name() == "Primer Database"
    assert page.inner is dialog and dialog.parent() is page
    assert not (dialog.windowFlags() & Qt.WindowType.Window)
    assert close.isHidden() and cancel.isHidden() and not other.isHidden()

    tabs = QTabWidget()
    tabs.addTab(page, page.tab_name())
    assert tabs.tabText(0) == "Primer Database"
    page.deleteLater()


def test_database_opens_as_a_tab(qapp, tmp_path):
    from repeatviewer.views.primer_database_dialog import PrimerDatabaseDialog
    from repeatviewer.views.tool_page import ToolDialogPage

    dialog = PrimerDatabaseDialog(tmp_path / "primers.sqlite")
    page = ToolDialogPage(dialog, "Primer Database")
    page.resize(900, 500)
    page.grab()                                   # paints inside a tab, not a window
    visible_close = [b for b in dialog.findChildren(QPushButton)
                     if b.text().strip() == "Close" and not b.isHidden()]
    assert visible_close == []
    page.deleteLater()
