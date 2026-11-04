"""Right rail for the tree preview: vertical Display tab matching MSA Display."""

from __future__ import annotations

from PySide6.QtCore import QSettings, Qt, Signal
from PySide6.QtGui import QColor, QPalette
from PySide6.QtWidgets import (
    QButtonGroup,
    QFrame,
    QHBoxLayout,
    QScrollArea,
    QSizePolicy,
    QStackedWidget,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.widgets.sequence_annotation_rail import (
    _RAIL_TAB_WIDTH,
    _RotatedTabButton,
    _h_divider,
    _page_title_with_rule,
)

_MODULE_WIDTH = 260
_COLLAPSED_RAIL_WIDTH = _RAIL_TAB_WIDTH
_SETTINGS_GROUP = "tree_display_rail"
_SETTINGS_COLLAPSED_KEY = "panel_collapsed"
_SETTINGS_TAB_KEY = "active_tab"
_LEGACY_OPTIONS_VISIBLE = "tree_view/options_panel"
# Same divider colour as SequenceMapPanel ↔ Annotation rail.
_RAIL_DIVIDER = QColor("#9a9a9a")


def _fill_window_surface(widget: QWidget) -> None:
    """Opaque Window fill — same as MSA SequenceAnnotationRail pages."""
    widget.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)
    widget.setAutoFillBackground(True)
    color = widget.palette().color(QPalette.ColorRole.Window)
    pal = widget.palette()
    pal.setColor(QPalette.ColorRole.Window, color)
    pal.setColor(QPalette.ColorRole.Base, color)
    widget.setPalette(pal)


class TreeDisplayRail(QWidget):
    """Vertical Display pill + options panel (MSA Annotation-rail chrome)."""

    panelCollapsedChanged = Signal(bool)
    tabChanged = Signal(int)

    TAB_DISPLAY = 0
    TAB_FIND = 1

    def __init__(
        self,
        content: QWidget,
        parent: QWidget | None = None,
        *,
        footer: QWidget | None = None,
        find_page: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.setObjectName("treeDisplayRail")
        self.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Expanding)
        self._module_width = _MODULE_WIDTH
        self._panel_collapsed = False
        self._active_tab = self.TAB_DISPLAY
        self._updating_tabs = False
        self.setFixedWidth(self._expanded_rail_width())
        _fill_window_surface(self)

        page = QWidget(self)
        page.setObjectName("TreeDisplayPage")
        _fill_window_surface(page)
        page_layout = QVBoxLayout(page)
        page_layout.setContentsMargins(0, 0, 0, 0)
        page_layout.setSpacing(0)

        title_host = QWidget(page)
        title_layout = QVBoxLayout(title_host)
        title_layout.setContentsMargins(10, 10, 10, 0)
        title_layout.setSpacing(0)
        display_title, _rule = _page_title_with_rule("Display", title_host)
        title_layout.addWidget(display_title)
        page_layout.addWidget(title_host)

        scroll = QScrollArea(page)
        scroll.setObjectName("TreeOptionsPanel")
        scroll.setWidget(content)
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        _fill_window_surface(scroll)
        _fill_window_surface(scroll.viewport())
        _fill_window_surface(content)
        scroll.setStyleSheet(
            "QScrollArea#TreeOptionsPanel { background: palette(window); border: none; }"
        )
        page_layout.addWidget(scroll, 1)
        self._options_panel = scroll

        if footer is not None:
            footer.setParent(page)
            _fill_window_surface(footer)
            page_layout.addWidget(footer, 0)

        self._stack = QStackedWidget(self)
        self._stack.setObjectName("TreeDisplayStack")
        self._stack.setFixedWidth(self._module_width)
        self._stack.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Expanding)
        _fill_window_surface(self._stack)
        page.setFixedWidth(self._module_width)
        self._stack.addWidget(page)
        if find_page is not None:
            find_page.setParent(self._stack)
            _fill_window_surface(find_page)
            find_page.setFixedWidth(self._module_width)
            self._stack.addWidget(find_page)

        # Divider between module column and vertical Display tab.
        self._tab_divider = QFrame(self)
        self._tab_divider.setFrameShape(QFrame.Shape.VLine)
        self._tab_divider.setFixedWidth(1)
        self._tab_divider.setStyleSheet("color: palette(mid);")

        self._tab_strip = QWidget(self)
        self._tab_strip.setFixedWidth(_RAIL_TAB_WIDTH)
        self._tab_strip.setObjectName("treeDisplayRailTabs")
        _fill_window_surface(self._tab_strip)
        tab_layout = QVBoxLayout(self._tab_strip)
        tab_layout.setContentsMargins(1, 2, 1, 2)
        tab_layout.setSpacing(2)

        self._tab_group = QButtonGroup(self)
        self._tab_group.setExclusive(False)
        self._tab_buttons: list[_RotatedTabButton] = []
        tabs = [("Display", "Display — layout, labels, and colours (click again to hide)")]
        if find_page is not None:
            tabs.append(("Find", "Find — tips by name or trait value (click again to hide)"))
        for index, (label, tip) in enumerate(tabs):
            btn = _RotatedTabButton(label, self._tab_strip)
            btn.setToolTip(tip)
            self._tab_group.addButton(btn, index)
            tab_layout.addWidget(btn, 0, Qt.AlignmentFlag.AlignHCenter)
            self._tab_buttons.append(btn)
            btn.clicked.connect(lambda _c=False, i=index: self._on_tab_clicked(i))
        self._display_btn = self._tab_buttons[0]
        tab_layout.addStretch(1)

        # Outer left edge of the rail (canvas | rail), same as MSA map divider.
        self._canvas_divider = QFrame(self)
        self._canvas_divider.setObjectName("treeMapRailDivider")
        self._canvas_divider.setFrameShape(QFrame.Shape.VLine)
        self._canvas_divider.setFixedWidth(1)
        self._canvas_divider.setAutoFillBackground(True)
        dp = self._canvas_divider.palette()
        dp.setColor(QPalette.ColorRole.Window, _RAIL_DIVIDER)
        self._canvas_divider.setPalette(dp)
        self._canvas_divider.setStyleSheet(
            f"QFrame#treeMapRailDivider {{ background-color: {_RAIL_DIVIDER.name()}; "
            f"border: none; max-width: 1px; }}"
        )

        root = QHBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)
        root.addWidget(self._canvas_divider, 0)
        root.addWidget(self._stack, 0)
        root.addWidget(self._tab_divider, 0)
        root.addWidget(self._tab_strip, 0)

        self._restore_collapsed()

    def options_scroll(self) -> QScrollArea:
        return self._options_panel

    def _expanded_rail_width(self) -> int:
        # module + tab strip + two 1px dividers (canvas|module and module|tabs)
        return int(self._module_width) + _RAIL_TAB_WIDTH + 2

    def is_panel_collapsed(self) -> bool:
        return bool(self._panel_collapsed)

    def set_panel_collapsed(self, collapsed: bool, *, persist: bool = True) -> None:
        collapsed = bool(collapsed)
        if collapsed == self._panel_collapsed:
            if collapsed:
                self._sync_tab_checked(False)
            return
        self._panel_collapsed = collapsed
        self._stack.setVisible(not collapsed)
        self._tab_divider.setVisible(not collapsed)
        # Keep the canvas|rail edge when collapsed so the tab strip stays framed.
        if collapsed:
            self.setFixedWidth(_COLLAPSED_RAIL_WIDTH + 1)  # tab + left divider
            self._sync_tab_checked(False)
        else:
            self.setFixedWidth(self._expanded_rail_width())
            self._sync_tab_checked(True)
        self.panelCollapsedChanged.emit(collapsed)
        if persist:
            self._save_collapsed()

    def current_tab(self) -> int:
        return self._active_tab

    def show_tab(self, index: int, *, persist: bool = True) -> None:
        """Show the Display (0) or Find (1) page and open the rail if it was hidden."""
        index = max(0, min(int(index), self._stack.count() - 1))
        self._active_tab = index
        self._stack.setCurrentIndex(index)
        if self._panel_collapsed:
            self.set_panel_collapsed(False, persist=persist)
        else:
            self._sync_tab_checked(True)
            if persist:
                self._save_collapsed()
        self.tabChanged.emit(index)

    def _sync_tab_checked(self, checked: bool) -> None:
        self._updating_tabs = True
        try:
            for i, btn in enumerate(self._tab_buttons):
                btn.setChecked(bool(checked) and i == self._active_tab)
                btn.update()
        finally:
            self._updating_tabs = False

    def _on_tab_clicked(self, index: int = 0) -> None:
        if self._updating_tabs:
            return
        if not self._panel_collapsed and index == self._active_tab:
            self.set_panel_collapsed(True)
            return
        self.show_tab(index)

    def _settings(self) -> QSettings:
        s = QSettings()
        s.beginGroup(_SETTINGS_GROUP)
        return s

    def _restore_collapsed(self) -> None:
        settings = self._settings()
        try:
            if settings.contains(_SETTINGS_COLLAPSED_KEY):
                raw = settings.value(_SETTINGS_COLLAPSED_KEY, False)
            else:
                legacy = QSettings().value(_LEGACY_OPTIONS_VISIBLE, True)
                visible = str(legacy).lower() in ("1", "true", "yes")
                raw = not visible
            collapsed = str(raw).lower() in ("1", "true", "yes")
            try:
                tab = int(settings.value(_SETTINGS_TAB_KEY, 0))
            except (TypeError, ValueError):
                tab = 0
        finally:
            settings.endGroup()
        self._active_tab = tab if 0 <= tab < self._stack.count() else 0
        self._stack.setCurrentIndex(self._active_tab)
        self.set_panel_collapsed(collapsed, persist=False)
        self._sync_tab_checked(not collapsed)

    def _save_collapsed(self) -> None:
        settings = self._settings()
        try:
            settings.setValue(_SETTINGS_COLLAPSED_KEY, bool(self._panel_collapsed))
            settings.setValue(_SETTINGS_TAB_KEY, int(self._active_tab))
        finally:
            settings.endGroup()
        QSettings().setValue("tree_view/options_panel", not self._panel_collapsed)
