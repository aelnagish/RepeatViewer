"""Save a drawing as PNG (at twice the screen resolution) or SVG (vector)."""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

from PySide6.QtCore import QRect, QSize, QSettings
from PySide6.QtGui import QColor, QImage, QPainter

IMAGE_FILTER = "PNG image (*.png);;SVG vector image (*.svg)"
PNG_SCALE = 2          # pixels per screen point in a saved PNG


def write_image(path: Path, w: int, h: int, render: Callable[[QPainter, int, int], None],
                background: QColor | None = None, title: str = "") -> Path:
    """Draw ``render(painter, w, h)`` into ``path`` (.svg: vector, else PNG)."""
    path = Path(path)
    if path.suffix.lower() == ".svg":
        from PySide6.QtSvg import QSvgGenerator

        gen = QSvgGenerator()
        gen.setFileName(str(path))
        gen.setSize(QSize(w, h))
        gen.setViewBox(QRect(0, 0, w, h))
        gen.setTitle(title or path.stem)
        gen.setDescription("RepeatViewer genome browser")
        p = QPainter(gen)
        try:
            if background is not None:
                p.fillRect(0, 0, w, h, background)
            render(p, w, h)
        finally:
            p.end()
        return path
    if path.suffix.lower() != ".png":
        path = path.with_suffix(".png")
    img = QImage(w * PNG_SCALE, h * PNG_SCALE, QImage.Format.Format_ARGB32_Premultiplied)
    img.setDevicePixelRatio(PNG_SCALE)
    img.fill(background or QColor("white"))
    p = QPainter(img)
    try:
        render(p, w, h)
    finally:
        p.end()
    if not img.save(str(path)):
        raise OSError(f"Cannot write {path}")
    return path


def ask_image_path(parent, default_name: str) -> Path | None:
    from PySide6.QtWidgets import QFileDialog

    folder = QSettings().value("genome_browser/image_dir", "") or str(Path.home())
    path, chosen = QFileDialog.getSaveFileName(parent, "Save image",
                                               str(Path(folder) / f"{default_name}.png"),
                                               IMAGE_FILTER)
    if not path:
        return None
    out = Path(path)
    if out.suffix.lower() not in (".png", ".svg"):
        out = out.with_suffix(".svg" if "SVG" in (chosen or "") else ".png")
    QSettings().setValue("genome_browser/image_dir", str(out.parent))
    return out


def save_widget_image(parent, widget, default_name: str,
                      render: Callable[[QPainter, int, int], None]) -> Path | None:
    """Ask where, then save ``widget``'s drawing there (PNG or SVG)."""
    from PySide6.QtWidgets import QMessageBox

    out = ask_image_path(parent, default_name)
    if out is None:
        return None
    pal = widget.palette()
    try:
        return write_image(out, widget.width(), widget.height(), render,
                           pal.color(pal.ColorRole.Base), default_name)
    except OSError as exc:
        QMessageBox.warning(parent, "Save image", str(exc))
        return None
