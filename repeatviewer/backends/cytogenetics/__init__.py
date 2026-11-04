"""Cytogenetics workspace backend (Qt-free).

Builds on the image store in ``backends/cyto_images.py``: every image is still a
``CytoImage`` record in ``.repeatviewer_cytogenetics/index.json``. Version 2 adds,
per image, the channel TIFFs, display settings, marked chromosomes and signals
(``model.ImageDetails``, kept in ``CytoImage.extra``). Slides, karyotypes and the
editable lists live in their own JSON files next to ``index.json``.
"""
