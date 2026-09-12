# Runtime hook for the Windows frozen app.
#
# Python 3.12's module repr treats a loader with origin=None as a namespace
# loader and reads loader._path. six._SixMetaPathImporter has no _path.
# PyInstaller's pyi_rth_pyside6 loads Shiboken first; Shiboken then inspects
# python-dateutil (via pandas) and inspect.getfile() repr()s the six.moves
# module, which crashes:
#
#   AttributeError: '_SixMetaPathImporter' object has no attribute '_path'
#
# Custom runtime hooks run before pyi_rth_pyside6.

from __future__ import annotations


def _patch_six_importer() -> None:
    try:
        import six
    except Exception:
        return

    cls = getattr(six, "_SixMetaPathImporter", None)
    if cls is not None and not hasattr(cls, "_path"):
        cls._path = []

    importer = getattr(six, "_importer", None)
    if importer is not None and not hasattr(importer, "_path"):
        try:
            importer._path = []
        except Exception:
            pass

    import sys

    for finder in sys.meta_path:
        if type(finder).__name__ == "_SixMetaPathImporter" and not hasattr(finder, "_path"):
            try:
                finder._path = []
            except Exception:
                pass


_patch_six_importer()
