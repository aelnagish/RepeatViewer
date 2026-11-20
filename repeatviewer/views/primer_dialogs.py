"""What is left of the old primer windows.

The Design and Map windows are gone: both live in the Primers panel now
(Design page) and as a tab (Map), so the settings of what is open sit in the
right panel and nothing floats over the sequence. This module keeps the two
things the rest of the program still asks for by this name:

* :func:`primer3_missing_message` – the install hint when primer3-py is absent
* :class:`PrimerDatabaseDialog`   – re-exported from its own module

Developer: Abdullah Yusuf
"""

from __future__ import annotations

from PySide6.QtWidgets import QMessageBox, QWidget

from repeatviewer.backends.primers.thermo import PRIMER3_MISSING, primer3_available
from repeatviewer.views.primer_database_dialog import PrimerDatabaseDialog


def primer3_missing_message(parent: QWidget | None) -> bool:
    """Show the install hint when primer3-py is missing; True when it is."""
    if primer3_available():
        return False
    QMessageBox.warning(parent, "Primers", PRIMER3_MISSING)
    return True


__all__ = [
    "PrimerDatabaseDialog",
    "primer3_missing_message",
]
