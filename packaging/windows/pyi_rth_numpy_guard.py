# Runtime hook: import numpy once, early, under a single module path.
# Mitigates Windows frozen-app failures where numpy's C extension is
# initialized twice ("cannot load module more than once per process").

from __future__ import annotations


def _guard() -> None:
    try:
        import numpy  # noqa: F401
    except Exception:
        # Let the real import site report the full traceback.
        return


_guard()
