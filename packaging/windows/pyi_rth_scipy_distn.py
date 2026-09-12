# Runtime hook for the Windows frozen app (Python 3.12.0 + PyInstaller).
#
# scipy.stats._distn_infrastructure ends with:
#
#   for obj in [s for s in dir() if s.startswith('_doc_')]:
#       exec('del ' + obj)
#   del obj
#
# In the frozen importer the list is empty, so `obj` is never bound and
# `del obj` raises NameError. That is a CPython 3.12.0 / PyInstaller issue
# (fixed in 3.12.1). Patch the source at import time instead of requiring
# a Python upgrade.

from __future__ import annotations

import sys


def _patch_source(source: str) -> str:
    target = "\ndel obj\n"
    replacement = "\nglobals().pop('obj', None)\n"
    if target in source:
        return source.replace(target, replacement, 1)
    return source


class _ScipyDistnLoader:
    def __init__(self, inner) -> None:
        self._inner = inner

    def __getattr__(self, name):
        return getattr(self._inner, name)

    def create_module(self, spec):
        return self._inner.create_module(spec)

    def exec_module(self, module):
        source = None
        try:
            source = self._inner.get_source(module.__name__)
        except Exception:
            source = None
        if source:
            patched = _patch_source(source)
            filename = getattr(self._inner, "path", module.__name__)
            exec(compile(patched, filename, "exec"), module.__dict__)
            return
        module.__dict__.setdefault("obj", None)
        self._inner.exec_module(module)


class _ScipyDistnFinder:
    _TARGET = "scipy.stats._distn_infrastructure"

    def find_spec(self, fullname, path=None, target=None):
        if fullname != self._TARGET:
            return None
        for finder in sys.meta_path:
            if finder is self:
                continue
            find = getattr(finder, "find_spec", None)
            if find is None:
                continue
            try:
                spec = find(fullname, path, target)
            except Exception:
                continue
            if spec is None or spec.loader is None:
                continue
            spec.loader = _ScipyDistnLoader(spec.loader)
            return spec
        return None


sys.meta_path.insert(0, _ScipyDistnFinder())
