"""Public API for TEviewer plugins."""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable, Iterable

from teviewer.project import TEviewerProject

if TYPE_CHECKING:
    from PySide6.QtWidgets import QWidget

ProgressCallback = Callable[[int, str], None]


@dataclass
class PluginBatchResult:
    """Same shape as ``run_batch`` so the main window can refresh the file list."""

    func_name: str
    total: int
    failures: list[str] = field(default_factory=list)
    created_outputs: list[str] = field(default_factory=list)
    source_files: list[str] = field(default_factory=list)

    def as_dict(self) -> dict:
        return {
            "func_name": self.func_name,
            "total": self.total,
            "failures": list(self.failures),
            "created_outputs": list(self.created_outputs),
            "source_files": list(self.source_files),
        }


class TEviewerPlugin(ABC):
    """
    Base class for third-party tools.

    Plugins receive the user's current file selection from the Files panel and
    run in a background thread (same as Blast / TEAid). Write outputs next to
    each input file unless you have a reason to use another folder.
    """

    #: Unique id (letters, numbers, underscores). Shown as ``plugin:<id>`` internally.
    id: str
    #: Toolbar / menu label.
    label: str
    tooltip: str = ""
    #: If True (default), only ``.fa`` / ``.fasta`` / ``.fas`` are passed in.
    fasta_only: bool = True
    #: When set, only these suffixes are accepted (lowercase, e.g. ``{".out"}``).
    file_suffixes: set[str] | None = None
    #: If True, the tool may run on multiple files in parallel (uses batch_workers).
    supports_parallel: bool = True
    #: If True, ``prepare_job`` shows a dialog where the user picks input files.
    picks_files_in_dialog: bool = False
    #: If False, the plugin can run without a pre-selected Files panel entry.
    needs_selected_files: bool = True

    @abstractmethod
    def validate(
        self,
        files: list[Path],
        project: TEviewerProject,
    ) -> str | None:
        """
        Return an error message to block the run, or None if OK.

        Use this for missing genome, missing external binaries, etc.
        """

    @abstractmethod
    def prepare_job(
        self,
        parent: QWidget,
        files: list[Path],
        project: TEviewerProject,
    ) -> bool:
        """
        Optional UI step before the background job (e.g. pick reference files).

        Return False to cancel. Store per-run options on ``self`` (e.g. ``self._job_options``).
        """
        self._job_options: dict[str, Any] = {}
        return True

    def run(
        self,
        file_path: Path,
        project: TEviewerProject,
        **job_options: Any,
    ) -> list[str]:
        """
        Run on one selected file.

        Returns paths (as strings) of **new files** created (for highlighting).
        Raise an exception if the run failed for this file.
        """

    def run_batch(
        self,
        files: Iterable[Path],
        project: TEviewerProject,
        progress_cb: ProgressCallback,
        *,
        max_workers: int = 1,
    ) -> PluginBatchResult:
        """
        Default batch runner (sequential or threaded). Override for custom batching.
        """
        from concurrent.futures import ThreadPoolExecutor, as_completed
        import threading

        file_list = [Path(p).expanduser().resolve() for p in files]
        total = len(file_list)
        failures: list[str] = []
        created: list[str] = []
        sources: list[str] = []
        lock = threading.Lock()
        completed = 0

        def _run_one(path: Path) -> None:
            nonlocal completed
            try:
                opts = getattr(self, "_job_options", {})
                outs = self.run(path, project, **opts)
            except Exception as exc:  # noqa: BLE001
                with lock:
                    failures.append(f"{path.name}: {exc}")
            else:
                with lock:
                    created.extend(str(Path(o).resolve()) for o in outs)
                    sources.append(str(path))
            with lock:
                completed += 1
                progress_cb(completed - 1, f"{self.label}: {path.name}")

        workers = max(1, min(max_workers, total)) if self.supports_parallel else 1
        if workers == 1:
            for path in file_list:
                _run_one(path)
        else:
            with ThreadPoolExecutor(max_workers=workers) as pool:
                futures = [pool.submit(_run_one, p) for p in file_list]
                for _ in as_completed(futures):
                    pass

        progress_cb(total, f"{self.label}: finished")
        return PluginBatchResult(
            func_name=f"plugin:{self.id}",
            total=total,
            failures=failures,
            created_outputs=created,
            source_files=sources,
        )
