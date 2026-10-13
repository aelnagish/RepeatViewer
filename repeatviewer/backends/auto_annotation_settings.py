"""User-configurable parameters for sequence-map auto annotation."""

from __future__ import annotations

from dataclasses import dataclass, field

from PySide6.QtCore import QSettings

DEFAULT_LTR_SIMILARITY_THRESHOLD = 0.6
DEFAULT_TIR_SIMILARITY_THRESHOLD = 0.5


DEFAULT_MIN_DISTANCE_TO_END = 1000


DEFAULT_ANNOTATION_THREADS = 4
# Shared by Feature detection + Similarity annotation (and Sequence Annotation Settings).
ANNOTATION_RUN_PREFIX = "annotation_run/"


def default_annotation_threads() -> int:
    """Default Threads for Feature detection / Similarity (same as other tools)."""
    return DEFAULT_ANNOTATION_THREADS


@dataclass
class TerminalRepeatSettings:
    min_length: int = 100
    max_length: int = 2000
    min_pair_distance: int = 100
    min_distance_to_end: int = DEFAULT_MIN_DISTANCE_TO_END
    similarity_threshold: float = DEFAULT_LTR_SIMILARITY_THRESHOLD


@dataclass
class TandemRepeatSettings:
    min_repeat_2: int = 5
    min_repeat_3: int = 4
    min_repeat_4: int = 3
    min_repeat_5: int = 3
    min_repeat_6: int = 3
    enable_approximate: bool = True
    atr_min_motif: int = 5
    atr_max_motif: int = 20
    atr_min_seed_repeat: int = 3
    atr_min_seed_length: int = 10
    atr_min_identity: float = 0.75
    atr_max_errors: int = 2
    atr_max_extend: int = 500


# Reading frames used by orfipy ORF detection (+1..+3 / -1..-3).
DEFAULT_ORF_FRAMES: tuple[int, ...] = (1, 2, 3, -1, -2, -3)


def normalize_orf_frames(frames: object) -> tuple[int, ...]:
    """Keep only valid frames, preserve +1..+3 then -1..-3 order."""
    allowed = set(DEFAULT_ORF_FRAMES)
    chosen: set[int] = set()
    if isinstance(frames, str):
        tokens = frames.replace(",", " ").split()
        for token in tokens:
            try:
                value = int(token)
            except ValueError:
                continue
            if value in allowed:
                chosen.add(value)
    elif isinstance(frames, (set, frozenset, list, tuple)):
        for item in frames:
            try:
                value = int(item)
            except (TypeError, ValueError):
                continue
            if value in allowed:
                chosen.add(value)
    if not chosen:
        return DEFAULT_ORF_FRAMES
    return tuple(frame for frame in DEFAULT_ORF_FRAMES if frame in chosen)


@dataclass
class OrfSettings:
    min_length: int = 30
    max_length: int = 10_000
    alternative_start_codons: str = ""
    # Which reading frames to keep (orfipy ORF_frame values).
    frames: tuple[int, ...] = DEFAULT_ORF_FRAMES

    def __post_init__(self) -> None:
        self.frames = normalize_orf_frames(self.frames)


@dataclass
class AutoAnnotationSettings:
    ltr: TerminalRepeatSettings | None = None
    tir: TerminalRepeatSettings | None = None
    tandem: TandemRepeatSettings | None = None
    orf: OrfSettings | None = None
    # CPU threads inside each sequence job (like BLAST Threads).
    num_threads: int = field(default_factory=default_annotation_threads)
    show_process_window: bool = True
    close_process_window_when_finished: bool = True

    def __post_init__(self) -> None:
        self.num_threads = max(1, int(self.num_threads or default_annotation_threads()))
        if self.ltr is None:
            self.ltr = TerminalRepeatSettings(
                min_length=100,
                max_length=2000,
                similarity_threshold=DEFAULT_LTR_SIMILARITY_THRESHOLD,
            )
        if self.tir is None:
            self.tir = TerminalRepeatSettings(
                min_length=20,
                max_length=500,
                min_pair_distance=50,
                similarity_threshold=DEFAULT_TIR_SIMILARITY_THRESHOLD,
            )
        if self.tandem is None:
            self.tandem = TandemRepeatSettings()
        if self.orf is None:
            self.orf = OrfSettings()


def default_auto_annotation_settings() -> AutoAnnotationSettings:
    return AutoAnnotationSettings()


def _get_int(s: QSettings, key: str, default: int) -> int:
    return int(s.value(key, default))


def _get_float(s: QSettings, key: str, default: float) -> float:
    return float(s.value(key, default))


def _get_bool(s: QSettings, key: str, default: bool) -> bool:
    value = s.value(key, default)
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.lower() in ("1", "true", "yes", "on")
    return bool(value)


def load_annotation_run_options() -> tuple[int, bool, bool]:
    """Return ``(num_threads, show_process_window, close_when_finished)``."""
    s = QSettings()
    threads = default_annotation_threads()
    show = True
    close = True
    if s.value(f"{ANNOTATION_RUN_PREFIX}num_threads") is not None:
        threads = _get_int(
            s, f"{ANNOTATION_RUN_PREFIX}num_threads", default_annotation_threads()
        )
    else:
        for prefix in ("auto_annotation/", "similarity_annotation/"):
            if s.value(f"{prefix}num_threads") is not None:
                threads = _get_int(s, f"{prefix}num_threads", threads)
                break
    if s.value(f"{ANNOTATION_RUN_PREFIX}show_process_window") is not None:
        show = _get_bool(s, f"{ANNOTATION_RUN_PREFIX}show_process_window", True)
    else:
        for prefix in ("auto_annotation/", "similarity_annotation/", "te_annotation/"):
            if s.value(f"{prefix}show_process_window") is not None:
                show = _get_bool(s, f"{prefix}show_process_window", True)
                break
    close_key = f"{ANNOTATION_RUN_PREFIX}close_process_window_when_finished"
    if s.value(close_key) is not None:
        close = _get_bool(s, close_key, True)
    else:
        for prefix in ("auto_annotation/", "similarity_annotation/", "te_annotation/"):
            key = f"{prefix}close_process_window_when_finished"
            if s.value(key) is not None:
                close = _get_bool(s, key, True)
                break
    return max(1, int(threads)), bool(show), bool(close)


def save_annotation_run_options(
    *,
    num_threads: int,
    show_process_window: bool,
    close_process_window_when_finished: bool,
) -> None:
    """Persist shared Threads / process-window flags for both annotation panels."""
    s = QSettings()
    threads = max(1, int(num_threads))
    show = bool(show_process_window)
    close = bool(close_process_window_when_finished)
    for prefix in (
        ANNOTATION_RUN_PREFIX,
        "auto_annotation/",
        "similarity_annotation/",
        "te_annotation/",
    ):
        s.setValue(f"{prefix}num_threads", threads)
        s.setValue(f"{prefix}show_process_window", show)
        s.setValue(f"{prefix}close_process_window_when_finished", close)


def load_auto_annotation_settings() -> AutoAnnotationSettings:
    s = QSettings()
    prefix = "auto_annotation/"
    ltr = TerminalRepeatSettings(
        min_length=_get_int(s, f"{prefix}ltr_min_length", 100),
        max_length=_get_int(s, f"{prefix}ltr_max_length", 2000),
        min_pair_distance=_get_int(s, f"{prefix}ltr_min_pair_distance", 100),
        min_distance_to_end=_get_int(
            s, f"{prefix}ltr_min_distance_to_end", DEFAULT_MIN_DISTANCE_TO_END
        ),
        similarity_threshold=_get_float(
            s, f"{prefix}ltr_similarity", DEFAULT_LTR_SIMILARITY_THRESHOLD
        ),
    )
    tir = TerminalRepeatSettings(
        min_length=_get_int(s, f"{prefix}tir_min_length", 20),
        max_length=_get_int(s, f"{prefix}tir_max_length", 500),
        min_pair_distance=_get_int(s, f"{prefix}tir_min_pair_distance", 50),
        min_distance_to_end=_get_int(
            s, f"{prefix}tir_min_distance_to_end", DEFAULT_MIN_DISTANCE_TO_END
        ),
        similarity_threshold=_get_float(
            s, f"{prefix}tir_similarity", DEFAULT_TIR_SIMILARITY_THRESHOLD
        ),
    )
    tandem = TandemRepeatSettings(
        min_repeat_2=_get_int(s, f"{prefix}tandem_min_repeat_2", 5),
        min_repeat_3=_get_int(s, f"{prefix}tandem_min_repeat_3", 4),
        min_repeat_4=_get_int(s, f"{prefix}tandem_min_repeat_4", 3),
        min_repeat_5=_get_int(s, f"{prefix}tandem_min_repeat_5", 3),
        min_repeat_6=_get_int(s, f"{prefix}tandem_min_repeat_6", 3),
        enable_approximate=_get_bool(s, f"{prefix}tandem_enable_approximate", True),
        atr_min_motif=_get_int(s, f"{prefix}tandem_atr_min_motif", 5),
        atr_max_motif=_get_int(s, f"{prefix}tandem_atr_max_motif", 20),
        atr_min_seed_repeat=_get_int(s, f"{prefix}tandem_atr_min_seed_repeat", 3),
        atr_min_seed_length=_get_int(s, f"{prefix}tandem_atr_min_seed_length", 10),
        atr_min_identity=_get_float(s, f"{prefix}tandem_atr_min_identity", 0.75),
        atr_max_errors=_get_int(s, f"{prefix}tandem_atr_max_errors", 2),
        atr_max_extend=_get_int(s, f"{prefix}tandem_atr_max_extend", 500),
    )
    orf = OrfSettings(
        min_length=_get_int(s, f"{prefix}orf_min_length", 30),
        max_length=_get_int(s, f"{prefix}orf_max_length", 10_000),
        alternative_start_codons=str(s.value(f"{prefix}orf_alt_starts", "") or ""),
        frames=normalize_orf_frames(s.value(f"{prefix}orf_frames", "")),
    )
    threads, show, close = load_annotation_run_options()
    return AutoAnnotationSettings(
        ltr=ltr,
        tir=tir,
        tandem=tandem,
        orf=orf,
        num_threads=threads,
        show_process_window=show,
        close_process_window_when_finished=close,
    )


def save_auto_annotation_settings(settings: AutoAnnotationSettings) -> None:
    s = QSettings()
    prefix = "auto_annotation/"
    ltr = settings.ltr or TerminalRepeatSettings()
    tir = settings.tir or TerminalRepeatSettings(
        min_length=20,
        similarity_threshold=DEFAULT_TIR_SIMILARITY_THRESHOLD,
    )
    tandem = settings.tandem or TandemRepeatSettings()
    orf = settings.orf or OrfSettings()
    s.setValue(f"{prefix}ltr_min_length", ltr.min_length)
    s.setValue(f"{prefix}ltr_max_length", ltr.max_length)
    s.setValue(f"{prefix}ltr_min_pair_distance", ltr.min_pair_distance)
    s.setValue(f"{prefix}ltr_min_distance_to_end", ltr.min_distance_to_end)
    s.setValue(f"{prefix}ltr_similarity", ltr.similarity_threshold)
    s.setValue(f"{prefix}tir_min_length", tir.min_length)
    s.setValue(f"{prefix}tir_max_length", tir.max_length)
    s.setValue(f"{prefix}tir_min_pair_distance", tir.min_pair_distance)
    s.setValue(f"{prefix}tir_min_distance_to_end", tir.min_distance_to_end)
    s.setValue(f"{prefix}tir_similarity", tir.similarity_threshold)
    s.setValue(f"{prefix}tandem_min_repeat_2", tandem.min_repeat_2)
    s.setValue(f"{prefix}tandem_min_repeat_3", tandem.min_repeat_3)
    s.setValue(f"{prefix}tandem_min_repeat_4", tandem.min_repeat_4)
    s.setValue(f"{prefix}tandem_min_repeat_5", tandem.min_repeat_5)
    s.setValue(f"{prefix}tandem_min_repeat_6", tandem.min_repeat_6)
    s.setValue(f"{prefix}tandem_enable_approximate", tandem.enable_approximate)
    s.setValue(f"{prefix}tandem_atr_min_motif", tandem.atr_min_motif)
    s.setValue(f"{prefix}tandem_atr_max_motif", tandem.atr_max_motif)
    s.setValue(f"{prefix}tandem_atr_min_seed_repeat", tandem.atr_min_seed_repeat)
    s.setValue(f"{prefix}tandem_atr_min_seed_length", tandem.atr_min_seed_length)
    s.setValue(f"{prefix}tandem_atr_min_identity", tandem.atr_min_identity)
    s.setValue(f"{prefix}tandem_atr_max_errors", tandem.atr_max_errors)
    s.setValue(f"{prefix}tandem_atr_max_extend", tandem.atr_max_extend)
    s.setValue(f"{prefix}orf_min_length", orf.min_length)
    s.setValue(f"{prefix}orf_max_length", orf.max_length)
    s.setValue(f"{prefix}orf_alt_starts", orf.alternative_start_codons.strip())
    s.setValue(
        f"{prefix}orf_frames",
        ",".join(str(frame) for frame in normalize_orf_frames(orf.frames)),
    )
    save_annotation_run_options(
        num_threads=settings.num_threads,
        show_process_window=settings.show_process_window,
        close_process_window_when_finished=settings.close_process_window_when_finished,
    )


def load_orf_prediction_settings() -> OrfSettings:
    """ORF settings for Sequence → ORF Prediction (separate from TE annotation)."""
    s = QSettings()
    prefix = "orf_prediction/"
    return OrfSettings(
        min_length=_get_int(s, f"{prefix}min_length", 30),
        max_length=_get_int(s, f"{prefix}max_length", 10_000),
        alternative_start_codons=str(s.value(f"{prefix}alt_starts", "") or ""),
        frames=normalize_orf_frames(s.value(f"{prefix}frames", "")),
    )


def save_orf_prediction_settings(settings: OrfSettings) -> None:
    s = QSettings()
    prefix = "orf_prediction/"
    s.setValue(f"{prefix}min_length", settings.min_length)
    s.setValue(f"{prefix}max_length", settings.max_length)
    s.setValue(f"{prefix}alt_starts", settings.alternative_start_codons.strip())
    s.setValue(
        f"{prefix}frames",
        ",".join(str(frame) for frame in normalize_orf_frames(settings.frames)),
    )
