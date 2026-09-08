"""Standalone TE analysis backends (vendored from TEtrimmer proof-curation tools).

Keep package import lightweight: individual backends such as the bundled
aligner should not load the full pandas/plotting curation stack.
"""


def run_batch(*args, **kwargs):
    from repeatviewer.backends.curation import run_batch as implementation

    return implementation(*args, **kwargs)


def run_on_sequence(*args, **kwargs):
    from repeatviewer.backends.curation import run_on_sequence as implementation

    return implementation(*args, **kwargs)


__all__ = ["run_batch", "run_on_sequence"]
