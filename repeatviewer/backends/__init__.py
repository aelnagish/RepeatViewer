"""Standalone TE analysis backends (vendored from TEtrimmer proof-curation tools).

Keep package import lightweight: individual backends such as the bundled
aligner should not load the full pandas/plotting curation stack.
"""

import os as _os

# NCBI BLAST+ sends a usage report to NCBI over the network every time a
# program starts, and waits for it: about 1 s per blastn/makeblastdb launch
# versus 0.05 s with reporting off. Feature detection launches BLAST once or
# twice per sequence, so that wait was most of its run time. Child processes
# inherit this; a user who sets the variable themselves keeps their value.
_os.environ.setdefault("BLAST_USAGE_REPORT", "false")


def run_batch(*args, **kwargs):
    from repeatviewer.backends.curation import run_batch as implementation

    return implementation(*args, **kwargs)


def run_on_sequence(*args, **kwargs):
    from repeatviewer.backends.curation import run_on_sequence as implementation

    return implementation(*args, **kwargs)


__all__ = ["run_batch", "run_on_sequence"]
