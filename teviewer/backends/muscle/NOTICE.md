# Bundled MUSCLE executable

TEviewer bundles MUSCLE 5.3 executables downloaded from the official release:

https://github.com/rcedgar/muscle/releases/tag/v5.3

MUSCLE was created by Robert C. Edgar. The author states that the MUSCLE
software, including its object code, source code, and documentation, is
donated to the public domain:

https://www.drive5.com/muscle/manual/license.html

Please cite the relevant MUSCLE publication when using results produced by
this software. The software is provided without warranty.

The official macOS executables dynamically link to the GCC runtime. To keep
the bundled aligner self-contained, the macOS folders also contain the
required runtime libraries from the Homebrew GCC 11.5.0 bottles. Their GPL,
LGPL, and GCC Runtime Library Exception notices are included as `GCC-COPYING`,
`GCC-COPYING.LIB`, and `GCC-COPYING.RUNTIME`.
