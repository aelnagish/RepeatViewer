# Bundled MUSCLE executable

RepeatViewer bundles MUSCLE 5.3 executables downloaded from the official release:

https://github.com/rcedgar/muscle/releases/tag/v5.3

MUSCLE was created by Robert C. Edgar. The author states that the MUSCLE
software, including its object code, source code, and documentation, is
donated to the public domain:

https://www.drive5.com/muscle/manual/license.html

Please cite the relevant MUSCLE publication when using results produced by
this software. The software is provided without warranty.

The official Apple Silicon macOS executable dynamically links to the GCC
runtime. To keep the bundled aligner self-contained, its folder also contains the
required runtime libraries from the Homebrew GCC 11.5.0 bottles. Their GPL,
LGPL, and GCC Runtime Library Exception notices are included as `GCC-COPYING`,
`GCC-COPYING.LIB`, and `GCC-COPYING.RUNTIME`.
