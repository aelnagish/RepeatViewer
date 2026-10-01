# minimap2 bundled binaries

These directories contain **minimap2 2.31 (r1302)** prebuilt binaries for RepeatViewer:

- **Linux** (`linux/bin/`): official x86_64 build from [lh3/minimap2 v2.31](https://github.com/lh3/minimap2/releases/download/v2.31/minimap2-2.31_x64-linux.tar.bz2)
- **Windows** (`windows/bin/`): unofficial community build from [win-ngs/minimap2-windows-build v2.31-r1302](https://github.com/win-ngs/minimap2-windows-build/releases/download/v2.31-r1302/minimap2-2.31-r1302-windows-x86_64-ucrt64.zip), including `minimap2.exe` and required DLLs (`libwinpthread-1.dll`, `zlib1.dll`)
- **macOS** (`macos/arm64/`, `macos/x86_64/`): Homebrew core bottles built from upstream v2.31

Only executables and any sibling libraries required to run them are vendored here.

minimap2 is distributed under the **MIT License**. Upstream project and documentation: https://lh3.github.io/minimap2
