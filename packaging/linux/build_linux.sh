#!/usr/bin/env bash
# Build RepeatViewer for Linux (folder + tar.gz).
#
# Usage (from repo root):
#   conda activate <env-with-deps>   # or: source venv/bin/activate
#   bash packaging/linux/build_linux.sh
#
# Output:
#   dist/RepeatViewer/RepeatViewer
#   dist/RepeatViewer-<version>-linux-<arch>.tar.gz
#
# Notes:
#   - Must run on Linux (matching the target CPU arch: x86_64 or aarch64).
#   - Optional conda tools: set File → Setting → External tools (conda).

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"

PYTHON="${PYTHON:-python3}"
if [[ -x "$ROOT/venv/bin/python" ]]; then
  PYTHON="$ROOT/venv/bin/python"
fi

ARCH="$(uname -m)"
VERSION="$("$PYTHON" -c "from repeatviewer import __version__; print(__version__)")"
ARCHIVE_NAME="RepeatViewer-${VERSION}-linux-${ARCH}.tar.gz"
DIST="$ROOT/dist"
BUILD="$ROOT/build"
APP_DIR="$DIST/RepeatViewer"
ARCHIVE_PATH="$DIST/$ARCHIVE_NAME"

echo "==> Using Python: $PYTHON"
echo "==> Version: $VERSION  arch: $ARCH"

if [[ "$(uname -s)" != "Linux" ]]; then
  echo "WARNING: This script is meant to run on Linux. Continuing anyway…" >&2
fi

echo "==> Ensuring build tools (pyinstaller)…"
"$PYTHON" -m pip install -q 'pyinstaller>=6.3'

echo "==> Building app icon (.png)…"
"$PYTHON" packaging/linux/packaging_utils.py || true

echo "==> Cleaning previous build…"
rm -rf "$BUILD/RepeatViewer" "$APP_DIR"
mkdir -p "$DIST"

echo "==> PyInstaller (this can take several minutes)…"
"$PYTHON" -m PyInstaller \
  --noconfirm \
  --clean \
  --distpath "$DIST" \
  --workpath "$BUILD" \
  packaging/linux/RepeatViewer.spec

EXE="$APP_DIR/RepeatViewer"
if [[ ! -x "$EXE" && ! -f "$EXE" ]]; then
  echo "ERROR: $EXE was not created" >&2
  exit 1
fi

echo "==> chmod +x on bundled tools…"
find "$APP_DIR" \
  \( -name 'blast*_linux' -o -name 'makeblastdb_linux' -o -name 'makeprofiledb_linux' \
     -o -name 'rpstblastn_linux' -o -name 'mafft' -o -name 'muscle' \
     -o -path '*/mafft/linux/mafftdir/libexec/*' -o -path '*/mafft/linux/mafftdir/bin/*' \
     -o -name 'RepeatViewer' \) \
  -type f -exec chmod +x {} + 2>/dev/null || true

echo "==> Writing README.txt and desktop entry…"
cat > "$APP_DIR/README.txt" <<'EOF'
RepeatViewer for Linux
======================

1. Unpack this folder anywhere.
2. Run:  ./RepeatViewer

Optional: desktop shortcut
--------------------------
Copy RepeatViewer.desktop to ~/.local/share/applications/ and edit Exec=/Icon=
to absolute paths under this folder.

Optional: RepeatMasker / FastGA / bedtools
-----------------------------------------
Install them in a conda env, then in RepeatViewer open
File → Setting → External tools (conda) and set the env name
(e.g. tetrimmer) or the full path to the env folder.

Bundled BLAST+, MAFFT, and MUSCLE are included for Linux.
EOF

ICON_SRC="$ROOT/packaging/linux/RepeatViewer.png"
if [[ -f "$ICON_SRC" ]]; then
  cp -f "$ICON_SRC" "$APP_DIR/RepeatViewer.png"
fi

cat > "$APP_DIR/RepeatViewer.desktop" <<EOF
[Desktop Entry]
Type=Application
Name=RepeatViewer
Comment=TE curation and visualization
Exec=REPLACE_ME/RepeatViewer
Icon=REPLACE_ME/RepeatViewer.png
Terminal=false
Categories=Science;Biology;
EOF

echo "==> Creating $ARCHIVE_NAME …"
rm -f "$ARCHIVE_PATH"
tar -C "$DIST" -czf "$ARCHIVE_PATH" RepeatViewer

echo
echo "Done."
echo "  App: $EXE"
echo "  Archive: $ARCHIVE_PATH"
echo
echo "Test: \"$EXE\""
echo "Share the .tar.gz; recipients unpack and run ./RepeatViewer"
