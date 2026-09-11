#!/usr/bin/env bash
# Build RepeatViewer.app and wrap it in a double-clickable .dmg (macOS).
#
# Usage (from repo root):
#   conda activate <env-with-deps>   # or: source venv/bin/activate
#   bash packaging/macos/build_macos_dmg.sh
#
# Output:
#   dist/RepeatViewer.app
#   dist/RepeatViewer-<version>-macos-<arch>.dmg
#
# Notes:
#   - Unsigned builds: users may need right-click → Open the first time.
#   - Apple notarization requires a Developer ID certificate (optional later).

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"

PYTHON="${PYTHON:-python3}"
if [[ -x "$ROOT/venv/bin/python" ]]; then
  PYTHON="$ROOT/venv/bin/python"
fi

ARCH="$(uname -m)"
VERSION="$("$PYTHON" -c "from repeatviewer import __version__; print(__version__)")"
DMG_NAME="RepeatViewer-${VERSION}-macos-${ARCH}.dmg"
DIST="$ROOT/dist"
BUILD="$ROOT/build"
APP="$DIST/RepeatViewer.app"
STAGE="$DIST/dmg_stage"
DMG_PATH="$DIST/$DMG_NAME"

echo "==> Using Python: $PYTHON"
echo "==> Version: $VERSION  arch: $ARCH"

echo "==> Ensuring build tools (pyinstaller)…"
"$PYTHON" -m pip install -q 'pyinstaller>=6.3'

echo "==> Building app icon (.icns)…"
"$PYTHON" packaging/macos/packaging_utils.py || true

echo "==> Cleaning previous build…"
rm -rf "$BUILD/RepeatViewer" "$DIST/RepeatViewer" "$DIST/RepeatViewer.app" "$STAGE"
mkdir -p "$DIST"

echo "==> PyInstaller (this can take several minutes)…"
"$PYTHON" -m PyInstaller \
  --noconfirm \
  --clean \
  --distpath "$DIST" \
  --workpath "$BUILD" \
  packaging/macos/RepeatViewer.spec

if [[ ! -d "$APP" ]]; then
  echo "ERROR: $APP was not created" >&2
  exit 1
fi

# Ensure bundled BLAST/MAFFT binaries are executable.
echo "==> chmod +x on bundled tools…"
find "$APP/Contents/Resources" \
  \( -name 'blast*_mac' -o -name 'makeblastdb_mac' -o -name 'makeprofiledb_mac' \
     -o -name 'rpstblastn_mac' -o -name 'mafft' -o -name 'muscle' \
     -o -path '*/mafft/macos/libexec/*' \) \
  -type f -exec chmod +x {} + 2>/dev/null || true

echo "==> Staging DMG contents…"
rm -rf "$STAGE"
mkdir -p "$STAGE"
# Copy app (not move) so dist/RepeatViewer.app remains for testing.
ditto "$APP" "$STAGE/RepeatViewer.app"
ln -s /Applications "$STAGE/Applications"

# Optional: drop a short readme for first-launch Gatekeeper tip.
cat > "$STAGE/README.txt" <<'EOF'
RepeatViewer for macOS
======================

1. Drag RepeatViewer into Applications.
2. Open it from Applications (or Launchpad).

If macOS says the app is from an unidentified developer:
  - Right-click (or Control-click) RepeatViewer → Open → Open.

This build is not Apple-notarized. Ask the maintainer for a signed
build if your organization requires it.
EOF

echo "==> Creating $DMG_NAME …"
rm -f "$DMG_PATH"
hdiutil create \
  -volname "RepeatViewer" \
  -srcfolder "$STAGE" \
  -ov \
  -format UDZO \
  "$DMG_PATH"

rm -rf "$STAGE"

echo
echo "Done."
echo "  App: $APP"
echo "  DMG: $DMG_PATH"
echo
echo "Test: open \"$APP\""
echo "Or mount the DMG and drag RepeatViewer to Applications."
