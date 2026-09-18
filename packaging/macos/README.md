# macOS packaging (DMG)

Double-clickable **RepeatViewer.app** inside a `.dmg`.

## Build

From the repo root, in an environment where `pip install -e .` already works:

```bash
bash packaging/macos/build_macos_dmg.sh
```

Outputs:

- `dist/RepeatViewer.app`
- `dist/RepeatViewer-<version>-macos-<arch>.dmg` (e.g. `arm64` or `x86_64`)

## Install (end users)

1. Open the DMG.
2. Drag **RepeatViewer** to **Applications**.
3. Launch from Applications.

Unsigned builds: first launch may require **right-click → Open**.

## What’s included

- PySide6 UI + Qt WebEngine (HTML previews)
- Bundled **macOS-only** BLAST+, MAFFT, MUSCLE (linux/windows tool trees are not packed)
- Flexidot (run in-process when frozen — no second Dock icon)
- AliView runtime (`aliview.jar` + launcher); installer scripts are omitted
- Annotation deps **orfipy** (ORF detection) and **pyhmmer** (HMM / Dfam search),
  collected in full so frozen imports work

## Signing / notarization (optional)

For Gatekeeper-clean distribution you need an Apple Developer ID:

```bash
codesign --deep --force --options runtime \
  --sign "Developer ID Application: …" \
  dist/RepeatViewer.app
# then notarize with notarytool and staple
```

Re-run DMG creation after signing.
