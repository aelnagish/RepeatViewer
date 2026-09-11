# Build RepeatViewer for Windows (folder + zip).
#
# Usage (from repo root, in PowerShell):
#   conda activate <env-with-deps>
#   powershell -ExecutionPolicy Bypass -File packaging/windows/build_windows.ps1
#
# Output:
#   dist/RepeatViewer/RepeatViewer.exe
#   dist/RepeatViewer-<version>-windows-x64.zip
#
# Notes:
#   - Must run on Windows (PyInstaller cannot cross-compile a real .exe from macOS).
#   - Optional conda tools (RepeatMasker, FastGA, …) come from WSL2 — set them in
#     File → Setting → External tools (conda / WSL2) after install.

$ErrorActionPreference = "Stop"

$Root = Resolve-Path (Join-Path $PSScriptRoot "..\..")
Set-Location $Root

$Python = if (Test-Path "$Root\venv\Scripts\python.exe") {
    "$Root\venv\Scripts\python.exe"
} elseif ($env:PYTHON) {
    $env:PYTHON
} else {
    "python"
}

$Version = & $Python -c "from repeatviewer import __version__; print(__version__)"
$ZipName = "RepeatViewer-$Version-windows-x64.zip"
$Dist = Join-Path $Root "dist"
$Build = Join-Path $Root "build"
$AppDir = Join-Path $Dist "RepeatViewer"
$ZipPath = Join-Path $Dist $ZipName

Write-Host "==> Using Python: $Python"
Write-Host "==> Version: $Version"

Write-Host "==> Ensuring build tools (pyinstaller)…"
# NumPy 2.4+ raises "cannot load module more than once per process" when a
# frozen app loads _multiarray_umath twice. Keep numpy <2.4 for Windows packs
# and use a PyInstaller with the updated numpy hook.
& $Python -m pip install -q --upgrade "pyinstaller>=6.14.1" "numpy>=1.24,<2.4"

Write-Host "==> Building app icon (.ico)…"
& $Python "packaging\windows\packaging_utils.py"
if ($LASTEXITCODE -ne 0) {
    Write-Warning "Icon generation failed; continuing without custom .ico"
}

Write-Host "==> Cleaning previous build…"
Remove-Item -Recurse -Force -ErrorAction SilentlyContinue `
    (Join-Path $Build "RepeatViewer"), $AppDir
New-Item -ItemType Directory -Force -Path $Dist | Out-Null

Write-Host "==> PyInstaller (this can take several minutes)…"
& $Python -m PyInstaller `
    --noconfirm `
    --clean `
    --distpath $Dist `
    --workpath $Build `
    "packaging\windows\RepeatViewer.spec"
if ($LASTEXITCODE -ne 0) {
    throw "PyInstaller failed with exit code $LASTEXITCODE"
}

$Exe = Join-Path $AppDir "RepeatViewer.exe"
if (-not (Test-Path $Exe)) {
    throw "ERROR: $Exe was not created"
}

Write-Host "==> Writing README.txt in app folder…"
@"
RepeatViewer for Windows
========================

1. Unzip this folder anywhere (or keep the zip contents as-is).
2. Double-click RepeatViewer.exe

Optional: RepeatMasker / FastGA / bedtools via WSL2
----------------------------------------------------
These tools are Linux binaries. Install them in a WSL2 conda env, then in
RepeatViewer open File → Setting → External tools (conda / WSL2) and set the
env name (e.g. tetrimmer) or a Linux path like ~/miniconda3/envs/tetrimmer.

Bundled BLAST+, MAFFT, and MUSCLE run natively on Windows and do not need WSL.
"@ | Set-Content -Encoding UTF8 (Join-Path $AppDir "README.txt")

Write-Host "==> Creating $ZipName …"
Remove-Item -Force -ErrorAction SilentlyContinue $ZipPath
Compress-Archive -Path $AppDir -DestinationPath $ZipPath -Force

Write-Host ""
Write-Host "Done."
Write-Host "  App: $Exe"
Write-Host "  Zip: $ZipPath"
Write-Host ""
Write-Host "Test: start `"$Exe`""
