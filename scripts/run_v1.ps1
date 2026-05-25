# =============================================================================
# run_v1.ps1
# -----------------------------------------------------------------------------
# Run Version 1: 2D axisymmetric silicon wafer pure thermal model.
# Usage (from project root):
#     powershell -ExecutionPolicy Bypass -File scripts\run_v1.ps1
# or:
#     .\scripts\run_v1.ps1
#
# Outputs go to results\v1\ :
#     d3plot / d3plot01 / ...  binary temperature field (open in LS-PrePost)
#     d3hsp                    solver echo, confirm parameters were ingested
#     messag                   warnings / errors
#     glstat                   global statistics
#     tprint                   ASCII nodal/element temperature time history
#
# NOTE: this script is ASCII-only on purpose -- Windows PowerShell parses .ps1
#       files using the system ANSI codepage unless a BOM is present, so any
#       non-ASCII byte will break the parser. Keep this file pure ASCII.
# =============================================================================

$ErrorActionPreference = "Stop"

# Project root (this script lives in scripts\)
$ScriptDir   = Split-Path -Parent $MyInvocation.MyCommand.Path
$ProjectRoot = (Resolve-Path (Join-Path $ScriptDir "..")).Path
Set-Location $ProjectRoot

$ModelRel = "models\v1_thermal_silicon.k"
$ModelAbs = Join-Path $ProjectRoot $ModelRel
if (-not (Test-Path $ModelAbs)) {
    Write-Host "[ERROR] keyword file not found: $ModelAbs" -ForegroundColor Red
    exit 1
}

# Output directory: results\v1\
$OutDir = Join-Path $ProjectRoot "results\v1"
New-Item -ItemType Directory -Force -Path $OutDir | Out-Null

# Copy the .k file into OutDir so LS-DYNA writes results next to it
$ModelLocal = Join-Path $OutDir "v1_thermal_silicon.k"
Copy-Item -Force $ModelAbs $ModelLocal

Write-Host "================================================================"
Write-Host " Version 1 -- 2D axisymmetric silicon thermal model" -ForegroundColor Cyan
Write-Host " model    : $ModelRel"
Write-Host " out dir  : $OutDir"
Write-Host "================================================================"

Push-Location $OutDir
try {
    $kFile = "v1_thermal_silicon.k"
    Write-Host "[INFO] launching LS-DYNA ..." -ForegroundColor Yellow
    $sw = [System.Diagnostics.Stopwatch]::StartNew()
    & lsdyna i=$kFile ncpu=4
    $sw.Stop()
    Write-Host ""
    Write-Host "[INFO] LS-DYNA exit code: $LASTEXITCODE   elapsed: $($sw.Elapsed.TotalSeconds.ToString('F2')) s" -ForegroundColor Yellow
}
finally {
    Pop-Location
}

if ($LASTEXITCODE -ne 0) {
    Write-Host ""
    Write-Host "[WARN] non-zero exit code -- inspect $OutDir\messag and $OutDir\d3hsp" -ForegroundColor Yellow
    exit $LASTEXITCODE
}

# Quick verdict: look for 'Normal termination' near the end of messag.
# LS-DYNA prints the banner with spaces between every letter, e.g.
#   ' N o r m a l    t e r m i n a t i o n'
# so we collapse whitespace before matching.
$messagPath = Join-Path $OutDir "messag"
$ok = $false
if (Test-Path $messagPath) {
    $tail = (Get-Content $messagPath -Tail 50) -join "`n"
    $compact = ($tail -replace '\s+','').ToLower()
    if ($compact -match 'normaltermination') { $ok = $true }
}

if ($ok) {
    Write-Host ""
    Write-Host "[OK] Normal termination -- results in $OutDir" -ForegroundColor Green
    Write-Host "next steps:" -ForegroundColor Green
    Write-Host "    1. open $OutDir\d3plot in LS-PrePost and inspect the T field"
    Write-Host "    2. run  .\.venv\Scripts\python.exe scripts\check_v1_outputs.py"
} else {
    Write-Host ""
    Write-Host "[WARN] did not see 'Normal termination' near end of messag -- inspect log." -ForegroundColor Yellow
}
