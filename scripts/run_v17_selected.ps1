# =============================================================================
# run_v17_selected.ps1 -- Launch V1.7 local 30 ps LS-DYNA cases
# =============================================================================
#
# Usage (from project root):
#     .\scripts\run_v17_selected.ps1 -DryRun
#     .\scripts\run_v17_selected.ps1 -OnlyCase ep_1p0uj -Ncpu 4
#
# =============================================================================

[CmdletBinding()]
param(
    [string]$RegistryCsv = "",
    [string]$OnlyCase    = "",
    [int]$Ncpu           = 4,
    [switch]$DryRun
)

$ErrorActionPreference = "Stop"

$scriptDir   = Split-Path -Parent $MyInvocation.MyCommand.Path
$projectRoot = Resolve-Path (Join-Path $scriptDir "..")
Push-Location $projectRoot
try {
    if (-not $RegistryCsv) {
        $RegistryCsv = Join-Path $projectRoot "models\v17_30ps_local\v17_case_registry.csv"
    }
    if (-not (Test-Path $RegistryCsv)) {
        Write-Error "Registry CSV not found: $RegistryCsv`nRun scripts\build_v17_30ps_local_mesh.py first."
        return
    }

    Write-Host "==========================================================" -ForegroundColor Cyan
    Write-Host " V1.7 LS-DYNA local 30 ps launcher" -ForegroundColor Cyan
    Write-Host " project root  : $projectRoot"
    Write-Host " registry      : $RegistryCsv"
    Write-Host " ncpu per case : $Ncpu"
    if ($DryRun) { Write-Host " DRY-RUN: nothing will actually be executed" -ForegroundColor Yellow }
    if ($OnlyCase) { Write-Host " filter        : only case '$OnlyCase'" -ForegroundColor Yellow }
    Write-Host "==========================================================" -ForegroundColor Cyan

    $rows = Import-Csv $RegistryCsv
    $toRun  = @()
    $toSkip = @()
    foreach ($r in $rows) {
        if ($OnlyCase -and ($r.name -ne $OnlyCase)) { continue }
        $runFlag = $r.run_lsdyna_final -match '^(true|1|yes)$'
        if ($runFlag) { $toRun += $r } else { $toSkip += $r }
    }

    Write-Host ""
    Write-Host "Cases to RUN  ($($toRun.Count)):" -ForegroundColor Green
    foreach ($r in $toRun) {
        Write-Host ("  + {0,-12} Ep={1,5} uJ  tau={2} ps  nodes={3}  dz_min={4} nm" -f `
            $r.name, $r.Ep_uJ, $r.tau_ps, $r.n_nodes, $r.dz_min_nm)
    }
    if ($toSkip.Count -gt 0) {
        Write-Host ""
        Write-Host "Cases to SKIP ($($toSkip.Count)):" -ForegroundColor DarkGray
        foreach ($r in $toSkip) {
            Write-Host ("  - {0,-12} {1}" -f $r.name, $r.notes) -ForegroundColor DarkGray
        }
    }
    Write-Host ""

    if ($DryRun) {
        Write-Host "Dry-run complete (no LS-DYNA invoked)." -ForegroundColor Yellow
        return
    }

    if ($toRun.Count -eq 0) {
        Write-Host "Nothing to run." -ForegroundColor Yellow
        return
    }

    $caseIdx = 0
    foreach ($r in $toRun) {
        $caseIdx++
        $name = $r.name
        $kAbs = Join-Path $projectRoot ($r.k_file -replace '/', '\')
        if (-not (Test-Path $kAbs)) {
            Write-Warning "Skipping '$name': .k missing at $kAbs"
            continue
        }

        $outDir = Join-Path $projectRoot "results\v17_30ps_local\$name"
        if (-not (Test-Path $outDir)) {
            New-Item -ItemType Directory -Path $outDir | Out-Null
        }

        $kInOut = Join-Path $outDir ("v17_{0}.k" -f $name)
        Copy-Item $kAbs $kInOut -Force

        Write-Host ""
        Write-Host "[$caseIdx/$($toRun.Count)] Running case '$name' ..." -ForegroundColor Cyan
        Write-Host "  input  : $kInOut"
        Write-Host "  output : $outDir"

        Push-Location $outDir
        try {
            $proc = Start-Process -FilePath "lsdyna" `
                -ArgumentList "i=$kInOut", "ncpu=$Ncpu" `
                -Wait -PassThru -NoNewWindow
            if ($proc.ExitCode -ne 0) {
                Write-Warning "lsdyna exit code $($proc.ExitCode) for case '$name'"
            }
        } finally {
            Pop-Location
        }

        $messag = Join-Path $outDir "messag"
        if (Test-Path $messag) {
            $raw = Get-Content $messag -Raw -ErrorAction SilentlyContinue
            $collapsed = ($raw -replace '\s', '').ToLower()
            if ($collapsed -match 'normaltermination') {
                Write-Host "  -> Normal termination detected." -ForegroundColor Green
            } else {
                Write-Host "  -> WARNING: Normal termination NOT found in messag." -ForegroundColor Yellow
            }
        } else {
            Write-Host "  -> WARNING: messag not found." -ForegroundColor Yellow
        }
    }

    Write-Host ""
    Write-Host "V1.7 sweep complete.  Run check_v17_outputs.py next." -ForegroundColor Green
} finally {
    Pop-Location
}
