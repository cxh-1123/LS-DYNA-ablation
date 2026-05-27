# =============================================================================
# run_v42_selected.ps1 -- Launch V4C recoil-drive sweep cases
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
        $RegistryCsv = Join-Path $projectRoot "models\v42_dynamic_recoil_sweep\v42_case_registry.csv"
    }
    if (-not (Test-Path $RegistryCsv)) {
        Write-Error "Registry CSV not found: $RegistryCsv`nRun scripts\build_v42_dynamic_recoil_sweep.py first."
        return
    }

    Write-Host "==========================================================" -ForegroundColor Cyan
    Write-Host " V4C LS-DYNA recoil-drive sweep launcher" -ForegroundColor Cyan
    Write-Host " project root  : $projectRoot"
    Write-Host " registry      : $RegistryCsv"
    Write-Host " ncpu per case : $Ncpu"
    if ($DryRun) { Write-Host " DRY-RUN: nothing will actually be executed" -ForegroundColor Yellow }
    if ($OnlyCase) { Write-Host " filter        : only case '$OnlyCase'" -ForegroundColor Yellow }
    Write-Host "==========================================================" -ForegroundColor Cyan

    $rows = Import-Csv $RegistryCsv
    $toRun = @()
    foreach ($r in $rows) {
        if ($OnlyCase -and ($r.name -ne $OnlyCase)) { continue }
        if ($r.run_lsdyna_final -match '^(true|1|yes)$') { $toRun += $r }
    }

    Write-Host ""
    Write-Host "Cases to RUN  ($($toRun.Count)):" -ForegroundColor Green
    foreach ($r in $toRun) {
        Write-Host ("  + {0,-25} recoil={1,9} m/s  removed={2,5}" -f `
            $r.name, $r.recoil_velocity_m_s, $r.removed_elements)
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

        $outDir = Join-Path $projectRoot "results\v42_dynamic_recoil_sweep\$name"
        if (-not (Test-Path $outDir)) {
            New-Item -ItemType Directory -Path $outDir | Out-Null
        }

        $kInOut = Join-Path $outDir ("v42_{0}.k" -f $name)
        Copy-Item $kAbs $kInOut -Force

        Write-Host ""
        Write-Host "[$caseIdx/$($toRun.Count)] Running V4C case '$name' ..." -ForegroundColor Cyan
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
    Write-Host "V4C recoil sweep complete." -ForegroundColor Green
} finally {
    Pop-Location
}
