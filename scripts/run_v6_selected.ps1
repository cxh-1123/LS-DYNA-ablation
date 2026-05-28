# =============================================================================
# run_v6_selected.ps1 -- Launch V6 LAMMPS TTM-MD pilot
# =============================================================================

[CmdletBinding()]
param(
    [string]$RegistryCsv = "",
    [string]$OnlyCase    = "",
    [string]$LammpsExe   = "",
    [switch]$DryRun
)

$ErrorActionPreference = "Stop"

$scriptDir   = Split-Path -Parent $MyInvocation.MyCommand.Path
$projectRoot = Resolve-Path (Join-Path $scriptDir "..")
Push-Location $projectRoot
try {
    if (-not $RegistryCsv) {
        $RegistryCsv = Join-Path $projectRoot "models\v6_lammps_ttm_md_pilot\v6_case_registry.csv"
    }
    if (-not (Test-Path $RegistryCsv)) {
        Write-Error "Registry not found: $RegistryCsv`nRun: python scripts\build_v6_lammps_ttm_input.py"
        return
    }

    $swPath = Join-Path $projectRoot "potentials\Si.sw"
    if (-not (Test-Path $swPath)) {
        Write-Warning "Missing potentials\Si.sw -- copy from your LAMMPS install (potentials/Si.sw)"
    }

    if (-not $LammpsExe) {
        foreach ($cand in @("lmp","lmp_mpi","lammps","lmp.exe")) {
            $p = Get-Command $cand -ErrorAction SilentlyContinue
            if ($p) { $LammpsExe = $p.Source; break }
        }
    }
    if (-not $LammpsExe -and -not $DryRun) {
        Write-Error "LAMMPS executable not found. Set -LammpsExe or add lmp to PATH."
        return
    }

    $rows = Import-Csv $RegistryCsv
    $toRun = @()
    foreach ($r in $rows) {
        if ($OnlyCase -and ($r.name -ne $OnlyCase)) { continue }
        if ($r.run_lammps -match '^(true|1|yes)$') { $toRun += $r }
    }

    Write-Host "==========================================================" -ForegroundColor Cyan
    Write-Host " V6 LAMMPS TTM-MD pilot launcher" -ForegroundColor Cyan
    Write-Host " project root : $projectRoot"
    Write-Host " registry     : $RegistryCsv"
    Write-Host " lammps       : $LammpsExe"
    if ($DryRun) { Write-Host " DRY-RUN" -ForegroundColor Yellow }
    Write-Host "==========================================================" -ForegroundColor Cyan

    foreach ($r in $toRun) {
        Write-Host ("  + {0,-20} Ep={1} uJ  t_end={2} ps" -f $r.name, $r.Ep_uJ, $r.t_end_ps)
    }
    if ($DryRun) { return }

    if (-not $env:OMP_NUM_THREADS) {
        $env:OMP_NUM_THREADS = [Math]::Max(1, [Environment]::ProcessorCount)
        Write-Host " OMP_NUM_THREADS = $env:OMP_NUM_THREADS"
    }

    foreach ($r in $toRun) {
        $name = $r.name
        $caseDir = Join-Path $projectRoot "models\v6_lammps_ttm_md_pilot\$name"
        $inFile  = Join-Path $caseDir "in_$name.lammps"
        if (-not (Test-Path $inFile)) {
            Write-Warning "Skip $name : missing $inFile"
            continue
        }

        $outDir = Join-Path $projectRoot "results\v6_lammps_ttm_md_pilot\$name"
        if (-not (Test-Path $outDir)) { New-Item -ItemType Directory -Path $outDir | Out-Null }

        Push-Location $caseDir
        try {
            Write-Host "`n[RUN] $name from $caseDir" -ForegroundColor Green
            & $LammpsExe -in "in_$name.lammps" -log (Join-Path $outDir "log_$name.lammps")
            $artifacts = Get-ChildItem -File snapshot_*.lammpstrj, traj.lammpstrj, final.data -ErrorAction SilentlyContinue
            if ($artifacts) {
                $artifacts | Copy-Item -Destination $outDir -Force
                Write-Host "[OK] copied $($artifacts.Count) file(s) -> $outDir" -ForegroundColor Green
            }
            if ($LASTEXITCODE -ne 0) {
                Write-Warning "LAMMPS exit code $LASTEXITCODE for $name"
            }
        } finally {
            Pop-Location
        }
    }
} finally {
    Pop-Location
}

Write-Host "`nV6 LAMMPS run pass complete." -ForegroundColor Cyan
