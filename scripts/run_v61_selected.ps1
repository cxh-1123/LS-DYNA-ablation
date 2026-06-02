# =============================================================================
# run_v61_selected.ps1 -- Launch V6B LAMMPS ttm/mod ablation scout cases
# =============================================================================

[CmdletBinding()]
param(
    [string]$RegistryCsv = "",
    [string]$OnlyCase    = "",
    [string]$LammpsExe   = "",
    [string]$PotentialFile = "",
    [int]$Ncpu           = 1,
    [switch]$DryRun
)

$ErrorActionPreference = "Stop"

$scriptDir   = Split-Path -Parent $MyInvocation.MyCommand.Path
$projectRoot = Resolve-Path (Join-Path $scriptDir "..")
Push-Location $projectRoot
try {
    if (-not $RegistryCsv) {
        $RegistryCsv = Join-Path $projectRoot "models\v61_lammps_ttm_mod_ablation_scan\v61_case_registry.csv"
    }
    if (-not (Test-Path $RegistryCsv)) {
        Write-Error "Registry CSV not found: $RegistryCsv`nRun scripts\build_v61_lammps_ttm_mod_ablation_scan.py first."
        return
    }

    if (-not $LammpsExe) {
        if ($env:LAMMPS_EXE) {
            $LammpsExe = $env:LAMMPS_EXE
        } else {
            $cmd = Get-Command lmp -ErrorAction SilentlyContinue
            if ($cmd) { $LammpsExe = $cmd.Source }
        }
    }
    if ((-not $LammpsExe) -and $DryRun) {
        $LammpsExe = "<set LAMMPS_EXE before real run>"
    }
    if ((-not $LammpsExe) -and (-not $DryRun)) {
        Write-Error "LAMMPS executable not found. Set `$env:LAMMPS_EXE or pass -LammpsExe."
        return
    }

    if (-not $PotentialFile) {
        $candidates = @(
            (Join-Path $projectRoot "potentials\Si.sw"),
            (Join-Path $projectRoot "models\v61_lammps_ttm_mod_ablation_scan\Si.sw")
        )
        if ($env:LAMMPS_POTENTIALS) {
            $candidates += (Join-Path $env:LAMMPS_POTENTIALS "Si.sw")
        }
        foreach ($p in $candidates) {
            if ($p -and (Test-Path $p)) { $PotentialFile = $p; break }
        }
    }

    $rows = Import-Csv $RegistryCsv
    $toRun = @()
    foreach ($r in $rows) {
        if ($OnlyCase -and ($r.name -ne $OnlyCase)) { continue }
        if ($r.run_lammps_hint -match '^(true|1|yes)$') { $toRun += $r }
    }

    Write-Host "==========================================================" -ForegroundColor Cyan
    Write-Host " V6B LAMMPS ttm/mod ablation scout launcher" -ForegroundColor Cyan
    Write-Host " project root  : $projectRoot"
    Write-Host " lammps exe    : $LammpsExe"
    Write-Host " Si.sw         : $PotentialFile"
    Write-Host " registry      : $RegistryCsv"
    Write-Host " ncpu          : $Ncpu"
    if ($DryRun) { Write-Host " DRY-RUN: nothing will actually be executed" -ForegroundColor Yellow }
    Write-Host "==========================================================" -ForegroundColor Cyan
    foreach ($r in $toRun) {
        Write-Host ("  + {0,-22} F={1,5} J/cm2  time={2,5} ps atoms={3}" -f $r.name, $r.fluence_J_cm2, $r.production_ps, $r.atoms_expected)
    }
    if ($DryRun) { return }

    $caseIdx = 0
    foreach ($r in $toRun) {
        $caseIdx++
        $name = $r.name
        $inAbs = Join-Path $projectRoot ($r.input_file -replace '/', '\')
        $ttmAbs = Join-Path $projectRoot ($r.ttm_mod_file -replace '/', '\')
        $outDir = Join-Path $projectRoot "results\v61_lammps_ttm_mod_ablation_scan\$name"
        if (-not (Test-Path $outDir)) { New-Item -ItemType Directory -Path $outDir | Out-Null }

        Copy-Item $inAbs (Join-Path $outDir (Split-Path -Leaf $inAbs)) -Force
        Copy-Item $ttmAbs (Join-Path $outDir (Split-Path -Leaf $ttmAbs)) -Force
        if ($PotentialFile -and (Test-Path $PotentialFile)) {
            Copy-Item $PotentialFile (Join-Path $outDir "Si.sw") -Force
        } else {
            Write-Warning "Si.sw was not found. Pass -PotentialFile."
        }

        Write-Host ""
        Write-Host "[$caseIdx/$($toRun.Count)] Running V6B case '$name' ..." -ForegroundColor Cyan
        Push-Location $outDir
        try {
            $inputName = Split-Path -Leaf $inAbs
            if ($Ncpu -gt 1) {
                $mpi = Get-Command mpiexec -ErrorAction SilentlyContinue
                if ($mpi) {
                    $proc = Start-Process -FilePath $mpi.Source -ArgumentList "-np", "$Ncpu", "$LammpsExe", "-in", "$inputName" -Wait -PassThru -NoNewWindow
                } else {
                    $proc = Start-Process -FilePath $LammpsExe -ArgumentList "-in", "$inputName" -Wait -PassThru -NoNewWindow
                }
            } else {
                $proc = Start-Process -FilePath $LammpsExe -ArgumentList "-in", "$inputName" -Wait -PassThru -NoNewWindow
            }
            if ($proc.ExitCode -ne 0) { Write-Warning "LAMMPS exit code $($proc.ExitCode) for case '$name'" }
        } finally {
            Pop-Location
        }

        $log = Join-Path $outDir "log.lammps"
        if (Test-Path $log) {
            $raw = Get-Content $log -Raw -ErrorAction SilentlyContinue
            if ($raw -match "Total wall time") {
                Write-Host "  -> LAMMPS run completed." -ForegroundColor Green
            } else {
                Write-Host "  -> WARNING: completion marker not found in log.lammps." -ForegroundColor Yellow
            }
        }
    }

    Write-Host ""
    Write-Host "V6B LAMMPS ttm/mod ablation scout run complete." -ForegroundColor Green
} finally {
    Pop-Location
}

