# =============================================================================
# run_v80_selected.ps1 -- Launch explicit V8 Si ablation cases
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

function Get-FinalThermo {
    param([string]$LogPath)
    $finalStep = $null
    $finalTime = $null
    if (-not (Test-Path $LogPath)) {
        return @{ Step = $null; Time = $null }
    }
    foreach ($line in Get-Content $LogPath) {
        $trim = $line.Trim()
        if ($trim -match '^\d+\s+[-+0-9.Ee]+') {
            $parts = $trim -split '\s+'
            if ($parts.Count -ge 2) {
                $finalStep = $parts[0]
                $finalTime = $parts[1]
            }
        }
    }
    return @{ Step = $finalStep; Time = $finalTime }
}

function Test-RunSucceeded {
    param(
        [int]$ExitCode,
        [string]$LogPath,
        [string]$ExpectedTime
    )
    if ($ExitCode -ne 0) { return $false }
    if (-not (Test-Path $LogPath)) { return $false }
    $raw = Get-Content $LogPath -Raw -ErrorAction SilentlyContinue
    if ($raw -match 'ERROR:') { return $false }
    if ($raw -notmatch 'Total wall time') { return $false }
    $thermo = Get-FinalThermo -LogPath $LogPath
    if ($null -eq $thermo.Time) { return $false }
    $actual = [double]$thermo.Time
    $expected = [double]$ExpectedTime
    if ([math]::Abs($actual - $expected) -gt [math]::Max(0.001, 0.01 * $expected)) {
        return $false
    }
    return $true
}

$scriptDir   = Split-Path -Parent $MyInvocation.MyCommand.Path
$projectRoot = Resolve-Path (Join-Path $scriptDir "..")
Push-Location $projectRoot
try {
    if (-not $RegistryCsv) {
        $RegistryCsv = Join-Path $projectRoot "models\v80_lammps_ttm_md_physical_model\v80_case_registry.csv"
    }
    if (-not (Test-Path $RegistryCsv)) {
        Write-Error "Registry CSV not found: $RegistryCsv`nRun scripts\build_v80_lammps_ttm_md_physical_model.py first."
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
            (Join-Path $projectRoot "models\v80_lammps_ttm_md_physical_model\Si.sw")
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
        if ($OnlyCase) {
            if ($r.name -eq $OnlyCase) { $toRun += $r }
        } elseif ($r.run_lammps_hint -match '^(true|1|yes)$') {
            $toRun += $r
        }
    }
    if ($toRun.Count -eq 0) {
        Write-Error "No V8 case selected. Check -OnlyCase or run_lammps_hint in $RegistryCsv."
        return
    }

    Write-Host "==========================================================" -ForegroundColor Cyan
    Write-Host " V8 Si ablation launcher" -ForegroundColor Cyan
    Write-Host " project root  : $projectRoot"
    Write-Host " lammps exe    : $LammpsExe"
    Write-Host " Si.sw         : $PotentialFile"
    Write-Host " registry      : $RegistryCsv"
    Write-Host " ncpu          : $Ncpu"
    if ($DryRun) { Write-Host " DRY-RUN: nothing will actually be executed" -ForegroundColor Yellow }
    Write-Host "==========================================================" -ForegroundColor Cyan
    foreach ($r in $toRun) {
        Write-Host ("  + {0,-52} mode={1,-30} boundary={2,-5} ttm={3,-5} time={4} ps steps={5}" -f $r.name, $r.case_mode, $r.boundary, $r.uses_fix_ttm_mod, $r.expected_final_time_ps, $r.run_steps)
    }
    if ($DryRun) { return }

    $failed = 0
    $caseIdx = 0
    foreach ($r in $toRun) {
        $caseIdx++
        $name = $r.name
        $inAbs = Join-Path $projectRoot ($r.input_file -replace '/', '\')
        $outDir = Join-Path $projectRoot "results\v80_lammps_ttm_md_physical_model\$name"
        if (-not (Test-Path $outDir)) { New-Item -ItemType Directory -Path $outDir | Out-Null }

        Copy-Item $inAbs (Join-Path $outDir (Split-Path -Leaf $inAbs)) -Force
        if ($r.ttm_mod_file) {
            $ttmAbs = Join-Path $projectRoot ($r.ttm_mod_file -replace '/', '\')
            Copy-Item $ttmAbs (Join-Path $outDir (Split-Path -Leaf $ttmAbs)) -Force
        }
        if ($PotentialFile -and (Test-Path $PotentialFile)) {
            Copy-Item $PotentialFile (Join-Path $outDir "Si.sw") -Force
        } else {
            Write-Warning "Si.sw was not found. Pass -PotentialFile."
        }

        Write-Host ""
        Write-Host "[$caseIdx/$($toRun.Count)] Running V8 case '$name' ..." -ForegroundColor Cyan
        $exitCode = 9999
        Push-Location $outDir
        try {
            $inputName = Split-Path -Leaf $inAbs
            if ($Ncpu -gt 1) {
                $mpi = Get-Command mpiexec -ErrorAction SilentlyContinue
                if ($mpi) {
                    $proc = Start-Process -FilePath $mpi.Source -ArgumentList "-np", "$Ncpu", "$LammpsExe", "-in", "$inputName" -Wait -PassThru -NoNewWindow
                } else {
                    Write-Warning "mpiexec not found; running single executable without MPI launcher."
                    $proc = Start-Process -FilePath $LammpsExe -ArgumentList "-in", "$inputName" -Wait -PassThru -NoNewWindow
                }
            } else {
                $proc = Start-Process -FilePath $LammpsExe -ArgumentList "-in", "$inputName" -Wait -PassThru -NoNewWindow
            }
            $exitCode = $proc.ExitCode
        } finally {
            Pop-Location
        }

        $log = Join-Path $outDir "log.lammps"
        $ok = Test-RunSucceeded -ExitCode $exitCode -LogPath $log -ExpectedTime $r.expected_final_time_ps
        $thermo = Get-FinalThermo -LogPath $log
        if ($ok) {
            Write-Host "  -> run complete" -ForegroundColor Green
            Write-Host "     case       : $name"
            Write-Host "     final Step : $($thermo.Step)"
            Write-Host "     final Time : $($thermo.Time) ps"
            Write-Host "     output dir : $outDir"
        } else {
            $failed++
            Write-Host "  -> failed" -ForegroundColor Red
            Write-Host "     case            : $name" -ForegroundColor Red
            Write-Host "     log file        : $log" -ForegroundColor Red
            Write-Host "     LAMMPS exit code: $exitCode" -ForegroundColor Red
            Write-Host "     final Step      : $($thermo.Step)" -ForegroundColor Red
            Write-Host "     final Time      : $($thermo.Time)" -ForegroundColor Red
        }
    }

    if ($failed -gt 0) {
        Write-Host ""
        Write-Host "V8 run finished with $failed failed case(s)." -ForegroundColor Red
        exit 1
    }
    Write-Host ""
    Write-Host "V8 selected run complete." -ForegroundColor Green
} finally {
    Pop-Location
}

