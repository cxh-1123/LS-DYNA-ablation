# =============================================================================
# run_v60_selected.ps1 -- Launch V6A LAMMPS TTM-MD pilot cases
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
        $RegistryCsv = Join-Path $projectRoot "models\v60_lammps_ttm_md_pilot\v60_case_registry.csv"
    }
    if (-not (Test-Path $RegistryCsv)) {
        Write-Error "Registry CSV not found: $RegistryCsv`nRun scripts\build_v60_lammps_ttm_md_pilot.py first."
        return
    }

    if (-not $LammpsExe) {
        if ($env:LAMMPS_EXE) {
            $LammpsExe = $env:LAMMPS_EXE
        } else {
            $candidates = @("lmp", "lmp_serial", "lmp_mpi", "lammps")
            foreach ($c in $candidates) {
                $cmd = Get-Command $c -ErrorAction SilentlyContinue
                if ($cmd) {
                    $LammpsExe = $cmd.Source
                    break
                }
            }
        }
    }
    if ((-not $LammpsExe) -and $DryRun) {
        $LammpsExe = "<set LAMMPS_EXE before real run>"
    }
    if ((-not $LammpsExe) -and (-not $DryRun)) {
        Write-Error "LAMMPS executable not found. Set `$env:LAMMPS_EXE to your lmp.exe path, or pass -LammpsExe."
        return
    }

    if (-not $PotentialFile) {
        $potentialCandidates = @(
            (Join-Path $projectRoot "potentials\Si.sw"),
            (Join-Path $projectRoot "models\v60_lammps_ttm_md_pilot\Si.sw")
        )
        if ($env:LAMMPS_POTENTIALS) {
            $potentialCandidates += (Join-Path $env:LAMMPS_POTENTIALS "Si.sw")
        }
        foreach ($p in $potentialCandidates) {
            if ($p -and (Test-Path $p)) {
                $PotentialFile = $p
                break
            }
        }
    }

    $rows = Import-Csv $RegistryCsv
    $toRun = @()
    foreach ($r in $rows) {
        if ($OnlyCase -and ($r.name -ne $OnlyCase)) { continue }
        if ($r.run_lammps_hint -match '^(true|1|yes)$') { $toRun += $r }
    }

    Write-Host "==========================================================" -ForegroundColor Cyan
    Write-Host " V6A LAMMPS TTM-MD pilot launcher" -ForegroundColor Cyan
    Write-Host " project root  : $projectRoot"
    Write-Host " lammps exe    : $LammpsExe"
    if ($PotentialFile) {
        Write-Host " Si.sw         : $PotentialFile"
    } else {
        Write-Host " Si.sw         : not found yet" -ForegroundColor Yellow
    }
    Write-Host " registry      : $RegistryCsv"
    Write-Host " ncpu          : $Ncpu"
    if ($DryRun) { Write-Host " DRY-RUN: nothing will actually be executed" -ForegroundColor Yellow }
    Write-Host "==========================================================" -ForegroundColor Cyan

    foreach ($r in $toRun) {
        Write-Host ("  + {0,-30} mode={1,-18} atoms={2}" -f $r.name, $r.mode, $r.atoms_expected)
    }

    if ($DryRun) { return }

    $caseIdx = 0
    foreach ($r in $toRun) {
        $caseIdx++
        $name = $r.name
        $inAbs = Join-Path $projectRoot ($r.input_file -replace '/', '\')
        $teAbs = Join-Path $projectRoot ($r.Te_init_file -replace '/', '\')
        if (-not (Test-Path $inAbs)) {
            Write-Warning "Skipping '$name': input missing at $inAbs"
            continue
        }

        $outDir = Join-Path $projectRoot "results\v60_lammps_ttm_md_pilot\$name"
        if (-not (Test-Path $outDir)) {
            New-Item -ItemType Directory -Path $outDir | Out-Null
        }

        Copy-Item $inAbs (Join-Path $outDir (Split-Path -Leaf $inAbs)) -Force
        if (Test-Path $teAbs) {
            Copy-Item $teAbs (Join-Path $outDir (Split-Path -Leaf $teAbs)) -Force
        }
        if ($PotentialFile -and (Test-Path $PotentialFile)) {
            Copy-Item $PotentialFile (Join-Path $outDir "Si.sw") -Force
        } else {
            Write-Warning "Si.sw was not found. LAMMPS may fail at pair_coeff. Pass -PotentialFile or set `$env:LAMMPS_POTENTIALS."
        }

        Write-Host ""
        Write-Host "[$caseIdx/$($toRun.Count)] Running V6A case '$name' ..." -ForegroundColor Cyan
        Push-Location $outDir
        try {
            $inputName = Split-Path -Leaf $inAbs
            if ($Ncpu -gt 1) {
                $mpi = Get-Command mpiexec -ErrorAction SilentlyContinue
                if ($mpi) {
                    $proc = Start-Process -FilePath $mpi.Source `
                        -ArgumentList "-np", "$Ncpu", "$LammpsExe", "-in", "$inputName" `
                        -Wait -PassThru -NoNewWindow
                } else {
                    Write-Warning "mpiexec not found; running single-process LAMMPS."
                    $proc = Start-Process -FilePath $LammpsExe `
                        -ArgumentList "-in", "$inputName" `
                        -Wait -PassThru -NoNewWindow
                }
            } else {
                $proc = Start-Process -FilePath $LammpsExe `
                    -ArgumentList "-in", "$inputName" `
                    -Wait -PassThru -NoNewWindow
            }
            if ($proc.ExitCode -ne 0) {
                Write-Warning "LAMMPS exit code $($proc.ExitCode) for case '$name'"
            }
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
    Write-Host "V6A LAMMPS TTM-MD pilot run complete." -ForegroundColor Green
} finally {
    Pop-Location
}
