param(
    [string]$LsPrePostExe = $env:LSPREPOST_EXE,
    [string]$Case = "all"
)

$ErrorActionPreference = "Stop"
$root = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
Set-Location $root

if (-not $LsPrePostExe) {
    $cmd = Get-Command lsprepost.exe -ErrorAction SilentlyContinue
    if ($cmd) {
        $LsPrePostExe = $cmd.Source
    }
}

if (-not $LsPrePostExe -or -not (Test-Path $LsPrePostExe)) {
    throw "Set `$env:LSPREPOST_EXE to your LS-PrePost executable path, or pass -LsPrePostExe `"C:\path\to\lsprepost.exe`"."
}

python scripts\make_v44_lsprepost_cfiles.py --case $Case

$cfiles = if ($Case -eq "all") {
    Get-ChildItem post\v44_lsprepost_real_field_extract -Recurse -Filter export_fields.cfile
} else {
    Get-ChildItem "post\v44_lsprepost_real_field_extract\$Case" -Filter export_fields.cfile
}

foreach ($cfile in $cfiles) {
    Write-Host "[RUN] $($cfile.FullName)"
    & $LsPrePostExe "c=$($cfile.FullName)" -nographics
}

python scripts\summarize_v44_lsprepost_fields.py --case $Case

