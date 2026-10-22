<#
Strict preflight for the frozen JE paper runner. Windows may report the image
as TPSOLV~1 when the native executable was launched through GetShortPathNameW.
This wrapper never kills another process and does not alter the frozen runner.
#>
[CmdletBinding()]
param(
    [ValidateSet('probe','angular','radial','joint','supg')]
    [string]$Batch = 'probe',
    [ValidateSet(1,2)]
    [int]$JointLevel = 1,
    [string]$OutputRoot = 'runs/JE_RECONSTRUCTION_20261018/stage4_paper_mesh',
    [switch]$Execute,
    [switch]$PreflightOnly
)
$ErrorActionPreference = 'Stop'
$solverRoot = Split-Path -Parent $PSScriptRoot
$activeSolvers = @(Get-Process | Where-Object {
    $_.ProcessName -match '^(TPSolverV2|TPSOLV~[0-9]+)$'
} | Select-Object Id,ProcessName,Path)
if ($activeSolvers.Count -gt 0) {
    $activeSolvers | ConvertTo-Json -Compress | Write-Output
    throw 'ANOTHER_NATIVE_SOLVER_ACTIVE: complete the existing run before starting a batch.'
}
if ($PreflightOnly) {
    Write-Output 'PREFLIGHT_PASS: no long-name or short-name TPSolver process is active.'
    exit 0
}
$runnerArgs = @(
    (Join-Path $PSScriptRoot 'run_je_reconstruction_mesh_20261018.py'),
    '--verification-scope','paper','--batch',$Batch,
    '--joint-level',"$JointLevel",'--output-root',$OutputRoot
)
if ($Execute) { $runnerArgs += '--execute' }
Push-Location -LiteralPath $solverRoot
$savedPythonUtf8 = $env:PYTHONUTF8
try {
    # Frozen audit JSON contains Chinese paths; propagate UTF-8 to children.
    $env:PYTHONUTF8 = '1'
    # Resolve a workspace junction so the frozen identity uses its original
    # physical script path; do not rewrite any manifest or input file.
    $physicalRunner = & python -X utf8 -c 'from pathlib import Path; import sys; print(Path(sys.argv[1]).resolve())' $runnerArgs[0]
    if ($LASTEXITCODE -ne 0) { throw 'Unable to resolve the runner physical path.' }
    $runnerArgs[0] = $physicalRunner
    & python -X utf8 @runnerArgs
    $runExitCode = $LASTEXITCODE
} finally {
    $env:PYTHONUTF8 = $savedPythonUtf8
    Pop-Location
}
exit $runExitCode
