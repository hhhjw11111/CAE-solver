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
try {
    & python @runnerArgs
    $runExitCode = $LASTEXITCODE
} finally {
    Pop-Location
}
exit $runExitCode
