param(
  [Parameter(Mandatory=$true)][string]$RunRoot,
  [Parameter(Mandatory=$true)][string[]]$CaseNames,
  [ValidateRange(1,16)][int]$Threads=4
)
$ErrorActionPreference='Stop'
$solverRoot=[System.IO.Path]::GetFullPath((Split-Path -Parent $PSScriptRoot))
$inputRoot=Join-Path $solverRoot 'cases\F2_H1_VOLTAGE_SWEEP_20260908'
$exe=Join-Path $solverRoot 'build\F1_FAST_THREADED_0144\TPSolverV2.exe'
$runner=Join-Path $solverRoot 'runs\F1_FAST_HOLE_L0_20260907\evidence\delivery_scripts\run_fast_case.ps1'
foreach($name in $CaseNames) {
  & $runner -RunRoot $RunRoot -CaseName $name -InputRoot $inputRoot -Executable $exe -Threads $Threads
}
Write-Output "F2_WORKER_PASS=$($CaseNames.Count)"
