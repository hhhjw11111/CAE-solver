param(
  [Parameter(Mandatory=$true)][string]$RunRoot,
  [Parameter(Mandatory=$true)][string]$Executable
)
$ErrorActionPreference='Stop'
$solverRoot=[System.IO.Path]::GetFullPath((Split-Path -Parent $PSScriptRoot))
$root=[System.IO.Path]::GetFullPath($RunRoot)
$inputRoot=Join-Path $solverRoot 'cases\F3_RADIUS_SWEEP_20260908'
$exe=[System.IO.Path]::GetFullPath($Executable)
$runner=Join-Path $PSScriptRoot 'run_fast_case.ps1'
$names=@(
  'F3_R02_LQ_L0','F3_R02_N0_L0','F3_R02_NT_L0',
  'F3_R04_LQ_L0','F3_R04_N0_L0','F3_R04_NT_L0'
)
foreach($name in $names) {
  $dest=Join-Path $root $name
  if (Test-Path -LiteralPath $dest) { throw "Existing F3 case directory: $dest" }
  & $runner -RunRoot $root -CaseName $name -InputRoot $inputRoot -Executable $exe -Threads 4
}
Write-Output 'F3_RADIUS_SWEEP=PASS_NOT_MESH_ACCEPTANCE'
