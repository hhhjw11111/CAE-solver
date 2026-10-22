param([Parameter(Mandatory=$true)][string]$RunRoot)
$ErrorActionPreference='Stop'
$root=[System.IO.Path]::GetFullPath($RunRoot)
$inputRoot=[System.IO.Path]::GetFullPath((Join-Path (Split-Path -Parent $PSScriptRoot) 'cases/F1_FAST_HOLE_L0'))
$exe=[System.IO.Path]::GetFullPath((Join-Path (Split-Path -Parent $PSScriptRoot) 'build/F1_FAST_THREADED_0144/TPSolverV2.exe'))
$runner=Join-Path $PSScriptRoot 'run_fast_case.ps1'
$names=@('FAST_H0_LQ_L0','FAST_H0_N0_L0','FAST_H0_NT_L0','FAST_H1_LQ_L0','FAST_H1_N0_L0','FAST_H1_NT_L0')
foreach($name in $names) {
  $dest=Join-Path $root $name
  if (Test-Path -LiteralPath $dest) { throw "Existing fast-case directory: $dest" }
  & $runner -RunRoot $root -CaseName $name -InputRoot $inputRoot -Executable $exe
}
Write-Output 'FAST_HOLE_MATRIX=PASS'
