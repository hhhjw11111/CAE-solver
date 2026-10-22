param(
  [string]$ExecutablePath = 'solver_v2\build\F1_MAIN_0142\TPSolverV2.exe'
)
$ErrorActionPreference='Stop'
$solverRoot=[System.IO.Path]::GetFullPath((Split-Path -Parent $PSScriptRoot))
$runRoot=Join-Path $solverRoot 'runs\F1_20260906'
$inputRoot=Join-Path $solverRoot 'cases\F1_20260906'
$exe=[System.IO.Path]::GetFullPath($ExecutablePath)
$frozenExe=Join-Path $runRoot 'evidence\TPSolverV2_main_0142.exe'
if ((Get-FileHash -LiteralPath $exe).Hash -ne (Get-FileHash -LiteralPath $frozenExe).Hash) {
  throw 'Executable differs from the frozen F1 evidence copy.'
}
$env:TPS_PARDISO_STORAGE='OOC'
$env:MKL_PARDISO_OOC_MAX_CORE_SIZE='2048'
$env:MKL_PARDISO_OOC_MAX_SWAP_SIZE='0'
$scratch=Join-Path $runRoot 'ooc_scratch'
New-Item -ItemType Directory -Path $scratch -Force | Out-Null
$env:MKL_PARDISO_OOC_PATH=$scratch
$env:MKL_PARDISO_OOC_FILE_NAME=Join-Path $scratch 'f1_factor'
function Run-One([string]$CaseName) {
  $dest=Join-Path $runRoot $CaseName
  if (Test-Path -LiteralPath $dest) {
    $metric=Join-Path $dest 'metrics.txt'
    if ((Test-Path -LiteralPath $metric) -and
        ((Get-Content -LiteralPath $metric) -contains 'case_acceptance=PASS_NOT_MESH_ACCEPTANCE')) {
      Write-Output "REUSE_COMPLETED=$CaseName"
      return
    }
    throw "Existing incomplete or failed case must be diagnosed, not overwritten: $CaseName"
  }
  New-Item -ItemType Directory -Path $dest -ErrorAction Stop | Out-Null
  Write-Output "START_CASE=$CaseName"
  & $exe run-f1 (Join-Path $inputRoot ($CaseName+'.nml')) $dest 2>&1 |
    Tee-Object -FilePath (Join-Path $dest 'stdout.txt')
  if ($LASTEXITCODE -ne 0) { throw "F1 case failed: $CaseName" }
}
& python (Join-Path $PSScriptRoot 'audit_f1_results.py') small SMALL_H1_LQ_L0 SMALL_H1_NT_L0_0142
if ($LASTEXITCODE -ne 0) { throw 'Small-load comparison has not passed.' }
foreach ($caseName in @('MESH_H1_NT_L0','MESH_H1_NT_L1','MESH_H1_NT_L2')) { Run-One $caseName }
& python (Join-Path $PSScriptRoot 'audit_f1_results.py') mesh
if ($LASTEXITCODE -ne 0) { throw 'The final mesh comparison has not passed; dependent main cases are not started.' }
foreach ($caseName in @('MAIN_H0_LQ_L2','MAIN_H0_N0_L2','MAIN_H0_NT_L2','MAIN_H1_LQ_L2','MAIN_H1_N0_L2')) {
  Run-One $caseName
}
$repeat=Join-Path $runRoot 'REPEAT_H1_NT_L2'
if (Test-Path -LiteralPath $repeat) { throw 'Representative repeat output already exists; inspect before proceeding.' }
New-Item -ItemType Directory -Path $repeat -ErrorAction Stop | Out-Null
& $exe run-f1 (Join-Path $inputRoot 'MESH_H1_NT_L2.nml') $repeat 2>&1 |
  Tee-Object -FilePath (Join-Path $repeat 'stdout.txt')
if ($LASTEXITCODE -ne 0) { throw 'Representative repeat failed.' }
foreach ($name in @('metrics.txt','nodes.csv','field_samples.csv')) {
  $original=Join-Path (Join-Path $runRoot 'MESH_H1_NT_L2') $name
  $repeated=Join-Path $repeat $name
  if ((Get-FileHash -LiteralPath $original).Hash -ne (Get-FileHash -LiteralPath $repeated).Hash) {
    throw "Representative repeat differs: $name"
  }
}
Write-Output 'F1_MATRIX_AND_REPRESENTATIVE_REPEAT=PASS'
