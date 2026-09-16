param([Parameter(Mandatory=$true)][string]$RunRoot)
$ErrorActionPreference='Stop'
$root=[System.IO.Path]::GetFullPath($RunRoot)
$caseRunner=Join-Path $PSScriptRoot 'run_f1_case.ps1'
$audit=Join-Path $PSScriptRoot 'audit_f1_results.py'
function Run-One([string]$Name,[string]$InputName='') {
  $dest=Join-Path $root $Name
  if (Test-Path -LiteralPath $dest) {
    $execution=Join-Path $dest 'execution.json'
    if (Test-Path -LiteralPath $execution) {
      $state=Get-Content -LiteralPath $execution -Raw | ConvertFrom-Json
      if ($state.status -eq 'RUNNING') {
        Write-Output "WAIT_ACTIVE_CASE=$Name PID=$($state.pid)"
        do {
          if (-not (Get-Process -Id $state.pid -ErrorAction SilentlyContinue)) { break }
          Start-Sleep -Seconds 15
          $state=Get-Content -LiteralPath $execution -Raw | ConvertFrom-Json
        } while ($state.status -eq 'RUNNING')
      }
    }
    $m=Join-Path $dest 'metrics.txt'
    if ((Test-Path -LiteralPath $m) -and
        ((Get-Content -LiteralPath $m) -contains 'case_acceptance=PASS_NOT_MESH_ACCEPTANCE')) {
      Write-Output "REUSE_COMPLETED=$Name"
      return
    }
    throw "Existing incomplete or failed directory: $dest"
  }
  & $caseRunner -RunRoot $root -CaseName $Name -InputName $InputName
}
& python $audit small SMALL_H1_LQ_L0 SMALL_H1_NT_L0_0142 --run-root $root
if ($LASTEXITCODE -ne 0) { throw 'Small-load audit failed.' }
foreach ($name in @('MESH_H1_NT_L0','MESH_H1_NT_L1','MESH_H1_NT_L2')) { Run-One $name }
& python (Join-Path $PSScriptRoot 'summarize_f1_release.py') --run-root $root --cases MESH_H1_NT_L0 MESH_H1_NT_L1 MESH_H1_NT_L2
if ($LASTEXITCODE -ne 0) { throw 'Independent FE reconstruction failed.' }
& python $audit mesh --run-root $root
if ($LASTEXITCODE -ne 0) { throw 'Mesh audit failed; dependent main matrix stopped.' }
@{MAIN_H1_NT_L2='MESH_H1_NT_L2';reason='Reuse accepted final mesh without another PDE solve'} |
  ConvertTo-Json | Set-Content -LiteralPath (Join-Path $root 'main_case_alias.json') -Encoding utf8
foreach ($name in @('MAIN_H0_LQ_L2','MAIN_H0_N0_L2','MAIN_H0_NT_L2','MAIN_H1_LQ_L2','MAIN_H1_N0_L2')) { Run-One $name }
Run-One 'REPEAT_H1_NT_L2' 'MESH_H1_NT_L2'
$records=foreach ($name in @('metrics.txt','nodes.csv','field_samples.csv')) {
  $a=(Get-FileHash -LiteralPath (Join-Path $root "MESH_H1_NT_L2/$name")).Hash
  $b=(Get-FileHash -LiteralPath (Join-Path $root "REPEAT_H1_NT_L2/$name")).Hash
  [ordered]@{file=$name;original_sha256=$a;repeat_sha256=$b;identical=($a -ceq $b)}
}
$passed=@($records | Where-Object { -not $_.identical }).Count -eq 0
@{passed=$passed;files=@($records)} | ConvertTo-Json -Depth 5 |
  Set-Content -LiteralPath (Join-Path $root 'repeat_comparison.json') -Encoding utf8
if (-not $passed) { throw 'Representative repeat is not byte-identical.' }
Write-Output 'F1_MATRIX_AND_REPRESENTATIVE_REPEAT=PASS'
