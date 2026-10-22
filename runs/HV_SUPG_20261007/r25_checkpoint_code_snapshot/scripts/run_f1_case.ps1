param(
  [Parameter(Mandatory=$true)][string]$RunRoot,
  [Parameter(Mandatory=$true)][string]$CaseName,
  [string]$InputName = ''
)
$ErrorActionPreference='Stop'
$root=[System.IO.Path]::GetFullPath($RunRoot)
$freeze=Get-Content -LiteralPath (Join-Path $root 'evidence/release.json') -Raw | ConvertFrom-Json
$exe=Join-Path $root 'evidence/TPSolverV2_main.exe'
if ((Get-FileHash -LiteralPath $exe).Hash -ne $freeze.executable_sha256) {
  throw 'Executable differs from release manifest.'
}
if (-not $InputName) { $InputName=$CaseName }
$inputPath=Join-Path $root "evidence/inputs/$InputName.nml"
$dest=Join-Path $root $CaseName
if (Test-Path -LiteralPath $dest) { throw "Do not overwrite existing case: $dest" }
New-Item -ItemType Directory -Path $dest | Out-Null
Copy-Item -LiteralPath $inputPath -Destination (Join-Path $dest 'input.nml')
$env:TPS_PARDISO_STORAGE=$freeze.storage
$env:MKL_PARDISO_OOC_MAX_CORE_SIZE=[string]$freeze.ooc_max_core_mb
$env:MKL_PARDISO_OOC_MAX_SWAP_SIZE=[string]$freeze.ooc_max_swap_mb
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$scratch=Join-Path $root "ooc_scratch/$CaseName"
New-Item -ItemType Directory -Path $scratch | Out-Null
$env:MKL_PARDISO_OOC_PATH=$scratch
$env:MKL_PARDISO_OOC_FILE_NAME=Join-Path $scratch 'factor'
$record=[ordered]@{
  case=$CaseName; input_name=$InputName; started_utc=[DateTime]::UtcNow.ToString('o')
  executable_sha256=$freeze.executable_sha256; solver_version=$freeze.solver_version
  input_sha256=(Get-FileHash -LiteralPath $inputPath).Hash
  storage=$freeze.storage; ooc_max_core_mb=$freeze.ooc_max_core_mb
  ooc_max_swap_mb=$freeze.ooc_max_swap_mb; scratch=$scratch
  status='STARTING'
}
$record | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath (Join-Path $dest 'execution.json') -Encoding utf8
$process=Start-Process -FilePath $exe -ArgumentList @('run-f1',('"'+$inputPath+'"'),('"'+$dest+'"')) `
  -WorkingDirectory $root -WindowStyle Hidden -PassThru `
  -RedirectStandardOutput (Join-Path $dest 'stdout.txt') -RedirectStandardError (Join-Path $dest 'stderr.txt')
$record.pid=$process.Id
$record.status='RUNNING'
$record | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath (Join-Path $dest 'execution.json') -Encoding utf8
Write-Output "START_CASE=$CaseName PID=$($process.Id)"
while (-not $process.WaitForExit(30000)) {
  $process.Refresh()
  $record.last_monitor_utc=[DateTime]::UtcNow.ToString('o')
  $record.elapsed_seconds=([DateTime]::UtcNow-[DateTime]::Parse($record.started_utc).ToUniversalTime()).TotalSeconds
  $record.working_set_bytes=$process.WorkingSet64
  $record.cpu_seconds=$process.TotalProcessorTime.TotalSeconds
  $record.stdout_bytes=(Get-Item -LiteralPath (Join-Path $dest 'stdout.txt')).Length
  $record.long_runtime_advisory=($record.elapsed_seconds -gt 43200)
  $record | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath (Join-Path $dest 'execution.json') -Encoding utf8
}
$process.WaitForExit()
$record.exit_code=$process.ExitCode
$record.finished_utc=[DateTime]::UtcNow.ToString('o')
$metric=Join-Path $dest 'metrics.txt'
$accepted=($process.ExitCode -eq 0) -and (Test-Path -LiteralPath $metric)
if ($accepted) { $accepted=(Get-Content -LiteralPath $metric) -contains 'case_acceptance=PASS_NOT_MESH_ACCEPTANCE' }
$record.status=if ($accepted) { 'PASS_NOT_MESH_ACCEPTANCE' } else { 'FAILED' }
$record | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath (Join-Path $dest 'execution.json') -Encoding utf8
$outputs=Get-ChildItem -LiteralPath $dest -File | Where-Object Name -ne 'output_hashes.json' | ForEach-Object {
  [ordered]@{name=$_.Name;bytes=$_.Length;sha256=(Get-FileHash -LiteralPath $_.FullName).Hash}
}
ConvertTo-Json -InputObject @($outputs) -Depth 4 | Set-Content -LiteralPath (Join-Path $dest 'output_hashes.json') -Encoding utf8
if (-not $accepted) { throw "F1 case failed: $CaseName; see preserved stdout, stderr and metrics." }
Write-Output "PASS_CASE=$CaseName"
