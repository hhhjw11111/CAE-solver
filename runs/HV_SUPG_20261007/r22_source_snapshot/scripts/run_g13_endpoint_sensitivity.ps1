param(
  [string]$ExecutablePath = ''
)

$ErrorActionPreference = 'Stop'
$solverRoot = [System.IO.Path]::GetFullPath((Split-Path -Parent $PSScriptRoot))
if ([string]::IsNullOrWhiteSpace($ExecutablePath)) {
  $ExecutablePath = Join-Path $solverRoot 'build\Debug\TPSolverV2.exe'
}
$executableFullPath = [System.IO.Path]::GetFullPath($ExecutablePath)
if (-not (Test-Path -LiteralPath $executableFullPath)) { throw "Executable not found: $executableFullPath" }

& (Join-Path $PSScriptRoot 'generate_g13_meshes.ps1') | Out-Null
$runId = [DateTime]::UtcNow.ToString('yyyyMMddTHHmmssZ') + '_G13_ENDPOINTS'
$runDirectory = Join-Path $solverRoot (Join-Path 'runs' $runId)
New-Item -ItemType Directory -Force -Path $runDirectory | Out-Null
$endpointResults=@{}
$allPass=$true
$allReproducible=$true

foreach($endpoint in @('R060','R120')) {
  $paths=@(
    (Join-Path $solverRoot ('cases\validation\G13_endpoint_mesh_sensitivity\{0}\L0\case.nml' -f $endpoint)),
    (Join-Path $solverRoot ('cases\validation\G13_endpoint_mesh_sensitivity\{0}\L1\case.nml' -f $endpoint)),
    (Join-Path $solverRoot ('cases\validation\G13_endpoint_mesh_sensitivity\{0}\L2\case.nml' -f $endpoint)))
  $runOne=& $executableFullPath verify G10 $paths[0] $paths[1] $paths[2] 2>&1
  if($LASTEXITCODE -ne 0){throw "G13 $endpoint first run failed."}
  $runTwo=& $executableFullPath verify G10 $paths[0] $paths[1] $paths[2] 2>&1
  if($LASTEXITCODE -ne 0){throw "G13 $endpoint second run failed."}
  $runOne | Set-Content (Join-Path $runDirectory ($endpoint + '_run_1.txt')) -Encoding utf8
  $runTwo | Set-Content (Join-Path $runDirectory ($endpoint + '_run_2.txt')) -Encoding utf8
  $passed=(($runOne -join "`n") -match 'verification_result=PASS')
  $repro=(($runOne -join "`n") -ceq ($runTwo -join "`n"))
  $endpointResults[$endpoint]=[pscustomobject]@{Passed=$passed;Reproducible=$repro}
  $allPass=$allPass -and $passed
  $allReproducible=$allReproducible -and $repro
  foreach($level in @('L0','L1','L2')) {
    $inputDirectory=Join-Path $runDirectory ($endpoint + '_' + $level)
    New-Item -ItemType Directory -Force -Path $inputDirectory | Out-Null
    Copy-Item (Join-Path $solverRoot ('cases\validation\G13_endpoint_mesh_sensitivity\{0}\{1}\case.nml' -f $endpoint,$level)) $inputDirectory
    Copy-Item (Join-Path $solverRoot ('cases\validation\G13_endpoint_mesh_sensitivity\{0}\{1}\mesh.dat' -f $endpoint,$level)) $inputDirectory
  }
}

$reportStatus=if($allPass -and $allReproducible){'VERIFIED'}else{'FAILED'}
$report=@"
## Material Passport

- Artifact ID: TPS-V2-G13-ENDPOINTS-$runId
- Artifact type: deterministic endpoint mesh-sensitivity verification
- Verification Status: $reportStatus
- Solver executable: $executableFullPath
- Frozen protocol: g13_endpoint_mesh_sensitivity_spec_v1.md

## Executed checks

- R060 three-level G10 convergence: $($endpointResults['R060'].Passed)
- R120 three-level G10 convergence: $($endpointResults['R120'].Passed)
- R060 exact-repeat output: $($endpointResults['R060'].Reproducible)
- R120 exact-repeat output: $($endpointResults['R120'].Reproducible)
- Overall verdict: $reportStatus

## Interpretation boundary

G13 checks endpoint mesh stability under the existing G10 thresholds. It does not establish experimental validity, a monotonic radius law, or publication-level material calibration.
"@
$report | Set-Content (Join-Path $runDirectory 'validation_report.md') -Encoding utf8
Copy-Item (Join-Path $solverRoot 'docs\g13_endpoint_mesh_sensitivity_spec_v1.md') (Join-Path $runDirectory 'g13_endpoint_mesh_sensitivity_spec_v1.md')
Get-ChildItem $runDirectory -Recurse -File | Get-FileHash -Algorithm SHA256 | ForEach-Object { '{0}  {1}' -f $_.Hash.ToLowerInvariant(),$_.Path.Substring($runDirectory.Length+1) } | Set-Content (Join-Path $runDirectory 'SHA256SUMS.txt') -Encoding ascii
Write-Output ('G13_RUN_DIRECTORY=' + $runDirectory)
Write-Output ('G13_VERIFICATION_STATUS=' + $reportStatus)
