param(
  [string]$RunDirectory = ''
)

$ErrorActionPreference = 'Stop'
$solverRoot = [System.IO.Path]::GetFullPath((Split-Path -Parent $PSScriptRoot))
if ([string]::IsNullOrWhiteSpace($RunDirectory)) {
  $RunDirectory = (Get-ChildItem (Join-Path $solverRoot 'runs') -Directory |
    Where-Object Name -like '*G13_ENDPOINTS' | Sort-Object Name -Descending |
    Select-Object -First 1).FullName
}
$runDirectory = [System.IO.Path]::GetFullPath($RunDirectory)
if (-not (Test-Path -LiteralPath $runDirectory)) { throw "G13 run directory not found: $runDirectory" }

$required=@('R060_run_1.txt','R060_run_2.txt','R120_run_1.txt','R120_run_2.txt')
foreach($name in $required){ if(-not (Test-Path (Join-Path $runDirectory $name))){throw "Missing G13 output: $name"} }

$endpointResults=@{}
$allPass=$true
$allReproducible=$true
foreach($endpoint in @('R060','R120')) {
  $one=Get-Content (Join-Path $runDirectory ($endpoint+'_run_1.txt')) -Raw
  $two=Get-Content (Join-Path $runDirectory ($endpoint+'_run_2.txt')) -Raw
  $passed=$one -match 'verification_result=PASS'
  $repro=$one -ceq $two
  $endpointResults[$endpoint]=[pscustomobject]@{Passed=$passed;Reproducible=$repro}
  $allPass=$allPass -and $passed
  $allReproducible=$allReproducible -and $repro
}

foreach($endpoint in @('R060','R120')) {
  foreach($level in @('L0','L1','L2')) {
    $inputDirectory=Join-Path $runDirectory ($endpoint + '_' + $level)
    New-Item -ItemType Directory -Force -Path $inputDirectory | Out-Null
    $sourceDirectory=Join-Path $solverRoot ('cases\validation\G13_endpoint_mesh_sensitivity\{0}\{1}' -f $endpoint,$level)
    Copy-Item (Join-Path $sourceDirectory 'case.nml') $inputDirectory -Force
    Copy-Item (Join-Path $sourceDirectory 'mesh.dat') $inputDirectory -Force
  }
}

$spec=Join-Path $solverRoot 'docs\g13_endpoint_mesh_sensitivity_spec_v1.md'
Copy-Item $spec (Join-Path $runDirectory 'g13_endpoint_mesh_sensitivity_spec_v1.md') -Force
$exe=Join-Path $solverRoot 'build\Debug\TPSolverV2.exe'
& $exe version 2>&1 | Set-Content (Join-Path $runDirectory 'version.txt') -Encoding utf8
Copy-Item (Join-Path $solverRoot 'build\Debug\build_manifest.json') (Join-Path $runDirectory 'build_manifest.json') -Force
$exeHash=(Get-FileHash -Algorithm SHA256 -LiteralPath $exe).Hash.ToLowerInvariant()
$status=if($allPass -and $allReproducible){'VERIFIED'}else{'FAILED'}
$report=@"
## Material Passport

- Artifact ID: TPS-V2-G13-ENDPOINTS-20260902T132814Z_G13_ENDPOINTS
- Artifact type: deterministic endpoint mesh-sensitivity verification
- Verification Status: $status
- Solver: TPSolverV2 0.12.0-g12-radius-matrix
- Frozen protocol: g13_endpoint_mesh_sensitivity_spec_v1.md
- Executable SHA-256: `$exeHash`

## Executed checks

- R060 three-level G10 convergence: $($endpointResults['R060'].Passed)
- R060 exact-repeat output: $($endpointResults['R060'].Reproducible)
- R120 three-level G10 convergence: $($endpointResults['R120'].Passed)
- R120 exact-repeat output: $($endpointResults['R120'].Reproducible)
- Overall verdict: $status

## Interpretation boundary

G13 reuses the frozen G10 convergence thresholds for the two G12 endpoint radii. It establishes endpoint mesh sensitivity under the current discretization, not experimental validity, a monotonic radius law, or publication-level material calibration. The G00–G12 cross-build gate already verifies the solver executable; this G13 artifact is the endpoint mesh-sensitivity supplement.
"@
$report | Set-Content (Join-Path $runDirectory 'validation_report.md') -Encoding utf8

$shaPath=Join-Path $runDirectory 'SHA256SUMS.txt'
if(Test-Path $shaPath){Remove-Item -LiteralPath $shaPath -Force}
Get-ChildItem $runDirectory -Recurse -File | Get-FileHash -Algorithm SHA256 |
  ForEach-Object { '{0}  {1}' -f $_.Hash.ToLowerInvariant(),$_.Path.Substring($runDirectory.Length+1) } |
  Set-Content $shaPath -Encoding ascii
Write-Output ('G13_VERIFICATION_STATUS=' + $status)
Write-Output ('G13_RUN_DIRECTORY=' + $runDirectory)
