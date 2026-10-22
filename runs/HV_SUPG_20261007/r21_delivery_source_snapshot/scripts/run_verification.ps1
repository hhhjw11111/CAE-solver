param(
  [string]$ExecutablePath = ''
)

$ErrorActionPreference = 'Stop'
$solverRoot = [System.IO.Path]::GetFullPath((Split-Path -Parent $PSScriptRoot))
if ([string]::IsNullOrWhiteSpace($ExecutablePath)) {
  $ExecutablePath = Join-Path $solverRoot 'build\Debug\TPSolverV2.exe'
}
$executableFullPath = [System.IO.Path]::GetFullPath($ExecutablePath)
if (-not $executableFullPath.StartsWith($solverRoot + [System.IO.Path]::DirectorySeparatorChar, [System.StringComparison]::OrdinalIgnoreCase)) {
  throw "Executable is outside solver_v2: $executableFullPath"
}
if (-not (Test-Path -LiteralPath $executableFullPath)) {
  throw "Executable not found: $executableFullPath"
}

$runId = [DateTime]::UtcNow.ToString('yyyyMMddTHHmmssZ') + '_G00_G12_A1J_L01'
$runDirectory = Join-Path $solverRoot (Join-Path 'runs' $runId)
New-Item -ItemType Directory -Force -Path $runDirectory | Out-Null

$versionOutput = & $executableFullPath version 2>&1
if ($LASTEXITCODE -ne 0) { throw 'Version command failed.' }
$versionOutput | Set-Content -LiteralPath (Join-Path $runDirectory 'version.txt') -Encoding utf8
$versionLine = $versionOutput | Where-Object { "$_" -match '^version=' } | Select-Object -First 1
if ([string]::IsNullOrWhiteSpace($versionLine)) { throw 'Solver did not report its version.' }
$solverVersion = ("$versionLine" -split '=', 2)[1].Trim()

$casePath = Join-Path $solverRoot 'cases\validation\V01_patch\case.nml'
$inspectOutput = & $executableFullPath inspect $casePath 2>&1
if ($LASTEXITCODE -ne 0) { throw 'Inspect command failed.' }
$inspectOutput | Set-Content -LiteralPath (Join-Path $runDirectory 'inspect.txt') -Encoding utf8

$g09PlainCasePath = Join-Path $solverRoot 'cases\validation\G09_center_hole\H00_plain\case.nml'
$g09HoleCasePath = Join-Path $solverRoot 'cases\validation\G09_center_hole\H02_r090\case.nml'
$g10CasePaths = @(
  (Join-Path $solverRoot 'cases\validation\G10_mesh_convergence\L0\case.nml'),
  (Join-Path $solverRoot 'cases\validation\G10_mesh_convergence\L1\case.nml'),
  (Join-Path $solverRoot 'cases\validation\G10_mesh_convergence\L2\case.nml')
)
$g11PlainCasePath = Join-Path $solverRoot 'cases\validation\G11_matched_control\H00_L2\case.nml'
$g11HoleCasePath = $g10CasePaths[2]
$g12PlainCasePath = $g11PlainCasePath
$g12CasePaths = @(
  (Join-Path $solverRoot 'cases\validation\G12_radius_matrix\R060\case.nml'),
  (Join-Path $solverRoot 'cases\validation\G12_radius_matrix\R090\case.nml'),
  (Join-Path $solverRoot 'cases\validation\G12_radius_matrix\R120\case.nml')
)
$verificationIds = @('V01','G02','G03','G04','L01','G05','G06','G07','G08','G09','G10','G11','G12','G12C','A1J')
$verificationResults = @{}
$cliOutputs = @{}
$allChecksPass = (($inspectOutput -join "`n") -match 'inspect_result=PASS')
$allReproducible = $true
foreach ($verificationId in $verificationIds) {
  if ($verificationId -eq 'G09') {
    $runOne = & $executableFullPath verify $verificationId $g09PlainCasePath $g09HoleCasePath 2>&1
  } elseif ($verificationId -eq 'G10') {
    $runOne = & $executableFullPath verify $verificationId $g10CasePaths[0] $g10CasePaths[1] $g10CasePaths[2] 2>&1
  } elseif ($verificationId -eq 'G11') {
    $runOne = & $executableFullPath verify $verificationId $g11PlainCasePath $g11HoleCasePath 2>&1
  } elseif ($verificationId -eq 'G12') {
    $runOne = & $executableFullPath verify $verificationId $g12PlainCasePath $g12CasePaths[0] $g12CasePaths[1] $g12CasePaths[2] 2>&1
  } else {
    $runOne = & $executableFullPath verify $verificationId 2>&1
  }
  if ($LASTEXITCODE -ne 0) { throw "First $verificationId verification failed." }
  $runOne | Set-Content -LiteralPath (Join-Path $runDirectory ($verificationId + '_run_1.txt')) -Encoding utf8

  if ($verificationId -eq 'G09') {
    $runTwo = & $executableFullPath verify $verificationId $g09PlainCasePath $g09HoleCasePath 2>&1
  } elseif ($verificationId -eq 'G10') {
    $runTwo = & $executableFullPath verify $verificationId $g10CasePaths[0] $g10CasePaths[1] $g10CasePaths[2] 2>&1
  } elseif ($verificationId -eq 'G11') {
    $runTwo = & $executableFullPath verify $verificationId $g11PlainCasePath $g11HoleCasePath 2>&1
  } elseif ($verificationId -eq 'G12') {
    $runTwo = & $executableFullPath verify $verificationId $g12PlainCasePath $g12CasePaths[0] $g12CasePaths[1] $g12CasePaths[2] 2>&1
  } else {
    $runTwo = & $executableFullPath verify $verificationId 2>&1
  }
  if ($LASTEXITCODE -ne 0) { throw "Second $verificationId verification failed." }
  $runTwo | Set-Content -LiteralPath (Join-Path $runDirectory ($verificationId + '_run_2.txt')) -Encoding utf8

  $checkPass = (($runOne -join "`n") -match 'verification_result=PASS')
  $checkReproducible = (($runOne -join "`n") -ceq ($runTwo -join "`n"))
  $cliOutputs[$verificationId] = ($runOne -join "`n")
  $verificationResults[$verificationId] = [pscustomobject]@{
    passed = $checkPass
    reproducible = $checkReproducible
  }
  $allChecksPass = $allChecksPass -and $checkPass
  $allReproducible = $allReproducible -and $checkReproducible
}

$productionCases = @(
  [pscustomobject]@{ id='G08_REFERENCE'; casePath=(Join-Path $solverRoot 'cases\validation\G08_reference_bar\case.nml') },
  [pscustomobject]@{ id='G09_H00'; casePath=$g09PlainCasePath },
  [pscustomobject]@{ id='G09_H02'; casePath=$g09HoleCasePath }
)
$productionResults = @{}
foreach ($productionCase in $productionCases) {
  $resultOne = Join-Path $runDirectory ('production_' + $productionCase.id + '_cli_1.dat')
  $resultTwo = Join-Path $runDirectory ('production_' + $productionCase.id + '_cli_2.dat')
  $stdoutOne = & $executableFullPath run $productionCase.casePath $resultOne 2>&1
  if ($LASTEXITCODE -ne 0) { throw "First $($productionCase.id) production run failed." }
  $stdoutOne | Set-Content -LiteralPath (Join-Path $runDirectory ('production_' + $productionCase.id + '_cli_1.txt')) -Encoding utf8
  $stdoutTwo = & $executableFullPath run $productionCase.casePath $resultTwo 2>&1
  if ($LASTEXITCODE -ne 0) { throw "Second $($productionCase.id) production run failed." }
  $stdoutTwo | Set-Content -LiteralPath (Join-Path $runDirectory ('production_' + $productionCase.id + '_cli_2.txt')) -Encoding utf8
  $casePassed = (($stdoutOne -join "`n") -match 'run_result=PASS')
  $caseReproducible = (($stdoutOne -join "`n") -ceq ($stdoutTwo -join "`n")) -and
    ((Get-FileHash -Algorithm SHA256 -LiteralPath $resultOne).Hash -ceq
     (Get-FileHash -Algorithm SHA256 -LiteralPath $resultTwo).Hash)
  $productionResults[$productionCase.id] = [pscustomobject]@{
    passed=$casePassed
    reproducible=$caseReproducible
    crossBuild=$false
    cliOutput=($stdoutOne -join "`n")
    resultOne=$resultOne
  }
  $allChecksPass = $allChecksPass -and $casePassed
  $allReproducible = $allReproducible -and $caseReproducible
}

$vswhere = 'C:\Program Files (x86)\Microsoft Visual Studio\Installer\vswhere.exe'
if (-not (Test-Path -LiteralPath $vswhere)) { throw 'vswhere.exe was not found.' }
$devenv = & $vswhere -latest -products * -find 'Common7\IDE\devenv.com' | Select-Object -First 1
if ([string]::IsNullOrWhiteSpace($devenv)) { throw 'devenv.com was not found.' }
$ideBuildOutput = & $devenv (Join-Path $solverRoot 'TPSolverV2.sln') /Rebuild 'Debug|x64' 2>&1
$ideBuildExitCode = $LASTEXITCODE
$ideBuildOutput | Set-Content -LiteralPath (Join-Path $runDirectory 'visual_studio_rebuild.log') -Encoding utf8
if ($ideBuildExitCode -ne 0) { throw "Visual Studio rebuild failed with exit code $ideBuildExitCode." }
$ideBuildClean = (($ideBuildOutput -join "`n") -match 'TPSolverV2 - 0 error\(s\), 0 warning\(s\)')

$ideExecutable = Join-Path $solverRoot 'x64\Debug\TPSolverV2.exe'
if (-not (Test-Path -LiteralPath $ideExecutable)) { throw 'Visual Studio executable was not produced.' }
$ideVersionOutput = & $ideExecutable version 2>&1
if ($LASTEXITCODE -ne 0) { throw 'Visual Studio version command failed.' }
$ideVersionOutput | Set-Content -LiteralPath (Join-Path $runDirectory 'IDE_version.txt') -Encoding utf8
$allCrossBuild = (($ideVersionOutput -join "`n") -ceq ($versionOutput -join "`n"))
foreach ($verificationId in $verificationIds) {
  if ($verificationId -eq 'G09') {
    $ideRun = & $ideExecutable verify $verificationId $g09PlainCasePath $g09HoleCasePath 2>&1
  } elseif ($verificationId -eq 'G10') {
    $ideRun = & $ideExecutable verify $verificationId $g10CasePaths[0] $g10CasePaths[1] $g10CasePaths[2] 2>&1
  } elseif ($verificationId -eq 'G11') {
    $ideRun = & $ideExecutable verify $verificationId $g11PlainCasePath $g11HoleCasePath 2>&1
  } elseif ($verificationId -eq 'G12') {
    $ideRun = & $ideExecutable verify $verificationId $g12PlainCasePath $g12CasePaths[0] $g12CasePaths[1] $g12CasePaths[2] 2>&1
  } else {
    $ideRun = & $ideExecutable verify $verificationId 2>&1
  }
  if ($LASTEXITCODE -ne 0) { throw "Visual Studio $verificationId verification failed." }
  $ideRun | Set-Content -LiteralPath (Join-Path $runDirectory ('IDE_' + $verificationId + '.txt')) -Encoding utf8
  $allCrossBuild = $allCrossBuild -and
    (($ideRun -join "`n") -match 'verification_result=PASS') -and
    (($ideRun -join "`n") -ceq $cliOutputs[$verificationId])
}
foreach ($productionCase in $productionCases) {
  $ideResult = Join-Path $runDirectory ('production_' + $productionCase.id + '_ide.dat')
  $ideOutput = & $ideExecutable run $productionCase.casePath $ideResult 2>&1
  if ($LASTEXITCODE -ne 0) { throw "Visual Studio $($productionCase.id) production run failed." }
  $ideOutput | Set-Content -LiteralPath (Join-Path $runDirectory ('IDE_production_' + $productionCase.id + '.txt')) -Encoding utf8
  $caseCrossBuild = (($ideOutput -join "`n") -ceq $productionResults[$productionCase.id].cliOutput) -and
    ((Get-FileHash -Algorithm SHA256 -LiteralPath $ideResult).Hash -ceq
     (Get-FileHash -Algorithm SHA256 -LiteralPath $productionResults[$productionCase.id].resultOne).Hash)
  $productionResults[$productionCase.id].crossBuild=$caseCrossBuild
  $allCrossBuild = $allCrossBuild -and $caseCrossBuild
}
$crossBuildHashes = @(
  ('cli_executable_sha256=' + (Get-FileHash -Algorithm SHA256 -LiteralPath $executableFullPath).Hash.ToLowerInvariant())
  ('ide_executable_sha256=' + (Get-FileHash -Algorithm SHA256 -LiteralPath $ideExecutable).Hash.ToLowerInvariant())
)
$crossBuildHashes | Set-Content -LiteralPath (Join-Path $runDirectory 'cross_build_hashes.txt') -Encoding ascii
$passed = $allChecksPass -and $allReproducible -and $allCrossBuild -and $ideBuildClean

Copy-Item -LiteralPath (Join-Path (Split-Path -Parent $executableFullPath) 'build_manifest.json') -Destination $runDirectory
$inputDirectory = Join-Path $runDirectory 'input_snapshot'
New-Item -ItemType Directory -Force -Path $inputDirectory | Out-Null
Copy-Item -LiteralPath $casePath -Destination $inputDirectory
Copy-Item -LiteralPath (Join-Path (Split-Path -Parent $casePath) 'mesh.dat') -Destination $inputDirectory
foreach ($productionCase in $productionCases) {
  $productionInputDirectory = Join-Path $inputDirectory $productionCase.id
  New-Item -ItemType Directory -Force -Path $productionInputDirectory | Out-Null
  Copy-Item -LiteralPath $productionCase.casePath -Destination $productionInputDirectory
  Copy-Item -LiteralPath (Join-Path (Split-Path -Parent $productionCase.casePath) 'mesh.dat') -Destination $productionInputDirectory
}
foreach ($g10CasePath in $g10CasePaths) {
  $g10Name = Split-Path -Leaf (Split-Path -Parent $g10CasePath)
  $g10InputDirectory = Join-Path $inputDirectory ('G10_' + $g10Name)
  New-Item -ItemType Directory -Force -Path $g10InputDirectory | Out-Null
  Copy-Item -LiteralPath $g10CasePath -Destination $g10InputDirectory
  Copy-Item -LiteralPath (Join-Path (Split-Path -Parent $g10CasePath) 'mesh.dat') -Destination $g10InputDirectory
}
$g11InputDirectory = Join-Path $inputDirectory 'G11_H00_L2'
New-Item -ItemType Directory -Force -Path $g11InputDirectory | Out-Null
Copy-Item -LiteralPath $g11PlainCasePath -Destination $g11InputDirectory
Copy-Item -LiteralPath (Join-Path (Split-Path -Parent $g11PlainCasePath) 'mesh.dat') -Destination $g11InputDirectory
foreach ($g12CasePath in $g12CasePaths) {
  $g12Name = Split-Path -Leaf (Split-Path -Parent $g12CasePath)
  $g12InputDirectory = Join-Path $inputDirectory ('G12_' + $g12Name)
  New-Item -ItemType Directory -Force -Path $g12InputDirectory | Out-Null
  Copy-Item -LiteralPath $g12CasePath -Destination $g12InputDirectory
  Copy-Item -LiteralPath (Join-Path (Split-Path -Parent $g12CasePath) 'mesh.dat') -Destination $g12InputDirectory
}
$specificationDirectory = Join-Path $runDirectory 'verification_spec'
New-Item -ItemType Directory -Force -Path $specificationDirectory | Out-Null
Copy-Item -LiteralPath (Join-Path $solverRoot 'docs\g05_m0_benchmark_spec_v1.md') -Destination $specificationDirectory
Copy-Item -LiteralPath (Join-Path $solverRoot 'docs\g06_v07_block_jacobian_spec_v1.md') -Destination $specificationDirectory
Copy-Item -LiteralPath (Join-Path $solverRoot 'docs\g07_true_residual_newton_spec_v1.md') -Destination $specificationDirectory
Copy-Item -LiteralPath (Join-Path $solverRoot 'docs\g08_continuation_production_spec_v1.md') -Destination $specificationDirectory
Copy-Item -LiteralPath (Join-Path $solverRoot 'docs\g08_continuation_production_spec_v2.md') -Destination $specificationDirectory
Copy-Item -LiteralPath (Join-Path $solverRoot 'docs\g09_center_hole_geometry_spec_v1.md') -Destination $specificationDirectory
Copy-Item -LiteralPath (Join-Path $solverRoot 'docs\g10_mesh_convergence_spec_v1.md') -Destination $specificationDirectory
Copy-Item -LiteralPath (Join-Path $solverRoot 'docs\g11_matched_control_spec_v1.md') -Destination $specificationDirectory
Copy-Item -LiteralPath (Join-Path $solverRoot 'docs\g11_matched_control_spec_v2.md') -Destination $specificationDirectory
Copy-Item -LiteralPath (Join-Path $solverRoot 'docs\g11_matched_control_findings_v1.md') -Destination $runDirectory
Copy-Item -LiteralPath (Join-Path $solverRoot 'docs\g11_matched_control_findings_v2.md') -Destination $runDirectory
Copy-Item -LiteralPath (Join-Path $solverRoot 'docs\g12_radius_matrix_spec_v1.md') -Destination $specificationDirectory
Copy-Item -LiteralPath (Join-Path $solverRoot 'docs\a1j_production_joule_flux_audit_spec_v1.md') -Destination $specificationDirectory
Copy-Item -LiteralPath (Join-Path $solverRoot 'docs\g12_radius_matrix_findings_v1.md') -Destination $runDirectory

$verificationStatus = if ($passed) { 'VERIFIED' } else { 'FAILED' }
$report = @"
## Material Passport

- Artifact ID: TPS-V2-G00-G12-A1J-L01-$runId
- Artifact type: deterministic code verification
- Verification Status: $verificationStatus
- Solver: TPSolverV2 $solverVersion
- Theory baseline: 2.1
- Architecture: TPS-V2-ARCH-001

## Executed checks

- G00 source/executable manifest: present
- Native input inspection: $((($inspectOutput -join "`n") -match 'inspect_result=PASS'))
- V01 HEX20 geometry and interpolation: $($verificationResults['V01'].passed)
- G02 material-point response and tangent: $($verificationResults['G02'].passed)
- G03 element residual and Jacobian: $($verificationResults['G03'].passed)
- G04 global CSR, scatter, Dirichlet and Jacobian: $($verificationResults['G04'].passed)
- L01 static PARDISO, scaling and manufactured global solve: $($verificationResults['L01'].passed)
- G05 M0 V02-V06 analytic solutions and convergence: $($verificationResults['G05'].passed)
- G05 error integration independent of solver quadrature: 5x5x5 Gauss
- G06 M2 full global 5x7 block-Jacobian audit: $($verificationResults['G06'].passed)
- G07 true-residual safeguarded Newton and discrete conservation: $($verificationResults['G07'].passed)
- G08 voltage/Joule continuation, cutback rollback and restart: $($verificationResults['G08'].passed)
- G09 matched plain/hole geometry, topology, Jacobian and continuation: $($verificationResults['G09'].passed)
- G10 three-level hole-edge mesh convergence with independent 5x5x3 metrics: $($verificationResults['G10'].passed)
- G11 matched-boundary H00/H02 comparison with 8x8 cut-cell integration and electron, hole, thermal and displacement metrics: $($verificationResults['G11'].passed)
- G12 60/90/120 nm radius matrix with matched H00 controls: $($verificationResults['G12'].passed)
- G12C per-radius H00 denominator and annulus guard: $($verificationResults['G12C'].passed)
- A1J production-path Joule heat, independent flux and negative-control audit: $($verificationResults['A1J'].passed)
- G08 reference production case execution: $($productionResults['G08_REFERENCE'].passed)
- G09 H00 production case execution: $($productionResults['G09_H00'].passed)
- G09 H02 production case execution: $($productionResults['G09_H02'].passed)
- All three production results exact-repeat: $allReproducible
- V01 exact-repeat comparison: $($verificationResults['V01'].reproducible)
- G02 exact-repeat comparison: $($verificationResults['G02'].reproducible)
- G03 exact-repeat comparison: $($verificationResults['G03'].reproducible)
- G04 exact-repeat comparison: $($verificationResults['G04'].reproducible)
- L01 exact-repeat comparison: $($verificationResults['L01'].reproducible)
- G05 exact-repeat comparison: $($verificationResults['G05'].reproducible)
- G06 exact-repeat comparison: $($verificationResults['G06'].reproducible)
- G07 exact-repeat comparison: $($verificationResults['G07'].reproducible)
- G08 exact-repeat comparison: $($verificationResults['G08'].reproducible)
- G09 exact-repeat comparison: $($verificationResults['G09'].reproducible)
- G10 exact-repeat comparison: $($verificationResults['G10'].reproducible)
- G11 exact-repeat comparison: $($verificationResults['G11'].reproducible)
- G12 exact-repeat comparison: $($verificationResults['G12'].reproducible)
- G12C exact-repeat comparison: $($verificationResults['G12C'].reproducible)
- A1J exact-repeat comparison: $($verificationResults['A1J'].reproducible)
- Visual Studio rebuild has 0 errors and 0 warnings: $ideBuildClean
- CLI and Visual Studio outputs match exactly: $allCrossBuild
- G08 reference CLI/Visual Studio production result match: $($productionResults['G08_REFERENCE'].crossBuild)
- G09 H00 CLI/Visual Studio production result match: $($productionResults['G09_H00'].crossBuild)
- G09 H02 CLI/Visual Studio production result match: $($productionResults['G09_H02'].crossBuild)

## Verdict

$verificationStatus
"@
$report | Set-Content -LiteralPath (Join-Path $runDirectory 'validation_report.md') -Encoding utf8

$artifactFiles = Get-ChildItem -LiteralPath $runDirectory -Recurse -File |
  Where-Object { $_.Name -ne 'SHA256SUMS.txt' } | Sort-Object FullName
$artifactHashes = foreach ($file in $artifactFiles) {
  '{0}  {1}' -f (Get-FileHash -Algorithm SHA256 -LiteralPath $file.FullName).Hash.ToLowerInvariant(),
    [System.IO.Path]::GetRelativePath($runDirectory, $file.FullName).Replace('\','/')
}
$artifactHashes | Set-Content -LiteralPath (Join-Path $runDirectory 'SHA256SUMS.txt') -Encoding ascii

if (-not $passed) {
  throw "G00-G12 plus A1J and L01 verification failed. Evidence: $runDirectory"
}
Write-Output "VERIFICATION_OK=$runDirectory"
