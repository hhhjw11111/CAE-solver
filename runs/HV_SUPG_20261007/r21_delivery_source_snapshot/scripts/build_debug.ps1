param(
  [switch]$Clean,
  [switch]$Optimize,
  [switch]$Threaded,
  [switch]$ParallelAssembly,
  [switch]$QuietArrayTemporaries,
  [string]$BuildSubdirectory = 'build\Debug'
)

$ErrorActionPreference = 'Stop'
$solverRoot = [System.IO.Path]::GetFullPath((Split-Path -Parent $PSScriptRoot))
$buildDirectory = [System.IO.Path]::GetFullPath((Join-Path $solverRoot $BuildSubdirectory))
$moduleDirectory = Join-Path $buildDirectory 'mod'
$executablePath = Join-Path $buildDirectory 'TPSolverV2.exe'
$setvars = 'C:\Program Files (x86)\Intel\oneAPI\setvars.bat'
$mklLibraryDirectory = 'C:\Program Files (x86)\Intel\oneAPI\mkl\latest\lib'
$compilerLibraryDirectory = 'C:\Program Files (x86)\Intel\oneAPI\compiler\latest\lib'

if (-not $buildDirectory.StartsWith($solverRoot + [System.IO.Path]::DirectorySeparatorChar, [System.StringComparison]::OrdinalIgnoreCase)) {
  throw "Unsafe build directory: $buildDirectory"
}
if (-not (Test-Path -LiteralPath $setvars)) {
  throw "Intel oneAPI environment script not found: $setvars"
}
$mklLibraries = if ($Threaded) {
  @('mkl_intel_lp64.lib','mkl_intel_thread.lib','mkl_core.lib')
} else {
  @('mkl_intel_lp64.lib','mkl_sequential.lib','mkl_core.lib')
}
foreach ($mklLibrary in $mklLibraries) {
  if (-not (Test-Path -LiteralPath (Join-Path $mklLibraryDirectory $mklLibrary))) {
    throw "Required static MKL library not found: $mklLibrary"
  }
}
$extraLibraries = @()
if ($Threaded) {
  $ompLibrary = Join-Path $compilerLibraryDirectory 'libiomp5md.lib'
  if (-not (Test-Path -LiteralPath $ompLibrary)) {
    throw "Threaded OpenMP runtime not found: $ompLibrary"
  }
  $extraLibraries += $ompLibrary
}
if ($Clean -and (Test-Path -LiteralPath $buildDirectory)) {
  Remove-Item -LiteralPath $buildDirectory -Recurse -Force
}
New-Item -ItemType Directory -Force -Path $moduleDirectory | Out-Null

$relativeSources = @(
  'src\core\kinds_m.f90',
  'src\core\constants_m.f90',
  'src\core\status_m.f90',
  'src\core\field_layout_m.f90',
  'src\mesh\mesh_types_m.f90',
  'src\mesh\boundary_topology_m.f90',
  'src\core\model_types_m.f90',
  'src\core\build_info_m.f90',
  'src\io\input_schema_m.f90',
  'src\physics\equilibrium_m.f90',
  'src\physics\reference_material_m.f90',
  'src\physics\zno_candidate_m.f90',
  'src\physics\transport_m.f90',
  'src\physics\recombination_m.f90',
  'src\physics\thermo_piezo_m.f90',
  'src\physics\material_point_m.f90',
  'src\mesh\hex20_shape_m.f90',
  'src\mesh\hex20_quadrature_m.f90',
  'src\mesh\dof_map_m.f90',
  'src\mesh\mesh_checks_m.f90',
  'src\io\native_input_m.f90',
  'src\fem\transport_supg_m.f90',
  'src\fem\element_kernel_m.f90',
  'src\fem\sparse_pattern_m.f90',
  'src\fem\dirichlet_m.f90',
  'src\fem\global_assembly_m.f90',
  'src\linear\pardiso_solver_m.f90',
  'src\nonlinear\scaling_m.f90',
  'src\nonlinear\newton_solver_m.f90',
  'src\nonlinear\temperature_control_m.f90',
  'src\post\balance_audit_m.f90',
  'src\post\independent_flux_audit_m.f90',
  'src\post\field_metrics_m.f90',
  'src\nonlinear\continuation_solver_m.f90',
  'src\problems\reference_bar_problem_m.f90',
  'src\post\result_writer_m.f90',
  'src\app\production_runner_m.f90',
  'src\app\f1_case_runner_m.f90',
  'src\verification\material_tangent_check_m.f90',
  'src\verification\element_tangent_check_m.f90',
  'src\verification\supg_check_m.f90',
  'src\verification\global_assembly_check_m.f90',
  'src\verification\global_block_jacobian_check_m.f90',
  'src\verification\linear_solver_check_m.f90',
  'src\verification\m0_benchmark_m.f90',
  'src\verification\a1_joint_check_m.f90',
  'src\verification\a2_hole_check_m.f90',
  'src\verification\analytic_benchmarks_m.f90',
  'src\verification\newton_solver_check_m.f90',
  'src\verification\continuation_solver_check_m.f90',
  'src\verification\temperature_control_check_m.f90',
  'src\verification\center_hole_check_m.f90',
  'src\verification\mesh_convergence_check_m.f90',
  'src\verification\matched_control_check_m.f90',
  'src\verification\radius_matrix_check_m.f90',
  'src\verification\production_joule_check_m.f90',
  'src\verification\verification_driver_m.f90',
  'src\app\main.f90'
)
$sourceArguments = ($relativeSources | ForEach-Object {
  '"' + (Join-Path $solverRoot $_) + '"'
}) -join ' '

$optimization = if ($Optimize) { '/O2' } else { '/Od' }
$runtimeChecks = if ($Optimize) { '/check:none' } else { '/check:all' }
# Array-temporary diagnostics are performance warnings, emitted once per node
# by formatted section I/O. Keep all correctness checks in the opt-in quiet run.
if ($QuietArrayTemporaries -and -not $Optimize) { $runtimeChecks += ' /check:noarg_temp_created' }
if($ParallelAssembly -and -not $Threaded) { throw 'ParallelAssembly requires Threaded for the OpenMP runtime.' }
$parallelOption = if($ParallelAssembly) { '/Qopenmp' } else { '' }
$compilerOptions = @(
  '/nologo',
  $optimization,
  '/debug:full',
  '/traceback',
  $runtimeChecks,
  $parallelOption,
  '/warn:all',
  '/warn:interfaces',
  '/stand:f18',
  '/standard-semantics',
  ('/module:"' + $moduleDirectory + '"'),
  ('/exe:"' + $executablePath + '"')
) -join ' '

$libraryArguments = ($mklLibraries | ForEach-Object {
  '"' + (Join-Path $mklLibraryDirectory $_) + '"'
}) -join ' '
$libraryArguments = (($libraryArguments + ' ' + (($extraLibraries | ForEach-Object { '"' + $_ + '"' }) -join ' ')).Trim())
$compilerCommand = "ifx $compilerOptions $sourceArguments $libraryArguments"
$command = 'call "' + $setvars + '" intel64 vs2022 >nul && ' + $compilerCommand

Push-Location $buildDirectory
try {
  & $env:ComSpec /d /s /c $command
  if ($LASTEXITCODE -ne 0) {
    throw "ifx build failed with exit code $LASTEXITCODE"
  }
} finally {
  Pop-Location
}

& (Join-Path $PSScriptRoot 'build_manifest.ps1') -ExecutablePath $executablePath
[ordered]@{ compiler='Intel ifx'; options=$compilerOptions; optimized=[bool]$Optimize; threaded=[bool]$Threaded; libraries=$mklLibraries+$extraLibraries; sources=$relativeSources } | ConvertTo-Json -Depth 4 | Set-Content -LiteralPath (Join-Path $buildDirectory 'compiler_options.json') -Encoding utf8
Write-Output "BUILD_OK=$executablePath"
