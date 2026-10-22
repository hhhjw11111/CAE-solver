$ErrorActionPreference = 'Stop'
$solverRoot = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot '../../..'))
$objects = @('kinds_m','constants_m','status_m','field_layout_m','mesh_types_m','model_types_m',
  'equilibrium_m','reference_material_m','transport_m','recombination_m','thermo_piezo_m',
  'material_point_m','hex20_shape_m','hex20_quadrature_m','dof_map_m',
  'transport_supg_m','element_kernel_m','material_tangent_check_m')
foreach ($revision in @('R13','R14')) {
  $buildPath = Join-Path $solverRoot "build/HV_SUPG_20261008_$revision"
  $outputPath = Join-Path $PSScriptRoot $revision
  New-Item -ItemType Directory -Path $outputPath -ErrorAction Stop | Out-Null
  $objectArgs = ($objects | ForEach-Object { '"' + (Join-Path $buildPath ($_ + '.obj')) + '"' }) -join ' '
  $command = 'call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64 vs2022 >nul && ifx /nologo /O2 /Qopenmp /standard-semantics /I:"' + (Join-Path $buildPath 'mod') + '" /exe:benchmark.exe "' + (Join-Path $PSScriptRoot 'compare_kernel.f90') + '" ' + $objectArgs
  Push-Location $outputPath
  try {
    & $env:ComSpec /d /s /c $command
    if ($LASTEXITCODE -ne 0) { throw "Benchmark link failed: $revision" }
    $env:OMP_NUM_THREADS='1';$env:MKL_NUM_THREADS='1'
    $env:PATH='C:\Program Files (x86)\Intel\oneAPI\compiler\latest\bin;' + $env:PATH
    & .\benchmark.exe > timing.txt
    if ($LASTEXITCODE -ne 0) { throw "Benchmark execution failed: $revision" }
  } finally { Pop-Location }
}
