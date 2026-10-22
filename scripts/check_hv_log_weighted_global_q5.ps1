param(
  [Parameter(Mandatory=$true)][string]$ProductionBuild,
  [Parameter(Mandatory=$true)][string]$OutputDirectory
)
$ErrorActionPreference='Stop'
$taskRoot=[IO.Path]::GetFullPath((Split-Path -Parent $PSScriptRoot))
$taskBuild=[IO.Path]::GetFullPath($ProductionBuild)
$taskOutput=[IO.Path]::GetFullPath($OutputDirectory)
if(Test-Path -LiteralPath $taskOutput) { throw 'A new output directory is required' }
$taskSource=Join-Path $taskRoot 'tests\hv_log_weighted_global_q5.f90'
$taskModule=Join-Path $taskRoot 'tests\log_weighted_probe_m.f90'
$taskObjects=@(Get-ChildItem -LiteralPath $taskBuild -Filter '*.obj' | Where-Object { $_.Name -ne 'main.obj' } | Sort-Object Name)
if($taskObjects.Count -lt 40) { throw 'Production object set is incomplete' }
New-Item -ItemType Directory -Path $taskOutput | Out-Null
New-Item -ItemType Directory -Path (Join-Path $taskOutput 'mod') | Out-Null
Copy-Item -LiteralPath $taskSource -Destination (Join-Path $taskOutput 'driver_source.f90')
Copy-Item -LiteralPath $taskModule -Destination (Join-Path $taskOutput 'module_source.f90')
$taskExe=Join-Path $taskOutput 'log_weighted_probe.exe'
$taskCompilerOptions='/nologo /Od /check:all /check:noarg_temp_created /traceback /standard-semantics /I"'+(Join-Path $taskBuild 'mod')+'" /module:"'+(Join-Path $taskOutput 'mod')+'" /exe:"'+$taskExe+'"'
$taskLibraries=@('mkl_intel_lp64.lib','mkl_intel_thread.lib','mkl_core.lib') | ForEach-Object { Join-Path 'C:\Program Files (x86)\Intel\oneAPI\mkl\latest\lib' $_ }
$taskLibraries+= 'C:\Program Files (x86)\Intel\oneAPI\compiler\latest\lib\libiomp5md.lib'
$taskQuotedObjects=($taskObjects | ForEach-Object { '"'+$_.FullName+'"' }) -join ' '
$taskQuotedLibraries=($taskLibraries | ForEach-Object { '"'+$_+'"' }) -join ' '
$taskResponse=Join-Path $taskOutput 'compile.rsp'
($taskCompilerOptions+' "'+$taskModule+'" "'+$taskSource+'" '+$taskQuotedObjects+' '+$taskQuotedLibraries) | Set-Content -LiteralPath $taskResponse -Encoding ascii
$taskCommand='call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64 vs2022 >nul && ifx @"'+$taskResponse+'"'
Push-Location $taskOutput
try {
  & $env:ComSpec /d /s /c $taskCommand *> (Join-Path $taskOutput 'build.txt')
  if($LASTEXITCODE -ne 0) { throw 'Log-weighted prototype harness build failed' }
} finally { Pop-Location }
$taskOldPath=$env:PATH
try {
  $env:PATH='C:\Program Files (x86)\Intel\oneAPI\compiler\latest\bin;'+$env:PATH
  & $taskExe *> (Join-Path $taskOutput 'stdout.txt')
  $taskExit=$LASTEXITCODE
} finally { $env:PATH=$taskOldPath }
$taskManifest=[ordered]@{
  source=$taskSource;source_sha256=(Get-FileHash -LiteralPath $taskSource -Algorithm SHA256).Hash.ToLowerInvariant()
  prototype_module=$taskModule;prototype_module_sha256=(Get-FileHash -LiteralPath $taskModule -Algorithm SHA256).Hash.ToLowerInvariant()
  production_build=$taskBuild;production_executable_sha256=(Get-FileHash -LiteralPath (Join-Path $taskBuild 'TPSolverV2.exe') -Algorithm SHA256).Hash.ToLowerInvariant()
  harness_executable_sha256=(Get-FileHash -LiteralPath $taskExe -Algorithm SHA256).Hash.ToLowerInvariant()
  options=$taskCompilerOptions;exit_code=$taskExit
  objects=@($taskObjects | ForEach-Object { [ordered]@{name=$_.Name;sha256=(Get-FileHash -LiteralPath $_.FullName -Algorithm SHA256).Hash.ToLowerInvariant()} })
}
$taskManifest | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath (Join-Path $taskOutput 'manifest.json') -Encoding utf8
Get-Content -LiteralPath (Join-Path $taskOutput 'stdout.txt') -Tail 4
if($taskExit -ne 0) { throw "Log-weighted prototype check failed: $taskExit" }
