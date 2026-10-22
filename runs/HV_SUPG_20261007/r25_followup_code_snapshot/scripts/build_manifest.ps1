param(
  [Parameter(Mandatory=$true)]
  [string]$ExecutablePath,
  [string]$OutputPath
)

$ErrorActionPreference = 'Stop'
$solverRoot = [System.IO.Path]::GetFullPath((Split-Path -Parent $PSScriptRoot))
$executableFullPath = [System.IO.Path]::GetFullPath($ExecutablePath)
if (-not $executableFullPath.StartsWith($solverRoot + [System.IO.Path]::DirectorySeparatorChar, [System.StringComparison]::OrdinalIgnoreCase)) {
  throw "Executable is outside solver_v2: $executableFullPath"
}
if (-not (Test-Path -LiteralPath $executableFullPath)) {
  throw "Executable not found: $executableFullPath"
}

# Verify the executable itself, source declaration and VERSION agree. A stale
# text label must never identify a different compiled solver as this version.
$taskVersion=(Get-Content -LiteralPath (Join-Path $solverRoot 'VERSION') -Raw).Trim()
$taskSourceText=Get-Content -LiteralPath (Join-Path $solverRoot 'src\core\build_info_m.f90') -Raw
$taskSourceMatch=[regex]::Match($taskSourceText,"SOLVER_VERSION\s*=\s*'([^']+)'")
if(-not $taskSourceMatch.Success -or $taskSourceMatch.Groups[1].Value -ne $taskVersion) {
  throw 'VERSION and build_info_m.f90 disagree'
}
$taskOldPath=$env:PATH
try {
  $env:PATH='C:\Program Files (x86)\Intel\oneAPI\compiler\latest\bin;'+$env:PATH
  $taskVersionLines=@(& $executableFullPath version)
  if($LASTEXITCODE -ne 0) { throw 'Cannot read executable version' }
} finally { $env:PATH=$taskOldPath }
$taskRuntime=@{}
foreach($taskLine in $taskVersionLines) {
  $taskParts=$taskLine -split '=',2
  if($taskParts.Count -eq 2) { $taskRuntime[$taskParts[0]]=$taskParts[1] }
}
if($taskRuntime['version'] -ne $taskVersion) { throw 'Executable and VERSION disagree' }

$files = @()
$files += Get-ChildItem -LiteralPath (Join-Path $solverRoot 'src') -Recurse -File -Filter '*.f90'
$files += Get-ChildItem -LiteralPath (Join-Path $solverRoot 'scripts') -File -Filter '*.ps1'
$files += Get-ChildItem -LiteralPath (Join-Path $solverRoot 'docs') -File -Filter '*.md'
$files += Get-Item -LiteralPath (Join-Path $solverRoot 'TPSolverV2.sln')
$files += Get-Item -LiteralPath (Join-Path $solverRoot 'TPSolverV2.vfproj')
$files += Get-Item -LiteralPath (Join-Path $solverRoot 'VERSION')
$files += Get-Item -LiteralPath $executableFullPath

$records = foreach ($file in ($files | Sort-Object FullName)) {
  $relative = [System.IO.Path]::GetRelativePath($solverRoot, $file.FullName)
  [pscustomobject]@{
    path = $relative.Replace('\','/')
    bytes = $file.Length
    modified_utc = $file.LastWriteTimeUtc.ToString('o')
    sha256 = (Get-FileHash -Algorithm SHA256 -LiteralPath $file.FullName).Hash.ToLowerInvariant()
  }
}

$manifest = [ordered]@{
  schema = 'TPSolverV2-build-manifest-v1'
  solver_version = $taskRuntime['version']
  version_verified_against_source_and_executable = $true
  theory_version = $taskRuntime['theory_version']
  architecture_id = $taskRuntime['architecture_id']
  generated_utc = [DateTime]::UtcNow.ToString('o')
  files = $records
}
$manifestPath = Join-Path (Split-Path -Parent $executableFullPath) 'build_manifest.json'
if($OutputPath) {
  $manifestPath=[IO.Path]::GetFullPath($OutputPath)
  if(Test-Path -LiteralPath $manifestPath) { throw 'Supplemental manifest output must be new' }
}
$manifest | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath $manifestPath -Encoding utf8
Write-Output "MANIFEST_OK=$manifestPath"
