param(
  [Parameter(Mandatory=$true)]
  [string]$ExecutablePath
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
  solver_version = (Get-Content -LiteralPath (Join-Path $solverRoot 'VERSION') -Raw).Trim()
  theory_version = '2.1'
  architecture_id = 'TPS-V2-ARCH-001'
  generated_utc = [DateTime]::UtcNow.ToString('o')
  files = $records
}
$manifestPath = Join-Path (Split-Path -Parent $executableFullPath) 'build_manifest.json'
$manifest | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath $manifestPath -Encoding utf8
Write-Output "MANIFEST_OK=$manifestPath"
