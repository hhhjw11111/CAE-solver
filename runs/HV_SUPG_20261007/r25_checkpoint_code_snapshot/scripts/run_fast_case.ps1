param(
  [Parameter(Mandatory=$true)][string]$RunRoot,
  [Parameter(Mandatory=$true)][string]$CaseName,
  [Parameter(Mandatory=$true)][string]$InputRoot,
  [Parameter(Mandatory=$true)][string]$Executable,
  [ValidateRange(1,16)][int]$Threads=4,
  [string]$InputName=''
)
$ErrorActionPreference='Stop'
$root=[System.IO.Path]::GetFullPath($RunRoot)
if (-not $InputName) { $InputName=$CaseName }
$input=[System.IO.Path]::GetFullPath((Join-Path $InputRoot ($InputName+'.nml')))
$exe=[System.IO.Path]::GetFullPath($Executable)
$dest=Join-Path $root $CaseName
if (Test-Path -LiteralPath $dest) { throw "Do not overwrite: $dest" }
$waitStarted=[DateTime]::UtcNow
while ((Get-CimInstance Win32_OperatingSystem).FreeVirtualMemory -lt 5242880) {
  if (([DateTime]::UtcNow-$waitStarted).TotalMinutes -gt 15) {
    throw "Insufficient commit headroom for $CaseName; no solver launched or result overwritten."
  }
  Start-Sleep -Seconds 10
}
New-Item -ItemType Directory -Path $dest | Out-Null
Copy-Item -LiteralPath $input -Destination (Join-Path $dest 'input.nml')
$env:TPS_PARDISO_STORAGE='IC'
$env:MKL_NUM_THREADS=[string]$Threads
$env:OMP_NUM_THREADS=[string]$Threads
$env:MKL_DYNAMIC='FALSE'
$env:KMP_AFFINITY='granularity=fine,compact,1,0'
$env:PATH='C:\Program Files (x86)\Intel\oneAPI\2024.2\bin;C:\Program Files (x86)\Intel\oneAPI\compiler\2024.2\bin;'+$env:PATH
$started=[DateTime]::UtcNow
$provenance=@{case=$CaseName;input_name=$InputName;threads=$Threads;storage='IC';
  started_utc=$started.ToString('o');executable_sha256=(Get-FileHash -LiteralPath $exe).Hash;
  input_sha256=(Get-FileHash -LiteralPath $input).Hash;MKL_DYNAMIC=$env:MKL_DYNAMIC;
  KMP_AFFINITY=$env:KMP_AFFINITY}
$provenance | ConvertTo-Json | Set-Content -LiteralPath (Join-Path $dest 'provenance.json') -Encoding utf8
$p=Start-Process -FilePath $exe -ArgumentList @('run-f1',"`"$input`"","`"$dest`"") `
  -WorkingDirectory $root -WindowStyle Hidden -PassThru `
  -RedirectStandardOutput (Join-Path $dest 'stdout.txt') -RedirectStandardError (Join-Path $dest 'stderr.txt')
while (-not $p.WaitForExit(15000)) {
  $p.Refresh()
  @{case=$CaseName;status='RUNNING';pid=$p.Id;elapsed_seconds=([DateTime]::UtcNow-$started).TotalSeconds;
    working_set_bytes=$p.WorkingSet64;cpu_seconds=$p.TotalProcessorTime.TotalSeconds} |
    ConvertTo-Json | Set-Content -LiteralPath (Join-Path $dest 'execution.json') -Encoding utf8
}
$p.WaitForExit()
$accepted=$p.ExitCode -eq 0 -and (Test-Path -LiteralPath (Join-Path $dest 'metrics.txt')) -and
  ((Get-Content -LiteralPath (Join-Path $dest 'metrics.txt')) -contains 'case_acceptance=PASS_NOT_MESH_ACCEPTANCE')
@{case=$CaseName;status=if($accepted){'PASS_NOT_MESH_ACCEPTANCE'}else{'FAILED'};pid=$p.Id;
  exit_code=$p.ExitCode;elapsed_seconds=([DateTime]::UtcNow-$started).TotalSeconds} |
  ConvertTo-Json | Set-Content -LiteralPath (Join-Path $dest 'execution.json') -Encoding utf8
if (-not $accepted) { throw "Fast case failed: $CaseName" }
Write-Output "PASS_FAST_CASE=$CaseName"
