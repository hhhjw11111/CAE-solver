# Resource housekeeping only: completed, new-study PARDISO factors are disposable.
# Scientific evidence, OOC configuration, and every pre-existing run are retained.
$ErrorActionPreference = 'Stop'
$taskStudyRoot = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot '../runs/JE_FINE_CROSS_20261016'))
$taskStudyRoot = (Resolve-Path -LiteralPath $taskStudyRoot).Path
$taskContract = Get-Content -LiteralPath (Join-Path $taskStudyRoot 'contract.json') -Raw | ConvertFrom-Json
$taskLockPath = Join-Path $taskStudyRoot 'collector_started.json'
if (Test-Path -LiteralPath $taskLockPath) { throw 'Collector was already started; inspect its status before restarting.' }
@{pid=$PID;scope='Only completed new-study regenerable PARDISO factors';started=(Get-Date -Format o)} | ConvertTo-Json | Set-Content -LiteralPath $taskLockPath -Encoding utf8
do {
    foreach ($taskJob in $taskContract.jobs) {
        $taskRun = [IO.Path]::GetFullPath((Join-Path $taskStudyRoot $taskJob.name))
        if (-not $taskRun.StartsWith($taskStudyRoot + '\', [StringComparison]::OrdinalIgnoreCase)) { throw 'Run escaped study root' }
        $taskCleanupRecord = Join-Path $taskRun 'temporary_factor_cleanup.json'
        if (-not (Test-Path -LiteralPath (Join-Path $taskRun 'summary.json'))) { continue }
        if (-not (Test-Path -LiteralPath (Join-Path $taskRun 'execution.json'))) { continue }
        if (Test-Path -LiteralPath $taskCleanupRecord) { continue }
        # execution.json is written only after Popen.wait() / kill()+wait() returned.
        $taskScratch = (Resolve-Path -LiteralPath (Join-Path $taskRun 'pardiso_ooc')).Path
        if ($taskScratch -ne [IO.Path]::GetFullPath((Join-Path $taskRun 'pardiso_ooc'))) { throw 'Unexpected scratch resolution' }
        if (-not $taskScratch.StartsWith($taskStudyRoot + '\', [StringComparison]::OrdinalIgnoreCase)) { throw 'Scratch escaped new study' }
        $taskFiles = @(Get-ChildItem -LiteralPath $taskScratch -File | Where-Object { $_.Name -match '^factors[0-9]*\.(ind|jal|jau|lnz|lpr|lup|sin|sle|unz|upr|uup)$' })
        $taskInventory = @($taskFiles | ForEach-Object { @{name=$_.Name;bytes=$_.Length} })
        foreach ($taskFile in $taskFiles) {
            if ($taskFile.DirectoryName -ne $taskScratch) { throw 'Unexpected factor location' }
            Remove-Item -LiteralPath $taskFile.FullName
        }
        @{scope='Regenerable temporary factors only; all scientific evidence and run_ooc.cfg retained';files=$taskInventory;completed=(Get-Date -Format o)} | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath $taskCleanupRecord -Encoding utf8
        Write-Output ('TEMPORARY_FACTORS_RELEASED ' + $taskJob.name)
    }
    $taskDone = Test-Path -LiteralPath (Join-Path $taskStudyRoot 'execution_complete.json')
    if (-not $taskDone) { Start-Sleep -Seconds 10 }
} while (-not $taskDone)
$taskAnalysis = Join-Path $taskStudyRoot 'analysis'
if (Test-Path -LiteralPath $taskAnalysis) { throw 'Analysis destination already exists; refusing overwrite' }
python (Join-Path $PSScriptRoot 'summarize_je_fine_cross_20261016.py') --output $taskAnalysis
if ($LASTEXITCODE -ne 0) { throw 'Evidence collection failed; preserve diagnostics' }
Write-Output 'FINAL_REPORT_READY'
