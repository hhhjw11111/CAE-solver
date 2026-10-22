param(
    [Parameter(Mandatory=$true)][string]$InputDocx,
    [Parameter(Mandatory=$true)][string]$OutputDocx
)
# Normalize a generated native-equation DOCX with Word. Never overwrite the source.
# This does not replace an actual MathType conversion acceptance test.
$ErrorActionPreference = 'Stop'
$inputPath = (Resolve-Path -LiteralPath $InputDocx).Path
$outputPath = [IO.Path]::GetFullPath($OutputDocx)
if ($inputPath -eq $outputPath -or (Test-Path -LiteralPath $outputPath)) {
    throw 'OutputDocx must be a new file; original and historical documents are preserved.'
}
if (-not (Test-Path -LiteralPath ([IO.Path]::GetDirectoryName($outputPath)))) {
    throw 'Output directory must already exist.'
}
$wordApp = $null
$wordDoc = $null
try {
    $wordApp = New-Object -ComObject Word.Application
    $wordApp.Visible = $false
    $wordApp.DisplayAlerts = 0
    $wordApp.AutomationSecurity = 3
    $wordDoc = $wordApp.Documents.Open($inputPath, $false, $true, $false)
    $before = $wordDoc.OMaths.Count
    $oldMode = $wordDoc.CompatibilityMode
    $wordDoc.SetCompatibilityMode(15)
    $wordDoc.SaveAs2($outputPath, 16)
    $after = $wordDoc.OMaths.Count
    if ($after -ne $before) { throw "Native equation count changed: $before -> $after" }
    [pscustomobject]@{
        Input = $inputPath; Output = $outputPath
        PreviousCompatibilityMode = $oldMode
        CompatibilityMode = $wordDoc.CompatibilityMode
        NativeEquations = $after
        MathTypeConversionTest = 'Required separately on a copy'
    } | ConvertTo-Json
} finally {
    if ($null -ne $wordDoc) {
        [object]$noSaveDocument = 0
        $wordDoc.Close([ref]$noSaveDocument)
        [void][Runtime.InteropServices.Marshal]::FinalReleaseComObject($wordDoc)
    }
    if ($null -ne $wordApp) {
        [object]$noSaveApplication = 0
        $wordApp.Quit([ref]$noSaveApplication)
        [void][Runtime.InteropServices.Marshal]::FinalReleaseComObject($wordApp)
    }
    [GC]::Collect()
    [GC]::WaitForPendingFinalizers()
}
