param([string]$Only = '')
$ErrorActionPreference = 'Stop'
$sprint = Split-Path $PSScriptRoot -Parent
$pandoc = 'C:\Program Files\Pandoc\pandoc.exe'
$xelatex = 'C:\Program Files\MiKTeX\miktex\bin\x64\xelatex.exe'
$documents = @('reports/study_report_final','manuscript/manuscript_working_draft')
if ($Only) { $documents = @($documents | Where-Object { $_.StartsWith($Only + '/') }) }
Push-Location $sprint
try {
  foreach ($document in $documents) {
    $stem = Split-Path $document -Leaf
    $log = Join-Path $sprint ('audit/' + $stem + '_render.log')
    $argsForPandoc = @("$document.md", '-o', "$document.pdf", '--resource-path=.;reports;manuscript', "--pdf-engine=$xelatex", '--pdf-engine-opt=--disable-installer', '--pdf-engine-opt=-no-shell-escape', '-V', 'mainfont=Cambria', '-V', 'CJKmainfont=Microsoft YaHei', '-V', 'sansfont=Microsoft YaHei', '-V', 'monofont=Consolas', '-V', 'mathfont=Cambria Math', '-V', 'geometry:margin=22mm', '-V', 'papersize=a4', '-V', 'fontsize=11pt', '-V', 'linestretch=1.2', '-V', 'colorlinks=true')
    if ($document.StartsWith('manuscript/')) { $argsForPandoc += @('--from=markdown-implicit_figures', '--lua-filter=scripts/paper_layout.lua', '--include-in-header=scripts/paper_layout_header.tex') }
    if ($document.StartsWith('reports/')) { $argsForPandoc += @('--from=markdown-implicit_figures', '--lua-filter=scripts/report_layout.lua', '--include-in-header=scripts/report_layout_header.tex') }
    & $pandoc @argsForPandoc 2>&1 | Tee-Object -FilePath $log
    if ($LASTEXITCODE -ne 0) { throw ('PDF_RENDER_FAILED ' + $document + ' exit=' + $LASTEXITCODE) }
    Write-Output ('PDF_RENDERED ' + $document + '.pdf')
  }
} finally { Pop-Location }
