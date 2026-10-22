param()

$ErrorActionPreference = 'Stop'
$solverRoot = [System.IO.Path]::GetFullPath((Split-Path -Parent $PSScriptRoot))
$culture = [System.Globalization.CultureInfo]::InvariantCulture
$length = 600.0e-9
$thickness = 20.0e-9

function Format-Real([double]$value) {
  if ([math]::Abs($value) -lt 1.0e-30) { $value = 0.0 }
  return $value.ToString('E16',$culture)
}

function Add-Node($nodes,$nodeMap,[string]$key,[double]$x,[double]$y,[double]$z) {
  if ($nodeMap.ContainsKey($key)) { return [int]$nodeMap[$key] }
  $identifier=$nodes.Count+1
  $nodes.Add([pscustomobject]@{id=$identifier;x=$x;y=$y;z=$z})
  $nodeMap[$key]=$identifier
  return $identifier
}

function Mapped-Coordinate([int]$index,[int]$elementCount) {
  $q=([double]$index-[double]$elementCount)/[double]$elementCount
  $value=0.5*$length+0.5*$length*[math]::Tan(0.25*[math]::PI*$q)
  if([math]::Abs($value) -lt 1.0e-21){$value=0.0}
  if([math]::Abs($value-$length) -lt 1.0e-21){$value=$length}
  return $value
}

function New-MatchedPlainMesh([string]$path) {
  $nx=8;$ny=8;$nz=1
  $nodes=[System.Collections.Generic.List[object]]::new()
  $elements=[System.Collections.Generic.List[object]]::new()
  $nodeMap=@{}
  $patterns=@(
    @(0,0,0),@(2,0,0),@(2,2,0),@(0,2,0),
    @(0,0,2),@(2,0,2),@(2,2,2),@(0,2,2),
    @(1,0,0),@(2,1,0),@(1,2,0),@(0,1,0),
    @(1,0,2),@(2,1,2),@(1,2,2),@(0,1,2),
    @(0,0,1),@(2,0,1),@(2,2,1),@(0,2,1))
  for($ez=0;$ez -lt $nz;$ez++) {
    for($ey=0;$ey -lt $ny;$ey++) {
      for($ex=0;$ex -lt $nx;$ex++) {
        $connectivity=[System.Collections.Generic.List[int]]::new()
        foreach($pattern in $patterns) {
          $ix=2*$ex+$pattern[0];$iy=2*$ey+$pattern[1];$iz=2*$ez+$pattern[2]
          $key='{0}:{1}:{2}' -f $ix,$iy,$iz
          $x=Mapped-Coordinate $ix $nx
          $y=Mapped-Coordinate $iy $ny
          $z=-0.5*$thickness+$thickness*[double]$iz/[double](2*$nz)
          $connectivity.Add((Add-Node $nodes $nodeMap $key $x $y $z))
        }
        $elements.Add($connectivity.ToArray())
      }
    }
  }
  $directory=[System.IO.Path]::GetDirectoryName($path)
  New-Item -ItemType Directory -Force -Path $directory | Out-Null
  $lines=[System.Collections.Generic.List[string]]::new()
  $lines.Add(('{0} {1}' -f $nodes.Count,$elements.Count))
  foreach($node in $nodes) {
    $lines.Add(('{0} {1} {2} {3}' -f $node.id,(Format-Real $node.x),(Format-Real $node.y),(Format-Real $node.z)))
  }
  for($element=0;$element -lt $elements.Count;$element++) {
    $lines.Add(('{0} 1 {1}' -f ($element+1),($elements[$element] -join ' ')))
  }
  [System.IO.File]::WriteAllLines($path,$lines,[System.Text.Encoding]::ASCII)
}

$output=Join-Path $solverRoot 'cases\validation\G11_matched_control\H00_L2\mesh.dat'
New-MatchedPlainMesh $output
Write-Output "G11_MESH_OK=$output"
