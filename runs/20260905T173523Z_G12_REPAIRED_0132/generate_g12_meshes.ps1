param()

$ErrorActionPreference = 'Stop'
$solverRoot = [System.IO.Path]::GetFullPath((Split-Path -Parent $PSScriptRoot))
$culture = [System.Globalization.CultureInfo]::InvariantCulture
$length = 600.0e-9
$thickness = 20.0e-9
$radialGradingPower = 2.0

function Format-Real([double]$value) {
  if ([math]::Abs($value) -lt 1.0e-30) { $value=0.0 }
  return $value.ToString('E16',$culture)
}
function Add-Node($nodes,$nodeMap,[string]$key,[double]$x,[double]$y,[double]$z) {
  if ($nodeMap.ContainsKey($key)) { return [int]$nodeMap[$key] }
  $id=$nodes.Count+1
  $nodes.Add([pscustomobject]@{id=$id;x=$x;y=$y;z=$z})
  $nodeMap[$key]=$id
  return $id
}
function New-HoleMesh([string]$path,[double]$holeRadius) {
  $nTheta=32;$nRadial=4;$nZ=1;$center=0.5*$length
  $nodes=[System.Collections.Generic.List[object]]::new()
  $elements=[System.Collections.Generic.List[object]]::new();$nodeMap=@{}
  $patterns=@(
    @(0,0,0),@(2,0,0),@(2,2,0),@(0,2,0),@(0,0,2),@(2,0,2),@(2,2,2),@(0,2,2),
    @(1,0,0),@(2,1,0),@(1,2,0),@(0,1,0),@(1,0,2),@(2,1,2),@(1,2,2),@(0,1,2),
    @(0,0,1),@(2,0,1),@(2,2,1),@(0,2,1))
  for($ez=0;$ez -lt $nZ;$ez++) { for($et=0;$et -lt $nTheta;$et++) { for($er=0;$er -lt $nRadial;$er++) {
    $connectivity=[System.Collections.Generic.List[int]]::new()
    foreach($pattern in $patterns) {
      $ir=2*$er+$pattern[0];$it=(2*$et+$pattern[1])%(2*$nTheta);$iz=2*$ez+$pattern[2]
      $key='{0}:{1}:{2}' -f $ir,$it,$iz
      $theta=2.0*[math]::PI*[double]$it/[double](2*$nTheta)
      $c=[math]::Cos($theta);$s=[math]::Sin($theta)
      $outerRadius=(0.5*$length)/[math]::Max([math]::Abs($c),[math]::Abs($s))
      $rhoLogical=[double]$ir/[double](2*$nRadial)
      if ($pattern[0] -eq 1) {
        $rhoInner=[double]($ir-1)/[double](2*$nRadial)
        $rhoOuter=[double]($ir+1)/[double](2*$nRadial)
        $radiusInner=$holeRadius + [math]::Pow($rhoInner,$radialGradingPower)*($outerRadius-$holeRadius)
        $radiusOuter=$holeRadius + [math]::Pow($rhoOuter,$radialGradingPower)*($outerRadius-$holeRadius)
        $radius=0.5*($radiusInner+$radiusOuter)
      } else {
        $rhoPhysical=[math]::Pow($rhoLogical,$radialGradingPower)
        $radius=$holeRadius+$rhoPhysical*($outerRadius-$holeRadius)
      }
      $x=$center+$radius*$c;$y=$center+$radius*$s
      if([math]::Abs($x) -lt 1.0e-21){$x=0.0};if([math]::Abs($x-$length) -lt 1.0e-21){$x=$length}
      if([math]::Abs($y) -lt 1.0e-21){$y=0.0};if([math]::Abs($y-$length) -lt 1.0e-21){$y=$length}
      $z=-0.5*$thickness+$thickness*[double]$iz/[double](2*$nZ)
      $connectivity.Add((Add-Node $nodes $nodeMap $key $x $y $z))
    }
    $elements.Add($connectivity.ToArray())
  }}}
  $directory=[System.IO.Path]::GetDirectoryName($path);New-Item -ItemType Directory -Force -Path $directory | Out-Null
  $lines=[System.Collections.Generic.List[string]]::new();$lines.Add(('{0} {1}' -f $nodes.Count,$elements.Count))
  foreach($node in $nodes){$lines.Add(('{0} {1} {2} {3}' -f $node.id,(Format-Real $node.x),(Format-Real $node.y),(Format-Real $node.z)))}
  for($element=0;$element -lt $elements.Count;$element++){$lines.Add(('{0} 1 {1}' -f ($element+1),($elements[$element] -join ' ')))}
  [System.IO.File]::WriteAllLines($path,$lines,[System.Text.Encoding]::ASCII)
}

foreach($spec in @(
  [pscustomobject]@{Name='R060';Radius=60.0e-9},
  [pscustomobject]@{Name='R090';Radius=90.0e-9},
  [pscustomobject]@{Name='R120';Radius=120.0e-9})) {
  New-HoleMesh (Join-Path $solverRoot ('cases\validation\G12_radius_matrix\{0}\mesh.dat' -f $spec.Name)) $spec.Radius
}
Write-Output 'G12_MESHES_OK=R060;R090;R120'
