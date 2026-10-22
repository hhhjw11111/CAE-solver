param()

$ErrorActionPreference = 'Stop'
$solverRoot = [System.IO.Path]::GetFullPath((Split-Path -Parent $PSScriptRoot))
$culture = [System.Globalization.CultureInfo]::InvariantCulture
$length = 600.0e-9
$thickness = 20.0e-9
$holeRadius = 90.0e-9
$radialGradingPower = 2.0

function Format-Real([double]$value) {
  if ([math]::Abs($value) -lt 1.0e-30) { $value = 0.0 }
  return $value.ToString('E16',$culture)
}

function Add-Node($nodes,$nodeMap,[string]$key,[double]$x,[double]$y,[double]$z) {
  if ($nodeMap.ContainsKey($key)) { return [int]$nodeMap[$key] }
  $identifier = $nodes.Count + 1
  $nodes.Add([pscustomobject]@{ id=$identifier; x=$x; y=$y; z=$z })
  $nodeMap[$key] = $identifier
  return $identifier
}

function Write-NativeMesh([string]$path,$nodes,$elements) {
  $directory = [System.IO.Path]::GetDirectoryName($path)
  New-Item -ItemType Directory -Force -Path $directory | Out-Null
  $lines = [System.Collections.Generic.List[string]]::new()
  $lines.Add(('{0} {1}' -f $nodes.Count,$elements.Count))
  foreach ($node in $nodes) {
    $lines.Add(('{0} {1} {2} {3}' -f $node.id,(Format-Real $node.x),(Format-Real $node.y),(Format-Real $node.z)))
  }
  for ($element=0; $element -lt $elements.Count; $element++) {
    $lines.Add(('{0} 1 {1}' -f ($element+1),($elements[$element] -join ' ')))
  }
  [System.IO.File]::WriteAllLines($path,$lines,[System.Text.Encoding]::ASCII)
}

function New-GradedHoleMesh([string]$path,[int]$nTheta,[int]$nRadial) {
  if (($nTheta % 4) -ne 0) { throw 'nTheta must be divisible by four.' }
  $nZ=1
  $center=0.5*$length
  $nodes=[System.Collections.Generic.List[object]]::new()
  $elements=[System.Collections.Generic.List[object]]::new()
  $nodeMap=@{}
  $patterns=@(
    @(0,0,0),@(2,0,0),@(2,2,0),@(0,2,0),
    @(0,0,2),@(2,0,2),@(2,2,2),@(0,2,2),
    @(1,0,0),@(2,1,0),@(1,2,0),@(0,1,0),
    @(1,0,2),@(2,1,2),@(1,2,2),@(0,1,2),
    @(0,0,1),@(2,0,1),@(2,2,1),@(0,2,1))
  for($ez=0;$ez -lt $nZ;$ez++) {
    for($et=0;$et -lt $nTheta;$et++) {
      for($er=0;$er -lt $nRadial;$er++) {
        $connectivity=[System.Collections.Generic.List[int]]::new()
        foreach($pattern in $patterns) {
          $ir=2*$er+$pattern[0]
          $it=(2*$et+$pattern[1])%(2*$nTheta)
          $iz=2*$ez+$pattern[2]
          $key='{0}:{1}:{2}' -f $ir,$it,$iz
          $theta=2.0*[math]::PI*[double]$it/[double](2*$nTheta)
          $cosine=[math]::Cos($theta); $sine=[math]::Sin($theta)
          $outerRadius=(0.5*$length)/[math]::Max([math]::Abs($cosine),[math]::Abs($sine))
          $rhoLogical=[double]$ir/[double](2*$nRadial)
          if ($pattern[0] -eq 1) {
            # A quadratic HEX20 radial midside must be the geometric midpoint
            # of its two corner nodes. Applying the grading law directly to
            # this node folds the first curved layer and can invert J.
            $rhoInner=[double]($ir-1)/[double](2*$nRadial)
            $rhoOuter=[double]($ir+1)/[double](2*$nRadial)
            $radiusInner=$holeRadius + [math]::Pow($rhoInner,$radialGradingPower)*($outerRadius-$holeRadius)
            $radiusOuter=$holeRadius + [math]::Pow($rhoOuter,$radialGradingPower)*($outerRadius-$holeRadius)
            $radius=0.5*($radiusInner+$radiusOuter)
          } else {
            $rhoPhysical=[math]::Pow($rhoLogical,$radialGradingPower)
            $radius=$holeRadius+$rhoPhysical*($outerRadius-$holeRadius)
          }
          $x=$center+$radius*$cosine; $y=$center+$radius*$sine
          if([math]::Abs($x) -lt 1.0e-21){$x=0.0}
          if([math]::Abs($x-$length) -lt 1.0e-21){$x=$length}
          if([math]::Abs($y) -lt 1.0e-21){$y=0.0}
          if([math]::Abs($y-$length) -lt 1.0e-21){$y=$length}
          $z=-0.5*$thickness+$thickness*[double]$iz/[double](2*$nZ)
          $connectivity.Add((Add-Node $nodes $nodeMap $key $x $y $z))
        }
        $elements.Add($connectivity.ToArray())
      }
    }
  }
  Write-NativeMesh $path $nodes $elements
}

$levels=@(
  [pscustomobject]@{Name='L0'; Theta=16; Radial=2},
  [pscustomobject]@{Name='L1'; Theta=24; Radial=3},
  [pscustomobject]@{Name='L2'; Theta=32; Radial=4})
$outputs=[System.Collections.Generic.List[string]]::new()
foreach($level in $levels) {
  $path=Join-Path $solverRoot ('cases\validation\G10_mesh_convergence\{0}\mesh.dat' -f $level.Name)
  New-GradedHoleMesh $path $level.Theta $level.Radial
  $outputs.Add($path)
}
Write-Output ('G10_MESHES_OK=' + ($outputs -join ';'))
