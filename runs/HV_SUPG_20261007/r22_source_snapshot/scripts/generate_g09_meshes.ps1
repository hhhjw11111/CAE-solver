param()

$ErrorActionPreference = 'Stop'
$solverRoot = [System.IO.Path]::GetFullPath((Split-Path -Parent $PSScriptRoot))
$culture = [System.Globalization.CultureInfo]::InvariantCulture
$length = 600.0e-9
$thickness = 20.0e-9
$holeRadius = 90.0e-9

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

function New-PlainMesh([string]$path) {
  $nx=4; $ny=4; $nz=1
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
          $ix=2*$ex+$pattern[0]; $iy=2*$ey+$pattern[1]; $iz=2*$ez+$pattern[2]
          $key='{0}:{1}:{2}' -f $ix,$iy,$iz
          $x=$length*[double]$ix/[double](2*$nx)
          $y=$length*[double]$iy/[double](2*$ny)
          $z=-0.5*$thickness+$thickness*[double]$iz/[double](2*$nz)
          $connectivity.Add((Add-Node $nodes $nodeMap $key $x $y $z))
        }
        $elements.Add($connectivity.ToArray())
      }
    }
  }
  Write-NativeMesh $path $nodes $elements
}

function New-HoleMesh([string]$path) {
  $nTheta=16; $nRadial=2; $nZ=1
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
          $rho=[double]$ir/[double](2*$nRadial)
          $radius=$holeRadius+$rho*($outerRadius-$holeRadius)
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

$plainPath=Join-Path $solverRoot 'cases\validation\G09_center_hole\H00_plain\mesh.dat'
$holePath=Join-Path $solverRoot 'cases\validation\G09_center_hole\H02_r090\mesh.dat'
New-PlainMesh $plainPath
New-HoleMesh $holePath
Write-Output "G09_MESHES_OK=$plainPath;$holePath"
