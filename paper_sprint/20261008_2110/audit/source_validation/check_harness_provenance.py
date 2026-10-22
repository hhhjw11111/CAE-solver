from pathlib import Path
import hashlib,json,subprocess
out=Path(__file__).absolute().parent;root=out.parents[3]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
arch=root/'runs/HV_SUPG_20261007/log_weighted_production_bar_r22'
m=json.loads((arch/'manifest.json').read_text(encoding='utf-8-sig'))
checks=[]
for p,want in [(arch/'log_weighted_probe.exe',m['harness_executable_sha256']),
               (arch/'driver_source.f90',m['source_sha256']),
               (arch/'module_source.f90',m['prototype_module_sha256']),
               (Path(m['production_build'])/'TPSolverV2.exe',m['production_executable_sha256'])]:
    checks.append({'path':str(p),'sha256':sha(p),'expected_sha256':want,'matches':sha(p)==want})
for obj in m['objects']:
    p=Path(m['production_build'])/obj['name'];checks.append({'path':str(p),'sha256':sha(p),'expected_sha256':obj['sha256'],'matches':sha(p)==obj['sha256']})
snapshot=root/'runs/HV_SUPG_20261007/r22_source_snapshot'
comparisons=[]
for f in ['src/physics/material_point_m.f90','src/physics/transport_m.f90','src/physics/thermo_piezo_m.f90','src/physics/zno_candidate_m.f90','src/fem/element_kernel_m.f90']:
    p=root/f;old=snapshot/f
    comparisons.append({'path':f,'current_sha256':sha(p),'r22_sha256':sha(old),'byte_identical':p.read_bytes()==old.read_bytes()})
    if p.read_bytes()!=old.read_bytes():
        d=subprocess.run(['git','diff','--no-index','--',str(old),str(p)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        (out/'r22_current_element_diff.patch').write_bytes(d.stdout)
result={'classification':'R22_PRODUCTION_KERNEL_HARNESS_PROVENANCE_NOT_CURRENT_R25_LOGWEIGHTED_REBUILD',
        'checked_items':checks,'all_manifest_items_match':all(x['matches'] for x in checks),
        'current_vs_r22':comparisons,'manifest_path':str(arch/'manifest.json'),'manifest_sha256':sha(arch/'manifest.json')}
(out/'harness_provenance.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'checked_items':len(checks),'all_match':result['all_manifest_items_match'],'current_vs_r22':comparisons},indent=2))
