"""Prepare, never solve, an approved NT/JE paper condition.

Default operation is a read-only plan. --prepare writes input.nml and an
INPUT_ONLY_NOT_SOLVED manifest in a fresh directory. No native launch exists.
The pure make_input function performs no file-system access.
"""
from pathlib import Path
import argparse
import ctypes
import hashlib
import json
import math
import re
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = ROOT/'runs/JE_6_8_HALF_V_20261013/input_V060.nml'
TEMPLATE_SHA256 = '3e9c15c38ea55c488c95594b36ba21f73603936d5453e035f996ce59f68e1f02'
TEMPLATE_TEXT_SHA256 = 'd23227689905841528fb411931cd52d5d7b5e96b5f1996a614ffb55b1cfc12db'
JE = 'E_DOT_TOTAL_CURRENT_DIAGNOSTIC_NOT_THERMODYNAMIC_CLOSURE'
DRIFT = 'DRIFT_CONDUCTIVITY_NOT_TOTAL_CURRENT_WORK'
CONDITIONS = {(3,v) for v in (6.,6.5,7.,7.5,8.)} | {(r,v) for r in (2,4) for v in (6.,8.)}
FIELDS = ('node','x_m','y_m','z_m','ux_m','uy_m','uz_m','phi_V','eta_n','eta_p','theta_K','n_m3','p_m3')


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream,'sha256').hexdigest()


def save(path, value):
    with Path(path).open('x',encoding='utf-8') as stream:
        json.dump(value,stream,indent=2,ensure_ascii=False,allow_nan=False)
        stream.write('\n')


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def _replace(text, key, value):
    pattern = r'(?<![A-Za-z0-9_])'+re.escape(key)+r"\s*=\s*(?:'[^']*'|[^,\n/]+)"
    result,count = re.subn(pattern,lambda _:key+'='+value,text)
    if count != 1:
        raise ValueError('Expected exactly one '+key)
    return result


def _quoted_path(value):
    value=str(value).replace('\\','/')
    if not value or len(value)>511 or any(c in value for c in "'\r\n"):
        raise ValueError('Invalid native path')
    value.encode('ascii')
    return "'"+value+"'"


def make_input(template_text, *, radius_um, voltage, mesh_path, seed_nodes_path, seed_voltage):
    """Pure conversion of the frozen V060 template; all other text is preserved."""
    text=template_text.replace('\r\n','\n')
    if hashlib.sha256(text.encode('utf-8')).hexdigest()!=TEMPLATE_TEXT_SHA256:
        raise ValueError('Frozen V060 template changed')
    if (radius_um,voltage) not in CONDITIONS:
        raise ValueError('Condition is outside the approved nine NT/JE cases')
    if not math.isfinite(seed_voltage) or not 0<seed_voltage<=voltage:
        raise ValueError('Require 0 < seed_voltage <= target_voltage')
    radius=float(f'{radius_um:g}e-6')
    replacements=dict(case_id=f"'JE_PAPER_R{int(radius_um):02d}_NT_V{round(voltage*10):03d}'",
        mesh_file=_quoted_path(mesh_path),hole_radius=f'{radius:.17g}',target_voltage=f'{-voltage:.17g}',
        target_voltage_factor='1.0',target_joule_factor='1.0',heat_source_current_work='.true.',
        restart_nodes_file=_quoted_path(seed_nodes_path),restart_voltage_factor=f'{seed_voltage/voltage:.17g}',
        restart_joule_factor='1.0')
    for key,value in replacements.items():
        text=_replace(text,key,value)
    # This includes SUPG=1, the original field scales (potential scale 8 V),
    # continuation factors, tolerances and arclength disabled. No rescaling.
    masked_before=template_text.replace('\r\n','\n');masked_after=text
    for key in replacements:
        masked_before=_replace(masked_before,key,'UNCHANGED_COMPARISON')
        masked_after=_replace(masked_after,key,'UNCHANGED_COMPARISON')
    if masked_before!=masked_after:
        raise ValueError('Unexpected edit outside preparation keys')
    text.encode('ascii')
    return text


def native_path(path):
    path=Path(path).resolve(strict=True)
    value=path.as_posix()
    try:
        value.encode('ascii')
    except UnicodeEncodeError:
        buffer=ctypes.create_unicode_buffer(32768)
        if not ctypes.windll.kernel32.GetShortPathNameW(str(path),buffer,len(buffer)):
            raise OSError('Cannot obtain existing ASCII path: '+str(path))
        value=buffer.value.replace('\\','/')
        value.encode('ascii')
    return value


def mesh_coordinates(path):
    with Path(path).open(encoding='utf-8-sig') as stream:
        nn,ne=map(int,stream.readline().split())
        xyz=np.empty((nn,3))
        for i in range(nn):
            row=stream.readline().split()
            if len(row)!=4 or int(row[0])!=i+1:
                raise ValueError('Invalid mesh node order')
            xyz[i]=list(map(float,row[1:]))
    if not np.isfinite(xyz).all():
        raise ValueError('Nonfinite mesh coordinates')
    return xyz,ne


def validate_mesh(mesh, manifest_path, radius_um):
    manifest=read_json(manifest_path)
    radius=float(f'{radius_um:g}e-6')
    required=dict(hole_radius_m=radius,ring_width_m=.5e-6,annulus_inner_radius_m=radius,
                  annulus_outer_radius_m=radius+.5e-6)
    for key,expected in required.items():
        value=manifest.get(key)
        if not isinstance(value,(int,float)) or not math.isclose(value,expected,rel_tol=1e-13,abs_tol=0.):
            raise ValueError('Mesh manifest geometry mismatch: '+key)
    if manifest.get('domain_symmetry')!='FULL' or manifest.get('mesh_sha256')!=sha(mesh):
        raise ValueError('Mesh hash or full-domain identity mismatch')
    if manifest.get('geometry_audit',{}).get('status')!='GEOMETRY_CHECKS_PASSED_NOT_SOLUTION_ACCEPTANCE':
        raise ValueError('Require the new geometry-only audit, not a copied mesh certificate')
    xyz,ne=mesh_coordinates(mesh)
    if manifest.get('nodes')!=len(xyz) or manifest.get('elements')!=ne:
        raise ValueError('Mesh manifest counts differ')
    radial=np.linalg.norm(xyz[:,:2],axis=1)
    if abs(float(radial.min())-radius)>1e-18:
        raise ValueError('Actual mesh hole radius mismatch')
    if not np.any(abs(radial-(radius+.5e-6))<1e-18):
        raise ValueError('Missing radius-specific ring interface')
    return xyz,manifest


def validate_seed_coordinates(seed_nodes, xyz, seed_voltage=None):
    data=np.genfromtxt(seed_nodes,delimiter=',',names=True)
    data=np.atleast_1d(data)
    if data.dtype.names!=FIELDS or len(data)!=len(xyz):
        raise ValueError('Seed header or node count mismatch')
    if not np.array_equal(data['node'],np.arange(1,len(xyz)+1)):
        raise ValueError('Seed node order mismatch')
    positions=np.column_stack([data[k] for k in ('x_m','y_m','z_m')])
    if not np.array_equal(positions,xyz):
        raise ValueError('Seed coordinates must exactly match the target mesh')
    if not all(np.isfinite(data[k]).all() for k in FIELDS):
        raise ValueError('Nonfinite seed fields')
    if (np.min(300.+data['theta_K'])<=0 or np.min(data['n_m3'])<=0 or np.min(data['p_m3'])<=0):
        raise ValueError('Nonpositive seed temperature or density')
    voltage_error=None
    if seed_voltage is not None:
        left=abs(xyz[:,0]-xyz[:,0].min())<1e-18
        right=abs(xyz[:,0]-xyz[:,0].max())<1e-18
        voltage_error=max(float(np.max(abs(data['phi_V'][left]))),
                          float(np.max(abs(data['phi_V'][right]+seed_voltage))))
        if voltage_error>1e-10:
            raise ValueError('Actual seed electrode voltage differs from its declared load')
    return dict(nodes=len(data),coordinates_exactly_equal=True,maximum_coordinate_difference_m=0.,
                finite_fields=True,positive_temperature_and_density=True,electrode_voltage_error_V=voltage_error)


def metrics(path):
    values={}
    for line in Path(path).read_text(encoding='utf-8-sig').splitlines():
        if '=' not in line or ',' in line:
            continue
        key,value=line.split('=',1)
        try:value=float(value.strip().replace('D','E'))
        except ValueError:value=value.strip()
        values[key.strip()]=value
    return values


def validate_source_metrics(m, radius_um, seed_voltage, allow_drift_seed=False):
    """Source identity only; its old PASS label never certifies the new JE input."""
    if m.get('model')!='NT' or not math.isclose(m.get('hole_radius_m',math.nan),float(f'{radius_um:g}e-6'),rel_tol=1e-13):
        raise ValueError('Seed source must be NT at the same physical hole radius')
    if m.get('voltage_V')!=seed_voltage:
        raise ValueError('Declared seed voltage differs from the saved source voltage')
    source_heat=m.get('heat_source')
    if source_heat==DRIFT:
        if not allow_drift_seed or seed_voltage!=3.:
            raise ValueError('Old drift seed requires explicit --allow-drift-seed and its actual 3 V source')
    elif source_heat!=JE:
        raise ValueError('Unknown seed heat-source branch')
    return source_heat


def build_plan(*, radius_um, voltage, mesh, mesh_manifest, seed_nodes, seed_run, seed_voltage,
               transfer_report=None, seed_source_mesh=None, allow_drift_seed=False):
    """Read and validate the provenance chain; never write or launch a solver."""
    paths={k:Path(v).resolve(strict=True) for k,v in dict(mesh=mesh,mesh_manifest=mesh_manifest,
           seed_nodes=seed_nodes,seed_run=seed_run).items()}
    mesh,mesh_manifest,seed_nodes,seed_run=(paths[k] for k in ('mesh','mesh_manifest','seed_nodes','seed_run'))
    if sha(TEMPLATE)!=TEMPLATE_SHA256:
        raise ValueError('Frozen input template file identity changed')
    xyz,manifest=validate_mesh(mesh,mesh_manifest,radius_um)
    seed_check=validate_seed_coordinates(seed_nodes,xyz,seed_voltage)
    source_nodes=seed_run/'nodes.csv';source_metrics=seed_run/'metrics.txt'
    m=metrics(source_metrics)
    source_heat=validate_source_metrics(m,radius_um,seed_voltage,allow_drift_seed)
    files=[mesh,mesh_manifest,seed_nodes,source_nodes,source_metrics,TEMPLATE,Path(__file__),
           ROOT/'src/app/f1_case_runner_m.f90',ROOT/'src/nonlinear/continuation_solver_m.f90']
    mapping=None
    if sha(seed_nodes)!=sha(source_nodes):
        report_path=Path(transfer_report or seed_nodes.parent/'transfer_report.json').resolve(strict=True)
        report=read_json(report_path)
        if report.get('status')!='INITIAL_GUESS_ONLY_NOT_A_SOLVED_FIELD':
            raise ValueError('Unidentified transferred seed')
        if report.get('source_nodes_sha256')!=sha(source_nodes) or report.get('target_mesh_sha256')!=sha(mesh):
            raise ValueError('Seed transfer source/target identity mismatch')
        if report.get('boundary_extension_nodes')!=0:
            raise ValueError('Preparation requires a transfer without boundary extrapolation')
        if seed_source_mesh is None:
            sm=read_json(seed_run/'manifest.json')
            seed_source_mesh=sm['mesh']
        source_mesh=Path(seed_source_mesh).resolve(strict=True)
        if report.get('source_mesh_sha256')!=sha(source_mesh):
            raise ValueError('Seed source mesh identity mismatch')
        source_xyz,_=mesh_coordinates(source_mesh)
        validate_seed_coordinates(source_nodes,source_xyz,seed_voltage)
        source_radius=float(np.linalg.norm(source_xyz[:,:2],axis=1).min())
        if abs(source_radius-float(f'{radius_um:g}e-6'))>1e-18:
            raise ValueError('Transfer came from a different physical hole')
        files += [report_path,source_mesh]
        mapping=dict(report=str(report_path),source_mesh=str(source_mesh),boundary_extension_nodes=0)
    else:
        validate_seed_coordinates(source_nodes,xyz,seed_voltage)
    text=make_input(TEMPLATE.read_text(encoding='utf-8-sig'),radius_um=radius_um,voltage=voltage,
        mesh_path=native_path(mesh),seed_nodes_path=native_path(seed_nodes),seed_voltage=seed_voltage)
    plan=dict(status='INPUT_ONLY_NOT_SOLVED',condition=dict(radius_um=radius_um,voltage_V=voltage,model='NT'),
        model_contract='JE_SIMPLIFIED_R0_V1',heat_source=JE,mesh=str(mesh),mesh_sha256=sha(mesh),
        mesh_manifest=str(mesh_manifest),mesh_counts=manifest['counts'],seed_nodes=str(seed_nodes),
        seed_source_run=str(seed_run),seed_source_heat=source_heat,seed_source_voltage_V=seed_voltage,
        restart_voltage_factor=seed_voltage/voltage,restart_joule_factor=1.,target_voltage=-voltage,
        target_voltage_factor=1.,target_joule_factor=1.,seed_validation=seed_check,transfer=mapping,
        seed_status='INITIAL_GUESS_REQUIRES_JE_REEQUILIBRATION',
        seed_is_JE_acceptance_evidence=False,source_acceptance_label=m.get('case_acceptance'),
        drift_seed_explicitly_enabled=allow_drift_seed,
        unchanged_numerics='All template text outside explicit identity/geometry/load/restart/JE keys retained; residual/update 1e-9, full heat 1, SUPG 1, potential field scale 8.',
        restart_semantics='Native continuation unconditionally solves the start state under JE at seed voltage before progressing; seed at target voltage is also re-equilibrated.',
        requires_spatial_and_supg_acceptance_before_production=True,
        native_solver_started=False,production_acceptance=False,
        input_sha256=hashlib.sha256(text.encode('ascii')).hexdigest(),input_text=text,
        hashes={str(p):sha(p) for p in files})
    return plan


def prepare(plan, directory):
    """Write a reviewed plan only to a fresh directory. Still never runs native."""
    for path,digest in plan['hashes'].items():
        if sha(Path(path))!=digest:
            raise ValueError('Preparation source changed after plan: '+path)
    directory=Path(directory).resolve()
    directory.mkdir(parents=True,exist_ok=False)
    text=plan['input_text']
    with (directory/'input.nml').open('x',encoding='ascii',newline='\n') as stream:
        stream.write(text)
    if sha(directory/'input.nml')!=plan['input_sha256']:
        raise ValueError('Written input identity mismatch')
    manifest={k:v for k,v in plan.items() if k!='input_text'}
    manifest['prepared_input']=str(directory/'input.nml')
    save(directory/'preparation_manifest.json',manifest)
    return manifest


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--radius-um',type=int,choices=(2,3,4),required=True)
    p.add_argument('--voltage',type=float,required=True)
    p.add_argument('--mesh',type=Path,required=True)
    p.add_argument('--mesh-manifest',type=Path)
    p.add_argument('--seed-nodes',type=Path,required=True)
    p.add_argument('--seed-run',type=Path,required=True,help='Original source run with nodes.csv and metrics.txt')
    p.add_argument('--seed-voltage',type=float,required=True)
    p.add_argument('--transfer-report',type=Path)
    p.add_argument('--seed-source-mesh',type=Path,help='Required for transfers from historical runs without a mesh manifest')
    p.add_argument('--allow-drift-seed',action='store_true',help='Allow same-radius old NT 3 V only as an initial guess')
    p.add_argument('--output',type=Path,help='Fresh preparation folder; required only with --prepare')
    p.add_argument('--prepare',action='store_true',help='Write input only; omission is a read-only plan')
    a=p.parse_args()
    if a.prepare and a.output is None:p.error('--prepare requires --output')
    plan=build_plan(radius_um=a.radius_um,voltage=a.voltage,mesh=a.mesh,
        mesh_manifest=a.mesh_manifest or a.mesh.parent/'manifest.json',seed_nodes=a.seed_nodes,
        seed_run=a.seed_run,seed_voltage=a.seed_voltage,transfer_report=a.transfer_report,
        seed_source_mesh=a.seed_source_mesh,allow_drift_seed=a.allow_drift_seed)
    result=prepare(plan,a.output) if a.prepare else plan
    print(json.dumps(result,indent=2,ensure_ascii=True,allow_nan=False))


if __name__=='__main__':main()
