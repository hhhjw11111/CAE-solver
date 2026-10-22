"""Recompute unchanged 6/8V gates and correct trace-metric power diagnostics.

Only completed native solves are reviewed. A conditional hole solve never enters
the numerical release decision. The historical directional audit's physical
integrals are cross-checked; its stabilization decomposition is superseded here.
"""
from pathlib import Path
import argparse
import json
import subprocess
import sys
import time
import traceback
import run_je_reconstruction_mesh_20261018 as r
import run_je_controlled_trace_batch_20261021 as batch
from audit_hv_mesh import audit as mesh_audit
from audit_hv_stabilization import audit as supg_audit
from audit_hv_field_shape import audit as local_audit

ROOT = batch.ROOT


def all_points():
    return [(i,v,1.) for i in (0,1,2) for v in (6,8)]+[(2,v,s) for s in (.5,2.) for v in (6,8)]


def power_view(level, voltage, scale):
    point = batch.point(level, voltage, scale)
    target = ROOT/'trace_power'/point.name
    if not (target/'report.json').exists():
        subprocess.run([sys.executable,'-X','utf8',str(r.ROOT/'scripts/audit_je_trace_power_20261021.py'),
                        '--run',str(point),'--output',str(target)],check=True)
    provenance = r.read(target/'provenance.json')
    for path, digest in provenance['source_hashes'].items():
        if r.sha(Path(path)) != digest:
            raise RuntimeError('Corrected power evidence changed: '+path)
    report = r.read(target/'report.json')
    original = r.read(point/'physical_power.json')
    for a,b in [('raw_half_JE_W','raw_half_JE_W'),
                ('raw_half_physical_exterior_power_W','raw_half_physical_exterior_power_W'),
                ('physical_exterior_gap_relative','physical_exterior_gap_relative')]:
        if abs(report[a]-original[b]) > 1e-12*max(abs(report[a]),1.):
            raise RuntimeError('Old physical integral and corrected metric integral disagree')
    return dict(voltage_V=voltage,counts=batch.COUNTS[level],supg_scale=scale,
                report=str(target/'report.json'),**{k:report[k] for k in (
                'physical_exterior_gap_relative','reconstructed_full_JE_W',
                'assembled_order3_free_residual_max_A','independent_vs_native_relative')})


def finalize():
    batch.validate_prerequisites()
    manifest = r.read(batch.OUT/'runner_manifest.json')
    resume = r.read(batch.OUT/'runner_manifest_v2.json')
    for path,digest in {**manifest['script_hashes'],**manifest['prerequisite_hashes'],**resume['script_hashes']}.items():
        if r.sha(Path(path)) != digest:
            raise RuntimeError('Runner dependency changed: '+path)
    baseline = r.read(ROOT/'baseline_contract.json')
    for path,digest in baseline['preserved_evidence_hashes'].items():
        if r.sha(Path(path)) != digest:
            raise RuntimeError('Preserved original evidence changed: '+path)
    mesh = [mesh_audit([batch.point(i,v) for i in (0,1,2)],[v]) for v in (6,8)]
    supg = [supg_audit([batch.point(2,v,s) for s in (.5,1.,2.)]) for v in (6,8)]
    local = []
    for i,v,s in all_points():
        point = batch.point(i,v,s)
        batch.checked_half(point)
        if i == 0 or i == 1:
            continue
        label = 'mesh_final_pair' if s == 1. else 'supg_same_mesh'
        source = batch.point(i-1,v) if s == 1. else batch.point(2,v)
        saved = r.read(point/(label+'_local_comparison_v2.json'))
        mapping_files = [Path(p) for p in saved['hashes'] if Path(p).name=='transfer_report.json']
        if len(mapping_files) != 1:
            raise RuntimeError('Ambiguous local comparison mapping')
        report = local_audit(mapping_files[0].parent, source, point)
        local.append(dict(voltage_V=v, supg_scale=s, passed=report['passed'],comparisons=report['comparisons']))
    powers = [power_view(*args) for args in all_points()]
    passed = (all(m['all_mesh_pass'] for m in mesh) and all(s['covered_metrics_passed'] for s in supg)
              and len(local)==6 and all(x['passed'] for x in local))
    out = ROOT/'final_review'
    out.mkdir(exist_ok=False)
    result = dict(scope='CONTROLLED_TRACE_NUMERICAL_REVIEW',original_thresholds_unchanged=True,
        original_core_and_model_unchanged=True,old_data_hashes_rechecked=True,
        numerical_setting='existing supg_directional_tau=false; original scales0.5/1/2',
        mesh=mesh,supg=supg,local=local,power=powers,
        all_original_endpoint_numerical_checks_passed=passed,
        production_allowed_by_endpoint_prerequisite=passed,production_started=False,
        power_gap_threshold=None,conditional_diagnostics_used_for_release=False,
        power_audit_scope='Physical JE/flux in run power_audit valid; directional S/reaction decomposition superseded by trace_power/',
        limitations=['Only r=3um symmetric NT endpoints; intermediate voltages, LQ and other radii need their own evidence',
                     'Current-work heat model is R=0 diagnostic, not thermodynamic closure or experimental calibration'])
    r.save(out/'report.json',result)
    lines = ['6/8 V受控数值改进复核','原模型、原热源、原门限、旧数据保持不变。',
             '完整度量SUPG的原生解析、低Pe、一致切线及既有checked-block复算通过。']
    for v,m,s in zip((6,8),mesh,supg):
        changes=m['voltages'][0]['changes']
        failed={k:x['L1_L2'] for k,x in changes.items() if x['L1_L2']>=x['limit']}
        fails={k:x['relative_span'] for k,x in s['comparisons'].items() if not x['passed']}
        lines.append(f'{v} V：联合网格={m["all_mesh_pass"]}，失败项={failed}；SUPG={s["covered_metrics_passed"]}，失败项={fails}。')
    lines.append('局部最终一对和SUPG比较：'+str(all(x['passed'] for x in local)))
    lines.append('端点全部原数值门限通过：'+str(passed))
    lines.append('后续生产仍须对应模型/孔径的验证；本报告不凭趋势或功率改善替代未通过项。')
    (out/'结论.txt').write_text('\n'.join(lines)+'\n',encoding='utf8')
    r.save(out/'artifact_hashes.json',{str(p.relative_to(out)):r.sha(p) for p in out.iterdir() if p.is_file()})
    print(json.dumps(dict(final_review=str(out),passed=passed)),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--watch',action='store_true')
    args=parser.parse_args()
    if args.watch:
        while not (batch.OUT/'final_audit_v2.json').exists():
            if (batch.OUT/'stopped_v2.json').exists():
                raise RuntimeError('Native controller stopped: '+str(r.read(batch.OUT/'stopped_v2.json')))
            for point in all_points():
                path=batch.point(*point)
                if (path/'completion_v2.json').exists() and not (ROOT/'trace_power'/path.name/'report.json').exists():
                    power_view(*point)
            time.sleep(15)
    finalize()


if __name__=='__main__':
    main()
