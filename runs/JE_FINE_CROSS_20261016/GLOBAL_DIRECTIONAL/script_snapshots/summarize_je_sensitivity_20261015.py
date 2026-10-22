"""Collect completed fixed-physics sensitivity evidence without changing gates."""
from pathlib import Path
import csv,json,re,shutil
import numpy as np
from run_hv_sweep import read_metrics,sha,classify
from audit_hv_stabilization import audit as stabilization_audit
from audit_hv_stabilization import METRICS as SUPG_METRICS,TOLERANCE as SUPG_TOLERANCE
from audit_hv_mesh import audit as mesh_audit
from audit_hv_field_shape import audit as field_audit,difference
from run_je_sensitivity_20261015 import replace
from diagnose_je_directional_pe_20261015 import run as pe_audit
from localize_je_supg_work_20261015 import run as localize_supg

ROOT=Path(__file__).resolve().parents[1]
STUDY=ROOT/'runs/JE_SENSITIVITY_20261015'


def write(path,obj):
    with path.open('x',encoding='utf-8') as f:json.dump(obj,f,indent=2,ensure_ascii=False,allow_nan=False)


def field_same_mesh(a,b):
    x=np.genfromtxt(a/'nodes.csv',delimiter=',',names=True)
    y=np.genfromtxt(b/'nodes.csv',delimiter=',',names=True)
    assert all(np.array_equal(x[k],y[k]) for k in ['node','x_m','y_m','z_m'])
    values={k:difference(x[k],y[k]) for k in ['phi_V','theta_K']}
    scale=max(float(np.max(abs(d[k]))) for d in [x,y] for k in ['ux_m','uy_m','uz_m'])
    values.update({k:difference(x[k],y[k],scale) for k in ['ux_m','uy_m','uz_m']})
    values.update({k+'_density_perturbation':difference(np.expm1(x[k]),np.expm1(y[k])) for k in ['eta_n','eta_p']})
    return {'runs':[a.name,b.name],'comparisons':values,'passed':all(v['passed'] for v in values.values()),
            'nodes_hashes':[sha(a/'nodes.csv'),sha(b/'nodes.csv')],'scope':'same-node field shape; 2% existing tolerance'}


def main():
    out=STUDY/'analysis';out.mkdir(exist_ok=False)
    bm=json.loads((ROOT/'build/ARC_JE_FINAL_20261012/build_manifest.json').read_text(encoding='utf-8-sig'))
    sources={f['path']:sha(ROOT/f['path'])==f['sha256'] for f in bm['files'] if f['path'].startswith('src/')}
    assert all(sources.values()),sources
    preserved={};replay={}
    for voltage in ['V060','V080']:
        prior=json.loads((ROOT/'runs/JE_POWER_AUDIT_20261015'/voltage/'provenance.json').read_text(encoding='utf-8-sig'))
        preserved[voltage]={s['path']:sha(Path(s['path']))==s['sha256'] for s in prior['sources']}
        assert all(preserved[voltage].values()),preserved[voltage]
        replay[voltage]=sha(ROOT/'runs/JE_6_8_HALF_V_20261013'/voltage/'nodes.csv')==sha(STUDY/('A16_S1_'+voltage)/'nodes.csv')
        assert replay[voltage]
    rows=[];metrics={};input_checks={}
    directories=[STUDY/f'A16_S{s:g}_V{v*10:03d}' for s in [.5,1.,2.] for v in [6,8]]
    directories += [STUDY/f'A{count}_S1_V{v*10:03d}' for count in [24,32] for v in [6,8]]
    for d in directories:
        manifest=json.loads((d/'manifest.json').read_text());execution=json.loads((d/'execution.json').read_text())
        state,m=classify(d,execution['exit_code']);metrics[d.name]=m
        assert sha(d/'input.nml')==execution['input_sha256']
        original=Path(manifest['input_template']).read_text(encoding='utf-8-sig')
        actual=(d/'input.nml').read_text()
        for key in manifest['changed_input_keys']:
            original=replace(original,key,"'NORMALIZED'");actual=replace(actual,key,"'NORMALIZED'")
        assert original==actual,d.name
        input_checks[d.name]={'only_declared_keys_differ':True,'changed_keys':manifest['changed_input_keys']}
        if state!='SOLVE_FLUX_PASS_MESH_PENDING':
            rows.append({'run':d.name,'voltage_V':manifest['voltages'][0],
                         'angular_elements':int(re.search(r'A(\d+)',d.name)[1]),'radial_elements':14,
                         'thickness_elements':14,'supg_scale':manifest['supg_scale'],
                         'OMP_threads':int(manifest['runtime_env']['OMP_NUM_THREADS']),
                         'MKL_threads':int(manifest['runtime_env']['MKL_NUM_THREADS']),
                         'state':state,'seconds':execution['seconds'],'timed_out':execution['timed_out']})
            continue
        p=json.loads((d/'power_audit/decomposition.json').read_text())
        if not (d/'directional_pe.json').exists():pe_audit(d)
        pe=json.loads((d/'directional_pe.json').read_text())
        if not (d/'SUPG_localization_v2.json').exists():localize_supg(d)
        local=json.loads((d/'SUPG_localization_v2.json').read_text())
        assert abs(local['weighted_SUPG_W']-p['parts']['SUPG_weighted_W'])<1e-10*m['joule_W']
        row={'run':d.name,'voltage_V':m['voltage_V'],'angular_elements':int(re.search(r'A(\d+)',d.name)[1]),
             'radial_elements':14,'thickness_elements':14,'supg_scale':m['supg_scale'],'state':state,
             'OMP_threads':int(manifest['runtime_env']['OMP_NUM_THREADS']),
             'MKL_threads':int(manifest['runtime_env']['MKL_NUM_THREADS']),
             'domain_volume_m3':m['total_volume_m3'],
             'directional_Pe_p95':pe['volume_weighted_p95'],
             'JE_mW':1e3*m['joule_W'],'boundary_power_mW':1e3*m['terminal_power_W'],
             'gap_mW':1e3*p['gap_W'],'gap_percent':100*p['gap_relative'],
             'electrode_gap_percent':100*p['electrode_gap_relative'],
             'SUPG_work_mW':1e3*p['parts']['SUPG_weighted_W'],'SUPG_share_percent':100*p['SUPG_fraction_of_gap'],
             'SUPG_absolute_local_work_mW':1e3*local['absolute_local_work_W'],
             'SUPG_signed_to_absolute_ratio':local['signed_to_absolute_ratio'],
             'theta_max_K':m['theta_max_K'],'electrode_current_mA':1e3*(m['electron_right_A']+m['hole_right_A']),
             'stress_vm_max_MPa':1e-6*m['stress_von_mises_max_Pa'],
             'heat_balance_percent':100*m['heat_balance_relative'],'electron_balance_percent':100*m['electron_balance_relative'],
             'native_residual':m['maximum_accepted_residual'],'native_update':m['maximum_accepted_update'],
             'reconstructed_free_residual_A':p['free_residual_max_A'],'seconds':execution['seconds'],
             'timed_out':execution['timed_out']}
        rows.append(row)
    with (out/'results.csv').open('x',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    supg={};mesh={};fields={};pair_metrics={}
    for v in [6,8]:
        runs=[STUDY/f'A16_S{s:g}_V{v*10:03d}' for s in [.5,1.,2.]]
        supg[v]=stabilization_audit(runs);write(out/f'SUPG_V{v:02d}.json',supg[v])
        a,b=(metrics[f'A16_S{s:g}_V{v*10:03d}'] for s in [.5,1.])
        pair_metrics[v]={k:{'relative_span':abs(a[k]-b[k])/max(abs(a[k]),abs(b[k]),1e-300),
                            'limit':SUPG_TOLERANCE,'values':[a[k],b[k]]} for k in SUPG_METRICS}
        for value in pair_metrics[v].values():value['passed']=value['relative_span']<value['limit']
        write(out/f'SUPG_pair_S0.5_1_V{v:02d}.json',{'scope':'Two completed coefficients; does not replace the three-coefficient audit','comparisons':pair_metrics[v]})
        for s in [.5,2.]:
            key=f'SUPG_S{s:g}_vs_1_V{v:02d}'
            other=STUDY/f'A16_S{s:g}_V{v*10:03d}'
            if json.loads((other/'summary.json').read_text())['state']=='SOLVE_FLUX_PASS_MESH_PENDING':
                fields[key]=field_same_mesh(other,STUDY/f'A16_S1_V{v*10:03d}')
            else:
                fields[key]={'passed':False,'comparisons':{},'status':'NO_ACCEPTED_COEFFICIENT_STATE'}
            write(out/(key+'.json'),fields[key])
        levels=[STUDY/f'A{a}_S1_V{v*10:03d}' for a in [16,24,32]]
        mesh[v]=mesh_audit(levels,[v]);write(out/f'angular_mesh_V{v:02d}.json',mesh[v])
        for a,b in [(16,24),(24,32)]:
            key=f'field_A{a}_A{b}_V{v:02d}'
            fields[key]=field_audit(STUDY/f'transfer_A{a}_A{b}_V{v*10:03d}',STUDY/f'A{a}_S1_V{v*10:03d}',STUDY/f'A{b}_S1_V{v*10:03d}')
            write(out/(key+'.json'),fields[key])
    result={'rows':rows,'SUPG':supg,'SUPG_pair_S0.5_1':pair_metrics,'angular_mesh':mesh,'fields':fields,'source_matches_build':sources,
            'input_checks':input_checks,'prior_evidence_preserved':preserved,'baseline_nodes_byte_identical':replay,'fully_validated':False,
            'limitations':['Angular-only refinement is not a full 3D mesh convergence study',
                           'J.E branch remains a current-work diagnostic, not a verified thermodynamic closure',
                           'Field shape audit uses nodal states, not continuity of gradient-derived fields',
                           'Existing solve/flux gates do not include an electrical-power-gap acceptance criterion']}
    write(out/'summary.json',result)
    lines=['J·E 分支：6 V / 8 V 网格与 SUPG 敏感性审计',
           '实验起始日期：2026-10-15。全部结果来自原生求解或明确标记的失败尝试。',
           '本轮只覆盖 6 V 和 8 V，不替代 6.5 / 7 / 7.5 V 在最终选定网格与方法上的重算。',
           '',
           '控制条件：相同可执行文件、材料、边界条件、J·E 热源、完整载荷；未修改源码或既有数值门限。',
           '新网格仅改变周向单元数 16/24/32，径向和厚度均保持 14。不能当作完整三维网格认证。',
           '计算资源：A16/A24 使用4线程，A32 单独运行并使用8线程；同一构建和数值判据，具体环境见每项 manifest.json。',
           '原网格比较 SUPG 系数 0.5/1/2。插值只作初值，每个新结果均需原生重新平衡。',
           '',
           '相对功率差额 = (∫Ω J·E dV - P边界)/∫Ω J·E dV。下表保留符号。',
           'P边界为全外边界 -∫∂Ω φJ·n dA；电极差额只积分两个电极。',
           '功率差额没有在本次新增 PASS 门限；原求解/通量 PASS 不等于电功率或空间精度认证。',
           '',
           '算例 | J·E [mW] | 全边界差额 [%] | 仅电极差额 [%] | 最大温升 [K] | 电极电流 [mA] | 峰值VM应力 [MPa]']
    by_name={r['run']:r for r in rows}
    overview=[]
    for v in [6,8]:
        series=[by_name[f'A{n}_S1_V{v*10:03d}'] for n in [16,24,32]]
        if all('gap_percent' in r for r in series):
            overview.append(f"固定 SUPG=1，{v} V 的 16→24→32 周向网格差额："+' → '.join(f"{r['gap_percent']:.6f}%" for r in series))
            overview.append(f"  Pe 的体积加权95分位：{series[0]['directional_Pe_p95']:.4f} → {series[2]['directional_Pe_p95']:.4f}；SUPG局部绝对贡献积分：{series[0]['SUPG_absolute_local_work_mW']:.8f} → {series[2]['SUPG_absolute_local_work_mW']:.8f} mW。")
    complete_rows=[r for r in rows if 'JE_mW' in r]
    overview.append(f"{len(complete_rows)} 个状态通过既有求解/通量检查（含2个基准重放）；另有 {len(rows)-len(complete_rows)} 个未获得合格完整解的尝试，全部保留。")
    overview.append('本轮 SUPG 三系数检查：'+str({v:s['covered_metrics_passed'] for v,s in supg.items()})+'；周向网格指标检查：'+str({v:g['all_mesh_pass'] for v,g in mesh.items()})+'。False 表示未通过或证据不完整。')
    if complete_rows:
        shares=[r['SUPG_share_percent'] for r in complete_rows]
        overview.append(f"成功状态中，电势加权 SUPG 项解释功率差额的 {min(shares):.3f}%～{max(shares):.3f}%。这是离散残差的分解，不是新增热源。")
    overview.append('建议：保持热源与默认设置不动，继续细网格上的SUPG交叉对照，以及径向/厚度方向的独立细化；现有结果不能作为完整空间精度或物理闭合认证。')
    lines[2:2]=['','主要观察：',*overview,'']
    for r in rows:
        if 'JE_mW' not in r:
            lines.append(r['run']+' | 未获得通过既有检查的结果：'+r['state']+f"；耗时={r['seconds']:.1f}s；触及运行时限={r['timed_out']}（不据此判定不存在解）")
        else:
            lines.append(f"{r['run']} | {r['JE_mW']:.8f} | {r['gap_percent']:.6f} | {r['electrode_gap_percent']:.6f} | {r['theta_max_K']:.6f} | {r['electrode_current_mA']:.8f} | {r['stress_vm_max_MPa']:.6f}")
    lines.extend(['','原有判据审计（没有修改门限）：'])
    for v in [6,8]:
        s=supg[v];g=mesh[v]
        bad=[k for k,c in s.get('comparisons',{}).items() if not c['passed']]
        lines.append(f"{v} V SUPG 三系数指标：{'通过' if s['covered_metrics_passed'] else '未通过'}；不合格指标={bad}；证据错误={s['errors']}")
        for k,c in pair_metrics[v].items():
            if not c['passed']:lines.append(f"  已完成的 0.5 与 1 对照：{k} 差异={100*c['relative_span']:.6f}% > 原 1% 门限；三系数记录不齐时仍保留这一实测差异。")
        lines.append(f"{v} V 周向三网格指标：{'通过' if g['all_mesh_pass'] else '未通过'}（只适用于周向序列；仍需局部场、其他网格方向和稳定化检查）")
        for r in g['voltages']:
            for k,c in r['changes'].items():
                lines.append(f"  {k}: 16→24={100*c['L0_L1']:.6f}%, 24→32={100*c['L1_L2']:.6f}%, 原门限={100*c['limit']:g}%")
    lines.extend(['','局部场比较：以两场的最大绝对幅值作统一归一化，原有最大差门限 2%。载流子比较 exp(eta)-1。'])
    for key,value in fields.items():
        lines.append(key+': '+('通过' if value['passed'] else '未通过'))
        for k,c in value.get('comparisons',{}).items():
            lines.append(f"  {k}: 最大差={100*c['maximum_relative']:.6f}%, RMS差={100*c['rms_relative']:.6f}%")
    lines.extend(['',
                  '功率分解：W5 - P外边界,5 = (W5-W3) + φᵀS3 + (P反力,3-P外边界,5) - φ自由ᵀR自由,3。',
                  '每个通过算例的 power_audit/ 保存 3阶/5阶独立重构、分解和边界分类。原生热源没有加入 SUPG 补偿项。',
                  f'构建一致性：{len(sources)} 个源码文件与所用构建的 SHA256 全部匹配。输入仅允许 case_id、网格、初值、初值载荷标记和 SUPG 系数不同。',
                  '旧证据保留：6 V / 8 V 各 12 项上一轮审计的来源哈希仍匹配；本次系数1原网格重放的节点文件与旧结果逐字节相同。',
                  'J·E 分支仍是电流做功诊断，不能由数值功率改善推断完整热电能量闭合或实验准确性。',
                  '节点场比较不认证梯度场连续性，峰值应力另行报告，不能由位移通过替代应力收敛。',
                  '',
                  '复现：从项目根目录运行 scripts/ 中原脚本，或从 solver_v2 下使用 scripts/ 路径；必须选择全新输出目录。',
                  '基础系数批次：python scripts/run_je_sensitivity_20261015.py --base-batch --study-root runs/<新目录>',
                  '网格：python scripts/generate_hv_surface_mesh.py --xy-level 0 --z-level 2 --angular-count 24 --output runs/<新网格目录>',
                  '精确重放本研究 A32/6V 的初值与资源配置（在 solver_v2 下执行，输出必须不存在）：',
                  'python scripts/run_je_sensitivity_20261015.py --voltage 6 --scale 1 --mesh runs/JE_SENSITIVITY_20261015/mesh_A32/H1_XY0_Z2.mesh --seed runs/JE_SENSITIVITY_20261015/transfer_A24_A32_V060/nodes.csv --output runs/JE_SENSITIVITY_REPLAY_NEW/A32_S1_V060 --threads 8 --timeout-seconds 7200',
                  '其余输入、网格、初值来源、可执行文件哈希、运行环境、退出码、审计脚本快照均保存在本研究目录。'])
    (out/'REPORT.txt').write_text('\n'.join(lines),encoding='utf-8')
    for script in ['run_je_sensitivity_20261015.py','audit_je_sensitivity_power_20261015.py',
                   'complete_je_sensitivity_20261015.py','localize_je_supg_work_20261015.py',
                   'diagnose_je_directional_pe_20261015.py',
                   'summarize_je_sensitivity_20261015.py','diagnose_je_power_gap_20261015.py',
                   'plot_je_8v_clouds_20261013.py','generate_hv_surface_mesh.py','generate_f1_main_cases.py',
                   'audit_hex20_boundary_geometry.py','transfer_hv_initial_guess.py','audit_hv_stabilization.py',
                   'audit_hv_mesh.py','audit_hv_field_shape.py','run_hv_sweep.py','run_hv_regression.py','hv_runtime.py']:
        shutil.copy2(ROOT/'scripts'/script,out/script)
    for name in ['build_manifest.json','compiler_options.json']:
        shutil.copy2(ROOT/'build/ARC_JE_FINAL_20261012'/name,out/name)
    shutil.copy2(STUDY/'study_contract.json',out/'study_contract.json')
    write(out/'evidence_hashes.json',{str(p.relative_to(STUDY)):sha(p) for p in STUDY.rglob('*')
                                      if p.is_file() and p.suffix in ['.json','.nml','.py','.txt','.csv','.mesh'] and 'pardiso_ooc' not in p.parts})
    print(json.dumps({'rows':rows,'SUPG_pass':{v:r['covered_metrics_passed'] for v,r in supg.items()},
                      'angular_mesh_pass':{v:r['all_mesh_pass'] for v,r in mesh.items()},
                      'field_pass':{k:r['passed'] for k,r in fields.items()}},indent=2),flush=True)


if __name__=='__main__':main()
