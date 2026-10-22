"""Read-only collection of fine-grid JE evidence; no tolerance changes."""
from pathlib import Path
import argparse,csv,json,shutil
import numpy as np
from run_je_fine_cross_20261016 import STUDY,OLD,ROOT,base,sha,save,accepted
from run_je_sensitivity_20261015 import replace
from run_hv_sweep import classify,read_metrics
from audit_hv_mesh import audit as mesh_audit
from audit_hv_stabilization import audit as supg_audit,METRICS,TOLERANCE
from audit_hv_field_shape import audit as field_audit
from summarize_je_sensitivity_20261015 import field_same_mesh

def locations(a,b):
    x=np.genfromtxt(a,delimiter=',',names=True);y=np.genfromtxt(b,delimiter=',',names=True)
    assert len(x)==len(y)
    result={}
    for eta,density in [('eta_n','n_m3'),('eta_p','p_m3')]:
        p,q=np.expm1(x[eta]),np.expm1(y[eta]);delta=abs(q-p);i=int(np.argmax(delta))
        result[eta]={'maximum_difference_node':int(y['node'][i]),
                     'xyz_um':[float(y[k][i])*1e6 for k in ['x_m','y_m','z_m']],
                     'source_perturbation_at_node':float(p[i]),'target_perturbation_at_node':float(q[i]),
                     'maximum_absolute_perturbation_difference':float(delta[i]),
                     'target_density_min_m3':float(y[density].min()),'target_density_max_m3':float(y[density].max())}
    return result

def main(output):
    output.mkdir(exist_ok=False)
    contract=json.loads((STUDY/'contract.json').read_text(encoding='utf-8'))
    preserved=json.loads((STUDY/'preserved_hashes.json').read_text(encoding='utf-8'))
    preservation={p:sha(Path(p))==h for p,h in preserved.items()}
    assert all(preservation.values()),'Old evidence or source was changed'
    baseline=[dict(name=base(v).name,v=v,scale=1.,counts=[32,14,14],kind='baseline',directory=str(base(v))) for v in [6,8]]
    rows=[];checks={};fields={};mesh={};supg={};pairs={};byname={}
    for j in baseline+contract['jobs']:
        d=Path(j.get('directory',str(STUDY/j['name'])));byname[j['name']]=d
        row=dict(run=j['name'],kind=j['kind'],voltage_V=j['v'],angular=j['counts'][0],radial=j['counts'][1],thickness=j['counts'][2],SUPG=j['scale'])
        if not (d/'summary.json').exists():
            row['state']='PENDING_OR_DEPENDENCY_BLOCKED';rows.append(row);continue
        ex=json.loads((d/'execution.json').read_text());mft=json.loads((d/'manifest.json').read_text())
        state,m=classify(d,ex['exit_code'])
        recorded_state=json.loads((d/'summary.json').read_text())['state']
        if recorded_state=='FAILED_FIELD_AUDIT':state=recorded_state
        row.update(state=state,seconds=ex['seconds'],timed_out=ex['timed_out'])
        assert ex['input_sha256']==sha(d/'input.nml') and ex['executable_sha256']==sha(Path(mft['executable']))
        assert mft['mesh_sha256']==sha(Path(mft['mesh']))
        a=Path(mft['input_template']).read_text(encoding='utf-8-sig');b=(d/'input.nml').read_text()
        for key in mft['changed_input_keys']:a=replace(a,key,"'NORMALIZED'");b=replace(b,key,"'NORMALIZED'")
        assert a==b
        checks[j['name']]={'input_only_declared_changes':True,'executable_matches':True,'mesh_matches':True}
        if state=='SOLVE_FLUX_PASS_MESH_PENDING':
            assert json.loads((d/'field_audit.json').read_text())['passed']
            row.update(JE_mW=1e3*m['joule_W'],boundary_mW=1e3*m['terminal_power_W'],theta_K=m['theta_max_K'],
                       electrode_mA=1e3*(m['electron_right_A']+m['hole_right_A']),VM_MPa=1e-6*m['stress_von_mises_max_Pa'],
                       native_residual=m['maximum_accepted_residual'],native_update=m['maximum_accepted_update'],
                       ring_n=m['ring_abs_n_relative_mean'],ring_p=m['ring_abs_p_relative_mean'])
            if (d/'power_audit/decomposition.json').exists():
                p=json.loads((d/'power_audit/decomposition.json').read_text())
                row.update(gap_percent=100*p['gap_relative'],electrode_gap_percent=100*p['electrode_gap_relative'],
                           SUPG_share_percent=100*p['SUPG_fraction_of_gap'],SUPG_mW=1e3*p['parts']['SUPG_weighted_W'],
                           free_residual_A=p['free_residual_max_A'])
            else:row['power_audit_status']='PENDING'
            if j['kind'] in ['radial','thickness']:
                source=Path(j['source']);mapping=Path(j['mapping'])
                f=field_audit(mapping,source,d);f['carrier_locations']=locations(mapping/'nodes.csv',d/'nodes.csv')
                fields[j['name']]=f
            elif j['kind']=='SUPG':
                f=field_same_mesh(base(j['v']),d);f['carrier_locations']=locations(base(j['v'])/'nodes.csv',d/'nodes.csv')
                fields[j['name']]=f
        rows.append(row)
    for v in [6,8]:
        sr=[STUDY/f'A32_R14_Z14_S0.5_V{v*10:03d}',base(v),STUDY/f'A32_R14_Z14_S2_V{v*10:03d}']
        supg[str(v)]=supg_audit(sr)
        for scale,d in [(.5,sr[0]),(2.,sr[2])]:
            if accepted(d):
                a,b=[read_metrics(p/'metrics.txt') for p in [base(v),d]]
                comparison={}
                for k in METRICS:
                    span=abs(a[k]-b[k])/max(abs(a[k]),abs(b[k]),1e-300)
                    comparison[k]=dict(values=[a[k],b[k]],relative_span=span,limit=TOLERANCE,passed=span<TOLERANCE)
                pairs[f'V{v}_S{scale:g}']=comparison
        for direction in ['radial','thickness']:
            levels=[base(v)]+[STUDY/f'A32_R{n if direction=="radial" else 14}_Z{n if direction=="thickness" else 14}_S1_V{v*10:03d}' for n in [16,18]]
            mesh[f'{direction}_V{v}']=mesh_audit(levels,[v])
    result=dict(rows=rows,fields=fields,mesh=mesh,SUPG=supg,SUPG_pairs=pairs,input_checks=checks,
                prior_files_unchanged=len(preservation),fully_validated=False,contract=contract)
    save(output/'summary.json',result)
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with (output/'results.csv').open('x',encoding='utf-8-sig',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=keys);writer.writeheader();writer.writerows(rows)
    lines=['J·E 新分支：细网格 SUPG 交叉与独立径向/厚度细化',
           '原模型、热源、源码、构建和验收阈值均保持不变；历史证据哈希已复核。',
           '基准为周向32、径向14、厚度14。径向14/16/18与厚度14/16/18为两项独立研究，不是联合细网格。',
           '本轮只检查6 V与8 V；不认证整个电压区间、周向收敛、应力收敛或完整热电闭合。',
           '', '功率差额=(J·E体积分−全外边界电功率)/(J·E体积分)，另列仅电极差额。未新增或修改功率验收门限。',
           '算例 | 状态 | 全边界差额% | 电极差额% | 温升K | 电流mA']
    for r in rows:
        vals=[f"{r[k]:.7g}" if k in r else '—' for k in ['gap_percent','electrode_gap_percent','theta_K','electrode_mA']]
        lines.append(r['run']+' | '+r['state']+' | '+' | '.join(vals))
    lines += ['', '三系数 SUPG 原有指标验收：']
    for key,value in supg.items():
        lines.append(f"{key} V: passed={value['covered_metrics_passed']}; errors={value['errors']}")
        for k,c in value['comparisons'].items():lines.append(f"  {k}: span={100*c['relative_span']:.6f}%, limit={100*c['limit']:g}%, passed={c['passed']}")
    lines += ['', '已完成的单系数与1的对照（不替代三系数验收）：']
    for name,comparison in pairs.items():
        lines.append(name)
        for k,c in comparison.items():lines.append(f"  {k}: span={100*c['relative_span']:.6f}%, limit={100*c['limit']:g}%, passed={c['passed']}")
    lines += ['', '独立方向网格指标验收（末两层比较）：']
    for key,value in mesh.items():
        lines.append(f"{key}: passed={value['all_mesh_pass']}; provenance_errors={value['provenance_errors']}")
        for vr in value['voltages']:
            lines.append('  status='+vr['status'])
            for k,c in vr['changes'].items():lines.append(f"  {k}: 14→16={100*c['L0_L1']:.6f}%, 16→18={100*c['L1_L2']:.6f}%, limit={100*c['limit']:g}%")
    lines += ['', '局部场比较：最大差/两场最大幅值；载流子采用exp(eta)-1，门限2%。映射是比较或初值，必须原生重新求解。']
    for key,f in fields.items():
        lines.append(key+f": passed={f['passed']}")
        for k,c in f['comparisons'].items():lines.append(f"  {k}: max={100*c['maximum_relative']:.6f}%, RMS={100*c['rms_relative']:.6f}%, passed={c['passed']}")
        for k,c in f['carrier_locations'].items():lines.append(f"  {k} 最大差位置[μm]={c['xyz_um']}; 目标浓度范围[m^-3]={c['target_density_min_m3']:.6g}..{c['target_density_max_m3']:.6g}")
    lines += ['',f'保全检查：{len(preservation)}个旧文件/源码/构建的SHA256均一致。',
              '超时是资源限制内未获得完整接受状态，不表示稳态不存在。新结果及失败记录保存在独立目录。',
              '主执行脚本：scripts/run_je_fine_cross_20261016.py；汇总脚本：scripts/summarize_je_fine_cross_20261016.py --output <全新目录>。',
              '输入来源、初值、网格、运行环境、退出码、场审计和独立功率分解均可逐例检查。']
    (output/'REPORT.txt').write_text('\n'.join(lines),encoding='utf-8')
    shutil.copy2(__file__,output/Path(__file__).name)
    save(output/'artifact_hashes.json',{str(p.relative_to(STUDY)):sha(p) for p in STUDY.rglob('*') if p.is_file() and p.suffix in ['.json','.csv','.nml','.mesh','.py','.txt'] and 'pardiso_ooc' not in p.parts})
    print(json.dumps({'rows':len(rows),'SUPG_pass':{k:v['covered_metrics_passed'] for k,v in supg.items()},'mesh_pass':{k:v['all_mesh_pass'] for k,v in mesh.items()},'fields_pass':{k:v['passed'] for k,v in fields.items()}},ensure_ascii=False),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args().output)
