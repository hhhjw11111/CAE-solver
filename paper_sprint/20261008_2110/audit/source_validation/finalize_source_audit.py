from pathlib import Path
import json, hashlib, subprocess, csv, re, os
out=Path(__file__).absolute().parent
root=out.parents[3]
val=root/'paper_sprint/20261008_2110/validation/runtime_retry'
audit=json.loads((out/'source_audit_runtime_retry.json').read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
refs=[
('src/core/field_layout_m.f90',5,15,'七场排序'),
('src/core/constants_m.f90',6,7,'q与kB'),
('src/physics/zno_candidate_m.f90',8,45,'全部ZnO候选常数与模型开关'),
('src/physics/transport_m.f90',40,58,'N0冻结及NT温导数'),
('src/physics/thermo_piezo_m.f90',11,24,'热压电本构'),
('src/physics/material_point_m.f90',107,149,'密度、载流子电流、热源'),
('src/physics/material_point_m.f90',151,204,'一致材料导数及LQ替换'),
('src/physics/recombination_m.f90',10,18,'R=0'),
('src/fem/element_kernel_m.f90',63,140,'积分、生产残差/切线和格式入口'),
('src/fem/element_kernel_m.f90',143,176,'log-weighted测试及其状态导数'),
('src/fem/element_kernel_m.f90',178,236,'普通弱式状态和测试算子'),
('src/app/f1_case_runner_m.f90',157,178,'LQ三阶段'),
('src/app/f1_case_runner_m.f90',315,329,'原验收门槛'),
('src/app/f1_case_runner_m.f90',543,595,'全域电极与A/B/C'),
('src/app/f1_case_runner_m.f90',598,657,'半域对称与刚体约束'),
('src/app/f1_case_runner_m.f90',659,680,'LQ固定字段'),
('src/io/native_input_m.f90',59,62,'新格式默认关闭'),
('src/verification/verification_driver_m.f90',41,76,'真实verify入口'),
('src/app/main.f90',29,90,'version/verify/run-f1命令分发'),
('tests/hv_log_weighted_production_bar.f90',89,113,'生产核条带采样误差定义'),
('tests/hv_log_weighted_production_bar.f90',123,153,'自定义缩约装配调用生产单元'),
('tests/hv_log_weighted_production_bar.f90',187,213,'实际Pe/网格及容差'),
]
lines=['# 精确源码索引','','路径以 solver_v2 为根；真实文件绝对路径同时写入 source_index.json。行号为本轮文件中的1基行号。','','| 路径与行 | 内容 | SHA-256 |','|---|---|---|']
items=[]
for f,a,b,desc in refs:
    p=root/f
    lines.append(f'| `{f}:{a}–{b}` | {desc} | `{sha(p)}` |')
    items.append({'path':str(p),'start':a,'end':b,'description':desc,'sha256':sha(p)})
(out/'source_index.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
(out/'source_index.json').write_text(json.dumps(items,ensure_ascii=False,indent=2),encoding='utf-8')
contract=(out/'source_contract.md').read_text(encoding='utf-8').split('\n\n## 源码精确引用')[0]
(out/'source_contract.md').write_text(contract+'\n\n## 源码精确引用\n\n以下文件全部含SHA-256，完整表也见 source_index.md。\n\n'+'\n'.join(lines[4:])+'\n',encoding='utf-8')
ifx=Path('C:/Program Files (x86)/Intel/oneAPI/compiler/latest/bin/ifx.exe')
p=subprocess.run([str(ifx),'/QV'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
(val/'ifx_version.txt').write_bytes(p.stdout)
mkl=Path('C:/Program Files (x86)/Intel/oneAPI/mkl/latest/include/mkl_version.h')
(val/'mkl_version_header.txt').write_bytes(mkl.read_bytes())
commands=json.loads((val/'command_manifest.json').read_text(encoding='utf-8'))
summary={'source_build_match_count':59,'source_build_mismatches':[],'current_executable_version':'0.15.9-log-inverse-checked','current_executable_sha256':sha(root/'build/HV_LOG_INVERSE_20261008_R25_CLEAN/TPSolverV2.exe'),'ifx_version':p.stdout.decode(errors='replace').strip(),'mkl_version':'2024.0.2; header build date 20240722','fresh_small_checks':[], 'spec_7_status':'PARTIAL_REUSED_EVIDENCE','gaps':['A=0 zero-field case','same-electron negative field case','relative L2 independent quadrature','quadrature refinement stability','saved full internal profile, free-row residual and separate endpoint-current metrics','global CSR assembly not exercised by custom-reduced bar driver'], 'runtime_first_attempt':{'exit_code':3221225781,'hex':'0xC0000135','classification':'MISSING_RUNTIME_SEARCH_PATH','logs_preserved':'../command_manifest.json'},'runtime_recovery':'Added installed Intel compiler/latest/bin to child PATH; no installation or production modification.'}
for c in commands:
    logfile=Path(c['stdout']['path']);text=logfile.read_text(errors='replace')
    summary['fresh_small_checks'].append({'command':c['command'],'exit_code':c['exit_code'],'verification_result':'PASS' if 'verification_result=PASS' in text else 'NOT_APPLICABLE_VERSION','log':str(logfile),'sha256':sha(logfile)})
bar=(val/'production_bar_replay.txt').read_text(errors='replace')
rows=[]
for match in re.finditer(r'production_Pe1000_carrier=(\d),errors=\s*([\d.E+-]+)\s+([\d.E+-]+)',bar):
    carrier,density,current=match.groups();rows.append({'kind':'PRODUCTION_KERNEL_ARCHIVED_HARNESS_R22_FRESH_REPLAY','carrier':carrier,'domain_Pe':1000,'A':2000,'elements':128,'sample_max_density_error':float(density),'sample_max_current_error':float(current),'metric_definition':'11 points per element, density normalized by 2; current by max(1,abs(exact)); not L2','source':str(val/'production_bar_replay.txt')})
summary['production_bar_replay']=rows
summary['production_bar_byte_identical_to_prior_stdout']=commands[-1]['byte_identical_to_prior_stdout']
(val/'verification_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
with (val/'production_bar_errors.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
notes=['# 本輪小型核验与部分完成边界','','当前R25_CLEAN 的 version、V01、V02、G02、F1M、G03 均退出0；各verify原输出为PASS。R22生产HEX20解析条带归档可执行文件本轮重新运行退出0。所有命令与原始日志列在 command_manifest.json。','','| 新鲜重放条带 | 最大采样浓度误差 | 最大采样电流误差 |','|---|---|---|']
for r in rows:notes.append(f"| carrier {r['carrier']}，Pe=1000，128元 | {r['sample_max_density_error']:.9g} | {r['sample_max_current_error']:.9g} |")
notes.extend(['','该Pe定义为mu E L/(2D)，附件A=mu E L/D=2Pe；不能称最大单元Pe=1000，实测最大单元Pe约23.255。初值内部eta=ln1.5，两端0/ln2，内部确实求解，未将解析浓度赋给全部节点。空穴使用电流扩散负号及镜像先验边界层网格，不等于电子反向电场验证。','', '原型48组Pe=0.1–1000采样表另存archived_bar_errors.csv，严格归属R21原型，不能冒称全部走R22生产实现。原型总体FAIL因16元低Pe新旧格式差7.66027e−7超过1e−7；32/64元新旧差5.45523e−8/3.34532e−9通过不抹去原失败。','', '附件第7节所需零场、电子负场、L2误差独立积分加密、完整剖面和自由行残差独立输出仍未完成。本轮复用/重放证明范围是生产材料点/单元核与该缩约条带，不是完整3D或全局CSR求解器证书。考虑原有大进程与内存限制，本轮不新编核心或启动三维作业。','', '首次未设置Intel DLL搜索路径导致7次调用均0xC0000135；修正仅在本任务子进程PATH中追加已安装compiler/latest/bin后恢复。首次日志在父validation目录保留。编译器只查询版本，未执行构建：ifx 2024.2.1 Build 20240711；MKL安装头文件为2024.0.2、20240722。现有R25构建采用/O2、/check:none、/Qopenmp，线程库；本轮运行OMP_NUM_THREADS=MKL_NUM_THREADS=1。'])
(val/'VERIFICATION_SCOPE.md').write_text('\n'.join(notes)+'\n',encoding='utf-8')
print(json.dumps({'checks':len(commands),'production_rows':rows,'byte_identical':summary['production_bar_byte_identical_to_prior_stdout']}))
