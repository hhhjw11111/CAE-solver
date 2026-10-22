"""Audit scheduling from native checkpoints; preserve failed hole postprocessing."""
import json,sys
from pathlib import Path
HERE=Path(__file__).absolute().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from run_hv_sweep import read_metrics,sha
test=HERE/'heat_restart_native_check_v2'
previous=json.loads((test/'report.json').read_text())
checkpoints=[]
for name in ('continuous','resumed'):
    path=sorted((test/name).glob('accepted_*_load.txt'))[-1]
    load=read_metrics(path)
    checkpoints.append(dict(path=str(path),sha256=sha(path),load=load))
native=all(c['load']['voltage_factor']==1. and c['load']['joule_factor']==1. and
           c['load']['residual']<=1e-9 and c['load']['update']<=1e-9 for c in checkpoints)
report=dict(scope='ONE_ELEMENT_HEAT_RESTART_ALGORITHM_CHECK_ONLY',
    native_checkpoint_gates_passed=native,field_agreement_passed=max(previous['scaled_field_difference'].values())<1e-8,
    maximum_scaled_field_difference=max(previous['scaled_field_difference'].values()),checkpoints=checkpoints,
    original_report_sha256=sha(test/'report.json'),original_report_preserved=True,
    full_case_or_flux_check_passed=False,exit_codes=previous['exit_codes'],
    postprocess_failure='Both plain-plate fixtures return Hole field metric input is invalid after f1_solve=COMPLETE.',
    high_voltage_or_mesh_validation=False)
report['heat_restart_algorithm_check_passed']=native and report['field_agreement_passed']
dest=HERE/'heat_restart_checkpoint_audit.json'
with dest.open('x') as f:json.dump(report,f,indent=2,allow_nan=False)
print(json.dumps({k:v for k,v in report.items() if k!='checkpoints'}))
