"""Report full-load fixed-mesh solutions separately from preliminary flux checks.

This does not relax native gates or certify spatial accuracy. A converged field
with a failed flux check stays explicitly failed, but is not confused with a
voltage point at which no complete solution was obtained.
"""
import argparse,csv,json,os
from pathlib import Path
from run_hv_sweep import sha
from summarize_hv_credibility import point_record

REQUESTED=[6.,6.5,7.,7.5,8.]

def preliminary_status(record):
    complete=(record.get('state') in ('SOLVE_FLUX_PASS_MESH_PENDING','SOLVED_FLUX_GATE_FAILED')
        and not record.get('evidence_errors')
        and record.get('field_audit_recomputed',{}).get('passed') is True)
    flux=complete and record['state']=='SOLVE_FLUX_PASS_MESH_PENDING'
    return dict(full_load_solution_available=complete,preliminary_flux_check_passed=flux,
        status=('FIXED_MESH_SOLUTION_BASIC_CHECKS_PASS' if flux else
                'FIXED_MESH_SOLUTION_FLUX_CHECK_FAILED' if complete else 'NO_ACCEPTABLE_COMPLETE_FIELD'),
        spatial_accuracy_certified=False)

def summarize(records):
    rows=[]
    for v in sorted(set(REQUESTED+[r['voltage_V'] for r in records])):
        candidates=[dict(r,**preliminary_status(r)) for r in records if r['voltage_V']==v]
        passed=[r for r in candidates if r['preliminary_flux_check_passed']]
        solved=[r for r in candidates if r['full_load_solution_available']]
        if candidates:row=(passed or solved or candidates)[-1]
        else:row=dict(voltage_V=v,status='NO_COMPLETE_RESULT_YET',full_load_solution_available=False,
                      preliminary_flux_check_passed=False,spatial_accuracy_certified=False)
        rows.append(row)
    requested=[r for r in rows if r['voltage_V'] in REQUESTED]
    return dict(scope='FIXED_MESH_FULL_LOAD_RESULTS_AND_PRELIMINARY_CHECKS_ONLY',points=rows,
        all_requested_full_load_solutions_available=all(r['full_load_solution_available'] for r in requested),
        all_requested_preliminary_checks_passed=all(r['preliminary_flux_check_passed'] for r in requested),
        spatial_accuracy_certified=False,complete_quantitative_0_8V=False)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-report',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    raw=args.source_report.read_bytes();source=json.loads(raw)
    records=[]
    for r in source['records']:
        directory=Path(os.path.abspath(r['directory']))
        records.append(point_record(r['voltage_V'],directory,directory/'manifest.json'))
    report=summarize(records)
    import hashlib
    report.update(source_report=str(args.source_report.absolute()),source_report_sha256=hashlib.sha256(raw).hexdigest(),
        reporter_sha256=sha(Path(__file__)),source_status=source.get('status'))
    args.output.mkdir(parents=True,exist_ok=False)
    (args.output/'preliminary_report.json').write_text(json.dumps(report,indent=2,allow_nan=False))
    with (args.output/'results_6_8V.csv').open('w',newline='',encoding='utf-8-sig') as stream:
        keys=['voltage_V','status','full_load_solution_available','preliminary_flux_check_passed',
              'theta_max_K','joule_half_W','electron_right_half_A','physical_exterior_heat_balance_relative','directory']
        writer=csv.DictWriter(stream,fieldnames=keys);writer.writeheader()
        for row in report['points']:
            metrics=row.get('metrics',{})
            writer.writerow({k:(metrics.get(k) if k in ['theta_max_K','physical_exterior_heat_balance_relative']
                else metrics.get('joule_W') if k=='joule_half_W' else metrics.get('electron_right_A') if k=='electron_right_half_A'
                else row.get(k)) for k in keys})
    print(json.dumps({k:report[k] for k in ['all_requested_full_load_solutions_available','all_requested_preliminary_checks_passed']}))

if __name__=='__main__':main()
