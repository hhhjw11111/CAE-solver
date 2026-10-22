"""Scale a near-uniform zero-heat checkpoint as an INITIAL GUESS ONLY.

The target native solve must re-equilibrate it and restore heat to one. This is
never a voltage result, interpolation evidence, or an accepted checkpoint.
"""
import argparse,csv,json,math
from pathlib import Path
from prepare_hv_zero_heat_guess import zero_heat_rows
from run_hv_sweep import read_metrics,sha


def voltage_guess(rows,load,target_voltage):
    source=abs(load.get('applied_voltage_V',0.))
    if (load.get('joule_factor')!=0. or not math.isfinite(source) or source<=0 or
        not math.isfinite(target_voltage) or not source<=target_voltage<=8.):
        raise ValueError('Require zero heat and 0 < source voltage <= target <= 8')
    for key in ('residual','update'):
        value=load.get(key,float('inf'))
        if not math.isfinite(value) or not 0<=value<=1e-9:
            raise ValueError('Source checkpoint has not met the original Newton gates')
    result,backgrounds=zero_heat_rows(rows)
    for row in rows:
        if any(abs(float(row[k]))>1e-10 for k in ('eta_n','eta_p','theta_K')):
            raise ValueError('Require a cold near-uniform source; heated fields are not scaled')
    phi=[float(r['phi_V']) for r in rows]
    if abs(min(phi)+source)>1e-8 or abs(max(phi))>1e-8:
        raise ValueError('Source electrode range does not match the checkpoint voltage')
    ratio=target_voltage/source
    for old,new in zip(rows,result):
        for key in ('phi_V','ux_m','uy_m','uz_m'):
            new[key]=f'{float(old[key])*ratio:.17e}'
    return result,backgrounds,ratio


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source-nodes',type=Path,required=True)
    p.add_argument('--source-load',type=Path,required=True)
    p.add_argument('--voltage',type=float,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    with a.source_nodes.open(newline='',encoding='utf-8-sig') as stream:
        reader=csv.DictReader(stream);columns=reader.fieldnames;rows=list(reader)
    result,backgrounds,ratio=voltage_guess(rows,read_metrics(a.source_load),a.voltage)
    a.output.mkdir(parents=True,exist_ok=False)
    dest=a.output/'nodes.csv'
    with dest.open('x',newline='',encoding='ascii') as stream:
        writer=csv.DictWriter(stream,fieldnames=columns);writer.writeheader();writer.writerows(result)
    report=dict(status='INITIAL_GUESS_ONLY_NOT_A_SOLUTION',intended_initial_heat_factor=0.,
        final_heat_factor_must_be=1.,target_voltage=a.voltage,scale_ratio=ratio,
        source_nodes=str(a.source_nodes.absolute()),source_load=str(a.source_load.absolute()),
        source_nodes_sha256=sha(a.source_nodes),source_load_sha256=sha(a.source_load),
        output_sha256=sha(dest),nodes=len(result),background_densities=backgrounds,
        required_action='NATIVE_REEQUILIBRATION_AT_TARGET_VOLTAGE_THEN_FULL_HEAT')
    (a.output/'initial_guess.json').write_text(json.dumps(report,indent=2,allow_nan=False))
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
