"""Make an explicitly non-equilibrium F1 seed for an independent restart check.

The smooth thermal perturbation vanishes on the outer x/y walls and respects
y-reflection symmetry. This changes an initial guess, never a physical input
or a solved result. A fresh native solve and all physical audits are required.
"""
import argparse
import csv
import json
import math
from pathlib import Path
from run_hv_sweep import sha


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source-nodes',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--temperature-amplitude-k',type=float,default=.01)
    a=p.parse_args()
    amplitude=a.temperature_amplitude_k
    if not math.isfinite(amplitude) or amplitude<=0:
        raise ValueError('A finite positive perturbation is required')
    with a.source_nodes.open(newline='',encoding='utf-8-sig') as f:
        reader=csv.DictReader(f);columns=reader.fieldnames;rows=list(reader)
    if not rows or any(not math.isfinite(float(value)) for row in rows for value in row.values()):
        raise ValueError('Empty or nonfinite initial field')
    xs=[float(row['x_m']) for row in rows];ys=[float(row['y_m']) for row in rows]
    width=1e-5
    if abs(min(xs)+width)>1e-15 or abs(max(xs)-width)>1e-15 or max(map(abs,ys))>width+1e-15:
        raise ValueError('This perturbation is restricted to the frozen F1 domain')
    changes=[]
    for row,x,y in zip(rows,xs,ys):
        delta=amplitude*max(0.,1-(x/width)**2)*max(0.,1-(y/width)**2)
        changes.append(delta)
        temperature=float(row['theta_K'])+delta
        if not math.isfinite(temperature):
            raise ValueError('Nonfinite perturbed temperature')
        row['theta_K']=f"{temperature:.17e}"
    if max(changes)==0:
        raise ValueError('No interior node was perturbed')
    a.output.mkdir(parents=True,exist_ok=False)
    output=a.output/'nodes.csv'
    with output.open('x',newline='',encoding='ascii') as f:
        writer=csv.DictWriter(f,fieldnames=columns);writer.writeheader();writer.writerows(rows)
    report=dict(status='PERTURBED_INITIAL_GUESS_ONLY_REQUIRES_NATIVE_RESOLVE',
                source_nodes=str(a.source_nodes.absolute()),source_nodes_sha256=sha(a.source_nodes),
                nodes_sha256=sha(output),temperature_amplitude_K=amplitude,
                maximum_actual_temperature_change_K=max(changes),
                formula='delta_theta = amplitude * max(0,1-(x/1e-5)^2) * max(0,1-(y/1e-5)^2)')
    (a.output/'perturbation.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report))


if __name__=='__main__':main()
