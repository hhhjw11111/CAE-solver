"""Prepare a zero-heat initial guess on the same mesh, never a solved field.

Preserve the source electric potential and geometry. Reset temperature rise,
displacements and log-density perturbations, and reconstruct background carrier
densities. This avoids importing a large heated deformation into the first
zero-heat load of voltage-then-heat continuation. All fields must be re-solved.
"""
from pathlib import Path
import argparse,csv,json,math
from run_hv_sweep import sha

RESET=['ux_m','uy_m','uz_m','theta_K','eta_n','eta_p']
REQUIRED=['node','x_m','y_m','z_m','phi_V','n_m3','p_m3']+RESET

def zero_heat_rows(rows):
    if not rows:raise ValueError('Empty source field')
    backgrounds=[];result=[]
    for index,row in enumerate(rows,1):
        if any(k not in row for k in REQUIRED):raise ValueError('Missing required field')
        if int(row['node'])!=index or not all(math.isfinite(float(row[k])) for k in REQUIRED):
            raise ValueError('Nonfinite or unordered source field')
        n0=float(row['n_m3'])*math.exp(-float(row['eta_n']))
        p0=float(row['p_m3'])*math.exp(-float(row['eta_p']))
        if not all(math.isfinite(v) and v>0 for v in [n0,p0]):raise ValueError('Invalid carrier background')
        if not backgrounds:backgrounds=[n0,p0]
        if not all(math.isclose(a,b,rel_tol=1e-10,abs_tol=0.) for a,b in zip([n0,p0],backgrounds)):
            raise ValueError('This initial-guess tool requires uniform background densities')
        new=dict(row)
        for key in RESET:new[key]='0.00000000000000000e+00'
        new['n_m3']=f'{backgrounds[0]:.17e}';new['p_m3']=f'{backgrounds[1]:.17e}'
        result.append(new)
    return result,backgrounds

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-nodes',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    with args.source_nodes.open(newline='',encoding='utf-8-sig') as stream:
        reader=csv.DictReader(stream);columns=reader.fieldnames;rows=list(reader)
    output,backgrounds=zero_heat_rows(rows)
    args.output.mkdir(parents=True,exist_ok=False)
    dest=args.output/'nodes.csv'
    with dest.open('x',newline='',encoding='ascii') as stream:
        writer=csv.DictWriter(stream,fieldnames=columns);writer.writeheader();writer.writerows(output)
    report=dict(status='INITIAL_GUESS_ONLY_NOT_A_SOLUTION',intended_initial_heat_factor=0.,
        source_nodes=str(args.source_nodes.absolute()),source_sha256=sha(args.source_nodes),
        output_sha256=sha(dest),reset_fields=RESET,preserved_fields=['node','x_m','y_m','z_m','phi_V'],
        background_n_m3=backgrounds[0],background_p_m3=backgrounds[1],nodes=len(output),
        final_heat_factor_must_be=1.)
    (args.output/'initial_guess.json').write_text(json.dumps(report,indent=2,allow_nan=False))
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
