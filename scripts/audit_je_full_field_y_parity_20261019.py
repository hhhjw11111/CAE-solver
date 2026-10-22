"""Read-only full-field parity prerequisite; never a half-domain certificate."""
from pathlib import Path
import argparse
import json
import numpy as np
import run_je_reconstruction_mesh_20261018 as runner
from audit_hex20_boundary_geometry import read_mesh
from generate_hv_half_y import verify_reflection, rigid_gauge, U
from audit_hv_half_y import LIMIT


def audit(run):
    manifest=runner.accepted(run)
    mesh=Path(manifest['mesh'])
    xyz,conn,materials=read_mesh(mesh)
    geometry=verify_reflection(xyz,conn,materials)
    keys=np.rint(xyz/1e-16).astype(np.int64)
    lookup={tuple(k):i for i,k in enumerate(keys)}
    reflected=keys.copy();reflected[:,1]*=-1
    ids=np.array([lookup[tuple(k)] for k in reflected])
    assert np.array_equal(ids[ids],np.arange(len(xyz)))
    nodes=np.genfromtxt(run/'nodes.csv',delimiter=',',names=True)
    assert np.array_equal(nodes['node'],np.arange(1,len(xyz)+1))
    assert np.max(abs(np.column_stack([nodes[k] for k in ('x_m','y_m','z_m')])-xyz))<1e-19
    u,gauge=rigid_gauge(xyz,np.column_stack([nodes[k] for k in U]))
    _,m=runner.classify(run,runner.read(run/'execution.json')['exit_code'])
    columns=[u[:,j] for j in range(3)]+[nodes[k] for k in ('phi_V','eta_n','eta_p','theta_K')]
    scales=[1e-10]*3+[max(1.,m['voltage_V']),.1,.1,35.]
    parity=[1,-1,1,1,1,1,1]
    fields={}
    for name,values,scale,sign in zip(U+['phi_V','eta_n','eta_p','theta_K'],columns,scales,parity):
        delta=abs(values-sign*values[ids])/scale
        index=int(np.argmax(delta))
        fields[name]=dict(maximum_scaled_difference=float(delta[index]),scale=scale,
            limit=LIMIT,passed=bool(delta[index]<LIMIT),sign_under_y_reflection=sign,
            worst_node=int(index+1),mirror_node=int(ids[index]+1))
    paths=[run/n for n in ('nodes.csv','input.nml','manifest.json','execution.json','metrics.txt','field_audit.json')]
    paths += [mesh,Path(__file__)]
    return dict(run=str(run),voltage_V=m['voltage_V'],geometry=geometry,gauge=gauge,fields=fields,
        full_field_parity_passed=all(v['passed'] for v in fields.values()),
        evidence_hashes={str(p):runner.sha(p) for p in paths})


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runs',nargs='+',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    results=[audit(p.resolve()) for p in args.runs]
    args.output.mkdir(parents=True,exist_ok=False)
    runner.save(args.output/'report.json',dict(
        scope='FULL_FIELD_PARITY_PREREQUISITE_ONLY',results=results,
        all_parity_passed=all(r['full_field_parity_passed'] for r in results),
        half_solver_started=False,half_equivalence_certified=False,production_acceptance=False,
        limitation='Does not test half-domain solve, perturbed recovery, vector/tensor flux parity, artificial cut flux or spatial/SUPG convergence. Half fallback remains conditional on actual full-domain resource insufficiency.'))
    print(json.dumps({r['run']:r['full_field_parity_passed'] for r in results}))


if __name__=='__main__':main()
