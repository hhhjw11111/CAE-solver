"""Map a solved field to another hole radius ONLY as an initial guess.

Each target point (rho,theta,z) is drawn from the source at the radially
stretched rho_s = a_s+(rho-a_t)(R-a_s)/(R-a_t), R(theta) the ray length to the
square boundary. Hole wall maps to hole wall and the outer boundary (electrodes
included) to itself. The target run must re-solve its own equations; nothing
here is a result for the target radius.
"""
from pathlib import Path
import argparse,csv,json
import numpy as np
from audit_hex20_boundary_geometry import read_mesh
from transfer_hv_initial_guess import FIELDS,transfer
from run_hv_sweep import sha

HALF_WIDTH=1e-5


def stretched_points(target,target_radius,source_radius):
    x,y,z=target[:,0],target[:,1],target[:,2]
    rho=np.hypot(x,y);theta=np.arctan2(y,x)
    ray=HALF_WIDTH/np.maximum(abs(np.cos(theta)),abs(np.sin(theta)))
    if np.any(rho<target_radius*(1-1e-9)) or np.any(rho>ray*(1+1e-9)):
        raise ValueError('Target node outside the plate/hole annulus')
    rho_s=source_radius+(rho-target_radius)*(ray-source_radius)/(ray-target_radius)
    return np.column_stack([rho_s*np.cos(theta),rho_s*np.sin(theta),z])


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source-mesh',type=Path,required=True);p.add_argument('--source-nodes',type=Path,required=True)
    p.add_argument('--source-radius',type=float,required=True)
    p.add_argument('--target-mesh',type=Path,required=True);p.add_argument('--target-radius',type=float,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    coords,conn,_=read_mesh(a.source_mesh);target,_,_=read_mesh(a.target_mesh)
    data=np.genfromtxt(a.source_nodes,delimiter=',',names=True)
    if len(data)!=len(coords) or not np.allclose(np.column_stack([data[k] for k in ['x_m','y_m','z_m']]),coords,atol=1e-16,rtol=0):
        raise ValueError('Source mesh/field mismatch')
    points=stretched_points(target,a.target_radius,a.source_radius)
    vals=np.column_stack([data[k] for k in FIELDS]);values,report=transfer(coords,conn,vals,points)
    n0=float(np.median(data['n_m3']/np.exp(data['eta_n'])))
    p0=float(np.median(data['p_m3']/np.exp(data['eta_p'])))
    a.output.mkdir(parents=True,exist_ok=False)
    with (a.output/'nodes.csv').open('x',newline='',encoding='ascii') as f:
        w=csv.writer(f);w.writerow(['node','x_m','y_m','z_m']+FIELDS+['n_m3','p_m3'])
        for i,(xyz,val) in enumerate(zip(target,values),1):
            w.writerow([i]+list(xyz)+list(val)+[n0*np.exp(val[4]),p0*np.exp(val[5])])
    report.update(status='INITIAL_GUESS_ONLY_NOT_A_SOLVED_FIELD',mapping='RADIAL_STRETCH_HOLE_TO_SQUARE_BOUNDARY',
        source_radius_m=a.source_radius,target_radius_m=a.target_radius,
        maximum_point_shift_m=float(np.max(np.linalg.norm(points-target,axis=1))),
        source_mesh_sha256=sha(a.source_mesh),source_nodes_sha256=sha(a.source_nodes),
        target_mesh_sha256=sha(a.target_mesh),script_sha256=sha(Path(__file__)))
    (a.output/'transfer_report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)


if __name__=='__main__':main()
