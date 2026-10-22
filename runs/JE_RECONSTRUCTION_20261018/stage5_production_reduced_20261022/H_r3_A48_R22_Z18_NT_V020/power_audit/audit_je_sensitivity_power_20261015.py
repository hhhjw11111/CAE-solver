"""Apply the previously independently checked power audit to each mesh and SUPG scale.

The old audit implements scale=1. The SUPG residual is linear in supg_scale
for a fixed saved state, so multiply only S; physical J, JE and flux unchanged.
"""
from pathlib import Path
import argparse,json,shutil
import numpy as np
import diagnose_je_power_gap_20261015 as audit


def run(directory,orders=(3,5)):
    directory=Path(directory).resolve()
    manifest=json.loads((directory/'manifest.json').read_text())
    summary=json.loads((directory/'summary.json').read_text())
    assert summary['state']=='SOLVE_FLUX_PASS_MESH_PENDING',summary
    scale=manifest['supg_scale'];mesh=Path(manifest['mesh'])
    assert audit.fe.sha256(mesh)==manifest['mesh_sha256']
    assert audit.fe.sha256(Path(manifest['executable']))==manifest['executable_sha256']
    met=audit.fe.read_metrics(directory/'metrics.txt')
    assert met['supg_scale']==scale
    assert met['heat_source']=='E_DOT_TOTAL_CURRENT_DIAGNOSTIC_NOT_THERMODYNAMIC_CLOSURE'
    out=directory/'power_audit';out.mkdir(exist_ok=False)
    original=audit.sample
    def scaled_sample(*args,**kwargs):
        result=original(*args,**kwargs)
        result['S']*=scale
        return result
    audit.sample=scaled_sample
    try:
        xyz,conn=audit.fe.read_mesh(mesh)
        nodes=np.loadtxt(directory/'nodes.csv',delimiter=',',skiprows=1)
        assert np.max(abs(nodes[:,1:4]-xyz))<1e-19
        state=nodes[:,4:11];_,lookup=audit.topology(xyz,conn)
        provenance={'supg_scale':scale,'mesh_sha256':audit.fe.sha256(mesh),
                    'nodes_sha256':audit.fe.sha256(directory/'nodes.csv'),
                    'verification':audit.basic_checks(),
                    'scripts':{str(p):audit.fe.sha256(p) for p in [Path(__file__),Path(audit.__file__),Path(audit.fe.__file__)]}}
        audit.fe.write_json(out/'provenance.json',provenance)
        results={order:audit.audit_order(xyz,conn,state,lookup,order,out) for order in orders}
        if 5 in results:
            for k,mk in [('JE_W','joule_W'),('all_exterior_power_W','terminal_power_W')]:
                assert abs(results[5][k]-met[mk])/abs(met[mk])<1e-8,(k,results[5][k],met[mk])
        if 3 in results and 5 in results:
            a,b=results[3],results[5]
            parts={'quadrature_W':b['JE_W']-a['JE_W'],
                   'SUPG_weighted_W':a['phi_weighted_SUPG_W'],
                   'reaction_minus_physical_boundary_W':a['reaction_terminal_power_W']-b['all_exterior_power_W'],
                   'free_residual_correction_W':-a['phi_weighted_free_residual_W']}
            report={'gap_W':b['gap_W'],'gap_relative':b['gap_relative'],
                    'electrode_gap_relative':(b['JE_W']-b['contact_only_power_W'])/b['JE_W'],
                    'parts':parts,'decomposition_remainder_W':b['gap_W']-sum(parts.values()),
                    'SUPG_fraction_of_gap':parts['SUPG_weighted_W']/b['gap_W'],
                    'free_residual_max_A':a['free_residual_max_A']}
            audit.fe.write_json(out/'decomposition.json',report)
            print(json.dumps({'run':directory.name,**report}),flush=True)
        shutil.copy2(__file__,out/Path(__file__).name)
    finally:
        audit.sample=original


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--runs',nargs='+',type=Path,required=True)
    a=p.parse_args()
    for directory in a.runs:run(directory)
