"""Additional exact FE midplane plots and provenance checks; no PDE solves."""
import postprocess_3v as pp
from pathlib import Path
import json, numpy as np, matplotlib.pyplot as plt, matplotlib.tri as tri

def main():
    audits=[];planes={}
    for case in pp.CASES:
        root=case['root'];exe=json.loads((root/'execution.json').read_text(encoding='utf-8-sig'));m=pp.metrics(root/'metrics.txt')
        limits={'maximum_stage_balance':1e-8,'maximum_accepted_residual':1e-9,'maximum_accepted_update':1e-9,'heat_balance_relative':.01,'electron_balance_relative':.01,'joule_quadrature_4_5_relative':.001}
        if case['half']:
            limits.update({k:.01 for k in ['hole_balance_diagnostic','physical_exterior_heat_balance_relative','physical_exterior_electron_balance_relative','physical_exterior_hole_balance_diagnostic','symmetry_heat_flux_relative','symmetry_electron_flux_relative','symmetry_hole_flux_diagnostic']})
        gates={k:bool(np.isfinite(m.get(k,np.nan)) and 0<=m[k]<v) for k,v in limits.items()}
        source_manifest=pp.F1/'evidence/build_manifest.json' if case['family']=='F1_L0' else case['source'].parent/'manifest.json'
        manifest=json.loads(source_manifest.read_text(encoding='utf-8-sig'));items=manifest['files'] if isinstance(manifest,dict) else manifest
        expected={i['path'].replace('\\','/'):i['sha256'].lower() for i in items}
        source_checks=[]
        for p in sorted(case['source'].rglob('*.f90')):
            key='src/'+p.relative_to(case['source']).as_posix()
            source_checks.append(dict(path=str(p.resolve()),matches_archive_manifest=pp.sha(p)==expected.get(key)))
        input_check=pp.sha(root/'input.nml')==exe['input_sha256'] if 'input_sha256' in exe else pp.sha(root/'input.nml')==pp.sha(pp.F1/f"evidence/inputs/{case['name']}.nml")
        record=dict(case=case['name'],execution_exit_code=exe['exit_code'],metric_gate_checks=gates,input_hash_matches_record_or_snapshot=input_check,source_files_checked=len(source_checks),source_hashes_match_archive=all(i['matches_archive_manifest'] for i in source_checks),source_checks=source_checks)
        audits.append(record)
        if not (exe['exit_code']==0 and all(gates.values()) and input_check and record['source_hashes_match_archive']):raise ValueError(record)
        if case['name'] not in ['FAST_H0_NT_L0','FAST_H1_NT_L0']:continue
        xyz,conn,_=pp.read_mesh(case['mesh']);d=np.genfromtxt(root/'nodes.csv',delimiter=',',names=True);state=np.column_stack([d[k] for k in ['ux_m','uy_m','uz_m','phi_V','eta_n','eta_p','theta_K']]);del d
        t=np.linspace(-1,1,7);points=np.array([[a,b,-1] for b in t for a in t]);rows=[];triangles=[];offset=0
        for eid,nn in enumerate(conn):
            co=xyz[nn]
            if abs(co[:,2].min())>1e-18:continue
            f=pp.fields(co[None,:],state[nn][None,:],points,'NT')
            for i in range(len(points)):
                rows.append(dict(x_m=float(f['xyz'][0,i,0]),y_m=float(f['xyz'][0,i,1]),z_m=float(f['xyz'][0,i,2]),element=eid+1,phi_V=float(f['val'][0,i,3]),E_V_m=float(np.linalg.norm(f['E'][0,i])),J_A_m2=float(np.linalg.norm(f['J'][0,i])),theta_K=float(f['val'][0,i,6]),Q_W_m3=float(f['pd'][0,i])))
            for j in range(6):
                for i in range(6):
                    a=offset+j*7+i;triangles.extend([[a,a+1,a+8],[a,a+8,a+7]])
            offset+=49
        pp.csvwrite(pp.OUT/f"data/{case['name']}_midplane.csv",rows)
        np.savetxt(pp.OUT/f"data/{case['name']}_midplane_triangles.csv",np.array(triangles),delimiter=',',fmt='%d',header='vertex0,vertex1,vertex2',comments='')
        planes[case['geom']]=(rows,np.array(triangles))
    pp.dump(pp.OUT/'audit/postprocess/provenance_gate_recheck.json',audits)
    fig,axs=plt.subplots(2,3,figsize=(7.2,4.6),layout='constrained')
    for col,(key,factor,label) in enumerate([('E_V_m',1e-5,'E (10^5 V/m)'),('J_A_m2',1e-7,'J (10^7 A/m²)'),('theta_K',1,'Temperature rise (K)')]):
        vmax=max(r[key]*factor for rows,_ in planes.values() for r in rows)
        for row,geom in enumerate(['H0','H1']):
            rows,tr=planes[geom];x=np.array([r['x_m'] for r in rows])*1e6;y=np.array([r['y_m'] for r in rows])*1e6;v=np.array([r[key] for r in rows])*factor
            art=axs[row,col].tripcolor(tri.Triangulation(x,y,tr),v,cmap='viridis' if col<2 else 'magma',vmin=0,vmax=vmax,shading='gouraud',rasterized=True)
            axs[row,col].set(aspect='equal',xlabel='x (µm)',ylabel='y (µm)',title=geom+' NT');fig.colorbar(art,ax=axs[row,col],shrink=.7,label=label)
    fig.suptitle('3 V F1 L0 exploratory, midplane z=0+; mesh pending\nExact FE samples within original element patches; common column scales',fontsize=9);pp.savefig(fig,'fig04_midplane_fields')
    pp.dump(pp.OUT/'audit/postprocess/supplemental_status.json',dict(provenance_rechecks='ALL_SELECTED_CASES_PASS_ARCHIVE_IDENTITY',plots='F1_NT_midplane_only',field_interpolation='HEX20 eta before exp; displayed linear subtriangles confined to each original patch'))
if __name__=='__main__':main()
