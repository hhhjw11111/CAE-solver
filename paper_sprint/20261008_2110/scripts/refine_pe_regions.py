"""Correct empty centroid-defined contact sets without redoing integration.
Contacts are the first complete HEX20 layer touching x=+/-10 um.
"""
import postprocess_3v as p
import numpy as np, csv, json

def main():
    regions=[]
    for case in p.CASES:
        fn=p.OUT/f"data/{case['name']}_element_peclet.csv"
        with fn.open(encoding='utf-8-sig',newline='') as f:rows=list(csv.DictReader(f))
        xyz,conn,_=p.read_mesh(case['mesh']);contact=np.max(abs(xyz[conn,0]),axis=1)>=1e-5-1e-16
        for i,row in enumerate(rows):
            x,y,z=map(float,[row['x_m'],row['y_m'],row['z_m']]);r='bulk'
            if np.hypot(x,y)<=3.5e-6:r='hole_neighbourhood' if case['geom']=='H1' else 'same_radius_control'
            if abs(z)>=.45e-6:r='surface_layers'
            if contact[i]:r='contact'
            row['region']=r
        p.csvwrite(fn,rows);pe=np.array([float(r['Pe_model_ratio']) for r in rows]);vol=np.array([float(r['volume_m3']) for r in rows]);reg=np.array([r['region'] for r in rows])
        for group in ['all',*sorted(set(reg))]:
            ix=np.ones(len(rows),bool) if group=='all' else reg==group
            regions.append(dict(case=case['name'],region=group,count=int(ix.sum()),volume_m3=float(vol[ix].sum()),Pe_median=float(np.median(pe[ix])),Pe_p95=float(np.quantile(pe[ix],.95)),Pe_max=float(pe[ix].max()),count_fraction_Pe_gt_1=float(np.mean(pe[ix]>1)),volume_fraction_Pe_gt_1=float(vol[ix&(pe>1)].sum()/vol[ix].sum()),count_fraction_Pe_gt_2=float(np.mean(pe[ix]>2)),volume_fraction_Pe_gt_2=float(vol[ix&(pe>2)].sum()/vol[ix].sum())))
    p.csvwrite(p.OUT/'data/peclet_regions.csv',regions)
    results=json.loads((p.OUT/'data/metrics_summary.json').read_text(encoding='utf-8'))['results'];p.figures(results,regions)
    p.dump(p.OUT/'audit/postprocess/region_definition_correction.json',dict(reason='Centroid |x|>=9um produced empty contact sets on F1 and R23 meshes; revised contact definition is boundary-touching first element layer, before interpreting regional conclusions.',final_definition='Contact first: element has at least one x=+/-10um node. Surface: remaining elements |z_centroid|>=0.45um. Near-hole: remaining r_centroid<=3.5um. Bulk: remainder.',physics_and_global_metrics_unchanged=True))
if __name__=='__main__':main()
