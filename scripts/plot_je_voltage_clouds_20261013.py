"""Compare saved J.E 6/7/8 V FEM fields; no new PDE solve or field smoothing."""
from pathlib import Path
import argparse
import csv
import json
import shutil
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Circle
from matplotlib.ticker import MaxNLocator
import plot_je_8v_clouds_20261013 as fe

QA = Path('C:/Users/Administrator/.codex/skills/nature-figure/scripts')
sys.path.insert(0, str(QA))
from audit_panel_alignment import require_matplotlib_panel_alignment

ROOT = fe.SOLVER / 'runs/JE_6_8_HALF_V_20261013'
OUT = ROOT / 'clouds_6_7_8V_template_20261013'
VOLTAGES = [6, 7, 8]
FACES = {'midplane': 0., 'top': .5e-6, 'bottom': -.5e-6}
# id, field, factor, plane, title, units, palette
SPECS = [
    ('01_joule_mid', 'Q', 1e-12, 'midplane', '焦耳热密度 · J·E', '10¹² W/m³', 'heat'),
    ('02_potential_mid', 'phi', 1., 'midplane', '电势', 'V', 'viridis'),
    ('03_electrons_mid', 'dn_percent', 1., 'midplane', '电子浓度相对变化', '(n/n0 − 1) × 100%', 'signed'),
    ('04_holes_mid', 'dp_percent', 1., 'midplane', '空穴浓度相对变化', '(p/p0 − 1) × 100%', 'signed'),
    ('05_total_stress_mid', 'vm', 1e-6, 'midplane', '耦合总等效应力（含热效应）', 'MPa', 'magma'),
    ('06_temperature_mid', 'theta', 1., 'midplane', '温升', 'K', 'magma'),
    ('07_electrons_bottom', 'dn_percent', 1., 'bottom', '下表面电子浓度相对变化', '(n/n0 − 1) × 100%', 'signed'),
    ('08_holes_top', 'dp_percent', 1., 'top', '上表面空穴浓度相对变化', '(p/p0 − 1) × 100%', 'signed'),
    ('09_thermal_term_mid', 'thermal_xx', 1e-6, 'midplane', '热膨胀本构应力项 −βxxΔT（非独立热应力解）', 'MPa', 'Blues_r'),
]


def prepare(render_only):
    OUT.mkdir(exist_ok=True)
    (OUT/'source_data').mkdir(exist_ok=True)
    (OUT/'qa').mkdir(exist_ok=True)
    contract = {
        'question': 'How do saved J.E-branch spatial fields change from 6 to 8 V?',
        'template_reuse': 'Style only: three horizontal square panels, circular hole; voltage replaces radius.',
        'data': 'Actual saved 6, 7, 8 V states on one r=3 um mesh; no voltage interpolation.',
        'geometry': '20 x 20 x 1 um; exact z=0, +/-0.5 um element faces.',
        'scales': 'Shared linear color scale across the three voltages for each figure; signed changes centered on zero.',
        'mapping': SPECS,
        'interpolation': 'HEX20 evaluation, 16 subdivisions/face; element-local triangulation. No cross-element averaging or cross-hole interpolation.',
        'uncertainty': 'Deterministic single-mesh solutions; no statistical intervals or mesh-convergence claim.',
        'stress': 'VM is full C:strain-e^T E-beta*dT. Direct -beta*dT is separately labeled, not an equilibrated thermal-only solution.',
        'carriers': {'n0_m3': fe.N0, 'p0_m3': fe.P0},
    }
    fe.write_json(OUT/'figure_contract.json', contract)
    mesh = fe.SOLVER/'cases/HV_SURFACE_20261008/XY0_Z2/H1_XY0_Z2.mesh'
    xyz, conn = fe.read_mesh(mesh)
    data, metrics, provenance = {}, {}, {'mesh': {'path': str(mesh), 'sha256': fe.sha256(mesh)}}
    for voltage in VOLTAGES:
        run = ROOT/f'V{voltage*10:03d}'
        met = fe.read_metrics(run/'metrics.txt')
        assert met['voltage_V'] == voltage and met['final_voltage_factor'] == met['final_joule_factor'] == 1
        assert met['heat_source'] == 'E_DOT_TOTAL_CURRENT_DIAGNOSTIC_NOT_THERMODYNAMIC_CLOSURE'
        assert met['case_acceptance'] == 'PASS_NOT_MESH_ACCEPTANCE'
        metrics[voltage] = met
        provenance[str(voltage)] = {name: {'path': str(run/name), 'sha256': fe.sha256(run/name)}
                                    for name in ['nodes.csv', 'metrics.txt', 'field_samples.csv']}
        if not render_only:
            nodes = np.loadtxt(run/'nodes.csv', delimiter=',', skiprows=1)
            assert np.isfinite(nodes).all() and np.max(abs(nodes[:,1:4]-xyz)) < 1e-19
            state = nodes[:,4:11]
            audit = fe.match_native_samples(xyz,conn,state,run,met)
            fe.write_json(OUT/'qa'/f'V{voltage}_native_audit.json', audit)
            for face, z in FACES.items():
                d = fe.plane(xyz,conn,state,z,16)
                np.savez_compressed(OUT/'source_data'/f'V{voltage}_{face}.npz', **d)
        data[voltage] = {face: dict(np.load(OUT/'source_data'/f'V{voltage}_{face}.npz')) for face in FACES}
    provenance['evaluator'] = {'path': str(Path(fe.__file__)), 'sha256': fe.sha256(Path(fe.__file__))}
    fe.write_json(OUT/'provenance.json', provenance)
    return data, metrics


def render(data):
    style = json.loads((fe.SOLVER/'scripts/styles/joule_map_height_frozen.json').read_text(encoding='utf-8-sig'))
    heat = LinearSegmentedColormap.from_list('reference_heat',style['palette_stops'],N=256)
    plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Microsoft YaHei', 'Arial', 'DejaVu Sans'],
                         'font.size': 9, 'axes.labelsize': 9, 'xtick.labelsize': 8, 'ytick.labelsize': 8,
                         'svg.fonttype': 'none', 'pdf.fonttype': 42})
    extrema, scales = [], {}
    for stem,key,factor,face,title,units,palette in SPECS:
        vals = [fe.values(data[v][face],key)*factor for v in VOLTAGES]
        lo, hi = min(float(x.min()) for x in vals), max(float(x.max()) for x in vals)
        if palette == 'signed':
            lo,hi = -max(abs(lo),abs(hi)),max(abs(lo),abs(hi))
        elif key in ['Q','vm','theta'] and lo >= 0:
            lo = 0.
        scales[stem] = [lo,hi]
        cmap = heat if palette == 'heat' else ('RdBu_r' if palette == 'signed' else palette)
        fig = plt.figure(figsize=(10.8,4.4))
        axes = [fig.add_axes([.058+i*.278,.18,.235,.235*10.8/4.4]) for i in range(3)]
        ztext = {'midplane': '中面 z = 0 μm', 'top': '上表面 z = +0.5 μm', 'bottom': '下表面 z = −0.5 μm'}[face]
        fig.text(.45,.957,title,ha='center',va='top',fontsize=12)
        fig.text(.45,.874,f'J·E 分支  |  r = 3 μm  |  {ztext}',ha='center',fontsize=8.5,color='#555555')
        for ax,voltage,values in zip(axes,VOLTAGES,vals):
            d = data[voltage][face]
            xy = d['xyz'][:,:2]*1e6
            tri = mtri.Triangulation(xy[:,0],xy[:,1],d['triangles'])
            artist = ax.tripcolor(tri,values,shading='gouraud',cmap=cmap,norm=Normalize(lo,hi),rasterized=True)
            ax.add_patch(Circle((0,0),3,fill=False,edgecolor='black',linewidth=.65))
            ax.set(xlim=(-10,10),ylim=(-10,10),aspect='equal',xlabel='x (μm)',ylabel='y (μm)',title=f'{voltage} V')
            ax.set_xticks([-10,-5,0,5,10]); ax.set_yticks([-10,-5,0,5,10])
            ax.tick_params(length=2,pad=6)
            for extreme,idx in [('min',np.argmin(values)),('max',np.argmax(values))]:
                extrema.append({'figure':stem,'voltage_V':voltage,'field':key,'plane':face,'extreme':extreme,
                                'value':float(values[idx]),'units':units,'x_um':float(xy[idx,0]),'y_um':float(xy[idx,1])})
        cax = fig.add_axes([.868,.18,.014,.235*10.8/4.4])
        cb = fig.colorbar(artist,cax=cax)
        cb.locator = MaxNLocator(nbins=5); cb.update_ticks()
        cb.ax.tick_params(labelsize=8,length=2)
        cb.set_label(units,fontsize=9,labelpad=7)
        fig.text(.45,.025,'三列共用色标；白色圆孔为非材料区域',ha='center',fontsize=8,color='#555555')
        require_matplotlib_panel_alignment(fig,axes=axes,panel_ids=['6V','7V','8V'],row_groups=[['6V','7V','8V']],
                                           json_out=OUT/'qa'/f'{stem}_alignment.json')
        fig.savefig(OUT/f'{stem}.png',dpi=330)
        fig.savefig(OUT/f'{stem}.pdf',dpi=330)
        fig.savefig(OUT/f'{stem}.svg',dpi=330)
        plt.close(fig)
        print('EXPORTED='+stem,flush=True)
    with (OUT/'plane_extrema.csv').open('w',newline='',encoding='utf-8-sig') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(extrema[0])); writer.writeheader(); writer.writerows(extrema)
    fe.write_json(OUT/'color_scales.json',scales)


def summarize(metrics):
    rows = []
    for v,m in metrics.items():
        rows.append({'voltage_V': v, 'joule_mW': m['joule_W']*1e3, 'theta_max_K': m['theta_max_K'],
                     'VM_max_MPa': m['stress_von_mises_max_Pa']*1e-6,
                     'terminal_power_mW': m['terminal_power_W']*1e3,
                     'heat_balance_percent': m['heat_balance_relative']*100,
                     'JE_terminal_gap_percent': abs(m['joule_W']-m['terminal_power_W'])/abs(m['joule_W'])*100})
    with (OUT/'global_metrics.csv').open('w',newline='',encoding='utf-8-sig') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    fe.write_json(OUT/'global_metrics.json',rows)
    shutil.copy2(__file__,OUT/Path(__file__).name)
    shutil.copy2(fe.__file__,OUT/Path(fe.__file__).name)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--render-only',action='store_true')
    args=parser.parse_args()
    data,metrics=prepare(args.render_only)
    render(data)
    summarize(metrics)
