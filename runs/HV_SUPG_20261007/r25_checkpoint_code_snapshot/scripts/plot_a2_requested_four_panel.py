"""Requested A2 centerline/wall potential and Joule comparison, frozen data."""
from pathlib import Path
import json
import csv
import shutil
import matplotlib.pyplot as plt
import plot_a2_skill_delivery as base

OUT=base.ROOT/'output/A2_解析解FEM_四宫格_黑色_20260909'
OUT.mkdir(parents=True,exist_ok=True)
base.OUT=OUT
base.REF_BLUE='#000000'
original_save=base.save_all

def save_requested(fig,stem):
    fig.set_size_inches(183/25.4,143/25.4)
    fig.subplots_adjust(left=.095,right=.98,bottom=.105,top=.95,wspace=.32,hspace=.40)
    titles=['(a) 中心线电势','(b) 孔壁电势','(c) 中心线焦耳热源','(d) 孔壁焦耳热源']
    for ax,title in zip(fig.axes,titles):
        ax.set_title(title,fontsize=9,pad=7)
    for ax in (fig.axes[0],fig.axes[2]):
        ax.axvspan(.35,.65,color='#f3f3f3',zorder=-10)
        ax.text(.5,.10,'孔洞',transform=ax.transAxes,ha='center',va='center',fontsize=7,color='.48')
        ax.set_xlabel(r'中心线位置 $(x+L/2)/L$')
    fig.text(.5,.014,r'$L=20\,\mu\mathrm{m}$，$a=3\,\mu\mathrm{m}$，$E_0=10^4\,\mathrm{V/m}$；$Q_\infty=\sigma E_0^2$',
             ha='center',fontsize=7)
    original_save(fig,'A2_analytic_FEM_four_panel')

base.save_all=save_requested
base.plot_reference_style_four_comparisons()
path=OUT/'A2_reference_four_panel_contract.json'
contract=json.loads(path.read_text(encoding='utf-8'))
contract.pop('reference_image',None)
contract['layout'].update({'figure_size_mm':[183,143],'panel_titles':True})
contract['layout'].update({'analytical':'black solid line with filled circular markers',
                           'FEM':'unconnected black open upward triangles'})
contract['geometry']={'plate_um':[20,20,1],'radius_um':3,'E0_V_m':10000,
                      'mesh':'A2 L2; 128 HEX20 elements, 1056 nodes',
                      'centerline':'y=z=0; normalized coordinate (x+L/2)/L, L=20 um',
                      'wall':'r=a; psi measured counterclockwise from +x'}
contract['boundaries']={'outer_four_sides':'phi=exact analytical potential',
                        'hole_and_top_bottom':'natural zero normal electric/carrier flux',
                        'all_nodes':'u=0, temperature rise=0'}
metrics=base.read_metrics(base.RUN/'A2_stdout.txt')
contract['independent_error_metrics']={k:v for k,v in metrics.items() if 'relative_L2' in k or 'wall_normal_E_RMS' in k}
contract['sources_sha256']={p.name:base.sha256(p) for p in [base.RUN/'A2_L2_nodes.csv',base.RUN/'A2_L2_wall.csv',base.RUN/'A2_stdout.txt']}
path.write_text(json.dumps(contract,ensure_ascii=False,indent=2),encoding='utf-8')
with (OUT/'A2_error_convergence.csv').open('w',newline='',encoding='utf-8-sig') as f:
    w=csv.writer(f);w.writerow(['level','phi_relative_L2','E_relative_L2','wall_normal_E_RMS_over_E0'])
    for level in ['L0','L1','L2']:
        w.writerow([level]+[metrics[f'a2_{level}_{key}'] for key in ['phi_relative_L2','E_relative_L2','wall_normal_E_RMS_over_E0']])
for name in ['plot_a2_requested_four_panel.py','plot_a2_skill_delivery.py']:
    shutil.copyfile(Path(__file__).parent/name,OUT/name)
for ext in ['svg','pdf','png','tiff']:
    p=OUT/f'A2_analytic_FEM_four_panel.{ext}'
    assert p.exists() and p.stat().st_size>1000
assert '<text' in (OUT/'A2_analytic_FEM_four_panel.svg').read_text(encoding='utf-8')
print(json.dumps(contract['displayed_comparison_maximum_absolute_differences'],indent=2))
print('Exports verified; SVG text is editable. No solver rerun or data rescaling.')
