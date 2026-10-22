"""Reproduce saved H1 L0 midplane Joule fields in the requested rainbow palette."""
from pathlib import Path
import sys, json, hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.colors import Normalize
import plot_f1_main_results as fe

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'output/H1_L0_焦耳热_JET_20260909'
RUN = ROOT / 'solver_v2/runs/F1_FAST_HOLE_L0_20260907'
SOURCE = ROOT / 'output/F1_汇报成果_20260907/04_数据表与图源'
MESH = ROOT / 'solver_v2/cases/F1_FAST_HOLE_L0/H1_L0.mesh'
OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({'font.size':10, 'axes.titlesize':12, 'svg.fonttype':'none', 'pdf.fonttype':42})

# Figure contract: locate the two hole-edge hotspots using verified FE fields.
# Quantitative grid; shared linear jet scale, 0--16 x 10^12 W/m^3.
# 3D height encodes Q, not physical deformation. No cross-hole interpolation.
xyz, conn, _ = fe.read_mesh(MESH)
assert len(xyz) == 10144 and len(conn) == 2240
fe.RUNS = RUN
data, audit = {}, {'nodes':len(xyz), 'elements':len(conn), 'models':{}}
for model in ['LQ','N0','NT']:
    case = f'FAST_H1_{model}_L0'
    node = np.genfromtxt(RUN/case/'nodes.csv', delimiter=',', names=True, usecols=range(13))
    np.testing.assert_allclose(np.column_stack([node[c] for c in ['x_m','y_m','z_m']]),xyz,rtol=0,atol=1e-19)
    state = np.column_stack([node[c] for c in ['ux_m','uy_m','uz_m','phi_V','eta_n','eta_p','theta_K']])
    error = fe.check_saved_sample(case, {'model':model}, xyz, conn, state)
    regenerated = fe.plane(xyz,conn,state,model,0.)
    source = SOURCE/f'{case}_midplane.npz'
    with np.load(source) as saved:
        for k in ['xyz','triangles','Q']:
            np.testing.assert_allclose(regenerated[k],saved[k],rtol=1e-12,atol=1e-20)
    data[model] = regenerated
    np.savez_compressed(OUT/f'H1_{model}_L0_midplane.npz',xyz=regenerated['xyz'],triangles=regenerated['triangles'],Q=regenerated['Q'])
    q=regenerated['Q']; i=q.argmax()
    audit['models'][model]={'solver_sample_relative_error':error,'Q_min_W_m3':float(q.min()),'Q_max_W_m3':float(q.max()),'max_xyz_um':(regenerated['xyz'][i]*1e6).tolist(),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest()}

norm=Normalize(0,16)
def tri(d):
    return mtri.Triangulation(d['xyz'][:,0]*1e6,d['xyz'][:,1]*1e6,d['triangles'])

def flat(ax,d,title):
    artist=ax.tripcolor(tri(d),d['Q']*1e-12,cmap='jet',norm=norm,shading='gouraud',rasterized=True)
    angle=np.linspace(0,2*np.pi,300)
    ax.plot(3*np.cos(angle),3*np.sin(angle),color='black',lw=.65)
    ax.plot([-10,10,10,-10,-10],[-10,-10,10,10,-10],color='black',lw=.7)
    ax.set(aspect='equal',xlim=(-10,10),ylim=(-10,10),xticks=[-10,-5,0,5,10],yticks=[-10,-5,0,5,10],xlabel=r'$x$ ($\mu$m)',ylabel=r'$y$ ($\mu$m)',title=title)
    return artist

def save(fig,name):
    for ext in ['png','svg','pdf']:
        fig.savefig(OUT/f'{name}.{ext}',dpi=320,bbox_inches='tight',facecolor='white')
    plt.close(fig)

fig,axs=plt.subplots(1,3,figsize=(12,4.4),layout='constrained')
for ax,(model,d) in zip(axs,data.items()):
    artist=flat(ax,d,'H1-'+model+' | L0')
fig.colorbar(artist,ax=axs,shrink=.76,pad=.025,label=r'$Q_J$ ($10^{12}$ W m$^{-3}$)',ticks=[0,4,8,12,16])
fig.suptitle(r'Joule heat density | $z=0$, $|\Delta V|=3$ V',fontsize=14)
save(fig,'H1_L0_three_models_jet')

fig=plt.figure(figsize=(11.5,5.2))
ax=fig.add_axes([.07,.17,.34,.72])
artist=flat(ax,data['NT'],'(a) Midplane field | H1-NT')
ax=fig.add_axes([.44,.10,.44,.83],projection='3d')
d=data['NT']; q=d['Q']*1e-12
ax.plot_trisurf(tri(d),q,cmap='jet',norm=norm,linewidth=0,antialiased=False,shade=False)
ax.set(xlim=(-10,10),ylim=(-10,10),zlim=(0,16),xticks=[-10,0,10],yticks=[-10,0,10],zticks=[0,8,16],xlabel=r'$x$ ($\mu$m)',ylabel=r'$y$ ($\mu$m)',zlabel=r'$Q_J$ ($10^{12}$ W m$^{-3}$)',title='(b) Field-height view | H1-NT')
ax.set_box_aspect((1,1,.65)); ax.view_init(elev=29,azim=-60)
for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:
    axis.pane.fill=False
    axis._axinfo['grid']['color']=(.7,.7,.7,.35)
cax=fig.add_axes([.94,.24,.015,.53])
fig.colorbar(artist,cax=cax,ticks=[0,4,8,12,16],label=r'$Q_J$ ($10^{12}$ W m$^{-3}$)')
fig.text(.5,.025,'L0: 2,240 HEX20 elements | Height represents heat density, not deformation.',ha='center',fontsize=10)
save(fig,'H1_NT_L0_jet_map_and_height')
audit.update({'color_scale':[0,16e12],'colormap':'jet (requested)','section':'z=0, trace from positive-z element','integrity':'Original HEX20 patch triangles; no smoothing between elements; empty hole; no clipping. 3D vertical coordinate is Q, not displacement.','scope':'Converged L0 solutions; not final mesh convergence acceptance.'})
(OUT/'data_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
print(json.dumps(audit,indent=2))
print(OUT)
