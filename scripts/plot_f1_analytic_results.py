"""Plot only saved F1 numerical evidence and the frozen analytical expressions."""
from pathlib import Path
import json,re
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'runs/F1_20260906'
OUT=ROOT.parent/'output/figures/F1_20260906'
OUT.mkdir(parents=True,exist_ok=True)
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],
 'font.size':8,'axes.titlesize':9,'axes.labelsize':8,'legend.fontsize':7,
 'svg.fonttype':'none','pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,
 'axes.linewidth':.7,'legend.frameon':False,'lines.linewidth':1.1})
BLUE='#236B8E'; ORANGE='#C96B35'; GREY='#54636E'

def metrics(path):
    pairs={}
    for line in path.read_text(encoding='utf-8-sig').splitlines():
        if '=' in line:
            k,v=line.split('=',1)
            try:pairs[k]=float(v)
            except ValueError: pass
    return pairs

def save(fig,name):
    for ext in ['svg','pdf','png']:
        fig.savefig(OUT/(name+'.'+ext),dpi=260,bbox_inches='tight',facecolor='white')
    plt.close(fig)

def a1():
    d=np.genfromtxt(DATA/'A1/A1_ne16_nodes.csv',delimiter=',',names=True)
    _,ix=np.unique(np.round(d['x_m'],18),return_index=True);d=d[ix]
    x=np.linspace(0,1e-6,501); L=1e-6;vt=1.380649e-23*300/1.602176634e-19
    eps=1e-8+1/1e11;g2=1.602176634e-19*(1e21+8e20)/vt/eps
    s=(1e-8+1e-5)/(eps*L);g=np.sqrt(g2)
    phi=s/g2*(np.cosh(g*(x-L/2))/np.cosh(g*L/2)-1)
    u=1e-5*x*x/(2*L)-phi/1e11
    fig,axes=plt.subplots(2,3,figsize=(7.2,4.5),layout='constrained')
    entries=[('phi_V',phi*1e3,1e3,'Potential (mV)'),
             ('eta_n',phi/vt*100,100,'Electron change (%)'),
             ('eta_p',-phi/vt*100,100,'Hole change (%)'),
             ('theta_K',x/L,1,'Temperature rise (K)'),
             ('u_m',u*1e12,1e12,'Axial displacement (pm)')]
    for ax,(col,exact,scale,ylabel) in zip(axes.flat,entries):
        ax.plot(x*1e6,exact,c=GREY,label='Analytical')
        ax.plot(d['x_m'][::2]*1e6,d[col][::2]*scale,'o',ms=3.5,mfc='white',mec=BLUE,label='HEX20 (16 elements)')
        ax.set(xlabel=r'$x$ ($\mu$m)',ylabel=ylabel);ax.tick_params(direction='out')
    axes.flat[0].legend(loc='best')
    m=metrics(DATA/'A1/A1_export_stdout.txt');ne=np.array([1,2,4,8,16])
    ax=axes.flat[5]
    for norm,marker,color in [('l2','o',BLUE),('h1','s',ORANGE)]:
        vals=[m[f'a1_ne{n}_phi_{norm}'] for n in ne]
        ax.loglog(ne,vals,marker=marker,c=color,label={'l2':r'$L^2$ error','h1':r'$H^1$ seminorm'}[norm])
    ax.set(xlabel='Elements along x',ylabel='Relative potential error',xticks=ne)
    ax.set_xticklabels(ne);ax.legend(loc='best')
    for label,ax in zip('abcdef',axes.flat):ax.text(-.17,1.05,label,transform=ax.transAxes,fontweight='bold')
    save(fig,'A1_joint_analytic')
    rows=[[n,m[f'a1_ne{n}_phi_l2'],m[f'a1_ne{n}_phi_h1'],m[f'a1_ne{n}_eta_n_l2'],m[f'a1_ne{n}_eta_p_l2']] for n in ne]
    np.savetxt(OUT/'A1_error_source.csv',rows,delimiter=',',header='elements,phi_L2,phi_H1,eta_n_L2,eta_p_L2',comments='')

def a2():
    d=np.genfromtxt(DATA/'A2/A2_L2_wall.csv',delimiter=',',names=True)
    angle=np.degrees(d['angle_rad'])%360;ix=np.argsort(angle);angle=angle[ix];d=d[ix]
    exactang=np.linspace(0,360,721)
    fig,grid=plt.subplots(2,2,figsize=(7.2,4.7),layout='constrained');axes=grid.ravel()
    nodes=np.genfromtxt(DATA/'A2/A2_L2_nodes.csv',delimiter=',',names=True)
    nodes=nodes[(np.abs(nodes['y_m'])<1e-15)&(np.abs(nodes['z_m'])<1e-15)]
    nodes.sort(order='x_m')
    for sign in [-1,1]:
        xx=np.linspace(3e-6,10e-6,200)*sign
        axes[0].plot(xx*1e6,-1e4*xx*(1+(3e-6/xx)**2),c=GREY,label='Analytical' if sign==1 else None)
    axes[0].plot(nodes['x_m']*1e6,nodes['phi_V'],'o',ms=3.5,mfc='white',mec=BLUE,label='HEX20 (L2)')
    axes[0].axvspan(-3,3,color='#E8EDF0',zorder=-2)
    axes[0].text(0,0,'Hole',ha='center',va='center',color=GREY)
    axes[0].set(xlabel=r'$x$ ($\mu$m), $y=z=0$',ylabel='Potential (V)');axes[0].legend(loc='upper right')
    for ax,col,exact,ylabel in zip(axes[1:3],['Et_over_E0','Q_over_Qinf'],
        [2*np.abs(np.sin(np.radians(exactang))),4*np.sin(np.radians(exactang))**2],
        [r'$|E_t|/E_0$',r'$Q_J/Q_\infty$']):
        ax.plot(exactang,exact,c=GREY,label='Analytical')
        ax.plot(angle[::10],d[col][::10],'o',ms=3,mfc='white',mec=BLUE,label='HEX20 (L2)')
        ax.set(xlabel='Hole-wall angle (deg)',ylabel=ylabel,xticks=[0,90,180,270,360])
    m=metrics(DATA/'A2/A2_export_stdout.txt')
    for col,label,color,marker in [('phi_relative_L2',r'$\phi$: relative $L^2$',BLUE,'o'),
         ('E_relative_L2',r'$E$: relative $L^2$',ORANGE,'s'),
         ('wall_normal_E_RMS_over_E0',r'$E_n$: wall RMS / $E_0$',GREY,'^')]:
        vals=[m[f'a2_L{i}_{col}'] for i in range(3)]
        axes[3].semilogy([0,1,2],vals,marker=marker,c=color,label=label)
        short={'phi_relative_L2':r'$\phi$: $L^2$','E_relative_L2':r'$E$: $L^2$',
               'wall_normal_E_RMS_over_E0':r'$E_n$: RMS'}[col]
        axes[3].text(2.07,vals[-1],short,va='center',color=color,fontsize=7)
    axes[3].set(xlabel='Mesh level',ylabel='Normalized error',xticks=[0,1,2],xlim=(-.1,2.7))
    for label,ax in zip('abcd',axes):ax.text(-.14,1.06,label,transform=ax.transAxes,fontweight='bold')
    save(fig,'A2_hole_analytic')
    names=['phi_relative_L2','E_relative_L2','wall_normal_E_RMS_over_E0','wall_Et_peak_over_E0','wall_Q_peak_over_Qinf']
    np.savetxt(OUT/'A2_error_source.csv',[[i]+[m[f'a2_L{i}_{n}'] for n in names] for i in range(3)],
        delimiter=',',header='level,'+','.join(names),comments='')

if __name__=='__main__':
    contract={'conclusion':'A1 and A2 meet the frozen analytical-field, order, and boundary-flux criteria; this does not validate all ZnO physics.',
      'archetype':'quantitative grid','backend':'saved Python preference',
      'evidence':'Saved 5-level A1 and 3-level A2 FEM outputs; closed-form curves; no generated numerical data.',
      'export':'183 mm nominal width; editable SVG/PDF text; PNG review; error CSVs',
      'statistics':'Deterministic discretization study, no biological replicates, error bars, or hypothesis tests.',
      'risks':'Nodal overlap alone is insufficient; independent volume errors and wall normal-flux residual are shown.',
      'captions':{
        'A1':'a-e: numerical nodal values (every second x node on the 16-element mesh) and analytical fields. f: independently integrated relative errors. Carrier errors match the potential errors to floating-point precision. Final L2 and H1 rates: 2.990 and 1.991.',
        'A2':'a: potential on y=z=0 with the hole excluded. b-c: L2 wall samples (every tenth record) versus the insulated-circle solution. d: independent volume and boundary errors on all three meshes. Peak factors 1.999857 and 3.999428. Synthetic isothermal test, not a ZnO heating prediction.'}}
    (OUT/'figure_contract.json').write_text(json.dumps(contract,indent=2),encoding='utf-8')
    a1();a2();print(OUT)
