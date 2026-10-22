"""Radius field maps (a=2/3/4 um, NT, 6 and 8 V) from the 2026-10-22 reduced production, in the
user-confirmed map/height layout (styles/joule_map_height_larger_surfaces.json, three radius columns).

Fields are evaluated element-locally from the saved HEX20 nodal states with the frozen evaluator of
plot_je_8v_clouds_20261013 (exponential carriers, i.e. NT only). Before plotting, every state is checked
by integrating the evaluated J.E source with 5x5x5 Gauss points and comparing the total with the solver's
joule_W, and the evaluated nodal theta maximum with theta_max_K. The Y-half solution is mirrored to y<0
for display only (the half-domain branch is symmetric by construction; full/half equivalence was verified).
Hole fields are not plotted (release 2026-10-22). Color and height limits are shared across radii.
"""
from pathlib import Path
import argparse
import csv
import json
import sys
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize
import plot_je_8v_clouds_20261013 as fe
import run_je_production_reduced_20261022 as reduced

QA = Path('C:/Users/Administrator/.codex/skills/nature-figure/scripts')
sys.path.insert(0, str(QA))
from audit_panel_alignment import require_matplotlib_panel_alignment  # noqa: E402

base, r = reduced.base, reduced.r
STYLE_PATH = fe.SOLVER/'scripts/styles/joule_map_height_larger_surfaces.json'
RADII, VOLTAGES = [2, 3, 4], [6, 8]
FACES = {'midplane': 0., 'top': .5e-6, 'bottom': -.5e-6}
SUBDIVISIONS = 4
SPECS = [
    dict(id='Q_mid', key='Q', factor=1e-12, face='midplane', title='J·E 热源密度（中面）',
         label=r'$Q=\mathbf{E}\cdot\mathbf{J}$ ($10^{12}$ W m$^{-3}$)', cmap='heat'),
    dict(id='theta_mid', key='theta', factor=1., face='midplane', title='温升（中面）',
         label=r'$\theta$ (K)', cmap='heat'),
    dict(id='phi_mid', key='phi', factor=1., face='midplane', title='电势（中面）',
         label=r'$\phi$ (V)', cmap='heat'),
    dict(id='electron_bottom', key='dn_percent', factor=1., face='bottom', title='电子相对浓度变化（下表面）',
         label=r'$\Delta n/n_0$ (%)', cmap='signed'),
    dict(id='electron_top', key='dn_percent', factor=1., face='top', title='电子相对浓度变化（上表面）',
         label=r'$\Delta n/n_0$ (%)', cmap='signed'),
    dict(id='vm_mid', key='vm', factor=1e-6, face='midplane', title='耦合总等效应力（中面，含热效应）',
         label=r'$\sigma_{\mathrm{VM}}$ (MPa)', cmap='heat'),
]


def run_dir(a, v):
    return base.OUTPUT/f'H_r{a}_A48_R22_Z18_NT_V{v*10:03d}'


def verify(a, v, xyz, conn, state):
    d = run_dir(a, v)
    if not base.accepted(d): raise RuntimeError('Not an accepted state: '+d.name)
    m = fe.read_metrics(d/'metrics.txt')
    gp, gw = np.polynomial.legendre.leggauss(5)
    points = np.array([[x, y, z] for z in gp for y in gp for x in gp])
    weights = np.array([a_*b_*c_ for c_ in gw for b_ in gw for a_ in gw])
    joule = 0.
    for nodes in conn:
        e = fe.evaluate(xyz[nodes], state[nodes], points)
        joule += float((e['det']*weights) @ e['Q'])
    checks = dict(joule_W=(joule, m['joule_W']), nodal_theta_max_K=(float(state[:, 6].max()), m['theta_max_K']))
    out = {k: dict(evaluated=x, solver=y, relative=abs(x-y)/abs(y)) for k, (x, y) in checks.items()}
    if out['joule_W']['relative'] > 1e-6: raise RuntimeError(f'Evaluator/solver power mismatch {d.name}: {out}')
    if out['nodal_theta_max_K']['relative'] > 1e-3: raise RuntimeError(f'theta mismatch {d.name}: {out}')
    return out


def mirror(d):
    """Append the y -> -y image of a half-domain plane (scalars only) for display."""
    n = len(d['xyz'])
    xyz = np.vstack([d['xyz'], d['xyz']*np.array([1, -1, 1])])
    tri = np.vstack([d['triangles'], d['triangles'][:, ::-1]+n])
    out = dict(xyz=xyz, triangles=tri)
    for key in ('Q', 'theta', 'phi', 'dn_percent', 'vm'):
        out[key] = np.concatenate([d[key], d[key]])
    return out


def nice_limits(lo, hi, signed):
    if signed:
        m = max(abs(lo), abs(hi)); step = 10**np.floor(np.log10(m))/2
        m = np.ceil(m/step)*step
        return [-m, m]
    span = hi-lo; step = 10**np.floor(np.log10(span))/2
    return [np.floor(lo/step)*step, np.ceil(hi/step)*step]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    out = a.output; (out/'qa').mkdir(parents=True, exist_ok=True)
    style = json.loads(STYLE_PATH.read_text(encoding='utf-8-sig'))
    planes, audit, sources = {}, {}, {}
    for rad in RADII:
        mesh = base.half_mesh(rad, base.FINEST)
        xyz, conn = fe.read_mesh(mesh)
        for v in VOLTAGES:
            d = run_dir(rad, v)
            if r.read(d/'manifest.json')['mesh_sha256'] != r.sha(mesh): raise RuntimeError('Mesh identity '+d.name)
            nodes = np.loadtxt(d/'nodes.csv', delimiter=',', skiprows=1)
            if np.max(np.abs(nodes[:, 1:4]-xyz)) > 1e-19: raise RuntimeError('Node coordinates '+d.name)
            state = nodes[:, 4:11]
            audit[d.name] = verify(rad, v, xyz, conn, state)
            sources[d.name] = {f: r.sha(d/f) for f in ('nodes.csv', 'metrics.txt', 'completion.json')}
            for face, z in FACES.items():
                # 4x4 element-local subdivision: the A48/R22 face count is ~7x the template mesh, so the
                # sampled point density matches the template's 16x16 on its coarser mesh.
                planes[rad, v, face] = mirror(fe.plane(xyz, conn, state, z, SUBDIVISIONS))
            print('VERIFIED', d.name, flush=True)
    heat = LinearSegmentedColormap.from_list('reference_heat', style['palette_stops'], N=style['colormap_levels'])
    signed = LinearSegmentedColormap.from_list('signed_blue_green_red', [
        (0, '#0000bd'), (.18, '#006bff'), (.32, '#00cfcd'), (.43, '#00db48'), (.5, '#30e526'),
        (.60, '#7ae517'), (.70, '#e5ed00'), (.82, '#ff9c00'), (1, '#ed0000')], N=256)
    mpl.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': style['fonts'], 'font.size': 9,
                         'axes.titlesize': 11, 'axes.labelsize': 9, 'xtick.labelsize': 9, 'ytick.labelsize': 9,
                         'axes.linewidth': .8, 'svg.fonttype': 'none', 'pdf.fonttype': 42, 'axes.unicode_minus': False,
                         'xtick.direction': 'out', 'ytick.direction': 'out', 'figure.facecolor': 'white',
                         'savefig.facecolor': 'white'})
    extrema, limits_used = [], {}
    for v in VOLTAGES:
        for s in SPECS:
            vals = {rad: planes[rad, v, s['face']][s['key']]*s['factor'] for rad in RADII}
            lo = min(float(x.min()) for x in vals.values()); hi = max(float(x.max()) for x in vals.values())
            sgn = s['cmap'] == 'signed' or lo < 0 < hi and s['key'] == 'Q'
            limits = nice_limits(lo, hi, sgn) if s['key'] != 'phi' else [-float(v), 0.]
            limits_used[f"{s['id']}_{v}V"] = dict(data_min=lo, data_max=hi, limits=limits)
            cmap = signed if sgn else heat
            norm = Normalize(*limits)
            fig = plt.figure(figsize=style['comparison_figsize_inches'])
            layout = style['comparison_layout'].copy(); rows = layout.pop('rows'); columns = layout.pop('columns')
            gs = fig.add_gridspec(rows, columns, **layout)
            upper, lower = [], []
            for j, rad in enumerate(RADII):
                d = planes[rad, v, s['face']]; values = vals[rad]
                xy = d['xyz'][:, :2]*1e6
                tri = mtri.Triangulation(xy[:, 0], xy[:, 1], d['triangles'])
                ax = fig.add_subplot(gs[0, j]); upper.append(ax)
                ax.tripcolor(tri, values, shading='gouraud', cmap=cmap, norm=norm, rasterized=True)
                ang = np.linspace(0, 2*np.pi, 400)
                ax.plot(rad*np.cos(ang), rad*np.sin(ang), color='black', lw=.7)
                ax.axhline(0, color='white', lw=.4, ls=':')
                ax.set(aspect='equal', xlim=style['xy_limits_um'], ylim=style['xy_limits_um'],
                       xticks=style['xy_ticks_um'], yticks=style['xy_ticks_um'],
                       xlabel=r'$x$ ($\mu$m)', ylabel=r'$y$ ($\mu$m)', title=f'a = {rad} μm')
                bx = fig.add_subplot(gs[1, j], projection='3d'); lower.append(bx)
                bx.set_position(style['lower_row_enlargement']['axes_positions'][j]); bx.patch.set_alpha(0)
                bx.plot_trisurf(tri, values, cmap=cmap, norm=norm, linewidth=0, antialiased=False, shade=False,
                                rasterized=True)
                view, enl = style['surface_view'], style['lower_row_enlargement']
                bx.set(xlim=(-10, 10), ylim=(-10, 10), zlim=limits, xticks=[-10, 0, 10], yticks=[-10, 0, 10],
                       xlabel=r'$x$ ($\mu$m)', ylabel=r'$y$ ($\mu$m)')
                bx.set_box_aspect(view['box_aspect'], zoom=enl['projection_zoom'])
                bx.view_init(elev=view['elev'], azim=view['azim'])
                bx.tick_params(labelsize=enl['tick_fontsize'], pad=1)
                for axis in (bx.xaxis, bx.yaxis, bx.zaxis):
                    axis.pane.fill = False; axis._axinfo['grid']['color'] = (.6, .6, .6, .3)
                half = xy[:, 1] >= 0
                for kind, idx in (('min', np.argmin(np.where(half, values, np.inf))),
                                  ('max', np.argmax(np.where(half, values, -np.inf)))):
                    extrema.append(dict(figure=s['id'], voltage_V=v, radius_um=rad, face=s['face'], extreme=kind,
                                        value=float(values[idx]), x_um=float(xy[idx, 0]), y_um=float(xy[idx, 1])))
            cax = fig.add_axes(style['colorbar_axes'])
            cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax)
            cb.set_label(s['label'], fontsize=9); cb.ax.tick_params(labelsize=9)
            fig.suptitle(f"{s['title']} · NT · {v} V", fontsize=15, y=.98)
            top_ids = [f'map_a{x}' for x in RADII]; bottom_ids = [f'surface_a{x}' for x in RADII]
            require_matplotlib_panel_alignment(fig, axes=upper+lower, panel_ids=top_ids+bottom_ids,
                row_groups=[top_ids, bottom_ids], column_groups=[[x, y] for x, y in zip(top_ids, bottom_ids)],
                exemptions=[{'panels': bottom_ids, 'checks': ['column'],
                             'reason': 'Template enlarges the lower 3D row; rows keep equal widths and gutters.'}],
                json_out=out/'qa'/f"{s['id']}_{v}V_alignment.json", strict=True)
            stem = out/f"map_{s['id']}_{v}V"
            for ext in style['export']['formats']:
                fig.savefig(stem.with_suffix('.'+ext), dpi=style['export']['dpi'], bbox_inches='tight', facecolor='white')
            plt.close(fig)
            print('EXPORTED', stem.name, flush=True)
    with (out/'plane_extrema.csv').open('w', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(extrema[0])); w.writeheader(); w.writerows(extrema)
    fe.write_json(out/'provenance.json', dict(
        sources=sources, evaluator_sha256=fe.sha256(Path(fe.__file__)), style_sha256=fe.sha256(STYLE_PATH),
        script_sha256=fe.sha256(Path(__file__)), evaluator_checks=audit, limits=limits_used,
        display='Y-half solution mirrored to y<0 for display only; white dotted line marks the symmetry plane.',
        excluded='Hole fields not plotted (release 2026-10-22); LQ maps not drawn (evaluator implements exponential carriers only).'))


if __name__ == '__main__':
    main()
