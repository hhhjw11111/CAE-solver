"""Source-backed HEX20 plane maps for the saved J.E 8 V candidate.

Reconstruct the current-work heat source, rather than reusing drift-only plots.
No PDE solve, global field smoothing, or extrapolation is performed.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import shutil
from pathlib import Path

import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.ticker import MaxNLocator


def find_solver_root():
    for candidate in Path(__file__).resolve().parents:
        if (candidate / "src/physics/zno_candidate_m.f90").exists():
            return candidate
        if (candidate / "solver_v2/src/physics/zno_candidate_m.f90").exists():
            return candidate / "solver_v2"
    raise FileNotFoundError("Cannot locate the current solver_v2 source tree")


SOLVER = find_solver_root()
NATURAL = np.array([
    [-1,-1,-1], [1,-1,-1], [1,1,-1], [-1,1,-1],
    [-1,-1,1], [1,-1,1], [1,1,1], [-1,1,1],
    [0,-1,-1], [1,0,-1], [0,1,-1], [-1,0,-1],
    [0,-1,1], [1,0,1], [0,1,1], [-1,0,1],
    [-1,-1,0], [1,-1,0], [1,1,0], [-1,1,0],
], dtype=float)
CHARGE = 1.602176634e-19
KB = 1.380649e-23
T0 = 300.0
N0 = 6e22
P0 = (2.9e-4)**2 / N0
MU_N0 = 0.0205
MU_P0 = 0.005
BETA_N = 2.0
BETA_P = 2.0
C = np.zeros((6, 6))
C[0,0] = C[1,1] = 209.7e9
C[2,2] = 210.9e9
C[0,1] = C[1,0] = 121.1e9
C[0,2] = C[2,0] = C[1,2] = C[2,1] = 105.1e9
C[3,3] = C[4,4] = 42.47e9
C[5,5] = 44.3e9
PIEZO = np.zeros((3,6))
PIEZO[0,4] = PIEZO[1,3] = -0.37
PIEZO[2,:3] = [-0.62,-0.62,0.96]
THERMAL_BETA = C @ np.array([4.31e-6,4.31e-6,2.49e-6,0,0,0])

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Microsoft YaHei", "Arial", "DejaVu Sans"],
    "font.size": 8, "axes.titlesize": 10, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7,
    "svg.fonttype": "none", "pdf.fonttype": 42,
    "axes.unicode_minus": False, "axes.linewidth": 0.65,
    "figure.facecolor": "white", "savefig.facecolor": "white",
    "xtick.direction": "in", "ytick.direction": "in",
})


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for part in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(part)
    return h.hexdigest()


def write_json(path, obj):
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False),
                    encoding="utf-8")


def read_metrics(path):
    result = {}
    for line in path.read_text(encoding="utf-8-sig").splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        try:
            result[key] = float(value)
        except ValueError:
            result[key] = value.strip()
    return result


def read_mesh(path):
    with path.open(encoding="utf-8-sig") as stream:
        nn, ne = map(int, stream.readline().split())
        xyz = np.empty((nn, 3))
        for i in range(nn):
            row = stream.readline().split()
            if int(row[0]) != i + 1:
                raise ValueError("Unexpected node numbering")
            xyz[i] = list(map(float, row[1:]))
        conn = np.empty((ne, 20), dtype=int)
        for i in range(ne):
            row = list(map(int, stream.readline().split()))
            if len(row) != 22 or row[0] != i + 1 or row[1] != 1:
                raise ValueError("Unexpected HEX20 element or material")
            conn[i] = np.array(row[2:]) - 1
    return xyz, conn


def shapes(points):
    pts = np.asarray(points)
    fac = 1 + pts[:,None,:] * NATURAL[None,:,:]
    out = np.empty((len(pts),20))
    linear = np.sum(pts[:,None,:] * NATURAL[None,:,:], axis=2) - 2
    out[:,:8] = np.prod(fac[:,:8,:], axis=2) * linear[:,:8] / 8
    for node in range(8,20):
        zero = int(np.flatnonzero(NATURAL[node] == 0)[0])
        others = [a for a in range(3) if a != zero]
        out[:,node] = ((1-pts[:,zero]**2) * fac[:,node,others[0]]
                      * fac[:,node,others[1]] / 4)
    return out


def derivatives(points):
    pts = np.asarray(points)
    fac = 1 + pts[:,None,:] * NATURAL[None,:,:]
    out = np.zeros((len(pts),3,20))
    linear = np.sum(pts[:,None,:] * NATURAL[None,:,:], axis=2) - 2
    for axis in range(3):
        others = [a for a in range(3) if a != axis]
        out[:,axis,:8] = (NATURAL[None,:8,axis] / 8
                         * fac[:,:8,others[0]] * fac[:,:8,others[1]]
                         * (linear[:,:8]+fac[:,:8,axis]))
    for node in range(8,20):
        zero = int(np.flatnonzero(NATURAL[node] == 0)[0])
        others = [a for a in range(3) if a != zero]
        out[:,zero,node] = (-pts[:,zero]/2 * fac[:,node,others[0]]
                           * fac[:,node,others[1]])
        for axis in others:
            remaining = next(a for a in others if a != axis)
            out[:,axis,node] = ((1-pts[:,zero]**2)/4 * NATURAL[node,axis]
                               * fac[:,node,remaining])
    return out


def evaluate(coords, state_e, points):
    shape = shapes(points)
    dn = derivatives(points)
    jac = np.einsum("pkn,nj->pkj", dn, coords)
    det = np.linalg.det(jac)
    if np.any(det <= 0):
        raise ValueError("Non-positive sampled geometry")
    grad = np.linalg.solve(jac, dn)
    vals = shape @ state_e
    dvals = np.einsum("pkn,ns->pks", grad, state_e)
    electric = -dvals[:,:,3]
    theta = vals[:,6]
    temp = T0 + theta
    if np.any(temp <= 0):
        raise ValueError("Non-positive absolute temperature")
    nr, pr = np.exp(vals[:,4]), np.exp(vals[:,5])
    mu_n = MU_N0 * (temp/T0)**(-BETA_N)
    mu_p = MU_P0 * (temp/T0)**(-BETA_P)
    diff_n, diff_p = KB * temp / CHARGE * mu_n, KB * temp / CHARGE * mu_p
    drift_n = CHARGE * N0 * nr[:,None] * mu_n[:,None] * electric
    drift_p = CHARGE * P0 * pr[:,None] * mu_p[:,None] * electric
    current_n = drift_n + CHARGE * N0 * nr[:,None] * diff_n[:,None] * dvals[:,:,4]
    current_p = drift_p - CHARGE * P0 * pr[:,None] * diff_p[:,None] * dvals[:,:,5]
    q_drift = np.sum(electric*(drift_n+drift_p), axis=1)
    q_je = np.sum(electric*(current_n+current_p), axis=1)
    strain = np.column_stack((
        dvals[:,0,0],dvals[:,1,1],dvals[:,2,2],
        dvals[:,2,1]+dvals[:,1,2],
        dvals[:,0,2]+dvals[:,2,0],
        dvals[:,1,0]+dvals[:,0,1],
    ))
    stress_thermal = -theta[:,None] * THERMAL_BETA
    stress = strain @ C.T - electric @ PIEZO + stress_thermal
    vm = np.sqrt(0.5*((stress[:,0]-stress[:,1])**2
                      +(stress[:,1]-stress[:,2])**2
                      +(stress[:,2]-stress[:,0])**2)
                 + 3*np.sum(stress[:,3:]**2, axis=1))
    return {
        "xyz": shape @ coords, "det": det, "phi": vals[:,3],
        "Evec": electric, "Q": q_je, "Q_drift_same_field": q_drift,
        "Q_diffusion_work": q_je-q_drift, "theta": theta,
        "nr": nr, "pr": pr, "dn_percent": 100*(nr-1),
        "dp_percent": 100*(pr-1), "eta_n": vals[:,4], "eta_p": vals[:,5],
        "n_m3": N0*nr, "p_m3": P0*pr, "u": vals[:,:3],
        "mu_ratio": mu_n/MU_N0, "stress": stress,
        "vm": vm, "stress_thermal": stress_thermal,
    }


def match_native_samples(xyz, conn, state, run, metrics):
    native = np.loadtxt(run/"field_samples.csv", delimiter=",", skiprows=1)
    if native.shape != (len(conn)*125,21):
        raise ValueError("Unexpected native sample count or columns")
    gp, gw = np.polynomial.legendre.leggauss(5)
    points = np.array([[r,s,t] for t in gp for s in gp for r in gp])
    weights = np.array([a*b*c for c in gw for b in gw for a in gw])
    peaks = np.max(np.abs(native),axis=0)
    error_abs = np.zeros(21)
    error_point = np.zeros(21)
    floor = np.maximum(peaks*1e-5,np.finfo(float).tiny)
    floor[:3] = 1e-11
    floor[4:7] = 10
    floor[9:11] = 1e-5
    floor[15:21] = 100
    joule, diffusion, volume = 0.0,0.0,0.0
    maximum_theta = maximum_vm = 0.0
    q_min, q_max = np.inf,-np.inf
    for element,nodes in enumerate(conn):
        d = evaluate(xyz[nodes],state[nodes],points)
        w = d["det"] * weights
        actual = np.column_stack((
            d["xyz"],w,d["Evec"],d["Q"],d["theta"],d["nr"],d["pr"],
            d["u"],d["mu_ratio"],d["stress"],
        ))
        reference = native[element*125:(element+1)*125]
        error_abs = np.maximum(error_abs,np.max(np.abs(actual-reference),axis=0))
        error_point = np.maximum(error_point,np.max(
            np.abs(actual-reference)/np.maximum(np.abs(reference),floor),axis=0))
        joule += float(w @ d["Q"])
        diffusion += float(w @ d["Q_diffusion_work"])
        volume += float(w.sum())
        maximum_theta = max(maximum_theta,float(d["theta"].max()))
        maximum_vm = max(maximum_vm,float(d["vm"].max()))
        q_min,q_max = min(q_min,float(d["Q"].min())),max(q_max,float(d["Q"].max()))
        if element % 600 == 0:
            print(f"NATIVE_SAMPLE_AUDIT={element}/{len(conn)}", flush=True)
    peak_relative = error_abs / np.maximum(peaks,np.finfo(float).tiny)
    checks = {
        "joule_W": (joule,metrics["joule_W"]),
        "diffusion_work_W": (diffusion,metrics["diffusion_work_W"]),
        "total_volume_m3": (volume,metrics["total_volume_m3"]),
        "theta_max_K": (maximum_theta,metrics["theta_max_K"]),
        "stress_von_mises_max_Pa": (maximum_vm,metrics["stress_von_mises_max_Pa"]),
    }
    checks = {key:{"reconstructed":a,"native":b,
                   "relative_error":abs(a-b)/max(abs(a),abs(b),np.finfo(float).tiny)}
              for key,(a,b) in checks.items()}
    if peak_relative.max()>5e-7 or max(v["relative_error"] for v in checks.values())>5e-7:
        raise ValueError(f"Native source mismatch: {peak_relative.max()}; {checks}")
    return {
        "status": "PASS_POSTPROCESSING_RECONSTRUCTION",
        "sample_count": len(native),
        "native_columns": (run/"field_samples.csv").open().readline().strip().split(","),
        "maximum_error_normalized_by_column_peak": peak_relative.tolist(),
        "maximum_pointwise_error_with_documented_floors": error_point.tolist(),
        "pointwise_normalization_floor": floor.tolist(),
        "metric_checks": checks,
        "native_sampling_Q_range_W_m3": [float(q_min),float(q_max)],
        "interpretation": "Checks plotting against the saved solver fields; not physical, energy-closure or mesh certification.",
    }


def plane(xyz,conn,state,z,subdivisions):
    pieces, triangles, ids = {},[],[]
    offset = 0
    t = np.linspace(-1,1,subdivisions+1)
    side = 1.0 if z > 0 else -1.0
    points = np.array([[r,s,side] for s in t for r in t])
    count = 0
    for element,nodes in enumerate(conn):
        coords = xyz[nodes]
        selected = (abs(coords[:,2].max()-z)<1e-19 if z>0
                    else abs(coords[:,2].min()-z)<1e-19)
        if not selected:
            continue
        d = evaluate(coords,state[nodes],points)
        if np.max(np.abs(d["xyz"][:,2]-z))>1e-18:
            raise ValueError("Plane is not represented exactly")
        for key,value in d.items():
            if key != "det":
                pieces.setdefault(key,[]).append(value)
        ids.extend([element+1]*len(points))
        width = subdivisions+1
        for j in range(subdivisions):
            for i in range(subdivisions):
                a = offset+j*width+i
                triangles.extend([[a,a+1,a+width+1],[a,a+width+1,a+width]])
        offset += len(points)
        count += 1
    if count == 0:
        raise ValueError("No element face intersects the plane")
    d = {key:np.concatenate(value) for key,value in pieces.items()}
    d.update({"triangles":np.array(triangles),"element_id":np.array(ids)})
    if any(not np.isfinite(value).all() for key,value in d.items() if key not in ("triangles","element_id")):
        raise ValueError("Non-finite plotted field")
    d["patch_count"] = np.array(count)
    return d


SPECS = [
    ("joule_density","Q",1e-12,"焦耳热密度 · J·E",r"$Q_{J\cdot E}$ ($10^{12}$ W m$^{-3}$)","heat"),
    ("potential","phi",1,"电势",r"$\phi$ (V)","viridis"),
    ("electron_density","dn_percent",1,"电子浓度变化",r"$(n/n_0-1)$ (%)","signed"),
    ("hole_density","dp_percent",1,"空穴浓度变化",r"$(p/p_0-1)$ (%)","signed"),
    ("coupled_stress","vm",1e-6,"耦合总应力 · von Mises",r"$\sigma_{\mathrm{VM}}$ (MPa)","magma"),
    ("thermal_constitutive_stress","thermal_xx",1e-6,"热膨胀本构项 · xx",r"$-\beta_{xx}\theta$ (MPa)","Blues_r"),
]
TEMPERATURE = ("temperature","theta",1,"温升",r"$\theta$ (K)","magma")


def values(d,key):
    return d["stress_thermal"][:,0] if key=="thermal_xx" else d[key]


def scales(planes):
    out = {}
    for name,key,factor,title,label,palette in SPECS+[TEMPERATURE]:
        lo = min(float((values(d,key)*factor).min()) for d in planes.values())
        hi = max(float((values(d,key)*factor).max()) for d in planes.values())
        if palette=="signed" or (palette=="heat" and lo<0):
            limit = max(abs(lo),abs(hi))
            lo,hi = -limit,limit
        if hi-lo < 1e-14:
            eps = max(abs(lo),1)*1e-8
            lo,hi = lo-eps,hi+eps
        out[name] = [lo,hi]
    return out


def palette(spec,bounds,heat):
    name = spec[-1]
    if name=="heat":
        return "RdBu_r" if bounds[0]<0 else heat
    return "RdBu_r" if name=="signed" else name


def cloud(fig,ax,d,spec,bounds,heat):
    name,key,factor,title,label,_ = spec
    xyz = d["xyz"]*1e6
    tri = mtri.Triangulation(xyz[:,0],xyz[:,1],d["triangles"])
    artist = ax.tripcolor(tri,values(d,key)*factor,shading="gouraud",
                          cmap=palette(spec,bounds,heat),norm=Normalize(*bounds),
                          rasterized=True)
    curve = np.linspace(0,2*np.pi,601)
    ax.plot(3*np.cos(curve),3*np.sin(curve),color="black",lw=.65)
    ax.set(aspect="equal",xlim=(-10,10),ylim=(-10,10),
           xticks=[-10,-5,0,5,10],yticks=[-10,-5,0,5,10],
           xlabel=r"$x$ ($\mu$m)",ylabel=r"$y$ ($\mu$m)",title=title)
    for edge in ax.spines.values():
        edge.set_linewidth(.65)
    cb = fig.colorbar(artist,ax=ax,fraction=.046,pad=.025)
    cb.set_label(label,fontsize=8)
    cb.ax.tick_params(labelsize=7,direction="in")
    cb.locator = MaxNLocator(5)
    cb.update_ticks()
    if max(abs(bounds[0]),abs(bounds[1]))<.01:
        cb.formatter.set_powerlimits((0,0))
        cb.update_ticks()
    return artist


def export(fig,stem,qa):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    # The export keeps the full canvas so final dimensions remain explicit.
    overflows = []
    for ax in fig.axes:
        box = ax.get_tightbbox(renderer)
        if box is not None and (box.x0 < -1 or box.y0 < -1
                                or box.x1 > fig.bbox.width+1
                                or box.y1 > fig.bbox.height+1):
            overflows.append([float(v) for v in box.extents])
    if overflows:
        raise ValueError(f"Figure labels outside canvas: {stem.name}: {overflows}")
    fig.savefig(stem.with_suffix(".png"),dpi=330,facecolor="white")
    fig.savefig(stem.with_suffix(".svg"),dpi=330,facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"),dpi=330,facecolor="white")
    if stem.name=="JE_V080_clouds_midplane":
        fig.savefig(stem.with_suffix(".tiff"),dpi=600,facecolor="white",
                    pil_kwargs={"compression":"tiff_lzw"})
    svg = stem.with_suffix(".svg").read_text(encoding="utf-8")
    if "<text" not in svg:
        raise ValueError("Vector text has not remained editable")
    qa.append({
        "figure":stem.name,"canvas_inches":list(fig.get_size_inches()),
        "label_overflows":overflows,"editable_svg_text":True,
        "png_dpi":330,"tiff_dpi":600 if stem.name=="JE_V080_clouds_midplane" else None,
    })
    plt.close(fig)
    print("EXPORTED="+stem.name,flush=True)


def render(out,planes,bounds,heat):
    exports = []
    for face,z in [("midplane",0),("top",.5),("bottom",-.5)]:
        local_bounds = scales({face:planes[face]})
        fig,axes = plt.subplots(2,3,figsize=(11.2,7.3),layout="constrained")
        fig.get_layout_engine().set(w_pad=.06,h_pad=.10,wspace=.06,hspace=.04)
        for ax,spec,letter in zip(axes.flat,SPECS,"abcdef"):
            cloud(fig,ax,planes[face],spec,local_bounds[spec[0]],heat)
            ax.text(0,1.085,letter,transform=ax.transAxes,weight="bold",fontsize=11)
        fig.suptitle(f"J·E 新分支 · 8 V · r = 3 μm · z = {z:g} μm",fontsize=13)
        fig.supxlabel("各面独立色标；跨截面对比请看共用色标图。热本构项不是单独平衡后的纯热应力。",fontsize=8)
        export(fig,out/("JE_V080_clouds_"+face),exports)
    for spec in SPECS+[TEMPERATURE]:
        fig,axes = plt.subplots(1,3,figsize=(11.8,4.3),layout="constrained")
        for ax,face,z in zip(axes,["bottom","midplane","top"],[-.5,0,.5]):
            cloud(fig,ax,planes[face],spec,bounds[spec[0]],heat)
            ax.set_title(f"z = {z:g} μm",fontsize=10)
        fig.suptitle("J·E · 8 V · r = 3 μm | "+spec[3],fontsize=13)
        export(fig,out/(spec[0]+"_three_planes"),exports)
    fig,axes = plt.subplots(2,3,figsize=(11.3,7.3),layout="constrained")
    for row,spec in enumerate(SPECS[2:4]):
        for col,(face,z) in enumerate([("bottom",-.5),("midplane",0),("top",.5)]):
            ax = axes[row,col]
            cloud(fig,ax,planes[face],spec,bounds[spec[0]],heat)
            ax.set_title(spec[3]+f" · z = {z:g} μm",fontsize=10)
    fig.suptitle("J·E · 8 V · 电子与空穴的浓度变化",fontsize=13)
    fig.supxlabel("相对背景浓度的变化；各载流子的三幅图共用色标。",fontsize=8)
    export(fig,out/"carriers_three_planes",exports)
    return exports


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    root = SOLVER/"runs/JE_6_8_HALF_V_20261013"
    parser.add_argument("--run",type=Path,default=root/"V080")
    parser.add_argument("--mesh",type=Path,default=SOLVER/"cases/HV_SURFACE_20261008/XY0_Z2/H1_XY0_Z2.mesh")
    parser.add_argument("--output",type=Path,default=root/"clouds_V080_20261013")
    parser.add_argument("--subdivisions",type=int,default=16)
    parser.add_argument("--render-only",action="store_true")
    args = parser.parse_args()
    run,mesh,out = args.run.resolve(),args.mesh.resolve(),args.output.resolve()
    if args.render_only:
        planes = {face:dict(np.load(out/"source_data"/(face+".npz")))
                  for face in ["midplane","top","bottom"]}
        metadata = json.loads((out/"metadata.json").read_text(encoding="utf-8"))
        bounds = scales(planes)
        metadata["shared_color_limits"] = bounds
        metadata["local_color_limits"] = {face:scales({face:d}) for face,d in planes.items()}
        write_json(out/"metadata.json",metadata)
        shutil.copy2(root/"cloud_figure_contract_20261013.json",out/"figure_contract.json")
    else:
        out.mkdir(parents=True,exist_ok=False)
        (out/"source_data").mkdir()
        contract = root/"cloud_figure_contract_20261013.json"
        shutil.copy2(contract,out/"figure_contract.json")
        metrics = read_metrics(run/"metrics.txt")
        if (metrics.get("model")!="NT"
            or metrics.get("heat_source")!="E_DOT_TOTAL_CURRENT_DIAGNOSTIC_NOT_THERMODYNAMIC_CLOSURE"
            or metrics.get("case_acceptance")!="PASS_NOT_MESH_ACCEPTANCE"
            or metrics["voltage_V"]!=8
            or metrics["final_voltage_factor"]!=1
            or metrics["final_joule_factor"]!=1):
            raise ValueError("Wrong case, model, heat source or load")
        xyz,conn = read_mesh(mesh)
        nodes = np.loadtxt(run/"nodes.csv",delimiter=",",skiprows=1)
        if len(nodes)!=len(xyz) or np.max(abs(nodes[:,1:4]-xyz))>1e-19:
            raise ValueError("Mesh and saved nodal field mismatch")
        if not np.isfinite(nodes).all():
            raise ValueError("Non-finite saved field")
        state = nodes[:,4:11]
        print("SOURCE_NODES="+str(len(nodes)),flush=True)
        audit = match_native_samples(xyz,conn,state,run,metrics)
        write_json(out/"native_sample_audit.json",audit)
        print("NATIVE_RECONSTRUCTION=PASS",flush=True)
        planes = {face:plane(xyz,conn,state,z,args.subdivisions)
                  for face,z in [("midplane",0),("top",.5e-6),("bottom",-.5e-6)]}
        bounds = scales(planes)
        rows = []
        for face,d in planes.items():
            np.savez_compressed(out/"source_data"/(face+".npz"),**d)
            table = np.column_stack((
                d["xyz"]*1e6,d["element_id"],d["Q"],d["phi"],d["n_m3"],d["p_m3"],
                d["nr"],d["pr"],d["dn_percent"],d["dp_percent"],d["vm"]*1e-6,
                d["stress"]*1e-6,d["stress_thermal"][:,0]*1e-6,d["theta"],
            ))
            np.savetxt(out/"source_data"/(face+".csv"),table,delimiter=",",
                       fmt="%.12e",comments="",
                       header="x_um,y_um,z_um,element_id,Q_JE_W_m3,phi_V,n_m3,p_m3,n_over_n0,p_over_p0,dn_percent,dp_percent,von_mises_MPa,sxx_MPa,syy_MPa,szz_MPa,syz_MPa,szx_MPa,sxy_MPa,thermal_constitutive_xx_MPa,theta_K")
            for spec in SPECS+[TEMPERATURE]:
                v = values(d,spec[1])*spec[2]
                for extreme,idx in [("min",int(np.argmin(v))),("max",int(np.argmax(v)))]:
                    rows.append({"plane":face,"field":spec[0],"extreme":extreme,
                                 "value":float(v[idx]),"x_um":float(d["xyz"][idx,0]*1e6),
                                 "y_um":float(d["xyz"][idx,1]*1e6),"z_um":float(d["xyz"][idx,2]*1e6)})
        with (out/"plane_extrema.csv").open("w",newline="",encoding="utf-8-sig") as stream:
            writer = csv.DictWriter(stream,fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        relevant = [run/"nodes.csv",run/"metrics.txt",run/"field_samples.csv",
                    run.parent/"input_V080.nml",mesh,
                    SOLVER/"src/physics/material_point_m.f90",
                    SOLVER/"src/physics/zno_candidate_m.f90",
                    SOLVER/"src/physics/transport_m.f90",
                    SOLVER/"src/physics/thermo_piezo_m.f90",
                    SOLVER/"src/mesh/hex20_shape_m.f90",
                    SOLVER/"build/ARC_JE_FINAL_20261012/build_manifest.json"]
        manifest = json.loads(relevant[-1].read_text(encoding="utf-8-sig"))
        manifest_hashes = {item["path"]:item["sha256"] for item in manifest["files"]}
        source_build_matches = {}
        for path in relevant:
            if path.is_relative_to(SOLVER):
                rel = path.relative_to(SOLVER).as_posix()
                if rel in manifest_hashes:
                    source_build_matches[rel] = sha256(path)==manifest_hashes[rel]
        if not all(source_build_matches.values()):
            raise ValueError("Current constitutive source differs from the run build")
        with (out/"source_sha256.csv").open("w",newline="",encoding="utf-8-sig") as stream:
            writer = csv.writer(stream)
            writer.writerow(["source","sha256"])
            for path in relevant:
                writer.writerow([path.relative_to(SOLVER).as_posix(),sha256(path)])
        metadata = {
            "case":"JE_GRID_V080","model":"NT","voltage_magnitude_V":8,
            "target_voltage_V":-8,"hole_radius_um":3,
            "mesh":"XY0_Z2","node_count":len(xyz),"element_count":len(conn),
            "n0_m3":N0,"p0_m3":P0,
            "thermal_beta_Pa_per_K":THERMAL_BETA.tolist(),
            "shared_color_limits":bounds,
            "local_color_limits":{face:scales({face:d}) for face,d in planes.items()},
            "subdivisions_per_original_element_edge":args.subdivisions,
            "patch_count_per_plane":{face:int(d["patch_count"]) for face,d in planes.items()},
            "points_per_plane":{face:len(d["xyz"]) for face,d in planes.items()},
            "current_source_matches_native_build":source_build_matches,
            "heat_source":"E dot (Jn+Jp), includes diffusion current work",
            "stress_definition":"C:strain - e^T E - beta*theta",
            "thermal_constitutive_xx_definition":"-beta_xx*theta; not the total stress caused by the thermal load after re-equilibration",
            "planes":"Exact HEX20 z=0, +/-0.5 um evaluations; one-sided material gradients, not thickness averages.",
            "integrity":"Independent element-patch triangles, no averaging between element gradients, no cross-hole interpolation.",
            "status":"Converged numerical candidate, not physically or spatially certified.",
            "electrical_work_relative_gap":abs(metrics["total_current_work_W"]-metrics["terminal_power_W"])/abs(metrics["total_current_work_W"]),
            "native_global_theta_max_K":metrics["theta_max_K"],
            "native_global_von_mises_max_MPa":metrics["stress_von_mises_max_Pa"]*1e-6,
        }
        write_json(out/"metadata.json",metadata)
    style = json.loads((SOLVER/"scripts/styles/joule_map_height_frozen.json").read_text(encoding="utf-8"))
    heat = LinearSegmentedColormap.from_list("reference_heat",style["palette_stops"],N=256)
    qa = render(out,planes,bounds,heat)
    write_json(out/"export_qa.json",{
        "backend":"Python matplotlib","exports":qa,
        "visual_review":"PENDING","physics_acceptance":"NOT_ESTABLISHED",
    })
    shutil.copy2(Path(__file__).resolve(),out/Path(__file__).name)
    readme = (
        "J·E 新分支 8 V 云图\n"
        "孔半径 3 μm；板尺寸 20 × 20 × 1 μm；NT；XY0_Z2。\n\n"
        "JE_V080_clouds_midplane/top/bottom：中面、顶面、底面六幅组合图。\n"
        "六幅组合图用各面的独立线性色标，便于看清弱中面载流子扰动。\n"
        "*_three_planes：每个变量的底面—中面—顶面比较；各变量跨截面共用线性色标。\n"
        "电子、空穴使用相对浓度变化百分比。绝对浓度也保存在 source_data/*.csv。\n"
        "背景电子浓度 n0 = 6e22 m^-3；背景空穴浓度 p0 = 1.4016666667e-30 m^-3。\n"
        "空穴的相对变化不能被解读为很大的绝对空穴密度。\n\n"
        "Q 是完整漂移—扩散电流的 J·E，并非原漂移电导热源。\n"
        "sigma_VM 是全部耦合应力的 von Mises 值，含弹性、压电、热贡献。\n"
        "热膨胀本构项 xx = -beta_xx*theta 是本构式的直接温度项，"
        "并非单独平衡后的纯热应力，不能与总应力峰值等同。\n\n"
        "绘图从原 HEX20 单元的保存节点场直接求值，"
        "积分点 Q、应力及其他场与 native field_samples.csv 对拍；"
        "392000 个点的审核见 native_sample_audit.json。\n"
        "每个三角形属于一个原单元，不跨孔洞插值，不对单元间梯度平滑。\n"
        "中面梯度使用正 z 侧单元；表面梯度从材料内部取单侧迹。\n"
        "色标包含所有绘图采样值；不使用分位截断或非线性色阶。\n\n"
        "数值解已收敛，仍有约 2.97% 的端口功率与体积分电功率差异，"
        "尚未完成新分支的网格/稳定化和物理闭合验证。图件不是验收证书。\n"
        "截面和表面极值可能大于体积 Gauss 点采样极值，二者采样位置不同。\n\n"
        "重现：在 solver_v2 中运行 python scripts/plot_je_8v_clouds_20261013.py "
        "--output runs/JE_6_8_HALF_V_20261013/<新的输出目录>。\n"
    )
    (out/"README.txt").write_text(readme,encoding="utf-8")
    print("OUTPUT="+str(out),flush=True)
    print(json.dumps({"color_limits":bounds,"source_status":"PASS"},ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
