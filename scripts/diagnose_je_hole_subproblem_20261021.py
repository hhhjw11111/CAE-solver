"""Conditional hole solves for numerical diagnosis, never seven-field results.

Keep the six saved non-hole fields fixed and independently assemble the existing
R=0 conventional hole current and its residual-based SUPG.  Only a constant
q*p0 is cancelled; no physical coefficient, heat source or native gate changes.
Order/metric comparisons must subsequently be solved by the native coupled code.
"""
from pathlib import Path
import argparse
import gc
import json
import sys
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'.codex_tmp/je_controlled_deps'))
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
import diagnose_je_power_gap_20261015 as fe
from run_hv_sweep import sha, read_metrics, classify


def save(path, value):
    with Path(path).open('x', encoding='utf-8') as stream:
        json.dump(value, stream, indent=2, ensure_ascii=False, allow_nan=False)


class HoleOperator:
    def __init__(self, xyz, conn, state, order=3, directional=True, scale=1.):
        self.xyz, self.conn = xyz, conn
        self.order, self.directional, self.scale = order, directional, scale
        self.fixed = abs(abs(xyz[:, 0])-1e-5) < 1e-15
        self.free = ~self.fixed
        self.chunks = []
        points, weights, n, dn, hn = fe.rule(order)
        self.n = n
        for start in range(0, len(conn), 64):
            ids = conn[start:start+64]
            x, xe = xyz[ids], state[ids]
            jac = np.einsum('qka,eaj->eqkj', dn, x)
            det = np.linalg.det(jac)
            if np.any(det <= 0):
                raise ValueError('Nonpositive element geometry')
            ji = np.linalg.inv(jac)
            dg = np.einsum('eqik,qka->eqia', ji, dn)
            hx = np.einsum('qabn,enk->eqabk', hn, x)
            corrected = hn[None, ...]-np.einsum('eqabk,eqkn->eqabn', hx, dg)
            natural_metric = np.einsum('eqia,eqib->eqab', ji, ji)
            lap = np.einsum('eqab,eqabn->eqn', natural_metric, corrected)
            val = np.einsum('qa,eaf->eqf', n, xe)
            grad = np.einsum('eqia,eaf->eqif', dg, xe)
            E, gt, T = -grad[:, :, :, 3], grad[:, :, :, 6], 300.+val[:, :, 6]
            if np.any(T <= 0):
                raise ValueError('Nonpositive saved temperature')
            mu = .005*(T/300.)**-2
            diff = fe.fe.KB*T/fe.fe.CHARGE*mu
            mut, dt = -2*mu/T, -diff/T
            v = mu[:, :, None]*E
            divv = -mu*np.einsum('eqa,ea->eq', lap, xe[:, :, 3])+mut*np.sum(E*gt, axis=2)
            adv = np.einsum('eqia,eqi->eqa', dg, v)
            B = adv-diff[:, :, None]*lap-dt[:, :, None]*np.einsum('eqia,eqi->eqa', dg, gt)
            G = np.einsum('eqia,eqja->eqij', ji, ji)
            trace = np.trace(G, axis1=2, axis2=3)
            speed2 = np.sum(v*v, axis=2)
            advmetric = np.einsum('eqi,eqij,eqj->eq', v, G, v)
            epsilon = 1e-12*diff*diff*trace
            gdir = ((advmetric+epsilon*trace/3)/(speed2+epsilon)) if directional else trace
            tau = scale/np.sqrt(advmetric+(4*diff*gdir)**2)
            self.chunks.append(dict(ids=ids, dn=dg, w=det*weights, diff=diff,
                                    adv=adv, B=B, divv=divv, tau=tau))
        self.rows = np.repeat(conn, 20, axis=1).ravel()
        self.cols = np.tile(conn, (1, 20)).ravel()

    def assemble(self, eta, tangent=True):
        residual = np.zeros(len(self.xyz))
        matrices = []
        for d in self.chunks:
            ids, dn = d['ids'], d['dn']
            en = eta[ids]
            ev = en@self.n.T
            if np.any(abs(ev) > 100):
                raise ValueError('Diagnostic log-density trial outside bounded range')
            rho = np.exp(ev)
            g = np.einsum('eqia,ea->eqi', dn, en)
            gd = np.einsum('eqia,eqi->eqa', dn, g)
            strong = d['divv']+np.einsum('eqa,ea->eq', d['B'], en)-d['diff']*np.sum(g*g, axis=2)
            flux_test = d['adv']-d['diff'][:, :, None]*gd
            test = flux_test-d['tau'][:, :, None]*d['adv']*strong[:, :, None]
            wrho = d['w']*rho
            re = np.einsum('eq,eqa->ea', wrho, test)
            np.add.at(residual, ids.ravel(), re.ravel())
            if tangent:
                ds = d['B']-2*d['diff'][:, :, None]*gd
                ke = np.einsum('eq,eqa,qb->eab', wrho, test, self.n)
                ke -= np.einsum('eq,eqia,eqib->eab', wrho*d['diff'], dn, dn)
                ke -= np.einsum('eq,eqa,eqb->eab', wrho*d['tau'], d['adv'], ds)
                matrices.append(ke.reshape(-1))
        if not tangent:
            return residual
        jac = coo_matrix((np.concatenate(matrices), (self.rows, self.cols)),
                         shape=(len(eta), len(eta))).tocsr()
        return residual, jac

    def solve(self, initial):
        eta = initial.copy()
        eta[self.fixed] = 0.
        history = []
        for iteration in range(20):
            r, j = self.assemble(eta)
            jf = j[self.free][:, self.free]
            row_scale = np.asarray(abs(jf).max(axis=1).toarray()).ravel()
            if np.any(row_scale <= 0):
                raise ValueError('Unscaled zero row')
            norm = float(np.max(abs(r[self.free])/row_scale))
            if norm < 1e-11:
                history.append(dict(iteration=iteration, scaled_residual=norm, update=0.))
                return eta, history
            delta = spsolve(jf.multiply(1/row_scale[:, None]).tocsc(), -r[self.free]/row_scale)
            if not np.isfinite(delta).all():
                raise ValueError('Nonfinite conditional Newton direction')
            step = min(1., .5/max(float(abs(delta).max()), 1e-300))
            for backtrack in range(24):
                trial = eta.copy()
                trial[self.free] += step*delta
                rt = self.assemble(trial, False)
                newnorm = float(np.max(abs(rt[self.free])/row_scale))
                if newnorm < norm:
                    break
                step *= .5
            else:
                raise RuntimeError('Conditional Newton line search failed')
            history.append(dict(iteration=iteration, scaled_residual=norm,
                                update=float(step*abs(delta).max()), step=step))
            print(json.dumps(history[-1]), flush=True)
            eta = trial
        raise RuntimeError('Conditional Newton did not converge')


def ring_mean(xyz, conn, eta):
    gp, gw = np.polynomial.legendre.leggauss(5)
    gz, wz = np.polynomial.legendre.leggauss(3)
    points = np.array([(x,y,z) for z in gz for y in gp for x in gp])
    weights = np.array([a*b*c for c in wz for b in gw for a in gw])
    n, dn = fe.fe.shapes(points), fe.fe.derivatives(points)
    volume = value = 0.
    for start in range(0, len(conn), 64):
        ids = conn[start:start+64]
        pos = np.einsum('qa,eai->eqi', n, xyz[ids])
        det = np.linalg.det(np.einsum('qka,eaj->eqkj', dn, xyz[ids]))
        radius = np.linalg.norm(pos[:, :, :2], axis=2)
        w = det*weights*((radius >= 3e-6)&(radius <= 3.5e-6))
        volume += float(w.sum())
        value += float(np.sum(w*abs(np.expm1(eta[ids]@n.T))))
    return value/volume


def check():
    xyz = fe.fe.NATURAL*np.array([1e-6, 2e-6, .5e-6])
    conn = np.arange(20)[None, :]
    state = np.zeros((20, 7))
    state[:, 3] = -xyz@np.array([1e4, -2e4, 3e4])
    eta = xyz@np.array([2e4, -1e4, 3e4])+.1
    tests = []
    for directional in (False, True):
        op = HoleOperator(xyz, conn, state, 3, directional)
        r, j = op.assemble(eta)
        direction = np.sin(np.arange(20)+1.)
        eps = 1e-6
        fd = (op.assemble(eta+eps*direction, False)-op.assemble(eta-eps*direction, False))/(2*eps)
        analytic = j@direction
        error = float(np.linalg.norm(fd-analytic)/np.linalg.norm(analytic))
        if error >= 1e-7:
            raise AssertionError(('Conditional tangent', error))
        exact = -state[:, 3]/(fe.fe.KB*300/fe.fe.CHARGE)
        rr = op.assemble(exact, False)
        if np.max(abs(rr)) > 1e-15:
            raise AssertionError(('Boltzmann zero-current residual', float(abs(rr).max())))
        tests.append(dict(directional=directional, tangent_relative_error=error,
                          boltzmann_max_residual=float(abs(rr).max())))
    return dict(status='LOCAL_TANGENT_AND_ZERO_CURRENT_CHECKS_PASSED', checks=tests)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--orders', nargs='+', type=int, choices=(3, 5, 7), default=[3, 5])
    p.add_argument('--metric', choices=('directional', 'trace'), default='directional')
    p.add_argument('--scale', type=float, default=1.)
    p.add_argument('--perturbation', type=float, default=0.)
    a = p.parse_args()
    a.run, a.output = a.run.resolve(), a.output.resolve()
    if not np.isfinite([a.scale, a.perturbation]).all() or a.scale <= 0:
        raise ValueError('Invalid diagnostic settings')
    manifest = json.loads((a.run/'manifest.json').read_text(encoding='utf-8'))
    execution = json.loads((a.run/'execution.json').read_text(encoding='utf-8'))
    status, metrics = classify(a.run, execution['exit_code'])
    if status != 'SOLVE_FLUX_PASS_MESH_PENDING':
        raise ValueError('Complete native solve required')
    if sha(Path(manifest['mesh'])) != manifest['mesh_sha256']:
        raise ValueError('Changed native geometry')
    xyz, conn = fe.fe.read_mesh(Path(manifest['mesh']))
    data = np.loadtxt(a.run/'nodes.csv', delimiter=',', skiprows=1)
    if np.max(abs(data[:, 1:4]-xyz)) >= 1e-19:
        raise ValueError('Native field geometry mismatch')
    state, reference = data[:, 4:11], data[:, 9].copy()
    a.output.mkdir(parents=True, exist_ok=False)
    save(a.output/'local_checks.json', check())
    files = [a.run/x for x in ('manifest.json', 'execution.json', 'nodes.csv', 'input.nml', 'metrics.txt')]
    files += [Path(manifest['mesh']), Path(__file__), Path(fe.__file__)]
    hashes = {str(f): sha(f) for f in files}
    reports = []
    for order in a.orders:
        begin = time.monotonic()
        op = HoleOperator(xyz, conn, state, order, a.metric == 'directional', a.scale)
        initial = reference+a.perturbation*np.maximum(0., 1-(xyz[:, 0]/1e-5)**2)
        result, history = op.solve(initial)
        denominator = max(float(np.max(abs(np.expm1(reference)))), float(np.max(abs(np.expm1(result)))))
        item = dict(order=order, metric=a.metric, scale=a.scale, iterations=history,
                    ring_abs_p_relative_mean=ring_mean(xyz, conn, result),
                    native_ring_abs_p_relative_mean=metrics['ring_abs_p_relative_mean'],
                    maximum_eta_change=float(np.max(abs(result-reference))),
                    local_density_perturbation_change=float(np.max(abs(np.expm1(result)-np.expm1(reference)))/denominator),
                    seconds=time.monotonic()-begin)
        np.savez_compressed(a.output/f'conditional_eta_p_order{order}.npz', eta_p=result)
        reports.append(item)
        save(a.output/f'order{order}.json', item)
        print(json.dumps(item), flush=True)
        del op
        gc.collect()
    if any(sha(Path(f)) != digest for f, digest in hashes.items()):
        raise RuntimeError('Source evidence changed during diagnostic')
    save(a.output/'report.json', dict(scope='CONDITIONAL_HOLE_SUBPROBLEM_DIAGNOSTIC_ONLY',
         six_other_fields_held_fixed=True, production_acceptance=False, native_gates_changed=False,
         run=str(a.run), voltage_V=metrics['voltage_V'], native_model_unchanged=True,
         source_hashes=hashes, results=reports,
         limitation='Not a seven-field solve or a mesh/SUPG certificate. Native coupled replay required.'))


if __name__ == '__main__':
    main()
