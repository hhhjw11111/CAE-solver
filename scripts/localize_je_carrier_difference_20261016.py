"""Locate existing mapped/same-grid carrier errors without changing any gate."""
from pathlib import Path
import argparse, json
import numpy as np
from run_hv_sweep import sha


def localize(reference, target):
    reference, target = Path(reference), Path(target)
    a = np.genfromtxt(reference, delimiter=',', names=True)
    b = np.genfromtxt(target, delimiter=',', names=True)
    assert len(a) == len(b) and np.array_equal(a['node'], b['node'])
    for k in ['x_m', 'y_m', 'z_m']:
        assert np.allclose(a[k], b[k], rtol=0, atol=1e-18)
    xyz = np.column_stack([b[k] for k in ['x_m', 'y_m', 'z_m']])
    r = np.hypot(xyz[:, 0], xyz[:, 1])
    surface = np.isclose(abs(xyz[:, 2]), .5e-6, rtol=0, atol=1e-15)
    ring = r <= 3.5e-6 + 1e-15
    regions = {'all': np.ones(len(b), dtype=bool), 'top_bottom_surface': surface,
               'off_surface': ~surface, 'hole_annulus_r_le_3p5um': ring,
               'outside_annulus': ~ring, 'annulus_surface': ring & surface,
               'annulus_off_surface': ring & ~surface}
    result = {'scope': 'LOCATION_DIAGNOSTIC_ONLY_NOT_AN_ADDITIONAL_ACCEPTANCE_TEST',
              'normalization': 'maximum absolute exp(eta)-1 over both complete fields, same denominator for every region',
              'regions_overlap': True, 'rms_is_node_weighted_not_volume_weighted': True,
              'reference': str(reference.resolve()), 'target': str(target.resolve()),
              'hashes': {str(p.resolve()): sha(p) for p in [reference, target]}, 'carriers': {}}
    for eta in ['eta_n', 'eta_p']:
        pa, pb = np.expm1(a[eta]), np.expm1(b[eta])
        denominator = max(float(abs(pa).max()), float(abs(pb).max()), 1e-300)
        error = abs(pb-pa)/denominator
        rr = {}
        for name, mask in regions.items():
            ids = np.flatnonzero(mask)
            if not len(ids):
                rr[name] = {'node_count': 0}
                continue
            j = ids[np.argmax(error[ids])]
            rr[name] = {'node_count': len(ids), 'maximum_percent': float(100*error[j]),
                        'rms_percent': float(100*np.sqrt(np.mean(error[ids]**2))),
                        'maximum_node': int(b['node'][j]), 'xyz_um': (xyz[j]*1e6).tolist()}
        result['carriers'][eta] = rr
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference', type=Path, required=True)
    p.add_argument('--target', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    r = localize(args.reference, args.target)
    r['script_sha256'] = sha(Path(__file__))
    with args.output.open('x', encoding='utf-8') as f:
        json.dump(r, f, indent=2, allow_nan=False)
    print(json.dumps(r['carriers'], indent=2))
