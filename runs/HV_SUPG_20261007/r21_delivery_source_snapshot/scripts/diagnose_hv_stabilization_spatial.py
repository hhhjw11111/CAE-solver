"""Locate differences between aligned quadrature samples; never certify a mesh."""
import argparse
import bisect
import csv
import itertools
import json
import math
from pathlib import Path
from run_hv_sweep import classify, sha

EDGES_NM = [0., 2., 5., 10., 20., 50., 100., 250., 500.000001]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runs', type=Path, nargs=2, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    manifests = []
    for path in args.runs:
        execution = json.loads((path/'execution.json').read_text())
        state, metrics = classify(path, execution['exit_code'])
        if state != 'SOLVE_FLUX_PASS_MESH_PENDING':
            raise ValueError('Both states must pass their existing solve and flux gates')
        manifests.append(json.loads((path/'manifest.json').read_text()))
    for key in ['mesh_sha256', 'executable_sha256', 'voltages', 'model']:
        if manifests[0][key] != manifests[1][key]:
            raise ValueError('States do not have the same mesh/executable/voltage/model')
    for manifest in manifests:
        if sha(Path(manifest['mesh'])) != manifest['mesh_sha256']:
            raise ValueError('Mesh hash mismatch')
    bins = {}
    totals = {}
    files = [p/'field_samples.csv' for p in args.runs]
    count = 0
    with files[0].open(newline='') as fa, files[1].open(newline='') as fb:
        ra, rb = csv.reader(fa), csv.reader(fb)
        if next(ra) != next(rb):
            raise ValueError('Field layouts differ')
        for a, b in itertools.zip_longest(ra, rb):
            if a is None or b is None or a[:4] != b[:4]:
                raise ValueError('Quadrature coordinates/weights/count differ')
            x, y, z, weight = map(float, a[:4])
            if weight <= 0:
                raise ValueError('Nonpositive volume weight')
            pa, pb, na, nb = float(a[10]), float(b[10]), float(a[9]), float(b[9])
            if not all(math.isfinite(v) and v > 0 for v in [pa, pb, na, nb]):
                raise ValueError('Invalid sampled concentration')
            depth_nm = (.5e-6-abs(z))*1e9
            index = bisect.bisect_right(EDGES_NM, depth_nm)-1
            if not 0 <= index < len(EDGES_NM)-1:
                raise ValueError('Sample outside the frozen thickness')
            zone = 'hole_annulus' if math.hypot(x,y) < 3.5e-6 else (
                'electrode_strip' if abs(x) > 9.5e-6 else 'bulk')
            contribution = [weight, weight*abs(pa-1), weight*abs(pb-1), weight*abs(pa-pb),
                            weight*abs(na-1), weight*abs(nb-1), weight*abs(na-nb)]
            for region in ['all', zone]:
                key = (region, 'bottom' if z < 0 else 'top', index)
                values = bins.setdefault(key, [0.]*len(contribution))
                total = totals.setdefault(region, [0.]*len(contribution))
                for k, value in enumerate(contribution):
                    values[k] += value
                    total[k] += value
            count += 1
    def metrics(values):
        return dict(volume_m3=values[0], mean_abs_p_relative=[values[1]/values[0], values[2]/values[0]],
                    mean_abs_p_difference=values[3]/values[0],
                    mean_abs_n_relative=[values[4]/values[0], values[5]/values[0]],
                    mean_abs_n_difference=values[6]/values[0])
    rows = []
    for (region, side, index), values in sorted(bins.items()):
        row = metrics(values)
        row.update(region=region, side=side, surface_depth_nm=EDGES_NM[index:index+2],
                   fraction_of_region_p_difference=values[3]/max(totals[region][3], 1e-300))
        rows.append(row)
    report = dict(scope='SPATIAL_DIFFERENCE_DIAGNOSTIC_NOT_MESH_ACCEPTANCE',
                  sample_count=count, supg_scales=[m['supg_scale'] for m in manifests],
                  sources=[dict(path=str(p.absolute()), sha256=sha(p)) for p in files],
                  region_totals={k:metrics(v) for k,v in totals.items()}, bins=rows)
    with args.output.open('x', encoding='utf-8') as f:
        json.dump(report, f, indent=2, allow_nan=False)
    print(json.dumps(report['region_totals'], indent=2))
    print('Largest annular contributions:', json.dumps(sorted(
        [r for r in rows if r['region']=='hole_annulus'],
        key=lambda r:r['fraction_of_region_p_difference'], reverse=True)[:5], indent=2))


if __name__ == '__main__':
    main()
