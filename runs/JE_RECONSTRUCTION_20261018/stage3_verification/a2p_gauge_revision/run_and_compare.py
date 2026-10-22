"""One preregistered gauge-placement experiment; retain every failed record."""
import csv
import json
import math
import sys
from pathlib import Path

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[3]
sys.path.insert(0, str(ROOT / 'scripts'))
from run_je_reconstruction_verification_20261018 import sha, write_json, check_build, run_check

EXE = ROOT / 'build/JE_VERIFY_20261018_GAUGE/TPSolverV2.exe'
OLD = OUT.parent / 'a2p_diagnostic/A2P_fields'
NREF, CHARGE, MU_SUM, HALF_WIDTH, E0 = 1e21, 1.602176634e-19, .025, 1e-5, 1e4
SIGMA = CHARGE * NREF * MU_SUM

def rows(path):
    with path.open(newline='', encoding='utf-8-sig') as stream:
        return [{k: float(v) for k, v in row.items()} for row in csv.DictReader(stream)]

def main():
    before = {str(p.relative_to(ROOT)): sha(p) for p in (ROOT / 'src').rglob('*.f90')}
    build = check_build(EXE)
    contract = json.loads((OUT.parents[1] / 'contract.json').read_text(encoding='utf-8'))
    core = [dict(path=r['path'], before=r['sha256'], after=sha(ROOT.parent / r['path']))
            for r in contract['core_sources_before']]
    provenance = dict(executable=str(EXE), executable_sha256=sha(EXE),
                      build_manifest_sha256=sha(EXE.parent / 'build_manifest.json'),
                      protocol_sha256=sha(OUT.parent / 'gauge_revision.json'),
                      sources=before, candidate_build=build, core=core,
                      core_unchanged=all(r['before'] == r['after'] for r in core))
    write_json(OUT / 'provenance.json', provenance)
    if build['source_mismatches'] or not build['executable_matches'] or not provenance['core_unchanged']:
        raise SystemExit('Provenance mismatch; no native run')
    native, _ = run_check(EXE, 'A2P', OUT, 'gauge', 1800, exports=True)
    write_json(OUT / 'native_result.json', native)
    comparisons = []
    for new_file in sorted((OUT / 'A2P_fields').glob('*_nodes.csv')):
        name = new_file.stem.removeprefix('A2P_').removesuffix('_nodes')
        mode, sign, level = name.split('_')
        e = {'epos': 1., 'eneg': -1., 'ezero': 0.}[sign]
        a, b = (e * E0 / 1e11, -E0) if mode == 'electric' else (1e-6, e * 1e-6 / 1e-10)
        wscale, pscale = max(abs(a) * HALF_WIDTH, 1e-12), max(abs(b) * HALF_WIDTH, .1)
        qscale = SIGMA * (b * b if b else E0 * E0)
        old_nodes, new_nodes = rows(OLD / new_file.name), rows(new_file)
        same_geometry = len(old_nodes) == len(new_nodes) and all(
            all(o[k] == n[k] for k in ('x_m', 'y_m', 'z_m')) for o, n in zip(old_nodes, new_nodes))
        if not same_geometry:
            raise RuntimeError(f'Unmatched node geometry: {name}')
        offset = new_nodes[0]['phi_V'] - old_nodes[0]['phi_V']
        changes = dict(
            w=max(abs(n['w_m'] - o['w_m']) / wscale for o, n in zip(old_nodes, new_nodes)),
            phi_mod_constant=max(abs(n['phi_V'] - o['phi_V'] - offset) / pscale
                                 for o, n in zip(old_nodes, new_nodes)),
            n=max(abs(math.expm1(n['eta_n']) - math.expm1(o['eta_n'])) for o, n in zip(old_nodes, new_nodes)),
            p=max(abs(math.expm1(n['eta_p']) - math.expm1(o['eta_p'])) for o, n in zip(old_nodes, new_nodes)))
        wall_name = new_file.name.replace('_nodes.csv', '_wall.csv')
        old_wall, new_wall = rows(OLD / wall_name), rows(OUT / 'A2P_fields' / wall_name)
        same_wall = len(old_wall) == len(new_wall) and all(
            all(o[k] == n[k] for k in ('angle_rad', 'x_m', 'y_m')) for o, n in zip(old_wall, new_wall))
        if not same_wall:
            raise RuntimeError(f'Unmatched wall sampling: {name}')
        for column in ('tz_scaled', 'Dn_scaled', 'Jn_scaled', 'Jp_scaled'):
            changes[column] = max(abs(n[column] - o[column]) for o, n in zip(old_wall, new_wall))
        changes['Q'] = max(abs(n['Q_W_m3'] - o['Q_W_m3']) / qscale for o, n in zip(old_wall, new_wall))
        comparisons.append(dict(case=name, geometry_identical=True, phi_constant_offset_V=offset,
                                maximum_normalized_changes=changes,
                                passed=all(v < 1e-9 for v in changes.values())))
    stable = all(sha(ROOT / p) == h for p, h in before.items())
    result = dict(native=native, comparisons=comparisons, comparison_count=len(comparisons),
                  equivalence_passed=len(comparisons) == 18 and all(c['passed'] for c in comparisons),
                  maximum_equivalence_change=max((max(c['maximum_normalized_changes'].values())
                                                  for c in comparisons), default=None),
                  source_still_same=stable, core_sources_unchanged=provenance['core_unchanged'])
    result['all_passed'] = native['passed'] and result['equivalence_passed'] and stable
    write_json(OUT / 'result.json', result)
    write_json(OUT / 'artifact_hashes.json', {str(p.relative_to(OUT)): sha(p)
               for p in OUT.rglob('*') if p.is_file() and p.name != 'artifact_hashes.json'})
    print(json.dumps({k: v for k, v in result.items() if k not in ('native', 'comparisons')}), flush=True)
    raise SystemExit(0 if result['all_passed'] else 1)

if __name__ == '__main__':
    main()
