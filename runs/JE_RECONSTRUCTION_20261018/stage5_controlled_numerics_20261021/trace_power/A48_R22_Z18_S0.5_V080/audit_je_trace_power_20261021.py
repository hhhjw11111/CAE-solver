"""Independent power identities using the explicitly selected trace SUPG metric.

The historical power auditor uses directional tau. Its physical-current JE and
boundary integrals remain valid for trace runs, but its S/reaction decomposition
must not be reused. This private copy changes only gdir=trace, records both code
hashes, and writes to a new directory. Never correct the heat by the power gap.
"""
from pathlib import Path
import argparse
import importlib.util
import json
import shutil
import numpy as np
import run_je_reconstruction_mesh_20261018 as r
import diagnose_je_power_gap_20261015 as original


def trace_module(output, scale):
    source = Path(original.__file__)
    text = source.read_text(encoding='utf-8')
    old = 'gdir=(advmetric+eps*trace/3)/(speed2+eps)'
    if text.count(old) != 1:
        raise RuntimeError('Historical power auditor changed; explicit review needed')
    copied = output/'trace_power_kernel.py'
    copied.write_text(text.replace(old, 'gdir=trace'), encoding='utf-8')
    spec = importlib.util.spec_from_file_location('_je_trace_power_'+str(scale), copied)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    sample = module.sample

    def scaled_sample(*args, **kwargs):
        result = sample(*args, **kwargs)
        result['S'] *= scale
        return result

    module.sample = scaled_sample
    return module, copied


def audit(directory, output):
    directory, output = Path(directory).resolve(), Path(output).resolve()
    manifest, metrics = r.read(directory/'manifest.json'), r.read(directory/'execution.json')
    status, metrics = r.classify(directory, metrics['exit_code'])
    if (status != r.PASS or metrics['supg_directional_tau'] != 'F'
            or manifest.get('supg_directional_tau') is not False
            or metrics['heat_source'] != 'E_DOT_TOTAL_CURRENT_DIAGNOSTIC_NOT_THERMODYNAMIC_CLOSURE'):
        raise RuntimeError('Completed trace-metric JE state required')
    for kind in ('mesh', 'executable'):
        if r.sha(Path(manifest[kind])) != manifest[kind+'_sha256']:
            raise RuntimeError('Changed '+kind)
    output.mkdir(parents=True, exist_ok=False)
    kernel, copied = trace_module(output, manifest['supg_scale'])
    geometry, conn = kernel.fe.read_mesh(Path(manifest['mesh']))
    nodes = np.loadtxt(directory/'nodes.csv', delimiter=',', skiprows=1)
    if np.max(abs(nodes[:,1:4]-geometry)) >= 1e-19:
        raise RuntimeError('Field and geometry differ')
    _, lookup = kernel.topology(geometry, conn)
    sources = [Path(__file__), Path(original.__file__), copied, Path(kernel.fe.__file__)]
    sources += [directory/f for f in ('nodes.csv', 'manifest.json', 'execution.json', 'metrics.txt', 'input.nml')]
    provenance = dict(scope='TRACE_SUPG_POWER_AUDIT', supg_directional_tau=False,
        supg_scale=manifest['supg_scale'], verification=kernel.basic_checks(),
        source_hashes={str(p):r.sha(p) for p in sources}, heat_changed=False,
        historical_directional_decomposition_applicable=False)
    r.save(output/'provenance.json', provenance)
    results = {order:kernel.audit_order(geometry, conn, nodes[:,4:11], lookup, order, output) for order in (3,5)}
    a, b = results[3], results[5]
    for name, native in (('JE_W','joule_W'), ('all_exterior_power_W','terminal_power_W')):
        error = abs(b[name]-metrics[native])/abs(metrics[native])
        if error >= 1e-8:
            raise RuntimeError('Independent trace power differs from native metric: '+name)
    cut = b['boundary_sets']['y-']['power_W']
    physical = b['all_exterior_power_W']-cut
    relative = abs(physical-metrics['physical_exterior_terminal_power_W'])/abs(physical)
    if relative >= 1e-8:
        raise RuntimeError('Physical exterior mismatch')
    parts = dict(quadrature_W=b['JE_W']-a['JE_W'], SUPG_weighted_W=a['phi_weighted_SUPG_W'],
        reaction_minus_all_boundary_W=a['reaction_terminal_power_W']-b['all_exterior_power_W'],
        free_residual_correction_W=-a['phi_weighted_free_residual_W'])
    r.save(output/'report.json', dict(scope='TRACE_POWER_DIAGNOSTIC_NO_ADDED_THRESHOLD',
        voltage_V=metrics['voltage_V'], raw_half_JE_W=b['JE_W'],
        raw_half_physical_exterior_power_W=physical, raw_half_symmetry_cut_power_W=cut,
        reconstructed_full_JE_W=2*b['JE_W'],
        physical_exterior_gap_relative=(b['JE_W']-physical)/b['JE_W'],
        independent_vs_native_relative=relative, original_comparison_limit=1e-8,
        assembled_order3_free_residual_max_A=a['free_residual_max_A'],
        assembled_order3_identity_remainder_W=a['algebraic_identity_remainder_W'],
        parts=parts, decomposition_remainder_W=b['gap_W']-sum(parts.values()),
        power_gap_acceptance_threshold=None, production_acceptance=False))
    shutil.copy2(__file__, output/Path(__file__).name)
    print(json.dumps(r.read(output/'report.json')), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    audit(args.run, args.output)
