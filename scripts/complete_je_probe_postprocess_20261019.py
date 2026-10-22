"""Complete only missing postprocessing records after the UTF-8 audit recovery."""
from pathlib import Path
import json
import sys
import run_je_reconstruction_mesh_20261018 as runner

root = Path(__file__).resolve().parents[1]
output = root / 'runs/JE_RECONSTRUCTION_20261018/stage4_paper_mesh'
directory = output / 'A40_R14_Z14_S1_V060'
evidence = root / 'runs/JE_RECONSTRUCTION_20261018/paper_scope/postprocess_recovery_20261019'
assert sys.flags.utf8_mode == 1
gate = runner.validate_protocol('paper')
manifest = runner.accepted(directory)
comparison = runner.read(directory / 'mapped_field_comparison.json')
for filename, digest in comparison['hashes'].items():
    assert runner.sha(Path(filename)) == digest, filename
provenance = runner.read(directory / 'power_audit/provenance.json')
assert provenance['nodes_sha256'] == runner.sha(directory / 'nodes.csv')
assert provenance['mesh_sha256'] == manifest['mesh_sha256']
for filename, digest in provenance['scripts'].items():
    assert runner.sha(Path(filename)) == digest, filename
_, metrics = runner.classify(directory, runner.read(directory / 'execution.json')['exit_code'])
order5 = runner.read(directory / 'power_audit/order5.json')
for key, metric in [('JE_W', 'joule_W'), ('all_exterior_power_W', 'terminal_power_W')]:
    assert abs(order5[key] - metrics[metric]) / abs(metrics[metric]) < 1e-8
gap = runner.read(directory / 'power_audit/decomposition.json')
assert abs(gap['gap_W'] - (order5['JE_W'] - order5['all_exterior_power_W'])) < 1e-15
source = str(runner.baseline(6))
runner.save(directory / 'completion.json', dict(
    power_audit_completed=True, field_comparison_passed=comparison['passed'],
    fully_validated=False, source_run=source, geometry_and_transfer_only_initial_guess=True,
    recovery='UTF8 postprocessing only; original native output and failed batch record preserved'))
runner.accepted(directory, require_new_completion=True)
runner.save(evidence / 'recovery.json', dict(
    status='POSTPROCESSING_COMPLETED_LOCAL_FIELD_FAIL_RETAINED',
    utf8_mode=sys.flags.utf8_mode, native_rerun=False, native_seconds=runner.read(directory / 'execution.json')['seconds'],
    workspace_alias=str(Path.cwd()), physical_root=str(root),
    source=source, target=str(directory), gap=gap, local_comparisons=comparison['comparisons'],
    verification_scope=gate['scope'], fully_validated=False,
    wrapper_sha256=runner.sha(root / 'scripts/invoke_je_paper_batch_20261019.ps1'),
    recovery_script_sha256=runner.sha(Path(__file__)),
    completed_files={str(p): runner.sha(p) for p in [directory / 'completion.json', directory / 'mapped_field_comparison.json',
        directory / 'power_audit/decomposition.json', directory / 'nodes.csv', directory / 'input.nml']}))
print(json.dumps(dict(completed=True, local_field_passed=comparison['passed'], gap_relative=gap['gap_relative'])))
