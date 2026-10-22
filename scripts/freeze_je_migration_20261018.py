"""Freeze the approved migration protocol before new native verification runs."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'runs/JE_RECONSTRUCTION_20261018'
ARCHIVE = ROOT.parent / 'solver_backups/HEAT_MIGRATION_20261018'

def sha(path):
    return hashlib.file_digest(path.open('rb'), 'sha256').hexdigest()

def main():
    OUT.mkdir(exist_ok=True)
    contract_path = OUT / 'contract.json'
    if contract_path.exists():
        raise SystemExit('Refusing to overwrite a frozen contract')
    archive = json.loads((ARCHIVE / 'verification.json').read_text(encoding='utf-8'))
    assert archive['state'] == 'COMPLETE' and archive['all_archive_members_sha256_verified']
    sources = json.loads((ARCHIVE / '02_source_build_and_protocols_manifest.json').read_text(encoding='utf-8'))
    core = [r for r in sources if r['path'].startswith('solver_v2/src/') and '/verification/' not in r['path']]
    gates = ['scripts/audit_hv_mesh.py','scripts/audit_hv_field_shape.py',
             'scripts/audit_hv_stabilization.py','scripts/audit_hv_half_y.py',
             'scripts/audit_je_sensitivity_power_20261015.py']
    data = dict(model='JE_SIMPLIFIED_R0_V1', user_order=[1,3,4,5], skipped_step=2,
        archive=str(ARCHIVE), archive_verification_sha256=sha(ARCHIVE/'verification.json'),
        protocol='docs/je_model_migration_20261018.md',
        protocol_sha256=sha(ROOT/'docs/je_model_migration_20261018.md'),
        model_scope='Signed E dot (Jn+Jp), R=0, simplified heat model; not full thermoelectric closure',
        legacy_default_changed=False, original_results_retained=True,
        core_sources_before=core,
        frozen_gate_scripts={p:sha(ROOT/p) for p in gates},
        resources=dict(native_concurrency=1,threads=8,ooc_initial_MB=6000,hard_timeout_seconds=7200),
        stage3=dict(ids=['A1','A2','A1JJE','JEGRAD','A2P','G02JE','F1MJE','G04JE','SUPG','SUPGJ','SUPGJT'],
                    legacy_cli='16 byte-identical comparisons',
                    jegrad_field_power_tolerance=1e-8,jegrad_quadrature_tolerance=1e-10,
                    a2p_field_L2=1e-3,a2p_gradient_L2=.02,a2p_hole_normal=.01),
        stage4=dict(first_probe='FULL A40 R14 Z14 SUPG1 6V',
                    angular=[32,40,48],radial=[18,20,22],voltages=[6,8],supg=[.5,1,2],
                    local_carrier_tolerance=.02,supg_scalar_ring_tolerance=.01,
                    power_gap_arbitrary_new_threshold=None,requires_joint_mesh=True),
        stage5_gate='Requires stage3 and stage4 acceptance; do not promote failed candidates',
        statuses=dict(stage1='FROZEN',stage3='PENDING',stage4='PENDING',stage5='GATED'))
    contract_path.write_text(json.dumps(data,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(contract_path)

if __name__ == '__main__':
    main()
