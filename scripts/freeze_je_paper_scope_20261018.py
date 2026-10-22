"""Record the explicitly approved minimal JE paper scope without editing prior evidence."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STUDY = ROOT / 'runs/JE_RECONSTRUCTION_20261018'


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def main():
    folder = STUDY / 'paper_scope'
    folder.mkdir(exist_ok=True)
    workflow = ROOT / 'docs/je_paper_minimal_workflow_20261018.md'
    conditions = [dict(radius_um=3, voltage_V=v, branch='NT') for v in (6, 6.5, 7, 7.5, 8)]
    conditions += [dict(radius_um=r, voltage_V=v, branch='NT') for r in (2, 4) for v in (6, 8)]
    scope = dict(
        scope='JE_PAPER_MINIMAL_20261018', user_approved=True,
        authorization='好的 按这个流程执行，原A1、A2和旧数据不动，先核对复用已有验证，再补齐6V、8V的网格和SUPG检查，通过后完成9个生产工况和论文图表。',
        client_date='2026-10-18',
        required_ids=['A1','A2','A1JJE','JEGRAD','G02JE','F1MJE','G04JE','SUPG','SUPGJ','SUPGJT'],
        required_legacy_ids=['V01','V02','V03','V04','V05','V06','G02','G03','G04','G05','G06','G07','G08','L01','F1M','A1J'],
        excluded_ids=['A2P'],
        exclusion_reason='User explicitly retained original A1/A2 and removed analytic extensions from the minimal paper scope; original A2P and original full-suite FAIL remain unchanged.',
        original_contract=str(STUDY/'contract.json'), original_contract_sha256=sha(STUDY/'contract.json'),
        approved_workflow=str(workflow), approved_workflow_sha256=sha(workflow),
        source_attempt=str(STUDY/'stage3_verification/attempt01'),
        archive=str(ROOT.parent/'solver_backups/HEAT_MIGRATION_20261018'),
        model='JE_SIMPLIFIED_R0_V1', heat_source='signed E dot (Jn+Jp)',
        legacy_default_changed=False, original_A1_A2_changed=False,
        required_thresholds_relaxed=False, physics_and_linear_core_changes_authorized=False,
        verification_method='Recompute identities and evidence checks; reuse unchanged passed runs, not a new numerical replay.',
        verification_output=str(folder/'verification_reuse'),
        mesh_output=str(STUDY/'stage4_paper_mesh'),
        production_conditions=conditions,
        production_gate='Requires applicable mesh, local-field, SUPG, native and independent-flux acceptance; each radius requires applicable spatial evidence.',
        resource_policy=dict(native_concurrency=1,threads=8,ooc_initial_MB=6000,
                             hard_timeout_seconds=7200,start_disk_GiB=16,stop_disk_GiB=4,
                             stop_available_RAM_MiB=512,low_RAM_consecutive_samples=3),
        preserve='All original inputs, scientific data, failed runs, protocols and builds; cleanup only owned exited new-run reproducible factors with inventory.',
        scope_writer_sha256=sha(Path(__file__)))
    with (folder/'contract.json').open('x',encoding='utf-8') as stream:
        json.dump(scope,stream,indent=2,ensure_ascii=False)
        stream.write('\n')
    print(folder/'contract.json')


if __name__ == '__main__':
    main()
