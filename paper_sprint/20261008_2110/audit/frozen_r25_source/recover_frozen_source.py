"""Recover only byte-identical startup sources from existing immutable snapshots.
Never reads current src/ as a recovery source; never rewrites existing destinations.
"""
from pathlib import Path
import csv, hashlib, json, shutil
from datetime import datetime, timezone

HERE=Path(__file__).resolve().parent
AUDIT=HERE.parent
SOLVER=HERE.parents[3]
SNAPSHOTS=[SOLVER/'runs/HV_SUPG_20261007'/s for s in [
    'r25_checkpoint_code_snapshot','r25_followup_code_snapshot',
    'r25_autocontinue_code_snapshot','r25_source_snapshot',
    'r25_mesh_execution_audit_snapshot','r23_source_snapshot']]

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    expected_file=AUDIT/'source_hashes_before.csv'
    with expected_file.open(encoding='utf-8-sig',newline='') as stream:expected=list(csv.DictReader(stream))
    destroot=HERE/'src'
    if destroot.exists():raise FileExistsError('Frozen destination already exists; no overwrite: '+str(destroot))
    destroot.mkdir()
    records=[]
    for row in expected:
        original=row['Path'];relative=original.replace('\\','/').split('/solver_v2/src/',1)[1]
        required=row['Hash'].lower();checked=[];match=None
        for root in SNAPSHOTS:
            candidate=root/'src'/relative
            if not candidate.is_file():continue
            candidate_hash=sha(candidate)
            checked.append(dict(path=str(candidate),sha256=candidate_hash,matched=candidate_hash==required))
            if candidate_hash==required and match is None:match=candidate
        record=dict(original_startup_path=original,relative_path='src/'+relative,expected_sha256=required,candidates_checked=checked)
        if match is None:record.update(status='MISSING_EXACT_HASH_MATCH',recovered_from=None)
        else:
            destination=destroot/relative;destination.parent.mkdir(parents=True,exist_ok=True)
            # Recheck immediately before copy to guard concurrent snapshot changes.
            if sha(match)!=required:raise RuntimeError('Recovery source changed: '+str(match))
            shutil.copy2(match,destination)
            actual=sha(destination)
            if actual!=required:raise RuntimeError('Copy hash mismatch: '+str(destination))
            record.update(status='RECOVERED_EXACT_HASH_MATCH',recovered_from=str(match),destination=str(destination),copied_sha256=actual,bytes=destination.stat().st_size)
        records.append(record)
    recovered=sum(r['status']=='RECOVERED_EXACT_HASH_MATCH' for r in records)
    manifest=dict(schema='TPSolverV2-startup-source-recovery-v1',created_utc=datetime.now(timezone.utc).isoformat(),
        criterion='Exact SHA-256 equality to source_hashes_before.csv; no use of current changing src/ files',
        expected_manifest=str(expected_file),expected_manifest_sha256=sha(expected_file),snapshot_search_roots=[str(s) for s in SNAPSHOTS],
        expected_count=len(records),recovered_count=recovered,missing_count=len(records)-recovered,
        status='COMPLETE_EXACT_STARTUP_SOURCE_RECOVERY' if recovered==len(records) else 'PARTIAL_EXACT_STARTUP_SOURCE_RECOVERY',
        compilation_performed=False,verification_performed=False,current_source_modified=False,git_refs_modified=False,
        limitation='Preserves startup source bytes only; not a claim of rebuilt executable equivalence or new numerical validation.',files=records)
    (HERE/'frozen_source_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:manifest[k] for k in ['expected_count','recovered_count','missing_count','status']}))
    if recovered!=len(records):print('MISSING',json.dumps([r['relative_path'] for r in records if r['status'].startswith('MISSING')]))
if __name__=='__main__':main()
