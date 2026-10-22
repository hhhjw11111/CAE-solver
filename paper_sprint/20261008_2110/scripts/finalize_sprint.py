"""Create a local-only evidence index, protection audit and compact delivery ZIP."""
from pathlib import Path
import csv, hashlib, json, re, shutil, subprocess, zipfile
from datetime import datetime, timezone, timedelta

root = Path(__file__).resolve().parents[1]
repo = root.parents[1]
alias = Path(r'C:\Users\Administrator\Desktop\codex-workspaces\Console3\solver_v2\paper_sprint\20261008_2110')
now = datetime.now(timezone(timedelta(hours=8))).isoformat()
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''): h.update(block)
    return h.hexdigest()
def write_json(rel,obj):
    (root/rel).write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
def read_csv(rel):
    with (root/rel).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
git=shutil.which('git') or r'C:\Program Files\Git\cmd\git.exe'
def run_git(*args):
    r=subprocess.run([git,'-C',str(repo),*args],capture_output=True)
    if r.returncode: raise RuntimeError(r.stderr.decode('utf-8',errors='replace'))
    return r.stdout

before=read_csv('audit/source_hashes_before.csv')
source_changes=[r['Path'] for r in before if not Path(r['Path']).exists() or sha(r['Path']).lower()!=r['Hash'].lower()]
after_diff=run_git('diff','--binary')
(root/'audit/user_changes_after.patch').write_bytes(after_diff)
(root/'audit/git_status_after.txt').write_bytes(run_git('status','--porcelain=v1','-uall'))
protection={'checked_at':now,'head':run_git('rev-parse','HEAD').decode().strip(),
 'branch':run_git('branch','--show-current').decode().strip(),
 'tag_object':run_git('rev-parse','baseline-pre-stabilization-20260915').decode().strip(),
 'tag_commit':run_git('rev-parse','baseline-pre-stabilization-20260915^{}').decode().strip(),
 'main_ref':run_git('rev-parse','refs/heads/main').decode().strip(),
 'source_files_checked':len(before),'source_changed_since_start':source_changes,
 'tracked_worktree_diff_unchanged':after_diff==(root/'audit/user_changes_before.patch').read_bytes(),
 'history_write_scope':'No task writes to original runs/cases/source; hash identities of reused originals in postprocess registry.',
 'new_branch_only':True,'committed':False,'pushed':False}
baseline='a63ad1db24a5d5813afa4fbbe7eaa6e7235950bc'
protection['baseline_refs_unchanged']=all(protection[x]==baseline for x in ['head','tag_commit','main_ref']) and protection['tag_object']=='d01d73dd263f81c27423469c27b09ab4d28fa040'
protection['status']=('BASELINE_PASS_CONCURRENT_CHANGES_RECORDED' if source_changes or not protection['tracked_worktree_diff_unchanged'] else 'PASS') if protection['baseline_refs_unchanged'] else 'REQUIRES_REVIEW'
protection['task_source_writes']=False
protection['concurrent_change_detail']='audit/concurrent_source_changes.json'
protection['validated_source_version']='0.15.9-log-inverse-checked / R25_CLEAN'
protection['current_version_at_freeze']=(repo/'VERSION').read_text(encoding='utf-8-sig').strip()
protection['frozen_original_source_manifest']='audit/frozen_r25_source/frozen_source_manifest.json'
write_json('audit/protection_check.json',protection)

cases=read_csv('data/case_registry.csv')
metrics=read_csv('data/metrics_summary.csv')
powerfields=['case','model','integral_scope','Pd_W','Pdiff_W','Pdiff_abs_W','P_JdotE_W','Pterm_all_boundary_W','Pterm_physical_exterior_W','heat_out_physical_exterior_W','power_identity_gap_relative']
write_json('data/power_audit.json',{'scope':'3V archived fields; volume integrals recomputed; boundary integrals reused','rows':[{k:r[k] for k in powerfields} for r in metrics]})
write_json('data/spatial_comparison.json',{'case_family':'3V R23 LOG_WEIGHTED y-half L1_to_L2','source':'data/spatial_failure_source.csv','status':'FAIL_LOCAL_CARRIER_2_PERCENT','rows':read_csv('data/spatial_failure_source.csv')})
verification=json.loads((root/'validation/runtime_retry/verification_summary.json').read_text(encoding='utf-8-sig'))
evidence=[
 ('E001','frozen startup R25 source/build identity','PASS','audit/source_validation/source_audit_runtime_retry.json','58 source + executable matched at startup; excludes later concurrent 0.16 edits; not clean worktree'),
 ('E002','historical per-case layered status','PARTIAL','data/case_registry.csv','H0 contemporaneous provenance PENDING; all spatial limitations retained'),
 ('E003','fresh small verification of frozen executable','PASS_WITH_SCOPE','validation/runtime_retry/verification_summary.json','five frozen R25 CLI checks plus R22 archived kernel replay; not full 1D production specification'),
 ('E004','mathematical reference','PASS_MATH_ONLY','validation/formula_reference/formula_reference_summary.json','not a production/FEM verification'),
 ('E005','3V volume postprocessing','PASS_POSTPROCESS_ONLY','data/metrics_summary.csv','L0 exploratory; boundary integrals archived'),
 ('E006','local carrier spatial acceptance','FAIL','data/spatial_failure_source.csv','R23 half domain; 2 percent gate'),
 ('E007','literature','VERIFIED_WITH_READ_SCOPE','audit/literature/reference_registry.csv','two papers; metadata/abstract/indexed excerpts only'),
 ('E008','baseline refs and concurrent edits',protection['status'],'audit/protection_check.json','HEAD/main/tag unchanged; external 0.16 source changes preserved; original 58 sources recovered by exact hash'),
 ('E009','remaining tasks','PENDING','OPEN_ITEMS.md','no unrun cases or unknown values promoted to PASS')]
with (root/'audit/evidence_registry.csv').open('w',encoding='utf-8',newline='') as f:
    w=csv.writer(f);w.writerow(['evidence_id','description','status','source','claim_boundary']);w.writerows(evidence)
write_json('audit/evidence_register.json',{'created_at':now,'task_scope':'report first; 3V H0/H1 radius3um; no added production solve','evidence':[dict(zip(['id','description','status','source','claim_boundary'],r)) for r in evidence],'cases':cases,'verification':verification,'claims_file':'audit/claim_matrix.md'})
source=Path(r'C:\Users\Administrator\Downloads\CAE_solver_Codex_action_report_20261008.md')
write_json('audit/action_attachment_identity.json',{'path':str(source),'sha256':sha(source),'bytes':source.stat().st_size,'role':'User-authorized task plan; not numerical evidence'})

required=['reports/study_report_final.md','reports/study_report_final.pdf','manuscript/manuscript_working_draft.md','manuscript/manuscript_working_draft.pdf','OPEN_ITEMS.md','CHANGELOG.md','analysis/postprocess_report.md','audit/manuscript_pdf_qa.md','audit/render_qa.md']
missing=[p for p in required if not (root/p).exists() or (root/p).stat().st_size==0]
if missing:raise RuntimeError('Missing delivery files: '+str(missing))
broken=[]
for rel in ['reports/study_report_final.md','manuscript/manuscript_working_draft.md']:
    p=root/rel; txt=p.read_text(encoding='utf-8')
    for target in re.findall(r'!\[[^\]]*\]\(([^)]+)\)',txt):
        dest=Path(target) if re.match(r'^[A-Za-z]:',target) else p.parent/target
        if not dest.exists():broken.append(str(dest))
write_json('audit/delivery_structure_check.json',{'missing_required':missing,'broken_document_images':broken,'required_files':required,'status':'PASS' if not missing and not broken else 'FAIL'})
if broken:raise RuntimeError('Broken images')
manifest={'started_at':'2026-10-08T21:10+08:00','start_time_precision':'minute','frozen_at':now,'timezone':'Asia/Shanghai','delivery_root':str(alias),'physical_root':str(root),'branch':protection['branch'],'head':protection['head'],'protected_baseline_commit':baseline,'new_3d_solves':0,'processed_3V_cases':len(cases),'figure_sets':len(list((root/'figures').glob('*.png'))),'figure_formats':['PNG','SVG','PDF'],'data_csv_count':len(list((root/'data').rglob('*.csv'))),'document_status':'WORKING_DRAFT','formats':{'markdown':'DELIVERED','pdf':'DELIVERED_VISUAL_QA','docx':'PENDING_RENDER_ENV'},'source_modified_by_this_task':False,'concurrent_changed_source_count':len(source_changes),'validated_version':protection['validated_source_version'],'current_version_not_validated_here':protection['current_version_at_freeze'],'initial_source_files':len(before),'readme':'README_FIRST.md','open_items':'OPEN_ITEMS.md','protection_status':protection['status'],'raw_data_archived_in_place':True,'no_external_upload':True,'verification_scope':'five R25 small CLI checks + R22 archived production kernel replay + pure mathematical references'}
write_json('audit/run_manifest.json',manifest)

files=[]
for p in sorted(root.rglob('*')):
    if not p.is_file():continue
    rel=p.relative_to(root).as_posix()
    if rel.startswith('audit/render_qa/') or p.suffix=='.zip' or '__pycache__' in rel or rel=='audit/delivery_files_sha256.csv':continue
    files.append((rel,p.stat().st_size,sha(p)))
with (root/'audit/delivery_files_sha256.csv').open('w',encoding='utf-8',newline='') as f:
    w=csv.writer(f);w.writerow(['relative_path','bytes','sha256']);w.writerows(files)
zpath=root/'CAE_solver_paper_sprint_20261008.zip'
with zipfile.ZipFile(zpath,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for rel,_,_ in files:
        if rel.endswith('.patch') or rel=='audit/report_render_probe.pdf':continue
        z.write(root/rel,rel)
    z.write(root/'audit/delivery_files_sha256.csv','audit/delivery_files_sha256.csv')
with zipfile.ZipFile(zpath,'r') as z:
    corrupt=z.testzip()
    if corrupt:raise RuntimeError('ZIP CRC failure: '+corrupt)
    packed_sha={r:hashlib.sha256(z.read(r)).hexdigest() for r in z.namelist()}
    mismatch=[r for r,_,h in files if r in packed_sha and packed_sha[r]!=h]
    if mismatch:raise RuntimeError('ZIP manifest hash mismatch: '+str(mismatch))
    packed_count=len(packed_sha)
print(json.dumps({'protection':protection,'manifest':manifest,'zip_crc_status':'PASS','zip_member_hash_status':'PASS','zip_file_count':packed_count,'zip_bytes':zpath.stat().st_size,'zip_sha256':sha(zpath)},ensure_ascii=False))
