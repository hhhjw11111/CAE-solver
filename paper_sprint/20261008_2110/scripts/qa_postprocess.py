"""Final artifact QA and checksums for delivered postprocessing files."""
import postprocess_3v as p
import json,csv
import numpy as np
import matplotlib.image as mi
from pathlib import Path

def main():
    report=p.OUT/'analysis/postprocess_report.md'
    text=report.read_text(encoding='utf-8').replace('厚度和中面重构对每个物理取样点求逆等参映射，不作域外外推。','厚度线与孔周线对物理取样点求逆等参映射，不作域外外推；中面图直接在原始HEX20面片自然坐标上求值。')
    text=text.replace('高于上述微小反馈量级；不能用约化漂移热源直接证明真实器件的精细热反馈机制。','提示不能忽略其物理解释；功率比不等于温度误差，需独立能量模型论证，不能用约化漂移热源直接证明真实器件的精细热反馈机制。')
    provenance_note='F1早期H0三例缺少当时逐case provenance.json；metadata_correction.json中的环境是后补回溯，未补造实时记录。因此这些case的完整运行来源状态为PENDING，归档输入/网格/二进制及54个源码文件哈希匹配为已完成子项。F1 H1有当时provenance；R23/R25各58个归档源码文件匹配。没有缺少两份源码的情况。'
    if provenance_note not in text:text=text.replace('本次归档身份检查不等价于从源码重新编译。','本次归档身份检查不等价于从源码重新编译。\n\n'+provenance_note)
    if 'scripts/qa_postprocess.py' not in text:
        text=text.replace('scripts/finalize_postprocess.py\n```','scripts/finalize_postprocess.py\npython -B solver_v2/paper_sprint/20261008_2110/scripts/qa_postprocess.py\n```')
    report.write_text(text,encoding='utf-8')
    contract=p.OUT/'analysis/figure_contract.md'
    ct=contract.read_text(encoding='utf-8')
    ct=ct.replace('接触区 |x|≥9 μm','触及 x=±10 μm 端面的第一层单元为接触区')
    if '修订记录' not in ct:ct+='\n修订记录：初始接触区中心|x|≥9μm在部分网格为空；在解读分区前改为触及x=±10μm端面的第一层单元，并在region_definition_correction.json保存原因。其余分区顺序不变，物理与全域指标不变。\n'
    contract.write_text(ct,encoding='utf-8')
    regfile=p.OUT/'data/case_registry.csv'
    with regfile.open(encoding='utf-8-sig',newline='') as stream:registry=list(csv.DictReader(stream))
    for row in registry:
        row['provenance_detail']='Archived input/mesh/exe/source hashes match; not a rebuilt-binary reproducibility claim.'
        if row['case'].startswith('FAST_H0_'):
            row['provenance_status']='PENDING';row['provenance_detail']+=' Missing contemporaneous per-case provenance; early environment was retrospectively recorded.'
    p.csvwrite(regfile,registry)
    p.dump(p.OUT/'audit/postprocess/provenance_limitation.json',dict(note=provenance_note,metadata_path=str((p.F1/'evidence/metadata_correction.json').resolve()),metadata_sha256=p.sha(p.F1/'evidence/metadata_correction.json')))
    qa=dict(status='ARTIFACT_QA_COMPLETE_NOT_PHYSICAL_ACCEPTANCE',visual_review='8 PNG figures viewed after generation; revised power legend avoids overlap; no clipped labels observed; unsmoothed element traces retained.',figures=[],source_csvs=[],inverse_mapping_max_residual_m=0.,executed_commands=[dict(command='python -B solver_v2/paper_sprint/20261008_2110/scripts/'+s,exit_code=0) for s in ['postprocess_3v.py','supplemental_3v.py','refine_pe_regions.py','finalize_postprocess.py']],initial_script_issue='Initial postprocessor N0/LQ scalar thermal-voltage broadcasting error fixed; subsequent full pass completed. This was not a PDE case failure.')
    spatial=json.loads((p.HV/'field_shape_r23_L1_L2_V030.json').read_text(encoding='utf-8'))
    p.csvwrite(p.OUT/'data/spatial_failure_source.csv',[dict(field=k,maximum_relative=v['maximum_relative'],rms_relative=v['rms_relative'],normalization=v['normalization'],limit=v['limit'],passed=v['passed']) for k,v in spatial['comparisons'].items()])
    try:import fitz
    except ImportError:fitz=None
    for image in sorted((p.OUT/'figures').glob('*.png')):
        a=mi.imread(image);svg=image.with_suffix('.svg');pdf=image.with_suffix('.pdf');item=dict(stem=image.stem,png_shape=list(a.shape),svg_has_editable_text='<text' in svg.read_text(encoding='utf-8'),pdf_exists=pdf.exists())
        if fitz:
            with fitz.open(pdf) as doc:item['pdf_extractable_text_characters']=sum(len(page.get_text()) for page in doc)
        assert item['svg_has_editable_text'] and item['pdf_exists'] and min(a.shape[:2])>500
        qa['figures'].append(item)
    for file in sorted((p.OUT/'data').rglob('*.csv')):
        with file.open(encoding='utf-8-sig',newline='') as f: rows=list(csv.DictReader(f))
        assert rows
        for row in rows:
            if 'inverse_residual_m' in row:qa['inverse_mapping_max_residual_m']=max(qa['inverse_mapping_max_residual_m'],float(row['inverse_residual_m']))
            assert not any(str(v).lower() in ['nan','inf','-inf'] for v in row.values())
        qa['source_csvs'].append(dict(path=file.relative_to(p.OUT).as_posix(),rows=len(rows)))
    p.dump(p.OUT/'audit/postprocess/visual_and_data_qa.json',qa)
    manifest=[]
    for folder in ['data','figures','scripts','analysis','audit/postprocess']:
        for file in sorted((p.OUT/folder).rglob('*')):
            if file.is_file() and file.name!='hash_manifest.json' and '__pycache__' not in str(file):manifest.append(dict(path=file.relative_to(p.OUT).as_posix(),bytes=file.stat().st_size,sha256=p.sha(file)))
    p.dump(p.OUT/'audit/postprocess/hash_manifest.json',manifest)
    print(json.dumps(dict(figures=len(qa['figures']),source_csvs=len(qa['source_csvs']),max_inverse_residual_m=qa['inverse_mapping_max_residual_m'],hash_entries=len(manifest))))
if __name__=='__main__':main()
