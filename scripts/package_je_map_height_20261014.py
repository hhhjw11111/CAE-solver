"""Audit and package the completed high-voltage plotting exports."""
import csv
import json
import os
import subprocess
import sys
import zipfile
from pathlib import Path
import fitz
import plot_je_map_height_20261014 as plot

out = plot.OUT
env = dict(os.environ, PYTHONIOENCODING='utf-8')
summary = {}
for spec in plot.SPECS:
    stem = spec['id']
    pdf = out/'figures'/f'{stem}.pdf'
    assert (out/'figures'/f'{stem}.svg').is_file()
    reports = {}
    for tool in ['audit_pdf_text', 'audit_figure_collisions']:
        result = subprocess.run([sys.executable, str(plot.QA/f'{tool}.py'), str(pdf), '--json'],
                                capture_output=True, text=True, encoding='utf-8', env=env)
        assert result.returncode == 0, (stem, tool, result.stderr, result.stdout)
        data = json.loads(result.stdout)
        (out/'qa'/f'{stem}_{tool}.json').write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding='utf-8')
        reports[tool] = data
    collision = reports['audit_figure_collisions']
    assert collision['summary']['fail'] == 0
    assert all(f['kind']=='text-image-edge' for f in collision['findings'])
    alignment=json.loads((out/'qa'/f'{stem}_alignment.json').read_text(encoding='utf-8'))
    summary[stem] = {
        'minimum_pdf_font_pt': reports['audit_pdf_text']['minimum_found_pt'],
        'collision_summary': collision['summary'],
        'alignment': alignment.get('summary'),
        'visual_review': 'All labels and holes visible. Raster bounding-box edge warnings reviewed against rendered PNG; no actual text occlusion.',
    }
    print(stem, summary[stem]['minimum_pdf_font_pt'], collision['summary'], flush=True)

result=subprocess.run([sys.executable,str(plot.QA/'validate_figure.py'),str(Path(plot.__file__)),'--json'],
                      capture_output=True,text=True,encoding='utf-8',env=env)
assert result.returncode==0, result.stdout+result.stderr
(out/'qa'/'source_preflight.json').write_text(result.stdout,encoding='utf-8')
(out/'qa'/'audit_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
rows=[]
for name in ['plane_extrema_01_joule_mid.csv','plane_extrema_05_total_stress_mid.csv']:
    with (out/name).open(encoding='utf-8-sig',newline='') as handle:
        rows.extend(csv.DictReader(handle))
assert len(rows)==48
with (out/'plane_extrema.csv').open('w',encoding='utf-8-sig',newline='') as handle:
    writer=csv.DictWriter(handle,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
plot.write_index()
combined=fitz.open()
for spec in plot.SPECS:
    with fitz.open(out/'figures'/f'{spec["id"]}.pdf') as source:
        assert len(source)==1
        combined.insert_pdf(source)
combined.set_metadata({'title':'6–8 V verified current-work fields', 'subject':'Single-mesh J dot E diagnostic branch; not thermodynamic closure'})
combined.save(out/'高电压场图合集.pdf',garbage=4,deflate=True)
combined.close()
files=list((out/'figures').glob('*'))+[out/name for name in ['说明.md','图像索引.html','verified_results.csv','plane_extrema.csv','figure_contract.json','provenance.json','高电压场图合集.pdf','plot_je_map_height_20261014.py','plot_je_8v_clouds_20261013.py']]
with zipfile.ZipFile(out/'高电压场图交付.zip','w',compression=zipfile.ZIP_DEFLATED) as archive:
    for path in files:
        archive.write(path,path.relative_to(out))
print('PACKAGED', len(files), 'files', flush=True)
