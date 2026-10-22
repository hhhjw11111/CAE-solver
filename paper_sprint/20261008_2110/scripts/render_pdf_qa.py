"""Rasterize local final PDFs with bundled Poppler, preserving all pages for QA."""
from pathlib import Path
import json
import subprocess
import sys
import hashlib
from PIL import Image, ImageOps, ImageDraw
from pypdf import PdfReader

root = Path(__file__).resolve().parents[1]
pdf = Path(sys.argv[1]).resolve()
pdf_sha256 = hashlib.sha256(pdf.read_bytes()).hexdigest()
out = root / 'audit' / 'render_qa' / pdf.stem / pdf_sha256[:12]
out.mkdir(parents=True, exist_ok=True)
poppler = Path(r'C:\Users\Administrator\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\poppler\Library\bin')
cmd = [str(poppler/'pdftoppm.exe'), '-r', '105', '-png', str(pdf), str(out/'page')]
cp = subprocess.run(cmd, capture_output=True, text=True)
(out/'rasterize.log').write_text(cp.stdout+'\n'+cp.stderr, encoding='utf-8')
if cp.returncode:
    raise SystemExit(cp.returncode)
pages = sorted(out.glob('page-*.png'))
info = {'pdf':str(pdf), 'pdf_sha256':pdf_sha256, 'output_dir':str(out), 'pages':len(PdfReader(pdf).pages), 'raster_pages':len(pages), 'command':cmd, 'exit_code':cp.returncode}
for start in range(0,len(pages),4):
    sheet = Image.new('RGB',(1280,1860),'#c8c8c8')
    draw = ImageDraw.Draw(sheet)
    for k,path in enumerate(pages[start:start+4]):
        img=Image.open(path).convert('RGB')
        img.thumbnail((620,895))
        x=(k%2)*640+10; y=(k//2)*930+24
        sheet.paste(img,(x,y))
        draw.text((x,y-18),f'{pdf.stem} / page {start+k+1}',fill='black')
    sheet.save(out/f'contact-{start//4+1:02d}.png')
(out/'raster_manifest.json').write_text(json.dumps(info,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(info,ensure_ascii=False))
