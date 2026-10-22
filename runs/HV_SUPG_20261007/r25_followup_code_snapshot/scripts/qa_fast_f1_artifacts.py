"""Render report pages and image contact sheets; check export completeness."""
from pathlib import Path
import argparse
import json
import fitz
import numpy as np
from PIL import Image, ImageDraw


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--delivery-root", type=Path, required=True)
    parser.add_argument("--qa-root", type=Path, required=True)
    args = parser.parse_args()
    root, qa = args.delivery_root.resolve(), args.qa_root.resolve()
    qa.mkdir(parents=True, exist_ok=True)
    pdf = fitz.open(root/"F1_六算例汇报.pdf")
    overview = Image.new("RGB", (960, 3*475), "white")
    for i, page in enumerate(pdf):
        pix = page.get_pixmap(matrix=fitz.Matrix(1.5, 1.5), alpha=False)
        path = qa/f"report_page_{i+1:02}.png"
        pix.save(path)
        with Image.open(path) as im:
            thumb = im.copy()
            thumb.thumbnail((310, 440))
            x, y = (i % 3)*320, (i//3)*475
            overview.paste(thumb, (x, y+20))
            ImageDraw.Draw(overview).text((x+8, y+3), f"Page {i+1}", fill="black")
    overview.save(qa/"report_overview.png")
    images = sorted(p for folder in ("01_折线对比图", "02_单幅场图", "03_组合对比图")
                    for p in (root/folder).glob("*.png"))
    checks = {}
    for i in range(0, len(images), 16):
        sheet = Image.new("RGB", (1280, 4*290), "white")
        draw = ImageDraw.Draw(sheet)
        for j, path in enumerate(images[i:i+16]):
            with Image.open(path) as im:
                pixels = np.asarray(im.convert("RGB"))
                checks[path.name] = {
                    "nonblank": bool(pixels.std() > 10),
                    "size": list(im.size),
                    "svg": path.with_suffix(".svg").exists(),
                    "pdf": path.with_suffix(".pdf").exists(),
                }
                thumb = im.convert("RGB")
                thumb.thumbnail((310, 255))
            x, y = (j % 4)*320, (j//4)*290
            sheet.paste(thumb, (x+(320-thumb.width)//2, y+20))
            draw.text((x+5, y+3), path.stem[:47], fill="black")
        sheet.save(qa/f"gallery_contact_{i//16+1:02}.png")
    result = {"image_count": len(images), "pdf_pages": len(pdf), "images": checks,
              "export_checks_passed": all(c["nonblank"] and c["svg"] and c["pdf"] for c in checks.values()),
              "visual_review": "PENDING_HUMAN_STYLE_INSPECTION"}
    (qa/"artifact_checks.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps({k: v for k, v in result.items() if k != "images"}))


if __name__ == "__main__":
    main()
