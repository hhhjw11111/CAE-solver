"""Render the frozen Chinese execution plan. No solver execution or data changes."""
from __future__ import annotations

import hashlib
import html
import json
import re
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    BaseDocTemplate, Frame, KeepTogether, PageBreak, PageTemplate,
    Paragraph, Spacer, Table, TableStyle,
)
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'output/pdf/算例与执行流程冻结方案_20260906.md'
DEST = SOURCE.with_suffix('.pdf')
QA = ROOT / 'tmp/pdfs/case_freeze_20260906'

pdfmetrics.registerFont(TTFont('Song', 'C:/Windows/Fonts/simsun.ttc', subfontIndex=0))
pdfmetrics.registerFont(TTFont('Hei', 'C:/Windows/Fonts/simhei.ttf'))
pdfmetrics.registerFont(TTFont('Latin', 'C:/Windows/Fonts/arial.ttf'))
pdfmetrics.registerFont(TTFont('Symbols', 'C:/Windows/Fonts/seguisym.ttf'))
pdfmetrics.registerFontFamily('Song', normal='Song', bold='Hei', italic='Song', boldItalic='Hei')

NAVY = colors.HexColor('#14364B')
TEAL = colors.HexColor('#187F88')
GREY = colors.HexColor('#52616B')
LIGHT = colors.HexColor('#EDF4F6')
LINE = colors.HexColor('#D8E2E7')
WIDTH, HEIGHT = A4
MARGIN = 47
AVAILABLE = WIDTH - 2 * MARGIN

STYLES = {
    'body': ParagraphStyle('body', fontName='Song', fontSize=10.0, leading=15.5,
                           wordWrap='CJK', textColor=NAVY, spaceAfter=6),
    'small': ParagraphStyle('small', fontName='Song', fontSize=8.4, leading=12,
                            wordWrap='CJK', textColor=GREY, spaceAfter=6),
    'title': ParagraphStyle('title', fontName='Hei', fontSize=21, leading=28,
                            textColor=NAVY, spaceBefore=3, spaceAfter=13),
    'subtitle': ParagraphStyle('subtitle', fontName='Hei', fontSize=16, leading=22,
                               textColor=TEAL, spaceAfter=15),
    'h3': ParagraphStyle('h3', fontName='Hei', fontSize=12.1, leading=18,
                         textColor=TEAL, spaceBefore=7, spaceAfter=6, keepWithNext=True),
    'cell': ParagraphStyle('cell', fontName='Song', fontSize=9.1, leading=13,
                           textColor=NAVY, wordWrap='CJK'),
    'headcell': ParagraphStyle('headcell', fontName='Hei', fontSize=9.2, leading=14,
                               textColor=colors.white, wordWrap='CJK'),
    'quote': ParagraphStyle('quote', fontName='Song', fontSize=10.2, leading=17,
                            textColor=NAVY, wordWrap='CJK'),
}


def markup(raw: str) -> str:
    sup_map = str.maketrans('⁰¹²³⁴⁵⁶⁷⁸⁹⁻', '0123456789-')
    raw = re.sub('[⁰¹²³⁴⁵⁶⁷⁸⁹⁻]+', lambda m: '<sup>' + m[0].translate(sup_map) + '</sup>', raw)
    raw = raw.replace('₀', '<sub>0</sub>')
    out = html.escape(raw, quote=False)
    for tag in ('sub', 'sup', 'br', 'br/'):
        out = out.replace('&lt;' + tag + '&gt;', '<' + tag + '>')
        out = out.replace('&lt;/' + tag + '&gt;', '</' + tag + '>')
    out = re.sub(r'\*\*(.+?)\*\*', r'<b>\1</b>', out)
    for char in set(out):
        if ord(char) > 127 and ord(char) not in pdfmetrics.getFont('Song').face.charToGlyph:
            fallback = next((name for name in ('Latin', 'Symbols')
                             if ord(char) in pdfmetrics.getFont(name).face.charToGlyph), None)
            if fallback is None:
                raise ValueError('Missing glyph: ' + repr(char))
            out = out.replace(char, '<font name="' + fallback + '">' + char + '</font>')
    return out


def p(text: str, style: str = 'body') -> Paragraph:
    return Paragraph(markup(text), STYLES[style])


def table(rows: list[list[str]]) -> Table:
    ncols = len(rows[0])
    if any(len(row) != ncols for row in rows):
        raise ValueError('Inconsistent table column count: ' + repr(rows))
    if ncols == 2:
        ratios = (0.27, 0.73)
    elif ncols == 3:
        ratios = (0.21, 0.42, 0.37)
    else:
        ratios = (0.19, 0.29, 0.28, 0.24)
    cells = [[p(cell, 'headcell' if index == 0 else 'cell') for cell in row]
             for index, row in enumerate(rows)]
    result = Table(cells, colWidths=[r * AVAILABLE for r in ratios], repeatRows=1,
                   hAlign='LEFT')
    result.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), NAVY),
        ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, LIGHT]),
        ('VALIGN', (0, 0), (-1, -1), 'TOP'),
        ('LEFTPADDING', (0, 0), (-1, -1), 7),
        ('RIGHTPADDING', (0, 0), (-1, -1), 7),
        ('TOPPADDING', (0, 0), (-1, -1), 5.5),
        ('BOTTOMPADDING', (0, 0), (-1, -1), 5.5),
        ('LINEBELOW', (0, -1), (-1, -1), .6, LINE),
        ('LINEBELOW', (0, 0), (-1, 0), .8, TEAL),
    ]))
    return result


def header_footer(canvas, doc):
    canvas.saveState()
    canvas.setFillColor(TEAL)
    canvas.rect(MARGIN, HEIGHT - 36, 26, 3, stroke=0, fill=1)
    canvas.setFont('Hei', 8)
    canvas.setFillColor(GREY)
    canvas.drawString(MARGIN + 35, HEIGHT - 35, 'TPSolverV2  /  算例与执行流程  /  F1')
    canvas.setStrokeColor(LINE)
    canvas.line(MARGIN, 39, WIDTH - MARGIN, 39)
    canvas.setFont('Song', 8)
    canvas.drawString(MARGIN, 25, '2026-09-06  |  研究执行方案，不等同于实验定标结论')
    canvas.setFont('Hei', 9)
    canvas.drawRightString(WIDTH - MARGIN, 25, f'{doc.page:02d}')
    canvas.restoreState()


def parse_story(source: str):
    story = []
    lines = source.splitlines()
    index = 0
    while index < len(lines):
        text = lines[index].strip()
        index += 1
        if not text:
            continue
        if text == '<!-- page -->':
            story.append(PageBreak())
            continue
        if text.startswith('|'):
            rows = []
            while True:
                row = [part.strip() for part in text.strip('|').split('|')]
                if not all(re.fullmatch(r':?-+:?', cell) for cell in row):
                    rows.append(row)
                if index >= len(lines) or not lines[index].strip().startswith('|'):
                    break
                text = lines[index].strip()
                index += 1
            story.extend([table(rows), Spacer(1, 7)])
            continue
        if text.startswith('### '):
            story.append(p(text[4:], 'h3'))
        elif text.startswith('## '):
            story.append(p(text[3:], 'subtitle'))
        elif text.startswith('# '):
            story.append(p(text[2:], 'title'))
        elif text.startswith('> '):
            quote = Table([[p(text[2:], 'quote')]], colWidths=[AVAILABLE])
            quote.setStyle(TableStyle([
                ('BACKGROUND', (0, 0), (-1, -1), LIGHT),
                ('LINEBEFORE', (0, 0), (0, 0), 2.3, TEAL),
                ('LEFTPADDING', (0, 0), (-1, -1), 12),
                ('RIGHTPADDING', (0, 0), (-1, -1), 12),
                ('TOPPADDING', (0, 0), (-1, -1), 9),
                ('BOTTOMPADDING', (0, 0), (-1, -1), 9),
            ]))
            story.extend([quote, Spacer(1, 9)])
        else:
            if text.startswith('- '):
                text = '• ' + text[2:]
            small = text.startswith(('• solver_v2/', '既有材料审计的主要文献索引', 'Material Passport'))
            story.append(p(text, 'small' if small else 'body'))
    return story


def main():
    QA.mkdir(parents=True, exist_ok=True)
    source = SOURCE.read_text(encoding='utf-8')
    doc = BaseDocTemplate(str(DEST), pagesize=A4, leftMargin=MARGIN, rightMargin=MARGIN,
                          topMargin=58, bottomMargin=51,
                          title='孔洞多场耦合研究：算例与执行流程冻结方案',
                          author='TPSolverV2 项目', subject='执行版 F1；两类解析解与六个主算例')
    frame = Frame(MARGIN, 51, AVAILABLE, HEIGHT - 109,
                  leftPadding=0, bottomPadding=0, rightPadding=0, topPadding=0)
    doc.addPageTemplates(PageTemplate(id='report', frames=[frame], onPage=header_footer))
    doc.build(parse_story(source))
    reader = PdfReader(DEST)
    contents = [page.extract_text() or '' for page in reader.pages]
    # markup() raises for unsupported glyphs; superscripts/subscripts are real markup.
    markup(source)
    missing = []
    metrics = {
        'pdf': str(DEST), 'page_count': len(reader.pages), 'expected_pages': 9,
        'page_char_counts': list(map(len, contents)),
        'source_sha256': hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        'pdf_sha256': hashlib.sha256(DEST.read_bytes()).hexdigest(),
        'missing_song_glyphs': missing,
        'visual_review': 'PENDING',
    }
    (QA / 'text_check.json').write_text(json.dumps(metrics, ensure_ascii=False, indent=2), encoding='utf-8')
    (QA / 'extracted_text.txt').write_text('\n\n'.join(contents), encoding='utf-8')
    print(json.dumps(metrics, ensure_ascii=True, indent=2))


if __name__ == '__main__':
    main()
