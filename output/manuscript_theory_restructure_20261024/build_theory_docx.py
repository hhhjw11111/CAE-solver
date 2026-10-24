"""Build the restructured theory section as a DOCX with native OMML equations.

Layout follows docs/manuscript_word_format.md (body section of the frozen journal
template): A4, margins 1247/1106 twip, line grid 312, 10.5 pt SimSun/Times New Roman.
Equations are converted from LaTeX to OMML with pandoc and placed in borderless
three-column rows (blank | equation | number) that are not split across pages.

Usage: python build_theory_docx.py <template.docx> <output.docx>
"""
import os
import re
import shutil
import struct
import subprocess
import sys
import tempfile
import zipfile
from xml.sax.saxutils import escape

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, 'article_src.txt')

TEXT_WIDTH = 11906 - 2 * 1106  # twip
NUM_COL = 850
EQ_COL = TEXT_WIDTH - 2 * NUM_COL

FONT_BODY = '<w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman" w:eastAsia="宋体" w:cs="Times New Roman"/>'
FONT_HEI = '<w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman" w:eastAsia="黑体" w:cs="Times New Roman"/>'
MATH_RPR = ('<w:rPr><w:rFonts w:ascii="Cambria Math" w:hAnsi="Cambria Math" w:eastAsia="宋体"/>'
            '<w:sz w:val="21"/><w:szCs w:val="21"/></w:rPr>')
GRID_SPACING = '<w:spacing w:line="34" w:lineRule="atLeast"/>'


# ---------------------------------------------------------------- parsing
def parse(path):
    blocks, eq = [], None
    for raw in open(path, encoding='utf8'):
        line = raw.rstrip('\n')
        if eq is not None:
            if line.strip() == '/EQ':
                blocks.append(eq)
                eq = None
            elif line.strip():
                eq['lines'].append(line.strip())
            continue
        if blocks and blocks[-1]['type'] == 'table' and not blocks[-1].get('closed'):
            if line.strip() == '/TABLE':
                blocks[-1]['closed'] = True
            elif line.strip():
                blocks[-1]['rows'].append(line.split('|'))
            continue
        if not line.strip() or line.startswith('#'):
            continue
        key, _, rest = line.partition(' ')
        if key in ('H1', 'H2', 'H3', 'P', 'N'):
            blocks.append({'type': key, 'text': rest})
        elif key == 'EQ':
            eq = {'type': 'eq', 'label': rest.strip(), 'lines': []}
        elif key == 'FIG':
            f, zh, en = rest.split('|')
            blocks.append({'type': 'fig', 'file': f, 'zh': zh, 'en': en})
        elif key == 'TABLE':
            zh, en = rest.split('|')
            blocks.append({'type': 'table', 'zh': zh, 'en': en, 'rows': []})
        elif key == 'REFS':
            blocks.append({'type': 'refs'})
        elif key == 'REF':
            k, entry = rest.split('|', 1)
            blocks.append({'type': 'ref', 'key': k, 'text': entry})
        else:
            raise ValueError('unknown line: ' + line)
    return blocks


INLINE = re.compile(r'(\$[^$]+\$|\[\[[^\]]+\]\]|\(@[A-Za-z0-9_]+\))')


# ---------------------------------------------------------------- math via pandoc
def convert_math(snippets):
    """snippets: list of (latex, display). Returns list of OMML strings."""
    md = []
    for tex, display in snippets:
        md.append(('$$%s$$' if display else '$%s$') % tex)
    tmp = tempfile.mkdtemp()
    src = os.path.join(tmp, 'm.md')
    out = os.path.join(tmp, 'm.docx')
    open(src, 'w', encoding='utf8').write('\n\n'.join(md) + '\n')
    subprocess.run(['pandoc', src, '-f', 'markdown', '-o', out], check=True)
    xml = zipfile.ZipFile(out).read('word/document.xml').decode('utf8')
    body = xml[xml.index('<w:body>'):]
    paras = re.findall(r'<w:p>.*?</w:p>', body, flags=re.S)
    if len(paras) != len(snippets):
        raise RuntimeError('pandoc paragraph count %d != %d' % (len(paras), len(snippets)))
    result = []
    for (tex, display), p in zip(snippets, paras):
        m = re.search(r'<m:oMath>.*</m:oMath>', p, flags=re.S)
        if not m or '<w:t' in p.replace('<m:t', ''):
            raise RuntimeError('math conversion failed: ' + tex)
        result.append(fix_omml(m.group(0)))
    shutil.rmtree(tmp)
    return result


def fix_omml(x):
    x = x.replace(' />', '/>')
    # Schema order (ECMA-376): m:rPr is nor | (scr?, sty?); m:dPr is begChr, sepChr, endChr, grow.
    x = re.sub(r'<m:nor/><m:sty m:val="[a-z]+"/>', '<m:nor/>', x)
    x = re.sub(r'(<m:sty m:val="[a-z]+"/>)(<m:scr m:val="[a-z-]+"/>)', r'\2\1', x)
    x = re.sub(r'(<m:endChr m:val="[^"]*"/>)(<m:sepChr m:val="[^"]*"/>)', r'\2\1', x)
    # Vectors, tensors and matrices: bold italic (template rule).
    x = x.replace('<m:sty m:val="b"/>', '<m:sty m:val="bi"/>')
    # Word's standard bar accent is the combining overline.
    x = x.replace('<m:chr m:val="‾"/>', '<m:chr m:val="̅"/>')
    # Explicit 10.5 pt Cambria Math runs; CJK text inside equations in SimSun.
    x = re.sub(r'(<m:r>(?:<m:rPr>.*?</m:rPr>)?)', lambda m: m.group(1) + MATH_RPR, x)
    # Merge adjacent runs with identical properties (pandoc emits one run per letter).
    merge = re.compile(r'<m:r>((?:<m:rPr>(?:(?!</m:rPr>).)*</m:rPr>)?<w:rPr>(?:(?!</w:rPr>).)*</w:rPr>)'
                       r'<m:t>([^<]*)</m:t></m:r><m:r>\1<m:t>([^<]*)</m:t></m:r>')
    while True:
        y = merge.sub(r'<m:r>\1<m:t>\2\3</m:t></m:r>', x)
        if y == x:
            break
        x = y
    x = re.sub(r'<m:t>([^<]*\s[^<]*)</m:t>', r'<m:t xml:space="preserve">\1</m:t>', x)
    return x


# ---------------------------------------------------------------- XML helpers
def run(text, rpr_extra='', font=FONT_BODY, size=21):
    if not text:
        return ''
    return ('<w:r><w:rPr>%s%s<w:sz w:val="%d"/><w:szCs w:val="%d"/></w:rPr>'
            '<w:t xml:space="preserve">%s</w:t></w:r>' % (font, rpr_extra, size, size, escape(text)))


def para(content, ppr):
    return '<w:p><w:pPr>%s</w:pPr>%s</w:p>' % (ppr, content)


class Builder:
    def __init__(self, blocks):
        self.blocks = blocks
        self.eqnum, self.cites = {}, {}
        n = 0
        for b in blocks:
            if b['type'] == 'eq':
                n += 1
                self.eqnum[b['label']] = n
        # citation order = first appearance in text
        for b in blocks:
            texts = [b.get('text', '')] if b['type'] in ('H1', 'H2', 'H3', 'P', 'N') else []
            if b['type'] == 'table':
                texts = [c for r in b['rows'] for c in r]
            for t in texts:
                for keys in re.findall(r'\[\[([^\]]+)\]\]', t):
                    for k in keys.split(','):
                        self.cites.setdefault(k.strip(), len(self.cites) + 1)
        refs = [b['key'] for b in blocks if b['type'] == 'ref']
        missing = set(self.cites) ^ set(refs)
        if missing:
            raise RuntimeError('citation/reference mismatch: %s' % missing)
        # collect math
        self.snippets = []
        for b in blocks:
            for t in self._texts(b):
                for tok in INLINE.findall(t):
                    if tok.startswith('$'):
                        self.snippets.append((tok[1:-1], False))
            if b['type'] == 'eq':
                for line in b['lines']:
                    self.snippets.append((line, True))
        self.omml = convert_math(self.snippets)
        self.mi = 0

    @staticmethod
    def _texts(b):
        if b['type'] in ('H1', 'H2', 'H3', 'P', 'N'):
            return [b['text']]
        if b['type'] == 'table':
            return [c for r in b['rows'] for c in r]
        return []

    def next_math(self):
        x = self.omml[self.mi]
        self.mi += 1
        return x

    def cite_label(self, keys):
        nums = sorted(self.cites[k.strip()] for k in keys.split(','))
        parts, i = [], 0
        while i < len(nums):
            j = i
            while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1:
                j += 1
            if j - i >= 2:
                parts.append('%d-%d' % (nums[i], nums[j]))
            else:
                parts.extend(str(v) for v in nums[i:j + 1])
            i = j + 1
        return '[' + ','.join(parts) + ']'

    def inline(self, text, font=FONT_BODY, size=21, extra=''):
        out = []
        for piece in INLINE.split(text):
            if not piece:
                continue
            if piece.startswith('$') and piece.endswith('$') and len(piece) > 1:
                out.append(self.next_math())
            elif piece.startswith('[['):
                out.append(run(self.cite_label(piece[2:-2]), extra + '<w:vertAlign w:val="superscript"/>', font, size))
            elif piece.startswith('(@'):
                out.append(run('(%d)' % self.eqnum[piece[2:-1]], extra, font, size))
            else:
                out.append(run(piece, extra, font, size))
        return ''.join(out)

    # blocks
    def heading(self, b):
        level = b['type']
        text = b['text']
        if level == 'H1':
            ppr = ('<w:keepNext/>' '<w:spacing w:beforeLines="50" w:before="156" w:afterLines="50" w:after="156" '
                   'w:line="34" w:lineRule="atLeast"/><w:outlineLvl w:val="0"/>')
            return para(run(text, '<w:b/>', FONT_HEI, 24), ppr)
        if level == 'H2':
            ppr = ('<w:keepNext/><w:spacing w:beforeLines="50" w:before="156" w:line="34" w:lineRule="atLeast"/>'
                   '<w:outlineLvl w:val="1"/>')
            return para(run(text, '<w:b/>', FONT_HEI, 21), ppr)
        num, _, rest = text.partition('　')
        ppr = '<w:keepNext/>' + GRID_SPACING + '<w:outlineLvl w:val="2"/>'
        return para(run(num + '　', '<w:b/>') + run(rest), ppr)

    def body(self, b):
        ind = '<w:ind w:firstLineChars="200" w:firstLine="420"/>' if b['type'] == 'P' else ''
        return para(self.inline(b['text']), GRID_SPACING + ind + '<w:jc w:val="both"/>')

    def equation(self, b):
        no_border = ''.join('<w:%s w:val="nil"/>' % s for s in
                            ('top', 'left', 'bottom', 'right', 'insideH', 'insideV'))
        cell_ppr = '<w:keepNext/><w:keepLines/>' + GRID_SPACING
        lines = []
        for _ in b['lines']:
            omml = self.next_math()
            lines.append('<w:p><w:pPr>%s<w:jc w:val="center"/></w:pPr><m:oMathPara><m:oMathParaPr>'
                         '<m:jc m:val="center"/></m:oMathParaPr>%s</m:oMathPara></w:p>' % (cell_ppr, omml))

        def cell(width, content, valign='center'):
            return ('<w:tc><w:tcPr><w:tcW w:w="%d" w:type="dxa"/><w:vAlign w:val="%s"/></w:tcPr>%s</w:tc>'
                    % (width, valign, content))
        empty = '<w:p><w:pPr>%s</w:pPr></w:p>' % GRID_SPACING
        number = para(run('(%d)' % self.eqnum[b['label']]), GRID_SPACING + '<w:jc w:val="right"/>')
        return ('<w:tbl><w:tblPr><w:tblW w:w="%d" w:type="dxa"/><w:jc w:val="center"/>'
                '<w:tblBorders>%s</w:tblBorders><w:tblLayout w:type="fixed"/>'
                '<w:tblCellMar><w:left w:w="0" w:type="dxa"/><w:right w:w="0" w:type="dxa"/></w:tblCellMar>'
                '<w:tblLook w:val="0000"/></w:tblPr><w:tblGrid><w:gridCol w:w="%d"/><w:gridCol w:w="%d"/>'
                '<w:gridCol w:w="%d"/></w:tblGrid><w:tr><w:trPr><w:cantSplit/></w:trPr>%s%s%s</w:tr></w:tbl>'
                % (TEXT_WIDTH, no_border, NUM_COL, EQ_COL, NUM_COL,
                   cell(NUM_COL, empty), cell(EQ_COL, ''.join(lines)), cell(NUM_COL, number)))

    def figure(self, b, rid, pid):
        path = os.path.join(HERE, b['file'])
        w, h = struct.unpack('>II', open(path, 'rb').read(32)[16:24])
        cx = int(3.6 * 914400)
        cy = int(cx * h / w)
        drawing = (
            '<w:r><w:drawing><wp:inline distT="0" distB="0" distL="0" distR="0"><wp:extent cx="%d" cy="%d"/>'
            '<wp:docPr id="%d" name="Figure %d"/><wp:cNvGraphicFramePr><a:graphicFrameLocks noChangeAspect="1"/>'
            '</wp:cNvGraphicFramePr><a:graphic><a:graphicData uri="http://schemas.openxmlformats.org/drawingml/2006/picture">'
            '<pic:pic><pic:nvPicPr><pic:cNvPr id="%d" name="%s"/><pic:cNvPicPr/></pic:nvPicPr><pic:blipFill>'
            '<a:blip r:embed="%s"/><a:stretch><a:fillRect/></a:stretch></pic:blipFill><pic:spPr><a:xfrm>'
            '<a:off x="0" y="0"/><a:ext cx="%d" cy="%d"/></a:xfrm><a:prstGeom prst="rect"><a:avLst/></a:prstGeom>'
            '</pic:spPr></pic:pic></a:graphicData></a:graphic></wp:inline></w:drawing></w:r>'
            % (cx, cy, pid, pid, pid, b['file'], rid, cx, cy))
        cap = GRID_SPACING + '<w:jc w:val="center"/>'
        return (para(drawing, '<w:keepNext/><w:spacing w:before="120"/><w:jc w:val="center"/>')
                + para(run(b['zh'], '', FONT_BODY, 18), '<w:keepNext/>' + cap)
                + para(run(b['en'], '', FONT_BODY, 18), '<w:spacing w:after="120" w:line="34" w:lineRule="atLeast"/>'
                       '<w:jc w:val="center"/>'))

    def table(self, b):
        cap = '<w:keepNext/>' + GRID_SPACING + '<w:jc w:val="center"/>'
        out = [para(run(b['zh'], '<w:b/>', FONT_HEI, 18), '<w:keepNext/><w:spacing w:before="120" w:line="34" '
                    'w:lineRule="atLeast"/><w:jc w:val="center"/>'),
               para(run(b['en'], '', FONT_BODY, 18), cap)]
        widths = [680, 1800, 1700, 2250, 1640, 1624]
        assert sum(widths) == TEXT_WIDTH
        rows = []
        for r, cells in enumerate(b['rows']):
            last = r == len(b['rows']) - 1
            tcs = []
            for c, text in enumerate(cells):
                borders = ''
                if r == 0:
                    borders = ('<w:top w:val="single" w:sz="12" w:space="0" w:color="000000"/>'
                               '<w:bottom w:val="single" w:sz="6" w:space="0" w:color="000000"/>')
                elif last:
                    borders = '<w:bottom w:val="single" w:sz="12" w:space="0" w:color="000000"/>'
                content = self.inline(text, FONT_HEI if r == 0 else FONT_BODY, 18)
                tcs.append('<w:tc><w:tcPr><w:tcW w:w="%d" w:type="dxa"/><w:tcBorders>%s</w:tcBorders>'
                           '<w:vAlign w:val="center"/></w:tcPr>%s</w:tc>'
                           % (widths[c], borders, para(content, '<w:keepNext/><w:spacing w:before="40" '
                                                       'w:after="40"/><w:jc w:val="center"/>')))
            hdr = '<w:tblHeader/>' if r == 0 else ''
            rows.append('<w:tr><w:trPr><w:cantSplit/>%s</w:trPr>%s</w:tr>' % (hdr, ''.join(tcs)))
        nil = ''.join('<w:%s w:val="nil"/>' % s for s in ('top', 'left', 'bottom', 'right', 'insideH', 'insideV'))
        out.append('<w:tbl><w:tblPr><w:tblW w:w="%d" w:type="dxa"/><w:jc w:val="center"/><w:tblBorders>%s'
                   '</w:tblBorders><w:tblLayout w:type="fixed"/><w:tblLook w:val="0000"/></w:tblPr><w:tblGrid>%s'
                   '</w:tblGrid>%s</w:tbl>' % (TEXT_WIDTH, nil, ''.join('<w:gridCol w:w="%d"/>' % w for w in widths),
                                               ''.join(rows)))
        out.append(para('', GRID_SPACING))
        return ''.join(out)

    def build(self):
        out, rels, pid = [], [], 1
        refs = sorted((b for b in self.blocks if b['type'] == 'ref'), key=lambda b: self.cites[b['key']])
        for b in self.blocks:
            t = b['type']
            if t in ('H1', 'H2', 'H3'):
                out.append(self.heading(b))
            elif t in ('P', 'N'):
                out.append(self.body(b))
            elif t == 'eq':
                out.append(self.equation(b))
            elif t == 'fig':
                rid = 'rIdImg%d' % pid
                rels.append((rid, b['file']))
                out.append(self.figure(b, rid, pid))
                pid += 1
            elif t == 'table':
                out.append(self.table(b))
            elif t == 'refs':
                out.append(para(run('参考文献', '<w:b/>', FONT_HEI, 21),
                                '<w:keepNext/><w:spacing w:before="312" w:line="34" w:lineRule="atLeast"/>'
                                '<w:jc w:val="center"/>'))
                for r in refs:
                    out.append(para(run('[%d] ' % self.cites[r['key']], '', FONT_BODY, 18)
                                    + run(r['text'], '', FONT_BODY, 18),
                                    GRID_SPACING + '<w:ind w:left="420" w:hanging="420"/><w:jc w:val="both"/>'))
        assert self.mi == len(self.omml)
        return ''.join(out), rels


# ---------------------------------------------------------------- package
NS = ('xmlns:wpc="http://schemas.microsoft.com/office/word/2010/wordprocessingCanvas" '
      'xmlns:mc="http://schemas.openxmlformats.org/markup-compatibility/2006" '
      'xmlns:o="urn:schemas-microsoft-com:office:office" '
      'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships" '
      'xmlns:m="http://schemas.openxmlformats.org/officeDocument/2006/math" '
      'xmlns:v="urn:schemas-microsoft-com:vml" '
      'xmlns:wp="http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing" '
      'xmlns:w10="urn:schemas-microsoft-com:office:word" '
      'xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main" '
      'xmlns:wne="http://schemas.microsoft.com/office/word/2006/wordml" '
      'xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" '
      'xmlns:pic="http://schemas.openxmlformats.org/drawingml/2006/picture"')

SECT = ('<w:sectPr><w:footerReference w:type="default" r:id="rIdFooter1"/><w:pgSz w:w="11906" w:h="16838"/>'
        '<w:pgMar w:top="1247" w:right="1106" w:bottom="1247" w:left="1106" w:header="851" w:footer="992" '
        'w:gutter="0"/><w:cols w:space="425"/><w:docGrid w:type="lines" w:linePitch="312"/></w:sectPr>')

FOOTER = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?><w:ftr %s><w:p><w:pPr><w:jc w:val="center"/>'
          '</w:pPr><w:r><w:rPr><w:sz w:val="18"/></w:rPr><w:fldChar w:fldCharType="begin"/></w:r><w:r><w:rPr>'
          '<w:sz w:val="18"/></w:rPr><w:instrText xml:space="preserve"> PAGE </w:instrText></w:r><w:r><w:rPr>'
          '<w:sz w:val="18"/></w:rPr><w:fldChar w:fldCharType="separate"/></w:r><w:r><w:rPr><w:sz w:val="18"/>'
          '</w:rPr><w:t>1</w:t></w:r><w:r><w:rPr><w:sz w:val="18"/></w:rPr><w:fldChar w:fldCharType="end"/></w:r>'
          '</w:p></w:ftr>' % NS)

SETTINGS = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?><w:settings %s><w:zoom w:percent="100"/>'
            '<w:bordersDoNotSurroundHeader/><w:bordersDoNotSurroundFooter/><w:defaultTabStop w:val="420"/>'
            '<w:drawingGridVerticalSpacing w:val="156"/><w:displayHorizontalDrawingGridEvery w:val="0"/>'
            '<w:displayVerticalDrawingGridEvery w:val="2"/><w:characterSpacingControl w:val="compressPunctuation"/>'
            '<w:compat><w:spaceForUL/><w:balanceSingleByteDoubleByteWidth/><w:doNotLeaveBackslashAlone/>'
            '<w:ulTrailSpace/><w:doNotExpandShiftReturn/><w:adjustLineHeightInTable/><w:useFELayout/>'
            '<w:compatSetting w:name="compatibilityMode" w:uri="http://schemas.microsoft.com/office/word" w:val="15"/>'
            '<w:compatSetting w:name="overrideTableStyleFontSizeAndJustification" '
            'w:uri="http://schemas.microsoft.com/office/word" w:val="1"/>'
            '<w:compatSetting w:name="enableOpenTypeFeatures" w:uri="http://schemas.microsoft.com/office/word" w:val="1"/>'
            '<w:compatSetting w:name="doNotFlipMirrorIndents" w:uri="http://schemas.microsoft.com/office/word" w:val="1"/>'
            '</w:compat><m:mathPr><m:mathFont m:val="Cambria Math"/><m:brkBin m:val="before"/>'
            '<m:brkBinSub m:val="--"/><m:smallFrac m:val="0"/><m:dispDef/><m:lMargin m:val="0"/>'
            '<m:rMargin m:val="0"/><m:defJc m:val="centerGroup"/><m:wrapIndent m:val="1440"/>'
            '<m:intLim m:val="subSup"/><m:naryLim m:val="undOvr"/></m:mathPr>'
            '<w:themeFontLang w:val="en-US" w:eastAsia="zh-CN"/><w:decimalSymbol w:val="."/>'
            '<w:listSeparator w:val=","/></w:settings>' % NS)


def write_docx(template, output, body, rels):
    t = zipfile.ZipFile(template)
    doc = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?><w:document %s><w:body>%s%s</w:body>'
           '</w:document>' % (NS, body, SECT))
    rel_xml = ['<Relationship Id="rIdStyles" Type="http://schemas.openxmlformats.org/officeDocument/2006/'
               'relationships/styles" Target="styles.xml"/>',
               '<Relationship Id="rIdSettings" Type="http://schemas.openxmlformats.org/officeDocument/2006/'
               'relationships/settings" Target="settings.xml"/>',
               '<Relationship Id="rIdFonts" Type="http://schemas.openxmlformats.org/officeDocument/2006/'
               'relationships/fontTable" Target="fontTable.xml"/>',
               '<Relationship Id="rIdTheme" Type="http://schemas.openxmlformats.org/officeDocument/2006/'
               'relationships/theme" Target="theme/theme1.xml"/>',
               '<Relationship Id="rIdFooter1" Type="http://schemas.openxmlformats.org/officeDocument/2006/'
               'relationships/footer" Target="footer1.xml"/>']
    for rid, f in rels:
        rel_xml.append('<Relationship Id="%s" Type="http://schemas.openxmlformats.org/officeDocument/2006/'
                       'relationships/image" Target="media/%s"/>' % (rid, f))
    ct = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Types xmlns="http://schemas.openxmlformats.org/'
          'package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.'
          'relationships+xml"/><Default Extension="xml" ContentType="application/xml"/><Default Extension="png" '
          'ContentType="image/png"/><Override PartName="/word/document.xml" ContentType="application/vnd.'
          'openxmlformats-officedocument.wordprocessingml.document.main+xml"/><Override PartName="/word/styles.xml" '
          'ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.styles+xml"/><Override '
          'PartName="/word/settings.xml" ContentType="application/vnd.openxmlformats-officedocument.'
          'wordprocessingml.settings+xml"/><Override PartName="/word/fontTable.xml" ContentType="application/vnd.'
          'openxmlformats-officedocument.wordprocessingml.fontTable+xml"/><Override PartName="/word/theme/theme1.xml" '
          'ContentType="application/vnd.openxmlformats-officedocument.theme+xml"/><Override PartName="/word/'
          'footer1.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.footer+xml"/>'
          '<Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-'
          'properties+xml"/></Types>')
    root_rels = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Relationships xmlns="http://schemas.'
                 'openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.'
                 'openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/>'
                 '<Relationship Id="rId2" Type="http://schemas.openxmlformats.org/package/2006/relationships/'
                 'metadata/core-properties" Target="docProps/core.xml"/></Relationships>')
    core = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?><cp:coreProperties xmlns:cp="http://schemas.'
            'openxmlformats.org/package/2006/metadata/core-properties" xmlns:dc="http://purl.org/dc/elements/1.1/" '
            'xmlns:dcterms="http://purl.org/dc/terms/" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">'
            '<dc:title>热—压电半导体多场耦合的三维有限元方法</dc:title></cp:coreProperties>')
    with zipfile.ZipFile(output, 'w', zipfile.ZIP_DEFLATED) as z:
        z.writestr('[Content_Types].xml', ct)
        z.writestr('_rels/.rels', root_rels)
        z.writestr('docProps/core.xml', core)
        z.writestr('word/document.xml', doc)
        z.writestr('word/_rels/document.xml.rels',
                   '<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Relationships xmlns="http://schemas.'
                   'openxmlformats.org/package/2006/relationships">%s</Relationships>' % ''.join(rel_xml))
        z.writestr('word/styles.xml', t.read('word/styles.xml'))
        z.writestr('word/fontTable.xml', t.read('word/fontTable.xml'))
        z.writestr('word/theme/theme1.xml', t.read('word/theme/theme1.xml'))
        z.writestr('word/settings.xml', SETTINGS)
        z.writestr('word/footer1.xml', FOOTER)
        for _, f in rels:
            z.write(os.path.join(HERE, f), 'word/media/' + f)


def main():
    template, output = sys.argv[1], sys.argv[2]
    b = Builder(parse(SRC))
    body, rels = b.build()
    write_docx(template, output, body, rels)
    print('equations:', len(b.eqnum), 'math objects:', len(b.omml), 'references:', len(b.cites))


if __name__ == '__main__':
    main()
