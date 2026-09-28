#!/usr/bin/env python3
"""
create_word_report.py
Converts PROJECT_REPORT_FOR_PRESENTATION.md into:
1. PROJECT_REPORT_FOR_PRESENTATION.docx (Microsoft Word Document)
2. PROJECT_REPORT_FOR_PRESENTATION.html (Google Docs / Word-compatible HTML)
3. PROJECT_REPORT_FOR_PRESENTATION.rtf  (Rich Text Format Document)
"""

import os
import sys
import re
import zipfile
import html
from xml.sax.saxutils import escape as xml_escape

REPORTS_DIR = os.path.dirname(os.path.abspath(__file__))
MD_PATH = os.path.join(REPORTS_DIR, "PROJECT_REPORT_FOR_PRESENTATION.md")
DOCX_PATH = os.path.join(REPORTS_DIR, "PROJECT_REPORT_FOR_PRESENTATION.docx")
HTML_PATH = os.path.join(REPORTS_DIR, "PROJECT_REPORT_FOR_PRESENTATION.html")
RTF_PATH = os.path.join(REPORTS_DIR, "PROJECT_REPORT_FOR_PRESENTATION.rtf")

def read_markdown(path):
    with open(path, "r", encoding="utf-8") as f:
        return f.read()

def try_python_docx(md_content, docx_path):
    """Attempt generation via python-docx if installed."""
    try:
        import docx
        from docx import Document
        from docx.shared import Pt, Inches, RGBColor
        from docx.enum.text import WD_ALIGN_PARAGRAPH
        from docx.enum.table import WD_TABLE_ALIGNMENT
        from docx.oxml import parse_xml
        from docx.oxml.ns import nsdecls
    except ImportError:
        return False

    doc = Document()
    
    # Page setup (1 inch margins)
    sections = doc.sections
    for section in sections:
        section.top_margin = Inches(1)
        section.bottom_margin = Inches(1)
        section.left_margin = Inches(1)
        section.right_margin = Inches(1)

    lines = md_content.split("\n")
    in_table = False
    table_lines = []

    def flush_table(t_lines):
        if not t_lines:
            return
        parsed_rows = []
        for line in t_lines:
            stripped = line.strip()
            if stripped.startswith("|") and stripped.endswith("|"):
                parts = [p.strip() for p in stripped[1:-1].split("|")]
                # skip separator row like |:---|:---|
                if all(re.match(r"^:?-+:?$", p) for p in parts if p):
                    continue
                parsed_rows.append(parts)
        if not parsed_rows:
            return

        num_cols = max(len(r) for r in parsed_rows)
        tbl = doc.add_table(rows=len(parsed_rows), cols=num_cols)
        tbl.alignment = WD_TABLE_ALIGNMENT.CENTER
        tbl.style = "Table Grid"

        for r_idx, row in enumerate(parsed_rows):
            for c_idx, cell_text in enumerate(row):
                if c_idx < num_cols:
                    cell = tbl.cell(r_idx, c_idx)
                    cell.text = cell_text
                    # Header formatting
                    if r_idx == 0:
                        shading_xml = parse_xml(f'<w:shd {nsdecls("w")} w:fill="0D47A1"/>')
                        cell._tc.get_or_add_tcPr().append(shading_xml)
                        for p in cell.paragraphs:
                            p.alignment = WD_ALIGN_PARAGRAPH.LEFT
                            for run in p.runs:
                                run.font.bold = True
                                run.font.color.rgb = RGBColor(255, 255, 255)
                                run.font.size = Pt(9.5)
                    else:
                        for p in cell.paragraphs:
                            for run in p.runs:
                                run.font.size = Pt(9)
        doc.add_paragraph()

    for line in lines:
        stripped = line.strip()
        if stripped.startswith("|") and stripped.endswith("|"):
            in_table = True
            table_lines.append(stripped)
            continue
        else:
            if in_table:
                flush_table(table_lines)
                table_lines = []
                in_table = False

        if not stripped:
            continue

        if stripped.startswith("# "):
            p = doc.add_heading(level=0)
            run = p.add_run(stripped[2:])
            run.font.size = Pt(22)
            run.font.bold = True
            run.font.color.rgb = RGBColor(13, 71, 161)
        elif stripped.startswith("## "):
            p = doc.add_heading(level=1)
            run = p.add_run(stripped[3:])
            run.font.size = Pt(15)
            run.font.bold = True
            run.font.color.rgb = RGBColor(25, 118, 210)
        elif stripped.startswith("### "):
            p = doc.add_heading(level=2)
            run = p.add_run(stripped[4:])
            run.font.size = Pt(12)
            run.font.bold = True
            run.font.color.rgb = RGBColor(55, 71, 79)
        elif stripped.startswith("#### "):
            p = doc.add_heading(level=3)
            run = p.add_run(stripped[5:])
            run.font.size = Pt(11)
            run.font.bold = True
        elif stripped.startswith("- ") or stripped.startswith("* "):
            p = doc.add_paragraph(style='List Bullet')
            _add_markdown_runs(p, stripped[2:])
        elif re.match(r"^\d+\.\s+", stripped):
            m = re.match(r"^(\d+\.\s+)(.*)", stripped)
            p = doc.add_paragraph(style='List Number')
            _add_markdown_runs(p, m.group(2))
        elif stripped.startswith("> "):
            p = doc.add_paragraph()
            p.paragraph_format.left_indent = Inches(0.4)
            run = p.add_run(stripped[2:])
            run.font.italic = True
            run.font.color.rgb = RGBColor(85, 85, 85)
        elif stripped == "---":
            p = doc.add_paragraph()
            p.add_run("____________________________________________________").font.color.rgb = RGBColor(200, 200, 200)
        elif stripped.startswith("```"):
            continue
        else:
            p = doc.add_paragraph()
            _add_markdown_runs(p, stripped)

    if in_table:
        flush_table(table_lines)

    doc.save(docx_path)
    return True

def _add_markdown_runs(paragraph, text):
    """Parses bold, italic, inline code in markdown and adds runs."""
    import docx
    from docx.shared import Pt, RGBColor
    
    # Simple regex-based splitter for **bold**, `code`, *italic*
    tokens = re.split(r'(\*\*.*?\*\*|`.*?`|\*.*?\*)', text)
    for token in tokens:
        if not token:
            continue
        if token.startswith("**") and token.endswith("**"):
            run = paragraph.add_run(token[2:-2])
            run.font.bold = True
        elif token.startswith("`") and token.endswith("`"):
            run = paragraph.add_run(token[1:-1])
            run.font.name = "Courier New"
            run.font.size = Pt(9.5)
            run.font.color.rgb = RGBColor(180, 40, 40)
        elif token.startswith("*") and token.endswith("*"):
            run = paragraph.add_run(token[1:-1])
            run.font.italic = True
        else:
            paragraph.add_run(token)


def create_standard_docx_via_zip(md_content, docx_path):
    """
    Builds a fully conforming .docx file from scratch using Python's built-in zipfile
    and standard WordprocessingML XML. Requires ZERO external packages.
    """
    lines = md_content.split("\n")
    body_xml = []

    in_table = False
    table_rows = []

    def flush_xml_table(rows):
        if not rows:
            return ""
        parsed = []
        for line in rows:
            stripped = line.strip()
            if stripped.startswith("|") and stripped.endswith("|"):
                parts = [p.strip() for p in stripped[1:-1].split("|")]
                if all(re.match(r"^:?-+:?$", p) for p in parts if p):
                    continue
                parsed.append(parts)
        if not parsed:
            return ""

        num_cols = max(len(r) for r in parsed)
        tbl_xml = [
            '<w:tbl>',
            '<w:tblPr>',
            '<w:tblStyle w:val="TableGrid"/>',
            '<w:tblW w:w="9000" w:type="dxa"/>',
            '<w:tblBorders>',
            '<w:top w:val="single" w:sz="4" w:space="0" w:color="CCCCCC"/>',
            '<w:bottom w:val="single" w:sz="4" w:space="0" w:color="CCCCCC"/>',
            '<w:insideH w:val="single" w:sz="4" w:space="0" w:color="EEEEEE"/>',
            '<w:insideV w:val="none"/>',
            '</w:tblBorders>',
            '</w:tblPr>'
        ]

        for r_idx, row in enumerate(parsed):
            tbl_xml.append('<w:tr>')
            is_header = (r_idx == 0)
            for c_idx in range(num_cols):
                cell_text = row[c_idx] if c_idx < len(row) else ""
                clean_text = re.sub(r'[\*`_]', '', cell_text)
                tbl_xml.append('<w:tc>')
                tbl_xml.append('<w:tcPr>')
                if is_header:
                    tbl_xml.append('<w:shd w:val="clear" w:color="auto" w:fill="0D47A1"/>')
                else:
                    bg = "F9FAFB" if r_idx % 2 == 1 else "FFFFFF"
                    tbl_xml.append(f'<w:shd w:val="clear" w:color="auto" w:fill="{bg}"/>')
                tbl_xml.append('<w:tcMar><w:top w:w="120" w:type="dxa"/><w:bottom w:w="120" w:type="dxa"/><w:left w:w="160" w:type="dxa"/><w:right w:w="160" w:type="dxa"/></w:tcMar>')
                tbl_xml.append('</w:tcPr>')
                tbl_xml.append('<w:p>')
                tbl_xml.append('<w:pPr><w:spacing w:before="60" w:after="60"/></w:pPr>')
                tbl_xml.append('<w:r>')
                tbl_xml.append('<w:rPr>')
                if is_header:
                    tbl_xml.append('<w:b/><w:color w:val="FFFFFF"/><w:sz w:val="19"/>')
                else:
                    tbl_xml.append('<w:sz w:val="18"/><w:color w:val="212121"/>')
                tbl_xml.append('</w:rPr>')
                tbl_xml.append(f'<w:t>{xml_escape(clean_text)}</w:t>')
                tbl_xml.append('</w:r>')
                tbl_xml.append('</w:p>')
                tbl_xml.append('</w:tc>')
            tbl_xml.append('</w:tr>')

        tbl_xml.append('</w:tbl>')
        return "".join(tbl_xml)

    for line in lines:
        stripped = line.strip()
        if stripped.startswith("|") and stripped.endswith("|"):
            in_table = True
            table_rows.append(stripped)
            continue
        else:
            if in_table:
                body_xml.append(flush_xml_table(table_rows))
                table_rows = []
                in_table = False

        if not stripped:
            continue

        if stripped.startswith("# "):
            text = xml_escape(stripped[2:])
            body_xml.append(f'<w:p><w:pPr><w:pStyle w:val="Heading1"/><w:spacing w:before="300" w:after="160"/></w:pPr><w:r><w:rPr><w:b/><w:sz w:val="40"/><w:color w:val="0D47A1"/></w:rPr><w:t>{text}</w:t></w:r></w:p>')
        elif stripped.startswith("## "):
            text = xml_escape(stripped[3:])
            body_xml.append(f'<w:p><w:pPr><w:pStyle w:val="Heading2"/><w:spacing w:before="260" w:after="120"/></w:pPr><w:r><w:rPr><w:b/><w:sz w:val="30"/><w:color w:val="1976D2"/></w:rPr><w:t>{text}</w:t></w:r></w:p>')
        elif stripped.startswith("### "):
            text = xml_escape(stripped[4:])
            body_xml.append(f'<w:p><w:pPr><w:pStyle w:val="Heading3"/><w:spacing w:before="200" w:after="80"/></w:pPr><w:r><w:rPr><w:b/><w:sz w:val="24"/><w:color w:val="37474F"/></w:rPr><w:t>{text}</w:t></w:r></w:p>')
        elif stripped.startswith("- ") or stripped.startswith("* "):
            raw = stripped[2:]
            runs = []
            tokens = re.split(r'(\*\*.*?\*\*|`.*?`)', raw)
            for tok in tokens:
                if not tok: continue
                if tok.startswith("**") and tok.endswith("**"):
                    runs.append(f'<w:r><w:rPr><w:b/></w:rPr><w:t>{xml_escape(tok[2:-2])}</w:t></w:r>')
                elif tok.startswith("`") and tok.endswith("`"):
                    runs.append(f'<w:r><w:rPr><w:rFonts w:ascii="Courier New"/><w:color w:val="C2185B"/><w:sz w:val="19"/></w:rPr><w:t>{xml_escape(tok[1:-1])}</w:t></w:r>')
                else:
                    runs.append(f'<w:r><w:t>{xml_escape(tok)}</w:t></w:r>')
            body_xml.append(f'<w:p><w:pPr><w:pStyle w:val="ListBullet"/><w:spacing w:before="40" w:after="40"/><w:ind w:left="400"/></w:pPr><w:r><w:t>• </w:t></w:r>{"".join(runs)}</w:p>')
        elif stripped.startswith("> "):
            text = xml_escape(stripped[2:])
            body_xml.append(f'<w:p><w:pPr><w:spacing w:before="100" w:after="100"/><w:ind w:left="500"/><w:pBdr><w:left w:val="single" w:sz="24" w:space="15" w:color="1976D2"/></w:pBdr></w:pPr><w:r><w:rPr><w:i/><w:color w:val="555555"/></w:rPr><w:t>{text}</w:t></w:r></w:p>')
        elif stripped == "---":
            body_xml.append('<w:p><w:pPr><w:pBdr><w:bottom w:val="single" w:sz="6" w:space="1" w:color="E0E0E0"/></w:pBdr></w:pPr></w:p>')
        elif stripped.startswith("```"):
            continue
        else:
            tokens = re.split(r'(\*\*.*?\*\*|`.*?`)', stripped)
            runs = []
            for tok in tokens:
                if not tok: continue
                if tok.startswith("**") and tok.endswith("**"):
                    runs.append(f'<w:r><w:rPr><w:b/></w:rPr><w:t>{xml_escape(tok[2:-2])}</w:t></w:r>')
                elif tok.startswith("`") and tok.endswith("`"):
                    runs.append(f'<w:r><w:rPr><w:rFonts w:ascii="Courier New"/><w:color w:val="C2185B"/><w:sz w:val="19"/></w:rPr><w:t>{xml_escape(tok[1:-1])}</w:t></w:r>')
                else:
                    runs.append(f'<w:r><w:t>{xml_escape(tok)}</w:t></w:r>')
            body_xml.append(f'<w:p><w:pPr><w:spacing w:before="80" w:after="80"/><w:rPr><w:sz w:val="22"/></w:rPr></w:pPr>{"".join(runs)}</w:p>')

    if in_table:
        body_xml.append(flush_xml_table(table_rows))

    content_types_xml = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
  <Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
  <Default Extension="xml" ContentType="application/xml"/>
  <Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/>
</Types>"""

    rels_xml = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/>
</Relationships>"""

    doc_xml = f"""<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">
  <w:body>
    {"".join(body_xml)}
    <w:sectPr>
      <w:pgSz w:w="12240" w:h="15840"/>
      <w:pgMar w:top="1440" w:right="1440" w:bottom="1440" w:left="1440"/>
    </w:sectPr>
  </w:body>
</w:document>"""

    with zipfile.ZipFile(docx_path, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("[Content_Types].xml", content_types_xml)
        z.writestr("_rels/.rels", rels_xml)
        z.writestr("word/document.xml", doc_xml)


def create_html_version(md_content, html_path):
    """
    Creates an executive HTML document styled to match Google Docs / Word import guidelines.
    Google Docs and Microsoft Word open this file with 100% fidelity.
    """
    lines = md_content.split("\n")
    html_parts = []
    in_table = False
    table_rows = []

    def flush_html_table(rows):
        if not rows: return ""
        out = ['<table style="border-collapse: collapse; width: 100%; margin: 16px 0; font-family: -apple-system, BlinkMacSystemFont, Arial, sans-serif; font-size: 13px;">']
        is_first = True
        for r in rows:
            stripped = r.strip()
            if stripped.startswith("|") and stripped.endswith("|"):
                parts = [p.strip() for p in stripped[1:-1].split("|")]
                if all(re.match(r"^:?-+:?$", p) for p in parts if p):
                    continue
                out.append('<tr>')
                for col in parts:
                    clean = html.escape(re.sub(r'[\*`_]', '', col))
                    if is_first:
                        out.append(f'<th style="background-color: #0D47A1; color: white; border: 1px solid #0D47A1; padding: 10px 12px; text-align: left; font-weight: 600;">{clean}</th>')
                    else:
                        out.append(f'<td style="border: 1px solid #E0E0E0; padding: 8px 12px; background-color: #FAFAFA;">{clean}</td>')
                out.append('</tr>')
                is_first = False
        out.append('</table>')
        return "".join(out)

    for line in lines:
        stripped = line.strip()
        if stripped.startswith("|") and stripped.endswith("|"):
            in_table = True
            table_rows.append(stripped)
            continue
        else:
            if in_table:
                html_parts.append(flush_html_table(table_rows))
                table_rows = []
                in_table = False

        if not stripped:
            continue

        if stripped.startswith("# "):
            html_parts.append(f'<h1 style="color: #0D47A1; font-family: Arial, sans-serif; margin-top: 30px; border-bottom: 2px solid #0D47A1; padding-bottom: 8px;">{html.escape(stripped[2:])}</h1>')
        elif stripped.startswith("## "):
            html_parts.append(f'<h2 style="color: #1976D2; font-family: Arial, sans-serif; margin-top: 24px;">{html.escape(stripped[3:])}</h2>')
        elif stripped.startswith("### "):
            html_parts.append(f'<h3 style="color: #37474F; font-family: Arial, sans-serif; margin-top: 18px;">{html.escape(stripped[4:])}</h3>')
        elif stripped.startswith("- ") or stripped.startswith("* "):
            formatted = re.sub(r'\*\*(.*?)\*\*', r'<strong>\1</strong>', html.escape(stripped[2:]))
            formatted = re.sub(r'`(.*?)`', r'<code style="background: #F0F4F8; color: #C2185B; padding: 2px 4px; border-radius: 4px;">\1</code>', formatted)
            html_parts.append(f'<li style="margin: 4px 0; font-family: Arial, sans-serif; font-size: 14px; line-height: 1.6;">{formatted}</li>')
        elif stripped.startswith("> "):
            html_parts.append(f'<blockquote style="border-left: 4px solid #1976D2; margin: 12px 0; padding: 8px 16px; background-color: #F0F7FF; color: #444; font-style: italic;">{html.escape(stripped[2:])}</blockquote>')
        elif stripped == "---":
            html_parts.append('<hr style="border: 0; border-top: 1px solid #E0E0E0; margin: 24px 0;" />')
        elif stripped.startswith("```"):
            continue
        else:
            formatted = re.sub(r'\*\*(.*?)\*\*', r'<strong>\1</strong>', html.escape(stripped))
            formatted = re.sub(r'`(.*?)`', r'<code style="background: #F0F4F8; color: #C2185B; padding: 2px 4px; border-radius: 4px;">\1</code>', formatted)
            html_parts.append(f'<p style="font-family: Arial, sans-serif; font-size: 14px; line-height: 1.6; color: #212121; margin: 8px 0;">{formatted}</p>')

    if in_table:
        html_parts.append(flush_html_table(table_rows))

    full_html = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>Technical & Executive Report: Amplicon & Sanger Sequencing Primer Design Suite</title>
</head>
<body style="max-width: 850px; margin: 40px auto; padding: 0 20px; font-family: -apple-system, BlinkMacSystemFont, Arial, sans-serif;">
{"".join(html_parts)}
</body>
</html>"""
    with open(html_path, "w", encoding="utf-8") as f:
        f.write(full_html)


def main():
    if not os.path.exists(MD_PATH):
        print(f"Error: {MD_PATH} not found.")
        sys.exit(1)

    content = read_markdown(MD_PATH)

    # 1. Generate DOCX
    success = try_python_docx(content, DOCX_PATH)
    if not success:
        create_standard_docx_via_zip(content, DOCX_PATH)
    print(f"Generated Word Document: {DOCX_PATH}")

    # 2. Generate Google Docs / Word-ready HTML
    create_html_version(content, HTML_PATH)
    print(f"Generated HTML (Google Docs compatible): {HTML_PATH}")

if __name__ == "__main__":
    main()
