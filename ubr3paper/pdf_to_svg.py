"""
Convert a PDF to SVG (one SVG per page).
Uses PyMuPDF's built-in get_svg_image(), which produces true vector SVG
(text stays as text, paths as paths).
"""
import os, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import pymupdf

PDF_PATH = r'C:\Users\User\Downloads\SCREEN.pdf'
OUT_DIR  = r'C:\Users\User\Downloads'

base = os.path.splitext(os.path.basename(PDF_PATH))[0]
doc  = pymupdf.open(PDF_PATH)

print(f"PDF: {PDF_PATH}")
print(f"Pages: {doc.page_count}")

outputs = []
for i, page in enumerate(doc):
    svg = page.get_svg_image(text_as_path=False)
    if doc.page_count == 1:
        out_path = os.path.join(OUT_DIR, f'{base}.svg')
    else:
        out_path = os.path.join(OUT_DIR, f'{base}_page{i+1}.svg')
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write(svg)
    size_kb = os.path.getsize(out_path) / 1024
    print(f"  page {i+1}: {out_path}  ({size_kb:.1f} KB)")
    outputs.append(out_path)

doc.close()
print(f"\nDone. {len(outputs)} SVG file(s) written.")
