"""
Create a PowerPoint presentation with all UBR3 analysis figures,
one figure per slide, organized into logical sections.
"""

import os
import glob
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.enum.text import PP_ALIGN
from pptx.dml.color import RGBColor
from PIL import Image

OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))

SLIDE_GROUPS = [
    ("Overview", [
        ("ubr3_animal_domain_architecture.png",
         "UBR3 Domain Architecture Across the Animal Kingdom"),
        ("ubr3_pub_combined_panels.png",
         "Combined Publication Figure (Alignments + Domains)"),
    ]),
    ("Conservation Analysis", [
        ("ubr3_conservation_line.png",
         "Conservation Score Along UBR3 Alignment"),
        ("ubr3_conservation_heatmap.png",
         "Conservation Heatmap"),
        ("ubr3_conserved_only_heatmap.png",
         "Conserved Regions Heatmap"),
        ("ubr3_conserved_regions_map.png",
         "Conserved Regions Map"),
        ("ubr3_top_regions_annotated.png",
         "Top Conserved Regions with Amino Acid Annotations"),
    ]),
    ("Domain Analysis", [
        ("ubr3_pub_domain_architecture.png",
         "Domain Architecture Diagram"),
        ("ubr3_domain_conservation.png",
         "Domain-Specific Conservation"),
        ("ubr3_identity_matrix.png",
         "Pairwise Sequence Identity Matrix"),
    ]),
    ("Publication Alignments", [
        ("ubr3_pub_ubr_box.png",
         "UBR Box Domain Alignment (aa 118-189)"),
        ("ubr3_pub_region_400_600.png",
         "Conserved Region Alignment (aa 400-600)"),
        ("ubr3_pub_ring_domain.png",
         "RING Domain Alignment (aa 1306-1364)"),
    ]),
    ("Full Alignment Blocks", [
        (f"ubr3_annotated_block_{i:02d}.png",
         f"Annotated Alignment Block {i}")
        for i in range(1, 13)
    ]),
    ("Individual Conserved Regions", [
        ("ubr3_region_01_h120-149.png",   "Conserved Region: aa 120-149"),
        ("ubr3_region_02_h156-174.png",   "Conserved Region: aa 156-174"),
        ("ubr3_region_03_h393-409.png",   "Conserved Region: aa 393-409"),
        ("ubr3_region_04_h415-430.png",   "Conserved Region: aa 415-430"),
        ("ubr3_region_05_h435-458.png",   "Conserved Region: aa 435-458"),
        ("ubr3_region_06_h469-483.png",   "Conserved Region: aa 469-483"),
        ("ubr3_region_07_h504-527.png",   "Conserved Region: aa 504-527"),
        ("ubr3_region_08_h547-587.png",   "Conserved Region: aa 547-587"),
        ("ubr3_region_09_h626-637.png",   "Conserved Region: aa 626-637"),
        ("ubr3_region_10_h684-699.png",   "Conserved Region: aa 684-699"),
        ("ubr3_region_11_h975-974.png",   "Conserved Region: near aa 975"),
        ("ubr3_region_12_h1015-1014.png", "Conserved Region: near aa 1015"),
        ("ubr3_region_13_h1048-1047.png", "Conserved Region: near aa 1048"),
        ("ubr3_region_14_h1050-1049.png", "Conserved Region: near aa 1050"),
        ("ubr3_region_15_h1065-1064.png", "Conserved Region: near aa 1065"),
        ("ubr3_region_16_h1119-1118.png", "Conserved Region: near aa 1119"),
        ("ubr3_region_17_h1252-1264.png", "Conserved Region: aa 1252-1264"),
        ("ubr3_region_18_h1590-1589.png", "Conserved Region: near aa 1590"),
        ("ubr3_region_19_h1615-1615.png", "Conserved Region: near aa 1615"),
        ("ubr3_region_20_h1623-1622.png", "Conserved Region: near aa 1623"),
    ]),
]

SLIDE_W = Inches(13.333)
SLIDE_H = Inches(7.5)
TITLE_H = Inches(0.85)
MARGIN = Inches(0.3)


def add_section_slide(prs, section_title):
    """Add a section divider slide."""
    slide = prs.slides.add_slide(prs.slide_layouts[6])  # blank
    slide.background.fill.solid()
    slide.background.fill.fore_color.rgb = RGBColor(0x1A, 0x23, 0x7E)

    txBox = slide.shapes.add_textbox(
        Inches(1), Inches(2.5), SLIDE_W - Inches(2), Inches(2.5))
    tf = txBox.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.text = section_title
    p.font.size = Pt(44)
    p.font.bold = True
    p.font.color.rgb = RGBColor(0xFF, 0xFF, 0xFF)
    p.alignment = PP_ALIGN.CENTER


def add_figure_slide(prs, img_path, title):
    """Add a slide with a title and the figure image centered and maximized."""
    slide = prs.slides.add_slide(prs.slide_layouts[6])  # blank

    txBox = slide.shapes.add_textbox(
        MARGIN, Inches(0.1), SLIDE_W - 2 * MARGIN, TITLE_H)
    tf = txBox.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.text = title
    p.font.size = Pt(22)
    p.font.bold = True
    p.font.color.rgb = RGBColor(0x1A, 0x23, 0x7E)
    p.alignment = PP_ALIGN.CENTER

    with Image.open(img_path) as img:
        img_w, img_h = img.size

    avail_w = SLIDE_W - 2 * MARGIN
    avail_h = SLIDE_H - TITLE_H - Inches(0.3)

    scale_w = avail_w / img_w
    scale_h = avail_h / img_h
    scale = min(scale_w, scale_h)

    final_w = int(img_w * scale)
    final_h = int(img_h * scale)

    left = (SLIDE_W - final_w) // 2
    top = TITLE_H + (avail_h - final_h) // 2

    slide.shapes.add_picture(img_path, left, top, final_w, final_h)


def main():
    print("=" * 60)
    print("Creating UBR3 Analysis PowerPoint Presentation")
    print("=" * 60)

    prs = Presentation()
    prs.slide_width = SLIDE_W
    prs.slide_height = SLIDE_H

    title_slide = prs.slides.add_slide(prs.slide_layouts[6])
    title_slide.background.fill.solid()
    title_slide.background.fill.fore_color.rgb = RGBColor(0x0D, 0x47, 0xA1)

    txBox = title_slide.shapes.add_textbox(
        Inches(1), Inches(2), SLIDE_W - Inches(2), Inches(2))
    tf = txBox.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.text = "UBR3 Sequence Conservation Analysis"
    p.font.size = Pt(40)
    p.font.bold = True
    p.font.color.rgb = RGBColor(0xFF, 0xFF, 0xFF)
    p.alignment = PP_ALIGN.CENTER

    p2 = tf.add_paragraph()
    p2.text = "Cross-species alignment and domain architecture"
    p2.font.size = Pt(22)
    p2.font.color.rgb = RGBColor(0xBB, 0xDE, 0xFB)
    p2.alignment = PP_ALIGN.CENTER
    p2.space_before = Pt(16)

    p3 = tf.add_paragraph()
    p3.text = "13 species from Mammalia to Cnidaria"
    p3.font.size = Pt(18)
    p3.font.color.rgb = RGBColor(0x90, 0xCA, 0xF9)
    p3.alignment = PP_ALIGN.CENTER
    p3.space_before = Pt(8)

    slide_count = 1
    for section_title, figures in SLIDE_GROUPS:
        add_section_slide(prs, section_title)
        slide_count += 1
        print(f"\n  Section: {section_title}")

        for fname, fig_title in figures:
            img_path = os.path.join(OUTPUT_DIR, fname)
            if not os.path.exists(img_path):
                print(f"    SKIP (not found): {fname}")
                continue
            add_figure_slide(prs, img_path, fig_title)
            slide_count += 1
            print(f"    Added: {fname}")

    out_path = os.path.join(OUTPUT_DIR, "UBR3_Analysis_Figures.pptx")
    prs.save(out_path)
    print(f"\n{'='*60}")
    print(f"Saved: {out_path}")
    print(f"Total slides: {slide_count}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
