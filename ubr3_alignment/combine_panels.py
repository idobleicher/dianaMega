"""
Combine the conservation plot and animal domain architecture
into a single composite figure with Panel A and Panel B labels.
"""

import os
from PIL import Image, ImageDraw, ImageFont

OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))

PANEL_A = os.path.join(OUTPUT_DIR, "ubr3_domain_conservation.png")
PANEL_B = os.path.join(OUTPUT_DIR, "ubr3_animal_domain_architecture.png")
OUTPUT  = os.path.join(OUTPUT_DIR, "ubr3_combined_conservation_architecture.png")


def main():
    img_a = Image.open(PANEL_A)
    img_b = Image.open(PANEL_B)

    target_w = max(img_a.width, img_b.width)

    def scale_to_width(img, tw):
        r = tw / img.width
        return img.resize((tw, int(img.height * r)), Image.LANCZOS)

    img_a = scale_to_width(img_a, target_w)
    img_b = scale_to_width(img_b, target_w)

    label_h = 60
    gap = 20
    total_h = label_h + img_a.height + gap + label_h + img_b.height + gap

    combined = Image.new("RGB", (target_w, total_h), "white")
    draw = ImageDraw.Draw(combined)

    try:
        font = ImageFont.truetype("arialbd.ttf", 42)
    except:
        try:
            font = ImageFont.truetype("arial.ttf", 42)
        except:
            font = ImageFont.load_default()

    y = 10
    draw.text((20, y), "A", fill="black", font=font)
    y += label_h
    combined.paste(img_a, (0, y))
    y += img_a.height + gap

    draw.text((20, y), "B", fill="black", font=font)
    y += label_h
    combined.paste(img_b, (0, y))

    combined.save(OUTPUT, dpi=(300, 300))
    print(f"Saved: {OUTPUT}")


if __name__ == "__main__":
    main()
