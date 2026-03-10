import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as mpatches
import logomaker
from collections import Counter
import os

EXCEL_PATH = r'c:\Users\User\Downloads\UBR3 Nt screen (1).xlsx'
OUTPUT_DIR = r'c:\Users\User\Desktop\תינוקת\ubr3enrichmentlogo\figures'

AMINO_ACIDS = list('ACDEFGHIKLMNPQRSTVWY')

mpl.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'font.size': 11,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.spines.top': False,
    'axes.spines.right': False,
})

df_hits = pd.read_excel(EXCEL_PATH, sheet_name='sub_high')

sequences_2pos = []
for _, row in df_hits.iterrows():
    seq = str(row['AA2']) + str(row['AA3'])
    if len(seq) == 2 and all(c in AMINO_ACIDS for c in seq):
        sequences_2pos.append(seq)

counts_matrix = logomaker.alignment_to_matrix(sequences_2pos, to_type='counts')
info_matrix = logomaker.transform_matrix(counts_matrix, from_type='counts', to_type='information')

AA_CATEGORIES = {
    'D': 'Acidic', 'E': 'Acidic',
    'R': 'Basic', 'K': 'Basic', 'H': 'Basic',
    'G': 'Nonpolar', 'A': 'Nonpolar', 'V': 'Nonpolar', 'L': 'Nonpolar',
    'I': 'Nonpolar', 'P': 'Nonpolar', 'F': 'Nonpolar', 'M': 'Nonpolar', 'W': 'Nonpolar',
    'S': 'Polar', 'T': 'Polar', 'C': 'Polar',
    'Y': 'Polar', 'N': 'Polar', 'Q': 'Polar',
}

# ---- Define color schemes ----

schemes = {

    # 1. Canonical (standard WebLogo-like, but warm-shifted)
    'canonical_warm': {
        'colors': {
            'Acidic':   '#CC0000',
            'Basic':    '#E8751A',
            'Nonpolar': '#333333',
            'Polar':    '#B8860B',
        },
        'title': 'Canonical (warm)',
    },

    # 2. Classic biochemistry (red/blue but muted)
    'classic': {
        'colors': {
            'Acidic':   '#D32F2F',
            'Basic':    '#1976D2',
            'Nonpolar': '#424242',
            'Polar':    '#388E3C',
        },
        'title': 'Classic Biochemistry',
    },

    # 3. All-red spectrum (dark to light)
    'red_spectrum': {
        'colors': {
            'Acidic':   '#6B0F1A',
            'Basic':    '#B71C1C',
            'Nonpolar': '#E53935',
            'Polar':    '#EF9A9A',
        },
        'title': 'Red Spectrum',
    },

    # 4. Sunset (deep red → coral → peach → gold)
    'sunset': {
        'colors': {
            'Acidic':   '#8B0000',
            'Basic':    '#CD5C5C',
            'Nonpolar': '#E8751A',
            'Polar':    '#DAA520',
        },
        'title': 'Sunset',
    },

    # 5. Rose-pink palette
    'rose': {
        'colors': {
            'Acidic':   '#880E4F',
            'Basic':    '#C2185B',
            'Nonpolar': '#E91E63',
            'Polar':    '#F48FB1',
        },
        'title': 'Rose',
    },

    # 6. Terracotta (earthy warm)
    'terracotta': {
        'colors': {
            'Acidic':   '#5D4037',
            'Basic':    '#A1887F',
            'Nonpolar': '#BF360C',
            'Polar':    '#E64A19',
        },
        'title': 'Terracotta',
    },

    # 7. Coral-blush
    'coral_blush': {
        'colors': {
            'Acidic':   '#B71C1C',
            'Basic':    '#FF5252',
            'Nonpolar': '#FF8A80',
            'Polar':    '#FFAB91',
        },
        'title': 'Coral & Blush',
    },

    # 8. Charge-focused (warm): negative=deep red, positive=amber, hydrophobic=brown, polar=peach
    'charge_warm': {
        'colors': {
            'Acidic':   '#960018',
            'Basic':    '#FF8F00',
            'Nonpolar': '#6D4C41',
            'Polar':    '#FFAB91',
        },
        'title': 'Charge-focused (warm)',
    },
}


def make_logo(info_mat, scheme_name, scheme, suffix=''):
    cat_colors = scheme['colors']
    aa_colors = {aa: cat_colors[cat] for aa, cat in AA_CATEGORIES.items()}

    fig, ax = plt.subplots(figsize=(2.2, 2.2))
    logomaker.Logo(info_mat, ax=ax, color_scheme=aa_colors, font_name='DejaVu Sans')
    ax.set_ylabel('Bits', fontsize=7)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Pos 2', 'Pos 3'], fontsize=7)
    ax.tick_params(axis='y', labelsize=6)
    ax.set_title(f'Sequence Logo — Best Hits', fontweight='bold', fontsize=8, pad=3)

    legend_patches = [mpatches.Patch(color=cat_colors[c], label=c) for c in
                      ['Acidic', 'Basic', 'Nonpolar', 'Polar']]
    ax.legend(handles=legend_patches, fontsize=5, frameon=False,
              loc='upper left', bbox_to_anchor=(1.02, 1.0),
              handlelength=0.8, handletextpad=0.3, borderpad=0.2)
    plt.tight_layout()

    fname = f'fig4_logo_{scheme_name}{suffix}'
    fig.savefig(os.path.join(OUTPUT_DIR, f'{fname}.png'))
    fig.savefig(os.path.join(OUTPUT_DIR, f'{fname}.pdf'))
    plt.close(fig)
    print(f'  Saved: {fname}')


for name, scheme in schemes.items():
    make_logo(info_matrix, name, scheme)

print('\nDone — all logo variants saved.')
