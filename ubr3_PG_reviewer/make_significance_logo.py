"""
Significance logo for positions 4-24: letter height is -log10(p) instead of
fold change, so a residue that is MORE SIGNIFICANT is TALLER regardless of how
large its fold change happens to be. This is the pLogo / kpLogo convention.

Why this exists: fold change and significance are different quantities. W at
position 12 has the largest fold change in the data (9.65x) but only reaches
p < 0.05, because a 9.65-fold change is exactly what you get when a residue
has equal raw counts in both groups -- possibly 1 peptide vs 1 peptide. R at
positions 7, 8 and 17 has a smaller fold change but a much smaller p, because
it rests on more observations. A fold-change logo makes W look like the main
result; a significance logo makes R the main result.

LIMITATION: AMINO ACIDS.csv stores significance only as * / ** / *** bins, not
as p-values, so heights here are the bin representatives 0.05 / 0.01 / 0.001
(-log10 = 1.30 / 2.00 / 3.00). Within a bin, residues are indistinguishable.
For exact heights the per-cell p-values are needed -- see the note printed at
the end about sheet 26 of UBR3_PG_substrate_tables.xlsx.

Source: data/AMINO_ACIDS_PG_motif.csv
Substrates (n = 20) vs. non-substrate controls with the same motif (n = 193).
"""
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import logomaker

HERE       = os.path.dirname(os.path.abspath(__file__))
SRC        = os.path.join(HERE, 'data', 'AMINO_ACIDS_PG_motif.csv')
OUTPUT_DIR = os.path.join(HERE, 'figures')
os.makedirs(OUTPUT_DIR, exist_ok=True)

AA  = ["K", "H", "R", "D", "E", "F", "Y", "W", "G", "P",
       "A", "M", "V", "L", "I", "S", "T", "C", "N", "Q"]
POS = list(range(4, 25))
N_SUB, N_CTRL = 20, 193

# star bin -> representative p-value -> letter height
P_OF_STARS = {0: np.nan, 1: 0.05, 2: 0.01, 3: 0.001}

CATEGORY_MEMBERS = {
    'Basic':       list('KRH'),
    'Acidic':      list('DE'),
    'Aromatic':    list('FWY'),
    'Aliphatic':   list('GPAM'),
    'Hydrophobic': list('LIV'),
    'Polar':       list('STCNQ'),
}
CAT_COLORS = {
    'Acidic':      '#7B241C',
    'Aromatic':    '#A0522D',
    'Basic':       '#C0392B',
    'Hydrophobic': '#D35400',
    'Aliphatic':   '#E8A33D',
    'Polar':       '#F1C40F',
}
LEGEND_ORDER = ['Basic', 'Acidic', 'Aromatic', 'Aliphatic', 'Hydrophobic', 'Polar']
AA_COLOR_SCHEME = {a: CAT_COLORS[c] for c, mem in CATEGORY_MEMBERS.items() for a in mem}

# ------------------------------------------------------------------ load
raw = pd.read_csv(SRC, header=None, dtype=str, keep_default_na=False)


def block(row0, conv):
    b = raw.iloc[row0:row0 + 20, 1:22].map(conv)
    b.index, b.columns = raw.iloc[row0:row0 + 20, 0].tolist(), POS
    return b


def to_num(x):
    x = str(x).strip()
    if x in ("N/A", "", "-", "#NUM!", "nan"):
        return np.nan
    try:
        return float(x)
    except ValueError:
        return np.nan


V = block(1, to_num)
S = block(23, lambda x: str(x).strip().count("*"))
assert list(V.index) == AA and list(S.index) == AA

# height = -log10(p), and only for residues that are ENRICHED (log2FC > 0)
H = S.map(lambda n: P_OF_STARS[int(n)]).astype(float)
H = -np.log10(H)
H = H.where(V > 0).fillna(0.0)

# ------------------------------------------------------------------ style
mpl.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Helvetica', 'Arial', 'Liberation Sans', 'DejaVu Sans'],
    'font.size': 8, 'axes.linewidth': 0.8,
    'svg.fonttype': 'none', 'pdf.fonttype': 42,
    'figure.dpi': 600, 'savefig.dpi': 600, 'savefig.bbox': 'tight',
    'axes.spines.top': False, 'axes.spines.right': False,
})

fig, ax = plt.subplots(figsize=(10.5, 2.9))
logomaker.Logo(H.T.set_axis(POS), ax=ax, color_scheme=AA_COLOR_SCHEME,
               font_name='Helvetica', vpad=0.0, width=0.95,
               stack_order='big_on_top')
for patch in ax.patches:
    patch.set_edgecolor('black')
    patch.set_linewidth(0.35)

for y, lab in [(1.30103, 'p = 0.05'), (2.0, 'p = 0.01'), (3.0, 'p = 0.001')]:
    ax.axhline(y=y, color='#999', linewidth=0.5, linestyle=(0, (3, 2)), zorder=1)
    ax.text(POS[-1] + 0.55, y, lab, fontsize=5, color='#666', va='center', ha='left')

ax.axhline(y=0, color='black', linewidth=0.5)
ax.yaxis.grid(False)
ax.xaxis.grid(False)
ax.set_axisbelow(True)
ax.set_ylabel('$-$log$_{10}$ $p$   (stacked)', fontsize=8)
ax.set_xlabel('Position', fontsize=8, labelpad=2)
ax.set_xticks(POS)
ax.set_xticklabels(POS, fontsize=6.5)
ax.tick_params(axis='both', labelsize=6.5, length=2.5, width=0.7, pad=2)
ax.set_title('Significance Logo — height is $-$log$_{10}$ $p$, not fold change',
             fontweight='bold', fontsize=8.5, pad=4)
ax.set_ylim(-0.05, float(H.sum(axis=0).max()) * 1.12)
ax.set_xlim(POS[0] - 0.6, POS[-1] + 0.6)
ax.text(0.5, -0.15,
        f'n = {N_SUB} P/G-D/E/T substrates  vs  n = {N_CTRL} non-substrate controls   '
        f'(heights are * / ** / *** bin representatives, uncorrected)',
        transform=ax.transAxes, ha='center', fontsize=5.8, color='#666', style='italic')

ax.legend(handles=[mpatches.Patch(color=CAT_COLORS[c],
                                  label=f'{c}  ({" ".join(CATEGORY_MEMBERS[c])})')
                   for c in LEGEND_ORDER],
          fontsize=5.5, frameon=False, loc='upper left', bbox_to_anchor=(1.10, 1.0),
          handlelength=0.8, handletextpad=0.3, borderpad=0.2)
plt.tight_layout(pad=0.3)

stem = 'Figure15_significance_logo_pos4_24'
for ext in ('png', 'pdf', 'svg'):
    p = os.path.join(OUTPUT_DIR, f'{stem}.{ext}')
    fig.savefig(p, format=ext)
    print(f"  saved: {os.path.basename(p)}  ({os.path.getsize(p)/1024:.1f} KB)")
plt.close(fig)

# ---------------------------------------------- fold change vs significance
print("\nfold change and significance rank residues differently:")
rows = [{'Position': p, 'Residue': a,
         'Fold change': round(float(2 ** V.loc[a, p]), 2),
         'p (bin)': P_OF_STARS[int(S.loc[a, p])],
         'Stars': '*' * int(S.loc[a, p])}
        for p in POS for a in AA if S.loc[a, p] >= 1 and V.loc[a, p] > 0]
t = pd.DataFrame(rows)
print("\n  top 5 by FOLD CHANGE:")
print(t.nlargest(5, 'Fold change').to_string(index=False))
print("\n  top 5 by SIGNIFICANCE:")
print(t.nsmallest(5, 'p (bin)').to_string(index=False))
