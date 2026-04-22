"""
Regenerate fig3_dipeptide_enrichment_bars, truncated to show only up through
the 'ER' dipeptide, saved as PNG + PDF + SVG (true vector).
"""
import os, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from collections import Counter

EXCEL_PATH = r'c:\Users\User\Downloads\UBR3 Nt screen (1).xlsx'
OUTPUT_DIR = r'c:\Users\User\Desktop\תינוקת\dianaMega\ubr3enrichmentlogo\figures'
os.makedirs(OUTPUT_DIR, exist_ok=True)

AMINO_ACIDS = list('ACDEFGHIKLMNPQRSTVWY')

mpl.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'font.size': 11,
    'axes.titlesize': 14,
    'axes.labelsize': 12,
    'svg.fonttype': 'none',
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.spines.top': False,
    'axes.spines.right': False,
})

df_all  = pd.read_excel(EXCEL_PATH, sheet_name='Nprot_5_analyzed')
df_hits = pd.read_excel(EXCEL_PATH, sheet_name='sub_high')

def get_dipeptide_freq(df):
    dipeptides = df['AA2'].astype(str) + df['AA3'].astype(str)
    counts = Counter(dipeptides)
    total = sum(counts.values())
    return {k: v / total for k, v in counts.items()}

lib_dipep = get_dipeptide_freq(df_all)
hit_dipep = get_dipeptide_freq(df_hits)

top_hit_dipeps = sorted(hit_dipep.keys(), key=lambda x: hit_dipep[x], reverse=True)[:20]
enrichments_dp = [
    (hit_dipep[dp] / lib_dipep[dp]) if lib_dipep.get(dp, 0) > 0 else 0
    for dp in top_hit_dipeps
]

sorted_idx     = np.argsort(enrichments_dp)[::-1]
sorted_dipeps  = [top_hit_dipeps[i]   for i in sorted_idx]
sorted_enrich  = [enrichments_dp[i]   for i in sorted_idx]

CUTOFF_DIPEPTIDE = 'EL'
if CUTOFF_DIPEPTIDE in sorted_dipeps:
    cut = sorted_dipeps.index(CUTOFF_DIPEPTIDE) + 1
    sorted_dipeps = sorted_dipeps[:cut]
    sorted_enrich = sorted_enrich[:cut]
    print(f"Truncated at position {cut} (up through '{CUTOFF_DIPEPTIDE}')")
else:
    print(f"'{CUTOFF_DIPEPTIDE}' not in top 20 — keeping full chart")
    cut = len(sorted_dipeps)

n = len(sorted_dipeps)
print(f"Bars shown: {n}")
for rank, (dp, e) in enumerate(zip(sorted_dipeps, sorted_enrich), 1):
    print(f"  {rank:2d}. {dp}  enrichment = {e:.2f}")

color_gradient = plt.cm.YlOrRd(np.linspace(0.3, 0.95, n))[::-1]

fig, ax = plt.subplots(figsize=(max(6, 0.5 * n + 2), 5))
ax.bar(range(n), sorted_enrich,
       color=[color_gradient[i] for i in range(n)],
       edgecolor='white', linewidth=0.5)
ax.axhline(y=1.0, color='#7B7D7D', linestyle='--', linewidth=0.8, alpha=0.7)
ax.set_xticks(range(n))
ax.set_xticklabels(sorted_dipeps, rotation=45, ha='right', fontsize=10)
ax.set_xlabel('Dipeptide (Position 2-3)')
ax.set_ylabel('Enrichment (hits / library)')
ax.set_title('Top Dipeptide Enrichment at Positions 2-3 (after Met)\nBest Hits vs Library',
             fontweight='bold')
plt.tight_layout()

base = 'fig3_dipeptide_enrichment_bars'
for ext in ('png', 'pdf', 'svg'):
    path = os.path.join(OUTPUT_DIR, f'{base}.{ext}')
    fig.savefig(path, format=ext)
    size_kb = os.path.getsize(path) / 1024
    print(f"  saved: {path}  ({size_kb:.1f} KB)")

plt.close(fig)
