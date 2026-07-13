#!/usr/bin/env python3
"""Figure: amino-acid composition of the 5 residues AFTER the ribosome skip
(downstream positions P1-P5), to show the UBR3 P-[D/E] motif is uncommon.
P1 = the skip proline (always P); P2 = the UBR3-relevant residue (D or E)."""
import csv, sys, collections
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import logomaker
import pandas as pd

CSV = sys.argv[1]; OUT = sys.argv[2]
AAs = list("ACDEFGHIKLMNPQRSTVWY")           # 20 standard, alphabetical
NPOS = 5

rows = [r for r in csv.DictReader(open(CSV, encoding='utf-8'))
        if len(r['downstream_product']) >= NPOS]
N = len(rows)

# counts[pos][aa]
counts = [collections.Counter() for _ in range(NPOS)]
for r in rows:
    for p in range(NPOS):
        aa = r['downstream_product'][p]
        if aa in AAs:
            counts[p][aa] += 1

# percentage matrix: rows=AA, cols=position
pct = np.zeros((20, NPOS))
for p in range(NPOS):
    tot = sum(counts[p].values())
    for i, aa in enumerate(AAs):
        pct[i, p] = 100.0 * counts[p][aa] / tot

de_p2 = pct[AAs.index('D'), 1] + pct[AAs.index('E'), 1]

# ---------------- figure ----------------
plt.rcParams.update({'font.size': 11, 'font.family': 'DejaVu Sans'})
fig = plt.figure(figsize=(15, 6.2))
gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.25, 1.0], wspace=0.38)
POSLAB = [f'P{i+1}' for i in range(NPOS)]

# Panel A — probability sequence logo
axA = fig.add_subplot(gs[0, 0])
prob = pd.DataFrame(0.0, index=range(NPOS), columns=AAs)
for p in range(NPOS):
    tot = sum(counts[p].values())
    for aa in AAs:
        prob.loc[p, aa] = counts[p][aa] / tot
logomaker.Logo(prob, ax=axA, color_scheme='chemistry', show_spines=False)
axA.set_xticks(range(NPOS)); axA.set_xticklabels(POSLAB)
axA.set_ylabel('probability'); axA.set_ylim(0, 1)
axA.set_title('A  Sequence logo of downstream P1–P5', loc='left', fontweight='bold')
axA.axvline(-0.5, color='#bbbbbb', lw=0)  # spacer
axA.annotate('skip →', xy=(-0.45, 1.02), fontsize=9, color='#555')

# Panel B — single-hue frequency heatmap
axB = fig.add_subplot(gs[0, 1])
im = axB.imshow(pct, aspect='auto', cmap='Blues', vmin=0, vmax=100)
axB.set_xticks(range(NPOS)); axB.set_xticklabels(POSLAB)
axB.set_yticks(range(20)); axB.set_yticklabels(AAs, fontsize=9)
axB.set_ylabel('amino acid')
axB.set_title('B  Frequency (%) by position', loc='left', fontweight='bold')
# annotate every cell lightly; bold the dominant ones
for i in range(20):
    for p in range(NPOS):
        v = pct[i, p]
        if v >= 5:
            axB.text(p, i, f'{v:.0f}', ha='center', va='center', fontsize=7,
                     color='white' if v > 55 else '#333')
# highlight the UBR3 cells: D and E at P2
for aa in ('D', 'E'):
    r0 = AAs.index(aa)
    axB.add_patch(Rectangle((0.5, r0-0.5), 1, 1, fill=False, edgecolor='#E8743B', lw=2.2))
cbar = fig.colorbar(im, ax=axB, fraction=0.046, pad=0.03)
cbar.set_label('% of sequences', fontsize=9)

# Panel C — position-2 composition, D/E highlighted
axC = fig.add_subplot(gs[0, 2])
order = np.argsort(-pct[:, 1])
labels = [AAs[i] for i in order]
vals = [pct[i, 1] for i in order]
colors = ['#E8743B' if AAs[i] in ('D', 'E') else '#b8bcc4' for i in order]
axC.bar(range(20), vals, color=colors, width=0.8)
axC.set_xticks(range(20)); axC.set_xticklabels(labels, fontsize=9)
axC.set_ylabel('% at position 2')
axC.set_title('C  Position 2 composition', loc='left', fontweight='bold')
axC.spines[['top', 'right']].set_visible(False)
axC.annotate(f'UBR3 motif  P2 = D or E\n= {de_p2:.1f}%  of all sequences',
             xy=(0.62, 0.82), xycoords='axes fraction', fontsize=10,
             color='#E8743B', fontweight='bold', ha='left')

fig.suptitle(f'Amino-acid composition of the 5 residues after the ribosome skip '
             f'(N = {N:,} resolved 2A sequences)', fontsize=13, y=1.0)
fig.text(0.5, -0.02, 'P1 is proline in 100% of sequences (the skip proline). '
         'The UBR3 P–[D/E] motif requires D or E at P2, present in only '
         f'{de_p2:.1f}% — i.e. it is uncommon.', ha='center', fontsize=9.5, color='#444')
fig.savefig(OUT, dpi=200, bbox_inches='tight', facecolor='white')
print('wrote', OUT, '| N =', N, '| P2 D/E % =', round(de_p2, 2))
print('P1:', counts[0].most_common(1))
print('P2 top:', counts[1].most_common(6))
