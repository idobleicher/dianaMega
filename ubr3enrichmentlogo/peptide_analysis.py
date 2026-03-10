import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as mpatches
import logomaker
from collections import Counter
from scipy import stats
import os

EXCEL_PATH = r'c:\Users\User\Downloads\UBR3 Nt screen (1).xlsx'
OUTPUT_DIR = r'c:\Users\User\Desktop\תינוקת\ubr3enrichmentlogo\figures'
os.makedirs(OUTPUT_DIR, exist_ok=True)

COLOR_MAIN = '#C0392B'
COLOR_SEC = '#E67E22'
COLOR_LIGHT = '#F5B041'
COLOR_PALE = '#FADBD8'
COLOR_ACCENT = '#D35400'

AMINO_ACIDS = list('ACDEFGHIKLMNPQRSTVWY')

AA_CATEGORIES = {
    'D': 'Acidic',   'E': 'Acidic',
    'R': 'Basic',    'K': 'Basic',    'H': 'Basic',
    'G': 'Nonpolar', 'A': 'Nonpolar', 'V': 'Nonpolar', 'L': 'Nonpolar',
    'I': 'Nonpolar', 'P': 'Nonpolar', 'F': 'Nonpolar', 'M': 'Nonpolar', 'W': 'Nonpolar',
    'S': 'Polar',    'T': 'Polar',    'C': 'Polar',
    'Y': 'Polar',    'N': 'Polar',    'Q': 'Polar',
}

CAT_COLORS = {
    'Acidic':   '#922B21',
    'Basic':    '#C0392B',
    'Nonpolar': '#D35400',
    'Polar':    '#F39C12',
}

AA_COLOR_SCHEME = {aa: CAT_COLORS[cat] for aa, cat in AA_CATEGORIES.items()}

mpl.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'font.size': 11,
    'axes.titlesize': 14,
    'axes.labelsize': 12,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.spines.top': False,
    'axes.spines.right': False,
})

df_all = pd.read_excel(EXCEL_PATH, sheet_name='Nprot_5_analyzed')
df_hits = pd.read_excel(EXCEL_PATH, sheet_name='sub_high')
print(f"Library: {len(df_all)} peptides | Best hits: {len(df_hits)} peptides")

# ===================== HELPERS =====================

def get_freq(series):
    counts = Counter(series.dropna())
    total = sum(counts.values())
    return {aa: counts.get(aa, 0) / total for aa in AMINO_ACIDS}

def get_counts(series):
    counts = Counter(series.dropna())
    return {aa: counts.get(aa, 0) for aa in AMINO_ACIDS}

def enrichment(hit_freq, lib_freq):
    return {aa: (hit_freq.get(aa, 0) / lib_freq[aa]) if lib_freq.get(aa, 0) > 0 else 0
            for aa in AMINO_ACIDS}

def save(fig, name):
    fig.savefig(os.path.join(OUTPUT_DIR, f'{name}.png'))
    fig.savefig(os.path.join(OUTPUT_DIR, f'{name}.pdf'))
    plt.close(fig)
    print(f"  Saved: {name}")

lib_aa2_freq = get_freq(df_all['AA2'])
lib_aa3_freq = get_freq(df_all['AA3'])
hit_aa2_freq = get_freq(df_hits['AA2'])
hit_aa3_freq = get_freq(df_hits['AA3'])
enrich_aa2 = enrichment(hit_aa2_freq, lib_aa2_freq)
enrich_aa3 = enrichment(hit_aa3_freq, lib_aa3_freq)

# ============================================================
# FIG 1: Enrichment bar chart — sorted by enrichment
# ============================================================
print("\n--- Fig 1: Enrichment bars ---")
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
for ax, enr, pos_label in zip(axes, [enrich_aa2, enrich_aa3],
                                ['Position 2 (after Met)', 'Position 3']):
    sorted_aa = sorted(enr.keys(), key=lambda x: enr[x], reverse=True)
    vals = [enr[aa] for aa in sorted_aa]
    colors = [COLOR_MAIN if v >= 1.5 else COLOR_SEC if v >= 1.0
              else COLOR_LIGHT if v >= 0.5 else COLOR_PALE for v in vals]
    ax.bar(sorted_aa, vals, color=colors, edgecolor='white', linewidth=0.5)
    ax.axhline(y=1.0, color='#7B7D7D', linestyle='--', linewidth=0.8, alpha=0.7)
    ax.set_xlabel('Amino Acid')
    ax.set_ylabel('Enrichment (hits / library)')
    ax.set_title(f'{pos_label} — Best Hits vs Library')
    ax.set_ylim(0, max(vals) * 1.15 if vals else 2)
plt.suptitle('Amino Acid Enrichment at N-terminal Positions (after Met)\nBest Hits vs Full Library',
             fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
save(fig, 'fig1_enrichment_pos2_pos3_bars')

# ============================================================
# FIG 1b: Enrichment HEATMAP — 20x2 (AA × position)
# ============================================================
print("--- Fig 1b: Enrichment heatmap ---")
enrich_matrix = pd.DataFrame({
    'Position 2': enrich_aa2,
    'Position 3': enrich_aa3,
})
sorted_by_pos2 = enrich_matrix.sort_values('Position 2', ascending=False)

lib_counts_2 = get_counts(df_all['AA2'])
lib_counts_3 = get_counts(df_all['AA3'])
hit_counts_2 = get_counts(df_hits['AA2'])
hit_counts_3 = get_counts(df_hits['AA3'])
n_lib = len(df_all)
n_hit = len(df_hits)

pval_dict = {}
for aa in AMINO_ACIDS:
    for pos_label, hc, lc in [('Position 2', hit_counts_2, lib_counts_2),
                               ('Position 3', hit_counts_3, lib_counts_3)]:
        table = [[hc[aa], n_hit - hc[aa]],
                 [lc[aa], n_lib - lc[aa]]]
        _, p = stats.fisher_exact(table, alternative='two-sided')
        pval_dict[(aa, pos_label)] = p

fig, ax = plt.subplots(figsize=(4.5, 8))
from matplotlib.colors import LinearSegmentedColormap
cmap_warm = LinearSegmentedColormap.from_list('warm', ['#FEF9E7', '#F5B041', '#E67E22', '#C0392B', '#78281F'])
im = ax.imshow(sorted_by_pos2.values, aspect='auto', cmap=cmap_warm, vmin=0,
               vmax=max(sorted_by_pos2.values.max(), 5))
ax.set_xticks([0, 1])
ax.set_xticklabels(['Position 2', 'Position 3'], fontsize=12)
ax.set_yticks(range(len(sorted_by_pos2)))
ax.set_yticklabels(sorted_by_pos2.index, fontsize=11)

col_labels = ['Position 2', 'Position 3']
for i, aa in enumerate(sorted_by_pos2.index):
    for j in range(2):
        v = sorted_by_pos2.values[i, j]
        p = pval_dict[(aa, col_labels[j])]
        stars = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
        color = 'white' if v > 2.5 else 'black'
        ax.text(j, i, f'{v:.1f}{stars}', ha='center', va='center', fontsize=9,
                color=color, fontweight='bold' if stars else 'normal')

cbar = plt.colorbar(im, ax=ax, shrink=0.6, label='Enrichment')
ax.set_title(f'Enrichment Heatmap (Best Hits / Library)\n'
             f'(n={n_hit} hits, n={n_lib:,} library)', fontweight='bold')
fig.text(0.5, 0.005, '* p<0.05  ** p<0.01  *** p<0.001 (Fisher exact test)',
         ha='center', fontsize=7.5, color='#555')
plt.tight_layout(rect=[0, 0.025, 1, 1])
save(fig, 'fig1b_enrichment_heatmap')

# ============================================================
# FIG 1c: Dipeptide enrichment HEATMAP (pos2 × pos3)
# ============================================================
print("--- Fig 1c: Dipeptide enrichment heatmap ---")
def get_dipeptide_freq(df):
    dipeptides = df['AA2'].astype(str) + df['AA3'].astype(str)
    counts = Counter(dipeptides)
    total = sum(counts.values())
    return {k: v / total for k, v in counts.items()}

lib_dipep = get_dipeptide_freq(df_all)
hit_dipep = get_dipeptide_freq(df_hits)

heatmap_data = np.zeros((20, 20))
for i, aa2 in enumerate(AMINO_ACIDS):
    for j, aa3 in enumerate(AMINO_ACIDS):
        dp = aa2 + aa3
        lib_f = lib_dipep.get(dp, 0)
        hit_f = hit_dipep.get(dp, 0)
        heatmap_data[i, j] = (hit_f / lib_f) if lib_f > 0 else 0

fig, ax = plt.subplots(figsize=(10, 9))
masked = np.ma.masked_where(heatmap_data == 0, heatmap_data)
im = ax.imshow(masked, cmap=cmap_warm, aspect='auto', vmin=0, vmax=min(heatmap_data.max(), 50))
ax.set_xticks(range(20))
ax.set_xticklabels(AMINO_ACIDS, fontsize=10)
ax.set_yticks(range(20))
ax.set_yticklabels(AMINO_ACIDS, fontsize=10)
ax.set_xlabel('Position 3 (AA3)', fontsize=12)
ax.set_ylabel('Position 2 (AA2)', fontsize=12)
for i in range(20):
    for j in range(20):
        v = heatmap_data[i, j]
        if v > 5:
            ax.text(j, i, f'{v:.0f}', ha='center', va='center', fontsize=7,
                    color='white' if v > 20 else 'black', fontweight='bold')
plt.colorbar(im, ax=ax, shrink=0.7, label='Enrichment (hits / library)')
ax.set_title('Dipeptide Enrichment Heatmap (Position 2 × Position 3)\nBest Hits vs Library',
             fontweight='bold')
plt.tight_layout()
save(fig, 'fig1c_dipeptide_heatmap')

# ============================================================
# FIG 2: Frequency comparison — SORTED BY HIT FREQUENCY
# ============================================================
print("\n--- Fig 2: Frequency comparison (sorted by frequency) ---")
fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
for ax, lib_f, hit_f, pos_label in zip(
    axes, [lib_aa2_freq, lib_aa3_freq], [hit_aa2_freq, hit_aa3_freq],
    ['Position 2 (after Met)', 'Position 3']
):
    sorted_aa = sorted(AMINO_ACIDS, key=lambda aa: hit_f.get(aa, 0), reverse=True)
    lib_vals = [lib_f.get(aa, 0) * 100 for aa in sorted_aa]
    hit_vals = [hit_f.get(aa, 0) * 100 for aa in sorted_aa]
    x = np.arange(len(sorted_aa))
    width = 0.35
    ax.bar(x - width/2, lib_vals, width, label='Library (all)', color=COLOR_LIGHT, edgecolor='white')
    ax.bar(x + width/2, hit_vals, width, label='Best Hits', color=COLOR_MAIN, edgecolor='white')
    ax.set_xlabel('Amino Acid')
    ax.set_ylabel('Frequency (%)')
    ax.set_title(f'{pos_label} — Frequency Comparison')
    ax.set_xticks(x)
    ax.set_xticklabels(sorted_aa)
    ax.legend(frameon=False)
plt.suptitle('Amino Acid Frequency at N-terminal Positions (after Met)\nLibrary vs Best Hits (sorted by hit frequency)',
             fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
save(fig, 'fig2_frequency_comparison_bars')

# ============================================================
# FIG 3: Top dipeptide enrichment bars (unchanged)
# ============================================================
print("\n--- Fig 3: Top dipeptide enrichment bars ---")
top_hit_dipeps = sorted(hit_dipep.keys(), key=lambda x: hit_dipep[x], reverse=True)[:20]
enrichments_dp = []
for dp in top_hit_dipeps:
    enrichments_dp.append((hit_dipep[dp] / lib_dipep[dp]) if lib_dipep.get(dp, 0) > 0 else 0)

n = len(top_hit_dipeps)
color_gradient = plt.cm.YlOrRd(np.linspace(0.3, 0.95, n))[::-1]
sorted_idx = np.argsort(enrichments_dp)[::-1]
sorted_dipeps = [top_hit_dipeps[i] for i in sorted_idx]
sorted_enrich = [enrichments_dp[i] for i in sorted_idx]

fig, ax = plt.subplots(figsize=(10, 5))
ax.bar(range(n), sorted_enrich, color=[color_gradient[i] for i in range(n)],
       edgecolor='white', linewidth=0.5)
ax.axhline(y=1.0, color='#7B7D7D', linestyle='--', linewidth=0.8, alpha=0.7)
ax.set_xticks(range(n))
ax.set_xticklabels(sorted_dipeps, rotation=45, ha='right', fontsize=10)
ax.set_xlabel('Dipeptide (Position 2-3)')
ax.set_ylabel('Enrichment (hits / library)')
ax.set_title('Top Dipeptide Enrichment at Positions 2-3 (after Met)\nBest Hits vs Library',
             fontweight='bold')
plt.tight_layout()
save(fig, 'fig3_dipeptide_enrichment_bars')

# ============================================================
# FIG 4: Logo plot — compact, category-colored, legend outside
# ============================================================
print("\n--- Fig 4: Fold-change logo (best hits vs library) ---")
fc_matrix = pd.DataFrame({0: enrich_aa2, 1: enrich_aa3}).T

fc_display = fc_matrix.copy()
fc_display[fc_display <= 1] = 0

fig, ax = plt.subplots(figsize=(4.5, 3.5))
logo = logomaker.Logo(fc_display, ax=ax, color_scheme=AA_COLOR_SCHEME,
                      font_name='DejaVu Sans', vpad=0.0, width=0.95,
                      stack_order='big_on_top')
for patch in ax.patches:
    patch.set_edgecolor('black')
    patch.set_linewidth(0.5)

ax.axhline(y=0, color='black', linewidth=0.5)
ax.yaxis.grid(True, linestyle='--', alpha=0.25, color='#7B7D7D')
ax.set_axisbelow(True)

ax.set_ylabel('Fold Change (hits / library)', fontsize=9)
ax.set_xticks([0, 1])
ax.set_xticklabels(['Pos 2', 'Pos 3'], fontsize=9)
ax.tick_params(axis='y', labelsize=7)
ax.set_title('Fold-Change Logo — Best Hits vs Library', fontweight='bold', fontsize=10, pad=6)

max_stack = max(sum(v for v in enr.values() if v > 1)
                for enr in [enrich_aa2, enrich_aa3])
ax.set_ylim(-0.3, max_stack + 2.0)

for pos_idx, enr in enumerate([enrich_aa2, enrich_aa3]):
    top_aa = max(enr, key=enr.get)
    top_fc = enr[top_aa]
    stack_top = sum(v for v in enr.values() if v > 1)
    ax.text(pos_idx, stack_top + 0.4, f'{top_aa} {top_fc:.1f}x',
            ha='center', fontsize=7.5, fontweight='bold', color='#333')

ax.text(0.5, -0.10, f'n = {len(df_hits)} best hits  vs  n = {len(df_all):,} library',
        transform=ax.transAxes, ha='center', fontsize=7, color='#666', style='italic')

legend_patches = [mpatches.Patch(color=CAT_COLORS[c], label=c) for c in
                  ['Acidic', 'Basic', 'Nonpolar', 'Polar']]
ax.legend(handles=legend_patches, fontsize=6, frameon=False,
          loc='upper left', bbox_to_anchor=(1.02, 1.0),
          handlelength=0.8, handletextpad=0.3, borderpad=0.2)
plt.tight_layout()
save(fig, 'fig4_logo_pos2_pos3_besthits')

# ============================================================
# FIG 4b: Logo plot WITH Met at position 1
# ============================================================
print("--- Fig 4b: Logo with Met at pos 1 ---")
sequences_3pos = []
for _, row in df_hits.iterrows():
    seq = 'M' + str(row['AA2']) + str(row['AA3'])
    if len(seq) == 3 and all(c in AMINO_ACIDS for c in seq):
        sequences_3pos.append(seq)

counts_3 = logomaker.alignment_to_matrix(sequences_3pos, to_type='counts')
info_3 = logomaker.transform_matrix(counts_3, from_type='counts', to_type='information')

fig, ax = plt.subplots(figsize=(2.6, 2.2))
logomaker.Logo(info_3, ax=ax, color_scheme=AA_COLOR_SCHEME, font_name='DejaVu Sans')
ax.set_ylabel('Bits', fontsize=7)
ax.set_xticks([0, 1, 2])
ax.set_xticklabels(['1 (Met)', '2', '3'], fontsize=7)
ax.tick_params(axis='y', labelsize=6)
ax.set_title('Sequence Logo — Best Hits', fontweight='bold', fontsize=8, pad=3)
legend_patches = [mpatches.Patch(color=CAT_COLORS[c], label=c) for c in
                  ['Acidic', 'Basic', 'Nonpolar', 'Polar']]
ax.legend(handles=legend_patches, fontsize=5, frameon=False,
          loc='upper left', bbox_to_anchor=(1.02, 1.0),
          handlelength=0.8, handletextpad=0.3, borderpad=0.2)
plt.tight_layout()
save(fig, 'fig4b_logo_with_Met')

# ============================================================
# FIG 5: Extended logo (positions 2-6), category-colored
# ============================================================
print("--- Fig 5: Extended logo ---")
sequences_5pos = []
for _, row in df_hits.iterrows():
    seq = ''.join([str(row[f'AA{i}']) for i in range(2, 7)])
    if len(seq) == 5 and all(c in AMINO_ACIDS for c in seq):
        sequences_5pos.append(seq)

counts_5 = logomaker.alignment_to_matrix(sequences_5pos, to_type='counts')
info_5 = logomaker.transform_matrix(counts_5, from_type='counts', to_type='information')

fig, ax = plt.subplots(figsize=(5, 2.2))
logomaker.Logo(info_5, ax=ax, color_scheme=AA_COLOR_SCHEME, font_name='DejaVu Sans')
ax.set_ylabel('Bits', fontsize=7)
ax.set_xticks(range(5))
ax.set_xticklabels([f'{i}' for i in range(2, 7)], fontsize=7)
ax.tick_params(axis='y', labelsize=6)
ax.set_xlabel('Position', fontsize=7)
ax.set_title('Sequence Logo — Best Hits (Pos 2-6)', fontweight='bold', fontsize=8, pad=3)
legend_patches = [mpatches.Patch(color=CAT_COLORS[c], label=c) for c in
                  ['Acidic', 'Basic', 'Nonpolar', 'Polar']]
ax.legend(handles=legend_patches, fontsize=5, frameon=False,
          loc='upper left', bbox_to_anchor=(1.02, 1.0),
          handlelength=0.8, handletextpad=0.3, borderpad=0.2)
plt.tight_layout()
save(fig, 'fig5_logo_extended_besthits')

# ============================================================
# FIG 5b: Extended logo WITH Met (positions 1-6)
# ============================================================
print("--- Fig 5b: Extended logo with Met ---")
sequences_6pos = []
for _, row in df_hits.iterrows():
    seq = 'M' + ''.join([str(row[f'AA{i}']) for i in range(2, 7)])
    if len(seq) == 6 and all(c in AMINO_ACIDS for c in seq):
        sequences_6pos.append(seq)

counts_6 = logomaker.alignment_to_matrix(sequences_6pos, to_type='counts')
info_6 = logomaker.transform_matrix(counts_6, from_type='counts', to_type='information')

fig, ax = plt.subplots(figsize=(5.5, 2.2))
logomaker.Logo(info_6, ax=ax, color_scheme=AA_COLOR_SCHEME, font_name='DejaVu Sans')
ax.set_ylabel('Bits', fontsize=7)
ax.set_xticks(range(6))
ax.set_xticklabels(['1 (Met)', '2', '3', '4', '5', '6'], fontsize=7)
ax.tick_params(axis='y', labelsize=6)
ax.set_xlabel('Position', fontsize=7)
ax.set_title('Sequence Logo — Best Hits (Pos 1-6)', fontweight='bold', fontsize=8, pad=3)
legend_patches = [mpatches.Patch(color=CAT_COLORS[c], label=c) for c in
                  ['Acidic', 'Basic', 'Nonpolar', 'Polar']]
ax.legend(handles=legend_patches, fontsize=5, frameon=False,
          loc='upper left', bbox_to_anchor=(1.02, 1.0),
          handlelength=0.8, handletextpad=0.3, borderpad=0.2)
plt.tight_layout()
save(fig, 'fig5b_logo_extended_with_Met')

# ============================================================
# FIG 9 (NEW): AA Category enrichment — grouped by property
# ============================================================
print("\n--- Fig 9: AA category enrichment ---")
cat_order = ['Acidic', 'Basic', 'Nonpolar', 'Polar']

def category_freq(freq_dict):
    cat_freq = {c: 0 for c in cat_order}
    for aa, f in freq_dict.items():
        cat_freq[AA_CATEGORIES.get(aa, 'Polar')] += f
    return cat_freq

lib2_cat = category_freq(lib_aa2_freq)
lib3_cat = category_freq(lib_aa3_freq)
hit2_cat = category_freq(hit_aa2_freq)
hit3_cat = category_freq(hit_aa3_freq)

fig, axes = plt.subplots(1, 2, figsize=(10, 5))
for ax, lib_c, hit_c, pos_label in zip(axes, [lib2_cat, lib3_cat], [hit2_cat, hit3_cat],
                                         ['Position 2', 'Position 3']):
    sorted_cats = sorted(cat_order, key=lambda c: hit_c[c], reverse=True)
    x = np.arange(len(sorted_cats))
    width = 0.35
    lib_vals = [lib_c[c] * 100 for c in sorted_cats]
    hit_vals = [hit_c[c] * 100 for c in sorted_cats]
    ax.bar(x - width/2, lib_vals, width, label='Library', color=COLOR_LIGHT, edgecolor='white')
    ax.bar(x + width/2, hit_vals, width, label='Best Hits', color=COLOR_MAIN, edgecolor='white')

    for i, cat in enumerate(sorted_cats):
        lib_v = lib_c[cat]
        hit_v = hit_c[cat]
        if lib_v > 0:
            fold = hit_v / lib_v
            ax.text(i + width/2, hit_c[cat] * 100 + 1, f'{fold:.1f}x',
                    ha='center', fontsize=9, fontweight='bold', color='#333')

    ax.set_xticks(x)
    ax.set_xticklabels(sorted_cats, fontsize=11)
    ax.set_ylabel('Frequency (%)')
    ax.set_title(f'{pos_label}', fontweight='bold')
    ax.legend(frameon=False, fontsize=9)

plt.suptitle('Amino Acid Category Enrichment (after Met)\nGrouped by Biochemical Property',
             fontsize=14, fontweight='bold', y=1.03)
plt.tight_layout()
save(fig, 'fig9_AA_category_enrichment')

# ============================================================
# FIG 10 (NEW): iceLogo-style — difference plot (hits − library)
# ============================================================
print("--- Fig 10: iceLogo-style difference plot ---")
fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
for ax, lib_f, hit_f, pos_label in zip(
    axes, [lib_aa2_freq, lib_aa3_freq], [hit_aa2_freq, hit_aa3_freq],
    ['Position 2', 'Position 3']
):
    diffs = {aa: (hit_f.get(aa, 0) - lib_f.get(aa, 0)) * 100 for aa in AMINO_ACIDS}
    sorted_aa = sorted(diffs.keys(), key=lambda a: diffs[a], reverse=True)
    vals = [diffs[aa] for aa in sorted_aa]
    colors = [COLOR_MAIN if v > 0 else '#AEB6BF' for v in vals]
    ax.bar(sorted_aa, vals, color=colors, edgecolor='white', linewidth=0.5)
    ax.axhline(y=0, color='black', linewidth=0.8)
    ax.set_xlabel('Amino Acid')
    ax.set_ylabel('Δ Frequency (%) [Hits − Library]')
    ax.set_title(f'{pos_label} — Enriched vs Depleted', fontweight='bold')

plt.suptitle('iceLogo-style: Frequency Difference (Best Hits − Library)\nPositive = enriched in hits, Negative = depleted',
             fontsize=14, fontweight='bold', y=1.03)
plt.tight_layout()
save(fig, 'fig10_iceLogo_difference')

# ============================================================
# FIG 11 (NEW): Statistical significance — Fisher exact test
# ============================================================
print("--- Fig 11: Statistical enrichment with p-values ---")
lib_aa2_counts = get_counts(df_all['AA2'])
lib_aa3_counts = get_counts(df_all['AA3'])
hit_aa2_counts = get_counts(df_hits['AA2'])
hit_aa3_counts = get_counts(df_hits['AA3'])
n_lib = len(df_all)
n_hits = len(df_hits)

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
for ax, lib_c, hit_c, enr, pos_label in zip(
    axes,
    [lib_aa2_counts, lib_aa3_counts],
    [hit_aa2_counts, hit_aa3_counts],
    [enrich_aa2, enrich_aa3],
    ['Position 2', 'Position 3']
):
    sorted_aa = sorted(AMINO_ACIDS, key=lambda a: enr[a], reverse=True)
    vals = [enr[aa] for aa in sorted_aa]
    pvals = []
    for aa in sorted_aa:
        table = [[hit_c[aa], n_hits - hit_c[aa]],
                 [lib_c[aa], n_lib - lib_c[aa]]]
        _, p = stats.fisher_exact(table, alternative='greater')
        pvals.append(p)

    colors = [COLOR_MAIN if p < 0.05 else COLOR_PALE for p in pvals]
    bars = ax.bar(sorted_aa, vals, color=colors, edgecolor='white', linewidth=0.5)
    ax.axhline(y=1.0, color='#7B7D7D', linestyle='--', linewidth=0.8, alpha=0.7)

    for i, (aa, p) in enumerate(zip(sorted_aa, pvals)):
        if p < 0.001:
            label = '***'
        elif p < 0.01:
            label = '**'
        elif p < 0.05:
            label = '*'
        else:
            label = ''
        if label:
            ax.text(i, vals[i] + 0.05, label, ha='center', fontsize=10, fontweight='bold', color='#333')

    ax.set_xlabel('Amino Acid')
    ax.set_ylabel('Enrichment (hits / library)')
    ax.set_title(f'{pos_label} — Statistical Enrichment\n(Fisher exact test, * p<0.05, ** p<0.01, *** p<0.001)',
                 fontweight='bold', fontsize=11)
    ax.set_ylim(0, max(vals) * 1.2)

    sig_patch = mpatches.Patch(color=COLOR_MAIN, label='p < 0.05')
    ns_patch = mpatches.Patch(color=COLOR_PALE, label='n.s.')
    ax.legend(handles=[sig_patch, ns_patch], frameon=False, fontsize=9)

plt.suptitle('Statistically Significant Amino Acid Enrichment at N-terminal Positions',
             fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
save(fig, 'fig11_statistical_enrichment')

# ============================================================
# FIG 12 (NEW): Probability weight matrix heatmap (pos 2-8)
# ============================================================
print("--- Fig 12: Probability weight matrix heatmap ---")
n_pos = 7
prob_matrix = np.zeros((20, n_pos))
for col_idx, pos in enumerate(range(2, 2 + n_pos)):
    col = f'AA{pos}'
    counts = Counter(df_hits[col].dropna())
    total = sum(counts.values())
    for aa_idx, aa in enumerate(AMINO_ACIDS):
        prob_matrix[aa_idx, col_idx] = counts.get(aa, 0) / total if total > 0 else 0

fig, ax = plt.subplots(figsize=(6, 9))
im = ax.imshow(prob_matrix, cmap=cmap_warm, aspect='auto', vmin=0, vmax=prob_matrix.max())
ax.set_xticks(range(n_pos))
ax.set_xticklabels([f'Pos {i}' for i in range(2, 2 + n_pos)], fontsize=11)
ax.set_yticks(range(20))
ax.set_yticklabels(AMINO_ACIDS, fontsize=11)
for i in range(20):
    for j in range(n_pos):
        v = prob_matrix[i, j]
        if v >= 0.04:
            color = 'white' if v > 0.15 else 'black'
            ax.text(j, i, f'{v:.0%}', ha='center', va='center', fontsize=8, color=color)
plt.colorbar(im, ax=ax, shrink=0.5, label='Probability')
ax.set_title('Position Weight Matrix — Best Hits\n(Probability of each AA at positions 2-8)',
             fontweight='bold')
ax.set_xlabel('Position (after Met)')
ax.set_ylabel('Amino Acid')
plt.tight_layout()
save(fig, 'fig12_PWM_heatmap')

# ============================================================
# FIG 13 (NEW): Consensus vs. anti-consensus — top enriched/depleted
# ============================================================
print("--- Fig 13: Consensus dipeptide top/bottom ---")
all_dipep_enrich = {}
for dp in set(list(lib_dipep.keys()) + list(hit_dipep.keys())):
    if lib_dipep.get(dp, 0) > 0 and hit_dipep.get(dp, 0) > 0:
        all_dipep_enrich[dp] = hit_dipep[dp] / lib_dipep[dp]
    elif hit_dipep.get(dp, 0) > 0:
        all_dipep_enrich[dp] = float('inf')

finite_enrich = {k: v for k, v in all_dipep_enrich.items() if np.isfinite(v)}
top10 = sorted(finite_enrich, key=finite_enrich.get, reverse=True)[:10]
bot10 = sorted(finite_enrich, key=finite_enrich.get)[:10]

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

ax = axes[0]
vals = [finite_enrich[dp] for dp in top10]
ax.barh(range(len(top10)), vals, color=plt.cm.YlOrRd(np.linspace(0.5, 0.95, 10))[::-1],
        edgecolor='white')
ax.set_yticks(range(len(top10)))
ax.set_yticklabels(top10, fontsize=11)
ax.set_xlabel('Enrichment')
ax.set_title('Top 10 Enriched Dipeptides\n(Consensus)', fontweight='bold')
ax.invert_yaxis()

ax = axes[1]
vals = [finite_enrich[dp] for dp in bot10]
ax.barh(range(len(bot10)), vals, color='#AEB6BF', edgecolor='white')
ax.set_yticks(range(len(bot10)))
ax.set_yticklabels(bot10, fontsize=11)
ax.set_xlabel('Enrichment')
ax.set_title('Top 10 Depleted Dipeptides\n(Anti-consensus)', fontweight='bold')
ax.axvline(x=1.0, color='#7B7D7D', linestyle='--', linewidth=0.8)
ax.invert_yaxis()

plt.suptitle('Consensus vs Anti-consensus Dipeptides at Positions 2-3',
             fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
save(fig, 'fig13_consensus_vs_anticonsensus')

# ============================================================
# FIGURES 6-8: PSI/AAVS (unchanged from before)
# ============================================================
print("\n--- Figs 6-8: PSI/AAVS stability analysis ---")
psi_aavs = df_all[['PSI-293a', 'PSI-293b']].mean(axis=1)
df_all['PSI_AAVS'] = psi_aavs
target_dipeps = ['PD', 'PE', 'GE']
position_data = {dp: [] for dp in target_dipeps}

for pos in range(1, 24):
    col1 = f'AA{pos}'
    col2 = f'AA{pos+1}'
    for dp in target_dipeps:
        aa1, aa2 = dp[0], dp[1]
        mask = (df_all[col1] == aa1) & (df_all[col2] == aa2)
        count = mask.sum()
        mean_psi = df_all.loc[mask, 'PSI_AAVS'].mean() if count > 0 else np.nan
        position_data[dp].append({
            'position': f'{pos}-{pos+1}', 'pos_start': pos,
            'count': count, 'mean_PSI': mean_psi
        })

dp_colors = {'PD': '#C0392B', 'PE': '#E67E22', 'GE': '#F39C12'}

fig, axes = plt.subplots(2, 1, figsize=(14, 9), gridspec_kw={'height_ratios': [2, 1]})
ax = axes[0]
x = np.arange(1, 24)
width = 0.25
offsets = {'PD': -width, 'PE': 0, 'GE': width}
for dp in target_dipeps:
    vals = [d['mean_PSI'] for d in position_data[dp]]
    ax.bar(x + offsets[dp], vals, width, label=dp, color=dp_colors[dp],
           edgecolor='white', linewidth=0.5, alpha=0.85)
ax.set_xlabel('Dipeptide Starting Position in 24-mer')
ax.set_ylabel('Mean PSI (AAVS)')
ax.set_title('Stability (PSI AAVS) by Dipeptide Position in 24-mers\nP-D, P-E, and G-E motifs',
             fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels([f'{i}-{i+1}' for i in range(1, 24)], rotation=45, ha='right', fontsize=9)
ax.legend(frameon=False, fontsize=11)
ax.axhline(y=psi_aavs.mean(), color='#7B7D7D', linestyle='--', linewidth=0.8, alpha=0.6)
ax.annotate('Global mean PSI', xy=(22, psi_aavs.mean()), fontsize=9, color='#7B7D7D', va='bottom')

ax2 = axes[1]
for dp in target_dipeps:
    counts = [d['count'] for d in position_data[dp]]
    ax2.bar(x + offsets[dp], counts, width, label=dp, color=dp_colors[dp],
            edgecolor='white', linewidth=0.5, alpha=0.85)
ax2.set_xlabel('Dipeptide Starting Position in 24-mer')
ax2.set_ylabel('Count')
ax2.set_title('Number of Peptides with Motif at Each Position', fontweight='bold')
ax2.set_xticks(x)
ax2.set_xticklabels([f'{i}-{i+1}' for i in range(1, 24)], rotation=45, ha='right', fontsize=9)
ax2.legend(frameon=False, fontsize=11)
plt.tight_layout()
save(fig, 'fig6_PSI_AAVS_dipeptide_position')

fig, ax = plt.subplots(figsize=(10, 5))
for dp in target_dipeps:
    psi_by_pos = [d for d in position_data[dp]
                  if not np.isnan(d['mean_PSI']) and d['count'] >= 5]
    positions = [d['pos_start'] for d in psi_by_pos]
    psis = [d['mean_PSI'] for d in psi_by_pos]
    ax.plot(positions, psis, 'o-', color=dp_colors[dp], label=dp, linewidth=2,
            markersize=7, alpha=0.85)
    for d in psi_by_pos:
        if d['pos_start'] in [2, 3]:
            ax.annotate(f"n={d['count']}\nPSI={d['mean_PSI']:.2f}",
                        xy=(d['pos_start'], d['mean_PSI']),
                        xytext=(d['pos_start'] + 0.5, d['mean_PSI'] + 0.08),
                        fontsize=8, color=dp_colors[dp],
                        arrowprops=dict(arrowstyle='->', color=dp_colors[dp], lw=0.8))
ax.axhline(y=psi_aavs.mean(), color='#7B7D7D', linestyle='--', linewidth=0.8, alpha=0.6)
ax.axvspan(1.5, 3.5, alpha=0.1, color=COLOR_MAIN, label='Positions 2-3 (highlighted)')
ax.set_xlabel('Dipeptide Starting Position')
ax.set_ylabel('Mean PSI (AAVS)')
ax.set_title('Position-dependent Stability Effect of P-D, P-E, G-E Motifs\n(filtered for n>=5)',
             fontweight='bold')
ax.legend(frameon=False)
ax.set_xticks(range(1, 24))
ax.set_xticklabels(range(1, 24))
plt.tight_layout()
save(fig, 'fig7_PSI_position_effect_line')

fig, ax = plt.subplots(figsize=(7, 5))
bar_data = []
for dp in target_dipeps:
    pos23_psis, other_psis = [], []
    for pos in range(1, 24):
        col1, col2 = f'AA{pos}', f'AA{pos+1}'
        aa1, aa2 = dp[0], dp[1]
        mask = (df_all[col1] == aa1) & (df_all[col2] == aa2)
        if mask.sum() > 0:
            if pos in [2, 3]:
                pos23_psis.extend(df_all.loc[mask, 'PSI_AAVS'].tolist())
            else:
                other_psis.extend(df_all.loc[mask, 'PSI_AAVS'].tolist())
    bar_data.append({
        'dipeptide': dp,
        'pos_2_3': np.mean(pos23_psis) if pos23_psis else np.nan,
        'other_pos': np.mean(other_psis) if other_psis else np.nan,
        'n_23': len(pos23_psis), 'n_other': len(other_psis)
    })

x = np.arange(len(target_dipeps))
width = 0.35
ax.bar(x - width/2, [d['pos_2_3'] for d in bar_data], width, label='Position 2-3',
       color=COLOR_MAIN, edgecolor='white')
ax.bar(x + width/2, [d['other_pos'] for d in bar_data], width, label='Other positions',
       color=COLOR_LIGHT, edgecolor='white')
for i, d in enumerate(bar_data):
    ax.text(i - width/2, d['pos_2_3'] + 0.02, f"n={d['n_23']}", ha='center', fontsize=9, color='#333')
    ax.text(i + width/2, d['other_pos'] + 0.02, f"n={d['n_other']}", ha='center', fontsize=9, color='#333')
ax.set_xlabel('Dipeptide Motif')
ax.set_ylabel('Mean PSI (AAVS)')
ax.set_title('PSI (AAVS) Stability: Position 2-3 vs Other Positions\nfor P-D, P-E, and G-E Motifs',
             fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels([d['dipeptide'] for d in bar_data], fontsize=12)
ax.legend(frameon=False)
plt.tight_layout()
save(fig, 'fig8_PSI_pos23_vs_other')

# ============================================================
# FIG 14: Enrichment heatmap across ALL positions (2-12) — hits vs library
# ============================================================
print("\n--- Fig 14: Multi-position enrichment heatmap ---")
n_pos_wide = 11
enrich_wide = np.zeros((20, n_pos_wide))
pval_wide = np.ones((20, n_pos_wide))
n_lib_total = len(df_all)
n_hit_total = len(df_hits)
for col_idx, pos in enumerate(range(2, 2 + n_pos_wide)):
    col = f'AA{pos}'
    lib_counts = Counter(df_all[col].dropna())
    lib_total = sum(lib_counts.values())
    hit_counts = Counter(df_hits[col].dropna())
    hit_total = sum(hit_counts.values())
    for aa_idx, aa in enumerate(AMINO_ACIDS):
        lib_f = lib_counts.get(aa, 0) / lib_total if lib_total > 0 else 0
        hit_f = hit_counts.get(aa, 0) / hit_total if hit_total > 0 else 0
        enrich_wide[aa_idx, col_idx] = (hit_f / lib_f) if lib_f > 0 else 0
        table = [[hit_counts.get(aa, 0), hit_total - hit_counts.get(aa, 0)],
                 [lib_counts.get(aa, 0), lib_total - lib_counts.get(aa, 0)]]
        _, pval_wide[aa_idx, col_idx] = stats.fisher_exact(table, alternative='two-sided')

fig, ax = plt.subplots(figsize=(9, 9))
im = ax.imshow(enrich_wide, cmap=cmap_warm, aspect='auto', vmin=0,
               vmax=np.clip(np.nanpercentile(enrich_wide[enrich_wide > 0], 95), 3, 6))
ax.set_xticks(range(n_pos_wide))
ax.set_xticklabels([f'{i}' for i in range(2, 2 + n_pos_wide)], fontsize=11)
ax.set_yticks(range(20))
ax.set_yticklabels(AMINO_ACIDS, fontsize=11)
for i in range(20):
    for j in range(n_pos_wide):
        v = enrich_wide[i, j]
        p = pval_wide[i, j]
        stars = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
        if v >= 1.5:
            color = 'white' if v > 3 else 'black'
            ax.text(j, i, f'{v:.1f}{stars}', ha='center', va='center', fontsize=7,
                    color=color, fontweight='bold')
        elif stars:
            ax.text(j, i, f'{v:.1f}{stars}', ha='center', va='center', fontsize=6,
                    color='black', fontweight='bold')
        elif v > 0:
            ax.text(j, i, f'{v:.1f}', ha='center', va='center', fontsize=6, color='#666')
plt.colorbar(im, ax=ax, shrink=0.5, label='Enrichment (hits / library)')
ax.set_xlabel('Position in peptide', fontsize=12)
ax.set_ylabel('Amino Acid', fontsize=12)
ax.set_title('Amino Acid Enrichment Heatmap (Best Hits vs Library)\nPositions 2-12',
             fontweight='bold')
ax.text(0.5, -0.05, '* p<0.05  ** p<0.01  *** p<0.001 (Fisher exact test)',
        transform=ax.transAxes, ha='center', fontsize=9, color='#555')
plt.tight_layout()
save(fig, 'fig14_enrichment_heatmap_pos2_12')

# ============================================================
# FIG 15: Log2 enrichment heatmap (pos 2-12) — diverging
# ============================================================
print("--- Fig 15: Log2 enrichment heatmap ---")
from matplotlib.colors import TwoSlopeNorm
log2_enrich = np.where(enrich_wide > 0, np.log2(enrich_wide), np.nan)
log2_enrich_masked = np.ma.masked_invalid(log2_enrich)

cmap_div = LinearSegmentedColormap.from_list('div_warm',
    ['#2E4057', '#AEB6BF', '#FDFEFE', '#F5B041', '#C0392B', '#78281F'])
vabs = np.nanmax(np.abs(log2_enrich[np.isfinite(log2_enrich)]))
vabs = min(vabs, 4)
norm = TwoSlopeNorm(vmin=-vabs, vcenter=0, vmax=vabs)

fig, ax = plt.subplots(figsize=(9, 9))
im = ax.imshow(log2_enrich_masked, cmap=cmap_div, aspect='auto', norm=norm)
ax.set_xticks(range(n_pos_wide))
ax.set_xticklabels([f'{i}' for i in range(2, 2 + n_pos_wide)], fontsize=11)
ax.set_yticks(range(20))
ax.set_yticklabels(AMINO_ACIDS, fontsize=11)
for i in range(20):
    for j in range(n_pos_wide):
        v = log2_enrich[i, j]
        if np.isfinite(v) and abs(v) >= 0.7:
            color = 'white' if abs(v) > 2 else 'black'
            ax.text(j, i, f'{v:.1f}', ha='center', va='center', fontsize=7,
                    color=color, fontweight='bold')
plt.colorbar(im, ax=ax, shrink=0.5, label='Log2 Enrichment')
ax.set_xlabel('Position in peptide', fontsize=12)
ax.set_ylabel('Amino Acid', fontsize=12)
ax.set_title('Log2 Enrichment Heatmap (Best Hits vs Library)\nPositions 2-12\nRed = enriched, Blue = depleted',
             fontweight='bold')
plt.tight_layout()
save(fig, 'fig15_log2_enrichment_heatmap')

# ============================================================
# FIG 16: AA category enrichment heatmap across positions 2-12
# ============================================================
print("--- Fig 16: Category enrichment heatmap across positions ---")
cat_names = ['Acidic', 'Basic', 'Nonpolar', 'Polar']
cat_enrich_matrix = np.zeros((len(cat_names), n_pos_wide))

for col_idx, pos in enumerate(range(2, 2 + n_pos_wide)):
    col = f'AA{pos}'
    lib_counts = Counter(df_all[col].dropna())
    lib_total = sum(lib_counts.values())
    hit_counts = Counter(df_hits[col].dropna())
    hit_total = sum(hit_counts.values())
    for cat_idx, cat in enumerate(cat_names):
        cat_aas = [aa for aa, c in AA_CATEGORIES.items() if c == cat]
        lib_f = sum(lib_counts.get(aa, 0) for aa in cat_aas) / lib_total if lib_total > 0 else 0
        hit_f = sum(hit_counts.get(aa, 0) for aa in cat_aas) / hit_total if hit_total > 0 else 0
        cat_enrich_matrix[cat_idx, col_idx] = (hit_f / lib_f) if lib_f > 0 else 0

fig, ax = plt.subplots(figsize=(9, 3.5))
im = ax.imshow(cat_enrich_matrix, cmap=cmap_warm, aspect='auto', vmin=0,
               vmax=max(cat_enrich_matrix.max(), 3.5))
ax.set_xticks(range(n_pos_wide))
ax.set_xticklabels([f'{i}' for i in range(2, 2 + n_pos_wide)], fontsize=11)
ax.set_yticks(range(len(cat_names)))
ax.set_yticklabels(cat_names, fontsize=12)
for i in range(len(cat_names)):
    for j in range(n_pos_wide):
        v = cat_enrich_matrix[i, j]
        color = 'white' if v > 2.5 else 'black'
        ax.text(j, i, f'{v:.1f}x', ha='center', va='center', fontsize=9, color=color,
                fontweight='bold' if v >= 1.5 else 'normal')
plt.colorbar(im, ax=ax, shrink=0.7, label='Enrichment')
ax.set_xlabel('Position in peptide', fontsize=12)
ax.set_title('Amino Acid Category Enrichment Across Positions 2-12\n(Best Hits vs Library)',
             fontweight='bold')
plt.tight_layout()
save(fig, 'fig16_category_enrichment_heatmap')

# ============================================================
# FIG 17: Frequency heatmap — Best Hits vs Library side-by-side
# ============================================================
print("--- Fig 17: Frequency heatmap side-by-side ---")
n_pos_freq = 6
freq_lib = np.zeros((20, n_pos_freq))
freq_hit = np.zeros((20, n_pos_freq))
for col_idx, pos in enumerate(range(2, 2 + n_pos_freq)):
    col = f'AA{pos}'
    lib_counts = Counter(df_all[col].dropna())
    lib_total = sum(lib_counts.values())
    hit_counts = Counter(df_hits[col].dropna())
    hit_total = sum(hit_counts.values())
    for aa_idx, aa in enumerate(AMINO_ACIDS):
        freq_lib[aa_idx, col_idx] = lib_counts.get(aa, 0) / lib_total if lib_total > 0 else 0
        freq_hit[aa_idx, col_idx] = hit_counts.get(aa, 0) / hit_total if hit_total > 0 else 0

fig, axes = plt.subplots(1, 2, figsize=(12, 8))
for ax, data, title in zip(axes, [freq_lib, freq_hit],
                            ['Library (16,514 peptides)', 'Best Hits (54 peptides)']):
    im = ax.imshow(data, cmap=cmap_warm, aspect='auto', vmin=0, vmax=0.27)
    ax.set_xticks(range(n_pos_freq))
    ax.set_xticklabels([f'{i}' for i in range(2, 2 + n_pos_freq)], fontsize=10)
    ax.set_yticks(range(20))
    ax.set_yticklabels(AMINO_ACIDS, fontsize=10)
    for i in range(20):
        for j in range(n_pos_freq):
            v = data[i, j]
            if v >= 0.03:
                color = 'white' if v > 0.15 else 'black'
                ax.text(j, i, f'{v:.0%}', ha='center', va='center', fontsize=7, color=color)
    ax.set_title(title, fontweight='bold')
    ax.set_xlabel('Position')
    ax.set_ylabel('Amino Acid')
fig.subplots_adjust(right=0.88)
cbar_ax = fig.add_axes([0.90, 0.25, 0.02, 0.5])
fig.colorbar(im, cax=cbar_ax, label='Frequency')
plt.suptitle('Amino Acid Frequency Heatmap: Library vs Best Hits\nPositions 2-7',
             fontweight='bold', fontsize=14, y=1.01)
plt.tight_layout(rect=[0, 0, 0.88, 0.97])
save(fig, 'fig17_frequency_heatmap_lib_vs_hits')

# ============================================================
# FIG 18: Dipeptide enrichment heatmap — focused on enriched AA2
# ============================================================
print("--- Fig 18: Focused dipeptide heatmap (top AA2 residues) ---")
top_aa2 = ['P', 'G', 'T', 'Q', 'E', 'D']
heatmap_focused = np.zeros((len(top_aa2), 20))
for i, aa2 in enumerate(top_aa2):
    for j, aa3 in enumerate(AMINO_ACIDS):
        dp = aa2 + aa3
        lib_f = lib_dipep.get(dp, 0)
        hit_f = hit_dipep.get(dp, 0)
        heatmap_focused[i, j] = (hit_f / lib_f) if lib_f > 0 else 0

fig, ax = plt.subplots(figsize=(12, 4.5))
im = ax.imshow(heatmap_focused, cmap=cmap_warm, aspect='auto', vmin=0,
               vmax=min(heatmap_focused.max(), 45))
ax.set_xticks(range(20))
ax.set_xticklabels(AMINO_ACIDS, fontsize=11)
ax.set_yticks(range(len(top_aa2)))
ax.set_yticklabels(top_aa2, fontsize=13, fontweight='bold')
for i in range(len(top_aa2)):
    for j in range(20):
        v = heatmap_focused[i, j]
        if v >= 3:
            color = 'white' if v > 15 else 'black'
            ax.text(j, i, f'{v:.0f}x', ha='center', va='center', fontsize=8,
                    color=color, fontweight='bold')
plt.colorbar(im, ax=ax, shrink=0.7, label='Enrichment (hits / library)')
ax.set_xlabel('Position 3 (AA3)', fontsize=12)
ax.set_ylabel('Position 2 (AA2)', fontsize=12)
ax.set_title('Dipeptide Enrichment — Top Enriched Position 2 Residues\n(P, G, T, Q, E, D) × All AA at Position 3',
             fontweight='bold')
plt.tight_layout()
save(fig, 'fig18_focused_dipeptide_heatmap')

# ============================================================
# FIG 19: Volcano plot — UBR3 screen (delta PSI vs significance)
# ============================================================
print("\n--- Fig 19: Volcano plot ---")
mean_ctrl = (df_all['PSI-293a'] + df_all['PSI-293b']) / 2
mean_ubr3 = (df_all['PSI-UBR3a'] + df_all['PSI-UBR3b']) / 2
delta_psi = mean_ubr3 - mean_ctrl

pvals = []
for _, row in df_all.iterrows():
    ctrl = [row['PSI-293a'], row['PSI-293b']]
    ubr3 = [row['PSI-UBR3a'], row['PSI-UBR3b']]
    _, p = stats.ttest_ind(ubr3, ctrl, equal_var=False)
    pvals.append(p if np.isfinite(p) and p > 0 else 1.0)

pvals = np.array(pvals)
neg_log_p = -np.log10(pvals)
neg_log_p = np.clip(neg_log_p, 0, 50)

DPSI_CUTOFF = 0.5
P_CUTOFF = 0.05
neg_log_p_cutoff = -np.log10(P_CUTOFF)

is_hit = df_all['HITubr3'].values
sig_up = (delta_psi >= DPSI_CUTOFF) & (pvals < P_CUTOFF) & ~is_hit
sig_down = (delta_psi <= -DPSI_CUTOFF) & (pvals < P_CUTOFF) & ~is_hit
ns = ~sig_up & ~sig_down & ~is_hit

fig, ax = plt.subplots(figsize=(16, 7))
ax.scatter(delta_psi[ns], neg_log_p[ns], s=4, c='#BDC3C7', alpha=0.3,
           edgecolors='none', rasterized=True, label='Not significant')
ax.scatter(delta_psi[sig_down], neg_log_p[sig_down], s=6, c='#2E86C1', alpha=0.5,
           edgecolors='none', rasterized=True, label=f'dPSI <= -{DPSI_CUTOFF}')
ax.scatter(delta_psi[sig_up], neg_log_p[sig_up], s=6, c=COLOR_SEC, alpha=0.5,
           edgecolors='none', rasterized=True, label=f'dPSI >= {DPSI_CUTOFF}')
ax.scatter(delta_psi[is_hit], neg_log_p[is_hit], s=25, c=COLOR_MAIN, alpha=0.8,
           edgecolors='black', linewidths=0.4, zorder=5, label=f'Best hits (n={is_hit.sum()})')

ax.set_ylim(-0.1, 3.0)

hit_idx = np.where(is_hit)[0]
from adjustText import adjust_text
texts = []
for idx in hit_idx:
    gene = df_all.iloc[idx]['Gene_ID']
    texts.append(ax.text(delta_psi.iloc[idx], neg_log_p[idx], gene,
                         fontsize=7.5, color='#333', fontweight='bold', zorder=6))
adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='#999', lw=0.4),
            force_text=(0.6, 0.9), force_points=(0.3, 0.3), expand=(1.3, 1.5))

ax.axvline(x=DPSI_CUTOFF, color='#333', linestyle='--', linewidth=0.8, alpha=0.6)
ax.axvline(x=-DPSI_CUTOFF, color='#333', linestyle='--', linewidth=0.8, alpha=0.6)
ax.axhline(y=neg_log_p_cutoff, color='#333', linestyle='--', linewidth=0.8, alpha=0.6)

ax.text(DPSI_CUTOFF + 0.02, 2.95, f'dPSI = {DPSI_CUTOFF}',
        fontsize=8, color='#333', va='top')
ax.text(-DPSI_CUTOFF - 0.02, 2.95, f'dPSI = -{DPSI_CUTOFF}',
        fontsize=8, color='#333', va='top', ha='right')
ax.text(ax.get_xlim()[1] * 0.98, neg_log_p_cutoff + 0.05, f'p = {P_CUTOFF}',
        fontsize=8, color='#333', va='bottom', ha='right')

ax.set_xlabel(r'$\Delta$ PSI (UBR3 $-$ Control)', fontsize=12)
ax.set_ylabel('-log10(p-value)', fontsize=12)
ax.set_title('Volcano Plot — UBR3 N-terminal Screen\n'
             f'(n={len(df_all):,} peptides, t-test UBR3 vs Control replicates)',
             fontweight='bold', fontsize=12)
ax.legend(fontsize=9, frameon=False, loc='upper left', markerscale=1.5)
plt.tight_layout()
save(fig, 'fig19_volcano_plot')

# FIG 20: Volcano plot — colored by Pos2 residue motif
# ============================================================
print("\n--- Fig 20: Volcano plot (colored by motif) ---")
from adjustText import adjust_text

pos2_colors = {
    'G': '#2ECC71',   # green
    'P': '#E74C3C',   # red
    'E': '#3498DB',   # blue
    'T': '#F39C12',   # amber
    'D': '#9B59B6',   # purple
    'A': '#1ABC9C',   # teal
    'S': '#E67E22',   # orange
    'Q': '#34495E',   # dark gray-blue
}
DEFAULT_COLOR = '#95A5A6'

fig, ax = plt.subplots(figsize=(16, 7))
ax.scatter(delta_psi[ns], neg_log_p[ns], s=4, c='#BDC3C7', alpha=0.3,
           edgecolors='none', rasterized=True, label='Not significant')
ax.scatter(delta_psi[sig_down], neg_log_p[sig_down], s=6, c='#BDC3C7', alpha=0.3,
           edgecolors='none', rasterized=True)
ax.scatter(delta_psi[sig_up & ~is_hit], neg_log_p[sig_up & ~is_hit], s=6, c='#BDC3C7',
           alpha=0.3, edgecolors='none', rasterized=True)

hit_idx = np.where(is_hit)[0]
plotted_pos2 = set()
for idx in hit_idx:
    seq = df_all.iloc[idx]['AA_seq']
    p2 = seq[1]
    c = pos2_colors.get(p2, DEFAULT_COLOR)
    lbl = f'Pos2 = {p2} (n={sum(1 for i in hit_idx if df_all.iloc[i]["AA_seq"][1]==p2)})' if p2 not in plotted_pos2 else None
    ax.scatter(delta_psi.iloc[idx], neg_log_p[idx], s=50, c=c, alpha=0.85,
               edgecolors='black', linewidths=0.5, zorder=5, label=lbl)
    plotted_pos2.add(p2)

ax.set_ylim(-0.1, 3.0)

texts = []
for idx in hit_idx:
    seq = df_all.iloc[idx]['AA_seq']
    gene = df_all.iloc[idx]['Gene_ID']
    dipeptide = seq[1:3]
    p2 = seq[1]
    c = pos2_colors.get(p2, DEFAULT_COLOR)
    texts.append(ax.text(delta_psi.iloc[idx], neg_log_p[idx],
                         f'{gene} ({dipeptide})',
                         fontsize=7.5, color=c, fontweight='bold', zorder=6))
adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='#999', lw=0.4),
            force_text=(0.6, 0.9), force_points=(0.3, 0.3), expand=(1.3, 1.5))

ax.axvline(x=DPSI_CUTOFF, color='#333', linestyle='--', linewidth=0.8, alpha=0.6)
ax.axvline(x=-DPSI_CUTOFF, color='#333', linestyle='--', linewidth=0.8, alpha=0.6)
ax.axhline(y=neg_log_p_cutoff, color='#333', linestyle='--', linewidth=0.8, alpha=0.6)

ax.text(DPSI_CUTOFF + 0.02, 2.95, f'dPSI = {DPSI_CUTOFF}',
        fontsize=8, color='#333', va='top')
ax.text(-DPSI_CUTOFF - 0.02, 2.95, f'dPSI = -{DPSI_CUTOFF}',
        fontsize=8, color='#333', va='top', ha='right')
ax.text(ax.get_xlim()[1] * 0.98, neg_log_p_cutoff + 0.05, f'p = {P_CUTOFF}',
        fontsize=8, color='#333', va='bottom', ha='right')

ax.set_xlabel(r'$\Delta$ PSI (UBR3 $-$ Control)', fontsize=12)
ax.set_ylabel('-log10(p-value)', fontsize=12)
ax.set_title('Volcano Plot — Best Hits Colored by Position 2 Residue\n'
             f'(n={len(df_all):,} peptides, label shows Gene (Pos2–Pos3 dipeptide))',
             fontweight='bold', fontsize=12)
ax.legend(fontsize=9, frameon=True, loc='upper left', markerscale=1.2,
          edgecolor='#ccc', fancybox=True)
plt.tight_layout()
save(fig, 'fig20_volcano_motif')

print(f"\n=== All figures saved to: {OUTPUT_DIR} ===")
