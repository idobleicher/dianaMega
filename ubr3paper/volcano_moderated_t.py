"""
Volcano plot for the UBR3 N-terminal screen using a MODERATED t-test
(limma-style empirical Bayes variance shrinkage, Smyth 2004) with
Benjamini-Hochberg FDR correction.

Implements:
  1. Empirical Bayes variance shrinkage (Smyth, Stat. Appl. Genet. Mol. Biol. 2004)
  2. Moderated t-statistic  t_g = d_g / ( s_tilde_g * sqrt(1/n1 + 1/n2) )
  3. Benjamini-Hochberg FDR
  4. As a cross-check, also computes a Z-score of dPSI against the full
     16,514-gene background distribution.

Output: volcano plot colored by Position-2 residue, with labels for the
best hits.  Saves PNG + PDF + SVG, and an Excel table of all results.
"""
import os, sys, io, warnings
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
warnings.filterwarnings('ignore')

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from scipy.special import digamma, polygamma
from statsmodels.stats.multitest import multipletests
from adjustText import adjust_text

EXCEL_PATH = r'c:\Users\User\Downloads\UBR3 Nt screen (1).xlsx'
OUT_DIR    = r'c:\Users\User\Desktop\תינוקת\dianaMega\ubr3enrichmentlogo\figures'
os.makedirs(OUT_DIR, exist_ok=True)

# ══════════════════════════════════════════════════════════════════════════
# 1.  LOAD DATA
# ══════════════════════════════════════════════════════════════════════════
df_all  = pd.read_excel(EXCEL_PATH, sheet_name='Nprot_5_analyzed')
print(f"Loaded {len(df_all):,} peptides")

ctrl_a = df_all['PSI-293a'].values.astype(float)
ctrl_b = df_all['PSI-293b'].values.astype(float)
ubr3_a = df_all['PSI-UBR3a'].values.astype(float)
ubr3_b = df_all['PSI-UBR3b'].values.astype(float)

mean_ctrl = (ctrl_a + ctrl_b) / 2
mean_ubr3 = (ubr3_a + ubr3_b) / 2
delta_psi = mean_ubr3 - mean_ctrl

# Per-gene POOLED variance with n1=n2=2, equal variance assumption
# s_g^2 = [ (x1-xbar)^2 + (x2-xbar)^2 + (y1-ybar)^2 + (y2-ybar)^2 ] / (n1+n2-2)
# For n=2 each group:  var(x) = (x1-x2)^2 / 2 ;  similarly for y
# Pooled:  s^2 = ( (x1-x2)^2/2 + (y1-y2)^2/2 ) / (1+1) * 1   -- careful
# Correct formula (df = 2):
#   s_g^2 = [ (n1-1)*var(x) + (n2-1)*var(y) ] / (n1+n2-2)
#         = [ var(x) + var(y) ] / 2
var_ctrl = np.var(np.column_stack([ctrl_a, ctrl_b]), axis=1, ddof=1)
var_ubr3 = np.var(np.column_stack([ubr3_a, ubr3_b]), axis=1, ddof=1)
s2_gene  = (var_ctrl + var_ubr3) / 2          # pooled variance per gene
df_gene  = 2                                   # (n1-1)+(n2-1) = 1+1 = 2
n1, n2   = 2, 2

print(f"Genes with zero variance (both replicates identical in both groups): "
      f"{(s2_gene == 0).sum()}")

# ══════════════════════════════════════════════════════════════════════════
# 2.  EMPIRICAL BAYES VARIANCE SHRINKAGE  (Smyth 2004)
# ══════════════════════════════════════════════════════════════════════════
def trigamma(x):
    return polygamma(1, x)

def trigamma_inverse(x):
    """
    Newton-Raphson inverse of the trigamma function.
    Returns y such that trigamma(y) = x.
    Follows Smyth's limma implementation (R code, trigammaInverse).
    """
    x = np.asarray(x, dtype=float)
    y = np.where(x > 1e7, 1.0 / np.sqrt(np.maximum(x, 1e-12)),
                 0.5 + 1.0 / x)
    for _ in range(200):
        tri = trigamma(y)
        tet = polygamma(2, y)          # tetragamma (negative)
        delta = tri * (1.0 - tri / x) / tet
        y_new = y + delta
        y_new = np.where(y_new <= 0, y / 2.0, y_new)
        if np.max(np.abs(y_new - y)) < 1e-8:
            y = y_new
            break
        y = y_new
    return y

def fit_fdist(s2, df):
    """
    Fit a scaled F-distribution to s2 values (Smyth 2004).
    Returns (s0_sq, d0) - the prior variance and prior degrees of freedom.
    """
    s2 = np.asarray(s2, dtype=float)
    s2 = s2[s2 > 0]                    # exclude zero-variance genes
    if len(s2) < 10:
        raise RuntimeError("Too few non-zero variances to fit EB prior")
    z = np.log(s2)
    e = z - digamma(df / 2.0) + np.log(df / 2.0)
    mean_e = np.mean(e)
    var_e  = np.var(e, ddof=1)

    target = var_e - trigamma(df / 2.0)
    if target <= 0:
        # Prior df is effectively infinite (all variances equal)
        d0    = np.inf
        s0_sq = float(np.exp(mean_e))
    else:
        d0    = float(2.0 * trigamma_inverse(target))
        s0_sq = float(np.exp(mean_e + digamma(d0 / 2.0) - np.log(d0 / 2.0)))
    return s0_sq, d0

# Replace zero variances with the minimum non-zero variance (avoid log(0))
s2_for_eb = s2_gene.copy()
min_nz = s2_for_eb[s2_for_eb > 0].min()
s2_for_eb[s2_for_eb == 0] = min_nz

s0_sq, d0 = fit_fdist(s2_for_eb, df_gene)
print(f"\nEmpirical Bayes hyperparameters:")
print(f"  Prior variance s0^2 = {s0_sq:.4f}")
print(f"  Prior d.f.     d0   = {d0:.2f}")
print(f"  Gene d.f.      dg   = {df_gene}")
print(f"  Total d.f.     dg+d0= {df_gene + d0:.2f}  (vs. Welch ~ 1-2)")

# Moderated variance (shrunk toward the prior)
s_tilde_sq = (d0 * s0_sq + df_gene * s2_for_eb) / (d0 + df_gene)

# Moderated t-statistic
se_mod   = np.sqrt(s_tilde_sq * (1.0 / n1 + 1.0 / n2))
t_mod    = delta_psi / se_mod
df_total = df_gene + d0
p_mod    = 2.0 * stats.t.sf(np.abs(t_mod), df_total)   # two-sided

# ══════════════════════════════════════════════════════════════════════════
# 3.  BENJAMINI-HOCHBERG FDR
# ══════════════════════════════════════════════════════════════════════════
_, fdr_mod, _, _ = multipletests(p_mod, method='fdr_bh')

print(f"\nAfter moderated t-test + BH FDR:")
for alpha in (0.05, 0.1, 0.2):
    print(f"  Genes with FDR < {alpha}: {(fdr_mod < alpha).sum()}  "
          f"(both directions combined)")

# ══════════════════════════════════════════════════════════════════════════
# 4.  Z-SCORE OF ΔPSI AGAINST THE POPULATION  (cross-check)
# ══════════════════════════════════════════════════════════════════════════
mad = stats.median_abs_deviation(delta_psi, scale='normal')
med = np.median(delta_psi)
z_dpsi = (delta_psi - med) / mad                       # robust z-score
print(f"\nZ-score (robust) of ΔPSI: median = {med:.4f},  MAD*1.4826 = {mad:.4f}")
for zcut in (1.96, 2.58, 3.0):
    print(f"  |Z| > {zcut}:  {int((np.abs(z_dpsi) > zcut).sum())} genes")

# ══════════════════════════════════════════════════════════════════════════
# 5.  PREPARE PLOTTING
# ══════════════════════════════════════════════════════════════════════════
is_hit = df_all['HITubr3'].values.astype(bool)
hit_idx = np.where(is_hit)[0]

DPSI_CUTOFF = 0.5
FDR_CUTOFF  = 0.05
neg_log_fdr = -np.log10(np.clip(fdr_mod, 1e-300, 1.0))
neg_log_fdr_cutoff = -np.log10(FDR_CUTOFF)
neg_log_fdr = np.clip(neg_log_fdr, 0, 50)

sig_up   = (delta_psi >=  DPSI_CUTOFF) & (fdr_mod < FDR_CUTOFF) & ~is_hit
sig_down = (delta_psi <= -DPSI_CUTOFF) & (fdr_mod < FDR_CUTOFF) & ~is_hit
ns       = ~sig_up & ~sig_down & ~is_hit

n_sig_up   = int(sig_up.sum())
n_sig_down = int(sig_down.sum())

print(f"\nWith cutoffs  |ΔPSI| ≥ {DPSI_CUTOFF}  and  FDR < {FDR_CUTOFF}:")
print(f"  Up-regulated (stabilized):   {n_sig_up}")
print(f"  Down-regulated (destabilized): {n_sig_down}")
print(f"  Labelled best hits (HITubr3): {int(is_hit.sum())}")

# ══════════════════════════════════════════════════════════════════════════
# 6.  STYLING
# ══════════════════════════════════════════════════════════════════════════
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'font.size': 11,
    'axes.titlesize': 14,
    'axes.labelsize': 12,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'svg.fonttype': 'none',
    'axes.spines.top': False,
    'axes.spines.right': False,
})

pos2_colors = {
    'G': '#2ECC71', 'P': '#E74C3C', 'E': '#3498DB', 'T': '#F39C12',
    'D': '#9B59B6', 'A': '#1ABC9C', 'S': '#E67E22', 'Q': '#34495E',
}
DEFAULT_COLOR = '#95A5A6'

pos2_order = ['P', 'G', 'E', 'T', 'Q', 'S', 'A', 'D']
pos2_counts = {}
for idx in hit_idx:
    p2 = df_all.iloc[idx]['AA_seq'][1]
    pos2_counts[p2] = pos2_counts.get(p2, 0) + 1

# ══════════════════════════════════════════════════════════════════════════
# 7.  VOLCANO PLOT — Moderated-t  +  BH FDR,  coloured by Pos-2 residue
# ══════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(17, 7))

ymax = min(np.nanmax(neg_log_fdr), 20) + 1
ax.axvspan(DPSI_CUTOFF,  2.0, ymin=neg_log_fdr_cutoff / ymax,
           alpha=0.03, color='#E67E22', zorder=0)
ax.axvspan(-2.0, -DPSI_CUTOFF, ymin=neg_log_fdr_cutoff / ymax,
           alpha=0.03, color='#2E86C1', zorder=0)

ax.scatter(delta_psi[ns], neg_log_fdr[ns], s=4, c='#BDC3C7', alpha=0.3,
           edgecolors='none', rasterized=True, label='_nolegend_')
ax.scatter(delta_psi[sig_down], neg_log_fdr[sig_down], s=6, c='#BDC3C7',
           alpha=0.3, edgecolors='none', rasterized=True, label='_nolegend_')
ax.scatter(delta_psi[sig_up & ~is_hit], neg_log_fdr[sig_up & ~is_hit], s=6,
           c='#BDC3C7', alpha=0.3, edgecolors='none', rasterized=True,
           label='_nolegend_')

for p2_key in pos2_order:
    idxs = [idx for idx in hit_idx if df_all.iloc[idx]['AA_seq'][1] == p2_key]
    if not idxs:
        continue
    c = pos2_colors.get(p2_key, DEFAULT_COLOR)
    for i, idx in enumerate(idxs):
        lbl = f'{p2_key} (n={len(idxs)})' if i == 0 else None
        ax.scatter(delta_psi[idx], neg_log_fdr[idx], s=50, c=c, alpha=0.9,
                   edgecolors='black', linewidths=0.5, zorder=5, label=lbl)

plotted_pos2 = set(pos2_order)
for idx in hit_idx:
    p2 = df_all.iloc[idx]['AA_seq'][1]
    if p2 in plotted_pos2:
        continue
    c = pos2_colors.get(p2, DEFAULT_COLOR)
    lbl = f'{p2} (n={pos2_counts.get(p2, 1)})'
    ax.scatter(delta_psi[idx], neg_log_fdr[idx], s=50, c=c, alpha=0.9,
               edgecolors='black', linewidths=0.5, zorder=5, label=lbl)
    plotted_pos2.add(p2)

ax.set_ylim(-0.3, ymax)

texts = []
for idx in hit_idx:
    seq  = df_all.iloc[idx]['AA_seq']
    gene = df_all.iloc[idx]['Gene_ID']
    dipeptide = seq[1:3]
    p2 = seq[1]
    c = pos2_colors.get(p2, DEFAULT_COLOR)
    texts.append(ax.text(delta_psi[idx], neg_log_fdr[idx],
                         f'{gene} ({dipeptide})',
                         fontsize=6.5, color=c, fontweight='bold', zorder=6))
adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='#999', lw=0.4),
            force_text=(1.0, 1.5), force_points=(0.5, 0.5), expand=(1.8, 2.2),
            only_move={'points': 'y', 'text': 'xy'})

ax.axvline(x=DPSI_CUTOFF,  color='#333', linestyle='--', linewidth=0.8, alpha=0.6)
ax.axvline(x=-DPSI_CUTOFF, color='#333', linestyle='--', linewidth=0.8, alpha=0.6)
ax.axhline(y=neg_log_fdr_cutoff, color='#333', linestyle='--', linewidth=0.8, alpha=0.6)

ax.text(DPSI_CUTOFF + 0.02, ymax * 0.97, f'ΔPSI = {DPSI_CUTOFF}',
        fontsize=8, color='#555', va='top')
ax.text(-DPSI_CUTOFF - 0.02, ymax * 0.97, f'ΔPSI = -{DPSI_CUTOFF}',
        fontsize=8, color='#555', va='top', ha='right')
ax.text(ax.get_xlim()[1] * 0.98, neg_log_fdr_cutoff + 0.1,
        f'FDR = {FDR_CUTOFF}', fontsize=8, color='#555',
        va='bottom', ha='right')

ax.set_xlabel(r'$\Delta$ PSI (UBR3 $-$ Control)', fontsize=12)
ax.set_ylabel(r'$-\log_{10}$(BH-adjusted p,  moderated t)', fontsize=12)
ax.set_title('Volcano Plot — Moderated t-test (limma-style empirical Bayes) + BH FDR\n'
             f'Best Hits coloured by Position 2 residue '
             f'(n = {len(df_all):,} peptides,  prior d.f. d0 = {d0:.1f})',
             fontweight='bold', fontsize=12)

# Summary stats box
textstr = (f"Moderated t-test (Smyth 2004)\n"
           f"  prior $s_0^2$ = {s0_sq:.3f}\n"
           f"  prior d.f. $d_0$ = {d0:.1f}\n"
           f"  moderated d.f. = $d_g + d_0$ = {df_gene + d0:.1f}\n"
           f"BH FDR:\n"
           f"  FDR<0.05: {(fdr_mod < 0.05).sum():,}  genes\n"
           f"  FDR<0.10: {(fdr_mod < 0.10).sum():,}  genes\n"
           f"Cutoffs: |ΔPSI| ≥ {DPSI_CUTOFF}, FDR < {FDR_CUTOFF}\n"
           f"  → {n_sig_up} up,  {n_sig_down} down")
ax.text(0.015, 0.97, textstr, transform=ax.transAxes, fontsize=8.5,
        va='top', ha='left',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='#FDFEFE',
                  edgecolor='#BDC3C7', alpha=0.95))

ax.legend(fontsize=8.5, frameon=True, loc='upper right', markerscale=1.2,
          edgecolor='#ccc', fancybox=True, ncol=1, title='Pos 2 residue',
          title_fontsize=9)

plt.tight_layout()
base = 'fig20c_volcano_moderated_t'
for ext in ('png', 'pdf', 'svg'):
    path = os.path.join(OUT_DIR, f'{base}.{ext}')
    fig.savefig(path, format=ext)
    print(f"  saved: {path}  ({os.path.getsize(path)/1024:.1f} KB)")
plt.close(fig)

# ══════════════════════════════════════════════════════════════════════════
# 8.  EXPORT RESULTS TABLE
# ══════════════════════════════════════════════════════════════════════════
results = pd.DataFrame({
    'Gene_ID':          df_all.get('Gene_ID', pd.Series([np.nan]*len(df_all))).values,
    'AA_seq':           df_all.get('AA_seq',  pd.Series([np.nan]*len(df_all))).values,
    'AA2':              df_all.get('AA2',     pd.Series([np.nan]*len(df_all))).values,
    'AA3':              df_all.get('AA3',     pd.Series([np.nan]*len(df_all))).values,
    'PSI-293a':         ctrl_a,
    'PSI-293b':         ctrl_b,
    'PSI-UBR3a':        ubr3_a,
    'PSI-UBR3b':        ubr3_b,
    'mean_ctrl':        mean_ctrl,
    'mean_UBR3':        mean_ubr3,
    'delta_PSI':        delta_psi,
    's2_gene':          s2_gene,
    's2_moderated':     s_tilde_sq,
    't_moderated':      t_mod,
    'p_moderated':      p_mod,
    'FDR_BH':           fdr_mod,
    'neg_log10_FDR':    -np.log10(np.clip(fdr_mod, 1e-300, 1.0)),
    'Z_robust_dPSI':    z_dpsi,
    'is_HITubr3':       is_hit,
})
results_sorted = results.sort_values('FDR_BH')

out_xlsx = os.path.join(OUT_DIR, 'volcano_moderated_t_results.xlsx')
with pd.ExcelWriter(out_xlsx, engine='xlsxwriter') as w:
    # README sheet
    readme = pd.DataFrame({
        'Column': [
            'delta_PSI', 's2_gene', 's2_moderated', 't_moderated',
            'p_moderated', 'FDR_BH', 'Z_robust_dPSI', 'is_HITubr3',
            '', 'Prior s0^2', 'Prior d0', 'Moderated d.f.',
        ],
        'Meaning': [
            'mean_UBR3 - mean_ctrl',
            'Pooled per-gene variance (equal-variance, n1=n2=2, d.f.=2)',
            'Empirical-Bayes-shrunk variance (toward prior s0^2)',
            'Moderated t = delta_PSI / sqrt(s2_moderated*(1/n1+1/n2))',
            'Two-sided p-value from t-distribution with (d_g+d_0) d.f.',
            'Benjamini-Hochberg FDR-adjusted p_moderated',
            'Robust z-score of delta_PSI  (x - median) / (1.4826*MAD)',
            'Original HITubr3 flag from input sheet',
            '',
            f'{s0_sq:.4f}',
            f'{d0:.2f}',
            f'd_g + d_0 = {df_gene + d0:.2f}',
        ],
    })
    readme.to_excel(w, sheet_name='README', index=False)
    results_sorted.to_excel(w, sheet_name='All_genes_sorted_by_FDR', index=False)
    results_sorted.query('FDR_BH < 0.05').to_excel(w, sheet_name='FDR_lt_0.05', index=False)
    results_sorted.query('FDR_BH < 0.10').to_excel(w, sheet_name='FDR_lt_0.10', index=False)
    results_sorted.query('is_HITubr3').to_excel(w, sheet_name='Original_best_hits', index=False)

print(f"\nResults table: {out_xlsx}")
print("\nDone.")
