import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
import os, sys, io, warnings
warnings.filterwarnings('ignore')

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

OUTPUT_DIR = r'c:\Users\User\Desktop\תינוקת\dianaMega\ubr3enrichmentlogo\dipeptide_analysis'
os.makedirs(OUTPUT_DIR, exist_ok=True)

SCREEN_EXCEL = r'c:\Users\User\Downloads\UBR3 Nt screen (1).xlsx'

COMBINED_GROUPS = {
    'PD+PE+PT': ['PD', 'PE', 'PT'],
    'GE+GD':    ['GE', 'GD'],
}

df_all = pd.read_excel(SCREEN_EXCEL, sheet_name='Nprot_5_analyzed')
df_hits = pd.read_excel(SCREEN_EXCEL, sheet_name='sub_high')

df_all['PSI_AAVS'] = (df_all['PSI-293a'] + df_all['PSI-293b']) / 2
df_hits['PSI_AAVS'] = (df_hits['PSI-293a'] + df_hits['PSI-293b']) / 2
df_all['dipeptide'] = df_all['AA2'].astype(str) + df_all['AA3'].astype(str)
df_hits['dipeptide'] = df_hits['AA2'].astype(str) + df_hits['AA3'].astype(str)

dp_data = {}
for grp_name, members in COMBINED_GROUPS.items():
    h = pd.concat([df_hits[df_hits['dipeptide'] == m] for m in members], ignore_index=True)
    a = pd.concat([df_all[df_all['dipeptide'] == m] for m in members], ignore_index=True)
    dp_data[grp_name] = {'hits': h, 'all': a}

print("=" * 80)
print("COMPREHENSIVE ALTERNATIVE TESTS FOR GE+GD AND PD+PE+PT")
print("=" * 80)

for dp in ['PD+PE+PT', 'GE+GD']:
    hits_psi = dp_data[dp]['hits']['PSI_AAVS'].dropna().values
    all_psi  = dp_data[dp]['all']['PSI_AAVS'].dropna().values
    
    print(f"\n{'='*80}")
    print(f"  {dp}")
    print(f"  Top Hits: n={len(hits_psi)}, mean={np.mean(hits_psi):.4f}, sd={np.std(hits_psi, ddof=1):.4f}")
    print(f"  All Screen: n={len(all_psi)}, mean={np.mean(all_psi):.4f}, sd={np.std(all_psi, ddof=1):.4f}")
    print(f"  Mean difference: {np.mean(hits_psi) - np.mean(all_psi):.4f}")
    print(f"{'='*80}")
    
    # 1. Current: Mann-Whitney U (one-tailed)
    mw_stat, mw_p1 = stats.mannwhitneyu(hits_psi, all_psi, alternative='less')
    print(f"\n  1. Mann-Whitney U (one-tailed, H1: hits < all):")
    print(f"     U = {mw_stat:.1f}, p = {mw_p1:.6f}")
    
    # 2. Mann-Whitney U (two-tailed)
    _, mw_p2 = stats.mannwhitneyu(hits_psi, all_psi, alternative='two-sided')
    print(f"\n  2. Mann-Whitney U (two-tailed):")
    print(f"     p = {mw_p2:.6f}")
    
    # 3. Welch's t-test (one-tailed)
    t_stat, t_p2 = stats.ttest_ind(hits_psi, all_psi, equal_var=False)
    t_p1 = t_p2 / 2 if t_stat < 0 else 1 - t_p2 / 2
    print(f"\n  3. Welch's t-test (one-tailed, H1: hits < all):")
    print(f"     t = {t_stat:.4f}, p = {t_p1:.6f}")
    
    # 4. Welch's t-test (two-tailed)
    print(f"\n  4. Welch's t-test (two-tailed):")
    print(f"     t = {t_stat:.4f}, p = {t_p2:.6f}")
    
    # 5. Permutation test (one-tailed)
    rng = np.random.default_rng(42)
    combined = np.concatenate([hits_psi, all_psi])
    n_hits = len(hits_psi)
    observed_diff = np.mean(hits_psi) - np.mean(all_psi)
    n_perm = 100000
    perm_diffs = np.zeros(n_perm)
    for i in range(n_perm):
        perm = rng.permutation(combined)
        perm_diffs[i] = np.mean(perm[:n_hits]) - np.mean(perm[n_hits:])
    perm_p_1t = np.mean(perm_diffs <= observed_diff)
    perm_p_2t = np.mean(np.abs(perm_diffs) >= abs(observed_diff))
    print(f"\n  5. Permutation test (100K permutations):")
    print(f"     Observed diff = {observed_diff:.4f}")
    print(f"     One-tailed p = {perm_p_1t:.6f}")
    print(f"     Two-tailed p = {perm_p_2t:.6f}")
    
    # 6. Bootstrap CI for mean difference
    n_boot = 100000
    boot_diffs = np.zeros(n_boot)
    for i in range(n_boot):
        boot_hits = rng.choice(hits_psi, size=len(hits_psi), replace=True)
        boot_all = rng.choice(all_psi, size=len(all_psi), replace=True)
        boot_diffs[i] = np.mean(boot_hits) - np.mean(boot_all)
    ci_lower = np.percentile(boot_diffs, 2.5)
    ci_upper = np.percentile(boot_diffs, 97.5)
    ci95_lower = np.percentile(boot_diffs, 5)
    print(f"\n  6. Bootstrap CI for mean difference (100K resamples):")
    print(f"     95% CI: [{ci_lower:.4f}, {ci_upper:.4f}]")
    print(f"     90% CI lower bound: {ci95_lower:.4f}")
    print(f"     {'CI excludes 0 -> significant at alpha=0.05' if ci_upper < 0 else 'CI includes 0'}")
    
    # 7. Brunner-Munzel test (better for unequal variances + unequal N)
    bm_stat, bm_p = stats.brunnermunzel(hits_psi, all_psi, alternative='less')
    print(f"\n  7. Brunner-Munzel test (one-tailed, H1: hits stochastically < all):")
    print(f"     W = {bm_stat:.4f}, p = {bm_p:.6f}")
    bm_stat2, bm_p2 = stats.brunnermunzel(hits_psi, all_psi, alternative='two-sided')
    print(f"     Two-tailed: W = {bm_stat2:.4f}, p = {bm_p2:.6f}")
    
    # 8. Kolmogorov-Smirnov (tests distributional difference)
    ks_stat, ks_p = stats.ks_2samp(hits_psi, all_psi, alternative='less')
    print(f"\n  8. Kolmogorov-Smirnov test (one-tailed):")
    print(f"     D = {ks_stat:.4f}, p = {ks_p:.6f}")
    ks_stat2, ks_p2 = stats.ks_2samp(hits_psi, all_psi, alternative='two-sided')
    print(f"     Two-tailed: D = {ks_stat2:.4f}, p = {ks_p2:.6f}")
    
    # 9. Effect sizes
    pooled_sd = np.sqrt((np.var(hits_psi, ddof=1) + np.var(all_psi, ddof=1)) / 2)
    cohen_d = (np.mean(all_psi) - np.mean(hits_psi)) / pooled_sd if pooled_sd > 0 else 0
    glass_d = (np.mean(all_psi) - np.mean(hits_psi)) / np.std(all_psi, ddof=1)
    
    # Hedges' g (bias-corrected Cohen's d for small samples)
    n1, n2 = len(hits_psi), len(all_psi)
    correction = 1 - 3 / (4 * (n1 + n2) - 9)
    hedges_g = cohen_d * correction
    
    # Cliff's delta (non-parametric effect size)
    n_greater = sum(1 for a in all_psi for h in hits_psi if a > h)
    n_less = sum(1 for a in all_psi for h in hits_psi if a < h)
    cliffs_delta = (n_greater - n_less) / (n1 * n2)
    
    # Rank-biserial r (effect size for Mann-Whitney)
    r_rb = 1 - (2 * mw_stat) / (n1 * n2)
    
    print(f"\n  9. Effect sizes:")
    print(f"     Cohen's d = {cohen_d:.4f}")
    print(f"     Hedges' g = {hedges_g:.4f} (bias-corrected)")
    print(f"     Glass's delta = {glass_d:.4f} (using control SD)")
    print(f"     Cliff's delta = {cliffs_delta:.4f}")
    print(f"     Rank-biserial r = {r_rb:.4f}")
    
    # 10. Normality checks
    if len(hits_psi) >= 3:
        sw_hits, sw_p_hits = stats.shapiro(hits_psi)
        print(f"\n  10. Normality (Shapiro-Wilk):")
        print(f"      Hits: W={sw_hits:.4f}, p={sw_p_hits:.4f} {'(normal)' if sw_p_hits > 0.05 else '(non-normal)'}")
    sw_all, sw_p_all = stats.shapiro(all_psi[:min(len(all_psi), 5000)])
    print(f"      All:  W={sw_all:.4f}, p={sw_p_all:.4f} {'(normal)' if sw_p_all > 0.05 else '(non-normal)'}")
    
    print(f"\n  {'='*60}")
    print(f"  SUMMARY for {dp}:")
    print(f"  {'='*60}")
    tests = [
        ("Mann-Whitney U (1-tail)", mw_p1),
        ("Welch t-test (1-tail)", t_p1),
        ("Permutation (1-tail)", perm_p_1t),
        ("Brunner-Munzel (1-tail)", bm_p),
        ("KS test (1-tail)", ks_p),
        ("Mann-Whitney U (2-tail)", mw_p2),
        ("Welch t-test (2-tail)", t_p2),
        ("Permutation (2-tail)", perm_p_2t),
        ("Brunner-Munzel (2-tail)", bm_p2),
        ("KS test (2-tail)", ks_p2),
    ]
    for name, p in sorted(tests, key=lambda x: x[1]):
        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
        marker = " <-- SIGNIFICANT" if p < 0.05 else ""
        print(f"    {name:35s}  p = {p:.6f}  ({sig}){marker}")

print("\n\nDone!")
