#!/usr/bin/env python
"""Composed multi-panel manuscript figures for the UBR3 N24mer Pro/Gly analysis.

Panels live in ubr3_panels.py and are shared with make_panels.py (the slide deck),
so the two outputs can never drift apart.
Renders 8 figures as 300-dpi PNG + vector PDF into figures/.
Light/print mode only - these are manuscript figures, not screen artifacts.
"""
import os
import warnings

import matplotlib.pyplot as plt
import pandas as pd

import ubr3_core as U
import ubr3_panels as P

warnings.filterwarnings('ignore')

INK, INK2 = P.INK, P.INK2


def save(fig, name):
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(U.FIGS, f'{name}.{ext}'))
    plt.close(fig)
    print('  wrote', name + '.png / .pdf')


def header(fig, x, y_title, y_sub, ttl, sub):
    fig.suptitle(ttl, x=x, y=y_title, ha='left', fontsize=14.5, fontweight='bold', color=INK)
    fig.text(x, y_sub, sub, ha='left', fontsize=9.2, color=INK2)


def grid2x2(C, keys, name, ttl, sub, hspace=0.46):
    fig = plt.figure(figsize=(13.6, 10.6))
    gs = fig.add_gridspec(2, 2, hspace=hspace, wspace=0.26,
                          left=0.07, right=0.975, top=0.855, bottom=0.065)
    for k, cell in zip(keys, [gs[0, 0], gs[0, 1], gs[1, 0], gs[1, 1]]):
        P.PANELS[k](fig.add_subplot(cell), C)
    header(fig, 0.07, 0.968, 0.938, ttl, sub)
    save(fig, name)


def figure1(C):
    grid2x2(C, ['1A', '1B', '1C', '1D'], 'Figure1_PG_landscape',
            'Figure 1 | Pro/Gly N-termini and the [P/G]-[E/D] motif in the N24mer stability screen',
            'GPS screen of 16,514 human N-terminal 24-mers in control-KO vs UBR3-KO cells. '
            'Position 2 becomes the N-terminal residue after initiator-Met excision.')


def figure2(C):
    grid2x2(C, ['2A', '2B', '2C', '2D'], 'Figure2_motif_necessary_not_sufficient',
            'Figure 2 | The [P/G]-[E/D] motif is necessary but not sufficient for UBR3 regulation',
            'Only 16 of the 179 library peptides carrying a Pro/Gly-initiated acidic motif are '
            'stabilised when UBR3 is lost - a 27-fold enrichment over background, but still a '
            'small minority.', hspace=0.50)


def figure3(C):
    fig = plt.figure(figsize=(15.4, 11.9))
    gs = fig.add_gridspec(3, 1, hspace=0.62, left=0.065, right=0.845, top=0.858, bottom=0.055)
    for k, row in zip(['3A', '3B', '3C'], range(3)):
        P.PANELS[k](fig.add_subplot(gs[row]), C)
    from matplotlib.patches import Patch
    handles = [Patch(color=P.CLASS_COLOR[c], label=f'{c}  ({", ".join(U.CLASS_MEMBERS[c])})')
               for c in U.CLASS_ORDER]
    fig.legend(handles=handles, loc='center left', bbox_to_anchor=(0.855, 0.5),
               title='Amino-acid class', title_fontsize=9.5, labelspacing=0.95, alignment='left')
    header(fig, 0.065, 0.968, 0.938,
           'Figure 3 | Sequence logos of the UBR3 substrate N-termini, coloured by amino-acid class',
           'The constraint is concentrated at positions 2 and 3 - Pro/Gly followed by Asp/Glu. '
           'Only one cell further along the peptide (Arg at position 7) survives FDR correction.')
    save(fig, 'Figure3_sequence_logos')


def figure4(C):
    fig = plt.figure(figsize=(14.2, 11.7))
    gs = fig.add_gridspec(3, 2, hspace=0.62, wspace=0.24, height_ratios=[1, 1, 1.05],
                          left=0.07, right=0.98, top=0.872, bottom=0.062)
    P.PANELS['4A'](fig.add_subplot(gs[0, :]), C)
    P.PANELS['4B'](fig.add_subplot(gs[1, :]), C)
    P.PANELS['4C'](fig.add_subplot(gs[2, 0]), C)
    P.PANELS['4D'](fig.add_subplot(gs[2, 1]), C)
    header(fig, 0.07, 0.982, 0.955,
           'Figure 4 | Chemical-class architecture of the substrate N-termini across all 24 positions',
           'Every residue of every peptide is assigned to one of six chemical classes. '
           'Only positions 2, 3 and 7 differ from the library after FDR correction.')
    save(fig, 'Figure4_AA_class_analysis')


def figure5(C):
    fig = plt.figure(figsize=(15.0, 10.1))
    gs = fig.add_gridspec(1, 2, wspace=0.24, left=0.075, right=0.965, top=0.825, bottom=0.09)
    P.PANELS['5A'](fig.add_subplot(gs[0]), C)
    P.PANELS['5B'](fig.add_subplot(gs[1]), C)
    header(fig, 0.075, 0.962, 0.918,
           'Figure 5 | Position-by-residue heatmaps of the 54 candidate UBR3 substrates',
           'The dominant signal is a single block: Pro and Gly at position 2 and Asp at position 3. '
           'The only other cell surviving FDR is Arg at position 7 (22% of substrates vs 7% of the '
           'library).')
    save(fig, 'Figure5_position_residue_heatmaps')


def figure6(C):
    grid2x2(C, ['6A', '6B', '6C', '6D'], 'Figure6_PSI_baseline_stability',
            f'Figure 6 | Baseline stability: peptides above and below PSI {C["CUT"]:g}',
            f'PSI is the read-weighted mean FACS bin (~1-4). The cut at PSI {C["CUT"]:g} is the '
            'antimode of the control-PSI density and splits the library almost evenly '
            f'({100 * C["lo"].mean():.0f}% / {100 * (~C["lo"]).mean():.0f}%). A degron substrate '
            'should start unstable and gain stability when its E3 ligase is lost - which is what '
            'the 54 substrates do.', hspace=0.50)


def figure7(C):
    grid2x2(C, ['7A', '7B', '7C', '7D'], 'Figure7_motif_vs_baseline_stability',
            'Figure 7 | The [P/G]-[E/D] effect is not a by-product of baseline instability',
            'Substrates start unstable, so the motif could in principle be enriched simply because '
            'Pro/Gly N-termini are unstable. Stratifying by control PSI shows it is not: the motif '
            'is evenly distributed across strata yet enriches for substrates within each one.',
            hspace=0.50)


def figure8(C):
    grid2x2(C, ['8A', '8B', '8C', '8D'], 'Figure8_PSI_vs_dPSI_and_cutoff_robustness',
            'Figure 8 | Why baseline stability is measured by control PSI, not by $\\Delta$PSI',
            'Control PSI is how stable a peptide is before UBR3 is removed; $\\Delta$PSI is how '
            'much it changes when UBR3 goes. The 54 substrates were selected on $\\Delta$PSI, so '
            'splitting on it would be circular - every $\\Delta$PSI cut contains all 54 by '
            'construction.', hspace=0.50)


def main():
    os.makedirs(U.FIGS, exist_ok=True)
    lib, hit = U.load()
    enr = pd.read_excel(os.path.join(U.HERE, 'UBR3_PG_substrate_tables.xlsx'),
                        '06_enrich_pos2_residue')
    print('rendering figures ...')
    print('  running 480 position x residue Fisher tests ...')
    sig = U.enrichment_test(list(hit.peptide_24mer), list(lib.peptide_24mer))
    C = P.context(lib, hit, sig, enr, tags=True, standalone=False)
    for fn in (figure1, figure2, figure3, figure4, figure5, figure6, figure7, figure8):
        fn(C)
    print('done ->', U.FIGS)


if __name__ == '__main__':
    main()
