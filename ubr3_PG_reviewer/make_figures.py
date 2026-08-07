#!/usr/bin/env python
"""Composed multi-panel manuscript figures for the UBR3 N24mer Pro/Gly analysis.

Panels live in ubr3_panels.py and are shared with make_panels.py (the slide deck),
so the two outputs can never drift apart.
Renders 7 figures as 300-dpi PNG + vector PDF into figures/.
Light/print mode only - these are manuscript figures, not screen artifacts.
"""
import os
import warnings

import matplotlib.pyplot as plt

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


def figure9(C):
    grid2x2(C, ['9A', '9B', '9C', '9D'], 'Figure9_stability_x_substrate_classification',
            f'Figure 9 | Motif-bearing peptides classified by baseline stability and substrate status',
            f'Every [P/G]-[E/D] peptide falls into one of four cells: unstable or stable at control '
            f'PSI {C["CUT"]:g}, and a UBR3 substrate or not. 11 of the 16 substrates start unstable; '
            'the remaining 5 were already stable and gained further stability on UBR3 loss. The '
            'unstable-vs-stable difference within the motif class is directional only '
            '(Fisher OR 2.85, p = 0.067) - it is the motif, not the stability class, that is '
            'strongly significant.',
            hspace=0.50)


def figure10(C):
    """Downstream of the motif: substrates vs non-substrates, both motif-bearing."""
    fig = plt.figure(figsize=(16.0, 10.4))
    gs = fig.add_gridspec(2, 3, hspace=0.55, wspace=0.30, height_ratios=[0.85, 1.0],
                          left=0.055, right=0.975, top=0.845, bottom=0.07)
    P.PANELS['10A'](fig.add_subplot(gs[0, :]), C)
    for k, col in zip(['10B', '10C', '10D'], range(3)):
        P.PANELS[k](fig.add_subplot(gs[1, col]), C)
    header(fig, 0.055, 0.968, 0.928,
           'Positions 4-24: what separates motif-bearing substrates from motif-bearing '
           'non-substrates',
           'Both groups carry M-[P/G]-[E/D], so positions 1-3 are identical by construction and only '
           'the downstream window can differ. No single residue or position separates them, but the '
           'substrates carry a markedly more positive net charge across the window.')
    save(fig, 'Figure10_downstream_of_motif')


def main():
    os.makedirs(U.FIGS, exist_ok=True)
    lib, hit = U.load()
    print('rendering figures ...')
    print('  running 480 position x residue Fisher tests ...')
    sig = U.enrichment_test(list(hit.peptide_24mer), list(lib.peptide_24mer))
    C = P.context(lib, hit, sig, tags=True, standalone=False)
    for fn in (figure2, figure3, figure4, figure5, figure6, figure9, figure10):
        fn(C)
    print('done ->', U.FIGS)


if __name__ == '__main__':
    main()
