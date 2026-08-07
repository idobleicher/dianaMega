#!/usr/bin/env python
"""Render every figure panel on its own, large, for the slide deck.

Uses the same panel functions as make_figures.py, so a standalone panel and its
counterpart inside the composed figure are guaranteed to be the same plot.
Output: panels/Fig<N><Letter>.png at 200 dpi.
"""
import os
import warnings

import matplotlib.pyplot as plt
import pandas as pd

import ubr3_core as U
import ubr3_panels as P

warnings.filterwarnings('ignore')

OUT = os.path.join(U.HERE, 'panels')


def main():
    os.makedirs(OUT, exist_ok=True)
    lib, hit = U.load()
    enr = pd.read_excel(os.path.join(U.HERE, 'UBR3_PG_substrate_tables.xlsx'),
                        '06_enrich_pos2_residue')
    sig = U.enrichment_test(list(hit.peptide_24mer), list(lib.peptide_24mer))
    # tags off (the slide title carries the panel letter), standalone on
    # (panels that relied on a figure-level legend draw their own)
    C = P.context(lib, hit, sig, enr, tags=False, standalone=True)

    print(f'rendering {len(P.PANELS)} standalone panels ...')
    for key, fn in P.PANELS.items():
        w, h = P.PANEL_SIZE[key]
        fig = plt.figure(figsize=(w, h))
        ax = fig.add_axes([0.10, 0.10, 0.86, 0.80])
        fn(ax, C)
        path = os.path.join(OUT, f'Fig{key}.png')
        fig.savefig(path, dpi=200, bbox_inches='tight', facecolor=P.SURFACE)
        plt.close(fig)
        print(f'  Fig{key}')
    print('done ->', OUT)


if __name__ == '__main__':
    main()
