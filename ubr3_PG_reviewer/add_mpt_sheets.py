#!/usr/bin/env python
"""Add the M-P-T motif (P at pos2, T at pos3) to UBR3_PG_substrate_tables.xlsx.

The source file (Supplemental Data 1.xlsx) is no longer on disk, so build_workbook.py
cannot be re-run. The MPT peptides are a subset of sheet 03 (every P/G-at-pos2
peptide), so this script takes them from there and writes two new sheets into the
existing workbook:
    27_MPT_motif_all    every MPT peptide in the library, substrates and non-substrates
    28_MPT_enrichment   substrate rate in MPT vs the rest of the library (Fisher)
    29_PGED_MPT_combined  sheets 02 and 27 in one list, with a motif_group column
"""
import os

import pandas as pd
from openpyxl import load_workbook
from openpyxl.formatting.rule import ColorScaleRule, FormulaRule
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
XLSX = os.path.join(HERE, 'UBR3_PG_substrate_tables.xlsx')


def fisher_rows(mpt, rest, n_lib, n_hits, rest_label):
    a, b = int((mpt.UBR3_substrate == 'yes').sum()), len(mpt)
    c, d = int((rest.UBR3_substrate == 'yes').sum()), len(rest)
    odds, p = stats.fisher_exact([[a, b - a], [c, d - c]])
    return {'comparison': f'MPT vs {rest_label}',
            'n_MPT_peptides': b, 'n_MPT_substrates': a,
            'MPT_substrate_rate_pct': round(100 * a / b, 3),
            'n_other_peptides': d, 'n_other_substrates': c,
            'other_substrate_rate_pct': round(100 * c / d, 3),
            'fold_enrichment': round((a / b) / (c / d), 2),
            'odds_ratio': round(odds, 3), 'p_value': float(f'{p:.3g}'),
            'significant_p<0.05': 'yes' if p < 0.05 else 'no'}


def main():
    pg = pd.read_excel(XLSX, '03_PG_library_all')
    s01 = pd.read_excel(XLSX, '01_UBR3_substrates')
    counts = pd.read_excel(XLSX, '05_counts_summary')
    n_lib = int(counts.loc[counts.Metric == 'total N24mer peptides', 'Value'].iloc[0])
    n_hits = len(s01)

    # ---------- sheet 27 : every MPT peptide ----------
    mpt = pg[(pg.pos2 == 'P') & (pg.pos3 == 'T')].drop(columns='rank_by_dPSI').sort_values('mean_dPSI', ascending=False)
    ann = s01.set_index('Gene_ID')[['protein_name', 'function']].to_dict('index')
    mpt['protein_name'] = [ann.get(g, {}).get('protein_name', '') for g in mpt.Gene_ID]
    mpt['function'] = [ann.get(g, {}).get('function', '') for g in mpt.Gene_ID]
    mpt.insert(0, 'rank_by_dPSI', range(1, len(mpt) + 1))

    # ---------- sheet 28 : is MPT enriched among substrates? ----------
    # the non-P/G peptides are not in sheet 03, so "rest of library" is built from counts
    mpt_hits = int((mpt.UBR3_substrate == 'yes').sum())
    rest_lib = pd.DataFrame({'UBR3_substrate': ['yes'] * (n_hits - mpt_hits)
                             + ['no'] * (n_lib - len(mpt) - (n_hits - mpt_hits))})
    rest_pg = pg[~pg.index.isin(mpt.index)]
    rest_p = pg[(pg.pos2 == 'P') & (pg.pos3 != 'T')]
    s28 = pd.DataFrame([fisher_rows(mpt, rest_lib, n_lib, n_hits, 'rest of library'),
                        fisher_rows(mpt, rest_pg, n_lib, n_hits, 'other P/G-at-pos2 peptides'),
                        fisher_rows(mpt, rest_p, n_lib, n_hits, 'other P-at-pos2 peptides')])

    # ---------- sheet 29 : [P/G]-[E/D] and MPT together ----------
    # the two motifs cannot overlap (pos3 is E/D in one, T in the other)
    pged = pd.read_excel(XLSX, '02_PGED_motif_all')
    s29 = pd.concat([pged.assign(motif_group='[P/G]-[E/D]'), mpt.assign(motif_group='MPT')],
                    ignore_index=True).drop(columns='rank_by_dPSI')
    s29 = s29.sort_values('mean_dPSI', ascending=False)
    s29.insert(0, 'rank_by_dPSI', range(1, len(s29) + 1))
    s29.insert(1, 'motif_group', s29.pop('motif_group'))

    wb = load_workbook(XLSX)
    for name in ('27_MPT_motif_all', '28_MPT_enrichment', '29_PGED_MPT_combined'):
        if name in wb.sheetnames:
            del wb[name]
    for name, df in (('27_MPT_motif_all', mpt), ('28_MPT_enrichment', s28),
                     ('29_PGED_MPT_combined', s29)):
        ws = wb.create_sheet(name)
        ws.append(list(df.columns))
        for row in df.itertuples(index=False):
            ws.append([None if pd.isna(v) else v for v in row])
        style(ws, list(df.columns))

    rm = wb['00_README']
    for item, text in (
            ('27_MPT_motif_all', f'All {len(mpt)} library peptides starting M-P-T (P at pos2, T at pos3), '
                                 f'UBR3 substrates and non-substrates together, ranked by mean dPSI. '
                                 f'{mpt_hits} are UBR3 substrates (rows shaded); filter UBR3_substrate.'),
            ('28_MPT_enrichment', 'Two-sided Fisher exact test: UBR3-substrate rate in MPT peptides vs the '
                                  'rest of the library, vs the other P/G-at-pos2 peptides, and vs the other '
                                  'P-at-pos2 peptides.'),
            ('29_PGED_MPT_combined', f'Sheets 02 and 27 in one list: all {len(s29)} [P/G]-[E/D] and MPT '
                                     f'peptides, substrates and non-substrates, ranked by mean dPSI. '
                                     f'motif_group says which motif each peptide carries; '
                                     f'{int((s29.UBR3_substrate == "yes").sum())} are UBR3 substrates '
                                     f'(rows shaded).')):
        hits = [r for r in range(1, rm.max_row + 1) if rm.cell(r, 1).value == item]
        r = hits[0] if hits else rm.max_row + 1
        rm.cell(r, 1, item).font = Font(bold=True, size=10)
        rm.cell(r, 2, text).alignment = Alignment(wrap_text=True, vertical='top')

    wb.save(XLSX)
    print(f'wrote {XLSX}')
    print(f'   27_MPT_motif_all  {len(mpt)} peptides, {mpt_hits} substrates')
    print(s28.to_string(index=False))
    print(f'   29_PGED_MPT_combined  {len(s29)} peptides, '
          f'{int((s29.UBR3_substrate == "yes").sum())} substrates')


def style(ws, cols):
    # same look as build_workbook.style
    hdr_fill = PatternFill('solid', fgColor='1F4E78')
    hdr_font = Font(bold=True, color='FFFFFF', size=10)
    row_fill = PatternFill(fill_type='solid', start_color='FFF8CBAD', end_color='FFF8CBAD')
    row_font = Font(bold=True, color='FF7A2E0E')
    ws.freeze_panes = 'B2'
    ws.auto_filter.ref = ws.dimensions
    for c in range(1, ws.max_column + 1):
        cell = ws.cell(1, c)
        cell.fill, cell.font = hdr_fill, hdr_font
        cell.alignment = Alignment(vertical='center', wrap_text=True, horizontal='center')
        vals = [len(str(ws.cell(r, c).value)) for r in range(2, ws.max_row + 1)
                if ws.cell(r, c).value is not None]
        head = len(str(ws.cell(1, c).value or ''))
        ws.column_dimensions[get_column_letter(c)].width = min(max((max(vals) if vals else 8), min(head, 18)) + 2, 46)
    ws.row_dimensions[1].height = 34
    last = get_column_letter(ws.max_column)
    if 'UBR3_substrate' in cols:
        L = get_column_letter(cols.index('UBR3_substrate') + 1)
        ws.conditional_formatting.add(
            f'A2:{last}{ws.max_row}',
            FormulaRule(formula=[f'${L}2="yes"'], fill=row_fill, font=row_font, stopIfTrue=False))
    for target in ('mean_dPSI', 'fold_enrichment'):
        if target in cols:
            L = get_column_letter(cols.index(target) + 1)
            ws.conditional_formatting.add(
                f'{L}2:{L}{ws.max_row}',
                ColorScaleRule(start_type='min', start_color='FFF5F0', end_type='max', end_color='2A78D6'))


if __name__ == '__main__':
    main()
