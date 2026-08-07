#!/usr/bin/env python
"""Shared loading / annotation layer for the UBR3 N24mer reviewer analysis.

Sheet (A) is the full N24mer GPS library; sheet (B) is the author-defined list of
candidate UBR3 substrates.  Everything downstream (workbook + figures) imports
`load()` from here so the two deliverables can never drift apart.
"""
import datetime
import json
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data')
FIGS = os.path.join(HERE, 'figures')
XLSX = r'C:\Users\User\Downloads\Supplemental Data 1.xlsx'

SHEET_LIB = '(A) N24mer library'
SHEET_HIT = ' (B) UBR3 substrates'

D1, D2 = 'dPSI_rep1', 'dPSI_rep2'
AAS = list('ACDEFGHIKLMNPQRSTVWY')

# PSI cut separating unstable from stable peptides at baseline.
# Data-driven: the control-PSI density has three robust modes (~1.52, ~2.25, ~3.49);
# the antimode between the unstable modes and the stable one sits at 2.59-2.63
# across kernel bandwidths 0.20-0.35, so 2.6 is that boundary rounded.
PSI_CUT = 2.6

# MetAP removes the initiator Met when residue 2 has a small side chain
# (Sherman's rule: radius of gyration <= 1.29 A).
MET_EXCISED = set('ACGPSTV')

# Six-class chemical grouping of the 20 amino acids, used for every
# class-level figure and for sheet 07.
AA_CLASS = {
    'D': 'Acidic', 'E': 'Acidic',
    'K': 'Basic', 'R': 'Basic', 'H': 'Basic',
    'S': 'Polar', 'T': 'Polar', 'N': 'Polar', 'Q': 'Polar', 'C': 'Polar',
    'A': 'Hydrophobic', 'V': 'Hydrophobic', 'L': 'Hydrophobic',
    'I': 'Hydrophobic', 'M': 'Hydrophobic',
    'F': 'Aromatic', 'W': 'Aromatic', 'Y': 'Aromatic',
    'G': 'Special', 'P': 'Special',
}
CLASS_ORDER = ['Acidic', 'Basic', 'Polar', 'Hydrophobic', 'Aromatic', 'Special']
CLASS_MEMBERS = {c: [a for a in AAS if AA_CLASS[a] == c] for c in CLASS_ORDER}

# Kyte-Doolittle hydropathy, for the per-position property track.
KD = {'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5, 'E': -3.5,
      'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8,
      'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2}
CHARGE = {**{a: 0 for a in AAS}, 'D': -1, 'E': -1, 'K': 1, 'R': 1, 'H': 0.1}


def _flatten(cols):
    return ['_'.join(str(x) for x in c if not str(x).startswith('Unnamed')).strip()
            for c in cols]


def _tidy(df):
    """Rename the 30 raw columns into short, stable names."""
    ren = {}
    for c in df.columns:
        n = c
        n = n.replace('Corrected Read Counts: ', '')
        n = n.replace('Control KO  rep #', 'Ctrl_rep').replace('Control KO rep #', 'Ctrl_rep')
        n = n.replace('UBR3 KO  rep #', 'UBR3KO_rep').replace('UBR3 KO rep #', 'UBR3KO_rep')
        n = n.replace('_PSI', '_PSI').replace('_Total Reads', '_TotalReads')
        n = n.replace('\u0394PSI_UBR3 KO_rep#1', D1).replace('\u0394PSI_UBR3 KO_rep#2', D2)
        n = n.replace('Amino Acid Sequence', 'peptide_24mer')
        n = n.replace('Nucleotide Sequence', 'nucleotide_seq')
        ren[c] = n.strip()
    return df.rename(columns=ren)


def fix_gene_symbols(df):
    """Repair gene symbols that Excel turned into dates in the source file.

    MARCH/SEPT family symbols (MARCH5, SEPT2, ...) are auto-converted to dates by
    Excel on entry. 36 rows of sheet (A) are affected. The month encodes the family
    and the day the member, so the symbol is recoverable; where Ensembl resolves the
    transcript we use its current symbol, otherwise we fall back to that rule.
    Column `gene_symbol_repaired` flags every row we touched.
    """
    fixmap = {}
    path = os.path.join(DATA, 'gene_symbol_fix.json')
    if os.path.exists(path):
        fixmap = {k: v for k, v in json.load(open(path, encoding='utf-8')).items() if v}

    def is_date(v):
        return isinstance(v, (datetime.datetime, datetime.date, pd.Timestamp))

    repaired, out = [], []
    for enst, g in zip(df.ENST_ID.astype(str), df.Gene_ID):
        if not is_date(g):
            out.append(g)
            repaired.append('no')
            continue
        sym = fixmap.get(enst)
        if not sym:                       # retired transcript - use the date rule
            fam = {3: 'MARCHF', 9: 'SEPTIN'}.get(g.month)
            sym = f'{fam}{g.day}' if fam else str(g.date())
        out.append(sym)
        repaired.append('yes')
    df['Gene_ID'] = out
    df['gene_symbol_repaired'] = repaired
    return df


def annotate(df):
    """Add N-terminal / motif / property columns to a peptide table."""
    s = df['peptide_24mer'].astype(str).str.strip().str.upper()
    df['peptide_24mer'] = s
    df['pos1'] = s.str[0]
    df['pos2'] = s.str[1]
    df['pos3'] = s.str[2]
    df['pos4'] = s.str[3]

    excised = (df.pos1 == 'M') & df.pos2.isin(MET_EXCISED)
    df['Met_excised_pred'] = np.where(excised, 'yes', 'no')
    df['Nterm_after_MetAP'] = np.where(excised, df.pos2, df.pos1)
    # residue immediately following the exposed N-terminus
    df['Nterm_plus1'] = np.where(excised, df.pos3, df.pos2)

    df['is_PG'] = df.pos2.isin(['P', 'G'])
    df['is_ED_at3'] = df.pos3.isin(['E', 'D'])
    df['is_PG_ED'] = df.is_PG & df.is_ED_at3
    df['motif_class'] = np.select(
        [df.is_PG_ED, df.is_PG & ~df.is_ED_at3],
        ['[P/G]-[E/D]', '[P/G]-other'], default='non-P/G')
    df['Nterm_motif'] = df.pos2 + '-' + df.pos3

    df['mean_dPSI'] = df[[D1, D2]].mean(1)
    df['min_dPSI'] = df[[D1, D2]].min(1)
    df['rep_agreement'] = (df[D1] - df[D2]).abs()

    psi = [c for c in df.columns if c.endswith('_PSI')]
    tot = [c for c in df.columns if c.endswith('_TotalReads')]
    df['mean_PSI_control'] = df[[c for c in psi if c.startswith('Ctrl')]].mean(1)
    df['mean_PSI_UBR3KO'] = df[[c for c in psi if c.startswith('UBR3KO')]].mean(1)
    df['min_total_reads'] = df[tot].min(1)

    # Baseline-stability strata. PSI runs ~1-4 (read-weighted mean FACS bin);
    # PSI = 3 splits the library into unstable and stable halves.
    lo_ctrl = df.mean_PSI_control < PSI_CUT
    lo_ko = df.mean_PSI_UBR3KO < PSI_CUT
    df['PSI_class_control'] = np.where(lo_ctrl, f'PSI < {PSI_CUT}', f'PSI >= {PSI_CUT}')
    df['PSI_class_UBR3KO'] = np.where(lo_ko, f'PSI < {PSI_CUT}', f'PSI >= {PSI_CUT}')
    df['PSI_transition'] = np.select(
        [lo_ctrl & ~lo_ko, ~lo_ctrl & lo_ko, lo_ctrl & lo_ko],
        ['crosses up (unstable -> stable)', 'crosses down (stable -> unstable)',
         'stays unstable'], default='stays stable')
    df['crosses_PSI3_up'] = np.where(lo_ctrl & ~lo_ko, 'yes', 'no')

    df['net_charge_24mer'] = [sum(CHARGE[a] for a in p) for p in s]
    df['hydropathy_GRAVY'] = [round(np.mean([KD[a] for a in p]), 3) for p in s]
    df['acidic_count_24mer'] = [sum(a in 'DE' for a in p) for p in s]
    df['basic_count_24mer'] = [sum(a in 'KRH' for a in p) for p in s]
    return df


def load():
    """Return (library, hits) with hits flagged inside library."""
    lib = pd.read_excel(XLSX, sheet_name=SHEET_LIB, header=[7, 8])
    hit = pd.read_excel(XLSX, sheet_name=SHEET_HIT, header=[0, 1])
    lib.columns = _flatten(lib.columns)
    hit.columns = list(lib.columns)
    lib, hit = _tidy(lib), _tidy(hit)
    lib, hit = fix_gene_symbols(lib), fix_gene_symbols(hit)
    lib, hit = annotate(lib), annotate(hit)

    hitset = set(hit.peptide_24mer)
    lib['is_UBR3_substrate'] = lib.peptide_24mer.isin(hitset)
    lib['UBR3_substrate'] = np.where(lib.is_UBR3_substrate, 'yes', 'no')

    ann = json.load(open(os.path.join(DATA, 'annotation.json'), encoding='utf-8'))
    for col in ['uniprot', 'protein_name', 'function', 'localization',
                'go_biological_process', 'protein_length']:
        hit[col] = [ann.get(str(g).strip(), {}).get(col, '') for g in hit.Gene_ID]
    hit.loc[hit.Gene_ID == 'CU633980.2', ['protein_name', 'function']] = [
        'Uncharacterized protein (clone-based Ensembl gene)',
        'No reviewed UniProt entry; uncharacterized ORF.']
    # a handful of entries carry no FUNCTION comment - fall back to GO, then to
    # an explicit "uncharacterised" note, so no annotation cell is left blank
    blank = hit.function.astype(str).str.strip().isin(['', 'nan'])
    hit.loc[blank, 'function'] = [
        f'No UniProt function comment. GO biological process: {go}.' if str(go).strip()
        else 'Poorly characterised; no UniProt function comment or GO biological process.'
        for go in hit.loc[blank, 'go_biological_process']]

    hit = hit.sort_values('mean_dPSI', ascending=False).reset_index(drop=True)
    hit.insert(0, 'rank', np.arange(1, len(hit) + 1))
    return lib, hit


def position_matrix(seqs, positions=24):
    """counts[position, amino acid] over a list of equal-length peptides."""
    m = pd.DataFrame(0, index=range(1, positions + 1), columns=AAS, dtype=float)
    for s in seqs:
        for i, a in enumerate(s[:positions], 1):
            if a in m.columns:
                m.loc[i, a] += 1
    return m


def freq_matrix(seqs, positions=24):
    m = position_matrix(seqs, positions)
    return m.div(m.sum(1), axis=0)


def class_matrix(seqs, positions=24):
    """Fraction of each chemical class at each position."""
    f = freq_matrix(seqs, positions)
    return pd.DataFrame({c: f[CLASS_MEMBERS[c]].sum(1) for c in CLASS_ORDER})


def enrichment_test(hit_seqs, lib_seqs, groups=None, positions=24):
    """Per position x residue (or x class) Fisher exact, BH-FDR over all cells.

    With only ~50 substrates most cells are pure sampling noise, so every figure
    that shows enrichment masks on `q` rather than plotting the raw ratio.

    Returns (log2_ratio, q_value, signed_neglog10_p) as position-indexed frames.
    """
    from scipy import stats as _st
    from statsmodels.stats.multitest import multipletests

    if groups is None:
        groups = {a: [a] for a in AAS}
    nh, nl = len(hit_seqs), len(lib_seqs)
    ch = position_matrix(hit_seqs, positions)
    cl = position_matrix(lib_seqs, positions)

    names = list(groups)
    lr = pd.DataFrame(0.0, index=ch.index, columns=names)
    pv = pd.DataFrame(1.0, index=ch.index, columns=names)
    for pos in ch.index:
        for g in names:
            a = int(ch.loc[pos, groups[g]].sum())
            c = int(cl.loc[pos, groups[g]].sum())
            _, p = _st.fisher_exact([[a, nh - a], [c, nl - c]])
            pv.loc[pos, g] = p
            # Laplace-smoothed ratio: bounded even when a cell is empty
            fh = (a + 1) / (nh + len(groups[g]))
            fl = (c + 1) / (nl + len(groups[g]))
            lr.loc[pos, g] = np.log2(fh / fl)

    flat = pv.values.ravel()
    q = multipletests(flat, method='fdr_bh')[1].reshape(pv.shape)
    qv = pd.DataFrame(q, index=pv.index, columns=pv.columns)
    signed = -np.log10(pv.clip(lower=1e-300)) * np.sign(lr)
    return lr, qv, signed
