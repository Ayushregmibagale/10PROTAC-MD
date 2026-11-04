#!/usr/bin/env python3
"""
00b_build_merged_manifest.py — build tiered structure manifest for 10PROTAC-MD
================================================================================
Merges two structure sources into a single manifest:

  Source A — KLIFS crystal structures (~/10PROTAC-MD/structures/)
    - Experimental PDB coordinates, curated
    - DFG state annotated (in / out / out-like / inter)
    - 497 downloaded files across 328 kinases (KLIFS naming)
    - Covers ~400+ HGNC kinases after name mapping

  Source B — AlphaFold2 structures (~/kinome_scan/data/structures/)
    - Computationally predicted, full kinome coverage (631 / 641 kinases)
    - Almost always DFG-in (ground-state prediction)
    - Named by UniProt ID (e.g. P00519.pdb)

Tiered selection per PROTAC:
  Foretinib (SJF8240, SJFδ):  KLIFS DFG-out / out-like ONLY
  Lapatinib (SJF1521):         KLIFS DFG-in first → AF2 fallback
  VS-4718   (GSK215):          KLIFS DFG-in first → AF2 fallback
  JH-I-25   (JNJ-1013):        KLIFS DFG-in first → AF2 fallback

Output:
  ~/10PROTAC-MD/structures/merged_manifest.tsv  — one row per HGNC kinase

Usage:
  python 00b_build_merged_manifest.py
"""

import pandas as pd
import re
from pathlib import Path

ROOT        = Path(__file__).resolve().parents[2]
STRUCT_DIR  = ROOT / 'structures'
KLIFS_MAN   = STRUCT_DIR / 'klifs_manifest.tsv'
OUT_TSV     = STRUCT_DIR / 'merged_manifest.tsv'

AF2_DIR     = Path.home() / 'kinome_scan' / 'data' / 'structures'
KIN_LIST    = Path.home() / 'kinome_scan' / 'data' / 'kinase_list.tsv'


# ── 1. Load KLIFS kinase metadata (HGNC names + UniProt IDs) ──────────────────

def load_klifs_kinase_map():
    """Return DataFrame: klifs_name → hgnc_gene, uniprot_id (Human only)."""
    from opencadd.databases.klifs import setup_remote
    session = setup_remote()
    kin = session.kinases.all_kinases()
    human = kin[kin['species.klifs'] == 'Human'].copy()
    human = human[['kinase.klifs_name', 'kinase.gene_name', 'kinase.uniprot']].copy()
    human.columns = ['klifs_name', 'hgnc_gene', 'uniprot_id']
    human = human[human['uniprot_id'] != 0].reset_index(drop=True)
    return human


# ── 2. Scan downloaded KLIFS PDB files and match to manifest ──────────────────

def scan_klifs_files():
    """
    Returns DataFrame: klifs_name, pdb_id, dfg_state, resolution, filepath.
    Parses filename pattern: HUMAN_{GENE}_{PDB}_{...}_complex.pdb
    """
    records = []
    for f in STRUCT_DIR.rglob('*.pdb'):
        m = re.match(r'HUMAN_([^_]+)_([^_]+)_.*complex\.pdb', f.name)
        if not m:
            continue
        klifs_name = m.group(1)
        pdb_id     = m.group(2).lower()
        records.append({'klifs_name': klifs_name, 'pdb_id': pdb_id, 'filepath': str(f)})
    return pd.DataFrame(records)


def match_files_to_manifest(files_df, manifest_df):
    """Join scanned files to manifest on (klifs_name ≈ gene, pdb_id) to get DFG state."""
    # manifest 'gene' column = klifs_name
    merged = files_df.merge(
        manifest_df[['gene', 'pdb_id', 'dfg_state', 'resolution', 'missing_residues', 'quality_score']],
        left_on=['klifs_name', 'pdb_id'],
        right_on=['gene', 'pdb_id'],
        how='left'
    )
    merged = merged.drop(columns=['gene'])
    return merged


def best_klifs_per_kinase_dfg(files_manifest_df):
    """
    For each (klifs_name, dfg_state) keep the single best file:
    prefer quality_score desc, then missing_residues asc, then resolution asc.
    Returns wide DataFrame: one row per klifs_name, columns per DFG state.
    """
    df = files_manifest_df.dropna(subset=['dfg_state']).copy()
    df['quality_score']    = pd.to_numeric(df['quality_score'],    errors='coerce').fillna(0)
    df['missing_residues'] = pd.to_numeric(df['missing_residues'], errors='coerce').fillna(999)
    df['resolution']       = pd.to_numeric(df['resolution'],       errors='coerce').fillna(99)

    df = df.sort_values(['quality_score', 'missing_residues', 'resolution'],
                        ascending=[False, True, True])
    best = df.drop_duplicates(subset=['klifs_name', 'dfg_state'], keep='first')

    rows = {}
    for _, row in best.iterrows():
        kn  = row['klifs_name']
        dfg = row['dfg_state']
        if kn not in rows:
            rows[kn] = {'klifs_name': kn}
        rows[kn][f'klifs_{dfg}_path']  = row['filepath']
        rows[kn][f'klifs_{dfg}_res']   = row['resolution']
        rows[kn][f'klifs_{dfg}_pdb']   = row['pdb_id']

    return pd.DataFrame(list(rows.values()))


# ── 3. Load AF2 structures from kinome_scan ───────────────────────────────────

def load_af2_map():
    """Return DataFrame: uniprot_id, hgnc_gene, af2_path (None if not downloaded)."""
    if not KIN_LIST.exists():
        print(f'  WARN: kinome_scan kinase_list not found at {KIN_LIST}')
        return pd.DataFrame(columns=['uniprot_id', 'hgnc_gene', 'af2_path'])
    kin = pd.read_csv(KIN_LIST, sep='\t')
    records = []
    for _, row in kin.iterrows():
        uid   = str(row['uniprot_id'])
        gene  = str(row['gene'])
        fpath = AF2_DIR / f'{uid}.pdb'
        records.append({
            'uniprot_id': uid,
            'hgnc_gene':  gene,
            'af2_path':   str(fpath) if fpath.exists() else None,
        })
    return pd.DataFrame(records)


# ── 4. Merge and assign per-PROTAC recommended structure ─────────────────────

PROTAC_STRATEGY = {
    'SJF8240':   {'dfg_prefer': ['out', 'out-like'], 'af2_fallback': False},
    'SJF_delta': {'dfg_prefer': ['out', 'out-like'], 'af2_fallback': False},
    'SJF1521':   {'dfg_prefer': ['in'],              'af2_fallback': True},
    'GSK215':    {'dfg_prefer': ['in', 'out', 'out-like', 'inter'], 'af2_fallback': True},
    'JNJ1013':   {'dfg_prefer': ['in'],              'af2_fallback': True},
}


def assign_structure(row, strategy):
    """Pick the best structure path for a given PROTAC strategy."""
    for dfg in strategy['dfg_prefer']:
        col = f'klifs_{dfg}_path'
        if col in row.index and pd.notna(row[col]):
            return row[col], f'klifs_{dfg}'
    if strategy['af2_fallback'] and pd.notna(row.get('af2_path')):
        return row['af2_path'], 'af2'
    return None, None


def build_merged(klifs_kinase_map, klifs_best, af2_map):
    """
    Join on UniProt ID:
      klifs_kinase_map: klifs_name → hgnc_gene + uniprot_id
      klifs_best:       klifs_name → per-DFG paths
      af2_map:          uniprot_id → hgnc_gene + af2_path
    """
    # Start from AF2 map (full kinome), merge in KLIFS via UniProt
    merged = af2_map.merge(
        klifs_kinase_map[['klifs_name', 'uniprot_id']],
        on='uniprot_id', how='left'
    )
    merged = merged.merge(klifs_best, on='klifs_name', how='left')

    # For kinases in KLIFS but not in AF2 list, add them
    klifs_only = klifs_kinase_map[~klifs_kinase_map['uniprot_id'].isin(af2_map['uniprot_id'])].copy()
    klifs_only = klifs_only.merge(klifs_best, on='klifs_name', how='left')
    klifs_only['af2_path'] = None
    klifs_only = klifs_only.rename(columns={'hgnc_gene': 'hgnc_gene'})

    merged = pd.concat([merged, klifs_only], ignore_index=True)
    merged = merged.drop_duplicates(subset=['uniprot_id'], keep='first')

    # Assign recommended structure per PROTAC
    for protac, strategy in PROTAC_STRATEGY.items():
        paths, sources = [], []
        for _, row in merged.iterrows():
            p, s = assign_structure(row, strategy)
            paths.append(p)
            sources.append(s)
        merged[f'{protac}_path']   = paths
        merged[f'{protac}_source'] = sources

    return merged


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    print('Loading KLIFS kinase metadata ...')
    klifs_kinase_map = load_klifs_kinase_map()
    print(f'  {len(klifs_kinase_map)} human kinases with UniProt IDs')

    print('Scanning downloaded KLIFS PDB files ...')
    files_df   = scan_klifs_files()
    manifest   = pd.read_csv(KLIFS_MAN, sep='\t')
    files_man  = match_files_to_manifest(files_df, manifest)
    klifs_best = best_klifs_per_kinase_dfg(files_man)
    n_dfg_out  = sum(1 for c in klifs_best.columns if 'out' in c and 'path' in c
                     for v in [klifs_best[c].notna().sum()] if v > 0)
    print(f'  {len(klifs_best)} kinases with ≥1 KLIFS structure downloaded')
    dfg_out_kin = klifs_best['klifs_out_path'].notna().sum() if 'klifs_out_path' in klifs_best else 0
    dfg_outlike = klifs_best['klifs_out-like_path'].notna().sum() if 'klifs_out-like_path' in klifs_best else 0
    print(f'  Kinases with DFG-out:      {dfg_out_kin}')
    print(f'  Kinases with DFG-out-like: {dfg_outlike}')
    print(f'  Kinases with DFG-in:       {klifs_best["klifs_in_path"].notna().sum() if "klifs_in_path" in klifs_best else 0}')

    print('Loading AF2 structures from kinome_scan ...')
    af2_map = load_af2_map()
    n_af2 = af2_map['af2_path'].notna().sum()
    print(f'  {n_af2} / {len(af2_map)} AF2 structures present on disk')

    print('Building merged manifest ...')
    merged = build_merged(klifs_kinase_map, klifs_best, af2_map)
    print(f'  Total kinases in merged manifest: {len(merged)}')

    # Summary per PROTAC
    print('\nStructure availability per PROTAC:')
    for protac in PROTAC_STRATEGY:
        col = f'{protac}_path'
        src = f'{protac}_source'
        if col not in merged.columns:
            continue
        total    = merged[col].notna().sum()
        n_klifs  = merged[merged[src].str.startswith('klifs', na=False)][col].notna().sum()
        n_af2    = merged[merged[src] == 'af2'][col].notna().sum()
        print(f'  {protac:12s}  {total:4d} kinases  '
              f'(KLIFS: {n_klifs}, AF2: {n_af2})')

    merged.to_csv(OUT_TSV, sep='\t', index=False)
    print(f'\nSaved → {OUT_TSV}')


if __name__ == '__main__':
    main()
