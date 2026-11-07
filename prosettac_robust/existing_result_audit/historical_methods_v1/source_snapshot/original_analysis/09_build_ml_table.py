#!/usr/bin/env python3
"""
Build per-POI ML feature table from all available data sources.

Sources:
  bridging_summary.tsv          — pRosettaC bridging fraction (19 POIs)
  bridging_distance_features.tsv — geometry distance features (19 POIs)
  random_screen/manifest.tsv    — which poses completed MD (8 POIs)
  analysis/random_screen/*/metrics.csv — P_stable per pose

For POIs without MD (CDK17, CDK2/4/6/DAPK1/IRAK1/TLK2), all rs_* features
are set to 0 and has_md=0. P_stable=0 for these is the scientific signal.

Artifact POIs (TNIK/MAPK14 only): prostac_bridging_frac invalid due to kin_min=104 overlap bug.
MAP4K5/SLK: features recomputed from chain-D cluster PDBs (03c/03d Aug 12) — CLEAN non-degraders.
CDK17: confirmed NON-DEGRADER (Aug 12) — was incorrectly labelled DEG.
Clean ML dataset: 17 POIs (5 DEG, 12 NON).

Output: analysis/features/ml_table.tsv (1 row per POI, 19 rows)
"""
import numpy as np
import pandas as pd
from pathlib import Path

ROOT        = Path.home() / 'protac-md'
ANA_ROOT    = ROOT / 'analysis'
RS_ROOT     = ROOT / 'random_screen'
OUT_DIR     = ANA_ROOT / 'features'
OUT_TSV     = OUT_DIR / 'ml_table.tsv'

BRIDGE_TSV   = ANA_ROOT / 'bridging_summary.tsv'
DIST_TSV     = ANA_ROOT / 'bridging_distance_features.tsv'
POCKET_TSV   = ANA_ROOT / 'vhl_pocket_occupancy.tsv'
INTER_TSV    = ANA_ROOT / 'interprotein_near_bridging.tsv'
LYS_CSV      = OUT_DIR / 'lysine_features.csv'
VWH_CSV      = OUT_DIR / 'vwh_poi_level.csv'
MANIFEST     = RS_ROOT / 'manifest.tsv'
RS_ANA_DIR   = ANA_ROOT / 'random_screen'
E2_POSE_TSV  = ANA_ROOT / 'e2_reachability_all_pois.tsv'      # per-pose (has 50A columns)
E2_SUM_TSV   = ANA_ROOT / 'e2_reachability_poi_summary.tsv'
LYS_DYN_TSV  = ANA_ROOT / 'lysine_sasa_poi_summary.tsv'       # from 17_all_poi_lysine_sasa
VKI_SUM_TSV  = ANA_ROOT / 'vhl_kin_interface_poi_summary.tsv' # from 18_vhl_kin_interface
PER_LYS_TSV  = ANA_ROOT / 'per_lysine_e2_summary.tsv'         # from 20_per_lysine_e2_distance

ARTIFACT_POIS = {'TNIK', 'MAPK14'}
# TNIK, MAPK14: pRosettaC artifact (kin_min=104 chain-A overlap); 0% from chain-D but excluded as science finding.
# MAP4K5, SLK: chain-D reruns COMPLETE (Aug 12); all features (03c/03d) recomputed from chain-D cluster PDBs.
#   Confirmed geometry-blocked non-degraders: bridge_frac=0%, CoM=68.7Å/67.2Å, n_near_vhl=1/0.

BRIDGE_THRESHOLD = 5.0   # Å

DEGRADER_POIS = {
    'DDR2', 'MET', 'RIPK2_full', 'MAPK14_full', 'EPHB2',  # high-affinity degraders
    'TNIK',                                                  # low-affinity off-target (artifact POI, excluded from ML)
    # CDK17 is a NON-DEGRADER (confirmed experimentally Aug 12). Previously mislabelled.
    # MAPK14 removed: same protein as MAPK14_full; stays in ARTIFACT_POIS.
}


# ── P_stable from random screen ───────────────────────────────────────────────

def _parse_rs_metrics(csv_path):
    rows = []
    try:
        for line in open(csv_path):
            if line.startswith('#') or line.startswith('time'):
                continue
            parts = line.strip().split(',')
            if len(parts) >= 3:
                try:
                    rows.append([float(parts[0]), float(parts[1]), float(parts[2])])
                except ValueError:
                    continue
    except Exception:
        return None
    if not rows:
        return None
    arr = np.array(rows)
    t, d_vhl, d_kin = arr[:, 0], arr[:, 1], arr[:, 2]
    bridging = (d_vhl < BRIDGE_THRESHOLD) & (d_kin < BRIDGE_THRESHOLD)
    bpct = float(bridging.mean() * 100)

    # Bridging events
    events, in_b, s = [], False, 0
    dt = float(t[1] - t[0]) if len(t) > 1 else 0.1
    for i, b in enumerate(bridging):
        if b and not in_b:  s = i; in_b = True
        elif not b and in_b: events.append((i - s) * dt); in_b = False
    if in_b: events.append((len(bridging) - s) * dt)

    return {
        'bridging_pct':    bpct,
        'n_events':        len(events),
        'mean_dur_ns':     float(np.mean(events)) if events else 0.0,
        'max_dur_ns':      float(np.max(events))  if events else 0.0,
        'vhl_pct':         float((d_vhl < BRIDGE_THRESHOLD).mean() * 100),
        'kin_pct':         float((d_kin < BRIDGE_THRESHOLD).mean() * 100),
        'traj_len_ns':     float(t[-1]),
    }


def extract_p_stable(poi):
    """Return dict of rs_* features for a POI; zeros if no MD."""
    empty = {
        'rs_bridging_mean': 0.0, 'rs_bridging_std': 0.0,
        'rs_bridging_max': 0.0,  'rs_bridging_min': 0.0,
        'rs_n_poses': 0, 'rs_frac_poses_gt50': 0.0, 'rs_frac_poses_zero': 1.0,
        'rs_mean_event_dur_ns': 0.0, 'rs_max_event_dur_ns': 0.0,
        'rs_mean_n_events_per_ns': 0.0,
        'rs_vhl_engagement_pct': 0.0, 'rs_kin_engagement_pct': 0.0,
        'has_md': 0,
    }

    if not MANIFEST.exists():
        return empty

    manifest = pd.read_csv(MANIFEST, sep='\t')
    done = manifest[(manifest['poi'] == poi) & (manifest['stage'] == 'ANALYSIS_DONE')]
    if done.empty:
        return empty

    bpcts, mean_durs, max_durs, nevts, vhl_pcts, kin_pcts = [], [], [], [], [], []

    for _, row in done.iterrows():
        csv = RS_ANA_DIR / poi / row['pose_id'] / 'metrics.csv'
        if not csv.exists():
            continue
        r = _parse_rs_metrics(csv)
        if r is None:
            continue
        bpcts.append(r['bridging_pct'])
        mean_durs.append(r['mean_dur_ns'])
        max_durs.append(r['max_dur_ns'])
        nevts.append(r['n_events'] / max(r['traj_len_ns'], 0.01))
        vhl_pcts.append(r['vhl_pct'])
        kin_pcts.append(r['kin_pct'])

    if not bpcts:
        return empty

    pcts = np.array(bpcts)
    return {
        'rs_bridging_mean':        float(pcts.mean()),
        'rs_bridging_std':         float(pcts.std()),
        'rs_bridging_max':         float(pcts.max()),
        'rs_bridging_min':         float(pcts.min()),
        'rs_n_poses':              len(pcts),
        'rs_frac_poses_gt50':      float((pcts > 50).mean()),
        'rs_frac_poses_zero':      float((pcts == 0).mean()),
        'rs_mean_event_dur_ns':    float(np.mean(mean_durs)),
        'rs_max_event_dur_ns':     float(np.mean(max_durs)),
        'rs_mean_n_events_per_ns': float(np.mean(nevts)),
        'rs_vhl_engagement_pct':   float(np.mean(vhl_pcts)),
        'rs_kin_engagement_pct':   float(np.mean(kin_pcts)),
        'has_md': 1,
    }


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # ── Load base tables ──
    bridge_df = pd.read_csv(BRIDGE_TSV, sep='\t')
    dist_df   = pd.read_csv(DIST_TSV,   sep='\t')

    # Merge on poi
    df = bridge_df.merge(dist_df.drop(columns=['category', 'n_poses'], errors='ignore'),
                         on='poi', how='left')

    # Labels
    df['experimental_degrader'] = df['poi'].isin(DEGRADER_POIS).astype(int)
    df['prostac_artifact']      = df['poi'].isin(ARTIFACT_POIS).astype(int)

    # Rename bridging fraction
    if 'bridge_pct' in df.columns:
        df['prostac_bridging_frac'] = df['bridge_pct'] / 100.0
        df = df.drop(columns=['bridge_pct'])

    # ── Add P_stable features ──
    print("\nExtracting P_stable from random screen...\n")
    rs_rows = []
    for poi in df['poi']:
        r = extract_p_stable(poi)
        r['poi'] = poi
        rs_rows.append(r)
        print(f"  [{poi:<16}] has_md={r['has_md']}  "
              f"rs_bridging_mean={r['rs_bridging_mean']:.2f}%  "
              f"n_poses={r['rs_n_poses']}")

    rs_df = pd.DataFrame(rs_rows)
    df = df.merge(rs_df, on='poi', how='left')

    # ── Add VHL pocket and inter-protein geometry features ──
    pocket_df = pd.read_csv(POCKET_TSV, sep='\t')[['poi', 'n_near_vhl', 'pct_pocket_engaged']]
    inter_df  = pd.read_csv(INTER_TSV,  sep='\t')[['poi', 'ptc_mid_dist_mean', 'com_dist_mean', 'min_ca_dist_mean']]
    df = df.merge(pocket_df, on='poi', how='left')
    df = df.merge(inter_df,  on='poi', how='left')

    # ── Add lysine features (per-POI best rank) ──
    if LYS_CSV.exists():
        lys_df = pd.read_csv(LYS_CSV)
        # Aggregate per POI: use rank1 if available, else mean
        lys_agg = (lys_df.sort_values('rank')
                         .groupby('poi')[['lys_nearest_e2_dist_A', 'lys_n_exposed_static',
                                         'lys_delta_sasa_A2', 'lys_e2_reach_prob_80A',
                                         'lys_n_reachable_80A', 'lys_sasa_static_mean_A2']]
                         .first()
                         .reset_index())
        lys_agg.columns = ['poi'] + [f'lk_{c}' for c in lys_agg.columns[1:]]
        df = df.merge(lys_agg, on='poi', how='left')
        print(f"Lysine features merged for {lys_agg['poi'].nunique()} POIs")

    # ── Add E2 reachability dynamic features (from 16_all_poi_e2_reachability) ──
    if E2_POSE_TSV.exists():
        e2p = pd.read_csv(E2_POSE_TSV, sep='\t')
        # Compute per-POI means including delta_50A (not in summary TSV)
        e2_rows = []
        for poi_name, grp in e2p.groupby('poi'):
            # Only average poses that have bridging data (reach_50A_bridge not NaN)
            has_b = grp['reach_50A_bridge'].notna()
            d50 = (grp.loc[has_b, 'reach_50A_bridge'] - grp.loc[has_b, 'reach_50A_nonbridge']).mean() \
                  if has_b.any() else float('nan')
            d80 = (grp.loc[has_b, 'reach_80A_bridge'] - grp.loc[has_b, 'reach_80A_nonbridge']).mean() \
                  if has_b.any() else float('nan')
            e2_rows.append({
                'poi':                poi_name,
                'e2_delta_reach_50A': d50,
                'e2_delta_reach_80A': d80,
                'e2_mean_min_dist_A': float(grp['min_e2_dist_mean'].mean()),
                'e2_reach_80A_bridge':float(grp.loc[has_b, 'reach_80A_bridge'].mean())
                                      if has_b.any() else float('nan'),
            })
        e2_df = pd.DataFrame(e2_rows)
        df = df.merge(e2_df, on='poi', how='left')
        n_e2 = e2_df['poi'].nunique()
        print(f"E2 reachability features merged for {n_e2} POIs")

        # Δ50A vs Δ80A comparison table
        e2s = pd.read_csv(E2_SUM_TSV, sep='\t') if E2_SUM_TSV.exists() else e2_df
        print(f"\n  Δ50A vs Δ80A discriminability:")
        print(f"  {'POI':<15} {'Δ50A':>8} {'Δ80A':>8} {'Label'}")
        for _, r in e2_df.sort_values('e2_delta_reach_50A', ascending=False).iterrows():
            lbl = 'DEG' if r['poi'] in DEGRADER_POIS else 'non'
            print(f"  {r['poi']:<15} {r['e2_delta_reach_50A']:>+8.4f} {r['e2_delta_reach_80A']:>+8.4f}  {lbl}")
        print()
    else:
        print(f"[SKIP] {E2_POSE_TSV.name} not found — run 16_all_poi_e2_reachability.py")

    # ── Add dynamic lysine burial features (from 17_all_poi_lysine_sasa) ──
    if LYS_DYN_TSV.exists():
        ld = pd.read_csv(LYS_DYN_TSV, sep='\t')[
            ['poi', 'mean_burial_bridge', 'mean_burial_nonbridge',
             'delta_burial', 'frac_buried_bridge', 'delta_frac_buried']]
        ld.columns = ['poi'] + [f'lys_dyn_{c}' for c in ld.columns[1:]]
        df = df.merge(ld, on='poi', how='left')
        print(f"Dynamic lysine burial features merged for {ld['poi'].nunique()} POIs")
    else:
        print(f"[SKIP] {LYS_DYN_TSV.name} not found — run 17_all_poi_lysine_sasa.py")

    # ── Add VHL-kinase interface contact features (from 18_vhl_kin_interface) ──
    if VKI_SUM_TSV.exists():
        vki = pd.read_csv(VKI_SUM_TSV, sep='\t')[
            ['poi', 'contacts_bridge', 'contacts_nonbridge',
             'delta_contacts', 'com_dist_bridge', 'delta_com_dist']]
        vki.columns = ['poi'] + [f'vki_{c}' for c in vki.columns[1:]]
        df = df.merge(vki, on='poi', how='left')
        print(f"VHL-kinase interface features merged for {vki['poi'].nunique()} POIs")
    else:
        print(f"[SKIP] {VKI_SUM_TSV.name} not found — run 18_vhl_kin_interface.py")

    # ── Add per-lysine E2 features (from 20_per_lysine_e2_distance) ──
    if PER_LYS_TSV.exists():
        plys = pd.read_csv(PER_LYS_TSV, sep='\t')
        plys_agg = plys.groupby('poi').agg(
            per_lys_min_delta_e2  = ('delta_e2',        'min'),   # most negative = bridging helps most
            per_lys_best_e2_bridge= ('mean_e2_bridge',  'min'),   # nearest lysine when bridging
            per_lys_best_reach_50A= ('reach_50A_bridge','max'),   # best 50A reach across lysines
        ).reset_index()
        df = df.merge(plys_agg, on='poi', how='left')
        print(f"Per-lysine E2 features merged for {plys_agg['poi'].nunique()} POIs")
    else:
        print(f"[SKIP] {PER_LYS_TSV.name} not found — run 20_per_lysine_e2_distance.py")

    # ── Add VHL warhead burial features ──
    if VWH_CSV.exists():
        vwh_df = pd.read_csv(VWH_CSV)
        # Normalise names to match ml_table (MAPK14→MAPK14_full, RIPK2→RIPK2_full)
        vwh_df['poi'] = vwh_df['poi'].replace({'MAPK14': 'MAPK14_full', 'RIPK2': 'RIPK2_full'})
        vwh_cols = ['poi', 'vwh_dist_bridge', 'vwh_dist_nonbridge', 'vwh_delta',
                    'vwh_dist_overall', 'linker_ext_bridge', 'linker_ext_overall']
        vwh_df = vwh_df[[c for c in vwh_cols if c in vwh_df.columns]]
        df = df.merge(vwh_df, on='poi', how='left')
        n_merged = df['vwh_dist_bridge'].notna().sum()
        print(f"VWH warhead features merged for {n_merged} POIs")

    # ── Column order ──
    id_cols   = ['poi', 'category', 'experimental_degrader', 'prostac_artifact', 'has_md']
    prc_cols  = ['prostac_bridging_frac', 'n_bridging', 'n_total', 'best_vhl_A', 'best_kin_A']
    geo_cols  = [c for c in df.columns if c.startswith('d_') or c.startswith('bridge_frac')
                 or c in ('kin_when_vhl_near', 'vhl_when_kin_near', 'r_dvhl_dkin', 'r_pval',
                           'n_vhl_near', 'n_kin_near')]
    rs_cols   = [c for c in df.columns if c.startswith('rs_')]
    keep = id_cols + [c for c in prc_cols if c in df.columns] + \
           [c for c in geo_cols if c in df.columns] + \
           [c for c in rs_cols if c in df.columns]
    extra = [c for c in df.columns if c not in keep]
    df = df[keep + extra]

    df.to_csv(OUT_TSV, sep='\t', index=False, float_format='%.4f')
    print(f"\nSaved: {OUT_TSV}")
    print(f"Shape: {df.shape[0]} rows × {df.shape[1]} columns")

    # Quick preview
    print(f"\n{'POI':<16} {'DEG':>4} {'ART':>4} {'MD':>3} "
          f"{'prc_frac%':>10} {'d_sum_min':>10} {'rs_bridge%':>10}")
    print('-' * 60)
    for _, r in df.sort_values('experimental_degrader', ascending=False).iterrows():
        pf = r.get('prostac_bridging_frac', float('nan'))
        ds = r.get('d_sum_min', float('nan'))
        rs = r.get('rs_bridging_mean', 0.0)
        print(f"{r['poi']:<16} {r['experimental_degrader']:>4} {r['prostac_artifact']:>4} "
              f"{r['has_md']:>3} {pf*100:>10.3f} {ds:>10.2f} {rs:>10.2f}")


if __name__ == '__main__':
    main()
