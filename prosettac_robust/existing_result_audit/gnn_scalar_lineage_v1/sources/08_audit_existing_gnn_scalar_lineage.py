#!/usr/bin/env python3
"""Reconcile current retained GNN scalars to declared source tables, read-only.

Arithmetic agreement is not proof of the exact historical training input: no
historical immutable training receipt is inferred from matching current files.
"""
import argparse,json,sys
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from protac_qc import ValidationError,sha256,atomic_json
COLS=['e2_mean_min_dist_A','e2_delta_reach_50A','e2_delta_reach_80A','lk_lys_nearest_e2_dist_A']


def expected_scalars(e2,lys):
    if e2.duplicated(['poi','pose_id']).any():raise ValidationError('Duplicate E2 pose records')
    records={};denominators={}
    for poi,group in e2.groupby('poi'):
        has_bridge=group['reach_50A_bridge'].notna()
        values={'e2_mean_min_dist_A':float(group['min_e2_dist_mean'].mean())}
        counts={'e2_mean_min_dist_A':int(group['min_e2_dist_mean'].notna().sum())}
        for threshold in [50,80]:
            diff=group.loc[has_bridge,f'reach_{threshold}A_bridge']-group.loc[has_bridge,f'reach_{threshold}A_nonbridge']
            values[f'e2_delta_reach_{threshold}A']=float(diff.mean())
            counts[f'e2_delta_reach_{threshold}A']=int(diff.notna().sum())
        records[poi]=values;denominators[poi]=counts
    for poi,group in lys.sort_values('rank').groupby('poi'):
        valid=group[group['lys_nearest_e2_dist_A'].notna()]
        value=float(valid.iloc[0]['lys_nearest_e2_dist_A']) if len(valid) else float('nan')
        records.setdefault(poi,{})['lk_lys_nearest_e2_dist_A']=value
        denominators.setdefault(poi,{})['lk_lys_nearest_e2_dist_A']=1 if len(valid) else 0
        denominators[poi]['lk_selected_rank']=float(valid.iloc[0]['rank']) if len(valid) else None
    return records,denominators


def audit(project,output,cohort_file):
    project=project.resolve();output=output.resolve()
    if output.exists() or output.is_relative_to(project) or output.is_relative_to(Path('/mnt/d/PROTAC degradation Prediction')):raise ValidationError('Fresh output outside original projects required')
    files={'feature_table':project/'analysis/features/ml_table.tsv','E2_per_pose':project/'analysis/e2_reachability_all_pois.tsv','static_lysine':project/'analysis/features/lysine_features.csv','feature_builder':project/'scripts/02_analysis/09_build_ml_table.py','GNN_runner':project/'scripts/02_analysis/04g_e2_scalar_gnn.py','E2_producer':project/'scripts/02_analysis/16_all_poi_e2_reachability.py','cohort':cohort_file}
    hashes={name:sha256(path) for name,path in files.items()};source_hash=sha256(Path(__file__))
    table=pd.read_csv(files['feature_table'],sep='\t');e2=pd.read_csv(files['E2_per_pose'],sep='\t');lys=pd.read_csv(files['static_lysine'])
    if table['poi'].duplicated().any():raise ValidationError('Duplicate targets in retained feature table')
    labels=json.loads(cohort_file.read_text())['protein_groups'];expected,counts=expected_scalars(e2,lys);rows=[];missing=[]
    table=table.set_index('poi')
    for poi in labels:
        if poi not in table.index:missing.append(poi);continue
        for column in COLS:
            actual=float(table.loc[poi,column]);wanted=expected.get(poi,{}).get(column,float('nan'))
            both_missing=bool(np.isnan(actual) and np.isnan(wanted))
            if not both_missing and (np.isinf(actual) or np.isinf(wanted)):raise ValidationError('Infinite scalar')
            delta=abs(actual-wanted) if np.isfinite(actual) and np.isfinite(wanted) else None
            rows.append({'poi':poi,'column':column,'saved_value':None if np.isnan(actual) else actual,'reconstructed_value':None if np.isnan(wanted) else wanted,'absolute_difference':delta,'matches_stored_precision':bool(both_missing or (delta is not None and delta<=5.01e-5)),'source_observations':counts.get(poi,{}).get(column,0),'selected_rank':counts.get(poi,{}).get('lk_selected_rank') if column.startswith('lk_') else None})
    if any(sha256(path)!=hashes[name] for name,path in files.items()) or sha256(Path(__file__))!=source_hash:raise ValidationError('Scalar audit source/input changed during execution')
    output.mkdir(parents=True);pd.DataFrame(rows).to_csv(output/'scalar_reconciliation.tsv',sep='\t',index=False,float_format='%.12g')
    mismatches=[r for r in rows if not r['matches_stored_precision']]
    report={'status':'CURRENT_RETAINED_SCALARS_RECONCILE' if not missing and not mismatches else 'SCALAR_LINEAGE_REVIEW_REQUIRED','cohort_targets':len(labels),'checked_scalars':len(rows),'matching_scalars':sum(r['matches_stored_precision'] for r in rows),'missing_targets':missing,'mismatches':mismatches,'source_file_sha256':hashes,'source_paths':{k:str(v) for k,v in files.items()},'audit_source_sha256':source_hash,'historical_exact_training_input_proven':False,'original_files_modified':False,'GNN_retrained':False,'new_MD_or_docking_executed':False,'actual_scalar_definitions':{'e2_mean_min_dist_A':'Equal-pose mean of each pose mean nearest E2-ensemble-point/kinase-Lys-NZ distance over producer frames; not conditioned on bridging','e2_delta_reach_50A':'Mean per-pose bridge-minus-nonbridge geometric ensemble-reach fraction; only poses with both finite conditional values contribute','e2_delta_reach_80A':'As for 50 Å, using 80 Å; inherited bridging mask refers to legacy Elongin-B-labelled proximity','lk_lys_nearest_e2_dist_A':'First nonmissing static-lysine value after sorting rank within POI; not nearest across all MD frames and no mean fallback'},'limitations':['Matching current table arithmetic does not establish immutable provenance for the historical training invocation.','Old bridging-conditioned descriptors use the historical scaffold-proximity mask.','Modeled E2 ensemble proximity to any kinase lysine is a geometric descriptor, not a catalytic or degradation probability.','Source producer truncates to minimum frame counts without timestamp verification; audit that separately before reusing the algorithm.'],'output_sha256':sha256(output/'scalar_reconciliation.tsv')}
    atomic_json(output/'summary.json',report);print(json.dumps({k:v for k,v in report.items() if k in ['status','cohort_targets','checked_scalars','matching_scalars','missing_targets','mismatches']},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['project','output','cohort-file']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();audit(a.project,a.output,a.cohort_file)
