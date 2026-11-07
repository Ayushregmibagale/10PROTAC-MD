"""Reconcile retained 19-target configurations, structures and GNN artifacts.

Never execute original modules or deserialize the old pickle graph packets.
"""
from pathlib import Path
import ast,csv,collections,json,sys,re
import numpy as np
from Bio.SeqUtils import seq1
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'workflow/scripts/lib'))
from protac_qc import ValidationError,sha256,atomic_json,pdb_atoms
from protein_roles import identify
PROJECT=Path('/home/ayushregmibagale/protac-md')


def config(path,name):
    tree=ast.parse(path.read_text());node=next(n.value for n in tree.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id==name for t in n.targets));result={}
    for key,value in zip(node.keys,node.values):
        entry={}
        for k,v in zip(value.keys,value.values):
            field=ast.literal_eval(k)
            try:entry[field]=ast.literal_eval(v)
            except ValueError:
                if isinstance(v,ast.Call) and isinstance(v.func,ast.Name) and v.func.id=='Path' and len(v.args)==1 and not v.keywords:entry[field]=ast.literal_eval(v.args[0])
                else:raise ValidationError('Nonliteral source config field '+field)
        result[ast.literal_eval(key)]=entry
    return result


def main(output):
    if output.exists():raise ValidationError('Fresh review output required')
    paths=[PROJECT/'scripts/02_analysis/01_screen_bridging_all_pois.py',PROJECT/'scripts/02_analysis/04a_build_ensemble_graphs.py',PROJECT/'scripts/02_analysis/04c_build_md_contact_graphs.py',PROJECT/'scripts/02_analysis/04d_dual_encoder_gnn.py',PROJECT/'scripts/02_analysis/04g_e2_scalar_gnn.py',PROJECT/'analysis/gnn/e2gnn_lopo_results.tsv',PROJECT/'analysis/gnn/kinome_predictions.tsv',PROJECT/'analysis/gnn/predict_kinome.log',PROJECT/'analysis/gnn/ensemble_graphs.pt',PROJECT/'analysis/gnn/md_contact_graphs.pt',PROJECT/'analysis/gnn/ensemble_graphs_np.npz',PROJECT/'analysis/gnn/md_contact_graphs_np.npz',ROOT/'sjf8240_review/current_graph_packet_inventory.json',ROOT/'sjf8240_review/lineage_summary.json',ROOT/'existing_result_audit/scientific_signoff_v1/protein_groups_19.json',ROOT/'existing_result_audit/scientific_signoff_v1/all19_sequences.json',ROOT/'existing_result_audit/scientific_signoff_v1/coverage_19_targets.json']
    hashes={str(p):sha256(p) for p in paths};static=config(paths[1],'POIS_CONFIG');dynamic=config(paths[2],'POI_CONFIG');screen=config(paths[0],'POIS')
    groups=json.loads(paths[14].read_text());entries={r['primaryAccession']:r['sequence']['value'] for r in json.loads(paths[15].read_text())['results']};coverage={r['poi']:r for r in json.loads(paths[16].read_text())};lineage=json.loads(paths[13].read_text())
    inventory=json.loads(paths[12].read_text())
    for p in paths[8:10]:
        if lineage['artifact_sha256']['analysis/gnn/'+p.name]['sha256']!=hashes[str(p)]:raise ValidationError('Current graph packet changed since safe metadata inventory')
    packets={r['file']:{g['poi']:g for g in r['graphs']} for r in inventory}
    with paths[5].open() as f:pred={r['poi']:r for r in csv.DictReader(f,delimiter='\t')}
    with paths[6].open() as f:future={r['poi']:r for r in csv.DictReader(f,delimiter='\t')}
    npz_targets={}
    for p in paths[10:12]:
        with np.load(p,allow_pickle=False) as a:npz_targets[p.name]=sorted({k.split('/')[0] for k in a.files})
    log_text=paths[7].read_text();log_md={key:value=='True' for key,value in re.findall(r'\[(SIK3_full|RPS6KB2_full)\s*\].*?has_md=(True|False)',log_text)}
    log_matches=bool(log_md) and all(log_md.get(key)==(row['has_md'].lower()=='true') for key,row in future.items())
    gnn_tree=ast.parse(paths[4].read_text());clean=next(n.value for n in ast.walk(gnn_tree) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='clean' for t in n.targets));clean_source=ast.unparse(clean)
    records=[]
    for poi,uid in sorted(groups.items()):
        key=poi+'_full' if poi in {'SIK3','RPS6KB2'} else poi;cfg=static[key];graph=packets['ensemble_graphs.pt'].get(key)
        if graph is None:raise ValidationError('Confirmed target lacks current static graph: '+key)
        expected=pred[poi]['label'] if poi in pred else '-1'
        if int(graph['label'])!=int(expected):raise ValidationError('Graph and retained label differ: '+key)
        directory=Path(cfg['path']) if 'path' in cfg else Path('/mnt/d/PROTAC degradation Prediction/pRosettaC/runs')/cfg['base']/poi/cfg.get('src_subdir','Patchdock_Results')
        pose=next(directory.rglob('combined_*.pdb'),None)
        if pose is None:raise ValidationError('No retained representative for '+key)
        before=sha256(pose);atoms=pdb_atoms(pose);ca=[a for a in atoms if a['record']=='ATOM' and a['name']=='CA'];old_vhl=[a for a in ca if a['chain']=='A' and 1<=a['resnum']<=103]
        target=[a for a in ca if (a['chain']=='D' if cfg['kin_chain']=='D' else a['chain']==cfg['kin_chain'] and a['resnum']>=cfg['kin_min'])]
        findings={}
        for name,selected,candidates in [('legacy_VHL',old_vhl,{k:entries[k] for k in ['P40337','Q15370','Q15369']}),('POI',target,{uid:entries[uid]})]:
            try:findings[name]=identify(''.join(seq1(a['resname'],undef_code='X') for a in selected),candidates)
            except ValidationError as exc:findings[name]={'status':'SEQUENCE_IDENTITY_REVIEW_REQUIRED','reason':str(exc),'CA_atoms':len(selected)}
        if sha256(pose)!=before:raise ValidationError('Original representative changed')
        sc_key=poi.lower()+'_full' if poi in {'SIK3','RPS6KB2'} else poi
        records.append({'poi':poi,'graph_key':key,'expected_uniprot':uid,'historical_label':int(expected),'evaluation_role':coverage[poi]['role'],'static_graph':graph,'dynamic_graph':packets['md_contact_graphs.pt'].get(key),'saved_E2_scalar_GNN_evaluation':pred.get(poi),'saved_prediction_only_call':future.get(key),'static_screen_directory_category':screen.get(sc_key,{}).get('category'),'directory_category_is_independent_label_evidence':False,'representative_pdb':str(pose),'representative_sha256':before,'representative_selection':'First retained configured combined pose without label/score optimization','protein_selection_findings':findings,'MD_graph_poses':coverage[poi]['MD_graph_pose_count'],'assay_evidence_status':coverage[poi]['assay_evidence_status'],'source_config_matches_saved_label':cfg.get('degrader') is None if expected=='-1' else int(cfg['degrader'])==int(expected)})
    if len(records)!=19 or any(sha256(p)!=hashes[str(p)] for p in paths):raise ValidationError('Roster or source drift')
    output.mkdir(parents=True);atomic_json(output/'per_target.json',records)
    report={'status':'COMPLETE_19_TARGET_RETAINED_ARTIFACT_RECONCILIATION','targets':19,'static_representatives_checked':19,'POI_sequence_identified':sum(r['protein_selection_findings']['POI'].get('uniprot')==r['expected_uniprot'] for r in records),'legacy_VHL_selects_ELOB':sum(r['protein_selection_findings']['legacy_VHL'].get('uniprot')=='Q15370' for r in records),'current_PT_static_graphs':len(packets['ensemble_graphs.pt']),'current_PT_dynamic_graphs':len(packets['md_contact_graphs.pt']),'NPZ_export_target_sets':npz_targets,'saved_evaluation_model':'04g E2-scalar dual-encoder GNN; 17 retained rows','prediction_only_declared_writer':'04d dual-encoder GNN, without the E2 scalar branch; current table has two rows','prediction_log_matches_current_table':log_matches,'prediction_log_has_MD_by_target':log_md,'prediction_log_mismatch_evidence':'Retained predict_kinome.log records has_md=False and rounded scores about 0.0001; current table records has_md=True and different scores. This log does not establish current-table as-run lineage.','historical_scalar_cohort_expression':clean_source,'historical_scalar_runner_observed_cohort_filter':'ARTIFACT_ONLY_NO_UNKNOWN_LABEL_FILTER' if clean_source=='[g for g in static_graphs if not g.artifact]' else 'REVIEW_REQUIRED','source_input_sha256':hashes,'audit_source_sha256':sha256(Path(__file__)),'records_sha256':sha256(output/'per_target.json'),'original_files_modified':False,'new_MD_or_GNN_training':False,'limitations':['One static representative per target, not every docking pose.','Current packets are expanded artifacts, not frozen copies of the historical evaluated inputs.','NPZ sidecars may be partial exports; they are not interchangeable with the full current PT packet.','Directory categories and code comments are not independent experimental label evidence.','Prediction-only calls are not extra held-out validation targets.']};atomic_json(output/'summary.json',report);print(json.dumps({k:v for k,v in report.items() if k not in ['source_input_sha256','limitations','NPZ_export_target_sets']},indent=2))

if __name__=='__main__':
    parser=__import__('argparse').ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);main(parser.parse_args().output)
