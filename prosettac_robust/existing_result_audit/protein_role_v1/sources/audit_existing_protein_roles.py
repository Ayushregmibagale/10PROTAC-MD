"""Identify proteins in existing compiled TPRs by cached canonical sequences.

Read-only retrospective audit; never modify old groups, topologies or trajectories.
"""
import argparse,json,sys,importlib.util
from pathlib import Path
import numpy as np
import MDAnalysis as mda
from Bio.SeqUtils import seq1
from Bio.Align import PairwiseAligner,substitution_matrices
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'workflow/scripts/lib'))
from protac_qc import ValidationError,sha256,atomic_json
spec=importlib.util.spec_from_file_location('existing_index_reader',ROOT/'workflow/scripts/02_analysis/07_audit_existing_bridge_metrics.py');reader=importlib.util.module_from_spec(spec);spec.loader.exec_module(reader)
ALIGNER=PairwiseAligner();ALIGNER.mode='local';ALIGNER.substitution_matrix=substitution_matrices.load('BLOSUM62');ALIGNER.open_gap_score=-10;ALIGNER.extend_gap_score=-.5


def identify(sequence,entries):
    matches=[]
    for uid,target in entries.items():
        alignment=ALIGNER.align(target,sequence)[0]
        pairs=[(i,j) for (a,b),(c,d) in zip(*alignment.aligned) for i,j in zip(range(a,b),range(c,d))]
        if not pairs:continue
        identity=sum(target[i]==sequence[j] for i,j in pairs)/len(pairs);coverage=len(pairs)/len(sequence)
        if identity>=.90 and coverage>=.90:matches.append({'uniprot':uid,'identity':identity,'observed_CA_coverage':coverage,'canonical_first':min(i for i,j in pairs)+1,'canonical_last':max(i for i,j in pairs)+1})
    if len(matches)!=1:raise ValidationError('Protein role not uniquely identified with >=90% identity and observed-sequence coverage')
    return matches[0]


def main(args):
    output=args.output.resolve();project=args.project.resolve()
    if output.exists() or output.is_relative_to(project) or output.is_relative_to(Path('/mnt/d/PROTAC degradation Prediction')):raise ValidationError('Fresh output outside originals required')
    roster=json.loads(args.pose_roster.read_text());cache=json.loads(args.sequences.read_text());all_entries={x['primaryAccession']:x['sequence']['value'] for x in cache['results']}
    cohort=json.loads((ROOT/'workflow/config/sjf8240_publication_cohort.json').read_text())['protein_groups']
    records=[];failures=[]
    for row in roster:
        poi=row['poi'];pose=row['pose_id'];tpr=project/'random_screen'/poi/pose/'prod/prod.tpr';index=project/'analysis/random_screen'/poi/pose/'index.ndx'
        try:
            hashes={'tpr':sha256(tpr),'index':sha256(index)};u=mda.Universe(str(tpr));groups=reader.index_groups(index)
            uid=cohort[poi.removesuffix('_full')];entries={k:all_entries[k] for k in ['P40337','Q15370','Q15369',uid]};proteins={}
            for moltype in sorted(set(u.atoms.moltypes)):
                if not str(moltype).startswith('Protein'):continue
                protein=u.select_atoms('moltype '+str(moltype));ca=protein.select_atoms('name CA');sequence=''.join(seq1(a.resname,undef_code='X') for a in ca)
                if not sequence:raise ValidationError('Protein molecule lacks Cα sequence')
                found=identify(sequence,entries);found.update(moltype=str(moltype),CA_atoms=len(ca),atoms=len(protein),atom_index_first_one_based=int(protein.indices.min())+1,atom_index_last_one_based=int(protein.indices.max())+1,contiguous_indices=bool(np.array_equal(protein.indices,np.arange(protein.indices.min(),protein.indices.max()+1))))
                proteins[str(moltype)]=found
            assigned={}
            for name,ids in groups.items():
                group=u.atoms[np.asarray(ids)-1];types=sorted(set(map(str,group.moltypes)))
                assigned[name]={'atoms':len(group),'hydrogen_mass_count':int(np.sum(group.masses<1.5)),'moltypes':types,'identified_protein':proteins.get(types[0]) if len(types)==1 else None,'resnames':sorted(set(map(str,group.resnames)))}
            vhl=[p for p in proteins.values() if p['uniprot']=='P40337'];poi_matches=[p for p in proteins.values() if p['uniprot']==uid]
            if len(vhl)!=1 or len(poi_matches)!=1:raise ValidationError('True VHL and POI molecules must each be unique')
            legacy=assigned['VHL']['identified_protein'];kinase=assigned['Kinase']['identified_protein']
            records.append({'poi':poi,'pose_id':pose,'source_sha256':hashes,'proteins':proteins,'legacy_index_groups':assigned,'legacy_label_VHL_is_true_VHL':bool(legacy and legacy['uniprot']=='P40337'),'legacy_Kinase_matches_expected_POI':bool(kinase and kinase['uniprot']==uid),'sequence_identified_VHL':vhl[0],'sequence_identified_POI':poi_matches[0]})
            if sha256(tpr)!=hashes['tpr'] or sha256(index)!=hashes['index']:raise ValidationError('Compiled topology/index changed during audit')
        except (OSError,ValueError,KeyError) as exc:failures.append({'poi':poi,'pose_id':pose,'reason':str(exc)})
    output.mkdir(parents=True)
    atomic_json(output/'per_pose_roles.json',records)
    summary={'status':'COMPLETE_RETROSPECTIVE_ROLE_AUDIT' if not failures else 'INCOMPLETE_ROLE_AUDIT','poses_requested':len(roster),'poses_identified':len(records),'legacy_VHL_label_matches':sum(r['legacy_label_VHL_is_true_VHL'] for r in records),'legacy_VHL_label_mismatches':sum(not r['legacy_label_VHL_is_true_VHL'] for r in records),'legacy_POI_matches':sum(r['legacy_Kinase_matches_expected_POI'] for r in records),'failures':failures,'original_files_modified':False,'new_MD_or_docking_executed':False,'sequences_sha256':sha256(args.sequences),'roster_sha256':sha256(args.pose_roster),'source_sha256':sha256(Path(__file__)),'records_sha256':sha256(output/'per_pose_roles.json'),'MDAnalysis_version':mda.__version__,'interpretation':'This identifies actual proteins in the existing compiled systems. Reproducing a legacy metric does not validate its biological role labels. Correcting analysis selections can reuse existing trajectories; no MD rerun is inferred.'}
    atomic_json(output/'summary.json',summary);print(json.dumps(summary,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['project','output','pose-roster','sequences']:p.add_argument('--'+name,type=Path,required=True)
    main(p.parse_args())
