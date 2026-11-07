"""Regression coverage for failures that could corrupt a long computation campaign."""
import csv
import importlib.util
import json
import shutil
import sys
from pathlib import Path
import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts/lib'))
from protac_qc import ValidationError, parse_scores, pose_features, validate_outputs, verify_receipt, atomic_json, sha256, read_params

def load(name,path):
    spec=importlib.util.spec_from_file_location(name,ROOT/path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module

batch=load('batch','scripts/01_setup/02_batch_prosettac.py')
attempt=load('attempt','scripts/01_setup/03_record_attempt.py')
placement=load('placement','scripts/01_setup/01_place_warhead_batch.py')
preflight=load('preflight','scripts/01_setup/01b_preflight.py')
supplier=load('supplier','scripts/01_setup/00c_verify_supplier_registry.py')
extract=load('extract','scripts/02_analysis/01_extract_bridging_features.py')
manifest=load('manifest','scripts/01_setup/00b_build_merged_manifest.py')

@pytest.fixture
def native():
    return ROOT/'reference_structures/validated/GSK215'

@pytest.fixture
def run(tmp_path,native):
    shutil.copy2(native/'selection.json',tmp_path/'selection.json')
    shutil.copy2(ROOT/'molecules/GSK215/full.smi',tmp_path/'full.smi')
    return tmp_path

def fill_outputs(run,native,score=-123.4):
    base=run/'Patchdock_Results'; base.mkdir()
    shutil.copy2(native/'native_complex.pdb',base/'combined_1_0001.pdb')
    (base/'score.sc').write_text(f'SEQUENCE: \nSCORE: total_score description\nSCORE: {score} combined_1_0001\n')
    atomic_json(run/'engine_status.json',{'status':'COMPLETE','score_role':'raw',
        'engine_sha256':{'synthetic_worker':'regression_fixture_not_scientific_output'},'seed':42,
        'stages':{stage:{'attempted':1,'succeeded':1,'rejected':0,'failed':0}
                  for stage in ['global_sampling','local_docking','linker_modeling','scoring']}})

def receipt(run):
    rows=validate_outputs(run,json.loads((run/'selection.json').read_text()))
    atomic_json(run/'completion.json',{'status':'COMPLETE','exit_code':0,
        'selection_sha256':sha256(run/'selection.json'),
        'input_sha256':{'full.smi':sha256(run/'full.smi')},
        'score_sha256':sha256(run/'Patchdock_Results/score.sc'),
        'engine_status_sha256':sha256(run/'engine_status.json'),
        'pose_sha256':{r['description']:r['pdb_sha256'] for r in rows}})

def test_supplier_formulas_and_stereochemistry(tmp_path):
    (tmp_path/'config').mkdir();shutil.copy2(ROOT/'config/molecules.json',tmp_path/'config/molecules.json')
    results=supplier.verify(tmp_path)
    assert sum(r['status']=='IDENTITY_VERIFIED' for r in results.values())==5
    assert results['SJF_delta']['catalog']=='7267'

def test_published_delta_control_retains_explicit_model_evidence():
    folder=ROOT/'prepared_inputs/SJF_delta/MAPK13_O15264_published_overlay'
    registry=json.loads((ROOT/'config/molecules.json').read_text())['compounds']
    result=preflight.validate_pair(folder,registry)
    assert set(result['issues'])=={'BENCHMARK_RELEASE_PENDING','MODELED_HEAD_REVIEW_PENDING'}
    info=json.loads((folder/'preparation.json').read_text())
    evidence=info['reference_model_evidence']
    assert evidence['status']=='MODELED_STARTING_CONTROL_NOT_EXPERIMENTAL_COMPLEX'
    assert len(evidence['fit_positions'])==83
    assert evidence['pocket_rmsd_A']<2.0
    assert evidence['head_min_contact_A']>=1.5
    assert evidence['foretinib_template'].startswith('EPHA2/5IA4')
    assert evidence['target_template'].startswith('MAPK13/4EYJ')
    assert not info['production_ready']

def test_klifs_pocket_correspondence_rejects_ambiguous_and_sparse_maps():
    delta=load('delta_pocket_control','scripts/01_setup/00j_prepare_delta_p38_control.py')
    rows=[{'index':i,'Xray_position':str(i)} for i in range(1,86)]
    ca={str(i):np.array([i,0.,0.]) for i in range(1,86)}
    with pytest.raises(ValidationError,match='Duplicate'):
        delta.pocket_pairs(rows+[rows[0]],rows,ca,ca)
    with pytest.raises(ValidationError,match='coverage'):
        delta.pocket_pairs(rows,rows,ca,dict(list(ca.items())[:60]))
    keys,x,y=delta.pocket_pairs(rows,rows,ca,ca)
    assert len(keys)==85
    np.testing.assert_array_equal(x,y)

def test_native_is_supplier_gsk215(native):
    from rdkit import Chem
    registry=json.loads((ROOT/'config/molecules.json').read_text())['compounds']
    assert Chem.MolToInchiKey(Chem.MolFromMolFile(str(native/'bound_ligand.sdf')))==registry['GSK215']['inchikey']

def test_head_split_preserves_rings_and_bound_coordinates(native):
    from rdkit import Chem
    full=Chem.MolFromMolFile(str(native/'bound_ligand.sdf'))
    metadata=json.loads((native/'heads.json').read_text())
    sets=[]
    for role,n in [('poi',36),('vhl',31)]:
        mol=Chem.MolFromMolFile(str(native/(role+'_head.sdf')))
        ids=metadata[role]['full_ccd_indices_zero_based'];sets.append(set(ids))
        assert mol.GetNumAtoms()==n
        assert full.HasSubstructMatch(mol,useChirality=True)
        np.testing.assert_allclose(mol.GetConformer().GetPositions(),full.GetConformer().GetPositions()[ids],atol=1e-4)
    assert not sets[0]&sets[1]
    assert sum(len(s) for s in sets)==67

def test_native_contact_control(native):
    feats=pose_features(native/'native_complex.pdb',json.loads((native/'selection.json').read_text()))
    assert 2.4<feats['d_vhl']<2.6 and 2.9<feats['d_kin']<3.1
    assert feats['bridging_5A'] and not feats['severe_overlap']

def test_selection_overlap_rejected(native):
    s=json.loads((native/'selection.json').read_text());s['poi']=s['vhl']
    with pytest.raises(ValidationError,match='Overlapping'):pose_features(native/'native_complex.pdb',s)

def test_absent_chain_not_zero_distance(native):
    s=json.loads((native/'selection.json').read_text());s['poi']['chains']=['Z']
    with pytest.raises(ValidationError,match='Empty'):pose_features(native/'native_complex.pdb',s)

def test_score_sequence_and_repeated_headers(run):
    p=run/'score.sc';p.write_text('SEQUENCE: \nSCORE: total_score description\nSCORE: -4 a\nSCORE: total_score description\nSCORE: -5 b\n')
    assert [r['total_score'] for r in parse_scores(p)]==[-4,-5]

@pytest.mark.parametrize('data',['SCORE: -5 a\n','SCORE: total_score description\nSCORE: NaN a\n',
    'SCORE: total_score description\nSCORE: -5 a\nSCORE: -6 a\n',
    'SCORE: total_score description\nSCORE: -5\n',
    'SCORE: total_score description\nSCORE: -5 ../escape\n','SEQUENCE: \n'])
def test_invalid_scores_never_become_negatives(run,data):
    p=run/'score.sc';p.write_text(data)
    with pytest.raises(ValidationError):parse_scores(p)

def test_missing_generated_pose_rejected(run,native):
    fill_outputs(run,native);(run/'Patchdock_Results/combined_1_0001.pdb').unlink()
    with pytest.raises(ValidationError,match='mismatch'):validate_outputs(run,json.loads((run/'selection.json').read_text()))

def test_placeholder_rejected(run,native):
    fill_outputs(run,native,score=-1)
    with pytest.raises(ValidationError,match='placeholder'):validate_outputs(run,json.loads((run/'selection.json').read_text()))

def test_score_file_alone_never_completes(run,native):
    fill_outputs(run,native)
    with pytest.raises(FileNotFoundError):verify_receipt(run,json.loads((run/'selection.json').read_text()))

def test_receipt_and_geometry_extract(run,native):
    fill_outputs(run,native);receipt(run)
    rows,features=extract.extract_one(run)
    assert features['n_poses']==1 and features['prostac_bridging_frac']==1 and features['bridge_pct']==100
    assert features['r_dvhl_dkin'] is None

@pytest.mark.parametrize('target',['full.smi','Patchdock_Results/score.sc','Patchdock_Results/combined_1_0001.pdb','selection.json'])
def test_post_completion_tampering_rejected(run,native,target):
    fill_outputs(run,native);receipt(run)
    with (run/target).open('a') as f:f.write('\n')
    with pytest.raises(ValidationError):verify_receipt(run,json.loads((run/'selection.json').read_text()))

def test_nonzero_exit_terminal_failure(run):
    r=attempt.record(run,[sys.executable,'-c','raise SystemExit(7)'],10)
    assert r['status']=='FAILED' and r['exit_code']==7
    with pytest.raises(ValidationError,match='already used'):attempt.record(run,[sys.executable,'-c','pass'],10)

def test_success_exit_without_outputs_is_failure(run):
    r=attempt.record(run,[sys.executable,'-c','pass'],10)
    assert r['status']=='FAILED'

def test_timeout_records_failure(run):
    r=attempt.record(run,[sys.executable,'-c','import time;time.sleep(30)'],0.05)
    assert r['status']=='FAILED' and 'timeout' in r['reason']

def test_manifest_filtered_union_and_change_invalidation(tmp_path):
    p=tmp_path/'manifest.tsv'
    a={'protac':'GSK215','gene':'PTK2','uid':'Q05397','input_hash':'a'}
    b={'protac':'GSK215','gene':'EGFR','uid':'P00533','input_hash':'b'}
    batch.build_manifest(p,[a,b]); rows=batch.build_manifest(p,[dict(a,input_hash='new')])
    assert len(rows)==2 and {r['input_hash'] for r in rows}=={'b','new'}
    assert all(r['status']=='BLOCKED' for r in rows)

def test_alignment_inverse_and_rmsd():
    rng=np.random.default_rng(19); mobile=rng.normal(size=(40,3))
    r=np.array([[0,-1,0],[1,0,0],[0,0,1]]);t=np.array([20,30,-11])
    ref=mobile@r.T+t; fitted,translation=placement.kabsch(mobile,ref)
    np.testing.assert_allclose((ref-translation)@fitted,mobile,atol=1e-12)

def test_native_placement_correspondence(native):
    _,ca=placement.domain_atoms(native/'protein.pdb','D',1,9999)
    from rdkit import Chem
    head=Chem.MolFromMolFile(str(native/'poi_head.sdf')).GetConformer().GetPositions()
    _,_,meta=placement.align_domain(ca,ca,head)
    assert meta['identity']==1 and meta['pocket_rmsd_A']<1e-10

def test_full_length_structure_rejected(native):
    with pytest.raises(ValidationError):placement.domain_atoms(native/'protein.pdb','Z',1,9999)

def test_unknown_chemotype_not_substituted(tmp_path):
    with pytest.raises(ValidationError,match='No validated'):placement.prepare(tmp_path,{'protac':'JNJ1013'})

def test_legacy_config_rejected(run):
    p=run/'Protac_params.txt';p.write_text('Structures: a b\nChains: ABC A\nSmiles: wrong.smi\n')
    with pytest.raises(ValidationError,match='Missing runner keys'):read_params(p)

def test_prepared_native_candidate_only_needs_benchmark_release():
    registry=json.loads((ROOT/'config/molecules.json').read_text())['compounds']
    folder=ROOT/'prepared_inputs/GSK215/PTK2_Q05397'
    if not folder.exists(): pytest.skip('Native candidate generated during installation')
    result=preflight.validate_pair(folder,registry)
    assert result['issues']==['BENCHMARK_RELEASE_PENDING']

def test_missing_config_is_discovered_and_blocked(tmp_path):
    folder=tmp_path/'GSK215/PTK2_Q05397';folder.mkdir(parents=True)
    rows=batch.get_all_inputs(tmp_path)
    assert len(rows)==1
    registry=json.loads((ROOT/'config/molecules.json').read_text())['compounds']
    assert preflight.validate_pair(folder,registry)['status']=='BLOCKED'

def test_metadata_join_respects_chain():
    import pandas as pd
    files=pd.DataFrame([{'klifs_name':'K','pdb_id':'1abc','chain':'B','klifs_id':None,'filepath':'input.pdb'}])
    metadata=pd.DataFrame([{'klifs_id':1,'gene':'K','pdb_id':'1abc','chain':'A','dfg_state':'in'},
                           {'klifs_id':2,'gene':'K','pdb_id':'1abc','chain':'B','dfg_state':'out'}])
    result=manifest.match_files_to_manifest(files,metadata)
    assert len(result)==1 and result.iloc[0]['dfg_state']=='out' and result.iloc[0]['klifs_id']==2

def test_ambiguous_chain_is_not_given_arbitrary_dfg():
    import pandas as pd
    files=pd.DataFrame([{'klifs_name':'K','pdb_id':'1abc','chain':'A','klifs_id':None,'filepath':'input.pdb'}])
    metadata=pd.DataFrame([{'klifs_id':1,'gene':'K','pdb_id':'1abc','chain':'A','dfg_state':'in'},
                           {'klifs_id':2,'gene':'K','pdb_id':'1abc','chain':'A','dfg_state':'out'}])
    result=manifest.match_files_to_manifest(files,metadata)
    assert result.iloc[0]['metadata_status']=='AMBIGUOUS_METADATA'
    assert pd.isna(result.iloc[0]['dfg_state'])
    files.loc[0,'klifs_id']=2
    exact=manifest.match_files_to_manifest(files,metadata)
    assert exact.iloc[0]['metadata_status']=='EXACT_ID' and exact.iloc[0]['dfg_state']=='out'

def test_tandem_domains_survive_manifest_merge():
    import pandas as pd
    kinase_map=pd.DataFrame([{'klifs_name':'RSK_N','uniprot_id':'U','hgnc_gene':'RSK'},
                             {'klifs_name':'RSK_C','uniprot_id':'U','hgnc_gene':'RSK'}])
    best=pd.DataFrame([{'klifs_name':'RSK_N'},{'klifs_name':'RSK_C'}])
    af2=pd.DataFrame([{'uniprot_id':'U','hgnc_gene':'RSK','af2_path':None}])
    result=manifest.build_merged(kinase_map,best,af2)
    assert len(result)==2 and set(result['domain_id'])=={'RSK_N','RSK_C'}

def test_source_sidecar_hash_detects_modified_structure(tmp_path,monkeypatch):
    import hashlib
    file=tmp_path/'123_out.pdb';file.write_text('original')
    file.with_suffix('.provenance.json').write_text(json.dumps({'gene':'K','pdb_id':'1abc','chain':'A',
         'klifs_id':123,'sha256':hashlib.sha256(b'original').hexdigest()}))
    monkeypatch.setattr(manifest,'STRUCT_DIR',tmp_path)
    assert manifest.scan_klifs_files().iloc[0]['klifs_id']==123
    file.write_text('changed')
    with pytest.raises(ValueError,match='hash mismatch'):manifest.scan_klifs_files()

def test_foretinib_and_vhl_exact_bound_fragments():
    from rdkit import Chem
    refs=ROOT/'reference_structures/validated'
    info=json.loads((refs/'SJF8240/heads.json').read_text())
    for role,parent,n in [('poi','foretinib_parent',37),('vhl','VH032_parent',30)]:
        head=Chem.MolFromMolFile(str(refs/'SJF8240'/(role+'_head.sdf')))
        original=Chem.MolFromMolFile(str(refs/parent/'bound_ligand.sdf'))
        assert head.GetNumAtoms()==n
        indices=info[role]['parent_indices_zero_based']
        np.testing.assert_allclose(head.GetConformer().GetPositions(),original.GetConformer().GetPositions()[indices],atol=1e-4)
        full=Chem.MolFromSmiles((ROOT/'molecules/SJF8240/full.smi').read_text().strip())
        assert full.HasSubstructMatch(head,useChirality=True)

def test_foretinib_candidate_is_not_released():
    folder=ROOT/'prepared_inputs/SJF8240/MET_P08581'
    if not folder.exists():pytest.skip('Foretinib control generated during installation')
    registry=json.loads((ROOT/'config/molecules.json').read_text())['compounds']
    assert preflight.validate_pair(folder,registry)['issues']==['BENCHMARK_RELEASE_PENDING']

def test_delta_has_a_distinct_vhl_attachment():
    from rdkit import Chem
    registry=json.loads((ROOT/'config/molecules.json').read_text())['compounds']
    delta=Chem.MolFromSmiles(registry['SJF_delta']['smiles'])
    head=Chem.MolFromMolFile(str(ROOT/'reference_structures/validated/SJF8240/vhl_head.sdf'))
    assert not delta.HasSubstructMatch(head,useChirality=True)

def test_child_worker_failure_not_hidden_by_surviving_pose(run,native):
    fill_outputs(run,native)
    p=run/'engine_status.json';report=json.loads(p.read_text())
    report['stages']['local_docking'].update(attempted=2,failed=1)
    atomic_json(p,report)
    with pytest.raises(ValidationError,match='Failed child'):validate_outputs(run,json.loads((run/'selection.json').read_text()))

def test_missing_engine_accounting_rejected(run,native):
    fill_outputs(run,native);(run/'engine_status.json').unlink()
    with pytest.raises(FileNotFoundError):validate_outputs(run,json.loads((run/'selection.json').read_text()))

def test_output_stereochemistry_inversion_rejected(run,native):
    from protac_qc import validate_ligand_identity
    fill_outputs(run,native)
    p=run/'Patchdock_Results/combined_1_0001.pdb'
    lines=[]
    for line in p.read_text().splitlines():
        if line.startswith('HETATM'):line=line[:30]+f'{-float(line[30:38]):8.3f}'+line[38:]
        lines.append(line)
    p.write_text('\n'.join(lines)+'\n')
    selection=json.loads((run/'selection.json').read_text())
    with pytest.raises(ValidationError,match='stereochemistry'):validate_ligand_identity(p,selection['protac'],run/'full.smi')

@pytest.mark.parametrize('name,role,parent,count',[('SJF1521','poi','SJF1521_parent',27),
    ('JNJ1013','poi','JNJ1013_parent',22),('SJF_delta','vhl','VH032_parent',28)])
def test_modified_head_preserves_bound_scaffold_and_supplier_stereo(name,role,parent,count):
    from rdkit import Chem
    dest=ROOT/'modeled_reference_candidates'/name
    meta=json.loads((dest/(role+'_model_provenance.json')).read_text())
    assert not meta['production_ready'] and not meta['experimental_coordinates_for_unmatched_atoms']
    assert meta['common_scaffold_atoms']==count and meta['seeds_valid']==12
    head=Chem.MolFromMolFile(str(dest/(role+'_head.sdf')))
    reference=Chem.MolFromMolFile(str(ROOT/'reference_structures/validated'/parent/'bound_ligand.sdf'))
    np.testing.assert_allclose(head.GetConformer().GetPositions()[meta['core_head_indices_zero_based']],
        reference.GetConformer().GetPositions()[meta['core_reference_indices_zero_based']],atol=1e-4)
    registry=json.loads((ROOT/'config/molecules.json').read_text())['compounds']
    supplier_mol=Chem.MolFromSmiles(registry[name]['smiles'])
    Chem.AssignStereochemistryFrom3D(head,replaceExistingTags=True)
    assert supplier_mol.HasSubstructMatch(head,useChirality=True)

@pytest.mark.parametrize('name,run_name',[('SJF1521','EGFR_P00533'),('JNJ1013','IRAK1_P51617'),('SJF_delta','MET_P08581')])
def test_modeled_controls_require_explicit_review(name,run_name):
    registry=json.loads((ROOT/'config/molecules.json').read_text())['compounds']
    result=preflight.validate_pair(ROOT/'prepared_inputs'/name/run_name,registry)
    assert set(result['issues'])=={'BENCHMARK_RELEASE_PENDING','MODELED_HEAD_REVIEW_PENDING'}


def test_no_valid_poses_is_not_a_biological_negative(run):
    output=run/'Patchdock_Results';output.mkdir()
    (output/'score.sc').write_text('SEQUENCE: \nSCORE: total_score description\n')
    atomic_json(run/'engine_status.json',{'status':'COMPLETE_NO_VALID_POSES','score_role':'raw',
        'engine_sha256':{'fixture':'fixture'},'seed':73,
        'stages':{'global_sampling':{'attempted':5,'succeeded':2,'rejected':3,'failed':0},
                  'local_docking':{'attempted':4,'succeeded':4,'rejected':0,'failed':0},
                  'linker_modeling':{'attempted':4,'succeeded':0,'rejected':4,'failed':0},
                  'scoring':{'attempted':0,'succeeded':0,'rejected':0,'failed':0}}})
    assert validate_outputs(run,json.loads((run/'selection.json').read_text()))==[]
    receipt(run)
    data=json.loads((run/'completion.json').read_text());data['status']='COMPLETE_NO_VALID_POSES'
    atomic_json(run/'completion.json',data)
    rows,features=extract.extract_one(run)
    assert rows==[] and features['linker_geometry_yield']==0
    assert features['bridge_pct'] is None and features['n_protein_poses']==4


def test_tandem_domain_manifest_rows_do_not_collapse(tmp_path):
    common={'protac':'GSK215','gene':'JAK1','uid':'P23458','local_dir':'fixture','input_hash':'abc'}
    pairs=[dict(common,run_name='JAK1_P23458__kinase1'),dict(common,run_name='JAK1_P23458__kinase2')]
    assert len(batch.build_manifest(tmp_path/'manifest.tsv',pairs))==2


def test_stage_denominator_cannot_drop_rejected_poses(run,native):
    fill_outputs(run,native)
    data=json.loads((run/'engine_status.json').read_text())
    data['stages']['local_docking']={'attempted':10,'succeeded':10,'failed':0,'rejected':0}
    atomic_json(run/'engine_status.json',data)
    with pytest.raises(ValidationError,match='denominator'):
        validate_outputs(run,json.loads((run/'selection.json').read_text()))


def test_md_protocol_preserves_checkpoint_velocities(tmp_path):
    md=load('prepare_md','scripts/01_setup/06_prepare_md.py')
    md.write_mdps(tmp_path,107,20)
    assert 'gen-vel = yes' in (tmp_path/'nvt.mdp').read_text()
    for stage in ['npt_restrained','npt','prod']:
        text=(tmp_path/(stage+'.mdp')).read_text()
        assert 'gen-vel = no' in text and 'continuation = yes' in text
        assert 'compressed-x-grps = SOLUTE' in text
    assert 'POSRES' not in (tmp_path/'prod.mdp').read_text()
    assert 'nsteps = 10000000' in (tmp_path/'prod.mdp').read_text()


def test_ce_correspondence_transfers_native_domain(native):
    from rdkit import Chem
    _,reference=placement.domain_atoms(native/'protein.pdb','D',-999,9999)
    mobile=[dict(a,xyz=np.asarray(a['xyz'])+np.array([19.,-7.,3.])) for a in reference]
    head=Chem.MolFromMolFile(str(native/'poi_head.sdf')).GetConformer().GetPositions()
    rotation,translation,info=placement.align_domain(mobile,reference,head,'ce')
    assert info['core_rmsd_A']<1e-6 and info['pocket_coverage']>0.99
    assert np.max(np.abs(np.array([a['xyz'] for a in mobile])@rotation.T+translation-
                         np.array([a['xyz'] for a in reference])))<1e-5


def test_baseline_local_worker_failure_is_propagated(tmp_path,monkeypatch):
    monkeypatch.chdir(tmp_path)
    import ast,subprocess
    source=ast.parse((ROOT/'baseline_engine/LOCAL.py').read_text())
    source.body=[n for n in source.body if not (isinstance(n,ast.ImportFrom) and n.level)]
    namespace={'Cluster':object};exec(compile(source,'LOCAL.py','exec'),namespace)
    with pytest.raises(subprocess.CalledProcessError):
        namespace['LOCAL'](n_workers=2).runBatchCommands([f'{sys.executable} -c "raise SystemExit(7)"'],mem=1)


def test_baseline_seed_and_rotation_are_repeatable(tmp_path):
    import subprocess,os
    module=ROOT/'baseline_engine/pyrosetta_adapter.py'
    command=[sys.executable,'-c',f'import importlib.util; s=importlib.util.spec_from_file_location("a",{str(module)!r});m=importlib.util.module_from_spec(s);s.loader.exec_module(m);print(m.publication_seed("stage"));print(m._random_rotation().tolist())']
    env=dict(os.environ,PROTAC_SEED='1921',PROTAC_ATTEMPT_ROOT=str(tmp_path))
    first=subprocess.check_output(command,cwd=tmp_path,env=env)
    assert first==subprocess.check_output(command,cwd=tmp_path,env=env)
    env['PROTAC_TASK_KEY']='different_worker_task'
    assert first!=subprocess.check_output(command,cwd=tmp_path,env=env)
    env.pop('PROTAC_TASK_KEY')
    env['PROTAC_SEED']='1922'
    assert first!=subprocess.check_output(command,cwd=tmp_path,env=env)


def test_sjf1521_parameter_alias_preserves_graph_and_qm_stereochemistry(tmp_path):
    params=load('params','scripts/01_setup/04_parameterize_ligands.py')
    source=ROOT/'ligand_params/SJF1521/neutral_reference'
    target=tmp_path/'parameters';shutil.copytree(source,target)
    meta=json.loads((target/'atom_mapping.json').read_text())
    assert meta['tool_atom_name_aliases']=={'Cl65':'CL65'}
    report=params.validate_parameters(target,'SJF1521',meta)
    assert report['atoms']==137 and report['output_stereochemistry_verified']


def test_baseline_parallel_score_writes_keep_all_rows(tmp_path,monkeypatch):
    from concurrent.futures import ThreadPoolExecutor
    module=load('baseline_adapter','baseline_engine/pyrosetta_adapter.py')
    monkeypatch.chdir(tmp_path)
    def write(worker):
        for i in range(20):module._append_score_sc(worker+i/100,f'pose_{worker}_{i}')
    with ThreadPoolExecutor(max_workers=4) as pool:list(pool.map(write,range(1,5)))
    assert len(parse_scores(tmp_path/'score.sc'))==80
    assert (tmp_path/'score.sc').read_text().count('total_score description')==1


def test_baseline_head_extraction_uses_names_not_pdb_order_or_hydrogens(tmp_path,monkeypatch):
    from rdkit import Chem
    params=load('official_params','baseline_engine/official_parameters.py')
    monkeypatch.setenv('PROTAC_ATTEMPT_ROOT',str(tmp_path))
    paths=[];lines=[];expected=[]
    for role,smiles in enumerate(['CO','CN']):
        molecule=Chem.MolFromSmiles(smiles);conf=Chem.Conformer(2)
        coords=np.array([[role*10.,0.,0.],[role*10.+1.4,0.,0.]])
        for i,xyz in enumerate(coords):conf.SetAtomPosition(i,xyz)
        molecule.AddConformer(conf);path=tmp_path/f'head{role}.sdf'
        w=Chem.SDWriter(str(path));w.write(molecule);w.close();paths.append(path);expected.extend(coords)
        mapping={str(i):a.GetSymbol()+str(i+1) for i,a in enumerate(molecule.GetAtoms())}
        atomic_json(tmp_path/f'PT{role}_parameter_provenance.json',{'source_sha256':sha256(path),'heavy_atoms':2,'all_atom_index_to_name':mapping})
        for i,atom in enumerate(molecule.GetAtoms()):
            x,y,z=coords[i]
            lines.append(f'HETATM{i+1:5d} {mapping[str(i)]:>4s} PT{role} X   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {atom.GetSymbol():>2s}  ')
        lines.append(f'HETATM{9:5d} {"H99":>4s} PT{role} X   1    {9.:8.3f}{9.:8.3f}{9.:8.3f}  1.00  0.00           H  ')
    pdb=tmp_path/'docked.pdb';pdb.write_text('\n'.join(reversed(lines))+'\n')
    combined=params.docked_heads(str(pdb),paths)
    assert np.allclose(combined.GetConformer().GetPositions(),expected)


def test_baseline_qc_rejects_overlapping_bridge_without_biological_label(run,native):
    qc=load('baseline_output_audit','baseline_engine/audit_outputs.py')
    fill_outputs(run,native)
    good=qc.audit(run)
    assert good['qc_valid_poses']==1 and good['production_release'] is False
    p=run/'Patchdock_Results/combined_1_0001.pdb';atoms=__import__('protac_qc').pdb_atoms(p)
    selector=json.loads((run/'selection.json').read_text())
    selected=[a for a in atoms if a['record']==selector['protac']['record'] and a['chain'] in selector['protac']['chains'] and a['element'] not in {'H','D'}]
    protein=next(a for a in atoms if a['record']=='ATOM' and a['chain'] in selector['poi']['chains'])
    shift=np.array(protein['xyz'])-np.array(selected[0]['xyz']);changed=[]
    for line in p.read_text().splitlines():
        if line.startswith('HETATM') and line[21] in selector['protac']['chains']:
            xyz=np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])])+shift
            line=line[:30]+''.join(f'{v:8.3f}' for v in xyz)+line[54:]
        changed.append(line)
    p.write_text('\n'.join(changed)+'\n')
    bad=qc.audit(run)
    assert bad['qc_valid_poses']==0 and bad['biological_degrader_label'] is None
    assert bad['conditional_valid_bridge_fraction_5A'] is None


def test_pose_archive_restores_exact_bytes_and_detects_tampering(tmp_path,native):
    archive=load('archive','scripts/02_analysis/03_archive_pose_ensemble.py')
    source=tmp_path/'source';source.mkdir();text=(native/'native_complex.pdb').read_text()
    first=source/'combined_1_1_0001.pdb';first.write_text(text)
    lines=text.splitlines()
    for i,line in enumerate(lines):
        if line.startswith('ATOM  '):
            lines[i]=line[:30]+f'{-0.0:8.3f}'+line[38:54]+'  0.99'+line[60:];break
    second=source/'pd.1_docking_0001.pdb';second.write_text('\n'.join(lines))
    packed=tmp_path/'packed';report=archive.pack(source,packed,1)
    assert report['source_deleted'] is False
    restored=tmp_path/'restored';archive.restore(packed,restored)
    assert first.read_bytes()==(restored/first.name).read_bytes()
    assert second.read_bytes()==(restored/second.name).read_bytes()
    chunk=next(packed.glob('*.npz'));chunk.write_bytes(chunk.read_bytes()+b'tampered')
    with pytest.raises(ValidationError,match='chunk changed'):
        archive.restore(packed,tmp_path/'tampered_restore')


def test_fixed_head_linker_preserves_native_geometry_and_rejects_impossible_placement(native):
    from rdkit import Chem
    model=load('fixed_head_linker','baseline_engine/constrained_linker.py')
    full=Chem.MolFromSmiles((ROOT/'molecules/GSK215/full.smi').read_text().strip())
    heads=[Chem.MolFromMolFile(str(native/n),removeHs=True) for n in ['vhl_head.sdf','poi_head.sdf']]
    targets=[h.GetConformer().GetPositions() for h in heads]
    mol,proof=model.constrained_linker(full,heads,targets,392,attempts=20)
    assert mol is not None and proof['accepted']
    for mapping,target in zip(proof['head_mappings'],targets):
        assert np.max(np.linalg.norm(mol.GetConformer().GetPositions()[mapping]-target,axis=1))<1e-5
    impossible=[targets[0],targets[1]+np.array([1000.,0.,0.])]
    mol,proof=model.constrained_linker(full,heads,impossible,392,attempts=2)
    assert mol is None and not proof['accepted']


def test_protonation_candidates_keep_supplier_graph_and_stereochemistry():
    from rdkit import Chem
    from rdkit.Chem.MolStandardize import rdMolStandardize
    report=json.loads((ROOT/'molecules/protonation_assessment_pH7p4/assessment.json').read_text())
    assert report['primary_MD_state_selected'] is False
    for compound,info in report['compounds'].items():
        for state in info['states'].values():
            parent=rdMolStandardize.Uncharger().uncharge(Chem.MolFromSmiles(state['smiles']))
            assert Chem.MolToInchiKey(parent)==info['supplier_inchikey']
            assert state['predicted_population'] is None
        assert all(s['tool_stats']['protonation']['fallback_used']==0 for s in info['scenarios'].values())


def test_native_control_manifest_uses_filename_when_old_metadata_lacks_identity(tmp_path):
    folder=tmp_path/'GSK215/PTK2_Q05397';folder.mkdir(parents=True)
    atomic_json(folder/'preparation.json',{'status':'PREPARATION_CANDIDATE'})
    rows=batch.get_all_inputs(tmp_path)
    assert (rows[0]['gene'],rows[0]['uid'])==('PTK2','Q05397')


def test_post_score_qc_preserves_raw_stage_denominators(run,native):
    fill_outputs(run,native)
    data=json.loads((run/'engine_status.json').read_text())
    for stage in ['local_docking','linker_modeling','scoring']:
        data['stages'][stage]={'attempted':2,'succeeded':2,'rejected':0,'failed':0}
    data['stages']['pose_qc']={'attempted':2,'succeeded':1,'rejected':1,'failed':0}
    atomic_json(run/'engine_status.json',data)
    assert len(validate_outputs(run,json.loads((run/'selection.json').read_text())))==1


def test_raw_baseline_score_tampering_is_detected(run,native):
    fill_outputs(run,native)
    raw=run/'engine_raw';(raw/'Patchdock_Results').mkdir(parents=True)
    shutil.copy2(run/'Patchdock_Results/score.sc',raw/'Patchdock_Results/score.sc')
    atomic_json(raw/'baseline_output_qc.json',{'fixture':'integrity test only'})
    data=json.loads((run/'engine_status.json').read_text())
    data.update(engine_kind='existing_pRosettaC_baseline_repaired',raw_score_sha256=sha256(raw/'Patchdock_Results/score.sc'),raw_output_qc_sha256=sha256(raw/'baseline_output_qc.json'))
    atomic_json(run/'engine_status.json',data)
    validate_outputs(run,json.loads((run/'selection.json').read_text()))
    (raw/'Patchdock_Results/score.sc').write_text('changed')
    with pytest.raises(ValidationError,match='Raw baseline score'):
        validate_outputs(run,json.loads((run/'selection.json').read_text()))


def test_attempt_sigterm_records_interrupt_and_terminates_child(run):
    import subprocess,signal,time,os
    command=[sys.executable,str(ROOT/'scripts/01_setup/03_record_attempt.py'),'--run-dir',str(run),
        '--timeout-seconds','120','--',sys.executable,'-c','import time;time.sleep(30)']
    process=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    deadline=time.monotonic()+10;child=None
    try:
        while time.monotonic()<deadline:
            if (run/'attempt.json').exists():
                data=json.loads((run/'attempt.json').read_text());child=data.get('process_id')
                if child:break
            time.sleep(0.05)
        assert child is not None
        process.send_signal(signal.SIGTERM);process.communicate(timeout=15)
        completed=json.loads((run/'completion.json').read_text())
        assert completed['status']=='INTERRUPTED'
        with pytest.raises(ProcessLookupError):os.kill(child,0)
    finally:
        if process.poll() is None:process.kill();process.wait()


def test_pocket_map_cache_requires_unique_indices_and_residues():
    cache=load('pocket_cache_checks','scripts/01_setup/00k_cache_klifs_pockets.py')
    rows=[{'index':i,'Xray_position':str(i),'KLIFS_position':f'test.{i}'} for i in range(1,86)]
    assert cache.validate(rows)==rows
    with pytest.raises(ValidationError,match='85'):
        cache.validate(rows+[rows[0]])
    rows[1]=dict(rows[1],Xray_position='1')
    with pytest.raises(ValidationError,match='one-to-one'):
        cache.validate(rows)


def test_canonical_pocket_consensus_never_resolves_conflicts_by_voting():
    correspondence=load('pocket_consensus_checks','scripts/01_setup/00l_build_pocket_correspondence.py')
    assert correspondence.consensus([{'1':10},{'1':10,'2':11}])=={'1':10,'2':11}
    with pytest.raises(ValidationError,match='Conflicting'):
        correspondence.consensus([{'1':10},{'1':11},{'1':10}])
    with pytest.raises(ValidationError,match='two pocket'):
        correspondence.consensus([{'1':10,'2':10}])
    assert correspondence.source_residue('123A')==(123,'A')
    assert correspondence.source_residue('_') is None


def test_explicit_pocket_fit_preserves_known_rigid_transform(native):
    atoms=[a for a in __import__('protac_qc').pdb_atoms(native/'protein.pdb') if a['name']=='CA']
    reference=atoms[:]
    from rdkit import Chem
    head=Chem.MolFromMolFile(str(native/'poi_head.sdf')).GetConformer().GetPositions()
    neighborhood=[i for i,a in enumerate(atoms) if np.linalg.norm(np.array(a['xyz'])-head,axis=1).min()<10]
    chosen=(neighborhood+[i for i in range(len(atoms)) if i not in neighborhood])[:85]
    positions={str(k+1):atoms[i]['resnum'] for k,i in enumerate(chosen)}
    shift=np.array([7.,-4.,2.])
    target=[dict(a,xyz=np.array(a['xyz'])+shift) for a in atoms]
    evidence={'target':{'positions':positions},'reference':{'positions':positions}}
    rotation,translation,meta=placement.align_klifs_pocket(target,reference,head,evidence)
    np.testing.assert_allclose(rotation,np.eye(3),atol=1e-10)
    np.testing.assert_allclose(translation,-shift,atol=1e-10)
    assert meta['pocket_rmsd_A']<1e-10 and meta['core_atoms']>=100
    bad=dict(positions);bad['2']=bad['1']
    with pytest.raises(ValidationError,match='one-to-one'):
        placement.align_klifs_pocket(target,reference,head,{'target':{'positions':bad},'reference':{'positions':positions}})
    sparse={k:v for k,v in positions.items() if int(k)<68}
    with pytest.raises(ValidationError,match='coverage'):
        placement.align_klifs_pocket(target,reference,head,{'target':{'positions':sparse},'reference':{'positions':positions}})


def test_md_graceful_stop_is_not_full_completion():
    md=load('md_stage_completion','scripts/01_setup/07_run_md_stage.py')
    partial='           Step           Time\n             10        0.02000\nFinished mdrun on rank 0\n'
    with pytest.raises(ValidationError,match='before requested'):
        md.validate_stage_log(partial,'nvt',None,'nsteps = 100\ndt = 0.002\n')
    assert md.validate_stage_log(partial,'nvt',10,'nsteps = 10\n')['last_logged_step']==10
    with pytest.raises(ValidationError,match='before requested'):
        md.validate_stage_log(partial,'prod',None,'init_step = 10\nnsteps = 10\n')
    invalid=partial+'LINCS WARNING\n'
    with pytest.raises(ValidationError,match='instability'):
        md.validate_stage_log(invalid,'nvt',None,'nsteps = 10\n')


def test_md_minimization_requires_positive_convergence_and_finite_force():
    md=load('md_em_completion','scripts/01_setup/07_run_md_stage.py')
    common='Step Time\n10 10.00000\nMaximum force = 100.0 on atom 12\nFinished mdrun\n'
    good=common+'Steepest Descents converged to Fmax < 500 in 10 steps\n'
    assert md.validate_stage_log(good,'em',None,'nsteps = 100\nemtol = 500\n')['force_converged']
    bad=common+'Steepest Descents did not converge to Fmax < 500 in 10 steps\n'
    with pytest.raises(ValidationError,match='force convergence'):
        md.validate_stage_log(bad,'em',None,'nsteps = 10\nemtol = 500\n')
    with pytest.raises(ValidationError,match='nonfinite'):
        md.validate_stage_log(good.replace('100.0','nan'),'em',None,'nsteps = 100\nemtol = 500\n')


def test_md_rejects_unbound_or_changed_tpr(tmp_path):
    md=load('md_tpr_binding','scripts/01_setup/07_run_md_stage.py')
    tpr=tmp_path/'nvt.tpr';tpr.write_bytes(b'first tpr')
    with pytest.raises(ValidationError,match='lacks a generation receipt'):
        md.verify_tpr_binding(tmp_path,'nvt',{'stage':'nvt'})
    atomic_json(tmp_path/'nvt.tpr_generation.json',{'intent':{'stage':'nvt'},'tpr_sha256':sha256(tpr)})
    md.verify_tpr_binding(tmp_path,'nvt',{'stage':'nvt'})
    tpr.write_bytes(b'second tpr')
    with pytest.raises(ValidationError,match='generation receipt'):
        md.verify_tpr_binding(tmp_path,'nvt',{'stage':'nvt'})


def test_md_preparation_requires_local_hashed_force_field(tmp_path):
    md=load('md_preparation_binding','scripts/01_setup/07_run_md_stage.py')
    atomic_json(tmp_path/'md_preparation.json',{'status':'TOPOLOGY_PREPARED','files':{}})
    with pytest.raises(ValidationError,match='Portable'):
        md.verify_preparation(tmp_path,'nvt')
    (tmp_path/'amber99sb-ildn.ff').mkdir();(tmp_path/'mdp').mkdir()
    for name in ['amber99sb-ildn.ff/forcefield.itp','topol.top','index.ndx','system.gro','mdp/nvt.mdp']:
        (tmp_path/name).write_text('fixture')
    files={str(p.relative_to(tmp_path)):sha256(p) for p in tmp_path.rglob('*') if p.is_file() and p.name!='md_preparation.json'}
    atomic_json(tmp_path/'md_preparation.json',{'status':'TOPOLOGY_PREPARED','files':files})
    md.verify_preparation(tmp_path,'nvt')
    (tmp_path/'amber99sb-ildn.ff/forcefield.itp').write_text('changed')
    with pytest.raises(ValidationError,match='altered'):
        md.verify_preparation(tmp_path,'nvt')


def test_md_system_lock_prevents_concurrent_stages(tmp_path,monkeypatch):
    import fcntl
    md=load('md_system_lock','scripts/01_setup/07_run_md_stage.py')
    monkeypatch.setattr(md,'_run',lambda *args:'locked execution')
    with (tmp_path/'.md_stage.lock').open('w') as held:
        fcntl.flock(held,fcntl.LOCK_EX|fcntl.LOCK_NB)
        with pytest.raises(ValidationError,match='Another MD stage'):
            md.run(tmp_path,'em',Path('/fixture/gmx'),1)
    assert md.run(tmp_path,'em',Path('/fixture/gmx'),1)=='locked execution'


def test_native_interface_benchmark_detects_decoys_and_ignores_global_translation(native,tmp_path):
    recovery=load('native_recovery','scripts/02_analysis/04_native_recovery.py')
    reference=native/'native_complex.pdb';selection=json.loads((native/'selection.json').read_text())
    original=reference.read_text().splitlines();shifted=[];decoy=[]
    for line in original:
        if line.startswith(('ATOM  ','HETATM')):
            xyz=np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])])
            translated=xyz+np.array([7.,-4.,2.])
            shifted.append(line[:30]+''.join(f'{v:8.3f}' for v in translated)+line[54:])
            if line.startswith('ATOM  ') and line[21] in selection['vhl']['chains']:
                xyz=xyz+np.array([30.,0.,0.])
            decoy.append(line[:30]+''.join(f'{v:8.3f}' for v in xyz)+line[54:])
        else:shifted.append(line);decoy.append(line)
    translated=tmp_path/'translated.pdb';translated.write_text('\n'.join(shifted)+'\n')
    wrong=tmp_path/'decoy.pdb';wrong.write_text('\n'.join(decoy)+'\n')
    good=recovery.measure(reference,translated,selection)
    assert good['vhl_CA_rmsd_after_POI_fit_A']<1e-8
    assert good['fraction_native_contacts']==1
    assert good['independent_global_recovery_demonstrated'] is False
    bad=recovery.measure(reference,wrong,selection)
    assert bad['vhl_CA_rmsd_after_POI_fit_A']>29
    assert bad['fraction_native_contacts']<0.1


def test_linker_embedding_has_explicit_internal_retry_budget(native,monkeypatch):
    from rdkit import Chem
    linker=load('linker_retry_budget','baseline_engine/constrained_linker.py')
    full=Chem.MolFromMolFile(str(native/'bound_ligand.sdf'))
    heads=[Chem.MolFromMolFile(str(native/(role+'_head.sdf'))) for role in ['vhl','poi']]
    seen=[]
    def reject(molecule,params):
        seen.append((params.maxIterations,params.trackFailures))
        return -1
    monkeypatch.setattr(linker.AllChem,'EmbedMolecule',reject)
    mol,report=linker.constrained_linker(full,heads,[h.GetConformer().GetPositions() for h in heads],17,attempts=2,embed_iterations=7)
    assert mol is None and seen==[(7,True),(7,True)]
    assert report['embedding_max_iterations_per_attempt']==7
    assert len(report['embedding_failure_counts_by_attempt'])==2
    with pytest.raises(ValidationError,match='bounded'):
        linker.constrained_linker(full,heads,[],17,embed_iterations=0)


def test_execution_snapshot_isolated_from_working_tree_and_detects_drift(tmp_path):
    campaign=load('campaign_snapshot','scripts/lib/campaign.py')
    root=tmp_path/'project';engine=root/'baseline_engine';engine.mkdir(parents=True)
    (engine/'worker.py').write_text('original worker\n')
    (engine/'docking.xml').write_text('<original/>\n')
    for relative in ['scripts/lib/protac_qc.py','scripts/lib/campaign.py','scripts/01_setup/02_batch_prosettac.py','scripts/01_setup/03_record_attempt.py']:
        p=root/relative;p.parent.mkdir(parents=True,exist_ok=True);p.write_text('fixture\n')
    expected=campaign.source_hashes(root);snapshot=tmp_path/'snapshot'
    campaign.snapshot_sources(root,snapshot,expected)
    (engine/'worker.py').write_text('changed worker\n')
    (engine/'docking.xml').write_text('<changed/>\n')
    assert (snapshot/'baseline_engine/worker.py').read_text()=='original worker\n'
    assert campaign.source_hashes(snapshot)==expected
    assert campaign.source_hashes(root)!=expected
    with pytest.raises(ValidationError,match='Sources changed'):
        campaign.snapshot_sources(root,tmp_path/'drifted_snapshot',expected)


def test_all_supplier_heads_offer_a_fixed_rotamer_frame():
    from rdkit import Chem
    linker=load('fixed_neighbor_selection','baseline_engine/constrained_linker.py')
    for compound in ['GSK215','SJF8240','SJF1521','JNJ1013','SJF_delta']:
        base=ROOT/('reference_structures/validated' if compound in {'GSK215','SJF8240'} else 'modeled_reference_candidates')/compound
        heads=[Chem.MolFromMolFile(str(base/(role+'_head.sdf'))) for role in ['vhl','poi']]
        full=Chem.MolFromSmiles((ROOT/'molecules'/compound/'full.smi').read_text().strip())
        fixed=set(linker.match_heads(full,heads)[0][0]);center=linker.fixed_head_neighbor(full,fixed)
        neighbors={a.GetIdx() for a in full.GetAtomWithIdx(center).GetNeighbors() if a.GetAtomicNum()>1}
        assert center in fixed and len(neighbors)>=2 and neighbors<=fixed
        assert full.GetAtomWithIdx(center).IsInRing()


def test_parameter_tool_rejects_moving_rotamer_orientation_frame(tmp_path):
    from rdkit import Chem
    params=load('rotamer_frame_validation','baseline_engine/official_parameters.py')
    molecule=Chem.MolFromSmiles('CCCO');conf=Chem.Conformer(4)
    for i,xyz in enumerate([[0.,0.,0.],[1.5,0.,0.],[2.2,1.3,0.],[3.6,1.3,0.]]):conf.SetAtomPosition(i,xyz)
    molecule.AddConformer(conf);other=Chem.Mol(molecule);other.GetConformer().SetAtomPosition(3,[3.6,1.3,1.])
    source=tmp_path/'ensemble.sdf'
    with Chem.SDWriter(str(source)) as writer:writer.write(molecule);writer.write(other)
    with pytest.raises(RuntimeError,match='not invariant'):
        params.generate(str(source),'PTC',str(tmp_path/'PT'),conformers=True,nbr=2)


def test_thermodynamic_diagnostics_reject_invalid_series_and_preserve_replica_unit():
    thermo=load('thermodynamic_diagnostics','scripts/02_analysis/05_md_thermodynamics.py')
    time=np.arange(0.,11.);values=np.column_stack([300.+time,1000.-time])
    report=thermo.summarize(time,values,2.,4.)
    assert report['retained_samples']==9 and len(report['block_means'])==2
    assert report['independent_replicates']==1 and report['convergence_certified'] is False
    assert sum(b['samples'] for b in report['block_means'])==9
    np.testing.assert_allclose(report['linear_slopes_per_ps'],[1.,-1.])
    with pytest.raises(ValidationError,match='non-increasing'):
        thermo.summarize([0.,0.,1.],np.ones((3,2)),0.,1.)
    with pytest.raises(ValidationError,match='Insufficient'):
        thermo.summarize(time,values,20.,4.)


def test_equilibration_contacts_stay_diagnostic():
    contact=load('equilibration_contact_scope','scripts/02_analysis/02_analyze_md_contacts.py')
    completed={'status':'COMPLETE','intent':{'stage':'npt'}}
    contact.validate_analysis_stage(completed,True)
    with pytest.raises(ValidationError,match='production'):
        contact.validate_analysis_stage(completed,False)
    with pytest.raises(ValidationError,match='completed'):
        contact.validate_analysis_stage({'status':'FAILED'},True)


def test_reduced_topology_requires_unpermuted_complete_solute():
    reduced=load('reduced_topology_order','scripts/02_analysis/06_reduce_md_topology.py')
    assert reduced.solute_indices('[ SOLUTE ]\n1 2 3\n[ System ]\n1 2 3 4\n',3)==[1,2,3]
    for text in ['[ SOLUTE ]\n1 3 2\n','[ SOLUTE ]\n1 2\n','[ SOLUTE ]\n1 2 3\n[ SOLUTE ]\n1 2 3\n']:
        with pytest.raises(ValidationError,match='ordered'):
            reduced.solute_indices(text,3)


def test_snapshot_receipt_verification_detects_code_tampering(tmp_path):
    from protac_qc import verify_execution_snapshot
    campaign=load('snapshot_receipt_campaign','scripts/lib/campaign.py')
    root=tmp_path/'original';engine=root/'baseline_engine';engine.mkdir(parents=True)
    (engine/'worker.py').write_text('original\n')
    for relative in ['scripts/lib/protac_qc.py','scripts/lib/campaign.py','scripts/01_setup/02_batch_prosettac.py','scripts/01_setup/03_record_attempt.py']:
        p=root/relative;p.parent.mkdir(parents=True,exist_ok=True);p.write_text('fixture\n')
    expected=campaign.source_hashes(root);run=tmp_path/'attempt';run.mkdir()
    campaign.snapshot_sources(root,run/'source_snapshot',expected)
    verify_execution_snapshot(run,{'engine_sha256':expected})
    (run/'source_snapshot/baseline_engine/worker.py').write_text('changed\n')
    with pytest.raises(ValidationError,match='snapshot differs'):
        verify_execution_snapshot(run,{'engine_sha256':expected})


def test_campaign_input_hash_ignores_snapshot_directory_but_detects_input_drift(tmp_path):
    from protac_qc import preparation_input_hash
    (tmp_path/'input.smi').write_text('CC\n');before=preparation_input_hash(tmp_path)
    (tmp_path/'source_snapshot').mkdir();(tmp_path/'source_snapshot/source.py').write_text('code\n')
    assert preparation_input_hash(tmp_path)==before
    (tmp_path/'input.smi').write_text('CO\n')
    assert preparation_input_hash(tmp_path)!=before


def test_campaign_keeps_one_code_version_when_working_tree_changes(tmp_path,monkeypatch):
    from types import SimpleNamespace
    from protac_qc import preparation_input_hash
    campaign=load('campaign_frozen_batch','scripts/lib/campaign.py')
    root=tmp_path/'project';engine=root/'baseline_engine';engine.mkdir(parents=True)
    (engine/'worker.py').write_text('original worker\n')
    for relative in ['scripts/lib/protac_qc.py','scripts/lib/campaign.py','scripts/01_setup/02_batch_prosettac.py','scripts/01_setup/03_record_attempt.py']:
        p=root/relative;p.parent.mkdir(parents=True,exist_ok=True);p.write_text('fixture\n')
    pairs=[]
    for name in ['one','two']:
        folder=tmp_path/name;folder.mkdir();(folder/'selection.json').write_text('{}\n')
        pairs.append({'protac':'GSK215','run_name':name,'gene':'PTK2','uid':'Q05397',
                      'local_dir':str(folder),'input_hash':preparation_input_hash(folder)})
    args=SimpleNamespace(python=Path(sys.executable),seed=1,protocol='global',n_dock=1,n_local=1,
        conformers=1,workers=1,max_global_trials=1,timeout_seconds=10,runs=tmp_path/'runs',rerun=False)
    monkeypatch.setattr(campaign,'validate_plan',lambda *a:{'mode':'TEST','minimum_free_bytes':0})
    monkeypatch.setattr(campaign.shutil,'disk_usage',lambda *a:SimpleNamespace(free=10_000_000_000))
    observed=[]
    def recorder(command,**kwargs):
        attempt=Path(command[command.index('--run-dir')+1])
        observed.append((attempt/'source_snapshot/baseline_engine/worker.py').read_text())
        (engine/'worker.py').write_text('changed working tree\n')
        atomic_json(attempt/'completion.json',{'status':'COMPLETE_NO_VALID_POSES'})
        return SimpleNamespace(returncode=0,stdout='',stderr='')
    monkeypatch.setattr(campaign.subprocess,'run',recorder)
    campaign.execute(root,pairs,args)
    assert observed==['original worker\n','original worker\n']
    assert len(list((args.runs/'source_snapshots').iterdir()))==1


def test_resource_plan_distinguishes_requested_and_generated_conformers(tmp_path):
    resources=load('resource_budget_evidence','scripts/01_setup/02a_resource_plan.py')
    inputs=tmp_path/'inputs';pair=inputs/'SJF1521/EGFR';pair.mkdir(parents=True)
    (pair/'preparation.json').write_text('{}\n')
    pilot=tmp_path/'pilot';raw=pilot/'Patchdock_Results';raw.mkdir(parents=True)
    (raw/'pd.1_docking_0001.pdb').write_text('synthetic resource fixture\n')
    ensemble=raw/'confs_1_1.sdf';ensemble.write_text('synthetic ensemble bytes\n')
    atomic_json(raw/'confs_1_1.sdf.sampling.json',{'requested_conformers':100,'generated_conformers':74})
    atomic_json(raw/'PT_1_1_parameter_provenance.json',{'input_conformers':74,'source_sha256':sha256(ensemble)})
    report=resources.plan(inputs,pilot,200,50,100,available_bytes=10**12)
    assert report['requested_full_conformers']==100 and report['actual_generated_conformers_by_local_pose']==[74]
    with pytest.raises(ValidationError,match='requested conformer budget'):
        resources.plan(inputs,pilot,200,50,200,available_bytes=10**12)
    (raw/'pd.1_docking_0002.pdb').write_text('synthetic rejected local pose\n')
    atomic_json(raw/'confs_1_2.sdf.sampling.json',{'requested_conformers':100,'generated_conformers':0})
    completed=resources.plan(inputs,pilot,200,50,100,available_bytes=10**12)
    assert sorted(completed['actual_generated_conformers_by_local_pose'])==[0,74]
    assert completed['zero_conformer_files_absent_after_upstream_empty_file_cleanup']==['confs_1_2.sdf']
    atomic_json(raw/'PT_1_1_parameter_provenance.json',{'input_conformers':75,'source_sha256':sha256(ensemble)})
    with pytest.raises(ValidationError,match='Parameter ensemble'):
        resources.plan(inputs,pilot,200,50,100,available_bytes=10**12)


def test_full_complex_qc_rejects_elongin_clashes_outside_vhl_binding_chain(native,tmp_path):
    from protac_qc import pdb_atoms
    selection=json.loads((native/'selection.json').read_text());selection['qc']={'min_nonbonded_heavy_distance_A':1.5}
    original=native/'native_complex.pdb';good=pose_features(original,selection)
    assert not good['severe_overlap'] and good['protein_clash_e3_chains']==['A','B','C']
    atoms=pdb_atoms(original)
    ligand=next(a for a in atoms if a['record']=='HETATM' and a['resname']=='PTC' and a['element'] not in {'H','D'})
    def modified(target,source,name):
        lines=original.read_text().splitlines();index=lines.index(target['line']);xyz=source['xyz']
        lines[index]=target['line'][:30]+''.join(f'{v:8.3f}' for v in xyz)+target['line'][54:]
        result=tmp_path/name;result.write_text('\n'.join(lines)+'\n');return pose_features(result,selection)
    elongin=next(a for a in atoms if a['record']=='ATOM' and a['chain']=='B' and a['name']=='CA')
    clash=modified(elongin,ligand,'ligand_elongin_clash.pdb')
    assert clash['d_vhl']>=1.5 and clash['d_kin']>=1.5
    assert clash['min_ligand_all_protein_A']==0 and clash['severe_overlap']
    poi=next(a for a in atoms if a['record']=='ATOM' and a['chain']=='D' and a['name']=='CA')
    clash=modified(elongin,poi,'protein_protein_clash.pdb')
    assert clash['min_e3_complex_poi_A']==0 and clash['severe_overlap']


def test_intraligand_clash_check_excludes_bonds_and_angles():
    from rdkit import Chem
    from protac_qc import check_ligand_nonbonded_geometry
    molecule=Chem.MolFromSmiles('CCCC');conf=Chem.Conformer(4)
    for i in range(4):conf.SetAtomPosition(i,[1.5*i,0.,0.])
    molecule.AddConformer(conf)
    assert check_ligand_nonbonded_geometry(molecule)==4.5
    molecule.GetConformer().SetAtomPosition(3,[0.,0.,0.])
    with pytest.raises(ValidationError,match='intraligand'):
        check_ligand_nonbonded_geometry(molecule)


def test_ce_alignment_compatibility_preserves_old_and_new_api_correspondences():
    from collections import namedtuple
    modern=namedtuple('ModernCE','length path')(3,([1,2,3],[4,5,6]))
    legacy=[[1,2,3],[4,5,6]]
    assert placement.ce_path_indices(modern)==placement.ce_path_indices(legacy)==((1,2,3),(4,5,6))
    with pytest.raises(ValidationError,match='Malformed'):
        placement.ce_path_indices([[1,2],[3]])


def test_distance_sampler_requires_exact_unambiguous_stereo_anchors():
    from rdkit import Chem
    helper=load('distance_head_mapping','baseline_engine/constrained_linker.py')
    folder=ROOT/'prepared_inputs/SJF8240/MET_P08581'
    heads=[Chem.SDMolSupplier(str(folder/n))[0] for n in ['vhl_head.sdf','poi_head.sdf']]
    full=Chem.MolFromSmiles((folder/'full.smi').read_text().split()[0])
    anchors=[int(x)-1 for x in read_params(folder/'Protac_params.txt')['Anchor atoms'].split()]
    assert helper.distance_anchor_indices(full,heads,anchors)==(49,36)
    with pytest.raises(ValidationError,match='outside'):
        helper.distance_anchor_indices(full,heads,[999,anchors[1]])
    with pytest.raises(ValidationError,match='Ambiguous'):
        helper.distance_anchor_indices(Chem.MolFromSmiles('CC.CC'),[Chem.MolFromSmiles('C')]*2,[0,0])
    for atom in full.GetAtoms():
        if atom.GetChiralTag()!=Chem.ChiralType.CHI_UNSPECIFIED:
            atom.InvertChirality();break
    with pytest.raises(ValidationError,match='stereochemical'):
        helper.distance_anchor_indices(full,heads,anchors)


def test_initial_distance_sampler_records_exact_mapping_and_single_thread(tmp_path,monkeypatch):
    from rdkit import Chem
    monkeypatch.syspath_prepend(str(ROOT/'baseline_engine'))
    monkeypatch.setenv('PROTAC_ATTEMPT_ROOT',str(tmp_path));monkeypatch.setenv('PROTAC_SEED','155921')
    monkeypatch.chdir(tmp_path)
    adapter=load('distance_sampler_adapter','baseline_engine/pyrosetta_adapter.py')
    folder=ROOT/'prepared_inputs/SJF8240/MET_P08581'
    params=read_params(folder/'Protac_params.txt')
    heads=[str(folder/n) for n in ['vhl_head.sdf','poi_head.sdf']]
    anchors=[int(x)-1 for x in params['Anchor atoms'].split()]
    bounds=adapter.sample_dist_free_conf(heads,anchors,str(folder/'full.smi'),n=3)
    ledger=json.loads((tmp_path/'anchor_distance_sampling.json').read_text())
    assert ledger['num_threads']==1 and ledger['requested_conformers']==3
    assert 0<ledger['generated_conformers']<=3
    assert ledger['full_molecule_anchor_indices_zero_based']==[49,36]
    assert list(bounds)==ledger['range_A'] and bounds[0]<bounds[1]
    assert sum(ledger['histogram'].values())==ledger['generated_conformers']


def test_gnn_contract_excludes_unknowns_and_checks_cohort_labels():
    from types import SimpleNamespace
    from gnn_contract import select_training_graphs
    graph=lambda p,y:SimpleNamespace(poi=p,y=np.array([y]),artifact=False)
    kept,excluded=select_training_graphs([graph('known',1),graph('prediction',-1)],{'known':1})
    assert [x.poi for x in kept]==['known'] and excluded==['prediction']
    with pytest.raises(ValidationError,match='label'):
        select_training_graphs([graph('known',-1)],{'known':1})
    with pytest.raises(ValidationError,match='Duplicate'):
        select_training_graphs([graph('known',1),graph('known',1)],{'known':1})


def test_gnn_scalar_normalization_uses_training_rows_only():
    import pandas as pd
    from gnn_contract import normalize_fold_scalars
    t=pd.DataFrame({'x':[1.,3.,100.]},index=['a','b','heldout'])
    assert normalize_fold_scalars('heldout',t,['a','b'],['x'])[0]==pytest.approx(98/np.sqrt(2))
    t.loc['heldout','x']=np.inf
    with pytest.raises(ValidationError,match='Infinite'):
        normalize_fold_scalars('heldout',t,['a','b'],['x'])


def test_gnn_periodic_contacts_and_joint_frequency_reject_independence_proxy():
    from gnn_contract import periodic_minima,observed_joint_contacts
    d=periodic_minima([[.2,0,0]],[[9.8,0,0]],[10,10,10,90,90,90])
    assert d[0]==pytest.approx(.4,abs=1e-5)
    left=[[True],[False]];right=[[False],[True]]
    assert observed_joint_contacts(left,right)[0,0]==0 # Marginal product would be 0.25.
    with pytest.raises(ValidationError,match='periodic box'):
        periodic_minima([[0,0,0]],[[1,1,1]],None)


def test_gnn_release_rejects_historical_chemistry_before_model_loading(tmp_path):
    from gnn_contract import verify_gnn_release
    release=tmp_path/'release.json';release.write_text(json.dumps({'status':'RELEASED_FOR_EXPLORATORY_GNN','compound':'SJF8240','ligand_inchikey':'MCJFROJTDORRBB-ZWUDBOSPSA-N'}))
    with pytest.raises(ValidationError,match='chemistry'):
        verify_gnn_release(release,{},'SCUZPTNDZHEHLB-ZWUDBOSPSA-N')


def test_compiled_ligand_guard_checks_connectivity_beyond_formula():
    import MDAnalysis as mda
    from rdkit import Chem
    from gnn_contract import require_verified_ligand_graph
    smiles=json.loads((ROOT/'config/molecules.json').read_text())['compounds']['SJF8240']['smiles']
    mol=Chem.MolFromSmiles(smiles);u=mda.Universe.empty(mol.GetNumAtoms(),n_residues=1,atom_resindex=np.zeros(mol.GetNumAtoms(),dtype=int))
    u.add_TopologyAttr('resnames',['PTC']);u.add_TopologyAttr('names',[f'A{i}' for i in range(mol.GetNumAtoms())]);u.add_TopologyAttr('masses',[a.GetMass() for a in mol.GetAtoms()])
    bonds=[(b.GetBeginAtomIdx(),b.GetEndAtomIdx()) for b in mol.GetBonds()];u.add_TopologyAttr('bonds',bonds)
    assert len(require_verified_ligand_graph(u,smiles))==79
    u.delete_bonds([bonds[0]]) # Atom identities/counts/formula are unchanged.
    with pytest.raises(ValidationError,match='connectivity'):
        require_verified_ligand_graph(u,smiles)


def test_static_graph_uses_conditional_distances_and_frozen_node_correspondence(tmp_path,monkeypatch):
    builder=load('static_sjf_graph','scripts/04_gnn/03_static_contact_graph.py')
    run=tmp_path/'run';poses=run/'Patchdock_Results';poses.mkdir(parents=True)
    smiles=json.loads((ROOT/'config/molecules.json').read_text())['compounds']['SJF8240']['smiles']
    (run/'full.smi').write_text(smiles+'\n');(run/'completion.json').write_text('{}')
    selection={'vhl':{'chains':['A']},'poi':{'chains':['D']},'protac':{'chains':['X'],'resname':'PTC'}}
    (run/'selection.json').write_text(json.dumps(selection))
    def atom(record,i,name,res,chain,x):
        return f'{record:<6}{i:5d} {name:^4s} {res:3s} {chain}   1    {x:8.3f}{0.:8.3f}{0.:8.3f}  1.00  0.00           C  '
    for pose,offset in [('one',1.),('two',3.)]:
        lines=[atom('ATOM',1,'CA','ALA','A',0.),atom('ATOM',2,'CA','ALA','D',10.)]
        lines.extend(atom('HETATM',i+3,f'C{i}','PTC','X',offset+i*.01) for i in range(79))
        (poses/(pose+'.pdb')).write_text('\n'.join(lines)+'\n')
    # Unit-test graph mathematics after the independently tested QC boundary.
    monkeypatch.setattr(builder,'verify_receipt',lambda *a:[{'description':'one'},{'description':'two'}])
    out=tmp_path/'graph';builder.build(run,'synthetic',-1,out)
    with np.load(out/'graph.npz',allow_pickle=False) as a:
        assert a['x'].shape==(81,9)
        def edge(i,j):return a['edge_attr'][np.where((a['edge_index'][0]==i)&(a['edge_index'][1]==j))[0][0]]
        assert edge(2,0)==pytest.approx([1.,.2,.1])
        assert edge(2,1)==pytest.approx([.5,.7,.7]) # Mean is conditional on contact, not 0.35 nm.


def test_gnn_release_prevents_aliases_of_same_protein_crossing_folds(tmp_path):
    from gnn_contract import verify_gnn_release
    p=tmp_path/'release.json';p.write_text(json.dumps({'status':'RELEASED_FOR_EXPLORATORY_GNN','compound':'SJF8240','ligand_inchikey':'SCUZPTNDZHEHLB-ZWUDBOSPSA-N','feature_schema':'SJF8240-periodic-contacts-v1','experimental_labels':{'MAPK14':1,'MAPK14_full':0},'protein_groups':{'MAPK14':'Q16539','MAPK14_full':'Q16539'}}))
    with pytest.raises(ValidationError,match='multiple times'):
        verify_gnn_release(p,{},'SCUZPTNDZHEHLB-ZWUDBOSPSA-N')


def test_md_graph_builder_uses_observed_joint_contacts_and_aligned_nodes(monkeypatch,tmp_path):
    module=load('periodic_sjf_md_graph','scripts/04_gnn/02_md_contact_graphs.py')
    monkeypatch.setattr(module,'HAS_TORCH',False)
    result={'vhl_resids':np.array([1]),'kin_resids':np.array([200]),'vhl_freq':np.array([.5]),'kin_freq':np.array([.5]),'vhl_mean_d':np.array([4.]),'kin_mean_d':np.array([4.]),'joint_freq':np.array([[0.]]),'n_frames':2,'n_ptc_atoms':79,'input_sha256':{}}
    monkeypatch.setattr(module,'compute_pose_contacts',lambda *a:result)
    cfg={'kin_min':0,'degrader':None,'screen_dir':tmp_path};poses=[{'pose_id':'one'},{'pose_id':'two'}]
    graph,metadata=module.build_poi_md_graph('synthetic',poses,1,cfg)
    assert graph['degrader']==-1 and graph['edge_src'].size==0 and graph['edge_attr'].shape==(0,3)
    assert len(metadata)==2 # Product of the marginals would incorrectly create an edge.
    count=[0]
    def changed(*args):
        count[0]+=1
        return dict(result,kin_resids=np.array([200 if count[0]==1 else 201]))
    monkeypatch.setattr(module,'compute_pose_contacts',changed)
    with pytest.raises(ValidationError,match='Unaligned'):
        module.build_poi_md_graph('synthetic',poses,1,cfg)


def test_gnn_graph_schema_rejects_legacy_packets_even_when_dimensions_match():
    from types import SimpleNamespace
    from gnn_contract import validate_feature_graph
    graph=SimpleNamespace(x=np.ones((2,9)),edge_index=np.empty((2,0),dtype=int),edge_attr=np.empty((0,3)))
    with pytest.raises(ValidationError,match='schema'):
        validate_feature_graph(graph,9)
    graph.feature_schema='SJF8240-periodic-contacts-v1';validate_feature_graph(graph,9)
    graph.edge_index=np.array([[0],[9]]);graph.edge_attr=np.ones((1,3))
    with pytest.raises(ValidationError,match='index'):
        validate_feature_graph(graph,9)


def test_gnn_runner_rejects_historical_scalar_table_without_versioned_chemistry(tmp_path):
    module=load('isolated_sjf_e2_runner','scripts/04_gnn/01_sjf8240_e2_lopo.py')
    p=tmp_path/'old_table.tsv';p.write_text('poi\tprostac_artifact\nMET\t0\n')
    with pytest.raises(ValidationError,match='scalar feature schema'):
        module.load_e2_table(p)


def test_sjf8240_sequence_repair_inputs_are_portable_and_hash_bound():
    inventory=load('portable_sjf_sequences','scripts/01_setup/00i_build_domain_inventory.py')
    cache=ROOT/'sequence_provenance';source=json.loads((cache/'SOURCE.json').read_text());p=cache/'uniprot_batches/sjf8240_reference_cohort.json'
    assert sha256(p)==source['selected_cache_sha256']
    entries=inventory.entries(cache);cohort=json.loads((ROOT/'config/sjf8240_publication_cohort.json').read_text())
    assert set(cohort['protein_groups'].values())|{'P40337','Q15370','Q15369'}==set(entries)
    assert len(entries)==20 and all(x['sequence']['value'] for x in entries.values())


def test_md_preparation_rejects_repaired_protein_from_a_different_pose(tmp_path):
    md=load('md_pose_repair_binding','scripts/01_setup/06_prepare_md.py')
    pose=tmp_path/'pose.pdb';protein=tmp_path/'protein.pdb';other=tmp_path/'other_pose.pdb'
    pose.write_text('Original observed protein/ligand frame\n');protein.write_text('Repaired protein in that frame\n');other.write_text('Different docking arrangement\n')
    repair={'status':'PROTEIN_MODEL_COMPLETED','source_sha256':sha256(pose),'output_sha256':sha256(protein)}
    md.verify_protein_repair(pose,protein,repair)
    with pytest.raises(ValidationError,match='different starting pose'):
        md.verify_protein_repair(other,protein,repair)


def test_distributed_parameter_bundles_contain_every_receipt_bound_artifact():
    reports=sorted((ROOT/'ligand_params').rglob('parameter_validation.json'))
    assert len(reports)>=6
    for path in reports:
        report=json.loads(path.read_text())
        assert report['status']=='PARAMETERS_VERIFIED'
        for name,digest in report['files'].items():
            artifact=path.parent/name
            assert artifact.is_file(),str(artifact)
            assert sha256(artifact)==digest,str(artifact)


def test_sequence_repair_rejects_new_loop_overlap_with_another_chain():
    repair=load('protein_interchain_geometry','scripts/01_setup/05_complete_protein.py')
    atoms=[{'chain':'A','element':'C','xyz':np.array([0.,0.,0.])},
           {'chain':'D','element':'C','xyz':np.array([2.5,0.,0.])}]
    assert repair.validate_interchain_contacts(atoms)=={'A:D':2.5}
    atoms.append({'chain':'D','element':'C','xyz':np.array([0.5,0.,0.])})
    with pytest.raises(ValidationError,match='interchain clash'):
        repair.validate_interchain_contacts(atoms)


def test_attempt_refuses_incomplete_source_snapshot_before_starting_worker(tmp_path):
    recorder=load('preexecution_snapshot_guard','scripts/01_setup/03_record_attempt.py')
    (tmp_path/'selection.json').write_text('{}')
    (tmp_path/'source_snapshot').mkdir()
    marker=tmp_path/'worker_started'
    with pytest.raises(ValidationError,match='source manifest required'):
        recorder.record(tmp_path,[sys.executable,'-c',"from pathlib import Path; Path('worker_started').touch()"],10)
    assert not marker.exists() and not (tmp_path/'attempt.json').exists()


@pytest.mark.parametrize('translation',[[0.,0.,0.],[1234.,-987.,456.]])
def test_default_linker_preserves_real_sjf_bound_heads_and_coordinate_frame(translation):
    linker=load('default_sjf_bound_head_regression','baseline_engine/constrained_linker.py')
    from rdkit import Chem
    fixture=ROOT/'tests/fixtures/sjf8240_bound_heads'
    full=Chem.MolFromSmiles((fixture/'full.smi').read_text().split()[0])
    heads=[Chem.SDMolSupplier(str(fixture/name),removeHs=True)[0] for name in ['vhl_head_H.sdf','poi_head_H.sdf']]
    xyz=np.asarray(json.loads((fixture/'head_coordinates.json').read_text())['10_1']['coordinates_A'])
    targets=[xyz[:30]+translation,xyz[30:]+translation]
    mol,report=linker.constrained_linker(full,heads,targets,155921,attempts=3)
    assert mol is not None and report['embedding_method']=='CENTERED_TORSION_DG_FIXED_HEADS_UFF'
    assert linker.stereo_key(mol)==Chem.MolToInchiKey(full)
    for indices,expected in zip(linker.match_heads(full,heads)[0],targets):
        assert np.max(np.linalg.norm(mol.GetConformer().GetPositions()[list(indices)]-expected,axis=1))<1e-5


def test_md_graph_producer_has_no_implicit_original_paths_or_biological_labels():
    graph=load('no_legacy_md_graph_defaults','scripts/04_gnn/02_md_contact_graphs.py')
    assert graph.POI_CONFIG=={} and graph.RANDOM_SCREEN is None
    with pytest.raises(ValidationError,match='Explicit rebuilt screen directory'):
        graph.build_poi_md_graph('MET',[],1,{'kin_min':0,'degrader':None})


def test_gromacs_energy_columns_follow_legends_instead_of_selection_order():
    energy=load('gmx_named_energy','scripts/lib/gromacs_energy.py')
    text='@ s0 legend "Potential"\n@ s1 legend "Temperature"\n0 -1000 299\n2 -990 301\n'
    time,columns=energy.named_energy_series(text,['Temperature','Potential'])
    assert np.array_equal(columns['Temperature'],[299,301])
    assert np.array_equal(columns['Potential'],[-1000,-990])
    assert np.array_equal(time,[0,2])


def test_gromacs_energy_refuses_unlabelled_or_nonphysical_temperature_columns():
    energy=load('gmx_energy_guard','scripts/lib/gromacs_energy.py')
    with pytest.raises(ValidationError,match='explicit indexed legend'):
        energy.named_energy_series('0 -1000 299\n2 -990 301\n',['Temperature'])
    with pytest.raises(ValidationError,match='Nonpositive Temperature'):
        energy.named_energy_series('@ s0 legend "Temperature"\n0 -1000\n2 -990\n',['Temperature'])


def test_existing_bridge_audit_uses_joint_frames_and_strict_historical_cutoff():
    audit=load('legacy_bridge_joint','scripts/02_analysis/07_audit_existing_bridge_metrics.py')
    time=np.array([0.,1.]);vhl=np.array([4.,6.]);poi=np.array([6.,4.])
    row=audit.joint_statistics(time,vhl,time,poi,5.)
    assert row['same_frame_bridge_fraction']==0. and row['product_of_marginals']==0.25
    boundary=audit.joint_statistics(time,np.array([5.,4.]),time,np.array([4.,4.]),5.)
    assert boundary['same_frame_bridge_fraction']==0.5 and boundary['exact_threshold_vhl_frames']==1


def test_existing_bridge_audit_refuses_timestamp_misalignment_or_truncation():
    audit=load('legacy_bridge_time_guard','scripts/02_analysis/07_audit_existing_bridge_metrics.py')
    with pytest.raises(ValidationError,match='timestamps'):
        audit.joint_statistics([0.,1.],[4.,4.],[0.,2.],[4.,4.],5.)
    with pytest.raises(ValidationError,match='never truncate'):
        audit.joint_statistics([0.,1.],[4.,4.],[0.],[4.],5.)


def test_existing_bridge_distance_reader_requires_explicit_units(tmp_path):
    audit=load('legacy_bridge_units','scripts/02_analysis/07_audit_existing_bridge_metrics.py')
    path=tmp_path/'distance.xvg'
    path.write_text('@ xaxis label "Time (ns)"\n@ yaxis label "Distance (nm)"\n@ s0 legend "Other-VHL"\n0 0.4\n1 0.5\n')
    _,distance=audit.distance_series(path,'Other-VHL');assert np.array_equal(distance,[4.,5.])
    path.write_text(path.read_text().replace('Distance (nm)','Distance (A)'))
    with pytest.raises(ValidationError,match='Explicit ns/nm'):
        audit.distance_series(path,'Other-VHL')


def test_existing_bridge_duplicate_repair_requires_explicit_identical_samples(tmp_path):
    audit=load('legacy_exact_duplicate_policy','scripts/02_analysis/07_audit_existing_bridge_metrics.py')
    path=tmp_path/'distance.xvg'
    header='@ xaxis label "Time (ns)"\n@ yaxis label "Distance (nm)"\n@ s0 legend "Other-VHL"\n'
    path.write_text(header+'0 0.4\n1 0.5\n1 0.5\n2 0.6\n')
    with pytest.raises(ValidationError,match='time is not increasing'):
        audit.distance_series(path,'Other-VHL')
    time,distance=audit.distance_series(path,'Other-VHL',allow_identical_duplicates=True)
    assert np.array_equal(time,[0,1,2]) and np.array_equal(distance,[4,5,6])
    path.write_text(header+'0 0.4\n1 0.5\n1 0.7\n2 0.6\n')
    with pytest.raises(ValidationError,match='Conflicting duplicate'):
        audit.distance_series(path,'Other-VHL',allow_identical_duplicates=True)


@pytest.mark.parametrize('refresh_offsets',[False,True])
def test_readonly_trajectory_reader_never_writes_offset_or_lock_sidecars(tmp_path,refresh_offsets):
    mda=pytest.importorskip('MDAnalysis')
    module=load('read_only_xtc','scripts/lib/read_only_trajectory.py')
    universe=mda.Universe.empty(2,trajectory=True);path=tmp_path/'existing.xtc'
    with mda.Writer(str(path),n_atoms=2) as writer:
        for frame in range(2):
            universe.atoms.positions=np.array([[frame,0.,0.],[frame+1.,0.,0.]],dtype=np.float32)
            universe.trajectory.ts.time=float(frame);universe.dimensions=[10.,10.,10.,90.,90.,90.];writer.write(universe.atoms)
    before={p.name for p in tmp_path.iterdir()};digest=sha256(path)
    with module.ReadOnlyXTCReader(str(path),refresh_offsets=refresh_offsets) as reader:
        assert len(reader)==2 and reader[1].time==1.
    assert {p.name for p in tmp_path.iterdir()}==before and sha256(path)==digest


def test_gnn_fold_checkpoint_replays_score_without_retraining(tmp_path):
    from types import SimpleNamespace
    module=load('fold_model_replay','scripts/04_gnn/01_sjf8240_e2_lopo.py')
    if not module.HAS_DEPS:pytest.skip('Requires the analysis GNN environment; docking environment has no PyG')
    import torch,pandas as pd
    from torch_geometric.data import Data
    torch.manual_seed(7919);torch.set_num_threads(1)
    args=SimpleNamespace(hidden=8,heads=2,dropout=0.)
    model=module.E2DualEncoderGNN(hidden=8,heads=2,dropout=0.).eval()
    graph=Data(x=torch.ones((3,9)),edge_index=torch.tensor([[0,1,1,2],[1,0,2,1]],dtype=torch.long),edge_attr=torch.ones((4,3)))
    table=pd.DataFrame([[1.,2.,3.,4.],[2.,3.,4.,5.]],index=['TRAIN_A','TRAIN_B'],columns=module.E2_COLS)
    e2=torch.zeros((1,4));score=module.predict(model,graph,None,e2,'cpu');path=tmp_path/'fold.pt'
    artifact=module.save_fold_checkpoint(path,model,'HELD_OUT',['TRAIN_A','TRAIN_B'],7919,table,e2,args,score)
    assert artifact['sha256']==sha256(path)
    saved=torch.load(path,map_location='cpu',weights_only=True)
    restored=module.E2DualEncoderGNN(**saved['architecture']);restored.load_state_dict(saved['model_state_dict'])
    assert module.predict(restored,graph,None,saved['held_out_normalized_scalars'],'cpu')==score
    assert saved['held_out_score']==score and saved['held_out_poi'] not in saved['training_pois']
    assert torch.equal(saved['scalar_reference_values'],torch.tensor(table.to_numpy(),dtype=torch.float64))
    with pytest.raises(ValidationError,match='already exists'):
        module.save_fold_checkpoint(path,model,'HELD_OUT',['TRAIN_A'],7919,table,e2,args,score)


def test_gnn_fold_checkpoint_rejects_held_out_target_in_training_metadata(tmp_path):
    from types import SimpleNamespace
    module=load('fold_checkpoint_leak_guard','scripts/04_gnn/01_sjf8240_e2_lopo.py')
    if not module.HAS_DEPS:pytest.skip('Requires the analysis GNN environment; docking environment has no PyG')
    import torch,pandas as pd
    args=SimpleNamespace(hidden=8,heads=2,dropout=0.)
    model=module.E2DualEncoderGNN(hidden=8,heads=2,dropout=0.)
    table=pd.DataFrame([[1.,2.,3.,4.]],index=['HELD_OUT'],columns=module.E2_COLS)
    with pytest.raises(ValidationError,match='leaked'):
        module.save_fold_checkpoint(tmp_path/'fold.pt',model,'HELD_OUT',['HELD_OUT'],1,table,torch.zeros((1,4)),args,.5)
    assert not (tmp_path/'fold.pt').exists()


def test_tiny_synthetic_lopo_retains_full_precision_scores_and_each_fold_model(tmp_path):
    from types import SimpleNamespace
    module=load('synthetic_fold_export','scripts/04_gnn/01_sjf8240_e2_lopo.py')
    if not module.HAS_DEPS:pytest.skip('Requires the analysis GNN environment; docking environment has no PyG')
    import torch,pandas as pd
    from torch_geometric.data import Data
    args=SimpleNamespace(seed=17,epochs=1,hidden=8,heads=2,dropout=0.,lr=.001,verbose=False,output=tmp_path)
    graphs=[Data(x=torch.ones((3,9)),edge_index=torch.tensor([[0,1,1,2],[1,0,2,1]],dtype=torch.long),edge_attr=torch.ones((4,3)),y=torch.tensor([i%2]),poi='SYNTHETIC_'+str(i)) for i in range(4)]
    table=pd.DataFrame([[i,i+1,i+2,i+3] for i in range(4)],index=[g.poi for g in graphs],columns=module.E2_COLS)
    predictions,_=module.run_lopo(graphs,{},table,args)
    assert len(predictions)==4 and len(list((tmp_path/'fold_checkpoints').glob('*.pt')))==4
    for row in predictions.to_dict('records'):
        path=tmp_path/row['fold_checkpoint'];assert sha256(path)==row['fold_checkpoint_sha256']
        saved=torch.load(path,map_location='cpu',weights_only=True)
        assert saved['held_out_score']==row['prob_degrader'] and saved['held_out_poi']==row['poi']
        assert row['poi'] not in saved['training_pois']
