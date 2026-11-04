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
