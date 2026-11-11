"""Shared strict file contracts for preparation, execution and analysis."""
import hashlib
import json
import math
import os
import tempfile
from pathlib import Path
import numpy as np

class ValidationError(ValueError):
    pass

def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()

def engine_source_hashes(root):
    root=Path(root)
    files=[f for f in (root/'baseline_engine').rglob('*') if f.is_file() and '__pycache__' not in f.parts
           and f.suffix in {'.py','.json','.smarts','.xml','.flags','.params'}]
    files.extend(root/relative for relative in ['scripts/lib/protac_qc.py','scripts/lib/campaign.py',
        'scripts/01_setup/02_batch_prosettac.py','scripts/01_setup/03_record_attempt.py'])
    return {str(f.relative_to(root)):sha256(f) for f in sorted(files)}

def preparation_input_hash(folder):
    files=sorted(p for p in Path(folder).iterdir() if p.is_file())
    return hashlib.sha256(json.dumps([(p.name,sha256(p)) for p in files]).encode()).hexdigest()

def verify_execution_snapshot(run,report):
    run=Path(run);snapshot=run/'source_snapshot'
    if not snapshot.exists():return  # Historical non-snapshot diagnostic receipts.
    manifest=json.loads((snapshot/'source_manifest.json').read_text())
    if manifest['source_sha256']!=report['engine_sha256'] or engine_source_hashes(snapshot)!=manifest['source_sha256']:
        raise ValidationError('Execution source snapshot differs from completed engine version')
    intent=run/'campaign_intent.json'
    if intent.exists() and json.loads(intent.read_text())['source_snapshot_manifest_sha256']!=sha256(snapshot/'source_manifest.json'):
        raise ValidationError('Execution source manifest differs from attempt intent')

def atomic_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix='.' + path.name, dir=path.parent)
    try:
        with os.fdopen(fd, 'w') as f:
            json.dump(data, f, indent=2, allow_nan=False)
            f.write('\n'); f.flush(); os.fsync(f.fileno())
        os.replace(name, path)
    finally:
        if os.path.exists(name): os.unlink(name)

def read_params(path):
    result = {}
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.lstrip().startswith('#'): continue
        if ':' not in line: raise ValidationError('Malformed params line')
        key, value = (x.strip() for x in line.split(':', 1))
        if key in result: raise ValidationError('Duplicate parameter: ' + key)
        result[key] = value
    required = {'Structures', 'Chains', 'Heads', 'Anchor atoms', 'Protac', 'Full', 'ClusterName'}
    if required - result.keys():
        raise ValidationError('Missing runner keys: ' + ', '.join(sorted(required - result.keys())))
    if result['Full'] not in {'True', 'False'} or result['ClusterName'] != 'LOCAL':
        raise ValidationError('Explicit Full boolean and ClusterName: LOCAL required')
    for key in ['Structures', 'Chains', 'Heads', 'Anchor atoms']:
        if len(result[key].split()) != 2: raise ValidationError(key + ' requires two entries')
    try:
        if min(map(int, result['Anchor atoms'].split())) < 1: raise ValueError()
    except ValueError as exc:
        raise ValidationError('Anchor atoms must be positive one-based indices') from exc
    return result

def pdb_atoms(path):
    """Read one model, explicit single-conformer atoms, no guessed chain identity."""
    atoms = []
    seen = set()
    for line in Path(path).read_text().splitlines():
        if line.startswith('ENDMDL'): break
        if not line.startswith(('ATOM  ', 'HETATM')): continue
        if line[16:17] not in (' ', 'A'): continue
        key = (line[21], line[22:27], line[12:16])
        if key in seen: raise ValidationError('Duplicate atom/alternate conformer: ' + str(key))
        seen.add(key)
        try:
            pos = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            if not np.isfinite(pos).all(): raise ValueError()
            resnum = int(line[22:26])
        except ValueError as exc: raise ValidationError('Invalid PDB coordinates/residue') from exc
        element = line[76:78].strip().upper() or line[12:16].strip().lstrip('0123456789')[0]
        atoms.append({'record': line[:6].strip(), 'chain': line[21], 'resnum': resnum,
                      'icode': line[26:27], 'resname': line[17:20].strip(),
                      'name': line[12:16].strip(), 'element': element, 'xyz': pos, 'line': line})
    if not atoms: raise ValidationError('No atoms: ' + str(path))
    return atoms

def selected_xyz(atoms, selector):
    # Residue restrictions are optional and explicit, never inferred from numbering.
    required = {'chains', 'record'}
    if not required.issubset(selector) or not selector['chains']:
        raise ValidationError('Atom selector needs chains and record')
    selected = [a['xyz'] for a in atoms if a['chain'] in selector['chains']
                and a['record'] == selector['record'] and a['element'] not in {'H','D'}
                and ('resname' not in selector or a['resname'] == selector['resname'])
                and ('resnums' not in selector or a['resnum'] in selector['resnums'])]
    if not selected: raise ValidationError('Empty atom selection: ' + str(selector))
    return np.asarray(selected)

def min_distance(a, b):
    from scipy.spatial import cKDTree
    return float(cKDTree(b).query(a)[0].min())

def parse_scores(path,allow_empty=False):
    """Rosetta SCORE records; ignore SEQUENCE and repeated identical headers."""
    header, rows, seen = None, [], set()
    for lineno, line in enumerate(Path(path).read_text().splitlines(), 1):
        fields = line.split()
        if not fields or fields[0] != 'SCORE:': continue
        fields = fields[1:]
        if 'description' in fields and 'total_score' in fields:
            if header is not None and fields != header: raise ValidationError('Conflicting score headers')
            header = fields; continue
        if header is None: raise ValidationError('Score row before header at line ' + str(lineno))
        if len(fields) != len(header): raise ValidationError('Truncated score row at line ' + str(lineno))
        row = dict(zip(header, fields)); pose = row['description']
        if Path(pose).name != pose or pose in seen: raise ValidationError('Unsafe/duplicate pose ID: ' + pose)
        try:
            score = float(row['total_score'])
            if not math.isfinite(score): raise ValueError()
        except ValueError as exc: raise ValidationError('Nonfinite/nonnumeric total_score') from exc
        row['total_score'] = score
        rows.append(row); seen.add(pose)
    if not rows and not (allow_empty and header is not None): raise ValidationError('No scored poses')
    return rows

def pose_features(pdb, selection):
    atoms = pdb_atoms(pdb)
    for k in ('vhl', 'poi', 'protac'):
        if k not in selection: raise ValidationError('Missing selection: ' + k)
    groups = {k: selected_xyz(atoms, selection[k]) for k in ('vhl','poi','protac')}
    # Overlapping selectors can create artificial zero-distance bridging.
    indices = {}
    for k in groups:
        sel = selection[k]
        indices[k] = {(a['chain'], a['resnum'], a['icode'], a['name']) for a in atoms
                      if a['chain'] in sel['chains'] and a['record'] == sel['record']
                      and ('resname' not in sel or a['resname'] == sel['resname'])
                      and ('resnums' not in sel or a['resnum'] in sel['resnums'])}
    if any(indices[a] & indices[b] for a,b in [('vhl','poi'),('vhl','protac'),('poi','protac')]):
        raise ValidationError('Overlapping atom selections')
    d_vhl = min_distance(groups['protac'], groups['vhl'])
    d_kin = min_distance(groups['protac'], groups['poi'])
    cutoff=float(selection.get('qc',{}).get('min_nonbonded_heavy_distance_A',1.0))
    if not math.isfinite(cutoff) or cutoff<=0:raise ValidationError('Invalid overlap cutoff')
    # Normalized ternaries contain protein ATOM records and an explicitly
    # selected ligand. Check the entire E3 complex, including elongins, rather
    # than only the VHL binding chain used by the contact feature.
    protein=[a for i,a in enumerate(atoms) if a['record']=='ATOM' and i not in indices['protac'] and a['element'] not in {'H','D'}]
    e3=[a for a in protein if a['chain'] not in selection['poi']['chains']]
    if not protein or not e3:raise ValidationError('Missing complete protein-complex selection')
    ligand_protein=min_distance(groups['protac'],np.array([a['xyz'] for a in protein]))
    e3_poi=min_distance(np.array([a['xyz'] for a in e3]),groups['poi'])
    return {'d_vhl': d_vhl, 'd_kin': d_kin, 'bridging_5A': d_vhl < 5 and d_kin < 5,
            'bridging_8A': d_vhl < 8 and d_kin < 8,
            'clash_cutoff_A':cutoff,'min_ligand_all_protein_A':ligand_protein,
            'min_e3_complex_poi_A':e3_poi,'protein_clash_e3_chains':sorted({a['chain'] for a in e3}),
            'severe_overlap': min(d_vhl,d_kin,ligand_protein,e3_poi) < cutoff}

def validate_outputs(run, selection):
    """Require every scored combined pose, fail on mismatches or sentinel scores."""
    run = Path(run); base = run / 'Patchdock_Results'
    score_path = base / 'score.sc'
    report=json.loads((run/'engine_status.json').read_text())
    no_valid=report.get('status')=='COMPLETE_NO_VALID_POSES'
    rows = parse_scores(score_path,allow_empty=no_valid)
    validate_stage_report(run,len(rows))
    if score_path.with_suffix('.sc.shift').exists() or (base/'score.sc.orig').exists():
        raise ValidationError('Shifted scores require explicit original-score reconciliation')
    expected = {r['description'] + '.pdb' for r in rows}
    actual = {p.name for p in base.glob('combined_*.pdb')}
    if expected != actual: raise ValidationError('Scored/generated combined pose mismatch')
    result = []
    for row in rows:
        if row['total_score'] == -1.0: raise ValidationError('Possible adapter placeholder score -1.0')
        path = base / (row['description'] + '.pdb')
        feats = pose_features(path, selection)
        validate_ligand_identity(path,selection['protac'],run/'full.smi')
        if feats['severe_overlap']: raise ValidationError('Severe ligand/protein overlap: ' + path.name)
        result.append(dict(row, **feats, pdb_sha256=sha256(path)))
    return result

def verify_receipt(run, selection):
    run = Path(run)
    receipt = json.loads((run/'completion.json').read_text())
    if receipt.get('status') not in {'COMPLETE','COMPLETE_NO_VALID_POSES'} or receipt.get('exit_code') != 0:
        raise ValidationError('No successful terminal receipt')
    if receipt.get('selection_sha256') != sha256(run/'selection.json'):
        raise ValidationError('Selection changed since completion')
    if not receipt.get('input_sha256'):
        raise ValidationError('Missing immutable input hashes')
    for name,digest in receipt['input_sha256'].items():
        if Path(name).name != name or sha256(run/name) != digest:
            raise ValidationError('Input changed since completion: ' + name)
    if json.loads((run/'selection.json').read_text()) != selection:
        raise ValidationError('Requested selection differs from receipt selection')
    rows = validate_outputs(run, selection)
    if receipt.get('engine_status_sha256')!=sha256(run/'engine_status.json'):
        raise ValidationError('Engine stage report changed since completion')
    verify_execution_snapshot(run,json.loads((run/'engine_status.json').read_text()))
    if receipt.get('score_sha256') != sha256(run/'Patchdock_Results/score.sc'):
        raise ValidationError('Scores changed since completion')
    hashes = {r['description']:r['pdb_sha256'] for r in rows}
    if hashes != receipt.get('pose_sha256'): raise ValidationError('Pose files changed since completion')
    return rows

def check_ligand_nonbonded_geometry(mol,cutoff=1.5):
    from rdkit import Chem
    topology=Chem.GetDistanceMatrix(mol);coordinates=mol.GetConformer().GetPositions()
    mask=np.triu(topology>2,k=1)  # Exclude bonded and angle-related atom pairs.
    if not mask.any():return None
    distances=np.linalg.norm(coordinates[:,None,:]-coordinates[None,:,:],axis=2)
    minimum=float(distances[mask].min())
    if minimum<cutoff:raise ValidationError('Severe intraligand nonbonded heavy-atom overlap')
    return minimum

def validate_ligand_identity(pdb, selector, smiles_path):
    """Check the complete heavy graph and 3D stereochemistry against full PROTAC."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    full=Chem.MolFromSmiles(Path(smiles_path).read_text().split()[0])
    if full is None:raise ValidationError('Invalid expected molecular graph')
    atoms=pdb_atoms(pdb)
    chosen=[a for a in atoms if a['chain'] in selector['chains'] and a['record']==selector['record']
            and a['element'] not in {'H','D'}
            and ('resname' not in selector or a['resname']==selector['resname'])
            and ('resnums' not in selector or a['resnum'] in selector['resnums'])]
    block='\n'.join(a['line'] for a in chosen)+'\nEND\n'
    mol=Chem.MolFromPDBBlock(block,removeHs=True,sanitize=False,proximityBonding=True)
    if mol is None or mol.GetNumAtoms()!=full.GetNumAtoms():
        raise ValidationError('Incomplete/wrong full-PROTAC heavy-atom count')
    try:
        mol=AllChem.AssignBondOrdersFromTemplate(full,mol)
        Chem.AssignStereochemistryFrom3D(mol,replaceExistingTags=True)
    except (ValueError,RuntimeError) as exc:raise ValidationError('Ligand connectivity cannot match full PROTAC') from exc
    if Chem.MolToInchiKey(mol)!=Chem.MolToInchiKey(full):
        raise ValidationError('Output ligand stereochemistry/identity differs from full PROTAC')
    minimum=check_ligand_nonbonded_geometry(mol)
    return {'min_intraligand_nonbonded_heavy_A':minimum}

def validate_stage_report(run, n_scored):
    """Do not hide failed child workers behind the surviving pose set."""
    report=json.loads((Path(run)/'engine_status.json').read_text())
    if report.get('status') not in {'COMPLETE','COMPLETE_NO_VALID_POSES'} or report.get('score_role')!='raw':
        raise ValidationError('Missing successful raw-score engine status')
    if not report.get('engine_sha256') or type(report.get('seed')) is not int:
        raise ValidationError('Engine provenance and explicit integer seed required')
    stages=report.get('stages',{})
    for stage in ['global_sampling','local_docking','linker_modeling','scoring']+(['pose_qc'] if 'pose_qc' in stages else []):
        counts=stages.get(stage,{})
        keys=['attempted','succeeded','rejected','failed']
        if any(type(counts.get(k)) is not int or counts[k]<0 for k in keys):
            raise ValidationError('Invalid/missing stage accounting: '+stage)
        if counts['attempted']!=sum(counts[k] for k in ['succeeded','rejected','failed']):
            raise ValidationError('Inconsistent stage accounting: '+stage)
        if counts['failed']:raise ValidationError('Failed child workers in stage: '+stage)
    if 'pose_qc' in stages:
        if stages['pose_qc']['attempted']!=stages['scoring']['succeeded'] or stages['pose_qc']['succeeded']!=n_scored:
            raise ValidationError('Post-score QC/retained pose denominator mismatch')
    elif stages['scoring']['succeeded']!=n_scored or stages['linker_modeling']['succeeded']!=n_scored:
        raise ValidationError('Stage report/scored pose count mismatch')
    if stages['linker_modeling']['attempted']!=stages['local_docking']['succeeded']:
        raise ValidationError('Local docking/linker denominator mismatch')
    if stages['scoring']['attempted']!=stages['linker_modeling']['succeeded']:
        raise ValidationError('Linker/scoring denominator mismatch')
    if (report['status']=='COMPLETE_NO_VALID_POSES')!=(n_scored==0):
        raise ValidationError('Terminal status inconsistent with valid pose count')
    if report.get('engine_kind')=='existing_pRosettaC_baseline_repaired':
        raw=Path(run)/'engine_raw'
        if sha256(raw/'Patchdock_Results/score.sc')!=report.get('raw_score_sha256'):
            raise ValidationError('Raw baseline score table changed')
        if sha256(raw/'baseline_output_qc.json')!=report.get('raw_output_qc_sha256'):
            raise ValidationError('Raw baseline QC report changed')
        for name,digest in report.get('sampling_sha256',{}).items():
            if Path(name).name!=name or sha256(raw/'Patchdock_Results'/name)!=digest:
                raise ValidationError('Baseline sampling ledger changed')
    return report
