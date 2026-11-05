#!/usr/bin/env python3
"""Enumerate candidate protonation states with parent graph/stereo integrity.

pH 7.4 is a provisional physiological assumption unless assay conditions are
specified. Enumeration does not provide measured populations or select a primary
MD state. Supplier-reference structures/parameters remain separately identified.
"""
import argparse,json,hashlib,importlib.metadata as metadata,sys
from pathlib import Path
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from dimorphite_dl.protonate.run import Protonate
import dimorphite_dl
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from protac_qc import ValidationError,atomic_json,sha256

def assess(root,output,ph=7.4):
    if not 0<=ph<=14:raise ValidationError('pH must lie within 0–14')
    if output.exists():raise ValidationError('Assessment directory already used')
    output.mkdir(parents=True)
    registry=json.loads((root/'config/molecules.json').read_text())['compounds']
    scenarios=[('assumed_pH_'+str(ph),ph,ph),('sensitivity_6p5_to_8p0',6.5,8.0)]
    report={'status':'PREDICTED_CANDIDATES_ONLY','primary_MD_state_selected':False,
        'pH_basis':'Explicit provisional physiological assumption; replace with assay conditions if supplied',
        'assumed_pH':ph,'precision_factor':1.0,'max_variants':128,
        'method':'Dimorphite-DL ionization-state enumeration; no measured population estimates',
        'versions':{n:metadata.version(n) for n in ['dimorphite-dl','rdkit','loguru']},
        'sources':['https://github.com/durrantlab/dimorphite_dl/releases/tag/v2.1.0',
                   'https://durrantlab.github.io/dimorphite_dl/api/protonate/run/'],
        'compounds':{},'tool_source_sha256':{}}
    package=Path(dimorphite_dl.__file__).parent
    for file in sorted(package.rglob('*')):
        if file.is_file() and '__pycache__' not in file.parts:
            report['tool_source_sha256'][str(file.relative_to(package))]=sha256(file)
    for compound,supplier in registry.items():
        source=root/'molecules'/compound/'full.smi';mol=Chem.MolFromSmiles(source.read_text().strip())
        if Chem.MolToInchiKey(mol)!=supplier['inchikey']:raise ValidationError('Supplier parent mismatch')
        source_signature={i+1:(a.GetSymbol(),a.GetIsotope()) for i,a in enumerate(mol.GetAtoms())}
        source_charge={i+1:a.GetFormalCharge() for i,a in enumerate(mol.GetAtoms())}
        for i,atom in enumerate(mol.GetAtoms()):atom.SetAtomMapNum(i+1)
        states={};scenario_reports={}
        for scenario,low,high in scenarios:
            engine=Protonate(Chem.MolToSmiles(mol),ph_min=low,ph_max=high,
                precision=1.0,max_variants=128,validate_output=True)
            raw=list(engine);stats=engine.get_stats()
            if stats['protonation']['fallback_used'] or stats['processor']['errors']:
                raise ValidationError('Enumeration fallback/error for '+compound)
            if len(raw)>=128:raise ValidationError('Enumeration may be truncated for '+compound)
            ids=[]
            for smiles in raw:
                state=Chem.MolFromSmiles(smiles)
                if state is None:raise ValidationError('Invalid predicted state')
                signature={a.GetAtomMapNum():(a.GetSymbol(),a.GetIsotope()) for a in state.GetAtoms()}
                if signature!=source_signature or len(signature)!=state.GetNumAtoms():
                    raise ValidationError('Heavy-atom identity/map changed during protonation')
                parent=rdMolStandardize.Uncharger().uncharge(Chem.Mol(state))
                if Chem.MolToInchiKey(parent)!=supplier['inchikey']:
                    raise ValidationError('Predicted state does not return the exact supplier graph/stereochemistry')
                changes=[{'supplier_atom_one_based':a.GetAtomMapNum(),'element':a.GetSymbol(),
                    'reference_formal_charge':source_charge[a.GetAtomMapNum()],
                    'candidate_formal_charge':a.GetFormalCharge()}
                    for a in state.GetAtoms() if a.GetFormalCharge()!=source_charge[a.GetAtomMapNum()]]
                mapped=Chem.MolToSmiles(state,isomericSmiles=True)
                clean=Chem.Mol(state)
                for atom in clean.GetAtoms():atom.SetAtomMapNum(0)
                canonical=Chem.MolToSmiles(clean,isomericSmiles=True)
                charge=sum(a.GetFormalCharge() for a in clean.GetAtoms())
                state_id='q'+('pos'+str(charge) if charge>0 else 'neg'+str(-charge) if charge<0 else '0')+'_'+hashlib.sha256(canonical.encode()).hexdigest()[:12]
                item={'state_id':state_id,'smiles':canonical,'mapped_smiles':mapped,
                    'inchikey':Chem.MolToInchiKey(clean),'parent_supplier_inchikey':supplier['inchikey'],
                    'formal_charge':charge,'heavy_atoms':clean.GetNumAtoms(),
                    'changed_charge_sites':changes,'parent_graph_and_stereo_verified':True,
                    'predicted_population':None,'production_selected':False}
                states[state_id]=item;ids.append(state_id)
            scenario_reports[scenario]={'pH_range':[low,high],'state_ids':sorted(set(ids)),'tool_stats':stats}
        report['compounds'][compound]={'catalog':supplier['catalog'],'supplier_inchikey':supplier['inchikey'],
            'supplier_smiles_sha256':sha256(source),'scenarios':scenario_reports,'states':states,
            'tautomer_scope':'Ionization enumeration only; exhaustive tautomer populations are not inferred'}
        folder=output/compound;folder.mkdir()
        for state_id,state in sorted(states.items()):(folder/(state_id+'.smi')).write_text(state['smiles']+'\n')
    atomic_json(output/'assessment.json',report)
    summary={c:{'states':len(d['states']),'formal_charges':sorted({s['formal_charge'] for s in d['states'].values()}),
        'by_scenario':{scenario:{'n_states':len(info['state_ids']),
            'formal_charges':sorted({d['states'][sid]['formal_charge'] for sid in info['state_ids']})}
            for scenario,info in d['scenarios'].items()}}
        for c,d in report['compounds'].items()}
    print(json.dumps({'status':report['status'],'assumed_pH':ph,'compounds':summary},indent=2))
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2]);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--ph',type=float,default=7.4);a=p.parse_args();assess(a.root,a.output,a.ph)
