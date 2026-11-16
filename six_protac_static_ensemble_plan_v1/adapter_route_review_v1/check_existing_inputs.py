"""Read-only audit of the proposed existing multi-compound input directories."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    rows = []
    required = ('Protac','Anchor atoms','Full')
    for compound in ('GSK215','JNJ1013','SJF1521','SJF8240','SJF_delta','MS143'):
        directory = args.inputs/compound
        files = sorted(directory.glob('*/Protac_params.txt'))
        missing, molecules, chains = Counter(), Counter(), Counter()
        provenance = hashlib.sha256()
        for path in files:
            parameters = dict((key.strip(),value.strip()) for line in path.read_text().splitlines()
                              if ':' in line for key,value in [line.split(':',1)])
            missing.update(key for key in required if key not in parameters)
            molecules.update([parameters.get('Protac',parameters.get('Smiles','MISSING'))])
            chains.update([parameters.get('Chains','MISSING')])
            provenance.update(str(path.relative_to(args.inputs)).encode()+b'\0')
            provenance.update(hashlib.sha256(path.read_bytes()).digest())
        rows.append(dict(compound=compound,directory_exists=directory.is_dir(),
                         parameter_files=len(files),missing_required_keys=dict(missing),
                         declared_molecule_paths=dict(molecules),chain_selections=dict(chains),
                         parameter_set_digest=provenance.hexdigest()))
    report = dict(scope='Parameter-file audit only; not molecule/pose validation',
                  input_root=str(args.inputs),compounds=rows)
    args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__ == '__main__':
    main()
