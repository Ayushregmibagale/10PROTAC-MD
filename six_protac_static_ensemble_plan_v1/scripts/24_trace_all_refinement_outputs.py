"""Retrospective all-pose refinement trace; never selects generation inputs."""
import argparse
import csv
import importlib.util
import json
import re
from pathlib import Path

import numpy as np


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('runs', 'references', 'output'):
        parser.add_argument('--' + key, type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    scripts = Path(__file__).parent
    reader = module('trace_reader', scripts / '08_diagnose_global_coverage.py')
    recovery = module('trace_recovery', scripts / '01_assess_reference_recovery.py')
    rows, summaries = [], []
    for ref in ('5T35_MZ1', '7PI4_GSK215'):
        run = args.runs / f'global_{ref}_seed42'
        status = json.loads((run / 'engine_status.json').read_text())
        raw = run / 'engine_raw' / 'Patchdock_Results'
        native = reader.protein_roles(args.references / ref / 'native_complex.pdb')

        def rmsd(path):
            pose = reader.protein_roles(path)
            if any(pose[k][0] != native[k][0] for k in ('vhl', 'poi')):
                raise ValueError(f'Residue correspondence differs: {path}')
            rotation, translation = recovery.transform(pose['vhl'][1], native['vhl'][1])
            return float(np.sqrt(np.mean(np.sum(
                (pose['poi'][1] @ rotation + translation - native['poi'][1]) ** 2, axis=1))))

        globals_ = {int(p.name.split('.')[1]): p for p in raw.glob('pd.*.pdb')
                    if re.fullmatch(r'pd\.\d+\.pdb', p.name)}
        initial = {index: rmsd(path) for index, path in sorted(globals_.items())}
        rr = []
        for path in sorted(raw.glob('pd.*_docking_*.pdb')):
            match = re.fullmatch(r'pd\.(\d+)_docking_(\d+)\.pdb', path.name)
            if match is None:
                if re.fullmatch(r'pd\.\d+_docking_\d+_ABCD\.pdb', path.name):
                    continue  # Linker-stage protein-only derivative, not a local output.
                raise ValueError(f'Unexpected local filename: {path}')
            index, local = map(int, match.groups())
            value = rmsd(path)
            ledger_path = raw / f'confs_{index}_{local}.sdf.sampling.json'
            ledger = json.loads(ledger_path.read_text())
            row = dict(reference=ref, global_index=index, local_index=local,
                       global_RMSD_A=initial[index], local_RMSD_A=value,
                       local_minus_global_RMSD_A=value-initial[index],
                       linker_status=ledger['status'],
                       global_sha256=recovery.sha(globals_[index]),
                       local_sha256=recovery.sha(path))
            rr.append(row)
        expected = status['stages']['local_docking']['succeeded']
        if len(rr) != expected or not rr or not initial:
            raise ValueError(f'Incomplete trace {ref}: {len(rr)} != {expected}')
        feasible = [r for r in rr if r['linker_status'] == 'FEASIBLE_CONFORMERS_GENERATED']
        summaries.append(dict(reference=ref, global_count=len(initial), local_count=len(rr),
            best_global_RMSD_A=min(initial.values()),
            best_local_RMSD_A=min(r['local_RMSD_A'] for r in rr),
            global_RMSD_le5A=sum(v <= 5 for v in initial.values()),
            local_RMSD_le5A=sum(r['local_RMSD_A'] <= 5 for r in rr),
            feasible_local_RMSD_le5A=sum(r['local_RMSD_A'] <= 5 for r in feasible),
            feasible_linker_count=len(feasible),
            median_RMSD_change_A=float(np.median([r['local_minus_global_RMSD_A'] for r in rr]))))
        rows.extend(rr)
    with (args.output / 'ALL_REFINEMENT_TRACE.tsv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    report = dict(status='ALL_SAVED_REFINEMENTS_EVALUATED', summary=summaries,
        source_sha256=recovery.sha(Path(__file__)), limitations=[
            'Native coordinates used for retrospective measurement only, never for selecting a new search.',
            'RMSD-only stage trace; final frozen gate also requires contact recovery and chemical QC.',
            'One seed; not a convergence or biological validation study.'])
    (args.output / 'ALL_REFINEMENT_TRACE.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
