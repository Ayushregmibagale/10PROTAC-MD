"""Stage task inputs for original pRosettaC, preserving source input hashes."""
import argparse
import json
from pathlib import Path
import shutil

from check_paths import sha

COMPOUNDS = ('GSK215','SJF8240','SJF1521','JNJ1013','SJF_delta','MS143')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--staged', type=Path, required=True,
                        help='Directory containing the existing six-drug INPUT_MANIFEST.json')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--full', choices=('True','False'), required=True,
                        help='Original protocol setting; False=500 global x10 local; True=1000x50')
    parser.add_argument('--engine-config', type=Path,
                        help='Override installation paths/resource settings when building a new frozen bundle')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    source = Path(__file__).resolve().parent
    for name in ('campaign.py','check_paths.py','run_original.py','submit_ernie.sh'):
        shutil.copy2(source/name, args.output/name)
    staged = json.loads((args.staged/'INPUT_MANIFEST.json').read_text())
    tasks = []
    for row in staged['controls']:
        if row['compound'] not in COMPOUNDS:
            raise ValueError('Unexpected compound')
        original = args.staged/row['directory']
        for name, digest in row['input_file_sha256'].items():
            if sha(original/name) != digest:
                raise ValueError('Staged source input changed: '+str(original/name))
        relative = Path('inputs')/row['compound']/row['target']
        target = args.output/relative
        shutil.copytree(original, target)
        params = target/'Protac_params.txt'
        lines = params.read_text().splitlines()
        if not any(line.startswith('Full:') for line in lines):
            raise ValueError('Missing original Full switch')
        params.write_text('\n'.join('Full: '+args.full if line.startswith('Full:') else
            'ClusterName: LOCAL' if line.startswith('ClusterName:') else line for line in lines)+'\n')
        tasks.append(dict(id=row['compound']+'__'+row['target'], compound=row['compound'],
            target=row['target'], role=row['role'], input=str(relative),
            source_files=row['input_file_sha256'],
            files={p.name:sha(p) for p in sorted(target.iterdir()) if p.is_file()}))
    if {t['compound'] for t in tasks} != set(COMPOUNDS):
        raise ValueError('Not all six compounds staged')
    sampling = dict(original_Full=args.full, global_results_limit=1000 if args.full=='True' else 500,
                    local_per_global=50 if args.full=='True' else 10,
                    constrained_conformers='Original GenConstConf default, nominal up to 100; no replacement sampler')
    manifest = dict(engine='Original LondonLab pRosettaC + PatchDock + standalone Rosetta',
        scope='Six existing technical inputs; not all kinase targets or a biological validation cohort',
        sampling=sampling, tasks=tasks, MD=False)
    (args.output/'TASKS.json').write_text(json.dumps(manifest,indent=2)+'\n')
    cfg=dict(source='/home/s1638308/PRosettaC',patchdock='/home/s1638308/PatchDock',
             rosetta='/home/s1638308/rosetta',openbabel='/home/s1638308/miniconda3/envs/protac/bin',
             python='/home/s1638308/miniconda3/envs/protac/bin/python3', workers=4,
             task_timeout_seconds=43200, min_free_bytes=50*1024**3)
    if args.engine_config:
        cfg.update(json.loads(args.engine_config.read_text()))
    (args.output/'engine_config.json').write_text(json.dumps(cfg,indent=2)+'\n')
    (args.output/'BUNDLE_MANIFEST.json').write_text(json.dumps(dict(
        scope='Executable scripts and input provenance; remote installation paths remain subject to check_paths.py',
        files={str(p.relative_to(args.output)):sha(p) for p in sorted(args.output.rglob('*')) if p.is_file()}),indent=2)+'\n')
    print(json.dumps(dict(bundle=str(args.output),tasks=len(tasks),sampling=sampling),indent=2))


if __name__ == '__main__':
    main()
