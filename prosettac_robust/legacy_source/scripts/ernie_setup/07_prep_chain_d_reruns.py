#!/usr/bin/env python3
"""
Prepare chain-D-fixed pRosettaC re-runs for 4 artifact POIs:
  TNIK, MAPK14, MAP4K5, SLK

Run this LOCALLY before transferring to ernie.

What it does:
  1. MAPK14 / MAP4K5 / SLK:
     - Creates <POI>_protein_clean_D.pdb with chain A renamed to D
     - Updates Protac_params.txt: Chains: ABC D, path → _D.pdb
  2. TNIK:
     - protein_clean.pdb is ALREADY chain D and Protac_params.txt is already correct
     - Just clears the stale Patchdock_Results and Results so ernie runs fresh
  3. For all 4: clears Results/ (so the ">50 PDBs → skip" guard doesn't fire on ernie)

After this script:
  Transfer:
    rsync -avz --progress \
      "/mnt/d/PROTAC degradation Prediction/pRosettaC/runs/low_affinity_degraders/TNIK/" \
      "/mnt/d/PROTAC degradation Prediction/pRosettaC/runs/low_affinity_degraders/MAPK14/" \
      "/mnt/d/PROTAC degradation Prediction/pRosettaC/runs/high_affinity_no_degradation/MAP4K5/" \
      "/mnt/d/PROTAC degradation Prediction/pRosettaC/runs/high_affinity_no_degradation/SLK/" \
      "/mnt/d/PROTAC degradation Prediction/data/structures/POI/high_affinity_degraders/MAPK14/" \
      "/mnt/d/PROTAC degradation Prediction/data/structures/POI/high_affinity_no_degradation/MAP4K5/" \
      "/mnt/d/PROTAC degradation Prediction/data/structures/POI/high_affinity_no_degradation/SLK/" \
      s1638308@ernie.bio.ed.ac.uk:~/protac/ ...
  Then run: bash ~/protac/scripts/ernie_setup/08_ernie_launch_chain_d.sh
"""

import shutil
from pathlib import Path

BASE    = Path("/mnt/d/PROTAC degradation Prediction")
RUNS    = BASE / "pRosettaC" / "runs"
STRUCTS = BASE / "data" / "structures" / "POI"

# ── Per-POI config ─────────────────────────────────────────────────────────────
POIS = {
    'TNIK': {
        'run_base':    'low_affinity_degraders',
        'struct_base': 'low_affinity_degraders',
        'already_D':   True,   # protein_clean.pdb + Protac_params already fixed
    },
    'MAPK14': {
        'run_base':    'low_affinity_degraders',
        'struct_base': 'high_affinity_degraders',  # Protac_params.txt points here
        'already_D':   False,
    },
    'MAP4K5': {
        'run_base':    'high_affinity_no_degradation',
        'struct_base': 'high_affinity_no_degradation',
        'already_D':   False,
    },
    'SLK': {
        'run_base':    'high_affinity_no_degradation',
        'struct_base': 'high_affinity_no_degradation',
        'already_D':   False,
    },
}


def rename_chain_to_D(src_pdb: Path, dst_pdb: Path):
    """Copy PDB renaming all chain A atoms to chain D. Preserves all formatting."""
    lines_out = []
    with open(src_pdb) as f:
        for line in f:
            rec = line[:6].strip()
            if rec in ('ATOM', 'HETATM') and line[21] == 'A':
                line = line[:21] + 'D' + line[22:]
            lines_out.append(line)
    with open(dst_pdb, 'w') as f:
        f.writelines(lines_out)
    print(f"  Written: {dst_pdb.name}")


def clear_patchdock_results(run_dir: Path):
    pd = run_dir / 'Patchdock_Results'
    if pd.exists():
        files = list(pd.glob('combined_*.pdb')) + list(pd.glob('PT_*.pdb'))
        n = 0
        for f in files:
            try:
                n += 1
            except (FileNotFoundError, PermissionError):
                pass
        print(f"  Cleared {n}/{len(files)} files from Patchdock_Results/")
    for fname in ('score.sc', 'score.sc.orig', 'log.txt', 'random_sampling.sdf'):
        f = run_dir / fname
        try:
            if f.exists():
                print(f"  Removed {fname}")
        except (FileNotFoundError, PermissionError):
            pass


def clear_results(run_dir: Path):
    res = run_dir / 'Results'
    if res.exists():
        pdbs = list(res.glob('**/*.pdb'))
        n = 0
        for f in pdbs:
            try:
                n += 1
            except (FileNotFoundError, PermissionError):
                pass
        for d in res.iterdir():
            if d.is_dir():
                try:
                except OSError:
                    pass
        print(f"  Cleared {n}/{len(pdbs)} PDBs from Results/")


def update_protac_params(params_path: Path, poi: str, new_pdb_path: str):
    """Update Protac_params.txt: set Chains: ABC D and update structure path."""
    text = params_path.read_text()
    lines = text.splitlines()
    new_lines = []
    for line in lines:
        if line.startswith('Chains:'):
            new_lines.append('Chains: ABC D')
        elif line.startswith('Structures:'):
            # Replace only the second path (kinase) — first path is VHL
            parts = line.split()
            # parts[0] = 'Structures:', parts[1] = VHL path, parts[2] = kinase path
            new_lines.append(f'Structures: {parts[1]} {new_pdb_path}')
        else:
            new_lines.append(line)
    params_path.write_text('\n'.join(new_lines) + '\n')
    print(f"  Updated: Protac_params.txt → Chains: ABC D, path → .../{Path(new_pdb_path).name}")


def main():
    print("=== Preparing chain-D pRosettaC re-runs for 4 artifact POIs ===\n")

    for poi, cfg in POIS.items():
        print(f"\n[{poi}]")
        run_dir    = RUNS / cfg['run_base'] / poi
        struct_dir = STRUCTS / cfg['struct_base'] / poi
        params     = run_dir / 'Protac_params.txt'

        if not cfg['already_D']:
            # 1. Create chain D version of protein_clean.pdb
            src_pdb = struct_dir / f'{poi}_protein_clean.pdb'
            dst_pdb = struct_dir / f'{poi}_protein_clean_D.pdb'
            if not src_pdb.exists():
                print(f"  ERROR: {src_pdb} not found — skipping")
                continue
            if dst_pdb.exists():
                print(f"  {dst_pdb.name} already exists — overwriting")
            rename_chain_to_D(src_pdb, dst_pdb)

            # 2. Update Protac_params.txt
            # The ernie path mirrors the local structure: swap the local mount prefix
            # /mnt/d/PROTAC degradation Prediction → /home/s1638308/protac
            ernie_pdb_path = str(dst_pdb).replace(
                '/mnt/d/PROTAC degradation Prediction',
                '/home/s1638308/protac'
            )
            update_protac_params(params, poi, ernie_pdb_path)

            # 3. Backup the old Side1.pdb and create a chain D version
            side1_src = run_dir / 'Side1.pdb'
            if side1_src.exists() and not (run_dir / 'Side1_A.pdb').exists():
                shutil.copy(side1_src, run_dir / 'Side1_A.pdb')
            rename_chain_to_D(src_pdb, run_dir / 'Side1.pdb')

        else:
            print(f"  protein_clean.pdb already chain D — no PDB fix needed")
            print(f"  Protac_params.txt already Chains: ABC D — no param fix needed")

        # 4. Clear stale Patchdock_Results (old chain A combined PDBs)
        clear_patchdock_results(run_dir)

        # 5. Clear Results (so ernie doesn't skip with ">50 PDBs" guard)
        clear_results(run_dir)

        print(f"  [{poi}] Ready for fresh ernie re-run")

    print("\n=== All 4 POIs prepared ===")
    print("\nNext steps:")
    print("1. Transfer run dirs and structure files to ernie:")
    print("   bash scripts/ernie_setup/transfer_chain_d.sh   (see below)")
    print("2. SSH to ernie and run:")
    print("   bash ~/protac/scripts/ernie_setup/08_ernie_launch_chain_d.sh")
    print("\nRsync commands (run from local WSL):")
    for poi, cfg in POIS.items():
        run_dir = RUNS / cfg['run_base'] / poi
        print(f"  rsync -avz --exclude 'Patchdock_Results/combined_*' \\")
        print(f"    '{run_dir}/' s1638308@ernie.bio.ed.ac.uk:~/protac/pRosettaC/runs/{cfg['run_base']}/{poi}/")
    print()
    for poi, cfg in POIS.items():
        if not cfg['already_D']:
            struct_dir = STRUCTS / cfg['struct_base'] / poi
            print(f"  rsync -avz '{struct_dir}/' s1638308@ernie.bio.ed.ac.uk:~/protac/data/structures/POI/{cfg['struct_base']}/{poi}/")


if __name__ == '__main__':
    main()
