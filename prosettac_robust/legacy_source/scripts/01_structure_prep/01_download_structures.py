"""
Download and clean PDB structures for all POIs and VHL.
Outputs cleaned PDBs into data/structures/POI/<category>/<TARGET>/
"""

import os
import subprocess
import urllib.request
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = PROJECT_ROOT / "data" / "structures"

POI_TABLE = [
    # (target, category, pdb_id, chain, priority)
    ("MET",     "high_affinity_degraders",     "3DKC", "A", 1),
    ("DDR2",    "high_affinity_degraders",     "3ZOS", "A", 1),
    ("RIPK2",   "high_affinity_degraders",     "5AR3", "A", 1),
    ("EPHB2",   "high_affinity_degraders",     "2HEN", "A", 1),
    ("MAPK14",  "high_affinity_degraders",     "3GCP", "A", 1),
    ("AXL",     "high_affinity_no_degradation","5U6B", "A", 2),
    ("SLK",     "high_affinity_no_degradation","2J51", "A", 2),
    ("ABL1",    "high_affinity_no_degradation","2HYY", "A", 2),
    ("EPHA2",   "high_affinity_no_degradation","3MBW", "A", 2),
    ("MAP4K5",  "high_affinity_no_degradation","6YSS", "A", 2),
    ("TNIK",    "low_affinity_degraders",      "4UXT", "A", 3),
    # MAPK14 reuses same PDB as above — symlink rather than re-download
    ("PIP4K2C", "low_affinity_degraders",      "2GK9", "A", 3),
    # CDK17 has no crystal structure — placeholder for AlphaFold model
    # ── Priority 4: low_affinity_no_degradation (SJF8240 panel expansion) ─────
    # All need chain A→D rename to avoid residue-number conflict with VHL 1-103
    ("CDK2",    "low_affinity_no_degradation", "1H1P", "A", 4),
    ("CDK6",    "low_affinity_no_degradation", "2EUF", "B", 4),
    ("CDK4",    "low_affinity_no_degradation", "2W9Z", "B", 4),
    ("DAPK1",   "low_affinity_no_degradation", "1JKS", "A", 4),
    ("TLK2",    "low_affinity_no_degradation", "5O0Y", "A", 4),
    ("IRAK1",   "low_affinity_no_degradation", "6BFN", "A", 4),
]

VHL_PDB = ("VHL", "4W9H", "B")   # chain B is VHL in 4W9H


def fetch_pdb(pdb_id: str, out_path: Path) -> bool:
    """Download PDB file. Returns True on success, False if not found (AF2/withdrawn)."""
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    print(f"  Fetching {pdb_id} → {out_path.name}")
    try:
        urllib.request.urlretrieve(url, out_path)
        return True
    except urllib.error.HTTPError as e:
        if e.code == 404:
            print(f"    SKIP: {pdb_id} not in RCSB (AlphaFold or withdrawn entry)")
            return False
        raise


def extract_chain(in_pdb: Path, chain: str, out_pdb: Path) -> None:
    """Keep only ATOM/HETATM records for the requested chain."""
    lines = []
    with open(in_pdb) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")) and line[21] == chain:
                lines.append(line)
            elif line.startswith(("TER", "END")):
                lines.append(line)
    with open(out_pdb, "w") as fh:
        fh.writelines(lines)
    print(f"    Chain {chain} extracted → {out_pdb.name}")


def rename_chain(in_pdb: Path, chain_from: str, chain_to: str, out_pdb: Path) -> None:
    """Rename chain identifier (col 22) from chain_from to chain_to in all ATOM/HETATM/TER lines."""
    lines = []
    with open(in_pdb) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM", "TER")) and len(line) > 21 and line[21] == chain_from:
                line = line[:21] + chain_to + line[22:]
            lines.append(line)
    with open(out_pdb, "w") as fh:
        fh.writelines(lines)
    print(f"    Chain {chain_from}→{chain_to} renamed → {out_pdb.name}")


def main():
    # VHL
    vhl_dir = DATA_DIR / "VHL"
    raw_vhl = vhl_dir / f"{VHL_PDB[1]}.pdb"
    if not raw_vhl.exists():
        fetch_pdb(VHL_PDB[1], raw_vhl)
    chain_vhl = vhl_dir / f"{VHL_PDB[1]}_chain{VHL_PDB[2]}.pdb"
    if not chain_vhl.exists():
        extract_chain(raw_vhl, VHL_PDB[2], chain_vhl)

    # POIs
    for target, category, pdb_id, chain, priority in POI_TABLE:
        out_dir = DATA_DIR / "POI" / category / target
        out_dir.mkdir(parents=True, exist_ok=True)

        raw_pdb = out_dir / f"{pdb_id}.pdb"
        if not raw_pdb.exists():
            ok = fetch_pdb(pdb_id, raw_pdb)
            if not ok:
                continue   # AF2 or withdrawn — no chain extraction or rename
        else:
            print(f"  {target} ({pdb_id}) already present, skipping download.")

        if not raw_pdb.exists():
            continue   # download was skipped (AF2 placeholder in table)

        chain_pdb = out_dir / f"{target}_{pdb_id}_chain{chain}.pdb"
        if not chain_pdb.exists():
            extract_chain(raw_pdb, chain, chain_pdb)

        # low_affinity_no_degradation: rename kinase chain A→D to avoid
        # residue-number conflict with VHL chain A (residues 1-103) in pRosettaC
        if category == "low_affinity_no_degradation":
            clean_pdb = out_dir / f"{target}_protein_clean.pdb"
            if not clean_pdb.exists():
                rename_chain(chain_pdb, chain, "D", clean_pdb)

    # MAPK14 low-affinity — symlink to high-affinity entry
    src = DATA_DIR / "POI" / "high_affinity_degraders" / "MAPK14"
    dst = DATA_DIR / "POI" / "low_affinity_degraders" / "MAPK14"
    dst.mkdir(parents=True, exist_ok=True)
    for f in src.iterdir():
        link = dst / f.name
        if not link.exists() and not link.is_symlink():
            link.symlink_to(f)
    print("  MAPK14 low_affinity_degraders symlinked to high_affinity entry.")

    print("\nDone. Next step: place warhead coordinates in each POI structure.")
    print("See scripts/structure_prep/place_warhead.py")


if __name__ == "__main__":
    main()
