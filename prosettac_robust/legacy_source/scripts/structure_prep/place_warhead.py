"""
Place SJF8240 kinase warhead into each POI binding site.

Strategy:
  1. REFERENCE: MET/3LQ8 contains 88Z (foretinib analog) — same pharmacophore
     as SJF8240 warhead. Extract 88Z 3D coords and use them directly as the
     placed warhead (remapped via MCS to SJF8240_warhead topology).
  2. ALL OTHER KINASES: sequence-guided Cα superposition onto 3LQ8 chain A
     using BioPython pairwise alignment + Kabsch superimposer, then apply
     the same rigid transform to warhead coordinates.

Output per target:
  data/structures/POI/<category>/<TARGET>/<TARGET>_warhead.pdb
  = protein ATOM records + warhead HETATM records (resname WRH, chain L)

Requires: RDKit, BioPython, numpy
"""

import sys, warnings
import numpy as np
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS
from rdkit.Chem.rdchem import RWMol
from Bio import PDB
from Bio.PDB import Superimposer
from Bio.Align import PairwiseAligner

warnings.filterwarnings("ignore")

PROJECT_ROOT = Path(__file__).resolve().parents[2]
POI_DIR      = PROJECT_ROOT / "data" / "structures" / "POI"

REF_PDB      = POI_DIR / "high_affinity_degraders" / "MET" / "MET_3LQ8_chainA.pdb"
REF_LIG_NAME = "88Z"

# SJF8240 warhead SMILES (foretinib core; OH = former PEG linker ether oxygen)
SMILES_WH = (
    "OC1=CC2=NC=CC(OC3=CC=C(C=C3F)N(C(=O)C3(CC3)C(N)=O)"
    "C3C=CC(F)=CC=3)=C2C=C1OC"
)
# 88Z SMILES (foretinib analog in 3LQ8; has morpholinopropoxy instead of OH)
SMILES_88Z = (
    "COc1cc2cc(Oc3ccc(NC(=O)C4(CC4)C(=O)Nc4ccc(F)cc4)c(F)c3)"
    "cnc2cc1OCCCN1CCOCC1"
)

# Three-letter → one-letter amino acid map
AA3TO1 = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C',
    'GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I',
    'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P',
    'SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V',
    'MSE':'M','HSD':'H','HSE':'H','HSP':'H',
}

POI_TABLE = [
    # (target, category, pdb_stem)
    ("MET",     "high_affinity_degraders",    "MET_3LQ8_chainA"),
    ("DDR2",    "high_affinity_degraders",    "DDR2_7AYM_chainA"),
    ("RIPK2",   "high_affinity_degraders",    "RIPK2_4C8B_chainA"),
    ("EPHB2",   "high_affinity_degraders",    "EPHB2_3ZFM_chainA"),
    ("MAPK14",  "high_affinity_degraders",    "MAPK14_3GCP_chainA"),
    ("AXL",     "high_affinity_no_degradation","AXL_5U6B_chainA"),
    ("SLK",     "high_affinity_no_degradation","SLK_2J51_chainA"),
    ("ABL1",    "high_affinity_no_degradation","ABL1_2HYY_chainA"),
    ("EPHA2",   "high_affinity_no_degradation","EPHA2_5I9Y_chainA"),
    ("MAP4K5",  "high_affinity_no_degradation","MAP4K5_kinase_domain"),
    ("TNIK",    "low_affinity_degraders",     "TNIK_5D7A_chainA"),
    ("PIP4K2C", "low_affinity_degraders",     "PIP4K2C_2GK9_chainA"),
    ("CDK17",   "low_affinity_degraders",     "CDK17_kinase_domain"),
]


# ─── helpers ─────────────────────────────────────────────────────────────────

def extract_ca_and_seq(pdb_path):
    """Return (residue_id_list, CA_atom_list, sequence_str) for chain A."""
    parser = PDB.PDBParser(QUIET=True)
    struct = parser.get_structure("x", str(pdb_path))
    res_ids, ca_atoms, seq = [], [], []
    for model in struct:
        for chain in model:
            for res in chain:
                if res.id[0] != " ":    # skip HETATM residues
                    continue
                if "CA" not in res:
                    continue
                aa = AA3TO1.get(res.resname, "X")
                res_ids.append(res.id[1])
                ca_atoms.append(res["CA"])
                seq.append(aa)
        break  # first model only
    return res_ids, ca_atoms, "".join(seq)


def sequence_guided_superposition(ref_ca, ref_seq, tgt_ca, tgt_seq,
                                   n_cycles=5, cutoff_sigma=2.0):
    """
    Pairwise align ref_seq vs tgt_seq, collect aligned Cα pairs, then
    iteratively reject outliers (>cutoff_sigma * σ per-pair RMSD) until
    convergence. Returns (Superimposer, n_pairs, rmsd).
    """
    aligner = PairwiseAligner()
    aligner.substitution_matrix = _blosum62()
    aligner.open_gap_score   = -10
    aligner.extend_gap_score = -0.5
    aligner.mode = "global"

    aln = next(iter(aligner.align(ref_seq, tgt_seq)))

    all_ref, all_tgt = [], []
    for (r0, r1), (t0, t1) in zip(aln.aligned[0], aln.aligned[1]):
        for rr, tt in zip(range(r0, r1), range(t0, t1)):
            if rr < len(ref_ca) and tt < len(tgt_ca):
                all_ref.append(ref_ca[rr])
                all_tgt.append(tgt_ca[tt])

    if len(all_ref) < 10:
        return None, 0, 999.9

    ref_pairs, tgt_pairs = all_ref[:], all_tgt[:]
    sup = Superimposer()

    for cycle in range(n_cycles):
        sup.set_atoms(ref_pairs, tgt_pairs)
        # Compute per-pair distance after superposition
        rot, tran = np.array(sup.rotran[0]), np.array(sup.rotran[1])
        dists = []
        for ra, ta in zip(ref_pairs, tgt_pairs):
            rv = np.array(ra.get_vector().get_array())
            tv = np.array(ta.get_vector().get_array())
            tv_rot = rot @ tv + tran
            dists.append(np.linalg.norm(rv - tv_rot))
        dists = np.array(dists)
        mean_d, std_d = dists.mean(), dists.std()
        cutoff = mean_d + cutoff_sigma * std_d
        keep = dists < cutoff
        if keep.sum() == len(ref_pairs):
            break   # converged
        ref_pairs = [a for a, k in zip(ref_pairs, keep) if k]
        tgt_pairs = [a for a, k in zip(tgt_pairs, keep) if k]
        if len(ref_pairs) < 10:
            break

    sup.set_atoms(ref_pairs, tgt_pairs)
    return sup, len(ref_pairs), sup.rms


def _blosum62():
    """Return Bio.Align substitution matrix BLOSUM62."""
    from Bio.Align import substitution_matrices
    return substitution_matrices.load("BLOSUM62")


def extract_ligand_as_rdmol(pdb_path, resname):
    """Read all HETATM atoms for resname into an RDKit mol with 3D coords."""
    lines = [l for l in Path(pdb_path).read_text().splitlines()
             if l.startswith("HETATM") and l[17:20].strip() == resname]
    if not lines:
        return None
    rw = RWMol()
    coords = []
    for l in lines:
        elem = l[76:78].strip() or l[12:14].strip()[0]
        elem = elem.capitalize()
        try:
            rw.AddAtom(Chem.Atom(elem))
        except Exception:
            rw.AddAtom(Chem.Atom("C"))
        x, y, z = float(l[30:38]), float(l[38:46]), float(l[46:54])
        coords.append((x, y, z))
    conf = Chem.Conformer(rw.GetNumAtoms())
    for i, (x, y, z) in enumerate(coords):
        conf.SetAtomPosition(i, (x, y, z))
    rw.AddConformer(conf, assignId=True)
    Chem.SanitizeMol(rw, catchErrors=True)
    return rw.GetMol()


def mcs_remap_coords(wh_smiles, ref_smiles, ref_mol_3d):
    """
    Build the warhead mol and assign 3D coords from ref via MCS.
    Returns warhead mol with coordinates in the binding-site frame.
    """
    wh_mol  = Chem.MolFromSmiles(wh_smiles)
    ref_mol = Chem.MolFromSmiles(ref_smiles)

    mcs = rdFMCS.FindMCS(
        [wh_mol, ref_mol], timeout=30,
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareOrderExact,
        ringMatchesRingOnly=True,
    )
    n_mcs = mcs.numAtoms
    print(f"  MCS: {n_mcs} atoms ({n_mcs/wh_mol.GetNumAtoms()*100:.0f}% warhead coverage)")

    core = Chem.MolFromSmarts(mcs.smartsString)
    wh_match  = wh_mol.GetSubstructMatch(core)
    ref_match = ref_mol.GetSubstructMatch(core)

    if not wh_match or not ref_match:
        return None

    ref_conf = ref_mol_3d.GetConformer()
    # Collect positions for MCS atoms from the 3D reference ligand
    # ref_mol_3d has no H — its atom order matches ref_mol (no-H SMILES mol)
    coord_map = {}
    for wh_idx, ref_idx in zip(wh_match, ref_match):
        if ref_idx < ref_mol_3d.GetNumAtoms():
            p = ref_conf.GetAtomPosition(ref_idx)
            coord_map[wh_idx] = np.array([p.x, p.y, p.z])

    # Embed the warhead with those constraints
    wh_h = Chem.AddHs(wh_mol)
    # Map heavy-atom indices in wh_mol to corresponding H-mol indices
    hmap = {}
    for i in range(wh_mol.GetNumAtoms()):
        hmap[i] = i   # AddHs preserves heavy-atom ordering

    from rdkit.Geometry import rdGeometry
    coord_map_pt = {k: rdGeometry.Point3D(*v) for k, v in coord_map.items()}

    # Try constrained embed on heavy-atom mol (no H) — simpler and more reliable
    wh_noH = Chem.RWMol(wh_mol)
    AllChem.EmbedMolecule(wh_noH, AllChem.ETKDGv3())
    # Override MCS atom positions with ref coords
    if wh_noH.GetNumConformers() > 0:
        conf = wh_noH.GetConformer()
        for idx, pt in coord_map_pt.items():
            conf.SetAtomPosition(idx, pt)
    else:
        # Build conformer entirely from coord map (place non-MCS atoms at centroid)
        conf = Chem.Conformer(wh_noH.GetNumAtoms())
        mcs_coords = np.array(list(coord_map.values()))
        centroid   = mcs_coords.mean(axis=0)
        for i in range(wh_noH.GetNumAtoms()):
            if i in coord_map:
                conf.SetAtomPosition(i, coord_map[i].tolist())
            else:
                conf.SetAtomPosition(i, centroid.tolist())
        wh_noH.AddConformer(conf, assignId=True)

    AllChem.MMFFOptimizeMolecule(wh_noH.GetMol())
    return wh_noH.GetMol()


def transform_warhead(wh_mol, rot, tran):
    """Map MET-frame warhead coordinates into the target kinase frame.

    BioPython Superimposer convention (confirmed empirically):
      set_atoms(fixed=MET, moving=target) → (rot, tran) s.t. rot.T @ target + tran ≈ MET
      (BioPython stores a right-multiplying matrix; col-vec equivalent: rot.T @ x + tran)
    Inverse (MET → target):  target_coord = rot @ (met_coord - tran)
    """
    conf = wh_mol.GetConformer()
    rw   = RWMol(wh_mol)
    new_conf = Chem.Conformer(rw.GetNumAtoms())
    for i in range(rw.GetNumAtoms()):
        p   = np.array(conf.GetAtomPosition(i))
        new = rot @ (p - tran)
        new_conf.SetAtomPosition(i, new.tolist())
    rw.RemoveAllConformers()
    rw.AddConformer(new_conf, assignId=True)
    return rw.GetMol()


def write_warhead_sdf(wh_mol, out_path):
    """Write warhead molecule (no H) as SDF for use in Protac_params.txt Heads."""
    from rdkit.Chem import SDWriter
    writer = SDWriter(str(out_path))
    writer.write(wh_mol)
    writer.close()


def write_warhead_sdf_with_H(wh_mol, out_path):
    """Write warhead molecule with explicit H as SDF for constraint_generation."""
    from rdkit.Chem import SDWriter, AddHs
    from rdkit.Chem.AllChem import EmbedMolecule, ETKDGv3
    mol_h = AddHs(wh_mol, addCoords=True)
    writer = SDWriter(str(out_path))
    writer.write(mol_h)
    writer.close()


def write_combined_pdb(protein_pdb, wh_mol, out_path, lig_resname="WRH"):
    protein_lines = [l for l in Path(protein_pdb).read_text().splitlines()
                     if l.startswith("ATOM")]
    conf = wh_mol.GetConformer()
    wh_lines = []
    for i, atom in enumerate(wh_mol.GetAtoms()):
        sym = atom.GetSymbol()
        pos = conf.GetAtomPosition(i)
        an  = f"{sym:<2}{i+1:<2}"[:4]
        wh_lines.append(
            f"HETATM{i+1:5d} {an:4s} {lig_resname} L{1:4d}    "
            f"{pos.x:8.3f}{pos.y:8.3f}{pos.z:8.3f}  1.00  0.00          {sym:>2}"
        )
    Path(out_path).write_text(
        "\n".join(protein_lines + wh_lines) + "\nEND\n"
    )


# ─── main ─────────────────────────────────────────────────────────────────────

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--force', nargs='*', metavar='TARGET',
                        help='Force re-generation for specific targets (e.g. --force EPHB2 DDR2). '
                             'If given with no targets, forces ALL.')
    args = parser.parse_args()
    force_set = set(args.force) if args.force is not None else (
        {'ALL'} if args.force is not None else set())
    force_all = args.force is not None and len(args.force) == 0

    print("=" * 60)
    print(" SJF8240 Warhead Placement Pipeline")
    print("=" * 60)

    # ── Step 1: Extract 88Z from 3LQ8 and remap to SJF8240 warhead ───────────
    print(f"\n[1] Extracting 88Z from {REF_PDB.name} and aligning warhead...")
    ref_lig = extract_ligand_as_rdmol(REF_PDB, REF_LIG_NAME)
    if ref_lig is None:
        sys.exit(f"ERROR: {REF_LIG_NAME} not found in {REF_PDB}")
    print(f"  88Z: {ref_lig.GetNumAtoms()} heavy atoms (from PDB coordinates)")

    placed_wh = mcs_remap_coords(SMILES_WH, SMILES_88Z, ref_lig)
    if placed_wh is None:
        sys.exit("ERROR: MCS-based warhead placement failed.")
    print(f"  Placed warhead: {placed_wh.GetNumAtoms()} atoms")

    # Write MET warhead PDB
    met_out = POI_DIR / "high_affinity_degraders" / "MET" / "MET_warhead.pdb"
    write_combined_pdb(REF_PDB, placed_wh, met_out)
    print(f"  -> MET_warhead.pdb")

    # ── Step 2: Load reference Cα + seq for superposition ────────────────────
    print(f"\n[2] Loading reference kinase domain (3LQ8)...")
    ref_ids, ref_ca, ref_seq = extract_ca_and_seq(REF_PDB)
    print(f"  {len(ref_seq)} Cα residues in reference")

    # ── Step 3: Sequence-guided superposition for all other targets ───────────
    print(f"\n[3] Placing warhead in all other POIs...")
    for target, category, pdb_stem in POI_TABLE:
        if target == "MET":
            continue

        out_pdb  = POI_DIR / category / target / f"{target}_warhead.pdb"
        poi_pdb  = POI_DIR / category / target / f"{pdb_stem}.pdb"

        force_this = force_all or target in force_set
        if out_pdb.exists() and not force_this:
            print(f"  [{target}] already exists, skipping.")
            continue
        if not poi_pdb.exists():
            print(f"  [{target}] PDB not found: {poi_pdb.name}")
            continue

        tgt_ids, tgt_ca, tgt_seq = extract_ca_and_seq(poi_pdb)
        if not tgt_ca:
            print(f"  [{target}] no Cα atoms found, skipping.")
            continue

        sup, n_pairs, rmsd = sequence_guided_superposition(
            ref_ca, ref_seq, tgt_ca, tgt_seq
        )
        if sup is None:
            print(f"  [{target}] alignment failed (<10 pairs), skipping.")
            continue

        flag = " *** CHECK ***" if rmsd > 3.0 else ""
        print(f"  [{target}] RMSD={rmsd:.2f} Å over {n_pairs} pairs{flag}")

        rot  = np.array(sup.rotran[0])
        tran = np.array(sup.rotran[1])
        new_wh = transform_warhead(placed_wh, rot, tran)
        write_combined_pdb(poi_pdb, new_wh, out_pdb)
        print(f"    -> {out_pdb.name}")

        # Write SDF files for pRosettaC Heads parameter
        sdf_noH = out_pdb.parent / f"{target}_warhead_placed.sdf"
        sdf_H   = out_pdb.parent / f"{target}_warhead_placed_H.sdf"
        write_warhead_sdf(new_wh, sdf_noH)
        write_warhead_sdf_with_H(new_wh, sdf_H)
        print(f"    -> {sdf_noH.name}  ({new_wh.GetNumAtoms()} atoms)")
        print(f"    -> {sdf_H.name}")

        # Sanity check: anchor atom 1 vs protein
        try:
            conf_chk = new_wh.GetConformer()
            a1 = np.array(conf_chk.GetAtomPosition(0))
            tgt_ca_xyz = np.array([a.get_vector().get_array() for a in tgt_ca])
            min_d = np.sqrt(((tgt_ca_xyz - a1)**2).sum(axis=1)).min()
            flag = ' *** FAR FROM PROTEIN — CHECK CAREFULLY ***' if min_d > 25 else ''
            print(f"    anchor-atom1 to nearest Cα: {min_d:.1f} Å{flag}")
        except Exception as e:
            print(f"    (sanity check error: {e})")

    print("\n" + "=" * 60)
    print("Done. Warhead PDBs written.")
    print("\nCRITICAL: Open all *_warhead.pdb files in PyMOL and verify:")
    print("  1. WRH residue sits inside the ATP-binding pocket")
    print("  2. No major steric clashes with the hinge region")
    print("  3. Targets flagged *** CHECK *** require manual inspection")
    print("     (RMSD > 3 Å indicates divergent kinase fold — may need redocking)")


if __name__ == "__main__":
    main()
