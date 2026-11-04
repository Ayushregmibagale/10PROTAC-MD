"""
Generate a Rosetta/PyRosetta-compatible .params file from an SDF using RDKit.
Lightweight alternative to Rosetta's molfile_to_params.py.

Produces:
  NAME / IO_STRING / TYPE / AA
  ATOM records   (atom_name, rosetta_type, mm_type, charge)
  BOND_TYPE records
  ICOOR_INTERNAL records  (required by PyRosetta to build the residue type)
  NBR_ATOM / NBR_RADIUS
  PDB_ROTAMERS placeholder

Usage:
  python sdf_to_params.py <input.sdf> <residue_name (3 chars)> <output.params>
"""

import sys
import numpy as np
from pathlib import Path
from collections import deque
from rdkit import Chem
from rdkit.Chem import AllChem, rdPartialCharges


# ── atom type mapping ─────────────────────────────────────────────────────────
# Only types that exist in PyRosetta 2026 fa_standard are used.
def rosetta_atom_type(atom, mol):
    sym  = atom.GetSymbol()
    arom = atom.GetIsAromatic()
    nH   = atom.GetTotalNumHs()

    if sym == 'C':
        if arom: return 'aroC'
        if nH == 3: return 'CH3'
        if nH == 2: return 'CH2'
        if nH == 1: return 'CH1'
        for nbr in atom.GetNeighbors():
            bt = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx()).GetBondTypeAsDouble()
            if bt >= 1.5:
                ns = nbr.GetSymbol()
                if ns == 'N': return 'CNH2'
                if ns in ('O', 'S'): return 'COO'
        return 'CH0'

    if sym == 'N':
        if arom: return 'Nhis'
        if nH >= 2: return 'NH2O'
        if nH == 1: return 'Nhis'
        return 'Nlys'

    if sym == 'O':
        if arom: return 'Oaro'
        if nH >= 1: return 'OH'
        for nbr in atom.GetNeighbors():
            bt = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx()).GetBondTypeAsDouble()
            if bt >= 1.5: return 'ONH2'
        return 'Oet2'

    if sym == 'S':
        if nH >= 1: return 'SH1'
        return 'S'

    if sym == 'F':  return 'F'
    if sym == 'Cl': return 'Cl'
    if sym == 'Br': return 'Br'
    if sym == 'I':  return 'I'
    if sym == 'H':  return 'Hpol'
    return 'VIRT'


def atom_name(atom, idx):
    sym = atom.GetSymbol()
    return (sym + str(idx + 1)).ljust(4)[:4]


# ── ICOOR computation ─────────────────────────────────────────────────────────

def _dihedral(p0, p1, p2, p3):
    """Dihedral angle (degrees) defined by atoms p0-p1-p2-p3."""
    b1 = p1 - p0
    b2 = p2 - p1
    b3 = p3 - p2
    n1 = np.cross(b1, b2);  n1 /= (np.linalg.norm(n1) + 1e-12)
    n2 = np.cross(b2, b3);  n2 /= (np.linalg.norm(n2) + 1e-12)
    m1 = np.cross(n1, b2 / (np.linalg.norm(b2) + 1e-12))
    return np.degrees(np.arctan2(np.dot(m1, n2), np.dot(n1, n2)))


def _angle(p0, p1, p2):
    """Bond angle (degrees) p0-p1-p2."""
    v1 = p0 - p1;  v1 /= (np.linalg.norm(v1) + 1e-12)
    v2 = p2 - p1;  v2 /= (np.linalg.norm(v2) + 1e-12)
    return np.degrees(np.arccos(np.clip(np.dot(v1, v2), -1, 1)))


def _build_spanning_tree(mol, root):
    """BFS spanning tree; returns (order, parent, grandparent, great_grandparent)."""
    n = mol.GetNumAtoms()
    adj = {i: [nb.GetIdx() for nb in mol.GetAtomWithIdx(i).GetNeighbors()]
           for i in range(n)}

    visited = {root}
    parent = {root: root}
    gp     = {root: root}
    ggp    = {root: root}
    order  = [root]
    queue  = deque([root])

    while queue:
        a = queue.popleft()
        for nb in adj[a]:
            if nb not in visited:
                visited.add(nb)
                parent[nb] = a
                gp[nb]     = parent[a]
                ggp[nb]    = gp[a]
                order.append(nb)
                queue.append(nb)

    return order, parent, gp, ggp


def generate_icoor_lines(mol_noH, atom_names, nbr_idx):
    """
    Compute ICOOR_INTERNAL records from 3D coordinates.
    Rosetta convention:
      phi   = dihedral(<atom>, <parent>, <angle_atom>, <torsion_atom>)
      theta = 180 - bond_angle(<angle_atom>, <parent>, <atom>)
      d     = bond_length(<parent>, <atom>)
    For the root atom all values are 0.
    """
    if mol_noH.GetNumConformers() == 0:
        return []

    conf = mol_noH.GetConformer()
    pos  = {i: np.array([conf.GetAtomPosition(i).x,
                          conf.GetAtomPosition(i).y,
                          conf.GetAtomPosition(i).z])
            for i in range(mol_noH.GetNumAtoms())}

    order, par, gpar, ggpar = _build_spanning_tree(mol_noH, nbr_idx)

    lines = []
    for rank, idx in enumerate(order):
        name     = atom_names[idx].strip()
        p_idx    = par[idx]
        gp_idx   = gpar[idx]
        ggp_idx  = ggpar[idx]
        p_name   = atom_names[p_idx].strip()
        gp_name  = atom_names[gp_idx].strip()
        ggp_name = atom_names[ggp_idx].strip()

        if rank == 0:
            phi, theta, d = 0.0, 0.0, 0.0
        elif rank == 1:
            d     = float(np.linalg.norm(pos[idx] - pos[p_idx]))
            theta = 0.0   # no angle_atom yet
            phi   = 0.0
        elif rank == 2:
            d     = float(np.linalg.norm(pos[idx] - pos[p_idx]))
            theta = 180.0 - _angle(pos[gp_idx], pos[p_idx], pos[idx])
            phi   = 0.0
        else:
            d     = float(np.linalg.norm(pos[idx] - pos[p_idx]))
            theta = 180.0 - _angle(pos[gp_idx], pos[p_idx], pos[idx])
            phi   = _dihedral(pos[idx], pos[p_idx], pos[gp_idx], pos[ggp_idx])

        lines.append(
            f"ICOOR_INTERNAL  {name:<4s}  {phi:12.6f}  {theta:12.6f}  {d:12.6f}"
            f"  {p_name:<4s}  {gp_name:<4s}  {ggp_name:<4s}"
        )

    return lines


# ── main generator ────────────────────────────────────────────────────────────

def generate_params(sdf_path: Path, res_name: str, out_path: Path):
    supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
    mol = next(supplier)
    if mol is None:
        mol = Chem.MolFromMolFile(str(sdf_path), removeHs=False)
    assert mol, f"Could not read {sdf_path}"

    mol_noH = Chem.RemoveHs(mol)
    AllChem.ComputeGasteigerCharges(mol_noH)

    lines = []
    lines.append(f"NAME {res_name}")
    lines.append(f"IO_STRING {res_name} X")
    lines.append(f"TYPE LIGAND")
    lines.append(f"AA UNK")

    # ATOM records
    atom_names = []
    for atom in mol_noH.GetAtoms():
        idx    = atom.GetIdx()
        name   = atom_name(atom, idx)
        atom_names.append(name)
        rtype  = rosetta_atom_type(atom, mol_noH)
        charge = float(atom.GetPropsAsDict().get("_GasteigerCharge", 0.0))
        if charge != charge:
            charge = 0.0
        lines.append(f"ATOM {name} {rtype} X {charge:6.3f}")

    # BOND_TYPE records
    for bond in mol_noH.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        btype = {
            Chem.BondType.SINGLE:   "1",
            Chem.BondType.DOUBLE:   "2",
            Chem.BondType.TRIPLE:   "3",
            Chem.BondType.AROMATIC: "ar",
        }.get(bond.GetBondType(), "1")
        lines.append(f"BOND_TYPE {atom_names[i]:4s} {atom_names[j]:4s} {btype}")

    # NBR_ATOM: heavy atom closest to geometric centre
    if mol_noH.GetNumConformers() > 0:
        conf = mol_noH.GetConformer()
        pos  = conf.GetPositions()
        centre = pos.mean(axis=0)
        dists  = np.linalg.norm(pos - centre, axis=1)
        nbr_idx    = int(dists.argmin())
        nbr_radius = float(dists.max()) + 1.0
    else:
        nbr_idx    = mol_noH.GetNumAtoms() // 2
        nbr_radius = 10.0

    lines.append(f"NBR_ATOM {atom_names[nbr_idx]:4s}")
    lines.append(f"NBR_RADIUS {nbr_radius:.3f}")

    # ICOOR_INTERNAL records — required by PyRosetta to build the residue type
    icoor = generate_icoor_lines(mol_noH, atom_names, nbr_idx)
    lines.extend(icoor)

    # PDB_ROTAMERS is omitted here; add it in callers that generate multi-conformer SDFs.

    out_path.write_text("\n".join(lines) + "\n")
    print(f"Params written: {out_path}")
    print(f"  {mol_noH.GetNumAtoms()} atoms, {mol_noH.GetNumBonds()} bonds, "
          f"{len(icoor)} ICOOR records")
    print(f"  NBR_ATOM: {atom_names[nbr_idx]}, NBR_RADIUS: {nbr_radius:.2f} Å")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python sdf_to_params.py <input.sdf> <RES_NAME> <output.params>")
        sys.exit(1)
    generate_params(Path(sys.argv[1]), sys.argv[2].upper(), Path(sys.argv[3]))
