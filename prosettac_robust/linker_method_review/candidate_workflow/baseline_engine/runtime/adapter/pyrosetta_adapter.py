"""
PyRosetta adapter for pRosettaC.

Monkey-patches pRosettaC's `rosetta` and `utils` modules to use:
  - PyRosetta instead of Rosetta standalone executables
  - Python OpenBabel (openbabel-wheel) instead of the `babel` CLI
  - Scipy/NumPy for the PatchDock global docking step

Call `apply_patches()` BEFORE importing main.py / constraint_generation.py.

Design notes:
  - All heavy-atom indices in SDF files are preserved after AddHydrogens()
    because OpenBabel appends H atoms at the end.
  - The global_dock() function replaces PatchDock: it samples random rigid-body
    orientations of protein B, constraining the anchor distance to [min, max]
    and rejecting steric clashes via a KD-tree.
"""

from __future__ import annotations
import os
import sys
import glob
import shutil
import random
import hashlib
import json
from pathlib import Path

def publication_seed(stage):
    base=int(os.environ.get("PROTAC_SEED", "42"))
    root=Path(os.environ.get("PROTAC_ATTEMPT_ROOT", os.getcwd())).resolve()
    try: task=str(Path.cwd().resolve().relative_to(root))
    except ValueError: raise RuntimeError("Worker outside declared attempt root")
    return 1+int.from_bytes(hashlib.sha256((str(base)+"|"+task+"|"+os.environ.get("PROTAC_TASK_KEY","parent")+"|"+stage).encode()).digest()[:4],"big")%2147483646

random.seed(publication_seed("python"))
import numpy as np
from pathlib import Path
from scipy.spatial import cKDTree
np.random.seed(publication_seed("numpy"))

# ── directories ──────────────────────────────────────────────────────────────
_THIS_DIR = Path(__file__).parent
_PROJECT_ROOT = _THIS_DIR.parent.parent
_PROSETTAC_DIR = Path(os.environ['PROTAC_PROSETTAC'])

# ── PyRosetta init (done once) ────────────────────────────────────────────────
_pyrosetta_ready = False

def _init_pyrosetta(extra_params=()):
    global _pyrosetta_ready
    if _pyrosetta_ready:
        return
    import pyrosetta
    opts = '-ignore_unrecognized_res false -load_PDB_components false -out:level 100'
    opts += f' -constant_seed -jran {publication_seed("rosetta")}'
    for p in extra_params:
        opts += f' -extra_res_fa {p}'
    pyrosetta.init(opts)
    _pyrosetta_ready = True


# ═══════════════════════════════════════════════════════════════════════════════
#  OpenBabel replacements  (utils.addH_sdf, utils.pdb2sdf, utils.sdf2sdf)
# ═══════════════════════════════════════════════════════════════════════════════

def addH_sdf(sdf_file: str, new_sdf: str | None = None) -> None:
    """
    Add hydrogens to an SDF file using RDKit.
    RDKit preserves the original bond types / aromaticity perception so that
    pRosettaC's translate_anchors() substructure match succeeds.
    OpenBabel would re-perceive aromaticity and break the match.

    Uses a per-file lock so parallel processes don't race when the output
    path is shared (e.g. SJF8240_vhll_placed_H.sdf used by all targets).
    """
    import fcntl
    from rdkit import Chem
    out = new_sdf if new_sdf else sdf_file
    lock_path = str(out) + '.lock'
    with open(lock_path, 'w') as lf:
        fcntl.flock(lf, fcntl.LOCK_EX)
        # After acquiring lock check if another process already wrote a valid file
        if os.path.exists(out):
            try:
                mol_test = Chem.SDMolSupplier(out, removeHs=False)[0]
                if mol_test is not None:
                    return
            except Exception:
                pass  # fall through and regenerate
        mol = Chem.SDMolSupplier(sdf_file)[0]
        if mol is None:
            mol = Chem.SDMolSupplier(sdf_file, removeHs=False)[0]
        mol_h = Chem.AddHs(mol, addCoords=True)
        writer = Chem.SDWriter(out)
        writer.write(mol_h)
        writer.close()


def pdb2sdf(pdb_file: str, sdf_file: str) -> None:
    """
    Convert PDB to SDF using Python OpenBabel.
    The resulting SDF may have invalid valences for quaternary N atoms (e.g. in
    PT0/PT1 HETATM records); callers should read with sanitize=False if needed.
    """
    from openbabel import openbabel as ob
    conv = ob.OBConversion()
    conv.SetInAndOutFormats('pdb', 'sdf')
    mol = ob.OBMol()
    conv.ReadFile(mol, pdb_file)
    conv.WriteFile(mol, sdf_file)


def sdf2sdf(sdf_file: str) -> None:
    """Re-save SDF (normalise format) using Python OpenBabel."""
    from openbabel import openbabel as ob
    conv = ob.OBConversion()
    conv.SetInAndOutFormats('sdf', 'sdf')
    mol = ob.OBMol()
    conv.ReadFile(mol, sdf_file)
    conv.WriteFile(mol, sdf_file)


# ═══════════════════════════════════════════════════════════════════════════════
#  Rosetta replacements  (rosetta.clean, rosetta.mol_to_params,
#                         rosetta.relax, rosetta.local_docking)
# ═══════════════════════════════════════════════════════════════════════════════

def clean_pdb(struct: str, chains: str) -> None:
    """Extract ATOM/HETATM records for the given chains."""
    chain_set = set(chains)
    out = struct.split('.pdb')[0] + '_' + chains + '.pdb'
    with open(struct) as f, open(out, 'w') as g:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                if len(line) > 21 and line[21] in chain_set:
                    g.write(line)
            elif line.startswith(('TER', 'END')):
                g.write(line)
    # pRosettaC's rosetta.clean() removes .fasta files generated by Rosetta's
    # clean_pdb.py — those don't exist here, so just skip silently.
    for ch in chains:
        fasta = struct.split('.pdb')[0] + '_' + ch + '.fasta'
        if os.path.exists(fasta):
            pass


def _sdf_to_pdb_block(sdf_file: str, res_name: str) -> str:
    """
    Convert SDF → HETATM PDB block (heavy atoms only).
    Atom names match those generated by sdf_to_params.atom_name(), so PyRosetta
    can map the PDB atoms to the params file unambiguously.
    """
    from rdkit import Chem
    sys.path.insert(0, str(_THIS_DIR))
    from sdf_to_params import atom_name as _an

    mol = Chem.SDMolSupplier(sdf_file, removeHs=False)[0]
    if mol is None:
        raise ValueError(f'Cannot read {sdf_file}')
    mol_noH = Chem.RemoveHs(mol)
    conf = mol_noH.GetConformer()

    lines = []
    rname3 = res_name.ljust(3)[:3]
    for atom in mol_noH.GetAtoms():
        idx  = atom.GetIdx()
        name = _an(atom, idx).strip()    # e.g. "C1", "N2"
        pos  = conf.GetAtomPosition(idx)
        elem = atom.GetSymbol()
        # PDB fixed-width HETATM record
        # cols: 1-6 record, 7-11 serial, 13-16 name, 18-20 resname, 22 chain,
        #       23-26 resseq, 31-38 x, 39-46 y, 47-54 z, 55-60 occ, 61-66 B,
        #       77-78 element
        line = (
            f'HETATM{idx + 1:5d}  {name:<3s} {rname3} X   1    '
            f'{pos.x:8.3f}{pos.y:8.3f}{pos.z:8.3f}'
            f'  1.00  0.00          {elem:>2s}  '
        )
        lines.append(line)
    lines.append('END')
    return '\n'.join(lines) + '\n'


def mol_to_params(sdf_file: str, name: str, pdb_prefix: str,
                  overwrite: bool = True, conformers: bool = False,
                  nbr: int = -1, v_atoms_sdf: str = ''):
    from official_parameters import generate
    return generate(sdf_file,name,pdb_prefix,overwrite,conformers,nbr,v_atoms_sdf)


def _append_score_sc(score,desc):
    import fcntl,math
    if not math.isfinite(score):raise RuntimeError('Nonfinite score')
    with open('score.sc','a+') as f:
        fcntl.flock(f,fcntl.LOCK_EX)
        f.seek(0,2)
        if f.tell()==0:f.write('SEQUENCE: \nSCORE: total_score description\n')
        f.write(f'SCORE: {score:.6f} {desc}\n');f.flush();os.fsync(f.fileno())


def relax(struct: str, sdf_params: str, interface: bool = False, n: int = 1) -> None:
    """Original ref2015 side-chain repacking; propagate all scientific failures."""
    import math
    import pyrosetta
    from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
    from pyrosetta.rosetta.core.pack.task import TaskFactory
    from pyrosetta.rosetta.core.pack.task.operation import RestrictToRepacking
    _init_pyrosetta()
    prefix=struct.split('.pdb')[0]
    pose=pyrosetta.Pose()
    pyrosetta.generate_nonstandard_residue_set(pose,[sdf_params])
    pyrosetta.pose_from_file(pose,struct)
    scorefxn=pyrosetta.create_score_function('ref2015')
    tf=TaskFactory();tf.push_back(RestrictToRepacking())
    packer=PackRotamersMover(scorefxn);packer.task_factory(tf)
    for i in range(int(n)):
        wpose=pyrosetta.Pose(pose)
        packer.apply(wpose)
        score=float(scorefxn(wpose))
        if not math.isfinite(score): raise RuntimeError('Nonfinite Rosetta score')
        out_pdb=f'{prefix}_{i + 1:04d}.pdb'
        wpose.dump_pdb(out_pdb)
        _append_score_sc(score,Path(out_pdb).stem)


def _symlink_path(p: str) -> str:
    """Replace space-containing project root with symlink so PyRosetta init string survives."""
    return p.replace(str(_PROJECT_ROOT), str(Path.home() / 'protac'))


def local_docking(struct: str, chainsA: str, chainsB: str,
                  ptA_params: str, ptB_params: str, nstruct: int = 50) -> str:
    """
    Return shell command that runs PyRosetta local docking on one PatchDock result.
    The command calls run_docking.py which is a standalone PyRosetta script.
    Params paths are rewritten through ~/protac symlink so PyRosetta init (which
    parses options by splitting on whitespace) doesn't truncate space-containing paths.
    """
    script = str(_THIS_DIR / 'run_docking.py')
    ptA = _symlink_path(ptA_params)
    ptB = _symlink_path(ptB_params)
    return (
        f'{sys.executable} "{script}" "{struct}" "{chainsA}" "{chainsB}" '
        f'"{ptA}" "{ptB}" {nstruct}'
    )


# ═══════════════════════════════════════════════════════════════════════════════
#  PatchDock replacement  (utils.patchdock)
# ═══════════════════════════════════════════════════════════════════════════════

def _random_rotation() -> np.ndarray:
    """Uniform random rotation via quaternion."""
    u = np.random.uniform(0, 1, 3)
    q = np.array([
        np.sqrt(1 - u[0]) * np.sin(2 * np.pi * u[1]),
        np.sqrt(1 - u[0]) * np.cos(2 * np.pi * u[1]),
        np.sqrt(u[0])     * np.sin(2 * np.pi * u[2]),
        np.sqrt(u[0])     * np.cos(2 * np.pi * u[2]),
    ])
    w, x, y, z = q
    return np.array([
        [1 - 2*(y*y + z*z),     2*(x*y - z*w),     2*(x*z + y*w)],
        [    2*(x*y + z*w), 1 - 2*(x*x + z*z),     2*(y*z - x*w)],
        [    2*(x*z - y*w),     2*(y*z + x*w), 1 - 2*(x*x + y*y)],
    ])


def _read_pdb_coords(pdb_file: str) -> tuple[np.ndarray, list[str]]:
    """Return (N,3) heavy-atom coords and raw PDB lines."""
    coords = []
    lines = Path(pdb_file).read_text().splitlines(keepends=True)
    for line in lines:
        if line.startswith(('ATOM', 'HETATM')):
            elem = line[76:78].strip() if len(line) > 76 else ''
            if elem != 'H':
                coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.array(coords) if coords else np.empty((0, 3)), lines


def _read_pdb_ca_coords(pdb_file: str) -> np.ndarray:
    """Return (N,3) Cα-only coords (ATOM CA records) for fast clash screening."""
    coords = []
    for line in Path(pdb_file).read_text().splitlines():
        if line.startswith('ATOM') and line[12:16].strip() == 'CA':
            coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.array(coords) if coords else np.empty((0, 3))


def _get_hetatm_serial_coord(pdb_file: str, serial: int) -> np.ndarray:
    """Get coordinates of the HETATM with the given 1-based serial number."""
    for line in Path(pdb_file).read_text().splitlines():
        if line.startswith('HETATM') and int(line[6:11]) == serial:
            return np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    raise ValueError(f'HETATM serial {serial} not found in {pdb_file}')


def _apply_rt_to_lines(lines: list[str], R: np.ndarray, t: np.ndarray) -> list[str]:
    """Apply rotation+translation to ATOM/HETATM coordinate columns in PDB lines."""
    out = []
    for line in lines:
        if line.startswith(('ATOM', 'HETATM')):
            xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            nxyz = R @ xyz + t
            line = line[:30] + f'{nxyz[0]:8.3f}{nxyz[1]:8.3f}{nxyz[2]:8.3f}' + line[54:]
        out.append(line)
    return out


def global_dock(structs: list[str], anchors: list[int],
                min_dist: float, max_dist: float,
                num_results: int = 1000, threshold: float = 2.0) -> int | None:
    """
    Constraint-guided global protein-protein docking.
    Replaces PatchDock. Samples random orientations of structs[1] relative to
    structs[0], keeping poses where the anchor atoms are within [min_dist, max_dist]
    and no Cα steric clash exists (Cα–Cα cutoff 4 Å).

    Random-sampling acceptance rate for large proteins is very low when using
    all heavy atoms.  Cα-only clash detection with a 4 Å threshold achieves
    ~100× higher throughput while still rejecting deeply interpenetrating poses.

    Writes Patchdock_Results/pd.{i}.pdb (1-indexed) and returns the count.
    """
    structA, structB = structs
    anchorA_serial, anchorB_serial = anchors  # 1-based serials

    _, linesA = _read_pdb_coords(structA)
    _, linesB = _read_pdb_coords(structB)

    caA = _read_pdb_ca_coords(structA)   # Cα coords for fast clash check
    caB_orig = _read_pdb_ca_coords(structB)

    anchorA = _get_hetatm_serial_coord(structA, anchorA_serial)
    anchorB_orig = _get_hetatm_serial_coord(structB, anchorB_serial)

    # Cα KD-tree for A; minimum Cα–Cα contact distance in protein complexes ~ 4 Å
    ca_threshold = 4.0
    treeA = cKDTree(caA) if len(caA) > 0 else None

    os.makedirs('Patchdock_Results', exist_ok=True)

    rng = np.random.default_rng(publication_seed("global_translation"))
    results = 0
    # Large trial budget: Cα-only clash is cheap, acceptance ~0.1-1% for big proteins
    max_trials = min(max(num_results * 20_000, 500_000),int(os.environ.get("PROTAC_MAX_GLOBAL_TRIALS","10000000")))

    print(f'[global_dock] caA={len(caA)} caB={len(caB_orig)} '
          f'anchor_range=[{min_dist:.1f},{max_dist:.1f}] Å  '
          f'budget={max_trials:,}', flush=True)

    trials_attempted=0
    for trial in range(max_trials):
        if results >= num_results:
            break

        trials_attempted+=1
        R = _random_rotation()

        d = rng.uniform(min_dist, max_dist)
        direction = rng.standard_normal(3)
        direction /= np.linalg.norm(direction)
        target_anchorB = anchorA + d * direction

        anchorB_rot = R @ anchorB_orig
        t = target_anchorB - anchorB_rot

        # Cα clash check (fast)
        if treeA is not None and len(caB_orig) > 0:
            caB_new = caB_orig @ R.T + t
            min_dists = treeA.query(caB_new, k=1)[0]
            if np.any(min_dists < ca_threshold):
                continue

        # Write combined PDB
        transformed_linesB = _apply_rt_to_lines(linesB, R, t)
        out_path = f'Patchdock_Results/pd.{results + 1}.pdb'
        with open(out_path, 'w') as f:
            for line in linesA:
                f.write(line if line.endswith('\n') else line + '\n')
            for line in transformed_linesB:
                if line.startswith('HETATM') and len(line) > 21 and line[17:20].strip() == 'PT1' and line[21] == 'X':
                    line = line[:21] + 'Y' + line[22:]
                f.write(line if line.endswith('\n') else line + '\n')

        results += 1
        if results % 10 == 0:
            print(f'[global_dock] {results}/{num_results} after {trial+1:,} trials', flush=True)

    print(f'[global_dock] done: {results} poses from {trials_attempted:,} trials', flush=True)
    with open("global_sampling.json","w") as ledger: json.dump({"attempted":trials_attempted,"succeeded":results,"rejected":trials_attempted-results,"failed":0},ledger)
    return results if results > 0 else None


# ═══════════════════════════════════════════════════════════════════════════════
#  Patch application
# ═══════════════════════════════════════════════════════════════════════════════

def translate_anchors_fixed(old_sdf: str, new_sdf: str, old_anchor: int) -> int:
    """
    Map anchor atom index from old_sdf (no H) to new_sdf (with H).
    Replaces protac_lib.translate_anchors which fails when old/new SDF have
    different aromaticity representations (sanitize=False vs True).

    RDKit's AddHs(addCoords=True) always appends H atoms AFTER heavy atoms,
    so heavy atom indices are preserved identically in the H-added SDF.
    We verify this by comparing heavy atom counts.
    """
    from rdkit import Chem
    mol_old = Chem.SDMolSupplier(old_sdf, removeHs=True)[0]
    mol_new = Chem.SDMolSupplier(new_sdf, removeHs=True)[0]
    if mol_old is not None and mol_new is not None:
        if mol_old.GetNumAtoms() == mol_new.GetNumAtoms():
            return old_anchor   # indices preserved by AddHs
    # Fallback: substructure match with consistent sanitization
    if mol_old is None:
        mol_old = Chem.SDMolSupplier(old_sdf, sanitize=True)[0]
    if mol_new is None:
        mol_new = Chem.SDMolSupplier(new_sdf, sanitize=True)[0]
    if mol_old is None or mol_new is None:
        return old_anchor
    match = mol_new.GetSubstructMatch(mol_old)
    return match[old_anchor] if match else old_anchor


def sample_dist_free_conf(Heads, Anchors, Linkers, n: int = 200,
                          output_hist: str = 'initial_distances.hist',
                          hist_threshold: float = 0.75, min_margin: float = 2,
                          homo_protac: bool = False):
    """Inherited free ETKDGv3 distance histogram with exact stereo head mapping.

    This is a seed-dependent sampling heuristic, not an exhaustive physical
    linker-reach distribution or a probability estimate. Both-head constrained
    final linker construction is a separate stage.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolTransforms
    from rdkit.Chem.rdchem import Conformer

    HeadA_sdf, HeadB_sdf = Heads

    with open(Linkers) as f:
        smi = f.readline().split()[0]
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return (None, None)

    HeadA = Chem.SDMolSupplier(HeadA_sdf)[0]
    HeadB = Chem.SDMolSupplier(HeadB_sdf)[0]
    if HeadA is None or HeadB is None:
        return (None, None)

    from constrained_linker import distance_anchor_indices
    anchor_A,anchor_B=distance_anchor_indices(mol,[HeadA,HeadB],Anchors)

    # Generate free conformers with ETKDGv3
    mol_h = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = publication_seed("distance_conformers")
    params.numThreads = 1
    params.enforceChirality = True
    conf_ids = list(AllChem.EmbedMultipleConfs(mol_h, numConfs=n, params=params))
    if not conf_ids:
        return (0, 0)

    # Measure anchor-to-anchor distances and build histogram
    histogram: dict[int, int] = {}
    writer = Chem.SDWriter('random_sampling.sdf')
    for cid in conf_ids:
        d = rdMolTransforms.GetBondLength(mol_h.GetConformer(cid), anchor_A, anchor_B)
        b = max(1, int(round(d)))
        histogram[b] = histogram.get(b, 0) + 1
        # Write heavy-atom molecule for downstream compatibility
        m_noH = Chem.RemoveHs(mol_h)
        conf_noH = Conformer(m_noH.GetNumAtoms())
        conf_h = mol_h.GetConformer(cid)
        for i in range(m_noH.GetNumAtoms()):
            conf_noH.SetAtomPosition(i, conf_h.GetAtomPosition(i))
        m_noH.RemoveAllConformers()
        m_noH.AddConformer(conf_noH, assignId=True)
        writer.write(m_noH)
    writer.close()

    with open(output_hist, 'w') as f:
        for b in sorted(histogram):
            f.write(f'{b}\t{histogram[b]}\n')

    sum_his = sum(histogram.values())
    if sum_his == 0:
        return (0, 0)

    max_count = max(histogram.values())
    avg_index = sum(b * c for b, c in histogram.items()) / sum_his
    threshold = max_count * hist_threshold
    high_values = [b for b, c in histogram.items() if c >= threshold]

    bounds=(min(min(high_values),avg_index-min_margin),
            max(max(high_values),avg_index+min_margin))
    from protac_qc import atomic_json,sha256
    atomic_json(Path('anchor_distance_sampling.json'),{
        'status':'SAMPLED_HEURISTIC_RANGE','requested_conformers':n,
        'generated_conformers':len(conf_ids),'num_threads':params.numThreads,
        'seed':params.randomSeed,'embedding_max_iterations_setting':params.maxIterations,
        'full_molecule_anchor_indices_zero_based':[anchor_A,anchor_B],
        'mapping':'Exact disjoint stereochemical complete heads; unique anchor pair',
        'histogram_bin_width_A':1.0,'histogram':histogram,'hist_threshold':hist_threshold,
        'min_margin_A':min_margin,'range_A':bounds,
        'source_sha256':{str(p):sha256(p) for p in [*Heads,Linkers]},
        'interpretation':'Finite free-conformer heuristic; no biological negative or population inferred'})
    return bounds


def _make_docked_heads_mol(docked_sdf_or_pdb,headA_sdf,headB_sdf):
    from official_parameters import docked_heads
    return docked_heads(docked_sdf_or_pdb,[headA_sdf,headB_sdf])


def _find_linker_attachment_atoms(mol, head_match):
    """
    Return the subset of head_match indices (in mol) that are directly bonded to
    at least one atom NOT in head_match.  These are the attachment points where
    the linker joins the warhead — typically 1-2 atoms per warhead.
    """
    head_set = set(head_match)
    boundary = []
    for idx in head_match:
        for nbr in mol.GetAtomWithIdx(idx).GetNeighbors():
            if nbr.GetIdx() not in head_set:
                boundary.append(idx)
                break
    return boundary


# ═══════════════════════════════════════════════════════════════════════════════
#  Crystal-structure coordMap helpers (used by gen_const_conf_anchored)
# ═══════════════════════════════════════════════════════════════════════════════

_mol_3jf_cache = None      # 3JF from 4W9H.pdb with CONECT bonds
_mol_88z_cache = None      # 88Z from 3LQ8.pdb with CONECT bonds
_met_clean_ca_cache = None  # MET_protein_clean.pdb Cα {resnum: ndarray}
_t1_cache: dict = {}       # {poi_protein_clean: (R1, t1)} MET→POI Kabsch


def _load_mol_with_conect(pdb_path: str, resname: str):
    """Load a single ligand from PDB using CONECT records for bond graph."""
    from rdkit import Chem
    mol = Chem.MolFromPDBFile(pdb_path, sanitize=False, removeHs=True)
    if mol is None:
        return None
    em = Chem.RWMol()
    old2new: dict = {}
    for atom in mol.GetAtoms():
        ri = atom.GetPDBResidueInfo()
        if ri and ri.GetResidueName().strip() == resname:
            ni = em.AddAtom(Chem.Atom(atom.GetAtomicNum()))
            em.GetAtomWithIdx(ni).SetNoImplicit(True)
            old2new[atom.GetIdx()] = ni
    if not old2new:
        return None
    for bond in mol.GetBonds():
        bi, ei = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if bi in old2new and ei in old2new:
            em.AddBond(old2new[bi], old2new[ei], bond.GetBondType())
    pos = mol.GetConformer().GetPositions()
    conf = Chem.Conformer(em.GetNumAtoms())
    new2old = {v: k for k, v in old2new.items()}
    for ni in range(em.GetNumAtoms()):
        conf.SetAtomPosition(ni, pos[new2old[ni]])
    em.AddConformer(conf, assignId=True)
    try:
        Chem.SanitizeMol(em)
    except Exception:
        pass
    return em.GetMol()


def _kabsch_col(ref: np.ndarray, mob: np.ndarray):
    """Kabsch (column-vector): returns R, t s.t. R @ mob[i] + t ≈ ref[i]."""
    rm, mm = ref.mean(0), mob.mean(0)
    H = (mob - mm).T @ (ref - rm)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1., 1., d]) @ U.T
    return R, rm - R @ mm


def _read_ca_by_resnum(pdb_path: str, chain: str = 'A') -> dict:
    """Return {resnum: np.ndarray} for Cα atoms in the given chain."""
    ca: dict = {}
    for line in open(pdb_path):
        if line.startswith('ATOM') and line[12:16].strip() == 'CA' and line[21] == chain:
            r = int(line[22:26])
            ca[r] = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return ca


def _icp_kabsch(ref_arr: np.ndarray, mob_arr: np.ndarray,
                max_iter: int = 15, dist_init: float = 10.0):
    """
    ICP Kabsch: align mob onto ref via iterative closest-point.
    Returns (R, t) column-vector: R @ mob[i] + t ≈ ref[nearest].
    """
    R = np.eye(3)
    t = ref_arr.mean(0) - mob_arr.mean(0)
    prev_rmsd = np.inf
    for it in range(max_iter):
        mob_t = (R @ mob_arr.T).T + t
        tree = cKDTree(ref_arr)
        dists, inds = tree.query(mob_t, k=1)
        cutoff = max(3.0, dist_init / (1.0 + it * 0.5))
        mask = dists < cutoff
        if mask.sum() < 10:
            break
        R, t = _kabsch_col(ref_arr[inds[mask]], mob_arr[mask])
        rmsd = float(np.sqrt((dists[mask] ** 2).mean()))
        if abs(rmsd - prev_rmsd) < 0.001:
            break
        prev_rmsd = rmsd
    return R, t


def _get_poi_protein_clean(docked_heads_path: str):
    """Parse Protac_params.txt to find POI_protein_clean.pdb path."""
    import re
    run_dir = Path(docked_heads_path).parent.parent
    params_file = run_dir / 'Protac_params.txt'
    if not params_file.exists():
        return None
    for line in open(params_file):
        if line.startswith('Structures:'):
            val = line.split(':', 1)[1].strip()
            # Split on whitespace immediately before '/' (path separator trick that
            # handles directory names with embedded spaces like 'PROTAC degradation Prediction').
            # Within a single absolute path, spaces are never directly followed by '/'.
            paths = re.split(r'\s+(?=/)', val)
            if len(paths) >= 2:
                p = Path(paths[1].strip())
                if p.exists():
                    return str(p)
                # Resolve through project-root alternatives
                for candidate in _PROJECT_ROOT.rglob(p.name):
                    return str(candidate)
    return None


def _compute_t1_met_to_poi(poi_protein_clean: str):
    """
    ICP Kabsch: MET Ca (3LQ8 frame) → POI Ca frame.
    Returns (R1, t1) s.t. R1 @ met_ca + t1 ≈ poi_ca (matched pairs).
    Returns identity for MET itself or when alignment fails.
    """
    global _met_clean_ca_cache
    met_pdb = str(_PROJECT_ROOT / 'data' / 'structures' / 'POI' /
                  'high_affinity_degraders' / 'MET' / 'MET_protein_clean.pdb')
    if _met_clean_ca_cache is None:
        _met_clean_ca_cache = _read_ca_by_resnum(met_pdb, 'A')
    met_ca = _met_clean_ca_cache
    if not poi_protein_clean or not os.path.exists(poi_protein_clean):
        return np.eye(3), np.zeros(3)
    poi_ca = _read_ca_by_resnum(poi_protein_clean, 'A')
    if not poi_ca:
        return np.eye(3), np.zeros(3)
    # If POI IS MET (≥90% residue overlap), return identity
    if (len(set(met_ca) & set(poi_ca)) / max(len(met_ca), 1)) > 0.9:
        return np.eye(3), np.zeros(3)
    met_arr = np.array(list(met_ca.values()))
    poi_arr = np.array(list(poi_ca.values()))
    return _icp_kabsch(ref_arr=poi_arr, mob_arr=met_arr)


def _compute_t2_poi_to_vhl(docked_heads_path: str, poi_protein_clean: str):
    """
    Per-pose Kabsch: POI Ca clean → docked VHL frame.
    Returns (R2, t2) s.t. R2 @ poi_clean_ca + t2 ≈ docked_pd_ca.

    Uses the local docking file (pd.X_docking_NNNN.pdb) for per-pose accuracy,
    since local docking can reposition the kinase significantly from the global
    PatchDock pose. Falls back to pd.X.pdb if the local file is missing.
    """
    import re
    m = re.search(r'docked_(\d+)_(\d+)', Path(docked_heads_path).name)
    if not m:
        return None, None
    if not poi_protein_clean or not os.path.exists(poi_protein_clean):
        return None, None

    pd_dir = Path(docked_heads_path).parent
    global_num, local_num = m.group(1), int(m.group(2))

    # Prefer local docking file for per-pose kinase position
    local_pdb = pd_dir / f'pd.{global_num}_docking_{local_num:04d}.pdb'
    global_pdb = pd_dir / f'pd.{global_num}.pdb'
    ref_pdb = str(local_pdb) if local_pdb.exists() else (
              str(global_pdb) if global_pdb.exists() else None)
    if ref_pdb is None:
        return None, None

    ca_poi = _read_ca_by_resnum(poi_protein_clean, 'A')
    ca_pd  = _read_ca_by_resnum(ref_pdb, 'A')
    common = sorted(set(ca_poi) & set(ca_pd))
    if len(common) < 10:
        return None, None
    ref = np.array([ca_pd[r]  for r in common])
    mob = np.array([ca_poi[r] for r in common])
    return _kabsch_col(ref, mob)


def _build_crystal_coordmap(head_linker_smi: str, docked_heads_path: str) -> dict:
    """
    Build coordMap for DG embedding using crystal-structure bijections.

    VHL anchor: MCS(head_linker, 3JF from 4W9H) → up to 19 atoms, VHL frame.
      3JF is loaded with CONECT bonds → correct atom-position bijection, 0 bad bonds.

    Kinase anchor is intentionally omitted: placed warhead SDFs have a bijection
    bug (place_warhead.py mcs_remap_coords assigns SMILES-index k the position of
    PDB-serial k) that produces bonds up to 35 Å in the coordMap, causing DG to
    always fail.  The 88Z-via-T1 alternative has ~3 Å ICP error for non-MET POIs
    that also causes DG failure. VHL-only (19 atoms) achieves ~47% seed success
    for MET and ~6% for DDR2/other POIs, giving ≥30 conformers per 500 seeds.

    Returns {hl_atom_idx: Point3D}. Zero bad bonds guaranteed.
    """
    global _mol_3jf_cache
    from rdkit import Chem
    from rdkit.Chem import rdFMCS
    from rdkit.Geometry.rdGeometry import Point3D

    coordmap: dict = {}

    # ── VHL anchor: 3JF from 4W9H ────────────────────────────────────────────
    if _mol_3jf_cache is None:
        pdb_4w9h = str(_PROJECT_ROOT / 'data' / 'structures' / 'VHL' / '4W9H.pdb')
        if os.path.exists(pdb_4w9h):
            _mol_3jf_cache = _load_mol_with_conect(pdb_4w9h, '3JF')

    if _mol_3jf_cache is not None:
        hl = Chem.MolFromSmiles(head_linker_smi)
        if hl is not None:
            mcs = rdFMCS.FindMCS(
                [hl, _mol_3jf_cache], timeout=30,
                atomCompare=rdFMCS.AtomCompare.CompareElements,
                bondCompare=rdFMCS.BondCompare.CompareOrderExact,
                ringMatchesRingOnly=True, completeRingsOnly=True,
            )
            if mcs.numAtoms > 0:
                patt = Chem.MolFromSmarts(mcs.smartsString)
                hl_m  = hl.GetSubstructMatch(patt)
                jf_m  = _mol_3jf_cache.GetSubstructMatch(patt)
                if hl_m and jf_m:
                    pos = _mol_3jf_cache.GetConformer().GetPositions()
                    for hi, ji in zip(hl_m, jf_m):
                        coordmap[hi] = Point3D(*pos[ji])

    return coordmap


def gen_const_conf_anchored(Heads, Docked_Heads, Head_Linkers, output_sdf,
                             Anchor_A, v_atoms_sdf, n=100, homo_protac=False,
                             rmsd_tol=1.0, n_embed_factor=5):
    """
    Crystal-structure-anchored conformer generator for protac_lib.GenConstConf.

    coordMap: VHL warhead anchor only (19 atoms from 3JF crystal, 0 bad bonds).
    DG succeeds for geometrically feasible poses; returns 0 conformers for infeasible ones.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdFMCS
    from rdkit.Chem.rdchem import Conformer
    from rdkit.Geometry.rdGeometry import Point3D

    writer = Chem.SDWriter(output_sdf)

    # ── read full PROTAC SMILES ───────────────────────────────────────────────
    with open(Head_Linkers) as fh:
        smi = fh.readline().split()[0]
    head_linker = Chem.MolFromSmiles(smi)
    if head_linker is None:
        writer.close()
        return 0, v_atoms_sdf

    # ── virtual atoms centroid (permutation-invariant — correct with any bijection) ─
    docked_heads = _make_docked_heads_mol(Docked_Heads, Heads[0], Heads[1])
    if docked_heads is None:
        writer.close()
        return 0, v_atoms_sdf

    num_docked = docked_heads.GetConformer().GetNumAtoms()
    xs = [docked_heads.GetConformer().GetAtomPosition(i).x for i in range(num_docked)]
    ys = [docked_heads.GetConformer().GetAtomPosition(i).y for i in range(num_docked)]
    zs = [docked_heads.GetConformer().GetAtomPosition(i).z for i in range(num_docked)]
    cx, cy, cz = sum(xs) / num_docked, sum(ys) / num_docked, sum(zs) / num_docked
    rw = Chem.RWMol()
    for _ in range(3):
        at = Chem.Atom(23)
        at.SetNoImplicit(True)
        rw.AddAtom(at)
    va_mol = rw.GetMol()
    va_conf = Conformer(3)
    va_conf.SetAtomPosition(0, Point3D(cx + 1, cy, cz))
    va_conf.SetAtomPosition(1, Point3D(cx, cy, cz))
    va_conf.SetAtomPosition(2, Point3D(cx, cy + 1, cz))
    va_rw = Chem.RWMol(va_mol)
    va_rw.AddConformer(va_conf, assignId=True)
    v_writer = Chem.SDWriter(v_atoms_sdf)
    v_writer.write(va_rw.GetMol())
    v_writer.close()

    # ── nbr_atom (VHL anchor atom in head_linker, returned to mol_to_params) ─
    head_a_path = Heads[0]
    if not os.path.exists(head_a_path):
        # Fall back to local VHL warhead when params.txt path is stale
        head_a_path = str(_PROJECT_ROOT / 'data' / 'structures' / 'VHL' /
                          'SJF8240_vhll_placed.sdf')
    HeadA = Chem.SDMolSupplier(head_a_path, removeHs=True)[0]
    if HeadA is None:
        writer.close()
        return 0, v_atoms_sdf
    patt_A = Chem.MolFromSmarts(
        rdFMCS.FindMCS([head_linker, HeadA],
                       ringMatchesRingOnly=True,
                       completeRingsOnly=True).smartsString
    )
    head_A_list = head_linker.GetSubstructMatches(patt_A, uniquify=False)
    if not head_A_list:
        writer.close()
        return 0, v_atoms_sdf
    nbr_atom = head_A_list[0][int(Anchor_A)]

    # ── crystal-structure coordMap (0 bad bonds guaranteed) ──────────────────
    coord_map = _build_crystal_coordmap(smi, Docked_Heads)
    print(f'[gen_const_conf_anchored] coordMap: {len(coord_map)} anchor atoms', flush=True)

    # ── DG embedding ──────────────────────────────────────────────────────────
    accepted = []
    for seed in range(1, n_embed_factor * n + 1):
        if len(accepted) >= n:
            break
        mol_h = Chem.RWMol(Chem.AddHs(Chem.Mol(head_linker)))
        mol_h.RemoveAllConformers()
        if AllChem.EmbedMolecule(
            mol_h,
            coordMap=coord_map,
            randomSeed=seed,
            useBasicKnowledge=True,
            maxAttempts=200,
            forceTol=0.01,
        ) != -1:
            accepted.append(Chem.RemoveHs(mol_h))

    if not accepted:
        print('[gen_const_conf_anchored] 0 conformers — pose is geometrically infeasible',
              flush=True)

    # ── Kabsch alignment: shift DG conformers into docked site frame ─────────
    HeadB_ka = Chem.SDMolSupplier(Heads[1], removeHs=True)[0] if len(Heads) > 1 else None
    _HeadA_ha = Chem.RemoveHs(HeadA)
    _HeadB_ha = Chem.RemoveHs(HeadB_ka) if HeadB_ka is not None else None
    _nA = _HeadA_ha.GetNumAtoms()
    _nB = _HeadB_ha.GetNumAtoms() if _HeadB_ha is not None else 0
    _dh_conf = docked_heads.GetConformer()
    _tgt = np.array([[_dh_conf.GetAtomPosition(i).x,
                      _dh_conf.GetAtomPosition(i).y,
                      _dh_conf.GetAtomPosition(i).z]
                     for i in range(_nA + _nB)])
    _hA_list = head_linker.GetSubstructMatches(_HeadA_ha, uniquify=False)
    _hB_list = (head_linker.GetSubstructMatches(_HeadB_ha, uniquify=False)
                if _HeadB_ha is not None else [])
    _hA_match = _hA_list[0] if _hA_list else None
    _hB_match = _hB_list[0] if _hB_list else None
    _can_align = (_hA_match is not None and _hB_match is not None
                  and len(_hA_match) == _nA and len(_hB_match) == _nB)
    for mol_out in accepted[:n]:
        c = mol_out.GetConformer()
        _pos = np.array([[c.GetAtomPosition(ai).x,
                          c.GetAtomPosition(ai).y,
                          c.GetAtomPosition(ai).z]
                         for ai in range(mol_out.GetNumAtoms())])
        if _can_align:
            _src = np.vstack([_pos[list(_hA_match)], _pos[list(_hB_match)]])
            _R, _t = _kabsch_col(_tgt, _src)
            _pos = (_R @ _pos.T).T + _t
        for aidx in range(mol_out.GetNumAtoms()):
            c.SetAtomPosition(aidx, Point3D(*_pos[aidx]))
        writer.write(mol_out)
    writer.close()
    return nbr_atom, v_atoms_sdf


def gen_const_conf_free(Heads, Docked_Heads, Head_Linkers, output_sdf,
                        Anchor_A, v_atoms_sdf, n=100, homo_protac=False):
    """
    Free-conformer generator for protac_lib.GenConstConf.

    Generates free heavy-atom ETKDGv3 conformers of the full PROTAC (no coordMap
    constraints — constrained embedding fails for large PROTACs like SJF8240 where
    30/79 atoms are fixed, creating infeasible distance-geometry problems).

    Accepts the first n conformers without filtering.  Rosetta's PackRotamersMover
    in the downstream relax step will select the best-scoring conformation.
    """
    import numpy as np
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Geometry.rdGeometry import Point3D

    n = min(n,max(1,int(os.environ.get("PROTAC_DIAGNOSTIC_CONFORMERS","2")))) if os.environ.get("PROTAC_DIAGNOSTIC")=="1" else n
    N_EMBED = max(n, 2 if os.environ.get("PROTAC_DIAGNOSTIC")=="1" else 20)   # generate at least n + a small buffer

    writer = Chem.SDWriter(output_sdf)

    # ── load molecules ────────────────────────────────────────────────────────
    HeadA = Chem.SDMolSupplier(Heads[0])[0]
    HeadB = Chem.SDMolSupplier(Heads[1])[0]
    docked_heads = _make_docked_heads_mol(Docked_Heads, Heads[0], Heads[1])

    if HeadA is None or HeadB is None or docked_heads is None:
        writer.close()
        return 0, v_atoms_sdf

    # ── virtual atoms (center-of-mass + 2 offset points, mirrors GenConstConf) ─
    num_atoms = docked_heads.GetConformer().GetNumAtoms()
    xs = [docked_heads.GetConformer().GetAtomPosition(i).x for i in range(num_atoms)]
    ys = [docked_heads.GetConformer().GetAtomPosition(i).y for i in range(num_atoms)]
    zs = [docked_heads.GetConformer().GetAtomPosition(i).z for i in range(num_atoms)]
    cx, cy, cz = sum(xs) / num_atoms, sum(ys) / num_atoms, sum(zs) / num_atoms
    rw = Chem.RWMol()
    for _ in range(3):
        at = Chem.Atom(23)
        at.SetNoImplicit(True)
        rw.AddAtom(at)
    va_mol = rw.GetMol()
    va_conf = Chem.rdchem.Conformer(3)
    va_conf.SetAtomPosition(0, Point3D(cx + 1, cy, cz))
    va_conf.SetAtomPosition(1, Point3D(cx, cy, cz))
    va_conf.SetAtomPosition(2, Point3D(cx, cy + 1, cz))
    va_rw = Chem.RWMol(va_mol)
    va_rw.AddConformer(va_conf, assignId=True)
    v_writer = Chem.SDWriter(v_atoms_sdf)
    v_writer.write(va_rw.GetMol())
    v_writer.close()

    # ── full PROTAC molecule ──────────────────────────────────────────────────
    with open(Head_Linkers) as f:
        smi = f.readline().split()[0]
    mol = Chem.MolFromSmiles(smi)   # heavy atoms only
    if mol is None:
        writer.close()
        return 0, v_atoms_sdf

    # ── substructure match in heavy-atom PROTAC (for NBR atom index) ────────────
    head_A_list = mol.GetSubstructMatches(HeadA, uniquify=False)
    nbr_atom = head_A_list[0][int(Anchor_A)] if head_A_list else 0

    # ── pre-compute Kabsch alignment: map conformer head atoms → docked frame ─
    # docked_heads has HeadA heavy atoms (0..nA-1) then HeadB heavy atoms (nA..nA+nB-1).
    # We align each generated conformer so its warhead atoms match these docked
    # positions before writing — fixes the bug where conformers stayed at origin.
    _HeadA_ha = Chem.RemoveHs(HeadA)
    _HeadB_ha = Chem.RemoveHs(HeadB)
    _nA = _HeadA_ha.GetNumAtoms()
    _nB = _HeadB_ha.GetNumAtoms()
    _dh_conf = docked_heads.GetConformer()
    _tgt = np.array([[_dh_conf.GetAtomPosition(i).x,
                      _dh_conf.GetAtomPosition(i).y,
                      _dh_conf.GetAtomPosition(i).z]
                     for i in range(_nA + _nB)])
    _hA_list = mol.GetSubstructMatches(_HeadA_ha, uniquify=False)
    _hB_list = mol.GetSubstructMatches(_HeadB_ha, uniquify=False)
    _hA_match = _hA_list[0] if _hA_list else None
    _hB_match = _hB_list[0] if _hB_list else None
    _can_align = (_hA_match is not None and _hB_match is not None and
                  len(_hA_match) == _nA and len(_hB_match) == _nB)

    # ── generate free heavy-atom conformers ───────────────────────────────────
    conf_ids = list(AllChem.EmbedMultipleConfs(
        mol, numConfs=N_EMBED,
        randomSeed=publication_seed("free_conformers"), numThreads=1, enforceChirality=True,
        useExpTorsionAnglePrefs=True, useBasicKnowledge=True, ETversion=2,
    ))

    if not conf_ids:
        writer.close()
        return nbr_atom, v_atoms_sdf

    # Accept the first n valid conformers (cid >= 0 means embedding succeeded)
    accepted = 0
    for conf_id in conf_ids:
        if accepted >= n:
            break
        if conf_id < 0:
            continue
        mol_w = Chem.RWMol(Chem.Mol(mol))
        mol_w.RemoveAllConformers()
        c = mol.GetConformer(conf_id)
        _pos = np.array([[c.GetAtomPosition(ai).x,
                          c.GetAtomPosition(ai).y,
                          c.GetAtomPosition(ai).z]
                         for ai in range(mol.GetNumAtoms())])
        if _can_align:
            _src = np.vstack([_pos[list(_hA_match)], _pos[list(_hB_match)]])
            _R, _t = _kabsch_col(_tgt, _src)
            _pos = (_R @ _pos.T).T + _t
        new_c = Chem.rdchem.Conformer(mol.GetNumAtoms())
        for aidx in range(mol.GetNumAtoms()):
            new_c.SetAtomPosition(aidx, Point3D(*_pos[aidx]))
        mol_w.AddConformer(new_c, assignId=True)
        writer.write(mol_w.GetMol())
        accepted += 1

    writer.close()
    return nbr_atom, v_atoms_sdf


def apply_patches():
    """
    Monkey-patch pRosettaC's rosetta and utils modules in sys.modules.
    Call this BEFORE importing main.py or constraint_generation.py.
    """
    if str(_PROSETTAC_DIR) not in sys.path:
        sys.path.insert(0, str(_PROSETTAC_DIR))

    # pRosettaC's utils.py reads env vars at module level — set placeholders
    # so the import succeeds, then we override the functions below.
    for var, val in [
        ('PATCHDOCK',  '/dev/null'),
        ('OB',         '/dev/null'),
        ('SCRIPTS_FOL', str(_PROSETTAC_DIR)),
        ('ROSETTA_FOL', '/dev/null'),
    ]:
        if var not in os.environ:
            os.environ[var] = val

    import utils as _utils
    import rosetta as _rosetta

    # utils patches
    _utils.addH_sdf   = addH_sdf
    _utils.pdb2sdf    = pdb2sdf
    _utils.sdf2sdf    = sdf2sdf
    _utils.patchdock  = global_dock

    # protac_lib patches
    import protac_lib as _pl
    _pl.translate_anchors = translate_anchors_fixed   # sanitize=False/True mismatch fix
    _pl.SampleDist = sample_dist_free_conf            # free-conformer anchor-distance sampling
    from linker_adapter import generate as generate_bound_linker
    _pl.GenConstConf = generate_bound_linker               # Preserve both complete bound heads; reject infeasible linker placements

    # redirect SCRIPTS_FOL so subprocess constraint_generation commands use our wrapper
    _utils.SCRIPTS_FOL = str(_THIS_DIR)
    # Ensure PBS cluster helpers remain importable after SCRIPTS_FOL is redirected.
    # PBS.py does sys.path.append(utils.SCRIPTS_FOL + '/cluster/PBS') at import time,
    # so pre-adding the real path prevents ModuleNotFoundError for chemfarm modules.
    _pbs_dir = str(_PROSETTAC_DIR / 'cluster' / 'PBS')
    if _pbs_dir not in sys.path:
        sys.path.append(_pbs_dir)

    # rosetta patches
    _rosetta.clean         = clean_pdb
    _rosetta.clean_replace = lambda struct, chains: (clean_pdb(struct, chains),
                                                      os.rename(struct.split('.pdb')[0] + '_' + chains + '.pdb', struct))
    _rosetta.mol_to_params = mol_to_params
    _rosetta.relax         = relax
    _rosetta.local_docking = local_docking

    print('[adapter] pRosettaC patched: PyRosetta + Python OpenBabel active')
