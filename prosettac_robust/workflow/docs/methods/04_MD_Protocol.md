# Repaired MD preparation and validation protocol

This protocol follows the working protac-md AMBER/GAFF2 baseline. It describes
actual generated inputs and the local validation; it does not certify trajectory
convergence or prescribe a universally adequate Nature simulation length.

## Chemistry and system preparation

Protein: AMBER99SB-ILDN. Ligand: GAFF2 with separately verified AM1-BCC charges.
Water: TIP3P. Temperature: 300 K. NaCl: 0.15 M with neutralization. Solvent box:
dodecahedral, 1.5 nm padding. The selected ligand microstate, exact atom map,
protein completion model and their hashes are recorded in md_preparation.json.
Neutral supplier graphs are reference states, not predicted physiological populations.
Protein protonation/terminal choices follow the recorded pdb2gmx preparation and
require assessment for the intended assay conditions; a ligand-only pH enumeration
does not establish all protein protonation states.

Observed protein heavy atoms remain unchanged during completion; rebuilt internal
loops are explicitly modeled. Ligand heavy atoms are placed using the complete
stereochemical graph and atom map; hydrogens are reconstructed in that frame.
No missing atoms or Gromacs warnings are silently bypassed. Force-field files are
copied locally, and the topology, restraints, index and stage MDPs are hash-bound.

## Stage definitions

| Stage | Generated baseline settings |
|---|---|
| EM | Steepest descent, at most 50,000 steps, emtol 500 kJ/mol/nm, emstep 0.01 nm |
| NVT | 200 ps, protein and ligand heavy-atom restraints, V-rescale thermostat |
| Restrained NPT | 500 ps, same restraints, isotropic C-rescale barostat |
| Unrestrained NPT | 500 ps, restraints removed, same thermostat/barostat |
| Production | Input-builder default 20 ns for the inherited screening protocol; final study length/replicas must be justified and frozen separately |

Dynamics: 2 fs timestep, hydrogen-bond constraints, LINCS order 4/2 iterations.
Electrostatics: PME, Coulomb cutoff 1.0 nm, PME order 4, Fourier spacing 0.12 nm.
Van der Waals cutoff: 1.0 nm, with EnerPres dispersion correction. Thermostat:
whole System, tau-t 1 ps. Barostat: 1 bar, tau-p 5 ps, compressibility 4.5e-5/bar.
Exact velocity/stochastic seeds are in the generated MDPs; independent replicas
must use independently specified seeds and separate system directories.

The inherited production length is not declared sufficient for convergence.
Replicate-level uncertainty, sampling adequacy, chemical-state sensitivity and
selection/exclusion rules must support the actual intended claims. Correlated
frames are not counted as independent simulations or experimental replicates.

## Output and analysis scope

The current protocol writes compressed SOLUTE trajectories every 10 ps and energy
records every 2 ps. Full-system final coordinates and checkpoints are retained.
This trajectory choice supports mapped protein/ligand analysis but does not retain
all time-dependent solvent coordinates; water-mechanism claims would require an
appropriate solvent-retaining output protocol and additional storage planning.
Full raw retention means all produced artifacts are preserved; it does not imply
that unrecorded solvent coordinates are available.

06_reduce_md_topology.py makes a verified topology for an ordered SOLUTE trajectory.
It checks the executed parent TPR and exact index/preparation hashes, rejects
permuted or incomplete solute groups, and records conversion provenance.
02_analyze_md_contacts.py uses mapped heavy atoms and triclinic periodic distances.
Equilibration analysis is explicitly diagnostic, never production evidence.
05_md_thermodynamics.py extracts receipt-bound energy data with explicit discard
and block windows; descriptive block means do not establish convergence or
independent-replica confidence intervals.

## Local control evidence

The neutral GSK215 validation system minimized to maximum force 472.168 kJ/mol/nm;
Gromacs reports convergence in 4,222 steps. It then completed all 200 ps NVT,
500 ps restrained NPT and 500 ps unrestrained NPT stages without accepted LINCS
warnings or truncated stages. Exact TPR and parent-output checks passed.

For the post-hoc diagnostic window 250–500 ps of unrestrained NPT, the recorded
energy samples average 299.989 K, pressure 8.964 bar and density 1004.720 kg/m^3.
Instantaneous pressure fluctuates substantially; these descriptive values are not
independent-replica statistics. Gromacs accumulator statistics over the entire
stage differ from means of the saved samples and are retained separately.

The 26 retained trajectory frames in that window maintain mapped contacts within
5 Å of both proteins. This is a short native-complex control, not proof of binding
free energy, long-term stability, degradation or conformational convergence.
No production trajectory or independent MD replica campaign has been launched.
