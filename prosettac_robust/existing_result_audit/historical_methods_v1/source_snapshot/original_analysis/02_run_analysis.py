#!/usr/bin/env python3
"""
run_analysis.py  —  PROTAC ternary complex MD analysis
======================================================
Analyses 100 ns GROMACS trajectories for pRosettaC MD candidates.
Runs GROMACS tools for all measurements (handles topology natively),
then plots results and exports a video.

Metrics per POI:
  1. PROTAC → VHL   minimum atomic distance (Å) — bridging criterion
  2. PROTAC → kinase minimum atomic distance (Å) — bridging criterion
  3. VHL + kinase Cα RMSD vs t=0 (Å, after alignment) — ternary stability
  4. PROTAC RMSD vs t=0 (Å, after alignment to VHL+kinase Cα) — PROTAC stability

Outputs per POI in ~/protac-md/analysis/<POI>/rank1/:
  index.ndx              — custom GROMACS index with VHL / Kinase / PROTAC groups
  prod_nojump.xtc        — PBC-corrected trajectory (centered, no jumps)
  stripped.xtc           — protein + PROTAC only (for PyMOL visualization)
  dist_protac_vhl.xvg    — PROTAC-VHL min distance over time
  dist_protac_kinase.xvg — PROTAC-kinase min distance over time
  rmsd_complex.xvg       — VHL+kinase Cα RMSD vs t=0
  rmsd_protac.xvg        — PROTAC RMSD vs t=0 (after backbone fit)
  metrics.csv            — all 4 metrics per frame, with summary stats header
  dashboard.png          — 2×2 panel plot of all 4 metrics
  dashboard_video.mp4    — animated dashboard video (100 ns)

Global summary in ~/protac-md/analysis/:
  summary_all_pois.png   — all POIs overlaid on same axes

Usage:
  python run_analysis.py              # all POIs with prod.xtc present
  python run_analysis.py DDR2 MET    # specific POIs only
"""

import sys, os, subprocess, textwrap
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation

ROOT     = Path.home() / 'protac-md'
SIM_ROOT = ROOT / 'simulations'
ANA_ROOT = ROOT / 'analysis'

ALL_POIS = ['DDR2', 'MET', 'RIPK2_full', 'MAPK14_full', 'EPHB2',
            'AXL', 'ABL1', 'EPHA2', 'TNIK', 'MAP4K5', 'SLK']
POI_LABELS = {
    'DDR2': 'DDR2', 'MET': 'MET',
    'RIPK2_full': 'RIPK2', 'MAPK14_full': 'MAPK14', 'EPHB2': 'EPHB2',
    'AXL': 'AXL', 'ABL1': 'ABL1', 'EPHA2': 'EPHA2', 'TNIK': 'TNIK',
    'MAP4K5': 'MAP4K5', 'SLK': 'SLK',
}

RANK = 1
REP  = 1

PASS_THRESHOLD     = 5.0   # Å  pRosettaC PASS bridging criterion
MARGINAL_THRESHOLD = 10.0  # Å  pRosettaC MARGINAL bridging criterion


# ─── GROMACS command runner ───────────────────────────────────────────────────

def gmx(args, stdin='', cwd=None):
    """
    Run gmx subcommand. args = string of flags (e.g. 'mindist -s ...')
    stdin = text piped to stdin (newlines supported).
    Returns CompletedProcess. Does not raise on failure.
    """
    GMX_BIN = Path.home() / '.conda/envs/gromacs/bin/gmx'
    gmx_exe = str(GMX_BIN) if GMX_BIN.exists() else 'gmx'
    argv = [gmx_exe] + args.split()
    r = subprocess.run(argv, input=stdin, text=True, capture_output=True, cwd=cwd)
    return r


# ─── XVG parser ──────────────────────────────────────────────────────────────

def read_xvg(path):
    """
    Parse GROMACS .xvg file.
    Returns (times, values) as 1-D and 2-D numpy arrays.
    time is in whatever units gmx wrote (ns if -tu ns was passed).
    values are in gmx native units (nm for distances/RMSD).
    """
    times, vals = [], []
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if not s or s[0] in ('#', '@'):
                continue
            p = s.split()
            times.append(float(p[0]))
            vals.append([float(x) for x in p[1:]])
    t = np.asarray(times)
    v = np.asarray(vals) if vals else np.zeros((0, 1))
    return t, v


# ─── Chain boundary detection ─────────────────────────────────────────────────

def find_chain_atom_ranges(protein_gro):
    """
    Detect chain boundaries in a pdb2gmx-produced GRO by looking for large
    jumps in residue numbers between consecutive atoms.

    Returns (chains, atoms) where chains is a list of dicts:
      {'start': int, 'end': int, 'first_res': int, 'last_res': int}
    and atoms is the raw GRO atom line list (for fallback kinase detection).

    Chain order in pdb2gmx output (for our split PDBs):
      Chain A (VHL, first_res=1),  Chain D (kinase, first_res=200+/553+/1053+),
      Chain B (ElonginB),           Chain C (ElonginC)
    """
    lines = open(protein_gro).readlines()
    n_atoms = int(lines[1])
    atoms = lines[2 : 2 + n_atoms]

    boundaries = [0]  # 0-based indices where new chain starts
    prev_res = None
    for i, l in enumerate(atoms):
        resnum = int(l[0:5])
        if prev_res is not None and abs(resnum - prev_res) > 30:
            boundaries.append(i)
        prev_res = resnum
    boundaries.append(n_atoms)

    chains = []
    for j in range(len(boundaries) - 1):
        s = boundaries[j]      # 0-based start
        e = boundaries[j+1]    # 0-based exclusive end
        chains.append({
            'start':     s + 1,              # 1-based inclusive
            'end':       e,                  # 1-based inclusive
            'first_res': int(atoms[s][0:5]),
            'last_res':  int(atoms[e-1][0:5]),
            'n_atoms':   e - s,
        })
    return chains, atoms


def identify_vhl_kinase(chains, atoms):
    """
    From the chain list, pick VHL (first_res == 1) and kinase (first_res > 103).
    Falls back to scanning atom-by-atom for residue 104 when the kinase starts
    immediately after VHL (no large residue jump, e.g. TNIK where poi_min=104).
    Returns (vhl_chain, kinase_chain) dicts.
    """
    vhl = next(c for c in chains if c['first_res'] == 1)
    try:
        kinase = next(c for c in chains if c['first_res'] > 103)
    except StopIteration:
        # kinase starts at residue 104 — no jump detected; scan for first res=104 atom
        kinase_start = None
        for i, l in enumerate(atoms):
            if int(l[0:5]) == 104:
                kinase_start = i + 1  # 1-based
                break
        if kinase_start is None:
            raise RuntimeError('Could not locate kinase (res>103) in GRO file')
        # find the segment that contains kinase_start and extend to its end
        seg = next(c for c in chains if c['start'] <= kinase_start <= c['end'])
        kinase = {
            'start':     kinase_start,
            'end':       seg['end'],
            'first_res': 104,
            'last_res':  seg['last_res'],
            'n_atoms':   seg['end'] - kinase_start + 1,
        }
        # Truncate VHL to exclude kinase atoms (vhl['end'] spanned the whole
        # combined block because there was no residue-number jump at res 103→104)
        vhl = {**vhl,
               'end':      kinase_start - 1,
               'last_res': 103,
               'n_atoms':  kinase_start - vhl['start']}
    return vhl, kinase


# ─── Custom GROMACS index ─────────────────────────────────────────────────────

def build_index(gro, tpr, protein_gro, ana_dir):
    """
    Build ~/protac-md/analysis/<POI>/rank1/index.ndx with groups:
      VHL, Kinase, VHL_CA, Kinase_CA, VHL_Kinase_CA, PROTAC, Protein_Other
    """
    ndx = ana_dir / 'index.ndx'
    if ndx.exists():
        return ndx

    chains, atoms = find_chain_atom_ranges(protein_gro)
    vhl, kinase = identify_vhl_kinase(chains, atoms)

    print(f'    VHL   : atoms {vhl["start"]:6d} – {vhl["end"]:6d} '
          f'(res {vhl["first_res"]}–{vhl["last_res"]}, {vhl["n_atoms"]} atoms)')
    print(f'    Kinase: atoms {kinase["start"]:6d} – {kinase["end"]:6d} '
          f'(res {kinase["first_res"]}–{kinase["last_res"]}, {kinase["n_atoms"]} atoms)')

    # Build make_ndx command string (interactive mode)
    # Groups created:
    #   default groups from .gro (0-20 typically)
    #   then we add VHL, Kinase, VHL_CA, Kinase_CA, VHL_Kinase_CA, Protein_Other
    # GROMACS make_ndx syntax:
    #   a 1-N   : select atoms 1 to N
    #   name X Y: rename group X to Y
    #   G1 | G2 : union
    #   G1 & t CA: intersect with CA atoms (atom type)
    #   t CA     : select all atoms named CA (atom name)

    n_default = 21   # default groups are 0..20 (21 groups total)
    cmds = textwrap.dedent(f"""\
        a {vhl['start']}-{vhl['end']}
        name {n_default} VHL
        a {kinase['start']}-{kinase['end']}
        name {n_default+1} Kinase
        a {vhl['start']}-{vhl['end']}
        name {n_default+2} VHL_all2
        {n_default+2} & 3
        name {n_default+3} VHL_CA
        a {kinase['start']}-{kinase['end']}
        name {n_default+4} Kinase_all2
        {n_default+4} & 3
        name {n_default+5} Kinase_CA
        {n_default+3} | {n_default+5}
        name {n_default+6} VHL_Kinase_CA
        1 | 12
        name {n_default+7} Protein_Other
        q
    """)

    r = gmx(f'make_ndx -f {gro} -o {ndx}', stdin=cmds)
    if not ndx.exists():
        raise RuntimeError(f'make_ndx failed:\n{r.stderr[-400:]}')
    return ndx


# ─── Trajectory preprocessing ─────────────────────────────────────────────────

def preprocess_traj(tpr, xtc_in, ndx, ana_dir):
    """
    Two steps:
    1. gmx trjconv -pbc nojump -center (center=Protein, output=System)
    2. gmx trjconv to strip water, output protein+PROTAC only xtc for PyMOL
    Returns (nojump_xtc, stripped_xtc).
    """
    xtc_nojump  = ana_dir / 'prod_nojump.xtc'
    xtc_stripped = ana_dir / 'stripped.xtc'

    if not xtc_nojump.exists():
        print('    Removing PBC jumps (center=Protein, output=System) ...')
        r = gmx(
            f'trjconv -s {tpr} -f {xtc_in} -o {xtc_nojump} '
            f'-pbc nojump -center -ur compact',
            stdin='Protein\nSystem\n'
        )
        if not xtc_nojump.exists():
            raise RuntimeError(f'trjconv nojump failed:\n{r.stderr[-300:]}')
    else:
        print('    prod_nojump.xtc exists — skipping')

    if not xtc_stripped.exists():
        print('    Stripping water/ions for PyMOL visualization ...')
        # "Protein_Other" = protein + PROTAC (group created in index)
        r = gmx(
            f'trjconv -s {tpr} -f {xtc_nojump} -o {xtc_stripped} '
            f'-n {ndx}',
            stdin='Protein_Other\n'
        )
    else:
        print('    stripped.xtc exists — skipping')

    return xtc_nojump, xtc_stripped


# ─── GROMACS analysis steps ───────────────────────────────────────────────────

def run_mindist(tpr, xtc, ndx, g1, g2, out):
    """Minimum atomic distance g1↔g2 over time. Output nm."""
    if out.exists():
        return
    r = gmx(f'mindist -s {tpr} -f {xtc} -n {ndx} -od {out} -tu ns',
            stdin=f'{g1}\n{g2}\n')
    if not out.exists():
        raise RuntimeError(f'gmx mindist failed (groups: {g1}, {g2}):\n{r.stderr[-400:]}')


def run_rms(tpr, xtc, ndx, fit_grp, rmsd_grp, out):
    """
    RMSD: fit trajectory to fit_grp, measure RMSD of rmsd_grp vs t=0.
    Reference = first frame of trajectory (from tpr). Output nm.
    """
    if out.exists():
        return
    r = gmx(f'rms -s {tpr} -f {xtc} -n {ndx} -o {out} -tu ns',
            stdin=f'{fit_grp}\n{rmsd_grp}\n')
    if not out.exists():
        raise RuntimeError(f'gmx rms failed (fit: {fit_grp}, rmsd: {rmsd_grp}):\n{r.stderr[-400:]}')


# ─── Main per-POI analysis ────────────────────────────────────────────────────

def analyse_poi(poi, rank=1, rep=1):
    """
    Full analysis pipeline for one POI. Returns data dict or None if skipped.
    """
    simdir  = SIM_ROOT / poi
    proddir = simdir / f'prod/rank{rank}/rep{rep}'
    tpr     = proddir / 'prod.tpr'
    xtc_in  = proddir / 'prod.xtc'
    gro     = simdir / f'topol/rank{rank}/complex_ions.gro'
    protein_gro = simdir / f'topol/rank{rank}/complex_protein.gro'
    ana_dir = ANA_ROOT / poi / f'rank{rank}'

    if not xtc_in.exists():
        print(f'[{poi}] prod.xtc not found — skipping')
        return None

    ana_dir.mkdir(parents=True, exist_ok=True)
    print(f'\n{"="*55}')
    print(f'  Analysing {poi} rank{rank}')
    print(f'{"="*55}')

    print('  [1/5] Building index ...')
    ndx = build_index(gro, tpr, protein_gro, ana_dir)

    print('  [2/5] Preprocessing trajectory ...')
    xtc, xtc_stripped = preprocess_traj(tpr, xtc_in, ndx, ana_dir)

    print('  [3/5] Computing PROTAC–VHL minimum distance ...')
    xvg_dv = ana_dir / 'dist_protac_vhl.xvg'
    run_mindist(tpr, xtc, ndx, 'Other', 'VHL', xvg_dv)

    print('  [3/5] Computing PROTAC–kinase minimum distance ...')
    xvg_dk = ana_dir / 'dist_protac_kinase.xvg'
    run_mindist(tpr, xtc, ndx, 'Other', 'Kinase', xvg_dk)

    print('  [4/5] Computing VHL+kinase Cα RMSD ...')
    xvg_rc = ana_dir / 'rmsd_complex.xvg'
    run_rms(tpr, xtc, ndx, 'VHL_Kinase_CA', 'VHL_Kinase_CA', xvg_rc)

    print('  [4/5] Computing PROTAC RMSD (after VHL+kinase fit) ...')
    xvg_rp = ana_dir / 'rmsd_protac.xvg'
    run_rms(tpr, xtc, ndx, 'VHL_Kinase_CA', 'Other', xvg_rp)

    # Parse XVG — gmx outputs nm; convert to Å (* 10)
    t_d, dv  = read_xvg(xvg_dv)
    _,   dk  = read_xvg(xvg_dk)
    t_r, rc  = read_xvg(xvg_rc)
    _,   rp  = read_xvg(xvg_rp)

    dv_A = dv[:, 0] * 10
    dk_A = dk[:, 0] * 10
    rc_A = rc[:, 0] * 10
    rp_A = rp[:, 0] * 10

    data = {
        'poi':     poi,
        'label':   POI_LABELS[poi],
        't_d':     t_d,
        'd_vhl':   dv_A,
        'd_kinase': dk_A,
        't_r':     t_r,
        'rmsd_c':  rc_A,
        'rmsd_p':  rp_A,
        'ana_dir': ana_dir,
        'xtc_stripped': xtc_stripped,
        'tpr':     tpr,
    }

    print('  [5/5] Saving CSV and plots ...')
    save_csv(data)
    plot_dashboard(data)
    make_dashboard_video(data)

    return data


# ─── Statistics ───────────────────────────────────────────────────────────────

def summarise(arr):
    return f'{arr.mean():.2f} ± {arr.std():.2f} Å  [min {arr.min():.2f}, max {arr.max():.2f}]'


def save_csv(data):
    n = min(len(data['t_d']), len(data['t_r']))
    out = data['ana_dir'] / 'metrics.csv'

    stats = (
        f'# {data["poi"]} MD metrics — 100 ns rank-1\n'
        f'# PROTAC-VHL:    {summarise(data["d_vhl"])}\n'
        f'# PROTAC-kinase: {summarise(data["d_kinase"])}\n'
        f'# Complex RMSD:  {summarise(data["rmsd_c"])}\n'
        f'# PROTAC RMSD:   {summarise(data["rmsd_p"])}\n'
        f'# Bridging (PROTAC-VHL<5Å AND PROTAC-kinase<5Å): '
        f'{((data["d_vhl"][:n]<5) & (data["d_kinase"][:n]<5)).mean()*100:.1f}% of frames\n'
        f'time_ns,d_protac_vhl_A,d_protac_kinase_A,rmsd_complex_A,rmsd_protac_A\n'
    )

    rows = np.column_stack([
        data['t_d'][:n],
        data['d_vhl'][:n],
        data['d_kinase'][:n],
        data['rmsd_c'][:n],
        data['rmsd_p'][:n],
    ])
    with open(out, 'w') as f:
        f.write(stats)
        np.savetxt(f, rows, delimiter=',', fmt='%.4f')

    print(f'    CSV: {out}')
    print(f'    PROTAC-VHL:    {summarise(data["d_vhl"])}')
    print(f'    PROTAC-kinase: {summarise(data["d_kinase"])}')
    print(f'    Complex RMSD:  {summarise(data["rmsd_c"])}')
    print(f'    PROTAC RMSD:   {summarise(data["rmsd_p"])}')
    n = min(len(data['d_vhl']), len(data['d_kinase']))
    pct = ((data['d_vhl'][:n] < 5) & (data['d_kinase'][:n] < 5)).mean() * 100
    print(f'    Bridging (both <5Å): {pct:.1f}% of trajectory')


# ─── Static dashboard plot ────────────────────────────────────────────────────

def _add_thresholds(ax):
    ax.axhline(PASS_THRESHOLD,     color='green',  ls='--', lw=1.2, alpha=0.7,
               label=f'PASS threshold ({PASS_THRESHOLD} Å)')
    ax.axhline(MARGINAL_THRESHOLD, color='orange', ls='--', lw=1.2, alpha=0.7,
               label=f'MARGINAL threshold ({MARGINAL_THRESHOLD} Å)')


def _rolling(t, y, window_ns=5.0):
    """Rolling mean. Returns same-length array."""
    if len(t) < 2:
        return y
    dt = (t[-1] - t[0]) / max(len(t) - 1, 1)
    w  = max(1, int(window_ns / dt))
    kernel = np.ones(w) / w
    return np.convolve(y, kernel, mode='same')


def plot_dashboard(data):
    """Static 2×2 dashboard saved as PNG."""
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
    fig.suptitle(
        f'{data["label"]}  —  100 ns MD (rank 1)\n'
        f'AMBER99SB-ILDN + GAFF2 + TIP3P | RTX 5090',
        fontsize=13, fontweight='bold'
    )

    t_d, t_r = data['t_d'], data['t_r']

    # ── PROTAC → VHL ──────────────────────────────────────────────────────────
    ax = axes[0, 0]
    ax.plot(t_d, data['d_vhl'], lw=0.5, color='steelblue', alpha=0.5)
    ax.plot(t_d, _rolling(t_d, data['d_vhl']), lw=2, color='navy',
            label='5 ns rolling mean')
    _add_thresholds(ax)
    ax.set(xlabel='Time (ns)', ylabel='Min distance (Å)',
           title='PROTAC → VHL (chain A)')
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=8)

    # ── PROTAC → kinase ───────────────────────────────────────────────────────
    ax = axes[0, 1]
    ax.plot(t_d, data['d_kinase'], lw=0.5, color='coral', alpha=0.5)
    ax.plot(t_d, _rolling(t_d, data['d_kinase']), lw=2, color='darkred',
            label='5 ns rolling mean')
    _add_thresholds(ax)
    ax.set(xlabel='Time (ns)', ylabel='Min distance (Å)',
           title='PROTAC → kinase (chain D)')
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=8)

    # ── VHL+kinase Cα RMSD ───────────────────────────────────────────────────
    ax = axes[1, 0]
    ax.plot(t_r, data['rmsd_c'], lw=1.2, color='seagreen')
    ax.set(xlabel='Time (ns)', ylabel='RMSD (Å)',
           title='VHL + kinase Cα RMSD (vs t=0)')
    ax.set_ylim(bottom=0)

    # ── PROTAC RMSD ──────────────────────────────────────────────────────────
    ax = axes[1, 1]
    ax.plot(t_r, data['rmsd_p'], lw=1.2, color='mediumpurple')
    ax.set(xlabel='Time (ns)', ylabel='RMSD (Å)',
           title='PROTAC RMSD (after VHL+kinase Cα fit)')
    ax.set_ylim(bottom=0)

    # Stats box at bottom
    n = min(len(data['d_vhl']), len(data['d_kinase']))
    pct_bridge = ((data['d_vhl'][:n] < 5) & (data['d_kinase'][:n] < 5)).mean() * 100
    txt = (
        f'Mean ± SD  |  PROTAC–VHL: {data["d_vhl"].mean():.1f}±{data["d_vhl"].std():.1f} Å  '
        f'|  PROTAC–kinase: {data["d_kinase"].mean():.1f}±{data["d_kinase"].std():.1f} Å  '
        f'|  Complex RMSD: {data["rmsd_c"].mean():.1f}±{data["rmsd_c"].std():.1f} Å  '
        f'|  PROTAC RMSD: {data["rmsd_p"].mean():.1f}±{data["rmsd_p"].std():.1f} Å  '
        f'|  Bridging (<5Å): {pct_bridge:.1f}%'
    )
    fig.text(0.5, -0.01, txt, ha='center', fontsize=9, style='italic',
             bbox=dict(boxstyle='round,pad=0.4', fc='lightyellow', alpha=0.9))

    out = data['ana_dir'] / 'dashboard.png'
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'    PNG: {out}')


# ─── Animated dashboard video ─────────────────────────────────────────────────

def make_dashboard_video(data, n_frames=300, fps=30):
    """
    Animated MP4 of the 4 metrics over time.
    Draws a cursor line that sweeps across all 4 panels simultaneously.
    n_frames: number of animation frames (downsampled from full trajectory)
    fps: output video frames per second
    """
    out = data['ana_dir'] / 'dashboard_video.mp4'
    if out.exists():
        print(f'    Video exists: {out}')
        return

    print(f'    Building animated video ({n_frames} frames @ {fps} fps) ...')

    t_d = data['t_d']
    t_r = data['t_r']

    # Choose time axis (use dist time for consistency)
    t_max = min(t_d[-1], t_r[-1])
    t_frames = np.linspace(0, t_max, n_frames)

    # Interpolate all series onto t_frames
    dv  = np.interp(t_frames, t_d, data['d_vhl'])
    dk  = np.interp(t_frames, t_d, data['d_kinase'])
    rc  = np.interp(t_frames, t_r, data['rmsd_c'])
    rp  = np.interp(t_frames, t_r, data['rmsd_p'])

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle(
        f'{data["label"]}  —  100 ns PROTAC Ternary Complex MD', fontsize=13,
        fontweight='bold'
    )

    configs = [
        (axes[0,0], dv,  'PROTAC → VHL (Å)',           'steelblue', True),
        (axes[0,1], dk,  'PROTAC → kinase (Å)',         'coral',     True),
        (axes[1,0], rc,  'VHL+kinase Cα RMSD (Å)',      'seagreen',  False),
        (axes[1,1], rp,  'PROTAC RMSD (Å)',             'mediumpurple', False),
    ]

    y_maxes = []
    lines, cursors, filled = [], [], []
    for ax, y, ylabel, col, thresholds in configs:
        ax.set_xlim(0, t_max)
        ymax = max(y.max() * 1.15, 12 if thresholds else 5)
        y_maxes.append(ymax)
        ax.set_ylim(0, ymax)
        ax.set_xlabel('Time (ns)', fontsize=9)
        ax.set_ylabel(ylabel, fontsize=9)
        ax.tick_params(labelsize=8)
        ax.set_facecolor('#f8f8f8')

        if thresholds:
            ax.axhline(PASS_THRESHOLD,     color='green',  ls='--', lw=1.0, alpha=0.6)
            ax.axhline(MARGINAL_THRESHOLD, color='orange', ls='--', lw=1.0, alpha=0.6)
            ax.text(t_max * 0.01, PASS_THRESHOLD + 0.3, 'PASS',
                    color='green', fontsize=7)
            ax.text(t_max * 0.01, MARGINAL_THRESHOLD + 0.3, 'MARGINAL',
                    color='orange', fontsize=7)

        # Full static trace (faint)
        ax.plot(t_frames, y, lw=0.8, color=col, alpha=0.25)

        # Active line (built frame by frame)
        line, = ax.plot([], [], lw=1.5, color=col)
        lines.append((line, y))

        # Cursor
        curs = ax.axvline(0, color='black', lw=1.2, ls='-', alpha=0.7)
        cursors.append(curs)

        # Fill under active portion
        fill = ax.fill_between([], [], color=col, alpha=0.15)
        filled.append((ax, fill, col))

    # Time label
    time_text = fig.text(0.5, 0.96, 'Time: 0.0 ns', ha='center',
                         fontsize=11, color='#333333')

    plt.tight_layout(rect=[0, 0, 1, 0.94])

    def init():
        for line, _ in lines:
            line.set_data([], [])
        for curs in cursors:
            curs.set_xdata([0])
        time_text.set_text('Time: 0.0 ns')
        return [l for l, _ in lines] + cursors + [time_text]

    def update(frame):
        t_now = t_frames[frame]
        mask  = t_frames[:frame+1]

        for (line, y_full), curs, (ax, old_fill, col) in zip(lines, cursors, filled):
            line.set_data(mask, y_full[:frame+1])
            curs.set_xdata([t_now])

            # Update fill (remove old, add new)
            old_fill.remove()
            new_fill = ax.fill_between(mask, y_full[:frame+1], color=col, alpha=0.15)
            filled[filled.index((ax, old_fill, col))] = (ax, new_fill, col)

        time_text.set_text(f'Time: {t_now:.1f} ns')
        return [l for l, _ in lines] + cursors + [time_text]

    ani = animation.FuncAnimation(
        fig, update, frames=n_frames, init_func=init,
        interval=1000/fps, blit=False
    )

    writer = animation.FFMpegWriter(fps=fps, bitrate=1800,
                                    extra_args=['-vcodec', 'libx264', '-pix_fmt', 'yuv420p'])
    ani.save(str(out), writer=writer, dpi=120)
    plt.close(fig)
    print(f'    Video: {out}')


# ─── Summary plot (all POIs) ─────────────────────────────────────────────────

def plot_summary(results):
    """All POIs overlaid on 4 shared panels."""
    if not results:
        return

    cmap = plt.cm.tab10
    colours = [cmap(i / 10) for i in range(len(results))]

    fig, axes = plt.subplots(2, 2, figsize=(14, 9), constrained_layout=True)
    fig.suptitle('PROTAC Ternary Complex MD — all POIs (rank 1)',
                 fontsize=14, fontweight='bold')

    for (data, col) in zip(results, colours):
        kw = dict(lw=1.8, color=col, label=data['label'])
        axes[0,0].plot(data['t_d'], data['d_vhl'],    **kw)
        axes[0,1].plot(data['t_d'], data['d_kinase'], **kw)
        axes[1,0].plot(data['t_r'], data['rmsd_c'],   **kw)
        axes[1,1].plot(data['t_r'], data['rmsd_p'],   **kw)

    for ax, title, ylabel in [
        (axes[0,0], 'PROTAC → VHL min distance',    'Min dist (Å)'),
        (axes[0,1], 'PROTAC → kinase min distance', 'Min dist (Å)'),
        (axes[1,0], 'VHL+kinase Cα RMSD',           'RMSD (Å)'),
        (axes[1,1], 'PROTAC RMSD',                  'RMSD (Å)'),
    ]:
        ax.set(title=title, xlabel='Time (ns)', ylabel=ylabel)
        ax.set_ylim(bottom=0)
        ax.legend(fontsize=9)

    for ax in axes[0, :]:
        ax.axhline(PASS_THRESHOLD,     color='green',  ls='--', lw=1, alpha=0.5)
        ax.axhline(MARGINAL_THRESHOLD, color='orange', ls='--', lw=1, alpha=0.5)

    out = ANA_ROOT / 'summary_all_pois.png'
    ANA_ROOT.mkdir(exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'\nSummary plot: {out}')


# ─── Entry point ─────────────────────────────────────────────────────────────

def main():
    # Usage: run_analysis.py [POI ...] [--rank N] [--rep N]
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('pois', nargs='*', default=ALL_POIS)
    parser.add_argument('--rank', type=int, default=1)
    parser.add_argument('--rep',  type=int, default=1)
    args = parser.parse_args()

    unknown = [p for p in args.pois if p not in ALL_POIS]
    if unknown:
        print(f'Unknown POIs: {unknown}. Valid: {ALL_POIS}')
        sys.exit(1)

    results = []
    for poi in args.pois:
        data = analyse_poi(poi, rank=args.rank, rep=args.rep)
        if data:
            results.append(data)

    if len(results) > 1:
        plot_summary(results)

    print('\n=== Analysis complete ===')
    for d in results:
        print(f'  {d["label"]:12s} → {d["ana_dir"]}')
    print()
    print('To view trajectory in PyMOL:')
    for d in results:
        print(f'  pymol {d["tpr"]}  {d["xtc_stripped"]}')


if __name__ == '__main__':
    main()
