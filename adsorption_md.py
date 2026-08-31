"""Unbiased adsorption MD: starting in bulk water, does the peptide bind the
bilayer, and how deep does it go?

This replaces the steered-MD protocol. Steered MD needs a correct work
integral and a long pull path, and the earlier runs had neither -- the
peptide was already buried in the hydrocarbon core, so the restraint moved
0.2 A and the "work" was positional jitter. Here there is no bias at all:
the peptide starts solvated above the membrane and either binds or does not.
The observables are geometric, so nothing depends on a force constant or a
pulling rate.

Reported per frame:
  z_rel      peptide COM z minus the near phosphate plane (negative = inserted)
  min_dist   closest approach of any peptide atom to any lipid heavy atom
  contacts   lipid heavy atoms within 4 A of the peptide, split POPC vs POPS
  n_wat      waters within 4 A of the peptide (desolvation on binding)

Usage:
  python adsorption_md.py <system> <replica> <gpu_index>
"""
import openmm as mm
import openmm.app as app
import openmm.unit as unit
import numpy as np
import json, sys, os, time, glob

system_name = sys.argv[1] if len(sys.argv) > 1 else 'PrP_23_93'
replica = int(sys.argv[2]) if len(sys.argv) > 2 else 0
gpu_index = sys.argv[3] if len(sys.argv) > 3 else '0'

PATHS = {
    'PrP_23_93': 'charmm_gui_prp2393/openmm',
    'KKRPKP_water': 'charmm_gui_kkrpkp_water/openmm',
    'NNRPNP_water': 'charmm_gui_nnrpnp_water/openmm',
}
rel = PATHS[system_name]
for base in ['/workspace/prion-neurotoxicity', '/Users/allan/Projects/cjd']:
    p = os.path.join(base, rel)
    if os.path.isdir(p):
        os.chdir(p)
        break
else:
    sys.exit(f"ERROR: cannot find system directory for {system_name} ({rel})")

OUT_DIR = '/workspace/adsorption' if os.path.isdir('/workspace') else \
          '/Users/allan/Projects/cjd/results_adsorption'
os.makedirs(OUT_DIR, exist_ok=True)
tag = f"{system_name}_rep{replica}"
SEED = 20000 + replica * 17 + abs(hash(system_name)) % 1000

PROD_NS = float(os.environ.get('PROD_NS', 50))
EQUIL_NS = float(os.environ.get('EQUIL_NS', 2))
DT = 0.002
SAMPLE_PS = 20.0   # frame every 20 ps

print(f"=== {tag} | GPU {gpu_index} | seed {SEED} ===", flush=True)

psf = app.CharmmPsfFile('step5_input.psf')
pdb = app.PDBFile('step5_input.pdb')
with open('sysinfo.dat') as f:
    dims = json.load(f)['dimensions']
psf.setBox(dims[0]*unit.angstroms, dims[1]*unit.angstroms, dims[2]*unit.angstroms)
print(f"Box: {dims[0]:.1f} x {dims[1]:.1f} x {dims[2]:.1f} A", flush=True)

toppar_dir = os.path.join(os.path.dirname(os.getcwd()), 'toppar')
pfiles = (sorted(glob.glob(os.path.join(toppar_dir, '*.rtf'))) +
          sorted(glob.glob(os.path.join(toppar_dir, '*.prm'))) +
          sorted(glob.glob(os.path.join(toppar_dir, '*.str'))) + ['toppar.str'])
params = app.CharmmParameterSet(*pfiles)

system = psf.createSystem(params, nonbondedMethod=app.PME,
                          nonbondedCutoff=1.2*unit.nanometers,
                          switchDistance=1.0*unit.nanometers,
                          constraints=app.HBonds)
system.addForce(mm.MonteCarloMembraneBarostat(
    1.0*unit.bar, 0.0*unit.bar*unit.nanometers, 303.15*unit.kelvin,
    mm.MonteCarloMembraneBarostat.XYIsotropic,
    mm.MonteCarloMembraneBarostat.ZFree))

# --- atom selections -------------------------------------------------------
PROTEIN_RES = {'ALA','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS','HSD','HSE','HSP',
               'ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'}
# CharmmPsfFile renames water on read: the PSF says TIP3/OH2, the OpenMM
# topology says HOH with atoms O/H1/H2. Accept both spellings.
WATER_RES = {'HOH', 'TIP3', 'WAT', 'SOL'}
WATER_O = {'O', 'OH2', 'OW'}
pep, popc, pops, phos, wat_o = [], [], [], [], []
for a in psf.topology.atoms():
    rn, nm = a.residue.name, a.name
    if rn in WATER_RES:
        if nm in WATER_O:
            wat_o.append(a.index)
        continue
    if nm.startswith('H'):
        continue
    if rn in PROTEIN_RES:
        pep.append(a.index)
    elif rn == 'POPC':
        popc.append(a.index)
        if nm == 'P':
            phos.append(a.index)
    elif rn == 'POPS':
        pops.append(a.index)
        if nm == 'P':
            phos.append(a.index)
pep = np.array(pep, dtype=int); popc = np.array(popc, dtype=int)
pops = np.array(pops, dtype=int); phos = np.array(phos, dtype=int)
wat_o = np.array(wat_o, dtype=int)
lipid = np.concatenate([popc, pops])
print(f"peptide {len(pep)} | POPC {len(popc)} | POPS {len(pops)} | "
      f"P {len(phos)} | water {len(wat_o)}", flush=True)
for name, sel in [('peptide', pep), ('POPC', popc), ('POPS', pops),
                  ('phosphate', phos), ('water', wat_o)]:
    if len(sel) == 0:
        sys.exit(f"ERROR: {name} selection is empty -- check residue naming")

platform = mm.Platform.getPlatformByName('CUDA')
props = {'Precision': 'mixed', 'DeviceIndex': gpu_index}
integ = mm.LangevinMiddleIntegrator(303.15*unit.kelvin, 1.0/unit.picosecond,
                                    DT*unit.picoseconds)
integ.setRandomNumberSeed(SEED)
sim = app.Simulation(psf.topology, system, integ, platform, props)
sim.context.setPositions(pdb.positions)

print("Minimizing...", flush=True)
sim.minimizeEnergy(maxIterations=10000)
sim.context.setVelocitiesToTemperature(303.15*unit.kelvin, SEED)
print(f"Equilibrating {EQUIL_NS} ns...", flush=True)
sim.step(int(EQUIL_NS * 1e6 / (DT * 1000)))


# is_pops marks which entries of `lipid` are the anionic lipid, so one
# distance matrix serves both the contact count and the POPC/POPS split.
is_pops = np.concatenate([np.zeros(len(popc), bool), np.ones(len(pops), bool)])


def _near(pos, sel, box, margin=1.2):
    """Indices of `sel` inside the peptide bounding box plus a margin.

    Computing every peptide-lipid distance is ~5M pairs per frame, which
    leaves the GPU idle waiting on numpy. The peptide is a compact object, so
    a bounding-box prefilter discards almost all candidates at a fraction of
    the cost and cannot change any distance below `margin`.
    """
    lo, hi = box
    if len(sel) == 0:
        return sel, np.zeros(0, dtype=bool)
    c = pos[sel]
    m = ((c[:, 0] > lo[0] - margin) & (c[:, 0] < hi[0] + margin) &
         (c[:, 1] > lo[1] - margin) & (c[:, 1] < hi[1] + margin) &
         (c[:, 2] > lo[2] - margin) & (c[:, 2] < hi[2] + margin))
    return sel[m], m


def frame(pos):
    """Geometric observables for one frame. pos in nm."""
    P = pos[pep]
    zc = float(P[:, 2].mean())
    zp = pos[phos, 2]
    # The peptide starts above the bilayer, so the relevant interface is the
    # upper leaflet: phosphates above the bilayer midplane.
    mid = float(zp.mean())
    upper = zp[zp > mid]
    z_up = float(upper.mean()) if len(upper) else mid

    box = (P.min(axis=0), P.max(axis=0))
    lsel, lmask = _near(pos, lipid, box)
    if len(lsel):
        d = np.linalg.norm(pos[lsel][None, :, :] - P[:, None, :], axis=2)
        dmin = float(d.min())
        close = (d < 0.4).any(axis=0)
        ps_mask = is_pops[lmask]
        n_ps = int(close[ps_mask].sum())
        n_pc = int(close[~ps_mask].sum())
    else:
        dmin, n_pc, n_ps = float('inf'), 0, 0

    wsel, _ = _near(pos, wat_o, box, margin=0.5)
    if len(wsel):
        dw = np.linalg.norm(pos[wsel][None, :, :] - P[:, None, :], axis=2)
        n_wat = int((dw.min(axis=0) < 0.4).sum())
    else:
        n_wat = 0

    return {'z_com': zc, 'z_rel': zc - z_up, 'min_dist': dmin,
            'contacts_popc': n_pc, 'contacts_pops': n_ps, 'n_wat_4A': n_wat}


n_steps = int(PROD_NS * 1e6 / (DT * 1000))
every = int(SAMPLE_PS / DT)
n_frames = n_steps // every
print(f"Production {PROD_NS} ns, {n_frames} frames", flush=True)

traj = []
t0 = time.time()
for i in range(n_frames):
    sim.step(every)
    pos = sim.context.getState(getPositions=True).getPositions(
        asNumpy=True).value_in_unit(unit.nanometers)
    f = frame(pos)
    f['t_ns'] = (i + 1) * SAMPLE_PS / 1000
    traj.append(f)
    if i % 50 == 0:
        el = time.time() - t0
        print(f"  {f['t_ns']:6.1f} ns  z_rel={f['z_rel']:+.2f} nm  "
              f"dmin={f['min_dist']:.2f} nm  PC={f['contacts_popc']:3d} "
              f"PS={f['contacts_pops']:3d}  "
              f"{f['t_ns']/(el/86400) if el else 0:.0f} ns/day", flush=True)

# --- summary ---------------------------------------------------------------
dmin = np.array([f['min_dist'] for f in traj])
zrel = np.array([f['z_rel'] for f in traj])
cpc = np.array([f['contacts_popc'] for f in traj])
cps = np.array([f['contacts_pops'] for f in traj])
tns = np.array([f['t_ns'] for f in traj])

BOUND = 0.4          # nm; any heavy-atom pair inside 4 A counts as contact
bound = dmin < BOUND
first = float(tns[bound][0]) if bound.any() else None
# POPS is 20% of the lipids; enrichment > 1 means the peptide prefers the
# anionic lipid over what random contact would give.
n_ps_atoms, n_pc_atoms = len(pops), len(popc)
tot = cpc + cps
enrich = float(np.mean((cps[tot > 0] / tot[tot > 0]) / (n_ps_atoms / (n_ps_atoms + n_pc_atoms)))) \
    if (tot > 0).any() else None

result = {
    'system': system_name, 'replica': replica, 'seed': SEED,
    'prod_ns': PROD_NS, 'equil_ns': EQUIL_NS,
    'n_peptide_atoms': len(pep),
    'bound_fraction': float(bound.mean()),
    'first_contact_ns': first,
    'min_dist_min_nm': float(dmin.min()),
    'min_dist_final_nm': float(dmin[-1]),
    'z_rel_min_nm': float(zrel.min()),
    'z_rel_final_nm': float(zrel[-1]),
    'z_rel_mean_last10ns': float(zrel[tns > tns.max() - 10].mean()),
    'contacts_popc_mean': float(cpc.mean()),
    'contacts_pops_mean': float(cps.mean()),
    'pops_enrichment': enrich,
    'pops_atom_fraction': n_ps_atoms / (n_ps_atoms + n_pc_atoms),
    'trajectory': traj,
}
with open(f"{OUT_DIR}/{tag}.json", 'w') as f:
    json.dump(result, f, indent=2)

print(f"\n{tag}")
print(f"  bound {bound.mean()*100:.0f}% of frames, first contact "
      f"{first if first is not None else 'never'} ns")
print(f"  deepest z_rel {zrel.min():+.2f} nm, final {zrel[-1]:+.2f} nm")
print(f"  contacts POPC {cpc.mean():.0f} / POPS {cps.mean():.0f}, "
      f"POPS enrichment {enrich if enrich else float('nan'):.2f}x")
print(f"  wrote {OUT_DIR}/{tag}.json  [{(time.time()-t0)/60:.0f} min]")
