"""Steered MD with independent replicas — the statistical version of
steered_md_real.py.

Two things this fixes relative to the single-trajectory run:

1. WORK INTEGRAL. The restraint moves by (z_end - z_start)/n_windows per
   window, not by PULL_RATE*DT (a nominal rate that never fed the linspace
   that actually generates the targets). And dW = -k*dz*dz0 carries a minus
   sign. The old expression was off by a constant factor of ~222 and had the
   sign folded into the interpretation rather than the math.

2. n=1. A single steered trajectory has no error bar, and the pulling work is
   a fluctuating quantity — comparing two single numbers says nothing about
   whether they differ. Each replica here gets its own velocity seed plus a
   decorrelation run before pulling starts, so the replicas are independent
   samples and the difference between systems can be tested.

Usage:
  python steered_replicas.py KKRPKP 0 0     # system, replica index, gpu index
"""
import openmm as mm
import openmm.app as app
import openmm.unit as unit
import numpy as np
import json, sys, os, time, glob

peptide_name = sys.argv[1] if len(sys.argv) > 1 else 'KKRPKP'
replica = int(sys.argv[2]) if len(sys.argv) > 2 else 0
gpu_index = sys.argv[3] if len(sys.argv) > 3 else '0'

PATHS = {
    'KKRPKP': 'charmm_gui_membrane/charmm-gui-8797249059/openmm',
    'NNRPNP': 'charmm_gui_nnrpnp/charmm-gui-8811352165/openmm',
}
rel = PATHS[peptide_name]
for base in ['/workspace/prion-neurotoxicity', '/Users/allan/Projects/cjd']:
    p = os.path.join(base, rel)
    if os.path.isdir(p):
        os.chdir(p)
        break
else:
    sys.exit(f"ERROR: cannot find system directory for {peptide_name}")

OUT_DIR = '/workspace/replicas' if os.path.isdir('/workspace') else \
          '/Users/allan/Projects/cjd/results_replicas'
os.makedirs(OUT_DIR, exist_ok=True)
tag = f"{peptide_name}_rep{replica}"
print(f"=== {tag} on GPU {gpu_index} ===", flush=True)

# Every replica gets a distinct seed stream; the offset keeps KKRPKP and
# NNRPNP replicas from sharing seeds, which would correlate the comparison.
SEED = 10000 + replica + (0 if peptide_name == 'KKRPKP' else 500)

psf = app.CharmmPsfFile('step5_input.psf')
pdb = app.PDBFile('step5_input.pdb')
with open('sysinfo.dat') as f:
    dims = json.load(f)['dimensions']
psf.setBox(dims[0]*unit.angstroms, dims[1]*unit.angstroms, dims[2]*unit.angstroms)

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

PROTEIN_RES = {'ALA','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS','HSD','HSE','HSP',
               'ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'}
pep_atoms = [a.index for a in psf.topology.atoms() if a.residue.name in PROTEIN_RES]
n_pep = len(pep_atoms)

K_PULL = 300.0    # kJ/mol/nm^2, total across the peptide
TOTAL_NS = 15.0   # matches the original protocol so results are comparable
DECORR_NS = 2.0   # free MD before pulling, to decorrelate replicas
DT = 0.002
COLLECT_EVERY = 500

platform = mm.Platform.getPlatformByName('CUDA')
props = {'Precision': 'mixed', 'DeviceIndex': gpu_index}

# --- Decorrelation: same starting coordinates, independent velocities ------
integ0 = mm.LangevinMiddleIntegrator(303.15*unit.kelvin, 1.0/unit.picosecond,
                                     DT*unit.picoseconds)
integ0.setRandomNumberSeed(SEED)
sim0 = app.Simulation(psf.topology, system, integ0, platform, props)
if os.path.exists('step7_production.rst'):
    sim0.loadState('step7_production.rst')
else:
    sim0.context.setPositions(pdb.positions)
    sim0.minimizeEnergy(maxIterations=5000)
sim0.context.setVelocitiesToTemperature(303.15*unit.kelvin, SEED)
print(f"Decorrelating {DECORR_NS} ns (seed {SEED})...", flush=True)
sim0.step(int(DECORR_NS * 1e6 / (DT * 1000)))

state = sim0.context.getState(getPositions=True, getVelocities=True)
pos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometers)
z_start = float(np.mean(pos[pep_atoms, 2]))
p_atoms = [a.index for a in psf.topology.atoms() if a.name == 'P']
mem_center = float(np.mean(pos[p_atoms, 2])) if p_atoms else 0.0
z_end = mem_center
print(f"z_start={z_start:.3f} nm, membrane center={mem_center:.3f} nm", flush=True)
del sim0

# --- Pulling ---------------------------------------------------------------
pull = mm.CustomExternalForce('0.5*k_smd*(z-z_smd_target)^2')
pull.addGlobalParameter('k_smd', K_PULL / n_pep)
pull.addGlobalParameter('z_smd_target', z_start)
for idx in pep_atoms:
    pull.addParticle(idx, [])
system.addForce(pull)

integ = mm.LangevinMiddleIntegrator(303.15*unit.kelvin, 1.0/unit.picosecond,
                                    DT*unit.picoseconds)
integ.setRandomNumberSeed(SEED + 1)
sim = app.Simulation(psf.topology, system, integ, platform, props)
sim.context.setState(state)

n_total = int(TOTAL_NS * 1e6 / (DT * 1000))
n_windows = n_total // COLLECT_EVERY
z_targets = np.linspace(z_start, z_end, n_windows)
dz0 = (z_end - z_start) / n_windows   # actual restraint displacement per window

print(f"Pulling {TOTAL_NS} ns over {n_windows} windows, dz0={dz0*1000:.4f} pm/window",
      flush=True)
t0 = time.time()
dz_list, traj = [], []
for si, zt in enumerate(z_targets):
    sim.context.setParameter('z_smd_target', zt)
    sim.step(COLLECT_EVERY)
    p = sim.context.getState(getPositions=True).getPositions(
        asNumpy=True).value_in_unit(unit.nanometers)
    pz = float(np.mean(p[pep_atoms, 2]))
    dz = pz - zt
    dz_list.append(dz)
    if si % 50 == 0:
        traj.append({'t': si*COLLECT_EVERY*DT/1000, 'zt': float(zt),
                     'z': pz, 'dz': dz})
    if si % 2000 == 0:
        ns = si*COLLECT_EVERY*DT/1000
        el = time.time() - t0
        print(f"  {ns:5.1f} ns  z={pz:.3f} dz={dz:+.3f}  "
              f"{ns/(el/86400) if el else 0:.0f} ns/day", flush=True)

# W = -k * sum(dz * dz0). Negative work means the peptide ran ahead of the
# restraint toward the membrane, i.e. insertion is favorable.
dz_arr = np.array(dz_list)
work_cum = -K_PULL * np.cumsum(dz_arr * dz0)
work_total = float(work_cum[-1])

inside = float(np.mean(dz_arr < 0))
result = {
    'system': peptide_name, 'replica': replica, 'seed': SEED,
    'work_kj_mol': work_total,
    'work_kcal_mol': work_total / 4.184,
    'dz_mean_nm': float(dz_arr.mean()),
    'dz_std_nm': float(dz_arr.std()),
    'frac_ahead_of_restraint': inside,
    'z_start': z_start, 'z_end': z_end, 'membrane_center': mem_center,
    'k_pull': K_PULL, 'total_ns': TOTAL_NS, 'decorr_ns': DECORR_NS,
    'n_windows': int(n_windows), 'dz0_nm': float(dz0),
    'kT_kj_mol': 8.314e-3 * 303.15,
    'trajectory': traj,
}
with open(f"{OUT_DIR}/{tag}.json", 'w') as f:
    json.dump(result, f, indent=2)

print(f"\n{tag}: W = {work_total:+.2f} kJ/mol "
      f"({work_total/(8.314e-3*303.15):+.1f} kT), "
      f"dz_mean = {dz_arr.mean():+.4f} nm, ahead {inside*100:.0f}% of the time")
print(f"Wrote {OUT_DIR}/{tag}.json  [{(time.time()-t0)/60:.1f} min]")
