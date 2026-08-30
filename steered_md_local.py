"""Steered MD on REAL POPC:POPS membrane — local M1 Pro (OpenCL/Metal)."""
import openmm as mm
import openmm.app as app
import openmm.unit as unit
import numpy as np
import json, sys, os, time

name = sys.argv[1] if len(sys.argv) > 1 else "KKRPKP"
base = "/Users/allan/Projects/cjd/charmm_gui_membrane/charmm-gui-8797249059/openmm"
os.chdir(base)

psf = app.CharmmPsfFile('step5_input.psf')
crd = app.CharmmCrdFile('step5_input.crd')
pdb = app.PDBFile('step5_input.pdb')
import glob
toppar_dir = os.path.join(os.path.dirname(base), 'toppar')
param_files = sorted(glob.glob(os.path.join(toppar_dir, '*.rtf'))) + \
              sorted(glob.glob(os.path.join(toppar_dir, '*.prm'))) + \
              sorted(glob.glob(os.path.join(toppar_dir, '*.str'))) + \
              ['toppar.str']
params = app.CharmmParameterSet(*param_files)

# Set periodic box from sysinfo.dat (JSON)
with open('sysinfo.dat') as f:
    sysinfo = json.load(f)
dims = sysinfo['dimensions']
psf.setBox(dims[0]*unit.angstroms, dims[1]*unit.angstroms, dims[2]*unit.angstroms)
print(f'Box: {dims[0]:.1f} x {dims[1]:.1f} x {dims[2]:.1f} A')

system = psf.createSystem(params,
    nonbondedMethod=app.PME,
    nonbondedCutoff=1.2*unit.nanometers,
    switchDistance=1.0*unit.nanometers,
    constraints=app.HBonds)

baro = mm.MonteCarloMembraneBarostat(
    1.0*unit.bar, 0.0*unit.bar*unit.nanometers, 303.15*unit.kelvin,
    mm.MonteCarloMembraneBarostat.XYIsotropic,
    mm.MonteCarloMembraneBarostat.ZFree)
system.addForce(baro)

# Find peptide atoms (PRO segments)
pep_atoms = []
for a in psf.topology.atoms():
    seg = a.residue.segment_id if hasattr(a.residue, 'segment_id') else ''
    chain = a.residue.chain.id
    if chain in ('A','B','C','D','E','F') and a.residue.name not in ('TIP3','POT','CLA','POPC','POPS'):
        pep_atoms.append(a.index)
if not pep_atoms:
    for a in psf.topology.atoms():
        if a.residue.name in ('LYS','ARG','PRO','ASN','GLY','TRP','THR','SER'):
            if a.index < 200:
                pep_atoms.append(a.index)
n_pep = len(pep_atoms)
print(f'{name}: {n_pep} peptide atoms, {system.getNumParticles()} total')

# Use CPU — OpenCL context creation fails in uv venv on M1
platform = mm.Platform.getPlatformByName('CPU')
prop = {}
import multiprocessing
n_threads = multiprocessing.cpu_count()
platform.setPropertyDefaultValue('Threads', str(n_threads))
print(f'Using CPU with {n_threads} threads')

DT = 0.002
integ_init = mm.LangevinMiddleIntegrator(303.15*unit.kelvin, 1.0/unit.picosecond, DT*unit.picoseconds)
sim = app.Simulation(psf.topology, system, integ_init, platform, prop)
sim.context.setPositions(pdb.positions)

print('Minimizing (5000 iter)...')
sim.minimizeEnergy(maxIterations=5000)
print('Equilibrating 0.5 ns...')
sim.context.setVelocitiesToTemperature(303.15*unit.kelvin)
sim.step(250000)
print('Equilibrated.')

# Get initial positions
state = sim.context.getState(getPositions=True)
pos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometers)
z_start = float(np.mean(pos[pep_atoms, 2]))

p_atoms = [a.index for a in psf.topology.atoms() if a.name == 'P']
if p_atoms:
    mem_center = float(np.mean(pos[p_atoms, 2]))
    print(f'Membrane center: {mem_center:.3f} nm')
else:
    mem_center = float(np.mean(pos[:, 2]))
    print(f'Using system center: {mem_center:.3f} nm')

z_end = mem_center
print(f'Pull: z={z_start:.3f} -> {z_end:.3f} nm')

# Add pulling force
PULL_RATE = 0.0003
K_PULL = 300.0
TOTAL_NS = 10.0

pull = mm.CustomExternalForce('0.5*k_smd_local*(z-z_smd_local_target)^2')
pull.addGlobalParameter('k_smd_local', K_PULL / n_pep)
pull.addGlobalParameter('z_smd_local_target', z_start)
for idx in pep_atoms:
    pull.addParticle(idx, [])
system.addForce(pull)

# Recreate simulation with pulling force
del sim
integ = mm.LangevinMiddleIntegrator(303.15*unit.kelvin, 1.0/unit.picosecond, DT*unit.picoseconds)
sim = app.Simulation(psf.topology, system, integ, platform, prop)
sim.context.setPositions(state.getPositions())
sim.context.setVelocitiesToTemperature(303.15*unit.kelvin)

n_total = int(TOTAL_NS * 1e6 / (DT * 1000))
collect_every = 500
z_targets = np.linspace(z_start, z_end, n_total // collect_every)

out_dir = f'/Users/allan/Projects/cjd/results_steered_real'
os.makedirs(out_dir, exist_ok=True)
z_data = []

print(f'Running {TOTAL_NS} ns steered MD on OpenCL/Metal...')
t0 = time.time()
for si, zt in enumerate(z_targets):
    sim.context.setParameter('z_smd_local_target', zt)
    sim.step(collect_every)
    st = sim.context.getState(getPositions=True)
    p = st.getPositions(asNumpy=True).value_in_unit(unit.nanometers)
    pz = float(np.mean(p[pep_atoms, 2]))
    dz = pz - zt
    z_data.append({'t': float(si*collect_every*DT/1000), 'zt': float(zt), 'z': float(pz), 'dz': float(dz)})
    if si % 200 == 0:
        elapsed = time.time() - t0
        ns = si * collect_every * DT / 1000
        spd = ns / (elapsed/86400) if elapsed > 0 else 0
        print(f'  t={ns:.1f}ns zt={zt:.3f} z={pz:.3f} dz={dz:+.3f} {spd:.0f}ns/day')

cum_work = np.cumsum([d['dz'] * K_PULL * PULL_RATE * DT for d in z_data])
results = {
    'system': f'{name} + POPC:POPS 80:20 (CHARMM-GUI, REAL membrane)',
    'platform': 'OpenCL/Metal (M1 Pro)',
    'membrane_center': mem_center, 'z_start': z_start, 'z_end': z_end,
    'pull_rate': PULL_RATE, 'k_pull': K_PULL, 'total_ns': TOTAL_NS,
    'total_work_kj': float(cum_work[-1]),
    'total_work_kcal': float(cum_work[-1] / 4.184),
    'n_peptide_atoms': n_pep,
    'trajectory': z_data[::10],
}
with open(f'{out_dir}/{name}_real_membrane.json', 'w') as f:
    json.dump(results, f, indent=2)

print(f'\nTotal work: {cum_work[-1]:.2f} kJ/mol ({cum_work[-1]/4.184:.2f} kcal/mol)')
print(f'Done in {(time.time()-t0)/60:.1f} min')
print(f'Saved to {out_dir}/{name}_real_membrane.json')
