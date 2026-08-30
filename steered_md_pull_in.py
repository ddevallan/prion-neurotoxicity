"""Steered MD: pull KKRPKP from WATER into REAL POPC:POPS membrane.
Uses the EXISTING CHARMM-GUI system — just repositions peptide above membrane first.
This gives the FREE ENERGY OF INSERTION."""
import openmm as mm
import openmm.app as app
import openmm.unit as unit
import numpy as np
import json, sys, os, time, glob

# Find openmm dir
for path in ['/workspace/charmm_gui/openmm',
             '/workspace/prion-neurotoxicity/charmm_gui_membrane/charmm-gui-8797249059/openmm',
             '/Users/allan/Projects/cjd/charmm_gui_membrane/charmm-gui-8797249059/openmm']:
    if os.path.isdir(path):
        os.chdir(path)
        break

psf = app.CharmmPsfFile('step5_input.psf')
pdb = app.PDBFile('step5_input.pdb')

with open('sysinfo.dat') as f:
    _sysinfo = json.load(f)
_dims = _sysinfo['dimensions']
psf.setBox(_dims[0]*unit.angstroms, _dims[1]*unit.angstroms, _dims[2]*unit.angstroms)
print(f'Box: {_dims[0]:.1f} x {_dims[1]:.1f} x {_dims[2]:.1f} A')

_toppar_dir = os.path.join(os.path.dirname(os.getcwd()), 'toppar')
if os.path.isdir(_toppar_dir):
    _pfiles = sorted(glob.glob(os.path.join(_toppar_dir, '*.rtf'))) + \
              sorted(glob.glob(os.path.join(_toppar_dir, '*.prm'))) + \
              sorted(glob.glob(os.path.join(_toppar_dir, '*.str'))) + ['toppar.str']
    params = app.CharmmParameterSet(*_pfiles)
    print(f'Loaded {len(_pfiles)} parameter files')
else:
    params = app.CharmmParameterSet('toppar.str')

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

PROTEIN_RES = {'ALA','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS','HSD','HSE','HSP',
               'ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'}
pep_atoms = [a.index for a in psf.topology.atoms() if a.residue.name in PROTEIN_RES]
n_pep = len(pep_atoms)

# Find membrane P atoms
p_atoms = [a.index for a in psf.topology.atoms() if a.name == 'P']
print(f'Peptide: {n_pep} atoms, P atoms: {len(p_atoms)}')

# Setup initial positions — move peptide 3 nm ABOVE membrane upper surface
try:
    platform = mm.Platform.getPlatformByName('CUDA')
    prop = {'Precision': 'mixed'}
    print('Using CUDA')
except Exception:
    platform = mm.Platform.getPlatformByName('CPU')
    prop = {}
    import multiprocessing
    platform.setPropertyDefaultValue('Threads', str(multiprocessing.cpu_count()))
    print(f'Using CPU ({multiprocessing.cpu_count()} threads)')

DT = 0.002
integ = mm.LangevinMiddleIntegrator(303.15*unit.kelvin, 1.0/unit.picosecond, DT*unit.picoseconds)
sim = app.Simulation(psf.topology, system, integ, platform, prop)
sim.context.setPositions(pdb.positions)

state = sim.context.getState(getPositions=True)
pos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometers)

# Find membrane surfaces
p_z = pos[p_atoms, 2]
mem_upper = float(np.max(p_z))
mem_center = float(np.mean(p_z))
mem_lower = float(np.min(p_z))
print(f'Membrane: lower={mem_lower:.2f}, center={mem_center:.2f}, upper={mem_upper:.2f} nm')

# Move peptide to 3 nm above upper surface
pep_com_z = float(np.mean(pos[pep_atoms, 2]))
target_z = mem_upper + 3.0
shift = target_z - pep_com_z
print(f'Moving peptide from z={pep_com_z:.2f} to z={target_z:.2f} (shift={shift:+.2f} nm)')
for idx in pep_atoms:
    pos[idx, 2] += shift
sim.context.setPositions(pos * unit.nanometers)

# Minimize and brief equilibration
print('Minimizing...')
sim.minimizeEnergy(maxIterations=5000)
sim.context.setVelocitiesToTemperature(303.15*unit.kelvin)
print('Equilibrating 0.2 ns...')
sim.step(100000)

# Get new starting z
state = sim.context.getState(getPositions=True)
pos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometers)
z_start = float(np.mean(pos[pep_atoms, 2]))
z_end = mem_upper  # pull to upper membrane surface
print(f'Pull: z={z_start:.2f} -> z={z_end:.2f} nm (distance: {z_start-z_end:.2f} nm)')

# Add pulling force
PULL_RATE = 0.0005  # nm/ps
K_PULL = 500.0
TOTAL_NS = (z_start - z_end) / (PULL_RATE * 1000) + 2.0  # time to reach + buffer
TOTAL_NS = min(TOTAL_NS, 20.0)  # cap at 20 ns

pull = mm.CustomExternalForce('0.5*k_pull_in*(z-z_pull_in_target)^2')
pull.addGlobalParameter('k_pull_in', K_PULL / n_pep)
pull.addGlobalParameter('z_pull_in_target', z_start)
for idx in pep_atoms:
    pull.addParticle(idx, [])
system.addForce(pull)

# Need new simulation with the pulling force
del sim
integ2 = mm.LangevinMiddleIntegrator(303.15*unit.kelvin, 1.0/unit.picosecond, DT*unit.picoseconds)
sim = app.Simulation(psf.topology, system, integ2, platform, prop)
sim.context.setPositions(state.getPositions())
sim.context.setVelocitiesToTemperature(303.15*unit.kelvin)

n_total = int(TOTAL_NS * 1e6 / (DT * 1000))
collect_every = 500
z_targets = np.linspace(z_start, z_end, n_total // collect_every)

out_dir = '/workspace/steered_pull_in' if os.path.isdir('/workspace') else '/Users/allan/Projects/cjd/results_steered_real'
os.makedirs(out_dir, exist_ok=True)
z_data = []

print(f'Steered MD: pulling from water into membrane')
print(f'  {TOTAL_NS:.1f} ns, rate={PULL_RATE} nm/ps, k={K_PULL}')
t0 = time.time()
for si, zt in enumerate(z_targets):
    sim.context.setParameter('z_pull_in_target', zt)
    sim.step(collect_every)
    st = sim.context.getState(getPositions=True)
    p = st.getPositions(asNumpy=True).value_in_unit(unit.nanometers)
    pz = float(np.mean(p[pep_atoms, 2]))
    dz = pz - zt
    force = K_PULL * dz
    z_data.append({'t':float(si*collect_every*DT/1000), 'zt':float(zt), 'z':float(pz), 'dz':float(dz), 'force':float(force)})
    if si % 200 == 0:
        elapsed = time.time()-t0
        ns = si*collect_every*DT/1000
        spd = ns/(elapsed/86400) if elapsed>0 else 0
        print(f'  t={ns:.1f}ns zt={zt:.3f} z={pz:.3f} dz={dz:+.3f} F={force:+.1f} {spd:.0f}ns/day')

# Calculate work = integral of force × displacement
cum_work = np.cumsum([d['force'] * PULL_RATE * DT / n_pep for d in z_data])

results = {
    'system': 'KKRPKP pulling from water INTO POPC:POPS 80:20 membrane',
    'type': 'insertion_free_energy',
    'membrane_upper': mem_upper, 'z_start': z_start, 'z_end': z_end,
    'pull_rate': PULL_RATE, 'k_pull': K_PULL, 'total_ns': TOTAL_NS,
    'total_work_kj': float(cum_work[-1]),
    'total_work_kcal': float(cum_work[-1]/4.184),
    'trajectory': z_data[::5],
}
with open(f'{out_dir}/KKRPKP_pull_in.json','w') as f:
    json.dump(results, f, indent=2)

print(f'\nInsertion work: {cum_work[-1]:.2f} kJ/mol ({cum_work[-1]/4.184:.2f} kcal/mol)')
print(f'Negative = FAVORABLE insertion')
print(f'Done in {(time.time()-t0)/60:.1f} min')
