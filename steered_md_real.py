"""Steered MD on REAL POPC:POPS membrane from CHARMM-GUI."""
import openmm as mm
import openmm.app as app
import openmm.unit as unit
import numpy as np
import json, sys, os, time

# Support both local and Vast.ai paths
for path in ['/workspace/charmm_gui/openmm',
             '/workspace/prion-neurotoxicity/charmm_gui_membrane/charmm-gui-8797249059/openmm']:
    if os.path.isdir(path):
        os.chdir(path)
        break
else:
    print("ERROR: Cannot find openmm directory")
    sys.exit(1)

psf = app.CharmmPsfFile('step5_input.psf')
pdb = app.PDBFile('step5_input.pdb')

# Set box vectors from sysinfo.dat
with open('sysinfo.dat') as f:
    _sysinfo = json.load(f)
_dims = _sysinfo['dimensions']
psf.setBox(_dims[0]*unit.angstroms, _dims[1]*unit.angstroms, _dims[2]*unit.angstroms)
print(f'Box: {_dims[0]:.1f} x {_dims[1]:.1f} x {_dims[2]:.1f} A')

import glob as _glob
_toppar_dir = os.path.join(os.path.dirname(os.getcwd()), 'toppar')
if os.path.isdir(_toppar_dir):
    _pfiles = sorted(_glob.glob(os.path.join(_toppar_dir, '*.rtf'))) + \
              sorted(_glob.glob(os.path.join(_toppar_dir, '*.prm'))) + \
              sorted(_glob.glob(os.path.join(_toppar_dir, '*.str'))) + ['toppar.str']
    params = app.CharmmParameterSet(*_pfiles)
    print(f'Loaded {len(_pfiles)} parameter files from toppar/')
else:
    params = app.CharmmParameterSet('toppar.str')
    print('Using toppar.str only')

system = psf.createSystem(params,
    nonbondedMethod=app.PME,
    nonbondedCutoff=1.2*unit.nanometers,
    switchDistance=1.0*unit.nanometers,
    constraints=app.HBonds)

# Barostat for membrane
baro = mm.MonteCarloMembraneBarostat(
    1.0*unit.bar, 0.0*unit.bar*unit.nanometers, 303.15*unit.kelvin,
    mm.MonteCarloMembraneBarostat.XYIsotropic,
    mm.MonteCarloMembraneBarostat.ZFree)
system.addForce(baro)

# Identify peptide atoms by residue name (protein amino acids)
PROTEIN_RES = {'ALA','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS','HSD','HSE','HSP',
               'ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'}
EXCLUDE_RES = {'TIP3','HOH','SOD','CLA','POT','POPC','POPS','NA','CL','K','WAT'}
pep_atoms = [a.index for a in psf.topology.atoms()
             if a.residue.name in PROTEIN_RES and a.residue.name not in EXCLUDE_RES]
n_pep = len(pep_atoms)
print(f'Peptide atoms: {n_pep}')

# Steered MD: pull peptide from current z toward membrane center (z=0)
PULL_RATE = 0.0003  # nm/ps (slower for real membrane)
K_PULL = 300.0      # kJ/mol/nm^2 per atom
TOTAL_NS = 15.0     # longer for real membrane
DT = 0.002

# Get initial peptide z
sim_init = app.Simulation(psf.topology, system,
    mm.LangevinMiddleIntegrator(303.15*unit.kelvin, 1.0/unit.picosecond, DT*unit.picoseconds),
    mm.Platform.getPlatformByName('CUDA'), {'Precision':'mixed'})

# Load or create equilibrated state
rst_file = 'step7_production.rst'
if os.path.exists(rst_file):
    sim_init.loadState(rst_file)
    print(f'Loaded equilibrated state from {rst_file}')
else:
    print('No restart file — running quick equilibration from PDB...')
    sim_init.context.setPositions(pdb.positions)
    sim_init.minimizeEnergy(maxIterations=5000)
    print('  Minimized. Running 1 ns NVT equilibration...')
    sim_init.context.setVelocitiesToTemperature(303.15*unit.kelvin)
    sim_init.step(500000)  # 1 ns at 2 fs
    print('  Equilibrated 1 ns.')

state = sim_init.context.getState(getPositions=True)
pos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometers)
z_start = float(np.mean(pos[pep_atoms, 2]))

# Find membrane center (average z of phosphorus atoms)
p_atoms = [a.index for a in psf.topology.atoms() if a.name == 'P']
if p_atoms:
    mem_center = float(np.mean(pos[p_atoms, 2]))
    mem_upper = float(np.max(pos[p_atoms, 2]))
    print(f'Membrane center: {mem_center:.3f} nm, upper: {mem_upper:.3f} nm')
else:
    mem_center = 0.0
    mem_upper = 2.0
    print('No P atoms found, using z=0 as center')

z_end = mem_center  # pull toward membrane center
print(f'Peptide start z: {z_start:.3f}, target: {z_end:.3f}')
del sim_init

# Add pulling force
pull = mm.CustomExternalForce('0.5*k_smd*(z-z_smd_target)^2')
pull.addGlobalParameter('k_smd', K_PULL / n_pep)
pull.addGlobalParameter('z_smd_target', z_start)
for idx in pep_atoms:
    pull.addParticle(idx, [])
system.addForce(pull)

integ = mm.LangevinMiddleIntegrator(303.15*unit.kelvin, 1.0/unit.picosecond, DT*unit.picoseconds)
sim = app.Simulation(psf.topology, system, integ,
    mm.Platform.getPlatformByName('CUDA'), {'Precision':'mixed'})

if os.path.exists(rst_file):
    sim.loadState(rst_file)
else:
    sim.context.setPositions(pdb.positions)
    sim.minimizeEnergy()

sim.context.setVelocitiesToTemperature(303.15*unit.kelvin)

n_total = int(TOTAL_NS * 1e6 / (DT * 1000))
collect_every = 500
z_targets = np.linspace(z_start, z_end, n_total // collect_every)

out_dir = '/workspace/steered_real'
os.makedirs(out_dir, exist_ok=True)
z_data = []

print(f'Running {TOTAL_NS} ns steered MD...')
t0 = time.time()
for si, zt in enumerate(z_targets):
    sim.context.setParameter('z_smd_target', zt)
    sim.step(collect_every)
    st = sim.context.getState(getPositions=True)
    p = st.getPositions(asNumpy=True).value_in_unit(unit.nanometers)
    pz = float(np.mean(p[pep_atoms, 2]))
    dz = pz - zt
    z_data.append({'t':float(si*collect_every*DT/1000), 'zt':float(zt), 'z':float(pz), 'dz':float(dz)})
    if si % 200 == 0:
        elapsed = time.time()-t0
        ns = si*collect_every*DT/1000
        spd = ns/(elapsed/86400) if elapsed>0 else 0
        print(f'  t={ns:.1f}ns zt={zt:.3f} z={pz:.3f} dz={dz:+.3f} {spd:.0f}ns/day')

cum_work = np.cumsum([d['dz']*K_PULL*PULL_RATE*DT for d in z_data])
results = {
    'system': 'KKRPKP + POPC:POPS 80:20 (CHARMM-GUI, real membrane)',
    'membrane_center': mem_center, 'z_start': z_start, 'z_end': z_end,
    'pull_rate': PULL_RATE, 'k_pull': K_PULL, 'total_ns': TOTAL_NS,
    'total_work_kj': float(cum_work[-1]),
    'total_work_kcal': float(cum_work[-1]/4.184),
    'trajectory': z_data[:100],  # subsample
}
with open(f'{out_dir}/results.json','w') as f:
    json.dump(results, f, indent=2)
print(f'\nTotal work: {cum_work[-1]:.2f} kJ/mol ({cum_work[-1]/4.184:.2f} kcal/mol)')
print(f'Done in {(time.time()-t0)/60:.1f} min')
print(f'Saved to {out_dir}/results.json')
