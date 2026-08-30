"""Steered MD: pull peptide toward membrane surface."""
import openmm as mm
import openmm.app as app
import openmm.unit as unit
import numpy as np
import json, sys, os, time

name = sys.argv[1]
pdb = app.PDBFile(f'/workspace/{name}_system.pdb')
with open(f'/workspace/{name}_system.xml') as f:
    system = mm.XmlSerializer.deserialize(f.read())

pep_atoms = [a.index for a in pdb.topology.atoms() if a.residue.name not in ('HOH','NA','CL')]
n_pep = len(pep_atoms)

# Steered MD parameters
PULL_RATE = 0.0005  # nm/ps = 0.5 nm/ns (slow pull)
K_PULL = 500.0      # kJ/mol/nm^2 (spring constant, per atom)
TOTAL_NS = 10.0     # total simulation time
DT = 0.002          # ps
Z_START = 2.5       # starting z (nm)
Z_END = 0.3         # target z (nm)

# Add pulling force: harmonic spring with moving center
pull = mm.CustomExternalForce('0.5*k_pull*(z-z_target)^2')
pull.addGlobalParameter('k_pull', K_PULL / n_pep)
pull.addGlobalParameter('z_target', Z_START)
for idx in pep_atoms:
    pull.addParticle(idx, [])
system.addForce(pull)

integ = mm.LangevinMiddleIntegrator(300*unit.kelvin, 1.0/unit.picosecond, DT*unit.picoseconds)
plat = mm.Platform.getPlatformByName('CUDA')
sim = app.Simulation(pdb.topology, system, integ, plat, {'Precision':'mixed'})
sim.context.setPositions(pdb.positions)
sim.context.setVelocitiesToTemperature(300*unit.kelvin)
sim.minimizeEnergy(maxIterations=2000)

n_total = int(TOTAL_NS * 1e6 / (DT * 1000))
collect_every = 500  # every 1 ps
n_steps_per_update = collect_every

z_target_per_step = np.linspace(Z_START, Z_END, n_total // n_steps_per_update)

results = {'name': name, 'pull_rate': PULL_RATE, 'k_pull': K_PULL}
z_data = []
force_data = []

print(f'Steered MD: {name}')
print(f'Pull from z={Z_START} to z={Z_END} nm over {TOTAL_NS} ns')
print(f'Rate: {PULL_RATE} nm/ps, k={K_PULL} kJ/mol/nm^2')

t0 = time.time()
for step_i, zt in enumerate(z_target_per_step):
    sim.context.setParameter('z_target', zt)
    sim.step(n_steps_per_update)

    state = sim.context.getState(getPositions=True, getEnergy=True)
    pos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometers)
    pep_z = float(np.mean(pos[pep_atoms, 2]))

    # Force = k * (z - z_target) -- work done by spring
    dz = pep_z - zt
    force = K_PULL * dz  # total force
    work_inst = 0.5 * K_PULL * dz**2  # instantaneous work

    z_data.append({'t': float(step_i * n_steps_per_update * DT / 1000),
                   'z_target': float(zt), 'z_actual': float(pep_z),
                   'force': float(force), 'work': float(work_inst)})

    if step_i % 200 == 0:
        elapsed = time.time() - t0
        ns_done = step_i * n_steps_per_update * DT / 1000
        speed = ns_done / (elapsed/86400) if elapsed > 0 else 0
        print(f'  t={ns_done:.1f}ns z_target={zt:.3f} z_actual={pep_z:.3f} dz={dz:+.3f} speed={speed:.0f}ns/day')

# Compute PMF from work
times = [d['t'] for d in z_data]
z_targets = [d['z_target'] for d in z_data]
z_actuals = [d['z_actual'] for d in z_data]
forces = [d['force'] for d in z_data]
works = [d['work'] for d in z_data]

# Cumulative work (Jarzynski-like, but single trajectory)
cum_work = np.cumsum([d['force'] * PULL_RATE * DT for d in z_data])

results['trajectory'] = z_data
results['cumulative_work_kj'] = cum_work.tolist()
results['total_work_kj'] = float(cum_work[-1])
results['total_work_kcal'] = float(cum_work[-1] / 4.184)

print(f'\nTotal work: {cum_work[-1]:.2f} kJ/mol ({cum_work[-1]/4.184:.2f} kcal/mol)')
print(f'Positive = unfavorable insertion; Negative = favorable insertion')

out_dir = f'/workspace/steered_{name}'
os.makedirs(out_dir, exist_ok=True)
with open(f'{out_dir}/results.json', 'w') as f:
    json.dump(results, f, indent=2)
print(f'Saved to {out_dir}/results.json')
print(f'Done in {(time.time()-t0)/60:.1f} min')
