"""Umbrella sampling: peptide z-distance to membrane surface."""
import openmm as mm
import openmm.app as app
import openmm.unit as unit
import numpy as np
import json, sys, os, time

name = sys.argv[1]
system_pdb = f'/workspace/{name}_system.pdb'
system_xml = f'/workspace/{name}_system.xml'

pdb = app.PDBFile(system_pdb)
with open(system_xml) as f:
    base_system = mm.XmlSerializer.deserialize(f.read())

N_WINDOWS = 20
Z_MIN, Z_MAX = 0.5, 3.0
K_UMB = 1000.0
NS_PER_WIN = 2.0
DT = 0.002

pep_atoms = [a.index for a in pdb.topology.atoms() if a.residue.name not in ('HOH','NA','CL')]
n_pep = len(pep_atoms)
z_centers = np.linspace(Z_MIN, Z_MAX, N_WINDOWS)

results_dir = f'/workspace/umbrella_{name}'
os.makedirs(results_dir, exist_ok=True)
all_z_data = {}

for wi, z0 in enumerate(z_centers):
    print(f'Window {wi+1}/{N_WINDOWS}: z0={z0:.3f}')
    system = mm.XmlSerializer.deserialize(mm.XmlSerializer.serialize(base_system))
    umb = mm.CustomExternalForce('0.5*k_umb_win*(z-z0_umb_win)^2')
    umb.addGlobalParameter('k_umb_win', K_UMB/n_pep)
    umb.addGlobalParameter('z0_umb_win', z0)
    for idx in pep_atoms:
        umb.addParticle(idx, [])
    system.addForce(umb)

    integ = mm.LangevinMiddleIntegrator(300*unit.kelvin, 1.0/unit.picosecond, DT*unit.picoseconds)
    plat = mm.Platform.getPlatformByName('CUDA')
    sim = app.Simulation(pdb.topology, system, integ, plat, {'Precision':'mixed'})
    sim.context.setPositions(pdb.positions)

    state = sim.context.getState(getPositions=True)
    pos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometers)
    shift = z0 - np.mean(pos[pep_atoms, 2])
    for idx in pep_atoms:
        pos[idx, 2] += shift
    sim.context.setPositions(pos*unit.nanometers)
    sim.context.setVelocitiesToTemperature(300*unit.kelvin)
    sim.minimizeEnergy(maxIterations=5000)

    n_equil = int(1.0e6 / (DT*1000))
    sim.step(n_equil)

    n_prod = int(NS_PER_WIN * 1e6 / (DT*1000))
    zvals = []
    t0 = time.time()
    done = 0
    while done < n_prod:
        sim.step(500)
        done += 500
        st = sim.context.getState(getPositions=True)
        p = st.getPositions(asNumpy=True).value_in_unit(unit.nanometers)
        zvals.append(float(np.mean(p[pep_atoms, 2])))

    print(f'  z={np.mean(zvals):.3f}+/-{np.std(zvals):.3f} ({time.time()-t0:.0f}s)')
    all_z_data[wi] = {'z0':float(z0), 'z_mean':float(np.mean(zvals)), 'z_std':float(np.std(zvals)), 'z_values':zvals, 'k':K_UMB}

# PMF
all_z = np.concatenate([d['z_values'] for d in all_z_data.values()])
z_edges = np.linspace(Z_MIN-0.5, Z_MAX+0.5, 51)
z_ctrs = 0.5*(z_edges[:-1]+z_edges[1:])
counts, _ = np.histogram(all_z, bins=z_edges)
counts = np.maximum(counts, 1)
kT = 2.479
pmf = -kT * np.log(counts.astype(float))
pmf -= np.min(pmf)

print('\nPMF:')
for z, g in zip(z_ctrs, pmf):
    if g < 30:
        print(f'  z={z:.2f} PMF={g:.2f} kJ/mol ({g/4.184:.2f} kcal/mol)')

output = {
    'name': name, 'n_windows': N_WINDOWS, 'k': K_UMB, 'ns_per_window': NS_PER_WIN,
    'windows': {str(k):{'z0':v['z0'],'z_mean':v['z_mean'],'z_std':v['z_std'],'n':len(v['z_values'])} for k,v in all_z_data.items()},
    'pmf': {'z':z_ctrs.tolist(), 'kj':pmf.tolist(), 'kcal':(pmf/4.184).tolist()},
}
with open(f'{results_dir}/results.json','w') as f:
    json.dump(output, f, indent=2)
print(f'Saved to {results_dir}/results.json')
