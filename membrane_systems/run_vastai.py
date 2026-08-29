"""Run membrane-mimetic simulation on Vast.ai GPU."""
import openmm as mm
import openmm.app as app
import openmm.unit as unit
import sys, os, time

name = sys.argv[1] if len(sys.argv) > 1 else "KKRPKP_membrane"
n_ns = int(sys.argv[2]) if len(sys.argv) > 2 else 300

pdb = app.PDBFile(f"{name}_system.pdb")
with open(f"{name}_system.xml") as f:
    system = mm.XmlSerializer.deserialize(f.read())

integrator = mm.LangevinMiddleIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 0.002*unit.picoseconds)

for pname in ['CUDA', 'OpenCL', 'CPU']:
    try:
        platform = mm.Platform.getPlatformByName(pname)
        props = {'Precision': 'mixed'} if pname == 'CUDA' else {}
        print(f"Platform: {pname}")
        break
    except Exception:
        continue

sim = app.Simulation(pdb.topology, system, integrator, platform, props if pname == 'CUDA' else {})
sim.context.setPositions(pdb.positions)

chk = f"{name}_checkpoint.chk"
if os.path.exists(chk):
    sim.loadCheckpoint(chk)
    print("Resumed from checkpoint")
else:
    sim.minimizeEnergy()
    sim.context.setVelocitiesToTemperature(300*unit.kelvin)

n_steps = int(n_ns * 500000)
sim.reporters.append(app.DCDReporter(f"{name}_traj.dcd", 5000))
sim.reporters.append(app.StateDataReporter(f"{name}_energy.csv", 5000,
    step=True, time=True, potentialEnergy=True, temperature=True, speed=True))
sim.reporters.append(app.CheckpointReporter(chk, 5000000))

print(f"Running {n_ns} ns ({n_steps:,} steps)...")
t0 = time.time()
chunk = 500000
done = 0
while done < n_steps:
    r = min(chunk, n_steps - done)
    sim.step(r)
    done += r
    ns_done = done * 0.002 / 1000
    s = sim.context.getState(getEnergy=True, getPositions=True)
    elapsed = time.time() - t0
    speed = ns_done / (elapsed / 86400) if elapsed > 0 else 0
    positions = s.getPositions(asNumpy=True).value_in_unit(unit.nanometers)
    pep_atoms = [a.index for a in pdb.topology.atoms() if a.residue.name not in ('HOH', 'NA', 'CL')]
    pep_z = float(positions[pep_atoms, 2].mean()) if pep_atoms else 0
    energy = s.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    print(f"  {ns_done:.0f}/{n_ns} ns | {speed:.0f} ns/day | E={energy:.0f} | z={pep_z:.2f}")

pos = sim.context.getState(getPositions=True).getPositions()
with open(f"{name}_final.pdb", "w") as f:
    app.PDBFile.writeFile(pdb.topology, pos, f)
print(f"Done in {(time.time()-t0)/3600:.1f}h")
