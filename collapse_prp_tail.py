"""Collapse PrP 23-93 from the extended AlphaFold conformation into a
compact globule, which is what an IDP actually looks like in water.

Implicit solvent (GBn2) so it is fast and the collapse is driven by
real solvation physics rather than an artificial restraint.
"""
import openmm as mm
import openmm.app as app
import openmm.unit as unit
import numpy as np
import time, os

BASE = None
for cand in ["/workspace/prion-neurotoxicity", "/Users/allan/Projects/cjd"]:
    if os.path.isdir(os.path.join(cand, "calibration_peptides")):
        BASE = cand
        break
SRC = os.path.join(BASE, "calibration_peptides/PrP_23_93.pdb")
OUT = os.path.join(BASE, "calibration_peptides/PrP_23_93_compact.pdb")

from pdbfixer import PDBFixer
fixer = PDBFixer(filename=SRC)
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(7.4)

tmp = "/tmp/prp_tail_fixed.pdb"
with open(tmp, 'w') as f:
    app.PDBFile.writeFile(fixer.topology, fixer.positions, f)

pdb = app.PDBFile(tmp)
ff = app.ForceField('amber14-all.xml', 'implicit/gbn2.xml')
system = ff.createSystem(pdb.topology,
                         nonbondedMethod=app.NoCutoff,
                         constraints=app.HBonds,
                         soluteDielectric=1.0,
                         solventDielectric=78.5)

# 350 K accelerates the collapse without unfolding anything that is folded
# (this region has no stable secondary structure to lose).
integrator = mm.LangevinMiddleIntegrator(350*unit.kelvin,
                                         1.0/unit.picosecond,
                                         0.002*unit.picoseconds)
try:
    platform = mm.Platform.getPlatformByName('CUDA')
    print("Platform: CUDA")
except Exception:
    platform = mm.Platform.getPlatformByName('CPU')
    platform.setPropertyDefaultValue('Threads', '4')
    print("Platform: CPU (4 threads)")

sim = app.Simulation(pdb.topology, system, integrator, platform)
sim.context.setPositions(pdb.positions)

n_atoms = system.getNumParticles()
print(f"PrP 23-93: {n_atoms} atoms, implicit solvent (GBn2)")
print("Minimizing...")
sim.minimizeEnergy(maxIterations=2000)
sim.context.setVelocitiesToTemperature(350*unit.kelvin)


def report(label):
    st = sim.context.getState(getPositions=True)
    pos = st.getPositions(asNumpy=True).value_in_unit(unit.angstroms)
    ext = pos.max(axis=0) - pos.min(axis=0)
    com = pos.mean(axis=0)
    rg = np.sqrt(((pos - com) ** 2).sum(axis=1).mean())
    print(f"  {label}: extent {ext[0]:.0f} x {ext[1]:.0f} x {ext[2]:.0f} A, Rg = {rg:.1f} A")
    return ext, rg, st


report("start")
t0 = time.time()
TOTAL_NS = 20.0
chunk = 250000  # 0.5 ns
done = 0
n_steps = int(TOTAL_NS * 1e6 / 2)

best_rg = 1e9
best_state = None

while done < n_steps:
    sim.step(min(chunk, n_steps - done))
    done += chunk
    ns = done * 0.002 / 1000
    ext, rg, st = report(f"{ns:.1f} ns")
    if rg < best_rg:
        best_rg = rg
        best_state = st
    # A 71-residue IDP has Rg ~ 22-26 A; stop once compact enough to fit the box
    if max(ext) < 60:
        print("  Compact enough — stopping early.")
        break

speed = (done * 0.002 / 1000) / ((time.time() - t0) / 86400)
print(f"Speed: {speed:.0f} ns/day")

state = best_state if best_state is not None else sim.context.getState(getPositions=True)
with open(OUT, 'w') as f:
    app.PDBFile.writeFile(pdb.topology, state.getPositions(), f)

pos = state.getPositions(asNumpy=True).value_in_unit(unit.angstroms)
ext = pos.max(axis=0) - pos.min(axis=0)
print(f"\nSaved {OUT}")
print(f"Final extent: {ext[0]:.1f} x {ext[1]:.1f} x {ext[2]:.1f} A (Rg {best_rg:.1f} A)")
print(f"Fits in 66 A box: {'YES' if max(ext) < 60 else 'NO — needs bigger patch'}")
