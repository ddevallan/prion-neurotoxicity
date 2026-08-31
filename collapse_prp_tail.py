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

# Hydrogens are pre-added (PrP_23_93_H.pdb) so this runs without pdbfixer,
# which conflicts with the openmm 8.1.1 required by openmm-cuda.
SRC_H = SRC.replace(".pdb", "_H.pdb")
pdb = app.PDBFile(SRC_H if os.path.exists(SRC_H) else SRC)
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


# Target Rg window. Flory scaling for a disordered chain, Rg ~ 2.54*N^0.522,
# gives ~24 A for 71 residues. Implicit solvent (GB) is known to over-compact
# IDPs, so we stop in the physiological window instead of running to the most
# compact structure, which would artificially bury the residues that should be
# available to contact the membrane.
N_RES = 71
RG_TARGET = 2.54 * N_RES ** 0.522
RG_LO, RG_HI = RG_TARGET * 0.85, RG_TARGET * 1.15
print(f"Flory Rg for {N_RES}-residue IDP: {RG_TARGET:.1f} A "
      f"(accept {RG_LO:.1f}-{RG_HI:.1f} A, need max extent < 60 A)")

report("start")
t0 = time.time()
TOTAL_NS = 20.0
chunk = 50000  # 0.1 ns — finer sampling so we can catch the target window
done = 0
n_steps = int(TOTAL_NS * 1e6 / 2)

chosen_state, chosen_rg = None, None
fallback_state, fallback_rg = None, 1e9

while done < n_steps:
    sim.step(min(chunk, n_steps - done))
    done += chunk
    ns = done * 0.002 / 1000
    ext, rg, st = report(f"{ns:.1f} ns")
    if rg < fallback_rg:
        fallback_state, fallback_rg = st, rg
    # Accept the first frame that is both compact enough to fit the membrane
    # patch and still within the physiological IDP size range.
    if max(ext) < 60 and RG_LO <= rg <= RG_HI:
        chosen_state, chosen_rg = st, rg
        print(f"  Rg {rg:.1f} A in target window and fits box — accepted.")
        break
    if rg < RG_LO:
        print(f"  Rg {rg:.1f} A dropped below window (over-collapsed) — stopping.")
        break

speed = (done * 0.002 / 1000) / ((time.time() - t0) / 86400)
print(f"Speed: {speed:.0f} ns/day")

if chosen_state is not None:
    state, final_rg, note = chosen_state, chosen_rg, "within physiological Rg window"
else:
    state, final_rg, note = fallback_state, fallback_rg, "FALLBACK: target window never hit"

with open(OUT, 'w') as f:
    app.PDBFile.writeFile(pdb.topology, state.getPositions(), f)

pos = state.getPositions(asNumpy=True).value_in_unit(unit.angstroms)
ext = pos.max(axis=0) - pos.min(axis=0)
print(f"\nSaved {OUT}  [{note}]")
print(f"Final extent: {ext[0]:.1f} x {ext[1]:.1f} x {ext[2]:.1f} A (Rg {final_rg:.1f} A)")
print(f"Flory reference Rg: {RG_TARGET:.1f} A")
print(f"Fits in 66 A box: {'YES' if max(ext) < 60 else 'NO — needs bigger patch'}")
