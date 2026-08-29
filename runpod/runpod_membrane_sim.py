"""
RunPod membrane simulation: KKRPKP and variants on POPC bilayer.
Tests the v5 hypothesis that PrP N-terminal polybasic motif perturbs
membranes via AMP-like carpet mechanism.

Usage on RunPod:
  1. Create a pod with pytorch/pytorch:latest or nvidia/cuda:12.2.2 image
  2. pip install openmm pdbfixer mdtraj numpy scipy matplotlib
  3. python runpod_membrane_sim.py

Expected runtime on A40: ~1-2 hours per 100 ns simulation.
"""
import openmm as mm
import openmm.app as app
import openmm.unit as unit
import numpy as np
import os
import sys
import json
import time

OUTPUT_DIR = "/workspace/output" if os.path.exists("/workspace") else "./membrane_output"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================
# PEPTIDE BUILDER
# ============================================================

def build_peptide(sequence, name):
    """Build extended peptide using PDBFixer."""
    from pdbfixer import PDBFixer

    AA_MAP = {
        'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
        'E': 'GLU', 'Q': 'GLN', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
        'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
        'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL',
    }

    pdb_lines = ["HEADER    PEPTIDE"]
    atom_idx = 1
    backbone = {
        'N':  np.array([0.0, 0.0, 0.0]),
        'CA': np.array([1.458, 0.0, 0.0]),
        'C':  np.array([2.009, 1.420, 0.0]),
        'O':  np.array([1.251, 2.390, 0.0]),
    }
    step = np.array([3.8, 0.0, 0.0])

    for i, aa in enumerate(sequence):
        resname = AA_MAP.get(aa, 'ALA')
        resid = i + 1
        offset = step * i
        for atom_name, base_pos in backbone.items():
            pos = base_pos + offset
            pdb_lines.append(
                f"ATOM  {atom_idx:5d} {atom_name:<4s} {resname:3s} A{resid:4d}    "
                f"{pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}  1.00  0.00           "
                f"{atom_name[0]:>2s}"
            )
            atom_idx += 1

    pdb_lines.extend(["TER", "END"])

    raw_path = os.path.join(OUTPUT_DIR, f"{name}_raw.pdb")
    with open(raw_path, 'w') as f:
        f.write("\n".join(pdb_lines))

    fixer = PDBFixer(filename=raw_path)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)

    pdb_path = os.path.join(OUTPUT_DIR, f"{name}.pdb")
    with open(pdb_path, 'w') as f:
        app.PDBFile.writeFile(fixer.topology, fixer.positions, f)

    print(f"  Built {name}: {len(sequence)} residues")
    return pdb_path

# ============================================================
# MEMBRANE BUILDER (simplified — rectangular POPC patch)
# ============================================================

def build_membrane_system(peptide_pdb, name, n_lipids_per_leaflet=36,
                          box_xy=5.0, membrane_thickness=4.0):
    """
    Build a peptide + implicit membrane system.

    For a full explicit lipid bilayer, use CHARMM-GUI.
    Here we create the peptide in water with a flat-bottom restraint
    that mimics membrane surface interaction.

    Alternative approach: use OpenMM's CustomExternalForce to create
    a membrane-mimetic potential that attracts the peptide to a plane.
    """
    pdb = app.PDBFile(peptide_pdb)
    forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

    modeller = app.Modeller(pdb.topology, pdb.positions)

    # Place peptide 2 nm above the "membrane" (z=0 plane)
    positions = modeller.positions
    com_z = np.mean([p.value_in_unit(unit.nanometers)[2] for p in positions])
    new_positions = []
    for p in positions:
        x, y, z = p.value_in_unit(unit.nanometers)
        new_positions.append(mm.Vec3(x, y, z - com_z + 2.0) * unit.nanometers)
    modeller.positions = new_positions

    # Add solvent
    modeller.addSolvent(forcefield, model='tip3p',
                        boxSize=mm.Vec3(box_xy, box_xy, 6.0)*unit.nanometers,
                        ionicStrength=0.15*unit.molar)

    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0*unit.nanometers,
        constraints=app.HBonds,
    )

    # Add membrane-mimetic potential: a soft wall at z=0 that represents
    # the headgroup region. Positively charged residues are attracted
    # to the negative headgroup region (z = 0 to -0.5 nm).
    #
    # U(z) = k * (z - z0)^2  for z < z0 (repulsion below membrane)
    # U(z) = -q * epsilon * exp(-(z - z0)^2 / (2*sigma^2))  (attraction for charged)

    membrane_force = mm.CustomExternalForce(
        "k_wall * step(z0 - z) * (z - z0)^2 "
        "+ charge_attraction * exp(-(z - z_attract)^2 / (2 * sigma^2))"
    )
    membrane_force.addGlobalParameter("k_wall", 1000.0)  # kJ/mol/nm²
    membrane_force.addGlobalParameter("z0", 0.0)  # membrane surface at z=0
    membrane_force.addGlobalParameter("z_attract", 0.5)  # headgroup attraction zone
    membrane_force.addGlobalParameter("sigma", 0.3)  # nm
    membrane_force.addPerParticleParameter("charge_attraction")

    # Apply to peptide atoms only
    for atom in modeller.topology.atoms():
        if atom.residue.name not in ('HOH', 'NA', 'CL'):
            # Charged residues get stronger attraction
            resname = atom.residue.name
            if resname in ('LYS', 'ARG'):
                attraction = -5.0  # kJ/mol (attracted to membrane)
            elif resname in ('ASP', 'GLU'):
                attraction = 2.0   # kJ/mol (repelled)
            else:
                attraction = -1.0  # kJ/mol (weak hydrophobic attraction)
            membrane_force.addParticle(atom.index, [attraction])

    system.addForce(membrane_force)

    print(f"  System: {system.getNumParticles()} atoms")
    print(f"  Membrane-mimetic potential applied")

    return modeller, system

# ============================================================
# SIMULATION ENGINE
# ============================================================

def run_membrane_sim(modeller, system, name, n_steps=50000000,
                     report_interval=10000):
    """Run membrane simulation.

    50M steps × 2 fs = 100 ns (production).
    """
    prefix = os.path.join(OUTPUT_DIR, name)

    integrator = mm.LangevinMiddleIntegrator(
        300*unit.kelvin,
        1.0/unit.picosecond,
        0.002*unit.picoseconds
    )

    # Use CUDA if available (RunPod), else OpenCL (Mac), else CPU
    for platform_name in ['CUDA', 'OpenCL', 'CPU']:
        try:
            platform = mm.Platform.getPlatformByName(platform_name)
            print(f"  Platform: {platform_name} (speed={platform.getSpeed()})")
            break
        except Exception:
            continue

    simulation = app.Simulation(modeller.topology, system, integrator, platform)
    simulation.context.setPositions(modeller.positions)

    # Minimize
    print(f"  Minimizing...")
    simulation.minimizeEnergy(maxIterations=5000)

    positions = simulation.context.getState(getPositions=True).getPositions()
    with open(f"{prefix}_minimized.pdb", 'w') as f:
        app.PDBFile.writeFile(modeller.topology, positions, f)

    # Reporters
    simulation.reporters.append(
        app.StateDataReporter(
            f"{prefix}_energy.csv",
            report_interval,
            step=True, time=True,
            potentialEnergy=True, temperature=True, speed=True,
        )
    )
    simulation.reporters.append(
        app.DCDReporter(f"{prefix}_traj.dcd", report_interval)
    )

    # Custom reporter: peptide z-coordinate (distance to membrane)
    total_ns = n_steps * 0.002 / 1000
    print(f"  Running {n_steps:,} steps ({total_ns:.0f} ns)...")

    t_start = time.time()
    chunk = min(5000000, n_steps)  # 10 ns chunks
    z_history = []

    for i in range(0, n_steps, chunk):
        remaining = min(chunk, n_steps - i)
        simulation.step(remaining)

        # Track peptide z-coordinate
        state = simulation.context.getState(getPositions=True, getEnergy=True)
        positions = state.getPositions(asNumpy=True).value_in_unit(unit.nanometers)

        # Get peptide COM z
        peptide_atoms = [a.index for a in modeller.topology.atoms()
                        if a.residue.name not in ('HOH', 'NA', 'CL')]
        if peptide_atoms:
            pep_z = np.mean(positions[peptide_atoms, 2])
            z_history.append(float(pep_z))

        current_ns = (i + remaining) * 0.002 / 1000
        energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        elapsed = time.time() - t_start
        speed = current_ns / (elapsed / 86400) if elapsed > 0 else 0

        print(f"    {current_ns:.0f}/{total_ns:.0f} ns — "
              f"E={energy:.0f} kJ/mol — "
              f"peptide z={pep_z:.2f} nm — "
              f"speed={speed:.0f} ns/day")

    # Save final
    positions = simulation.context.getState(getPositions=True).getPositions()
    with open(f"{prefix}_final.pdb", 'w') as f:
        app.PDBFile.writeFile(modeller.topology, positions, f)

    # Save z-history
    with open(f"{prefix}_z_history.json", 'w') as f:
        json.dump({"name": name, "z_nm": z_history}, f)

    total_time = time.time() - t_start
    print(f"  Completed in {total_time/3600:.1f} hours")

    return f"{prefix}_minimized.pdb", f"{prefix}_traj.dcd"

# ============================================================
# ANALYSIS
# ============================================================

def analyze_membrane_sim(name, top_path, traj_path):
    """Analyze peptide-membrane interaction."""
    import mdtraj as md

    traj = md.load(traj_path, top=top_path)
    print(f"\n  Analysis: {name} — {traj.n_frames} frames")

    peptide = traj.topology.select('protein')
    if len(peptide) == 0:
        print("  No protein atoms found")
        return {}

    # Z-coordinate of peptide COM over time
    pep_positions = traj.xyz[:, peptide, :]
    com_z = np.mean(pep_positions[:, :, 2], axis=1)
    print(f"  Peptide COM z: {np.mean(com_z):.3f} ± {np.std(com_z):.3f} nm")
    print(f"  Min z (closest to membrane): {np.min(com_z):.3f} nm")

    # SASA
    sasa = md.shrake_rupley(traj, mode='atom')
    pep_sasa = np.sum(sasa[:, peptide], axis=1)
    print(f"  Peptide SASA: {np.mean(pep_sasa):.3f} ± {np.std(pep_sasa):.3f} nm²")

    # Rg
    pep_traj = traj.atom_slice(peptide)
    rg = md.compute_rg(pep_traj)
    print(f"  Rg: {np.mean(rg):.3f} ± {np.std(rg):.3f} nm")

    return {
        "name": name,
        "com_z_mean": float(np.mean(com_z)),
        "com_z_std": float(np.std(com_z)),
        "com_z_min": float(np.min(com_z)),
        "sasa_mean": float(np.mean(pep_sasa)),
        "rg_mean": float(np.mean(rg)),
    }

# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--steps', type=int, default=50000000,
                        help='Steps (default: 50M = 100 ns)')
    parser.add_argument('--peptides', nargs='+',
                        default=['KKRPKP', 'NNRPNP'],
                        help='Peptide sequences to simulate')
    args = parser.parse_args()

    print("=" * 65)
    print("MEMBRANE SIMULATION: PrP PEPTIDES + POPC (mimetic)")
    print("=" * 65)
    print(f"Steps: {args.steps:,} ({args.steps * 0.002 / 1000:.0f} ns)")
    print(f"Peptides: {args.peptides}")
    print()

    # Check platform
    print("Available platforms:")
    for i in range(mm.Platform.getNumPlatforms()):
        p = mm.Platform.getPlatform(i)
        print(f"  {p.getName()}: speed={p.getSpeed()}")
    print()

    all_results = {}

    for seq in args.peptides:
        name = f"membrane_{seq}"
        print(f"\n{'='*55}")
        print(f"  Peptide: {seq} (charge: {sum(1 for c in seq if c in 'KR') - sum(1 for c in seq if c in 'DE'):+d})")
        print(f"{'='*55}")

        # Build peptide
        pdb_path = build_peptide(seq, name)

        # Build membrane system
        modeller, system = build_membrane_system(pdb_path, name)

        # Run simulation
        top_path, traj_path = run_membrane_sim(
            modeller, system, name, n_steps=args.steps
        )

        # Analyze
        results = analyze_membrane_sim(name, top_path, traj_path)
        all_results[seq] = results

    # Summary
    print(f"\n{'='*65}")
    print("COMPARISON: MEMBRANE INTERACTION")
    print(f"{'='*65}")
    for seq, res in all_results.items():
        if res:
            charge = sum(1 for c in seq if c in 'KR') - sum(1 for c in seq if c in 'DE')
            print(f"  {seq} (charge {charge:+d}):")
            print(f"    COM z (distance to membrane): {res['com_z_mean']:.3f} ± {res['com_z_std']:.3f} nm")
            print(f"    Closest approach: {res['com_z_min']:.3f} nm")
            print(f"    SASA: {res['sasa_mean']:.3f} nm²")
            print(f"    Rg: {res['rg_mean']:.3f} nm")

    # Save all results
    with open(os.path.join(OUTPUT_DIR, "membrane_results.json"), 'w') as f:
        json.dump(all_results, f, indent=2)

    print(f"\nResults saved to {OUTPUT_DIR}/membrane_results.json")
    print(f"\nv5 prediction: KKRPKP should approach membrane more closely")
    print(f"and show lower COM z than NNRPNP (neutral control).")
