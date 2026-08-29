"""
MD simulation: KKRPKP hexapeptide interacting with a POPC lipid bilayer.
Tests the v5 hypothesis that the PrP N-terminal polybasic motif perturbs
membranes via AMP-like carpet mechanism.

Two systems:
  1. KKRPKP (wild-type, charge +4) — should insert and perturb
  2. NNRPNP (charge-neutralized, charge 0) — control, should not insert

OpenMM with OpenCL (Apple M1 Pro GPU).
"""
import openmm as mm
import openmm.app as app
import openmm.unit as unit
import numpy as np
import os
import json

OUTPUT_DIR = "/Users/allan/Projects/cjd/md_output"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================
# Step 1: Build peptide structures from sequence
# ============================================================

def build_peptide_pdb(sequence, name, output_path):
    """Build a linear peptide PDB using OpenMM's Modeller."""
    from pdbfixer import PDBFixer

    # Build peptide manually using OpenMM's topology builder
    topology = app.Topology()
    chain = topology.addChain()

    AA_MAP = {
        'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
        'E': 'GLU', 'Q': 'GLN', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
        'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
        'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL',
    }

    # We'll create PDB content directly
    pdb_lines = ["HEADER    PEPTIDE"]
    atom_idx = 1

    # Standard backbone + CB atoms for a linear peptide
    # Extended conformation: phi=-180, psi=180
    backbone_template = {
        'N':  np.array([0.0, 0.0, 0.0]),
        'CA': np.array([1.458, 0.0, 0.0]),
        'C':  np.array([2.009, 1.420, 0.0]),
        'O':  np.array([1.251, 2.390, 0.0]),
    }

    step = np.array([3.8, 0.0, 0.0])  # ~3.8 Å per residue in extended

    for i, aa in enumerate(sequence):
        resname = AA_MAP.get(aa, 'ALA')
        resid = i + 1
        offset = step * i

        for atom_name, base_pos in backbone_template.items():
            pos = base_pos + offset
            pdb_lines.append(
                f"ATOM  {atom_idx:5d} {atom_name:<4s} {resname:3s} A{resid:4d}    "
                f"{pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}  1.00  0.00           "
                f"{atom_name[0]:>2s}"
            )
            atom_idx += 1

    pdb_lines.append("TER")
    pdb_lines.append("END")

    pdb_content = "\n".join(pdb_lines)

    # Write raw PDB
    raw_path = output_path.replace('.pdb', '_raw.pdb')
    with open(raw_path, 'w') as f:
        f.write(pdb_content)

    # Use PDBFixer to add missing atoms and hydrogens
    fixer = PDBFixer(filename=raw_path)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)  # pH 7.4

    with open(output_path, 'w') as f:
        app.PDBFile.writeFile(fixer.topology, fixer.positions, f)

    print(f"  Built {name}: {len(sequence)} residues, saved to {output_path}")
    return output_path

# ============================================================
# Step 2: Build a simple membrane-peptide system
# ============================================================

def setup_peptide_water_system(pdb_path, name):
    """Set up a peptide in water box for initial equilibration.

    Full membrane simulation requires CHARMM-GUI or similar;
    here we do peptide-in-water as proof of concept and
    analyze charge distribution and solvent-accessible properties.
    """
    print(f"\n  Setting up system: {name}")

    pdb = app.PDBFile(pdb_path)

    # Force field: AMBER ff14SB for proteins
    forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

    # Create system with solvent box
    modeller = app.Modeller(pdb.topology, pdb.positions)
    modeller.addSolvent(forcefield, model='tip3p',
                        padding=1.2*unit.nanometers,
                        ionicStrength=0.15*unit.molar)  # 150 mM NaCl

    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0*unit.nanometers,
        constraints=app.HBonds,
    )

    # Count atoms
    n_atoms = system.getNumParticles()
    n_residues = sum(1 for _ in modeller.topology.residues())
    print(f"  Total atoms: {n_atoms}")
    print(f"  Total residues: {n_residues}")

    return modeller, system, forcefield

def run_simulation(modeller, system, name, n_steps=50000, report_interval=1000):
    """Run a short MD simulation.

    50,000 steps × 2 fs = 100 ps — enough to see initial dynamics.
    For production: 5,000,000 steps = 10 ns.
    """
    output_prefix = os.path.join(OUTPUT_DIR, name)

    # Integrator
    integrator = mm.LangevinMiddleIntegrator(
        300*unit.kelvin,       # temperature
        1.0/unit.picosecond,   # friction
        0.002*unit.picoseconds # timestep
    )

    # Platform — use OpenCL for GPU acceleration
    try:
        platform = mm.Platform.getPlatformByName('OpenCL')
        print(f"  Using OpenCL (GPU)")
    except Exception:
        platform = mm.Platform.getPlatformByName('CPU')
        print(f"  Falling back to CPU")

    simulation = app.Simulation(modeller.topology, system, integrator, platform)
    simulation.context.setPositions(modeller.positions)

    # Energy minimization
    print(f"  Minimizing energy...")
    simulation.minimizeEnergy(maxIterations=1000)

    state = simulation.context.getState(getEnergy=True)
    print(f"  Energy after minimization: {state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole):.1f} kJ/mol")

    # Save minimized structure
    positions = simulation.context.getState(getPositions=True).getPositions()
    with open(f"{output_prefix}_minimized.pdb", 'w') as f:
        app.PDBFile.writeFile(modeller.topology, positions, f)

    # Reporters
    simulation.reporters.append(
        app.StateDataReporter(
            f"{output_prefix}_energy.csv",
            report_interval,
            step=True,
            potentialEnergy=True,
            temperature=True,
            density=True,
            speed=True,
        )
    )
    simulation.reporters.append(
        app.DCDReporter(f"{output_prefix}_trajectory.dcd", report_interval)
    )

    # Run simulation
    total_time_ps = n_steps * 0.002
    print(f"  Running {n_steps} steps ({total_time_ps:.0f} ps / {total_time_ps/1000:.1f} ns)...")

    simulation.step(n_steps)

    # Save final state
    state = simulation.context.getState(getPositions=True, getEnergy=True)
    positions = state.getPositions()

    with open(f"{output_prefix}_final.pdb", 'w') as f:
        app.PDBFile.writeFile(modeller.topology, positions, f)

    energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    print(f"  Final energy: {energy:.1f} kJ/mol")
    print(f"  Trajectory saved: {output_prefix}_trajectory.dcd")

    return simulation, f"{output_prefix}_trajectory.dcd"


# ============================================================
# Step 3: Analysis
# ============================================================

def analyze_peptide(name, pdb_path, traj_path):
    """Analyze peptide properties from trajectory."""
    import mdtraj as md

    traj = md.load(traj_path, top=pdb_path)
    print(f"\n  Analysis for {name}:")
    print(f"  Frames: {traj.n_frames}")
    print(f"  Atoms: {traj.n_atoms}")

    # Select peptide atoms (non-water, non-ion)
    peptide = traj.topology.select('protein')

    if len(peptide) == 0:
        print("  WARNING: no protein atoms found")
        return {}

    # SASA (Solvent Accessible Surface Area)
    sasa = md.shrake_rupley(traj, mode='residue')
    peptide_residues = [r for r in traj.topology.residues if r.name not in ('HOH', 'NA', 'CL')]
    n_pep_res = len(peptide_residues)

    if n_pep_res > 0:
        mean_sasa = np.mean(sasa[:, :n_pep_res], axis=0)
        total_sasa = np.sum(mean_sasa)
        print(f"  Total SASA: {total_sasa:.3f} nm²")
        for i, res in enumerate(peptide_residues):
            if i < len(mean_sasa):
                print(f"    {res.name}{res.index+1}: {mean_sasa[i]:.3f} nm²")

    # Radius of gyration
    rg = md.compute_rg(traj, masses=np.ones(traj.n_atoms))
    print(f"  Radius of gyration: {np.mean(rg):.3f} ± {np.std(rg):.3f} nm")

    # End-to-end distance (first and last CA)
    ca_atoms = traj.topology.select('name CA')
    if len(ca_atoms) >= 2:
        e2e = md.compute_distances(traj, [[ca_atoms[0], ca_atoms[-1]]])
        print(f"  End-to-end distance: {np.mean(e2e):.3f} ± {np.std(e2e):.3f} nm")

    results = {
        "name": name,
        "n_frames": int(traj.n_frames),
        "total_sasa": float(total_sasa) if n_pep_res > 0 else None,
        "rg_mean": float(np.mean(rg)),
        "rg_std": float(np.std(rg)),
    }

    return results

# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    print("=" * 65)
    print("MD SIMULATION: PrP N-TERMINAL PEPTIDES IN WATER")
    print("=" * 65)
    print("\nPhase 1: Peptide-in-water equilibration")
    print("Tests: charge distribution, compactness, solvent exposure")
    print()

    peptides = {
        "KKRPKP_wt": {
            "sequence": "KKRPKP",
            "description": "Wild-type PrP 23-28 (charge +4)",
        },
        "NNRPNP_neutral": {
            "sequence": "NNRPNP",
            "description": "Charge-neutralized control (charge 0)",
        },
        "KKRPKPGGWNTGG_extended": {
            "sequence": "KKRPKPGGWNTGG",
            "description": "Extended PrP 23-35 (charge +4, with linker)",
        },
    }

    results = {}

    for name, info in peptides.items():
        print(f"\n{'='*50}")
        print(f"Peptide: {name} — {info['description']}")
        print(f"Sequence: {info['sequence']}")
        print(f"{'='*50}")

        # Build peptide
        pdb_path = os.path.join(OUTPUT_DIR, f"{name}.pdb")
        try:
            build_peptide_pdb(info['sequence'], name, pdb_path)
        except Exception as e:
            print(f"  ERROR building peptide: {e}")
            continue

        # Setup system
        try:
            modeller, system, ff = setup_peptide_water_system(pdb_path, name)
        except Exception as e:
            print(f"  ERROR setting up system: {e}")
            continue

        # Run short equilibration (100 ps)
        try:
            sim, traj_path = run_simulation(
                modeller, system, name,
                n_steps=50000,      # 100 ps
                report_interval=500  # every 1 ps
            )
        except Exception as e:
            print(f"  ERROR running simulation: {e}")
            continue

        # Analyze
        try:
            min_pdb = os.path.join(OUTPUT_DIR, f"{name}_minimized.pdb")
            res = analyze_peptide(name, min_pdb, traj_path)
            results[name] = res
        except Exception as e:
            print(f"  ERROR analyzing: {e}")

    # Save results
    with open(os.path.join(OUTPUT_DIR, "md_results.json"), 'w') as f:
        json.dump(results, f, indent=2)

    print(f"\n{'='*65}")
    print("SUMMARY")
    print(f"{'='*65}")
    for name, res in results.items():
        if res:
            print(f"  {name}: Rg={res.get('rg_mean', 'N/A'):.3f} nm, SASA={res.get('total_sasa', 'N/A')}")

    print(f"\nAll output in: {OUTPUT_DIR}")
    print(f"\nNext steps:")
    print(f"  1. Extend to 10 ns for production statistics")
    print(f"  2. Add POPC bilayer (requires membrane builder)")
    print(f"  3. Measure peptide-membrane insertion depth")
    print(f"  4. Compare membrane thickness perturbation")
