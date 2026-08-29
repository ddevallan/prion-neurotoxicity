"""
Production MD script for CHARMM-GUI membrane systems on Vast.ai/RunPod.
Runs after CHARMM-GUI equilibration is complete.

Usage:
  python run_production.py --system system_KKRPKP --ns 300 --replicas 3

Requires: CHARMM-GUI generated files in the system directory.
"""
import openmm as mm
import openmm.app as app
import openmm.unit as unit
import numpy as np
import os
import sys
import time
import json
import argparse
import glob

def find_charmm_gui_files(system_dir):
    """Locate CHARMM-GUI output files."""
    openmm_dir = os.path.join(system_dir, "openmm")
    if not os.path.isdir(openmm_dir):
        openmm_dir = system_dir

    # Find the last equilibration step
    equil_files = sorted(glob.glob(os.path.join(openmm_dir, "step6*")))

    # Find topology and coordinate files
    psf = glob.glob(os.path.join(system_dir, "*.psf"))
    pdb = glob.glob(os.path.join(openmm_dir, "step5_input.pdb"))
    if not pdb:
        pdb = glob.glob(os.path.join(system_dir, "step5_input.pdb"))

    # Find force field parameters
    toppar = os.path.join(system_dir, "toppar")
    if not os.path.isdir(toppar):
        toppar = os.path.join(system_dir, "toppar_c36m")

    return {
        "openmm_dir": openmm_dir,
        "psf": psf[0] if psf else None,
        "pdb": pdb[0] if pdb else None,
        "toppar": toppar if os.path.isdir(toppar) else None,
        "equil_files": equil_files,
    }

def setup_from_charmm_gui(system_dir):
    """
    Load system from CHARMM-GUI equilibrated files.
    Expects the standard CHARMM-GUI output structure.
    """
    files = find_charmm_gui_files(system_dir)
    print(f"  System dir: {system_dir}")
    print(f"  PSF: {files['psf']}")
    print(f"  PDB: {files['pdb']}")
    print(f"  Toppar: {files['toppar']}")

    if not files['psf'] or not files['pdb']:
        print("ERROR: Missing PSF or PDB files from CHARMM-GUI")
        print("Expected structure:")
        print("  system_dir/")
        print("    ├── *.psf")
        print("    ├── step5_input.pdb")
        print("    ├── toppar/")
        print("    └── openmm/")
        sys.exit(1)

    # Load PSF and PDB
    psf = app.CharmmPsfFile(files['psf'])
    pdb = app.PDBFile(files['pdb'])

    # Load CHARMM parameters
    params = app.CharmmParameterSet(
        *glob.glob(os.path.join(files['toppar'], "*.prm")),
        *glob.glob(os.path.join(files['toppar'], "*.rtf")),
        *glob.glob(os.path.join(files['toppar'], "*.str")),
    )

    # Create system
    system = psf.createSystem(
        params,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.2*unit.nanometers,
        switchDistance=1.0*unit.nanometers,
        constraints=app.HBonds,
    )

    # Add barostat (semi-isotropic for membrane)
    barostat = mm.MonteCarloMembraneBarostat(
        1.0*unit.bar,
        0.0*unit.bar*unit.nanometers,
        300*unit.kelvin,
        mm.MonteCarloMembraneBarostat.XYIsotropic,
        mm.MonteCarloMembraneBarostat.ZFree,
    )
    system.addForce(barostat)

    return psf, pdb, system

def setup_from_amber(pdb_path):
    """
    Alternative: set up from PDB + AMBER forcefield (simpler, no CHARMM-GUI needed).
    Uses AMBER14/lipid17 which supports POPC.
    """
    pdb = app.PDBFile(pdb_path)
    forcefield = app.ForceField(
        'amber14-all.xml',
        'amber14/tip3pfb.xml',
        # Note: lipid17 for POPC would need separate parameter loading
    )

    modeller = app.Modeller(pdb.topology, pdb.positions)
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0*unit.nanometers,
        constraints=app.HBonds,
    )

    return modeller, system

def run_production(system, topology, positions, name, n_ns=300,
                   output_dir="./output", checkpoint_interval_ns=10):
    """
    Run production MD with checkpointing.

    Saves:
    - Trajectory (DCD, every 10 ps)
    - Energy (CSV, every 1 ps)
    - Checkpoints (every checkpoint_interval_ns)
    """
    os.makedirs(output_dir, exist_ok=True)
    prefix = os.path.join(output_dir, name)

    n_steps = int(n_ns * 1e6 / 2)  # 2 fs timestep
    report_interval = 5000          # 10 ps
    energy_interval = 500           # 1 ps
    checkpoint_steps = int(checkpoint_interval_ns * 1e6 / 2)

    # Integrator
    integrator = mm.LangevinMiddleIntegrator(
        310*unit.kelvin,          # physiological temperature
        1.0/unit.picosecond,
        0.002*unit.picoseconds
    )

    # Platform selection
    for pname in ['CUDA', 'OpenCL', 'CPU']:
        try:
            platform = mm.Platform.getPlatformByName(pname)
            properties = {}
            if pname == 'CUDA':
                properties = {'Precision': 'mixed'}
            print(f"  Platform: {pname}")
            break
        except Exception:
            continue

    simulation = app.Simulation(topology, system, integrator, platform,
                                properties if pname == 'CUDA' else {})
    simulation.context.setPositions(positions)

    # Check for checkpoint to resume from
    chk_path = f"{prefix}_checkpoint.chk"
    if os.path.exists(chk_path):
        simulation.loadCheckpoint(chk_path)
        current_step = simulation.context.getState().getStepCount()
        print(f"  Resumed from checkpoint at step {current_step}")
    else:
        simulation.context.setVelocitiesToTemperature(310*unit.kelvin)

    # Reporters
    simulation.reporters.append(
        app.DCDReporter(f"{prefix}_traj.dcd", report_interval)
    )
    simulation.reporters.append(
        app.StateDataReporter(
            f"{prefix}_energy.csv",
            energy_interval,
            step=True, time=True,
            potentialEnergy=True, kineticEnergy=True,
            temperature=True, volume=True,
            speed=True, separator=','
        )
    )
    simulation.reporters.append(
        app.CheckpointReporter(chk_path, checkpoint_steps)
    )

    # Run
    print(f"  Production: {n_ns} ns ({n_steps:,} steps)")
    t_start = time.time()

    chunk_steps = 500000  # 1 ns
    steps_done = 0

    while steps_done < n_steps:
        remaining = min(chunk_steps, n_steps - steps_done)
        simulation.step(remaining)
        steps_done += remaining

        ns_done = steps_done * 0.002 / 1000
        elapsed = time.time() - t_start
        speed = ns_done / (elapsed / 86400) if elapsed > 0 else 0

        state = simulation.context.getState(getEnergy=True)
        temp = state.getKineticEnergy() * 2 / (
            3 * system.getNumParticles() * unit.MOLAR_GAS_CONSTANT_R
        )

        if steps_done % (chunk_steps * 10) == 0 or steps_done >= n_steps:
            print(f"    {ns_done:.0f}/{n_ns} ns — "
                  f"speed={speed:.0f} ns/day — "
                  f"E={state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole):.0f} kJ/mol")

    # Save final
    state = simulation.context.getState(getPositions=True)
    with open(f"{prefix}_final.pdb", 'w') as f:
        app.PDBFile.writeFile(topology, state.getPositions(), f)

    total_h = (time.time() - t_start) / 3600
    print(f"  Completed {name}: {n_ns} ns in {total_h:.1f} hours")

    return {
        "name": name,
        "ns": n_ns,
        "hours": round(total_h, 2),
        "trajectory": f"{prefix}_traj.dcd",
        "topology": f"{prefix}_final.pdb",
    }


def main():
    parser = argparse.ArgumentParser(description="Membrane MD production run")
    parser.add_argument('--system', required=True,
                        help='Path to CHARMM-GUI system directory')
    parser.add_argument('--ns', type=int, default=300,
                        help='Production time in ns (default: 300)')
    parser.add_argument('--replicas', type=int, default=3,
                        help='Number of replicas (default: 3)')
    parser.add_argument('--output', default='./output',
                        help='Output directory')
    parser.add_argument('--name', default=None,
                        help='System name (default: from directory)')
    args = parser.parse_args()

    name = args.name or os.path.basename(args.system.rstrip('/'))
    print("=" * 65)
    print(f"MEMBRANE MD PRODUCTION: {name}")
    print(f"  {args.ns} ns × {args.replicas} replicas = {args.ns * args.replicas} ns total")
    print("=" * 65)

    # Load system
    psf, pdb, system = setup_from_charmm_gui(args.system)

    # Run replicas
    results = []
    for rep in range(1, args.replicas + 1):
        rep_name = f"{name}_rep{rep}"
        print(f"\n--- Replica {rep}/{args.replicas} ---")
        res = run_production(
            system, psf.topology, pdb.positions,
            rep_name, n_ns=args.ns, output_dir=args.output
        )
        results.append(res)

    # Save summary
    with open(os.path.join(args.output, f"{name}_summary.json"), 'w') as f:
        json.dump(results, f, indent=2)

    total_hours = sum(r['hours'] for r in results)
    print(f"\n{'='*65}")
    print(f"TOTAL: {total_hours:.1f} hours for {args.ns * args.replicas} ns")
    print(f"{'='*65}")


if __name__ == "__main__":
    main()
