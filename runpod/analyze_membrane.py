"""
Post-simulation analysis of peptide-membrane interactions.
Run locally after downloading trajectories from Vast.ai.

Measures the key v5 predictions:
1. Peptide insertion depth
2. Membrane thinning
3. Lipid order perturbation
4. Contact analysis
"""
import mdtraj as md
import numpy as np
import json
import os
import sys
import argparse

def load_trajectory(top_path, traj_path, stride=10):
    """Load trajectory with optional stride for memory efficiency."""
    print(f"  Loading: {traj_path}")
    traj = md.load(traj_path, top=top_path, stride=stride)
    print(f"  Frames: {traj.n_frames}, Atoms: {traj.n_atoms}")
    return traj

def peptide_insertion_depth(traj):
    """
    Calculate peptide COM z-position relative to membrane center.
    Negative z = inserted into membrane.
    """
    peptide = traj.topology.select('protein')
    if len(peptide) == 0:
        return None

    # Peptide COM z
    pep_z = np.mean(traj.xyz[:, peptide, 2], axis=1)

    # Phosphorus atoms mark membrane surface (P atoms in POPC headgroups)
    phosphorus = traj.topology.select('name P')
    if len(phosphorus) == 0:
        # No membrane — use box center as reference
        box_z = traj.unitcell_lengths[:, 2] / 2 if traj.unitcell_lengths is not None else np.zeros(traj.n_frames)
        return {
            "pep_z_mean": float(np.mean(pep_z)),
            "pep_z_std": float(np.std(pep_z)),
            "pep_z_min": float(np.min(pep_z)),
            "reference": "box_center",
            "note": "No membrane detected (water-only simulation)"
        }

    # Upper and lower leaflet P atoms
    p_z = traj.xyz[:, phosphorus, 2]
    membrane_center = np.mean(p_z, axis=1)
    upper_p = np.max(p_z, axis=1)
    lower_p = np.min(p_z, axis=1)

    # Relative to membrane center
    relative_z = pep_z - membrane_center

    # Relative to upper leaflet surface
    depth = pep_z - upper_p  # negative = inserted below surface

    return {
        "pep_z_mean": float(np.mean(relative_z)),
        "pep_z_std": float(np.std(relative_z)),
        "pep_z_min": float(np.min(relative_z)),
        "depth_mean": float(np.mean(depth)),
        "depth_min": float(np.min(depth)),
        "membrane_thickness": float(np.mean(upper_p - lower_p)),
        "reference": "membrane_center",
        "timeseries": {
            "z": [float(x) for x in relative_z[::10]],
            "depth": [float(x) for x in depth[::10]],
        }
    }

def membrane_thickness_map(traj, grid_points=20):
    """
    Calculate local membrane thickness as a 2D map.
    Reveals thinning near the peptide.
    """
    phosphorus = traj.topology.select('name P')
    if len(phosphorus) < 4:
        return None

    # Use last 50% of trajectory (equilibrated)
    n_start = traj.n_frames // 2
    p_pos = traj.xyz[n_start:, phosphorus, :]

    # Determine upper vs lower leaflet
    mean_z = np.mean(p_pos[:, :, 2])
    upper_mask = np.mean(p_pos[:, :, 2], axis=0) > mean_z
    lower_mask = ~upper_mask

    if not np.any(upper_mask) or not np.any(lower_mask):
        return None

    # Grid
    box_x = np.mean(traj.unitcell_lengths[n_start:, 0])
    box_y = np.mean(traj.unitcell_lengths[n_start:, 1])

    x_edges = np.linspace(0, box_x, grid_points + 1)
    y_edges = np.linspace(0, box_y, grid_points + 1)

    thickness_map = np.zeros((grid_points, grid_points))
    counts = np.zeros((grid_points, grid_points))

    for frame in range(len(p_pos)):
        for i in range(len(phosphorus)):
            x, y, z = p_pos[frame, i]
            xi = min(int(x / box_x * grid_points), grid_points - 1)
            yi = min(int(y / box_y * grid_points), grid_points - 1)

            if upper_mask[i]:
                thickness_map[xi, yi] += z
                counts[xi, yi] += 1
            else:
                thickness_map[xi, yi] -= z
                counts[xi, yi] += 1

    # Normalize
    valid = counts > 0
    thickness_map[valid] /= counts[valid]

    return {
        "thickness_map": thickness_map.tolist(),
        "mean_thickness": float(np.mean(thickness_map[valid])),
        "min_thickness": float(np.min(thickness_map[valid])),
        "max_thickness": float(np.max(thickness_map[valid])),
    }

def contact_analysis(traj, cutoff_nm=0.4):
    """Count peptide-lipid contacts over time."""
    peptide = traj.topology.select('protein')
    lipid = traj.topology.select('resname POPC')

    if len(peptide) == 0 or len(lipid) == 0:
        return {"note": "No peptide-lipid contacts (water-only simulation)"}

    headgroup = traj.topology.select('resname POPC and (name P or name N or name O*)')
    acyl = traj.topology.select('resname POPC and name C2* C3*')

    # Count contacts
    n_frames = traj.n_frames
    hg_contacts = np.zeros(n_frames)
    acyl_contacts = np.zeros(n_frames)

    for i in range(n_frames):
        for pa in peptide:
            for ha in headgroup:
                d = np.linalg.norm(traj.xyz[i, pa] - traj.xyz[i, ha])
                if d < cutoff_nm:
                    hg_contacts[i] += 1
            for aa in acyl[:50]:
                d = np.linalg.norm(traj.xyz[i, pa] - traj.xyz[i, aa])
                if d < cutoff_nm:
                    acyl_contacts[i] += 1

    return {
        "headgroup_contacts_mean": float(np.mean(hg_contacts)),
        "acyl_contacts_mean": float(np.mean(acyl_contacts)),
        "headgroup_contacts_max": float(np.max(hg_contacts)),
    }

def sasa_analysis(traj):
    """SASA of peptide over time."""
    peptide = traj.topology.select('protein')
    peptide_residues = [r for r in traj.topology.residues
                       if r.name not in ('HOH', 'NA', 'CL', 'K', 'POPC')]
    n_res = len(peptide_residues)

    if n_res == 0:
        return {}

    sasa = md.shrake_rupley(traj, mode='residue')
    total_sasa = np.sum(sasa[:, :n_res], axis=1)

    per_residue = {}
    mean_sasa = np.mean(sasa[:, :n_res], axis=0)
    for i, res in enumerate(peptide_residues):
        if i < len(mean_sasa):
            per_residue[f"{res.name}{res.index+1}"] = float(mean_sasa[i])

    return {
        "total_sasa_mean": float(np.mean(total_sasa)),
        "total_sasa_std": float(np.std(total_sasa)),
        "per_residue": per_residue,
        "timeseries": [float(x) for x in total_sasa[::10]],
    }

def full_analysis(top_path, traj_path, name, output_dir="./analysis"):
    """Run all analyses."""
    os.makedirs(output_dir, exist_ok=True)

    print(f"\n{'='*55}")
    print(f"  ANALYSIS: {name}")
    print(f"{'='*55}")

    traj = load_trajectory(top_path, traj_path, stride=5)

    results = {"name": name}

    print("\n  1. Insertion depth...")
    results["insertion"] = peptide_insertion_depth(traj)
    if results["insertion"]:
        print(f"     z = {results['insertion']['pep_z_mean']:.3f} ± {results['insertion']['pep_z_std']:.3f} nm")

    print("  2. SASA...")
    results["sasa"] = sasa_analysis(traj)
    if results["sasa"]:
        print(f"     SASA = {results['sasa']['total_sasa_mean']:.3f} ± {results['sasa']['total_sasa_std']:.3f} nm²")

    print("  3. Rg...")
    peptide = traj.topology.select('protein')
    if len(peptide) > 0:
        pep_traj = traj.atom_slice(peptide)
        rg = md.compute_rg(pep_traj)
        results["rg"] = {
            "mean": float(np.mean(rg)),
            "std": float(np.std(rg)),
        }
        print(f"     Rg = {np.mean(rg):.3f} ± {np.std(rg):.3f} nm")

    print("  4. Secondary structure...")
    if len(peptide) > 0:
        dssp = md.compute_dssp(traj)
        peptide_residues = [r for r in traj.topology.residues
                           if r.name not in ('HOH', 'NA', 'CL', 'K', 'POPC')]
        n_res = len(peptide_residues)
        ss_fractions = {}
        for ss_type in ['H', 'E', 'C']:
            frac = np.mean(dssp[:, :n_res] == ss_type)
            ss_fractions[{'H': 'helix', 'E': 'sheet', 'C': 'coil'}[ss_type]] = float(frac)
        results["secondary_structure"] = ss_fractions
        print(f"     Helix: {ss_fractions.get('helix', 0)*100:.1f}%, "
              f"Sheet: {ss_fractions.get('sheet', 0)*100:.1f}%, "
              f"Coil: {ss_fractions.get('coil', 0)*100:.1f}%")

    # Save
    out_path = os.path.join(output_dir, f"{name}_analysis.json")
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\n  Saved to {out_path}")

    return results

def compare_systems(results_list, output_dir="./analysis"):
    """Compare multiple systems side by side."""
    print(f"\n{'='*65}")
    print("COMPARISON")
    print(f"{'='*65}")

    print(f"\n{'System':<25} {'SASA':>10} {'Rg':>10} {'z-depth':>10}")
    print("-" * 55)

    for r in results_list:
        name = r.get('name', '?')
        sasa = r.get('sasa', {}).get('total_sasa_mean', 0)
        rg = r.get('rg', {}).get('mean', 0)
        z = r.get('insertion', {}).get('pep_z_mean', 0)
        print(f"{name:<25} {sasa:>10.3f} {rg:>10.3f} {z:>10.3f}")

    # Save comparison
    out_path = os.path.join(output_dir, "comparison.json")
    with open(out_path, 'w') as f:
        json.dump(results_list, f, indent=2, default=str)
    print(f"\nComparison saved to {out_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--top', required=True, nargs='+',
                        help='Topology file(s) (PDB)')
    parser.add_argument('--traj', required=True, nargs='+',
                        help='Trajectory file(s) (DCD)')
    parser.add_argument('--names', nargs='+', default=None,
                        help='System names')
    parser.add_argument('--output', default='./analysis')
    args = parser.parse_args()

    if len(args.top) != len(args.traj):
        print("ERROR: must provide same number of --top and --traj files")
        sys.exit(1)

    names = args.names or [f"system_{i}" for i in range(len(args.top))]
    all_results = []

    for top, traj, name in zip(args.top, args.traj, names):
        res = full_analysis(top, traj, name, args.output)
        all_results.append(res)

    if len(all_results) > 1:
        compare_systems(all_results, args.output)
