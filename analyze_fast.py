"""Fast parallel analysis — skip SASA (too slow for 40k atoms × 1250 frames)."""
import mdtraj as md
import numpy as np
import json
import sys
import os
from multiprocessing import Pool

RESULTS_DIR = "/Users/allan/Projects/cjd/results_mimetic"

SYSTEMS = [
    ("KKRPKP_membrane", "KKRPKP (+4, wild-type)", 4, "KKRPKP"),
    ("NNRPNP_membrane", "NNRPNP (0, neutral)", 0, "NNRPNP"),
    ("PrP_23_35_membrane", "PrP 23-35 (+4, Trp)", 4, "KKRPKPGGWNTGG"),
]

def analyze(args):
    name, label, charge, sequence = args
    traj_path = os.path.join(RESULTS_DIR, f"{name}_traj.dcd")
    top_path = os.path.join(RESULTS_DIR, f"{name}_system.pdb")

    traj = md.load(traj_path, top=top_path)
    peptide = traj.topology.select('protein')
    pep_res = [r for r in traj.topology.residues if r.name not in ('HOH','NA','CL')]
    n_res = len(pep_res)

    # Z-coordinate (FAST)
    pep_z = np.mean(traj.xyz[:, peptide, 2], axis=1)

    # Per-residue z (FAST)
    per_res_z = {}
    for i, res in enumerate(pep_res):
        res_atoms = traj.topology.select(f'resid {res.index}')
        if len(res_atoms) > 0:
            per_res_z[f"{res.name}{i+1}"] = round(float(np.mean(traj.xyz[:, res_atoms, 2])), 3)

    # Rg (FAST)
    pep_traj = traj.atom_slice(peptide)
    rg = md.compute_rg(pep_traj)

    # E2E (FAST)
    ca = traj.topology.select('name CA')
    e2e = md.compute_distances(traj, [[ca[0], ca[-1]]]) if len(ca) >= 2 else None

    # DSSP (FAST)
    dssp = md.compute_dssp(traj)

    # SASA on PEPTIDE ONLY (much faster than whole system)
    sasa = md.shrake_rupley(pep_traj, mode='residue')
    total_sasa = np.sum(sasa[:, :n_res], axis=1)
    mean_per_res_sasa = np.mean(sasa[:, :n_res], axis=0)

    # Min z per residue (closest approach)
    min_z_per_res = {}
    for i, res in enumerate(pep_res):
        res_atoms = traj.topology.select(f'resid {res.index}')
        if len(res_atoms) > 0:
            min_z_per_res[f"{res.name}{i+1}"] = round(float(np.min(traj.xyz[:, res_atoms, 2])), 3)

    return {
        "name": name, "label": label, "charge": charge, "sequence": sequence,
        "n_frames": int(traj.n_frames),
        "z": {
            "mean": round(float(np.mean(pep_z)), 4),
            "std": round(float(np.std(pep_z)), 4),
            "min": round(float(np.min(pep_z)), 4),
            "initial": round(float(pep_z[0]), 4),
            "final": round(float(pep_z[-1]), 4),
            "delta": round(float(pep_z[-1] - pep_z[0]), 4),
        },
        "sasa": {
            "mean": round(float(np.mean(total_sasa)), 3),
            "std": round(float(np.std(total_sasa)), 3),
            "per_residue": {f"{pep_res[i].name}{i+1}": round(float(mean_per_res_sasa[i]), 3) for i in range(min(n_res, len(mean_per_res_sasa)))},
        },
        "rg": {"mean": round(float(np.mean(rg)), 3), "std": round(float(np.std(rg)), 3)},
        "e2e": {"mean": round(float(np.mean(e2e)), 3)} if e2e is not None else None,
        "ss": {k: round(float(np.mean(dssp[:,:n_res] == v))*100, 1) for v,k in [('H','helix'),('E','sheet'),('C','coil')]},
        "per_residue_z": per_res_z,
        "per_residue_min_z": min_z_per_res,
    }

if __name__ == "__main__":
    print("Analyzing 3 systems in parallel...")
    with Pool(3) as pool:
        results_list = pool.map(analyze, SYSTEMS)

    all_results = {r["name"]: r for r in results_list}

    # Comparison table
    print(f"\n{'System':<30} {'Chg':>4} {'z_mean':>8} {'z_min':>8} {'Δz':>8} {'SASA':>8} {'Rg':>8}")
    print("-" * 80)
    for r in results_list:
        print(f"{r['label']:<30} {r['charge']:>+4d} "
              f"{r['z']['mean']:>8.3f} {r['z']['min']:>8.3f} "
              f"{r['z']['delta']:>+8.3f} {r['sasa']['mean']:>8.3f} "
              f"{r['rg']['mean']:>8.3f}")

    # Per-residue z
    print(f"\nPer-residue mean z (nm):")
    for r in results_list:
        print(f"  {r['label']}:")
        for res, z in r['per_residue_z'].items():
            min_z = r['per_residue_min_z'].get(res, 0)
            print(f"    {res}: mean={z:.3f}, min={min_z:.3f}")

    # Per-residue SASA
    print(f"\nPer-residue SASA (nm²):")
    for r in results_list:
        print(f"  {r['label']}:")
        for res, s in r['sasa']['per_residue'].items():
            print(f"    {res}: {s:.3f}")

    # v5 predictions
    wt = all_results["KKRPKP_membrane"]
    ctrl = all_results["NNRPNP_membrane"]
    ext = all_results["PrP_23_35_membrane"]

    print(f"\n{'='*60}")
    print("v5 PREDICTIONS vs RESULTS")
    print(f"{'='*60}")
    print(f"  1. KKRPKP closer to membrane than NNRPNP?")
    print(f"     KK z={wt['z']['mean']:.3f}, NN z={ctrl['z']['mean']:.3f}")
    print(f"     {'CONFIRMED' if wt['z']['mean'] < ctrl['z']['mean'] else 'REFUTED'}")
    print(f"  2. KKRPKP approached more (delta z more negative)?")
    print(f"     KK Δz={wt['z']['delta']:+.3f}, NN Δz={ctrl['z']['delta']:+.3f}")
    print(f"     {'CONFIRMED' if wt['z']['delta'] < ctrl['z']['delta'] else 'REFUTED'}")
    print(f"  3. PrP 23-35 (Trp) closer than KKRPKP?")
    print(f"     P35 z_min={ext['z']['min']:.3f}, KK z_min={wt['z']['min']:.3f}")
    print(f"     {'CONFIRMED (Trp anchors)' if ext['z']['min'] < wt['z']['min'] else 'NOT CONFIRMED'}")

    with open(os.path.join(RESULTS_DIR, "analysis", "comparison.json"), "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\nSaved to results_mimetic/analysis/comparison.json")
