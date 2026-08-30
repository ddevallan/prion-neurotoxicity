"""
Analyze mimetic membrane simulation results.
Compares KKRPKP (+4) vs NNRPNP (0) vs PrP 23-35 (+4 with Trp).

Measures:
1. Distance to membrane surface (z-coordinate over time)
2. SASA evolution
3. Radius of gyration
4. Per-residue analysis
5. Approach velocity (how fast does peptide move toward membrane)
"""
import mdtraj as md
import numpy as np
import json
import os
import sys

RESULTS_DIR = "/Users/allan/Projects/cjd/results_mimetic"
OUTPUT_DIR = "/Users/allan/Projects/cjd/results_mimetic/analysis"
os.makedirs(OUTPUT_DIR, exist_ok=True)

SYSTEMS = {
    "KKRPKP_membrane": {
        "label": "KKRPKP (+4, wild-type)",
        "charge": 4,
        "sequence": "KKRPKP",
    },
    "NNRPNP_membrane": {
        "label": "NNRPNP (0, neutral control)",
        "charge": 0,
        "sequence": "NNRPNP",
    },
    "PrP_23_35_membrane": {
        "label": "PrP 23-35 (+4, with Trp)",
        "charge": 4,
        "sequence": "KKRPKPGGWNTGG",
    },
}

def analyze_system(name, info):
    traj_path = os.path.join(RESULTS_DIR, f"{name}_traj.dcd")
    top_path = os.path.join(RESULTS_DIR, f"{name}_system.pdb")

    if not os.path.exists(traj_path) or not os.path.exists(top_path):
        print(f"  SKIP {name}: files not found")
        return None

    print(f"\n{'='*55}")
    print(f"  {info['label']}")
    print(f"  Sequence: {info['sequence']}, Charge: {info['charge']:+d}")
    print(f"{'='*55}")

    traj = md.load(traj_path, top=top_path)
    print(f"  Frames: {traj.n_frames}, Atoms: {traj.n_atoms}")

    peptide = traj.topology.select('protein')
    pep_residues = [r for r in traj.topology.residues
                    if r.name not in ('HOH', 'NA', 'CL')]
    n_res = len(pep_residues)
    print(f"  Peptide: {n_res} residues, {len(peptide)} atoms")

    results = {
        "name": name,
        "label": info["label"],
        "charge": info["charge"],
        "sequence": info["sequence"],
        "n_frames": int(traj.n_frames),
    }

    # 1. Z-coordinate (distance to membrane surface at z=0)
    print("  Computing z-distance to membrane...")
    pep_z = np.mean(traj.xyz[:, peptide, 2], axis=1)
    results["z"] = {
        "mean": round(float(np.mean(pep_z)), 4),
        "std": round(float(np.std(pep_z)), 4),
        "min": round(float(np.min(pep_z)), 4),
        "max": round(float(np.max(pep_z)), 4),
        "final": round(float(pep_z[-1]), 4),
        "initial": round(float(pep_z[0]), 4),
        "timeseries": [round(float(z), 4) for z in pep_z[::max(1, len(pep_z)//100)]],
    }
    print(f"    Initial z: {pep_z[0]:.3f} nm")
    print(f"    Final z:   {pep_z[-1]:.3f} nm")
    print(f"    Min z:     {np.min(pep_z):.3f} nm (closest approach)")
    print(f"    Mean z:    {np.mean(pep_z):.3f} +/- {np.std(pep_z):.3f} nm")

    # Approach: did it move toward membrane (z=0)?
    delta_z = pep_z[-1] - pep_z[0]
    results["z"]["delta"] = round(float(delta_z), 4)
    print(f"    Delta z:   {delta_z:+.3f} nm ({'APPROACHED' if delta_z < 0 else 'MOVED AWAY'})")

    # 2. SASA
    print("  Computing SASA...")
    sasa = md.shrake_rupley(traj, mode='residue')
    total_sasa = np.sum(sasa[:, :n_res], axis=1)
    mean_per_res = np.mean(sasa[:, :n_res], axis=0)

    results["sasa"] = {
        "mean": round(float(np.mean(total_sasa)), 3),
        "std": round(float(np.std(total_sasa)), 3),
        "per_residue": {
            f"{pep_residues[i].name}{i+1}": round(float(mean_per_res[i]), 3)
            for i in range(min(n_res, len(mean_per_res)))
        },
    }
    print(f"    SASA: {np.mean(total_sasa):.3f} +/- {np.std(total_sasa):.3f} nm²")
    for i, r in enumerate(pep_residues):
        if i < len(mean_per_res):
            print(f"      {r.name}{i+1}: {mean_per_res[i]:.3f} nm²")

    # 3. Rg
    print("  Computing Rg...")
    rg = md.compute_rg(traj, masses=np.ones(traj.n_atoms))
    results["rg"] = {
        "mean": round(float(np.mean(rg)), 3),
        "std": round(float(np.std(rg)), 3),
    }
    print(f"    Rg: {np.mean(rg):.3f} +/- {np.std(rg):.3f} nm")

    # 4. End-to-end distance
    ca = traj.topology.select('name CA')
    if len(ca) >= 2:
        e2e = md.compute_distances(traj, [[ca[0], ca[-1]]])
        results["e2e"] = {
            "mean": round(float(np.mean(e2e)), 3),
            "std": round(float(np.std(e2e)), 3),
        }
        print(f"    E2E: {np.mean(e2e):.3f} +/- {np.std(e2e):.3f} nm")

    # 5. Secondary structure
    print("  Computing DSSP...")
    dssp = md.compute_dssp(traj)
    ss = {}
    for code, label in [('H', 'helix'), ('E', 'sheet'), ('C', 'coil')]:
        frac = float(np.mean(dssp[:, :n_res] == code)) * 100
        ss[label] = round(frac, 1)
    results["ss"] = ss
    print(f"    SS: helix={ss['helix']:.1f}%, sheet={ss['sheet']:.1f}%, coil={ss['coil']:.1f}%")

    # 6. Z-position per residue (which residues are closest to membrane?)
    print("  Per-residue z-position...")
    per_res_z = {}
    for i, res in enumerate(pep_residues):
        res_atoms = traj.topology.select(f'resid {res.index}')
        if len(res_atoms) > 0:
            res_z = np.mean(traj.xyz[:, res_atoms, 2])
            per_res_z[f"{res.name}{i+1}"] = round(float(res_z), 3)
            print(f"      {res.name}{i+1}: z = {res_z:.3f} nm")
    results["per_residue_z"] = per_res_z

    return results

# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    print("=" * 65)
    print("MIMETIC MEMBRANE SIMULATION ANALYSIS")
    print("=" * 65)

    all_results = {}

    for name, info in SYSTEMS.items():
        res = analyze_system(name, info)
        if res:
            all_results[name] = res

    if len(all_results) < 2:
        print("\nNot enough systems to compare. Check that files are downloaded.")
        sys.exit(1)

    # ============================================================
    # COMPARISON
    # ============================================================

    print(f"\n{'='*65}")
    print("COMPARISON: KKRPKP vs NNRPNP vs PrP 23-35")
    print(f"{'='*65}")

    print(f"\n{'System':<30} {'Charge':>7} {'z_mean':>8} {'z_min':>8} {'delta_z':>8} {'SASA':>8} {'Rg':>8}")
    print("-" * 80)

    for name, res in all_results.items():
        print(f"{res['label']:<30} {res['charge']:>+7d} "
              f"{res['z']['mean']:>8.3f} {res['z']['min']:>8.3f} "
              f"{res['z']['delta']:>+8.3f} {res['sasa']['mean']:>8.3f} "
              f"{res['rg']['mean']:>8.3f}")

    # v5 predictions
    wt = all_results.get("KKRPKP_membrane", {})
    ctrl = all_results.get("NNRPNP_membrane", {})
    ext = all_results.get("PrP_23_35_membrane", {})

    print(f"\nv5 MODEL PREDICTIONS vs RESULTS:")
    print(f"-" * 60)

    if wt and ctrl:
        wt_closer = wt['z']['mean'] < ctrl['z']['mean']
        print(f"  1. KKRPKP closer to membrane than NNRPNP?")
        print(f"     KKRPKP z={wt['z']['mean']:.3f}, NNRPNP z={ctrl['z']['mean']:.3f}")
        print(f"     {'YES — CONFIRMED' if wt_closer else 'NO — REFUTED'}")

        wt_approached = wt['z']['delta'] < ctrl['z']['delta']
        print(f"  2. KKRPKP approached membrane more?")
        print(f"     KKRPKP delta={wt['z']['delta']:+.3f}, NNRPNP delta={ctrl['z']['delta']:+.3f}")
        print(f"     {'YES — CONFIRMED' if wt_approached else 'NO — REFUTED'}")

    if ext and wt:
        trp_closer = ext['z']['min'] < wt['z']['min']
        print(f"  3. PrP 23-35 (with Trp) closer than KKRPKP alone?")
        print(f"     PrP_23-35 z_min={ext['z']['min']:.3f}, KKRPKP z_min={wt['z']['min']:.3f}")
        print(f"     {'YES — Trp anchors confirmed' if trp_closer else 'NO — Trp anchoring not observed'}")

    # Save
    out_path = os.path.join(OUTPUT_DIR, "mimetic_comparison.json")
    with open(out_path, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\nResults saved to {out_path}")
