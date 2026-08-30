"""Analyze a single mimetic system. Run in parallel for each peptide."""
import mdtraj as md
import numpy as np
import json
import sys
import os

name = sys.argv[1]
label = sys.argv[2]
charge = int(sys.argv[3])
sequence = sys.argv[4]

RESULTS_DIR = "/Users/allan/Projects/cjd/results_mimetic"
traj_path = os.path.join(RESULTS_DIR, f"{name}_traj.dcd")
top_path = os.path.join(RESULTS_DIR, f"{name}_system.pdb")

print(f"Loading {name}...")
traj = md.load(traj_path, top=top_path)
peptide = traj.topology.select('protein')
pep_res = [r for r in traj.topology.residues if r.name not in ('HOH','NA','CL')]
n_res = len(pep_res)
print(f"  {traj.n_frames} frames, {n_res} residues")

# Z-coordinate
pep_z = np.mean(traj.xyz[:, peptide, 2], axis=1)
print(f"  z: {pep_z[0]:.3f} -> {pep_z[-1]:.3f} (delta {pep_z[-1]-pep_z[0]:+.3f})")

# SASA
print(f"  Computing SASA...")
sasa = md.shrake_rupley(traj, mode='residue')
total_sasa = np.sum(sasa[:, :n_res], axis=1)
mean_per_res = np.mean(sasa[:, :n_res], axis=0)

# Rg
rg = md.compute_rg(traj, masses=np.ones(traj.n_atoms))

# E2E
ca = traj.topology.select('name CA')
e2e = md.compute_distances(traj, [[ca[0], ca[-1]]]) if len(ca) >= 2 else None

# DSSP
dssp = md.compute_dssp(traj)

# Per-residue z
per_res_z = {}
for i, res in enumerate(pep_res):
    res_atoms = traj.topology.select(f'resid {res.index}')
    if len(res_atoms) > 0:
        per_res_z[f"{res.name}{i+1}"] = round(float(np.mean(traj.xyz[:, res_atoms, 2])), 3)

results = {
    "name": name, "label": label, "charge": charge, "sequence": sequence,
    "n_frames": int(traj.n_frames),
    "z": {
        "mean": round(float(np.mean(pep_z)), 4),
        "std": round(float(np.std(pep_z)), 4),
        "min": round(float(np.min(pep_z)), 4),
        "initial": round(float(pep_z[0]), 4),
        "final": round(float(pep_z[-1]), 4),
        "delta": round(float(pep_z[-1] - pep_z[0]), 4),
        "timeseries": [round(float(z), 4) for z in pep_z[::max(1, len(pep_z)//100)]],
    },
    "sasa": {
        "mean": round(float(np.mean(total_sasa)), 3),
        "std": round(float(np.std(total_sasa)), 3),
        "per_residue": {f"{pep_res[i].name}{i+1}": round(float(mean_per_res[i]), 3) for i in range(min(n_res, len(mean_per_res)))},
    },
    "rg": {"mean": round(float(np.mean(rg)), 3), "std": round(float(np.std(rg)), 3)},
    "e2e": {"mean": round(float(np.mean(e2e)), 3), "std": round(float(np.std(e2e)), 3)} if e2e is not None else None,
    "ss": {k: round(float(np.mean(dssp[:,:n_res] == v))*100, 1) for v,k in [('H','helix'),('E','sheet'),('C','coil')]},
    "per_residue_z": per_res_z,
}

out = os.path.join(RESULTS_DIR, f"{name}_analysis.json")
with open(out, "w") as f:
    json.dump(results, f, indent=2)
print(f"  DONE -> {out}")
