"""Analyze AlphaFold PrP structure for auto-inhibitory N-terminal contacts."""
import numpy as np
import json
import os

PDB_PATH = "/Users/allan/Projects/cjd/results_alphafold/prp_alphafold.pdb"
PLDDT_PATH = "/Users/allan/Projects/cjd/results_alphafold/prp_plddt.json"
OUT_PATH = "/Users/allan/Projects/cjd/results_alphafold/alphafold_analysis.json"

# v5 domain definitions (1-indexed residue numbers in UniProt sequence)
SIGNAL = (1, 22)
NTERMINAL = (23, 93)
LINKER = (94, 111)
HYDROPHOBIC = (112, 133)
GLOBULAR = (128, 225)
GPI_SIGNAL = (226, 253)

CONTACT_CUTOFF = 0.8  # nm

# Parse PDB
def parse_pdb(path):
    atoms = []
    for line in open(path):
        if line.startswith("ATOM"):
            atom = {
                'serial': int(line[6:11]),
                'name': line[12:16].strip(),
                'resname': line[17:20].strip(),
                'chain': line[21],
                'resid': int(line[22:26]),
                'x': float(line[30:38]) / 10.0,  # A -> nm
                'y': float(line[38:46]) / 10.0,
                'z': float(line[46:54]) / 10.0,
                'bfactor': float(line[60:66]),  # pLDDT in AlphaFold
            }
            atoms.append(atom)
    return atoms

atoms = parse_pdb(PDB_PATH)
print(f"Loaded {len(atoms)} atoms, residues {atoms[0]['resid']}-{atoms[-1]['resid']}")

# Get CA atoms for distance calculations
ca_atoms = [a for a in atoms if a['name'] == 'CA']
print(f"CA atoms: {len(ca_atoms)}")

# 1. pLDDT by region
def region_plddt(ca_list, start, end):
    vals = [a['bfactor'] for a in ca_list if start <= a['resid'] <= end]
    return np.mean(vals) if vals else 0

regions = {
    'Signal peptide (1-22)': SIGNAL,
    'N-terminal tail (23-93)': NTERMINAL,
    'Linker (94-111)': LINKER,
    'Hydrophobic gatekeeper (112-133)': HYDROPHOBIC,
    'Globular domain (128-225)': GLOBULAR,
    'GPI signal (226-253)': GPI_SIGNAL,
}

print("\n=== pLDDT by v5 region ===")
plddt_by_region = {}
for name, (s, e) in regions.items():
    p = region_plddt(ca_atoms, s, e)
    plddt_by_region[name] = round(p, 1)
    conf = "HIGH" if p > 70 else "LOW" if p > 50 else "VERY LOW (disordered)"
    print(f"  {name}: {p:.1f} ({conf})")

# 2. Inter-domain distances (N-terminal to globular)
print("\n=== N-terminal <-> Globular domain distances ===")
nterm_ca = [a for a in ca_atoms if NTERMINAL[0] <= a['resid'] <= NTERMINAL[1]]
glob_ca = [a for a in ca_atoms if GLOBULAR[0] <= a['resid'] <= GLOBULAR[1]]

contacts = []
min_dist = 999
for na in nterm_ca:
    for ga in glob_ca:
        d = np.sqrt((na['x']-ga['x'])**2 + (na['y']-ga['y'])**2 + (na['z']-ga['z'])**2)
        if d < CONTACT_CUTOFF:
            contacts.append({
                'nterm_res': f"{na['resname']}{na['resid']}",
                'glob_res': f"{ga['resname']}{ga['resid']}",
                'distance_nm': round(d, 3),
            })
        if d < min_dist:
            min_dist = d
            closest = (f"{na['resname']}{na['resid']}", f"{ga['resname']}{ga['resid']}")

print(f"  N-terminal residues: {len(nterm_ca)}")
print(f"  Globular residues: {len(glob_ca)}")
print(f"  Contacts within {CONTACT_CUTOFF} nm: {len(contacts)}")
print(f"  Closest pair: {closest[0]} - {closest[1]} at {min_dist:.3f} nm")

if contacts:
    print(f"  Contact pairs:")
    for c in sorted(contacts, key=lambda x: x['distance_nm'])[:15]:
        print(f"    {c['nterm_res']:>8} - {c['glob_res']:<8} : {c['distance_nm']:.3f} nm")

# 3. Residue 127 (G127V) location
print("\n=== G127V position ===")
res127_atoms = [a for a in atoms if a['resid'] == 127]
res127_ca = [a for a in ca_atoms if a['resid'] == 127]
if res127_ca:
    r127 = res127_ca[0]
    print(f"  Residue 127: {r127['resname']} at ({r127['x']:.2f}, {r127['y']:.2f}, {r127['z']:.2f}) nm")
    print(f"  pLDDT: {r127['bfactor']:.1f}")

    # Distance to N-terminal contacts
    dists_to_nterm = []
    for na in nterm_ca:
        d = np.sqrt((r127['x']-na['x'])**2 + (r127['y']-na['y'])**2 + (r127['z']-na['z'])**2)
        dists_to_nterm.append((f"{na['resname']}{na['resid']}", d))
    dists_to_nterm.sort(key=lambda x: x[1])
    print(f"  Closest N-terminal residues to G127:")
    for name, d in dists_to_nterm[:5]:
        print(f"    {name}: {d:.3f} nm")

    # Is G127 at the N-term/globular interface?
    near_nterm = any(d < 1.0 for _, d in dists_to_nterm)
    near_glob = any(
        np.sqrt((r127['x']-g['x'])**2 + (r127['y']-g['y'])**2 + (r127['z']-g['z'])**2) < 1.0
        for g in glob_ca
    )
    print(f"  G127 near N-terminal (<1 nm): {near_nterm}")
    print(f"  G127 near globular (<1 nm): {near_glob}")
    if near_nterm and near_glob:
        print(f"  -> G127 IS at the N-terminal/globular interface!")
    elif near_glob:
        print(f"  -> G127 is at the globular domain surface (gatekeeper region)")

# 4. Per-residue pLDDT for N-terminal
print("\n=== Per-residue pLDDT (N-terminal 23-93) ===")
nterm_plddt = []
for a in ca_atoms:
    if NTERMINAL[0] <= a['resid'] <= NTERMINAL[1]:
        nterm_plddt.append({'resid': a['resid'], 'resname': a['resname'], 'plddt': a['bfactor']})
        status = "■" if a['bfactor'] > 70 else "□" if a['bfactor'] > 50 else "·"
        print(f"  {a['resname']}{a['resid']:>4}: {a['bfactor']:5.1f} {status}")

# 5. Center of mass distances between regions
def com(atom_list):
    x = np.mean([a['x'] for a in atom_list])
    y = np.mean([a['y'] for a in atom_list])
    z = np.mean([a['z'] for a in atom_list])
    return np.array([x, y, z])

nterm_com = com(nterm_ca)
glob_com = com(glob_ca)
hydro_ca = [a for a in ca_atoms if HYDROPHOBIC[0] <= a['resid'] <= HYDROPHOBIC[1]]
hydro_com = com(hydro_ca) if hydro_ca else np.zeros(3)

com_dist_ng = np.linalg.norm(nterm_com - glob_com)
com_dist_nh = np.linalg.norm(nterm_com - hydro_com)

print(f"\n=== Center-of-mass distances ===")
print(f"  N-terminal COM <-> Globular COM: {com_dist_ng:.3f} nm")
print(f"  N-terminal COM <-> Hydrophobic COM: {com_dist_nh:.3f} nm")

# 6. Save results
results = {
    'pdb_file': PDB_PATH,
    'n_atoms': len(atoms),
    'n_residues': len(ca_atoms),
    'plddt_by_region': plddt_by_region,
    'n_contacts_nterm_glob': len(contacts),
    'contact_cutoff_nm': CONTACT_CUTOFF,
    'contacts': contacts[:20],
    'closest_pair': {'nterm': closest[0], 'glob': closest[1], 'dist_nm': round(min_dist, 3)} if contacts or min_dist < 999 else None,
    'g127': {
        'resname': res127_ca[0]['resname'] if res127_ca else None,
        'plddt': res127_ca[0]['bfactor'] if res127_ca else None,
        'near_nterm': near_nterm if res127_ca else None,
        'near_glob': near_glob if res127_ca else None,
        'at_interface': (near_nterm and near_glob) if res127_ca else None,
    },
    'com_distances': {
        'nterm_glob_nm': round(com_dist_ng, 3),
        'nterm_hydrophobic_nm': round(com_dist_nh, 3),
    },
    'interpretation': '',
}

# Interpretation
if len(contacts) > 0:
    results['interpretation'] = (
        f"AlphaFold predicts {len(contacts)} contacts between N-terminal and globular domain "
        f"(cutoff {CONTACT_CUTOFF} nm), supporting the auto-inhibition model. "
        f"The closest pair is {closest[0]}-{closest[1]} at {min_dist:.3f} nm."
    )
else:
    results['interpretation'] = (
        f"AlphaFold does NOT predict close contacts between N-terminal and globular domain "
        f"(closest pair at {min_dist:.3f} nm, cutoff {CONTACT_CUTOFF} nm). "
        f"This is expected: AlphaFold struggles with disordered regions (N-terminal pLDDT={plddt_by_region.get('N-terminal tail (23-93)', 0):.0f}). "
        f"The auto-inhibition is a dynamic interaction seen by NMR, not a static contact that AlphaFold captures."
    )

with open(OUT_PATH, 'w') as f:
    json.dump(results, f, indent=2)

print(f"\n=== SUMMARY ===")
print(results['interpretation'])
print(f"\nResults saved to {OUT_PATH}")
