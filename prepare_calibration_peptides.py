"""Prepare PDBs for calibration and PrP 23-93 membrane simulations.

Peptides:
  - PrP 23-93: extracted from AlphaFold (physiologically relevant conformation)
  - LL-37: human cathelicidin, the benchmark AMP for calibration
  - Melittin: second benchmark (bee venom, most-studied membrane peptide)
"""
import numpy as np
import os
from pdbfixer import PDBFixer
import openmm.app as app

OUT = "/Users/allan/Projects/cjd/calibration_peptides"
os.makedirs(OUT, exist_ok=True)

AA3 = {'A':'ALA','R':'ARG','N':'ASN','D':'ASP','C':'CYS','E':'GLU','Q':'GLN',
       'G':'GLY','H':'HIS','I':'ILE','L':'LEU','K':'LYS','M':'MET','F':'PHE',
       'P':'PRO','S':'SER','T':'THR','W':'TRP','Y':'TYR','V':'VAL'}


def extract_range_from_pdb(src, first, last, out_path):
    """Extract a residue range from a PDB, keeping only ATOM records."""
    kept = []
    for line in open(src):
        if line.startswith('ATOM'):
            resid = int(line[22:26])
            if first <= resid <= last:
                kept.append(line)
    with open(out_path, 'w') as f:
        f.writelines(kept)
        f.write('TER\nEND\n')

    coords = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in kept])
    extent = coords.max(axis=0) - coords.min(axis=0)
    n_res = len({int(l[22:26]) for l in kept})
    print(f"  {os.path.basename(out_path)}: {n_res} residues, {len(kept)} atoms")
    print(f"    extent: {extent[0]:.1f} x {extent[1]:.1f} x {extent[2]:.1f} A")
    return extent


def build_helix(seq, name):
    """Build peptide as an alpha-helix (correct for LL-37 and melittin,
    which are helical when membrane-bound)."""
    # Ideal alpha-helix: rise 1.5 A/residue, 100 deg rotation, radius 2.3 A
    lines = ["HEADER    PEPTIDE"]
    idx = 1
    for i, aa in enumerate(seq):
        resn = AA3.get(aa, 'ALA')
        phase = np.radians(100.0 * i)
        z = 1.5 * i
        # Backbone atoms placed around the helical axis
        for atom, r, dphi, dz in [('N', 1.5, -0.35, -0.5),
                                  ('CA', 2.3, 0.0, 0.0),
                                  ('C', 2.0, 0.35, 0.5),
                                  ('O', 2.4, 0.55, 0.7)]:
            x = r * np.cos(phase + dphi)
            y = r * np.sin(phase + dphi)
            lines.append(
                f"ATOM  {idx:5d} {atom:<4s} {resn:3s} A{i+1:4d}    "
                f"{x:8.3f}{y:8.3f}{z+dz:8.3f}  1.00  0.00           {atom[0]:>2s}")
            idx += 1
    lines += ["TER", "END"]

    raw = os.path.join(OUT, f"{name}_raw.pdb")
    with open(raw, 'w') as f:
        f.write("\n".join(lines))

    fixer = PDBFixer(filename=raw)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)
    out = os.path.join(OUT, f"{name}.pdb")
    with open(out, 'w') as f:
        app.PDBFile.writeFile(fixer.topology, fixer.positions, f)

    pos = np.array(fixer.positions.value_in_unit(
        __import__('openmm.unit', fromlist=['angstroms']).angstroms))
    extent = pos.max(axis=0) - pos.min(axis=0)
    print(f"  {name}.pdb: {len(seq)} residues, {len(pos)} atoms")
    print(f"    extent: {extent[0]:.1f} x {extent[1]:.1f} x {extent[2]:.1f} A")
    return out


print("=" * 60)
print("PREPARING CALIBRATION PEPTIDES")
print("=" * 60)

print("\n1. PrP 23-93 (from AlphaFold, disordered N-terminal):")
extract_range_from_pdb(
    "/Users/allan/Projects/cjd/results_alphafold/prp_alphafold.pdb",
    23, 93, os.path.join(OUT, "PrP_23_93.pdb"))

print("\n2. LL-37 (human cathelicidin, benchmark AMP, alpha-helix):")
build_helix("LLGDFFRKSKEKIGKEFKRIVQRIKDFLRNLVPRTES", "LL37")

print("\n3. Melittin (bee venom, benchmark, alpha-helix):")
build_helix("GIGAVLKVLTTGLPALISWIKRKRQQ", "Melittin")

print(f"\nAll PDBs in {OUT}/")
print("\nNOTE: extents above determine the CHARMM-GUI box size needed.")
print("A peptide longer than ~60 A in any dimension needs a bigger membrane patch.")
