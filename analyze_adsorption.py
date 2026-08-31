"""Characterise the binding mode from the adsorption trajectories.

z_rel is measured against the phosphorus plane, but phosphorus is not the
membrane surface: choline nitrogens and the glycerol/serine headgroups sit
above it. A positive z_rel therefore does not by itself mean the peptide
failed to engage the bilayer, and a peptide's centre of mass sits further out
the larger the peptide is. This script puts the numbers on a scale that can be
read: where the headgroup layer ends, and how the contact count and
desolvation evolve.
"""
import json, glob, os
import numpy as np
import openmm.app as app
import openmm.unit as unit

BASE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(BASE, 'results_adsorption')

SYSTEMS = {
    'KKRPKP_water': 'charmm_gui_kkrpkp_water',
    'NNRPNP_water': 'charmm_gui_nnrpnp_water',
    'PrP_23_93': 'charmm_gui_prp2393',
}


def membrane_reference(sysdir):
    """Height of the headgroup layer above the phosphorus plane, and the
    peptide's own radius, both in nm. Measured on the starting structure."""
    pdb = app.PDBFile(os.path.join(BASE, sysdir, 'openmm/step5_input.pdb'))
    pos = pdb.getPositions(asNumpy=True).value_in_unit(unit.nanometers)
    PR = {'ALA','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS','HSD','HSE','HSP',
          'ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'}
    P, outer, pep = [], [], []
    for a, c in zip(pdb.topology.atoms(), pos):
        rn, nm = a.residue.name, a.name
        if rn in ('POPC', 'POPS'):
            if nm == 'P':
                P.append(c)
            # Choline N (POPC) and serine N/carboxyl (POPS) cap the headgroup.
            if nm in ('N', 'C13', 'C14', 'C15', 'O13A', 'O13B'):
                outer.append(c)
        elif rn in PR and not nm.startswith('H'):
            pep.append(c)
    P = np.array(P); outer = np.array(outer); pep = np.array(pep)
    zp = P[:, 2]
    up = zp[zp > zp.mean()]
    zo = outer[:, 2]
    up_o = zo[zo > zo.mean()]
    pep_radius = float(np.linalg.norm(pep - pep.mean(axis=0), axis=1).max())
    return {
        'p_plane': float(up.mean()),
        'headgroup_top': float(up_o.mean()),
        'headgroup_thickness': float(up_o.mean() - up.mean()),
        'peptide_radius': pep_radius,
    }


print("=" * 74)
print("BINDING MODE")
print("=" * 74)

out = {}
for sysname, sysdir in SYSTEMS.items():
    files = sorted(glob.glob(f"{RES}/{sysname}_rep*.json"))
    if not files:
        continue
    ref = membrane_reference(sysdir)
    reps = [json.load(open(f)) for f in files]

    # Per-replica time series, averaged over the last 10 ns.
    #
    # Insertion depth is deliberately NOT estimated as (COM z - peptide
    # radius): that assumes a sphere with its farthest atom pointing straight
    # down. PrP 23-93 is an elongated 55 x 37 x 49 A globule, so the estimate
    # put it 23 A below the phosphate plane while still 2.6 A from the nearest
    # lipid atom, which is impossible. Per-atom depth is not recoverable from
    # the saved observables, so what is reported here is what was measured:
    # approach (COM z over time), contact (closest approach, contact count) and
    # desolvation (waters lost).
    dz, contacts, waters, mindist = [], [], [], []
    for r in reps:
        tr = r['trajectory']
        t = np.array([f['t_ns'] for f in tr])
        late = t > t.max() - 10
        zrel = np.array([f['z_rel'] for f in tr])
        dz.append(zrel[late].mean() - zrel[t < 1].mean())
        contacts.append(np.array([f['contacts_popc'] + f['contacts_pops']
                                  for f in tr])[late].mean())
        waters.append(np.array([f['n_wat_4A'] for f in tr])[late].mean())
        mindist.append(np.array([f['min_dist'] for f in tr])[late].mean())

    dz = np.array(dz); contacts = np.array(contacts)
    waters = np.array(waters); mindist = np.array(mindist)
    w_early = np.array([np.array([f['n_wat_4A'] for f in r['trajectory']
                                  if f['t_ns'] < 1]).mean() for r in reps])

    # Early vs late contact count: does engagement grow after first contact?
    early_c, late_c = [], []
    for r in reps:
        tr = r['trajectory']
        t = np.array([f['t_ns'] for f in tr])
        c = np.array([f['contacts_popc'] + f['contacts_pops'] for f in tr])
        early_c.append(c[t < 5].mean())
        late_c.append(c[t > t.max() - 10].mean())
    early_c = np.array(early_c); late_c = np.array(late_c)

    def pm(a, f="{:+.2f}"):
        s = a.std(ddof=1)/np.sqrt(len(a)) if len(a) > 1 else 0.0
        return f.format(a.mean()) + f" +/- {s:.2f}"

    print(f"\n{sysname}  (n={len(reps)})")
    print(f"  closest approach       : {pm(mindist, '{:.3f}')} nm")
    print(f"  COM moved toward bilayer: {pm(dz)} nm over the run")
    print(f"  lipid contacts (late)  : {pm(contacts, '{:.0f}')} atoms")
    print(f"  contacts early -> late : {early_c.mean():.0f} -> {late_c.mean():.0f}")
    lost = 100 * (1 - waters.mean() / max(w_early.mean(), 1e-9))
    print(f"  waters within 4 A      : {w_early.mean():.0f} -> {waters.mean():.0f} "
          f"({lost:.0f}% desolvated)")

    out[sysname] = {
        'min_dist_nm': float(mindist.mean()),
        'com_displacement_nm': float(dz.mean()),
        'contacts_late': float(contacts.mean()),
        'contacts_early': float(early_c.mean()),
        'waters_4A_early': float(w_early.mean()),
        'waters_4A_late': float(waters.mean()),
    }

with open(f"{RES}/BINDING_MODE.json", 'w') as f:
    json.dump(out, f, indent=2)
print(f"\nWrote {RES}/BINDING_MODE.json")
