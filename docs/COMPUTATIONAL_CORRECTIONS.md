# Computational corrections — what was wrong and what replaced it

This document records defects found in the molecular dynamics work that
produced `results_steered_real/COMPARISON_DEFINITIVE.json`, and the corrected
protocol that replaces it. It is written to be read *before* any of the MD
numbers elsewhere in this repository, because it retracts some of them.

**Status of the retracted result: the KKRPKP-vs-NNRPNP steered-MD comparison
does not support the conclusion that was drawn from it.** The claim was that a
charged peptide is drawn into the membrane while a neutral one is pushed out.
That claim may still be true — it is being retested — but the simulation
offered as evidence could not have shown it.

---

## Defect 1 — the peptide was never in water

The single most serious problem. In both CHARMM-GUI systems, the peptide was
placed at the **centre of the bilayer**, buried in the hydrocarbon core:

| | peptide mean z | bilayer centre | upper phosphate plane | waters within 6 Å |
|---|---|---|---|---|
| KKRPKP (old) | 42.7 Å | 42.5 Å | 61.1 Å | **0** |
| NNRPNP (old) | 45.5 Å | 45.7 Å | 64.6 Å | **0** |

The lipid tail termini span z = 39–46 Å, so the peptide sat squarely in the
apolar core with **zero hydration**. For KKRPKP that means three lysines and
an arginine — four full charges — desolvated in hydrocarbon, which costs on
the order of hundreds of kJ/mol. It is not a state the peptide ever occupies.

The cause is a CHARMM-GUI default: leaving "Translate Molecule along Z axis"
unchecked centres the uploaded molecule on the bilayer midplane. That is the
right default for a transmembrane protein and the wrong one for a peripheral
peptide.

## Defect 2 — the steered MD barely moved

The pull ran from `z_start` to the membrane centre. Since the peptide already
*was* at the membrane centre, that is a displacement of **0.2 Å over 15 ns**.
No insertion was sampled. What the run actually measured was whether a
peptide trapped in an unphysical state drifts a little up or a little down
while weakly restrained — positional jitter, not affinity.

## Defect 3 — the work integral was wrong

In `steered_md_real.py:157`:

```python
cum_work = np.cumsum([d['dz']*K_PULL*PULL_RATE*DT for d in z_data])
```

Two errors. The restraint's displacement per window is
`(z_end - z_start)/n_windows`, but the code used `PULL_RATE*DT` — a nominal
rate that never fed the `np.linspace` actually generating the targets. For the
run as configured the two differ by a factor of ~222. And `dW = -k·dz·dz₀`
carries a minus sign that was missing, left implicit in how the sign was read
off afterwards.

## Defect 4 — n = 1

Steered-MD work is a fluctuating quantity. Two single trajectories cannot be
compared without a spread. The reported values were **−0.15 and +0.07 kJ/mol**
against kT = 2.5 kJ/mol at 303 K — that is, quantities 20-40× *smaller* than
thermal noise were presented as a difference in sign. (Correcting Defect 3
scales them to roughly −33 and +15 kJ/mol, comfortably above kT, but that only
fixes the units; it does not create an error bar.)

## Defect 5 — the control was mislabelled

NNRPNP was described as the neutral (charge 0) control. It retains its
arginine, so it carries **+1**, not 0. The contrast is +4 vs +1 — a
three-charge difference, still a real contrast, but not the one stated.

---

## The corrected protocol

Rather than repair the steered-MD setup, the experiment was redesigned.
Steered MD requires a correct work integral and a long pull path; both were
sources of error here. **Unbiased adsorption MD** removes the failure mode
entirely: the peptide starts solvated above the bilayer and either binds or
does not, and every observable is geometric. Nothing depends on a force
constant or a pulling rate.

### Systems (rebuilt)

All three place the peptide in bulk water above the bilayer, POPC:POPS 80:20,
0.15 M KCl, CHARMM36m, 303.15 K, NPT.

| system | CHARMM-GUI job | atoms | box (Å) | peptide above phosphates |
|---|---|---|---|---|
| KKRPKP (+4) | 8813850296 | 34,962 | 65.9 × 65.9 × 87.1 | **+12.7 Å** |
| NNRPNP (+1) | 8813960988 | 35,253 | 65.9 × 65.9 × 87.1 | **+12 Å** |
| PrP 23-93 | 8813640350 | 79,492 | 84.9 × 84.9 × 117.1 | **+10.5 Å** |

The two hexapeptide systems use the identical 66 Å patch as the original
runs, so only the peptide's starting position differs. CHARMM-GUI's own
message — *"For this system, insertion method can not be used"* — is an
independent confirmation that the peptide no longer intersects the bilayer.

### Observables

- **bound fraction** — frames with any peptide-lipid heavy-atom pair < 4 Å
- **first contact** — time to binding
- **z_rel** — peptide COM relative to the upper phosphate plane (negative = inserted)
- **desolvation** — waters within 4 Å of the peptide
- **POPS enrichment** — fraction of lipid contacts made to the anionic lipid,
  divided by the 20% baseline. This is an *internal* control: it reads charge
  preference out of a single simulation, without needing the sequence-matched
  peptide.

### Statistics

4 independent replicas per system, each with its own velocity seed. Pooled to
mean ± SEM, with a Welch test between systems.

---

## Two build details worth recording

Both cost real time and are not obvious.

**Atom ordering.** CHARMM-GUI's PDB reader identifies residues by the N-CA-C-O
backbone sequence. OpenMM writes `N, CA, C, CB, O, ...` — O after CB. Uploading
an OpenMM-written PDB made CHARMM-GUI split a continuous 71-residue chain into
three fragments (1-63, 64, 65-71), which would have been simulated as three
peptides with spurious charged termini. Rewriting atoms in canonical order
fixed it.

**Rotation origin.** CHARMM-GUI centres an uploaded molecule in x and y but
leaves z alone, and Positioning rotations are about the origin. Peptides cut
from the old membrane systems had z ≈ 31-54 Å, so a 90° rotation about X
mapped that offset into y and displaced the peptide 43 Å sideways — beside the
patch rather than above it (`yextent` came back as 54.3 Å for a 23 Å peptide).
Centring and pre-orienting locally before upload removes the problem.

**Downloading large systems.** The `download.tgz` link truncates at
41,876,032 bytes for a system of this size. The individual files are served as
static paths and can be fetched directly:
`https://charmm-gui.org/uploaded_pdb/<jobid>/step5_assembly.psf`. Box
dimensions are in `step5_assembly.str`.

---

## What still stands

Nothing in the literature review, the kinetic model, the sequence analyses
(charge, Wimley-White, AMP scoring, pH-dependent aggregation), or the
therapeutic assessment depends on the retracted MD. Those results are
independent of it. What is retracted is specifically the claim that the
membrane simulations demonstrated charge-driven insertion.
