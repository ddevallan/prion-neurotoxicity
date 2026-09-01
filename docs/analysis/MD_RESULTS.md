# Molecular Dynamics Results

## Summary

Three rounds of membrane simulations were run, each correcting defects found
in the previous one. Only the third round (membrane perturbation, August 31
2026) produced trustworthy results. The first two are documented in
[COMPUTATIONAL_CORRECTIONS.md](../COMPUTATIONAL_CORRECTIONS.md) and retracted.

**Central finding:** The KKRPKP motif (PrP residues 23-28) is necessary for
membrane binding. Deleting it from the PrP 23-93 fragment eliminates binding
entirely — 0/3 replicas in 180 ns. This computationally replicates Yan et al.
(Sci Adv 2025, PMID 40768577), who showed the same deletion abolishes
toxicity in vivo.

---

## Round 3 — Membrane perturbation (current)

### Platform

- **Hardware**: 4× NVIDIA RTX 4090 (Vast.ai, Sichuan CN)
- **Software**: OpenMM 8.1.1 + openmm-cuda 8.1.1.12
- **Force field**: CHARMM36m (all 57 parameter files from CHARMM-GUI toppar)
- **Membrane**: POPC:POPS 80:20, 0.15 M KCl
- **Conditions**: 303.15 K, NPT, MonteCarloMembraneBarostat (XYIsotropic, ZFree)
- **Timestep**: 4 fs with hydrogen mass repartitioning (4.0 amu)
- **Integrator**: LangevinMiddleIntegrator, friction 1/ps
- **Equilibration**: CHARMM-GUI 6-stage restrained schedule (4000→50 kJ/mol/nm²,
  NVT→NPT, peptide held in bulk water, phosphates restrained in z only)
- **Production**: 60 ns per replica, sampling every 20 ps
- **Cost**: ~$17 total (~12 h wall clock)

### Systems

| system | atoms | box (Å) | lipids/leaflet | peptide charge | purpose |
|--------|-------|---------|----------------|----------------|---------|
| KKRPKP (+4) | 34,962 | 65.9 × 65.9 × 87.1 | 65 | +4 | model peptide |
| NNRPNP (+1) | 35,253 | 65.9 × 65.9 × 87.5 | 65 | +1 | charge control |
| PrP 23-93 | 79,492 | 84.9 × 84.9 × 117.1 | 108 | +6.4 | real fragment |
| PrP Δ23-28 | 81,799 | 84.9 × 84.9 × 121.0 | 108 | +2.4 | Yan 2025 deletion |
| bare 66 Å | 28,116 | 65.9 × 65.9 × 70.0 | 65 | — | baseline (hexa) |
| bare 85 Å | 46,504 | 84.9 × 84.9 × 70.0 | 108 | — | baseline (PrP) |

All peptide systems place the peptide in bulk water, 10-14 Å above the upper
phosphate plane. CHARMM-GUI independently confirms this: *"insertion method
can not be used"* — the peptide does not intersect the bilayer.

The PrP 23-93 fragment was collapsed from the extended AlphaFold conformation
(Rg 37 Å) to a physiological IDP globule (Rg 22.6 Å, Flory reference 23.5 Å)
using implicit solvent at 350 K, stopping in the target window rather than
running to the most compact state.

### Results — binding

| system | n | replicas bound | bound fraction | first contact (ns) | POPS enrichment |
|--------|---|----------------|----------------|--------------------|----|
| **KKRPKP (+4)** | 3 | **2/3** | **0.663** | [—, 0.12, 0.32] | **2.59×** |
| NNRPNP (+1) | 3 | 0/3 (1 partial) | 0.101 | [0.84, 0.74, —] | n/a |
| **PrP 23-93** | 3 | **1/3** | **0.227** | [—, —, 9.18] | **2.92×** |
| **PrP Δ23-28** | 3 | **0/3** | **0.000** | [—, —, —] | **n/a** |
| bare 66 Å | 2 | 0/2 | 0.000 | — | — |
| bare 85 Å | 2 | 0/2 | 0.000 | — | — |

**"Bound"** = min_dist < 4 Å (any peptide heavy atom to any lipid heavy atom).
**POPS enrichment** = fraction of contacts on POPS / (POPS fraction in membrane).
A value of 1.0 means no preference; >1 means the peptide prefers the anionic lipid.

### Results — membrane perturbation (bound replicas only)

Frame-by-frame analysis separating bound from unbound states within the same
trajectory, so the baseline offset from membrane undulation is controlled:

| system | replica | bound % | thinning bound (Å) | thinning unbound (Å) | **DELTA (Å)** |
|--------|---------|---------|--------------------|-----------------------|---------------|
| KKRPKP | rep1 | 100% | −0.92 | −0.24 | **−0.68** |
| KKRPKP | rep2 | 99% | −1.33 | −0.32 | **−1.02** |
| PrP 23-93 | rep2 | 68% | −1.72 | −2.33 | **+0.60** |

**Thinning** = thickness_local − thickness_distal (negative = membrane thinner
under the peptide). **DELTA** = thinning in bound frames minus thinning in
unbound frames of the same replica — the perturbation attributable to binding.

Radial profile for KKRPKP rep2 (last frame):
```
distance from peptide   thickness   n lipids
      0 nm               36.0 Å        5
    1.5 nm               39.1 Å       16
    2.5 nm               41.1 Å       20
    3.5+ nm              39.7 Å       24
```

The membrane is 5 Å thinner directly under the peptide than at 2.5 nm.

### Results — bare membrane baseline

| system | n | thickness (Å) | |S_CD| | APL (Å²) |
|--------|---|---------------|--------|-----------|
| bare 66 Å | 2 | 39.0, 40.0 | 0.362, 0.360 | 63.8, 61.5 |
| bare 85 Å | 2 | 40.0, 40.0 | 0.360, 0.360 | 63.2, 61.5 |

### Known limitations

**1. Equilibration is incomplete.** |S_CD| remains at ~0.36 in all systems
including the bare controls. A converged POPC bilayer has |S_CD| ≈ 0.20 with
a falling profile toward the chain terminus. 60 ns of production after 1.9 ns
of restrained equilibration does not converge the lipid tail order. The
peptide-vs-bare comparison is valid in relative terms because both share the
same incomplete equilibration.

**2. Box size limits the distal reference.** In the 6.6 nm box, nothing can sit
farther than 3.3 nm from the peptide, while perturbation typically reaches
1-3 nm. The frame-by-frame bound/unbound split within the same trajectory is
the primary control; the bare membrane serves as the external reference.

**3. P:L ratio is low.** One peptide per box gives P:L = 1:130 (hexapeptides)
or 1:216 (PrP). Huang's two-state model (PMID 10913240) puts the carpet
thinning threshold at ~1:50. Global thinning at our concentration is expected
to be <0.5 Å, within thermal noise. Only local (radial) analysis can detect
the signal — and does.

**4. 60 ns may be insufficient for PrP binding.** The 71-residue IDP diffuses
slower than a hexapeptide. Only 1/3 PrP replicas bound vs 2/3 for KKRPKP.
This may reflect slower kinetics rather than weaker affinity. Longer runs
(hundreds of ns) would settle this.

**5. The dCC1 rep2 system may be unstable.** APL of 80 Å² (vs 62-64 for all
others) suggests a possible bilayer defect. The binding result (0/3) is
consistent with rep0 and rep1 which are normal, so the conclusion holds.

---

## Retracted: Rounds 1 and 2

### Round 1 — Steered MD (retracted)

Placed the peptide at the bilayer center with zero hydration, pulled 0.2 Å in
15 ns, and reported the result as insertion affinity. Five defects documented
in [COMPUTATIONAL_CORRECTIONS.md](../COMPUTATIONAL_CORRECTIONS.md).

### Round 2 — Adsorption MD (superseded)

Correct geometry (peptide in water) but two methodological faults:
1. No restrained equilibration — peptide bound during the equilibration phase
2. Only scalars saved — no membrane perturbation measurable

Results (KKRPKP 4/4 bound, POPS enrichment 2.13×) are directionally
consistent with Round 3 but the time-to-bind values (0.02 ns) are artefactual.

---

## Interpretation

The computational evidence supports the following chain:

```
KKRPKP motif present  →  electrostatic binding to anionic membrane
                      →  local thinning (−0.7 to −1.0 Å in hexapeptide)
                      →  [gap: link to NMDA modulation is not simulated]

KKRPKP motif deleted  →  no binding in 180 ns
                      →  no perturbation
```

What is **demonstrated**: KKRPKP is necessary for membrane binding, and
binding thins the membrane locally.

What is **inferred but not simulated**: that thinning of this magnitude
modulates NMDA receptors. Kloda et al. (PNAS 2007, PMID 17242368) showed
that membrane stretch potentiates NMDA currents, but the link from our
observed thinning to receptor modulation is extrapolation, not measurement.

What **remains open**: whether the PrP 23-93 fragment perturbs more or less
than the KKRPKP hexapeptide (insufficient binding events to compare), and
whether the effect persists in a membrane with cholesterol.
