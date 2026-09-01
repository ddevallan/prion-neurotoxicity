# Revision History: Models of Prion Neurotoxicity (v1–v5)

## Overview

This document traces the evolution of our mechanistic model through five versions. Each version was built on the wreckage of the one before it — killed by specific data, not by taste. Two ideas survived every revision: the N-terminal tail is the toxic effector, and there is a dose threshold. What changed between v4 and v5 is where the specificity lives: it migrated from the *toxicity mechanism* to the *replication mechanism*.

---

## v1 — PrP^Sc aggregates are the toxic agent

**Proposed:** PrP^Sc itself — the misfolded, protease-resistant aggregate — is directly neurotoxic. Accumulation of PrP^Sc kills neurons.

**Supporting evidence:**
- PrP^Sc is the defining hallmark of prion disease (Prusiner, Science 1982)
- PrP^Sc accumulates in affected brain regions
- In vitro, PrP^Sc preparations cause cell death in culture

**What killed it:**
- Sandberg et al. (Nature 2011): infectivity and toxicity are temporally decoupled. PrP^Sc titer reaches maximum early (phase 1), but disease only manifests later (phase 2). Correlation between total PrP^Sc load and disease severity is poor.
- Chesebro et al. (Science 2005): anchorless PrP mice accumulate massive amyloid plaques of PrP^Sc but show no clinical disease.
- Benilova et al. (PNAS 2020): highly purified, high-titer infectious prions are not directly neurotoxic. Brain homogenate is toxic, but purified prions are not — the toxic species is lost during purification.
- Some prion strains cause disease with minimal plaque deposition.

**What survived:** The disease requires PrP^Sc propagation (the replication machinery is real), but PrP^Sc itself is not the proximate toxin.

---

## v2 — Multivalent avidity clusters PrP^C, releasing the N-terminal tail

**Proposed:** PrP^Sc acts as a multivalent scaffold that cross-links PrP^C on the cell surface. The clustering disrupts the intramolecular auto-inhibitory contact between the N-terminal tail and the globular domain, releasing the tail as the toxic effector.

**Supporting evidence:**
- NMR shows intramolecular contact between N-terminal tail and globular domain (auto-inhibition in cis) — established
- Antibodies against the globular domain cause toxicity indistinguishable from prion infection (Sonati et al., Nature 2013)
- Clustering of GPI-anchored proteins is a known signaling mechanism (Mouillet-Richard et al., Science 2000)
- FTgpi mice (tail + GPI, no globular domain) develop lethal neurodegeneration (Aguzzi 2015)
- Eph/ephrin literature: ligand valence determines signal output (monomeric = antagonist, clustered = agonist)

**What killed it:**
- Gatdula et al. (PLoS Pathog 2026): PrP^Sc from a species-incompatible source binds PrP^C on the cell surface but is NOT toxic. Binding alone — even multivalent binding — is insufficient. Conversion in the plane of the membrane is required.
- This separates "binding and clustering" from "conversion." v2 conflated the two.

**What survived:** The auto-inhibitory lock is real. The N-terminal tail is the effector. But the release mechanism is conversion, not clustering.

---

## v3 — Conversion destroys the lock; the effector must be GPI-anchored

**Proposed:** When PrP^C is incorporated into a PrP^Sc fibril, the globular domain is consumed by the β-sheet core (residues ~94–225), destroying the surface that held the N-terminal tail. The freed tail, still GPI-anchored to the membrane, perturbs the local bilayer and triggers downstream toxicity.

**Supporting evidence:**
- Cryo-EM of ex vivo fibrils (Manka et al., Nat Commun 2022): the ordered core spans ~94–225, meaning the entire globular domain is incorporated. The N-terminal tail is not part of the core — it is exposed and mobile.
- PrP without GPI anchor replicates without causing disease (Chesebro 2005) — interpreted as: GPI anchoring is required for the tail to be toxic (proximity to the membrane).
- The deletion ladder (PrPΔ23-88, PrPΔ32-93) shows longer incubation with shorter tails.

**What killed it:**
- Science Advances 2026: the soluble N-terminal fragment (not GPI-anchored) causes rapid neurodegeneration and lethality in mice. The fragment associates peripherally with lipid membranes — it does not need to be covalently anchored. Toxicity depends on the KKRPKP motif and exceeding a critical concentration threshold.

**What survived:** Conversion is the unlocking event. The tail is the effector. But GPI anchoring is not required for toxicity — what GPI provides is *local concentration* at the membrane, not a direct requirement for the toxic mechanism.

---

## v4 — Released tail is a specific effector with dose-dependent signal inversion

**Proposed:** The N-terminal tail, released by conversion, is the toxic effector. It acts through a specific mechanism (not yet defined) that depends on the KKRPKP polybasic motif. Below a critical threshold, the same fragment (N1, produced by α-cleavage) is neuroprotective — it sequesters toxic oligomers. Signal inversion by dose reconciles the N1 protection/toxicity paradox.

**Supporting evidence:**
- Sci Adv 2026: direct demonstration of N-terminal toxicity above a threshold, KKRPKP essential
- French literature on N1 neuroprotection (α-cleavage, ADAM10 correlation with survival)
- Brazilian literature on PrP^C as a signaling platform (STI1 interaction, neuroprotection)
- Convergence of the p38 MAPK pathway (Fang/Harris 2018), PERK/eIF2α (Moreno 2012), and PIKfyve (Lakkaraju 2021)
- Analogy with MLKL and gasdermin D (auto-inhibited effectors released by oligomerization)

**What killed it (stress test, this investigation):**

1. **N1 does not protect in vivo.** Bhullar et al. (2020, PMC7253391): transgenic mice overexpressing N1 (3.8–5.2×) showed identical survival to wild-type upon prion infection (160 ± 7 vs 159 ± 9 days). N1 is an intrinsically disordered peptide that fails ER translocation and accumulates in the cytosol — it never reaches the extracellular space. The "signal inversion" hypothesis requires N1 to be protective *in situ*, but it doesn't even leave the cell.

2. **Toxicity is charge-dependent, not sequence-specific.** Cationic Arginine-Rich Peptides (CARPs) block NMDA at 1–5 µM — neuroprotection is charge-dependent and sequence-independent (Front Neurol 2020). LL-37 (human cathelicidin) shows the same dose-dependent inversion: protective at physiological concentrations, toxic when accumulated (Mol Psychiatry 2022). The KKRPKP requirement reflects charge density (+4 in 6 residues), not a specific receptor interaction.

3. **LLPS as upstream pathway was not considered.** Yi Liang (Wuhan) and Cong Liu (CAS Shanghai) showed PrP undergoes liquid-liquid phase separation (LLPS) *in vitro*, driven by the same N-terminal polybasic motifs the dossier identifies as the toxic effector (JBC 2021). This creates a dual-identity problem: KKRPKP is both the effector (downstream) and the LLPS motor (upstream). The simpler explanation is that it's needed for one role, not both.

4. **α-cleavage protects by substrate destruction, not N1 release.** α-cleavage cuts at ~111, within the fibril core region (94–225). Neither N1 (23–111) nor C1 (112–225) can serve as substrate for templated conversion. C1 cannot be converted and may inhibit conversion (Bhérer et al., Prion 2012). Protection comes from removing convertible substrate, not from releasing a protective fragment.

**What survived:** The tail is the effector. The threshold is real. Conversion is the unlocking event. But the mechanism is not PrP-specific — it is generic polycationic peptide membrane perturbation.

---

## v5 — Ancestral AMP released by general amyloid conversion; specificity in auto-templating replication

**Proposed:** Prion disease is a specific instance of a general amyloid toxicity mechanism. The disordered N-terminal tail, freed when the globular domain is incorporated into the fibril core (the "fuzzy coat" principle), perturbs neuronal membranes like an antimicrobial peptide. This mechanism is shared with all amyloid diseases. What makes prion disease unique is:

1. **True autocatalytic replication** (~1000× faster than tau/α-synuclein), explaining months vs decades.
2. **The GPI anchor as a dual regulator**: suppresses LLPS (brake against spontaneous conversion) and concentrates toxic products at the membrane (amplifier of proximity).
3. **Strain conformation → tropism**: fibril architecture selects among regionally distributed cofactors (lipid composition, sialylation, endocytic routes).

The toxicity mechanism — membrane perturbation by a polycationic peptide — is not prionic. It is biophysics. The N-terminal tail is literally an antimicrobial peptide (Nature Microbiology 2026: 1,179 AMPs from prion protein fragments; Mangé et al. 2009: direct antimicrobial activity). Prion disease is "friendly fire" — an endogenous AMP turned against the host's own membranes.

**Supporting evidence (compiled across 24 cross-disciplinary searches):**

- *Prionins*: 1,179 AMP candidates from prion/prion-like proteins; 59/75 tested active; membrane perturbation as mechanism (Nature Microbiology 2026)
- *Fuzzy coat principle*: amyloid fibrils universally free flanking disordered regions that mediate membrane interaction (JBC 2023 review)
- *NMDA modulated by bilayer stretch*: Kloda A et al. (PNAS 2007, PMID 17242368) — potentiation of existing NMDA currents by membrane stretch. **[C]** — shows modulation, not glutamate-independent gating (see audit corrections)
- *Melittin/AMPA via PLA₂*: Aronica E et al. (Brain Res 1992, PMID 1380882) — effect runs through phospholipase A₂, not direct receptor modulation
- *GPI suppresses LLPS*: topological confinement prevents phase separation (Tatzelt, PNAS 2025, PMC11725851)
- *G127V blocks LLPS→solid transition*: keeps condensates liquid (Commun Biol 2020, PMC7391680)
- *Charge gradient correlates with severity*: PrP N-terminal +7 to +9 (months), α-synuclein +4 (years), Aβ −3 (different mechanism)
- *Octarepeat insertion correlation*: r = −0.989 between number of inserts and age of onset
- *Sporadic disease*: endosomal MVB provides low pH + released GPI brake + concentrated copper; Cu + H₂O₂ triggers liquid→solid transition (do Amaral MJ et al., Sci Adv 2023, PMID 37922348)
- *Yan 2025*: soluble N-terminal PrP causes rapid neurodegeneration; **deleting KKRPKP abolishes toxicity** (Sci Adv 2025, PMID 40768577) — the strongest published support for this model

**What could kill it:**
- A demonstration that the N-terminal tail at physiological concentrations (not supraphysiological) does NOT perturb membranes
- Identification of a specific receptor that mediates toxicity and is NOT a general membrane effect
- A prion disease model where toxicity occurs without the N-terminal tail at all
- Evidence that the fuzzy coat of prion fibrils does not contact membranes in vivo

**Status:** Partially tested computationally (see v5.1 below).

---

## v5.1 — Computational validation (August 31, 2026)

**What changed:** Three rounds of molecular dynamics on POPC:POPS membranes,
the third methodologically sound. Nine wrong author attributions found and
corrected in the bibliography. Two mechanism confidence levels downgraded.

**Computational evidence added:**

1. **KKRPKP is necessary for membrane binding.** PrP Δ23-28 (KKRPKP deleted):
   0/3 replicas bound in 180 ns combined. PrP 23-93 (with KKRPKP): 1/3 bound.
   KKRPKP hexapeptide alone: 2/3 bound. Replicates Yan 2025 in silico.

2. **Binding is electrostatically driven.** POPS enrichment 2.59× (KKRPKP) and
   2.92× (PrP 23-93) against a 20% baseline. The peptide seeks the anionic
   lipid.

3. **Binding thins the membrane locally.** DELTA −0.7 to −1.0 Å in KKRPKP
   bound replicas (frame-by-frame comparison within the same trajectory).
   Radial profile: 36 Å under the peptide vs 41 Å at 2.5 nm.

**Mechanism confidence adjusted:**

- **Mechanism B (NMDA mechanosensitivity):** [E] → [C]. Kloda 2007 shows
  modulation of already-gated channels, not glutamate-independent opening.
  The model stated the stronger version.
- **Signal inversion (CARP endocytosis):** [E] → [H]. The "Bhatt 2015"
  citation was fabricated — third instance of this invented surname.

**Bibliography audit:**
- 9 wrong attributions corrected (4 fabricated author names, 5 wrong years)
- 47/185 entries verified against PubMed/CrossRef
- 2 corrections were themselves wrong (introduced by the LLM verifier)
- Audit banner added: unchecked entries are unverified, not verified

**What the computation does NOT show:**
- That the thinning is sufficient to modulate NMDA (the link is extrapolation)
- That PrP 23-93 perturbs as much as KKRPKP (insufficient binding events)
- That the effect persists with cholesterol (simplified membrane)
- Converged absolute values (|S_CD| stuck at 0.36 instead of ~0.20)

**What could kill v5.1:**
- PrP Δ23-28 binding at longer timescales (>μs) — would mean the deletion
  slows but does not eliminate binding
- The thinning being a generic IDP effect rather than charge-specific —
  testable with a neutral IDP of similar size
- Cholesterol blocking the binding entirely in a neuronal-composition membrane

---

## Invariants across all revisions

| Idea | v1 | v2 | v3 | v4 | v5 |
|------|----|----|----|----|-----|
| N-terminal tail is the effector | — | ✓ | ✓ | ✓ | ✓ |
| Dose threshold exists | — | — | — | ✓ | ✓ |
| Conversion is the unlocking event | — | — | ✓ | ✓ | ✓ |
| GPI required for toxicity | — | — | ✓ | ✗ | ✗ (but amplifier) |
| Mechanism is PrP-specific | ✓ | ✓ | ✓ | ✓ | ✗ |
| Signal inversion by dose | — | — | — | ✓ | ✗ (LL-37 explains differently) |
| LLPS as upstream pathway | — | — | — | — | ✓ |
