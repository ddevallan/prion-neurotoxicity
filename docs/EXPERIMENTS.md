# Proposed Experiments

Experiments prioritized by their ability to confirm, refine, or refute the v5 model of prion neurotoxicity. Each entry states the design, the v5 prediction, and the difficulty.

---

## Priority 1 — Decide the model

These experiments test the core claims of v5. A negative result on any of E1–E4 would require significant model revision.

### E1. Memantine in RML-infected mice (survival)

**Design.** Administer memantine (oral, 10–30 mg/kg/day) to C57BL/6 mice inoculated intracerebrally with RML prions. Begin dosing at three timepoints: day 30 (pre-symptomatic), day 80 (early symptomatic), day 110 (late symptomatic). Monitor survival, clinical scoring, and neuropathology (spongiosis, synapse density, dendritic beading).

**v5 prediction.** 15–30% survival extension when started pre-symptomatically; diminishing effect with later start; improvement in synapse preservation at all timepoints. Falls within the 4–19% range Minikel documents for downstream interventions if started late.

**Rationale.** Memantine blocks the NMDA channel regardless of the activation stimulus — including mechanosensitive opening from AMP-driven membrane perturbation. Müller et al. (Eur J Pharmacol 1993) showed memantine protects neurons from PrPSc in vitro. The Prnp G92N knockin (JCI 2025) showed memantine rescues dendritic beading and extends survival. **No one has tested memantine in prion-infected mice for survival — a 33-year gap.**

**Difficulty.** Trivial. FDA-approved drug, standard RML mouse model, standard endpoints.

---

### E2. Quantify free N-terminal fragment throughout infection

**Design.** Collect serial brain homogenates from RML-infected mice at 7-day intervals from inoculation to terminal stage. Quantify the free N-terminal fragment (residues 23–90) by sandwich ELISA using an N-terminal-specific capture antibody (e.g., anti-PrP 23–30) and a detection antibody against octarepeats. Normalize to total PrPC. Run in parallel on mice expressing different levels of PrPC (wild-type, hemizygous, Tga20 overexpressor).

**v5 prediction.** The free N-terminal fragment rises at the phase 1 → phase 2 transition (Sandberg 2011), not earlier. The accumulation rate is proportional to PrPC concentration. The fragment coincides temporally with the protease-sensitive species described by Sandberg et al. (Nat Commun 2014).

**Rationale.** Tests the identification PrPL = released N-terminal tails (currently at 40% confidence, 5 competing candidates). If confirmed, directly validates the v5 causal chain: conversion → tail release → toxicity.

**Difficulty.** Easy. Requires only an N-terminal antibody and serial material. Both exist in multiple labs.

---

### E3. Low-MW heparan sulfate analog intrathecal in infected mice

**Design.** Intrathecal infusion of a low-molecular-weight 2-O-sulfated heparan sulfate analog (e.g., HS mimetic with high affinity for KKRPKP, <5 kDa for CNS distribution) in RML-infected mice. Start at day 80. Compare with intrathecal vehicle control and systemic heparan sulfate (which should not cross BBB). Add an arm combining the HS analog with ION717-like ASO (if available) to test synergy.

**v5 prediction.** Partial survival extension from the HS analog alone (neutralizes the freed N-terminal AMP). The effect should be additive with ASO (different mechanisms: HS neutralizes existing AMP, ASO prevents new AMP generation). Retroactively reinterprets pentosan polysulfate (PPS) — which showed preclinical efficacy but was attributed to anti-conversion activity; PPS inhibits toxicity in cells but *stimulates* conversion cell-free (CureFFI analysis), consistent with AMP neutralization rather than conversion blockade.

**Rationale.** Heparan sulfate binds PrP N-terminal at three sites (residues 23–52, 53–93, 110–128; JBC 2002) with nanomolar affinity. 2-O-sulfated HS neutralizes LL-37 by the same electrostatic mechanism (Hayashida 2025; PMC12806164). The intrathecal route avoids PPS's delivery problems (intraventricular infusion).

**Difficulty.** Moderate. Requires HS analog chemistry (available from glycobiology labs) and intrathecal infusion setup (same as ION717/PRiSM trials).

---

### E4. Memantine + ASO combination in RML mice

**Design.** Four-arm study in RML-infected mice: (1) vehicle; (2) ASO alone (intrathecal, started day 80); (3) memantine alone (oral, started day 80); (4) ASO + memantine (both started day 80). Primary endpoint: survival. Secondary: PrP levels in CSF, neuropathology, synapse counts.

**v5 prediction.** Combination arm survives longest. Memantine alone gives +15–30%. ASO alone gives >50% extension (extrapolating from PRiSM mouse data: 64% post-symptomatic). Combination should be additive or super-additive: memantine protects synapses during the weeks the ASO needs to reduce PrPC.

**Rationale.** The v5 "bridge" strategy. No combination therapy has been tested in prion disease despite a 2026 review in *Expert Opinion on Drug Discovery* explicitly calling for multi-target approaches. Anti-PrPSc drug combinations generate strain-specific resistance (Berry et al., PLoS Pathog 2020); v5's combination targets host-encoded steps — no evolutionary escape possible.

**Difficulty.** Moderate. Two agents, standard model, 4 arms × 10 mice = 40 animals minimum.

---

## Priority 2 — Mechanism

These refine v5's mechanistic claims. They don't decide the model but sharpen its specificity.

### E5. N1 fragment in RT-QuIC / PMCA

**Design.** Use recombinant N1 fragment (PrP 23–111) as substrate in RT-QuIC and PMCA, seeded with RML brain homogenate. Test whether N1 supports conversion. Compare with full-length PrPC and C1 fragment (PrP 112–231) as substrates.

**v5 prediction.** N1 should NOT form infectious fibrils. N1 may undergo LLPS (the N-terminal domain is sufficient for phase separation; Tatzelt PNAS 2025) but the product should not seed prion conversion, because α-cleavage cuts within the fibril core (~residue 111) and neither N1 nor C1 retains the full core sequence (94–225). This reframes α-cleavage as substrate destruction.

**Rationale.** Bhullar et al. (2020, PMC7253391) showed TgN1 mice have no survival benefit — N1 fails ER translocation and never reaches the extracellular space. The protection attributed to N1 should instead be attributed to substrate depletion (producing non-convertible C1).

**Difficulty.** Easy. RT-QuIC is routine; recombinant fragments are standard.

---

### E6. Endogenous AMP profiling throughout prion infection

**Design.** Transcriptomic (RNA-seq) and proteomic profiling of endogenous antimicrobial peptides — including CRAMP (mouse cathelicidin, homolog of human LL-37), α/β-defensins, and other innate immune effectors — in RML-infected mouse brain at serial timepoints. Correlate with microglial activation markers and clinical stage.

**v5 prediction.** CRAMP and defensins are upregulated during prion neuroinflammation, creating a positive feedback loop: prion conversion → N-terminal AMP release → membrane perturbation → neuroinflammation → microglial activation → endogenous AMP upregulation → additional membrane perturbation. If true, anti-inflammatory interventions that suppress AMP production would have a secondary neuroprotective effect beyond their anti-inflammatory action.

**Rationale.** No study has measured AMP levels in prion-infected brains — a complete blind spot. LL-37 is upregulated in AD brains (Bhatt et al., Biochem Pharmacol 2015) and causes AD phenotypes in mice and primates (Mol Psychiatry 2022). CRAMP is induced in astrocytes, microglia, and neurons during neuroinflammation (JCI 2025, PMC11785927). Prion neuroinflammation uses the same pathways (NALP3 inflammasome, TNF-α, IL-6) known to upregulate AMPs.

**Difficulty.** Easy — transcriptomic data on prion-infected brains likely already exists; re-analysis for AMP genes may be sufficient.

---

### E7. Synthetic charge-gradient peptides in dendritic spine assay

**Design.** Synthesize a series of 6-residue peptides with increasing net charge: NNRPNP (+0), KNRPNP (+1), KNRPKP (+2), KKRPNP (+3), KKRPKP (+4). Apply to hippocampal neuron cultures at sub-lytic concentrations. Measure dendritic spine density, NMDA-dependent Ca²⁺ influx, and p38 MAPK activation. Include LL-37 and melittin as positive controls.

**v5 prediction.** Dose-response scales with charge, not sequence. +4 peptides cause spine retraction similar to PrPSc. +0 peptides are inert. The charge threshold for toxicity is approximately +2 to +3, matching the AMP literature. Memantine co-application rescues the phenotype.

**Rationale.** CARP literature shows NMDA blockade is charge-dependent and sequence-independent (Front Neurol 2020). If synthetic peptides reproduce PrPSc-like synaptotoxicity purely by charge, the "PrP-specific mechanism" claim is definitively refuted — it's generic polycationic membrane perturbation.

**Difficulty.** Moderate. Peptide synthesis is trivial; the dendritic spine assay (Fang et al., PLoS Pathog 2016) requires electrophysiology expertise.

---

### E8. PrP Csat as function of age, PTMs, and lipid composition

**Design.** Measure the saturation concentration (Csat) for PrP LLPS in vitro under varying conditions: (1) PrP age (fresh vs oxidatively damaged); (2) with/without Cu²⁺ at physiological ratios; (3) Cu²⁺ + H₂O₂; (4) in the presence of lipid vesicles of different compositions (POPC, POPC/cholesterol, brain-region-specific lipid extracts). Use turbidity and fluorescence microscopy to detect phase separation.

**v5 prediction.** Csat decreases with oxidative damage, Cu²⁺ + H₂O₂, and aging-related PTMs. LLPS is suppressed by GPI-anchored membrane context (Tatzelt PNAS 2025) but enhanced in acidic, concentrated compartments (MVB-like conditions). Brain-region-specific lipids may differentially modulate Csat, contributing to strain tropism.

**Rationale.** Tests the sporadic disease mechanism: age-related Csat decline → spontaneous LLPS → conversion. Evans et al. (Sci Adv 2023) showed Cu + H₂O₂ triggers the liquid-to-solid transition. No one has measured Csat as a function of age or PTM state.

**Difficulty.** Moderate. Requires recombinant PrP, DLS/turbidity measurements, fluorescence microscopy. All standard biophysics.

---

## Priority 3 — Computational (in progress)

### E9. MD simulation: KKRPKP vs NNRPNP on POPC membrane

**Status.** In progress. OpenMM 8.6 installed with OpenCL acceleration on Apple M1 Pro. Peptide-in-water equilibrations complete (100 ps). Extended 10 ns simulation running. Full membrane simulation planned on Vast.ai (A40 GPU, ~$10 for complete study).

**Result so far.** KKRPKP (charge +4) has 56% more solvent-accessible surface area than NNRPNP (charge 0) at the same radius of gyration (1.79 nm). Charged sidechains project outward — the conformation an AMP needs for carpet-model membrane insertion.

---

### E10. Cross-species PrP N-terminal charge analysis

**Status.** Complete.

**Result.** Net charge is +6.5 across ALL mammals tested (human, mouse, hamster, sheep, cow, elk, bank vole, dog, rabbit, horse) — susceptible and resistant species alike. The N-terminal sequence is highly conserved. Resistance (dog, rabbit, horse) does NOT come from charge differences in the tail but from the globular domain and hydrophobic gatekeeper. This confirms v5: toxicity mechanism is generic (same AMP), specificity is in conversion.

---

### E11. Kinetic model of therapeutic window

**Status.** Complete. Calibrated to RML in wild-type mouse (~148 day baseline).

**Key results.**
- Memantine alone (day 80): +28 days (+19%) — matches Minikel's 4–19% for downstream interventions.
- ASO alone (day 110): survival.
- Bridge strategy (memantine day 90 + ASO day 104): survival — memantine buys time while ASO ramps.
- Full combination (day 110): survival.
- Late rescue (everything day 120): survival.

---

### E12. Mutation mapping to v5 functional domains

**Status.** Complete.

**Key results.**
- Octarepeat insertions: r = −0.989 correlation with onset age (1-OPRI → 62 yr, 9-OPRI → 20 yr). Each insert adds ~2 charges → more potent AMP + more LLPS.
- Hydrophobic gatekeeper mutations (A117V, G131V) → GSS (slow, chronic) — consistent with partial gating.
- Globular domain mutations → widest phenotype range (CJD rapid to FFI to GSS slow).
- Both protective mutations (G127V, E219K) are in the gatekeeper or globular domain — none in the N-terminal tail.
- D178N + M129 → FFI; D178N + V129 → CJD — the gatekeeper modulates a globular mutation's phenotype.
