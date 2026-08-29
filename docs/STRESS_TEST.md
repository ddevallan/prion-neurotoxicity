# Stress Test: Attacking the v4 Model

This document records the eight attacks formulated against the v4 model of prion neurotoxicity, plus one emergent challenge (LLPS) discovered during cross-disciplinary search. Each entry states the attack, the evidence found, the resolution or current status, and the impact on the model's evolution to v5.

---

## Attack 1 — The deletion ladder confound is worse than acknowledged

**The attack:** The v4 dossier cites the deletion ladder (PrPΔ23-88 → 161 days, PrPΔ32-93 → 232–313 days, wild-type → 50 days) as evidence that shorter tail = less toxicity. But it acknowledges in passing that the deletions also impair conversion. This confound is not minor — it is fatal to the interpretation. If less conversion means less PrP^Sc, the slower disease could reflect reduced replication, not reduced AMP potency. The ladder may be measuring conversion efficiency, not effector dose.

**Evidence found:**
- Minikel (CureFFI 2015): PrPΔ23-88 mice still die of prion disease at 161 days. If the tail were THE sole effector, deletion should confer much stronger protection. Downstream pathway interventions give only 4–19% survival extension vs 2–4× for anti-propagation strategies.
- The MHM2Δ23-88 / Prnp^0/0 mouse is resistant to 730 days — but this is a different species chimera on a knockout background, confounding species barrier with tail deletion.

**Resolution:** Partially resolved. The confound remains — the ladder cannot cleanly distinguish "less AMP" from "less conversion." However, the octarepeat insertion data (see mutation analysis) provides a cleaner test: more insertions → earlier onset with r = −0.989, and insertions ADD charge without changing the globular domain. Under v5, the ladder confound is expected: the N-terminal drives both LLPS (upstream conversion) and AMP toxicity (downstream), so deleting it impairs both.

**Status:** Open confound, but less damaging under v5 than under v4.

---

## Attack 2 — α-cleavage protects by substrate destruction, not N1 release

**The attack:** The v4 model assumes N1 (released by α-cleavage) is neuroprotective at low doses. But α-cleavage produces two fragments: N1 (23–111, soluble) and C1 (112–225, GPI-anchored). If C1 cannot be converted to PrP^Sc (it lacks the N-terminal portion of the fibril core), then α-cleavage protects by *removing convertible substrate from the membrane*, not by releasing a protective N1 fragment. This is simpler and doesn't require signal inversion.

**Evidence found:**
- **Bhullar et al. 2020 (PMC7253391):** TgN1 mice overexpressing N1 at 3.8–5.2× showed NO protection against prion disease. Survival: 160 ± 7 vs 159 ± 9 days (identical to wild-type). N1 fails ER translocation, retains its signal peptide, and accumulates in the cytosol. It never reaches the extracellular space where it would need to act.
- **C1 cannot convert:** C1 lacks the N-terminal portion of the fibril core (the ordered region starts at ~94). Bhérer et al. (Prion 2012) reported C1 may actively inhibit prion propagation.
- **C1 is abundant:** 10–50% of total PrP^C in normal brain is C1.

**Resolution:** Resolved. The "N1 paradox" dissolves: N1 is not protective in vivo (Bhullar 2020 proved this), and α-cleavage protects by substrate destruction. The signal-inversion hypothesis is unnecessary. Under v5, LL-37's dose-dependent effects (protective at low dose, toxic at high dose) provide an alternative framework that doesn't require N1 to be the protective species — the mechanism is generic cationic-peptide biphasic response (low dose → NMDA endocytosis; high dose → membrane perturbation dominates).

**Impact:** Fatal to v4's reconciliation of the N1 paradox. One of the cleanest kills.

---

## Attack 3 — PrP^L at 40% confidence is the load-bearing assumption

**The attack:** The v4 dossier hypothesizes that PrP^L (the protease-sensitive species Sandberg described rising during phase 2) is accumulated released N-terminal tails. This is assigned 40% confidence but is structurally load-bearing — without it, the model has no temporal mechanism linking Sandberg's kinetics to neurodegeneration.

**Evidence found:**
- **Sandberg/Collinge 2014:** PrP^L was deliberately left undefined. The authors called it "PK-sensitive disease-related PrP" and said it could be "a single defined species of high specific toxicity, or an ensemble of diverse species." No structural characterization, no mass spectrometry, no specific antibody.
- **Five competing candidates:** (1) β-rich oligomers / PrP* (~20 monomers, soluble); (2) sarkosyl-sensitive non-prion assemblies (Benilova 2020: sarkosyl destroys toxicity without reducing infectivity); (3) syntaxin-6-stabilized on-pathway intermediates (eLife 2024); (4) PK-resistant multimers (contradicting Sandberg); (5) released N-terminal tails (the v4 hypothesis, zero direct evidence).
- **Benilova 2020:** The strongest constraint — the toxic species is sarkosyl-sensitive. This is consistent with loose membrane-associated fragments but equally consistent with β-oligomers or transient intermediates. Does not discriminate.

**Resolution:** Partially resolved. Under v5, the identity of PrP^L is less critical because the model doesn't depend on a single specific toxic species. The toxicity comes from the generic AMP activity of the freed N-terminal, and PrP^L could be any ensemble that includes freed tails, LLPS intermediates, or both. The temporal link (why phase 2 correlates with disease) is maintained: conversion at steady state keeps generating subproducts at a rate proportional to PrP^C.

**Status:** PrP^L identity remains unresolved. Experiment E2 (dose the free N-terminal fragment along the infection) would close this gap.

---

## Attack 4 — The toxicity may be nonspecific polycationic membrane disruption

**The attack:** KKRPKP is a polybasic sequence. If the N-terminal fragment was tested at supraphysiological levels, the toxicity could be nonspecific — the same thing any cationic peptide does at high enough concentration. The v4 dossier notes the MLKL "antimicrobial peptide-like mechanism" analogy as if it were confirmation, but it is actually the most dangerous attack: if it is *literally* an antimicrobial peptide, there is no prion-specific toxicity mechanism.

**Evidence found:**
- **Mangé et al. (PLoS One 2009):** PrP N-terminal peptides containing KKRPKP function as cell-penetrating peptides, cause membrane perturbation similar to melittin and LL-37, and kill bacteria nonspecifically. The antimicrobial activity is mediated entirely by the N-terminal region.
- **De la Fuente Lab (Nature Microbiology 2026, "prionins"):** Deep learning scan of 19.3 million fragments from 2,897 prion and prion-like proteins found 1,179 candidate AMPs. 59/75 synthesized inhibited bacterial pathogens. 53 perturbed membranes. Many adopted ordered conformations only upon membrane contact. This confirms the identity: PrP N-terminal fragments ARE functional AMPs by definition.
- **CARPs (Front Neurol 2020):** Poly-arginine R18 at 1–5 µM reduces neuronal death after NMDA exposure. Neuroprotection is charge-dependent and sequence-independent. Any sufficiently cationic peptide blocks NMDA.

**Resolution:** Confirmed as correct — and this confirmation *strengthened* the model by transforming it into v5. The toxicity IS nonspecific polycationic membrane perturbation. This is not a weakness but the core insight of v5: the specificity of prion disease is in the *replication* (auto-templating, GPI concentration, strain architecture), not in the *toxicity* (which is generic AMP biophysics). The attack killed v4's claim of specificity but birthed v5.

**Impact:** Transformed the model. The most productive attack.

---

## Attack 5 — The specification sheet is incomplete

**The attack:** The v4 specification sheet (R1–R10) omits at least three constraints any valid model must satisfy: (R11) strain tropism — different strains target different brain regions in the same host; (R12) glycosylation ratio varies by strain and correlates with clinical phenotype; (R13) astrocytes propagate prions without degenerating.

**Evidence found:**
- **Tuzi et al. (J Virol 2008):** PrP^Sc containing only unglycosylated form maintained strain-specific patterns. Glycosylation ratio is a readout, not a driver of tropism. R12 refined.
- **Makarava et al. (IJMS 2020):** Sialylation is region-dependent and modulates replication rate. Lower sialylation = faster replication = greater vulnerability.
- **Barria et al. (JBC 2018):** Region-specific PMCA reproduces in vivo tropism. Correcting for PrP^C levels reduces but does not eliminate the pattern (r = 0.65 → 0.40). Non-PrP cofactors contribute.
- **Lipid raft atlas (Nat Commun 2024):** 84% of 419 lipid species correlate with cell-type markers. Sphingomyelin depletion → 4× more PrP^Sc. Cholesterol does NOT protect against AMP insertion at raft boundaries (PubMed 22885355).
- **Bhérer et al. (Biology 2024):** Strains differ in neuron-vs-astrocyte tropism. Astrocytes propagate without degenerating.

**Resolution:** Resolved. v5 explains tropism as a three-way interaction: fibril conformation (strain) × regionally distributed cofactors (lipid composition, sialylation) × local PrP^C concentration and endocytic routes. The toxicity mechanism is the same everywhere — what varies is WHERE and HOW FAST conversion happens. Astrocytes may resist because their membrane composition is less susceptible to AMP perturbation (lower sphingolipid, different raft structure).

**Status:** R11–R13 added to the v5 specification sheet. Astrocyte resistance mechanism is a hypothesis, not demonstrated.

---

## Attack 6 — The MLKL analogy is dangerously general

**The attack:** In MLKL, auto-inhibition is released by phosphorylation by RIPK3. In PrP, by fibril incorporation. These are completely different mechanisms. "Protein with an auto-inhibited domain" describes thousands of proteins. The analogy gives false confidence of understanding when what exists is a metaphor.

**Evidence found:**
- **NINJ1 (Nature 2024):** A much closer structural analogy. Auto-inhibited homodimer where the hydrophilic membrane-rupturing face is buried. Re-oligomerizes into filaments exposing the effector face → membrane rupture. Release mechanism is change of oligomeric state (dissociation → re-oligomerization into a different assembly), not phosphorylation or cleavage. Directly cytotoxic.
- **Fuzzy coat principle (JBC 2023 review):** Amyloid fibrils *universally* consist of an ordered core flanked by disordered regions that mediate membrane binding and disruption. This is the exact structural logic of the PrP model, generalized. The core locks into place, freeing the flanks for new interactions.
- **Geometric logic:** The intramolecular contact surface overlaps with the intermolecular oligomerization interface. In PrP, cryo-EM shows residues 94–225 forming the core with every surface engaged — the surface that contacted the N-terminal tail is physically unavailable in the fibril.

**Resolution:** Refined. MLKL remains a valid analogy but NINJ1 is structurally closer. The fuzzy coat principle generalizes the mechanism to all amyloids: fibril core formation frees flanking IDRs for membrane interactions. The analogy is not dangerously general — it identifies a *real structural principle* common to NINJ1, MLKL, gasdermin D, and amyloid fuzzy coats.

**Status:** Resolved. The NINJ1 analogy and fuzzy coat principle replace MLKL as the primary structural parallel.

---

## Attack 7 — The model does not explain sporadic disease

**The attack:** 85% of prion disease is sporadic. If conversion triggers the mechanism, what causes the first conversion? If spontaneous misfolding is rare enough to explain ~1/million/year incidence, the initial concentration of released tails would be infinitesimal — far below any threshold.

**Evidence found:**
- **Endosomal MVB as nucleation site (J Cell Sci 2015):** Late endosomes/multivesicular bodies are the major intracellular site of prion conversion. pH 4.5–5.5 destabilizes PrP. On intraluminal vesicles, PrP's ectodomain is exposed to the endosomal lumen — the GPI topological brake is released. Concentration is massively increased by geometric compression.
- **Copper + oxidative stress (Evans et al., Sci Adv 2023):** Copper drives PrP phase separation at physiological ratios. Copper alone keeps condensates liquid (proposed function: copper buffering at synapses). Copper + H₂O₂ triggers the liquid-to-solid transition → amyloid with dityrosine cross-links.
- **GPI anchor suppresses LLPS (Tatzelt, PNAS 2025):** Membrane confinement prevents LLPS. When released from membrane, PrP spontaneously transitions to aggregates.
- **Lipid raft degradation with aging:** Progressive reduction in sphingolipids and cholesterol weakens the GPI brake.
- **No somatic mutations (2024, 205 sCJD cases):** Deep sequencing found no enrichment of somatic PRNP mutations vs controls.

**Resolution:** Resolved. Sporadic disease arises from the convergence of three age-related vulnerabilities: (1) endosomal transit provides low pH, high concentration, and weakened GPI constraint (constitutive low-level risk); (2) copper/ROS dysregulation converts protective liquid condensates to pathological seeds; (3) lipid raft degradation weakens the GPI brake. The ~1/million/year incidence reflects the probability that all three factors align in the same cell, at the same moment. No single factor is sufficient; the convergence is the rare event.

**Status:** Resolved as a mechanistic framework. Untested prediction: PrP Csat should decrease with age, PTM state, and membrane composition (experiment E8).

---

## Attack 8 — Gain of function vs loss of function

**The attack:** The model is entirely gain-of-function (released tail is toxic). But PrP^C has a function (whatever it may be). Depletion of PrP^C reverses early spongiosis (Mallucci 2003) — consistent with the model, but equally consistent with loss of function (the PrP^C being sequestered into fibrils can no longer perform its normal role). The knockout being healthy does not resolve this because compensation by Doppel and Shadoo may mask the deficit.

**Evidence found:**
- **PrP^C as AMP (Mangé 2009, Nature Microbiology 2026):** PrP^C has direct antimicrobial activity, localized to the N-terminal. This IS a physiological function. Disease = the antimicrobial function turned against the host.
- **PrP^C as copper buffer (Evans, Sci Adv 2023):** Copper-PrP condensates buffer synaptic copper via LLPS. This is another physiological function. Disease disrupts it.
- **Knockout mice are healthy but show subtle phenotypes:** impaired olfaction, altered circadian rhythms, mildly increased susceptibility to seizures. Compensated but not silent.

**Resolution:** Partially resolved. Under v5, the gain-of-function (AMP release) and loss-of-function (loss of copper buffering, loss of antimicrobial defense) operate simultaneously. The gain-of-function dominates clinically (it's what kills neurons), but loss-of-function contributes to vulnerability (loss of copper homeostasis may promote further condensate solidification). The disease is both: the AMP is released (gain) because the protein that was performing a job (copper buffering, antimicrobial patrol) is consumed by fibrils (loss).

**Status:** Resolved conceptually. The "friendly fire" framing unifies gain and loss: the same function (membrane perturbation for antimicrobial defense) is the gain-of-function toxin when activated in the wrong context.

---

## Emergent Attack — LLPS as alternative upstream pathway

**Origin:** Not part of the initial eight attacks. Emerged from the Chinese prion research search (Yi Liang, Wuhan; Cong Liu, CAS Shanghai).

**The challenge:** The N-terminal polybasic domain (including KKRPKP) drives LLPS, which initiates conformational conversion. This creates a dual-identity problem: in the v4 model, KKRPKP is the toxic effector (downstream of conversion); in the LLPS literature, it is the motor of phase separation (upstream of conversion). These are opposite causal roles for the same sequence.

**Evidence found:**
- **PrP undergoes LLPS in vitro, driven by N-terminal polybasic motifs (JBC 2021, PMC8289115)**
- **LLPS → amyloid conversion is a validated general mechanism** across α-synuclein, tau, Aβ, FUS, TDP-43 (now textbook: Alberti & Dormann, Dev Cell 2020)
- **Two distinct thresholds:** one for LLPS onset (Csat), one for liquid→solid transition
- **Dog PrP Arg177/Asp159 slow LLPS and inhibit amyloid formation of human PrP (JBC 2023)** — resistance via LLPS modulation
- **G127V modulates LLPS:** promotes droplet fluidity, prevents liquid→solid transition (J Neurochem 2023; Commun Biol 2020)
- **GPI anchor suppresses LLPS (Tatzelt, PNAS 2025):** topological confinement prevents phase separation

**Resolution:** Integrated into v5. The N-terminal does BOTH — it drives LLPS upstream and acts as AMP downstream. This is not a contradiction but a dual role: the same sequence that causes the protein to condense (polybasic motifs driving multivalent interactions) is the same sequence that perturbs membranes when freed (polybasic motifs inserting into anionic membrane surface). The LLPS pathway adds a mechanistic layer that v4 lacked: it explains how conversion begins (phase separation → nucleation → templating), how G127V protects (keeps condensates liquid), and how sporadic disease initiates (age-related Csat decline).

**Impact:** The LLPS integration was the largest structural change from v4 to v5. It added the upstream pathway, connected G127V mechanistically, and provided the sporadic disease explanation.
