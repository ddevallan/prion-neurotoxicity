# Model v5: Neurotoxicity in Prion Disease

> Checkpoint — August 28, 2026. Compiled from exploratory investigation with 24 cross-disciplinary searches across prion biology, antimicrobial peptides, LLPS/condensate physics, membrane biophysics, necroptosis/pyroptosis effectors, and clinical therapeutics. Not peer-reviewed.

## 1. The Question

Formulated as an equation:

**ALTERED PRION PROTEIN + X = DAMAGED NEURON**

The goal is to identify **X** by reverse engineering: start from the observed damage and reconstruct the causal chain back to the molecular property that initiates it. This is an ill-posed inverse problem — many different causes feed into the same death pathways (proteostasis failure, unfolded protein response, oxidative stress, synaptic dysfunction), all shared with Alzheimer's and Parkinson's.

---

## 2. The v5 Model

### Core statement

Prion disease is a specific instance of a **general amyloid toxicity mechanism** — the disordered "fuzzy coat" released by fibril formation perturbs membranes as an antimicrobial peptide (AMP). The PrP N-terminal tail (residues 23–93, net charge +7 to +9) is literally an AMP: fragments from prion and prion-like proteins were confirmed as functional antimicrobial peptides by computational screening and synthesis (Nature Microbiology 2026, "prionins" — 59/75 tested inhibited pathogens, mechanism = membrane disruption). **[E]**

### What makes prion disease specific

The specificity of prion disease resides exclusively in three properties that no other amyloidosis shares at the same intensity:

1. **Truly autocatalytic self-templating replication** — prion doubling time is ~2–5 days in vivo, roughly 1,000× faster than tau (~5 years in human brain) or α-synuclein. This exponential fragmentation-elongation mechanism, not mere seeded aggregation, explains why prion disease kills in months while Alzheimer's and Parkinson's take decades. **[E]**
2. **The GPI anchor as a dual-function regulator** — it suppresses LLPS (brake against spontaneous conversion) and simultaneously concentrates conversion products at the neuronal membrane (amplifier of proximity). No other amyloidogenic protein has this built-in membrane localization coupled with LLPS suppression. **[E]** for each role individually; **[H]** for the dual-function framing.
3. **Strain conformation determines tropism** — the fibril architecture selects among regionally distributed cofactors (lipid composition, sialylation patterns, endocytic routes), producing the strain-specific neuropathology that is unique to prion disease. **[E]**

### What is NOT specific to prions

- The toxicity mechanism (AMP-like membrane perturbation) — shared with all amyloid fuzzy coats. **[H]**
- The dose-dependent signal inversion — documented for LL-37 and other cationic peptides. **[E]** for LL-37; **[H]** for PrP.
- LLPS as a nucleation pathway — confirmed for α-synuclein, tau, Aβ, FUS, TDP-43. **[E]**
- The downstream cascades (NMDA → Ca²⁺ → p38 MAPK; PERK → eIF2α) — common neurodegeneration pathways.

---

## 3. The Causal Cascade

### Initiation

**PrP^C is stable on the membrane** → GPI anchor suppresses LLPS through topological confinement (Tatzelt Lab, PNAS 2025). **[E]**

↓

**Initiation event**: external seed (infection, iatrogenic) **OR** spontaneous LLPS in the late endosome / multivesicular body (MVB), where: pH 4.5–5.5 destabilizes PrP; intraluminal vesicle topology releases the GPI brake; geometric compression concentrates protein; Cu²⁺ + H₂O₂ triggers the liquid-to-solid transition via dityrosine cross-links (Evans et al., Sci Adv 2023). **[E]** for endosomal conversion; **[H]** for the three-factor convergence as the sporadic initiator.

### Conversion

↓

**N-terminal domain drives LLPS** (polybasic motifs including KKRPKP). The same region is "required and sufficient" for PrP phase separation. **[E]**

↓

**Hydrophobic gatekeeper domain (112–133) controls the liquid → solid transition**. G127V keeps condensates liquid and prevents fibril formation — it is a natural "gate lock." α-cleavage at ~111 cuts within this region, destroying convertible substrate. **[E]** for G127V mechanism; **[E]** for α-cleavage site.

↓

**Fibril formation** incorporates the globular domain + hydrophobic domain into the ordered core (~94–225, PIRIBS architecture). The intramolecular contact surface between the N-terminal tail and the globular domain is now occupied by inter-rung β-sheet contacts. **[E]**

↓

**N-terminal tail is released** as the fibril's "fuzzy coat" — exposed, flexible, soluble, and positioned in situ at the membrane. **[E]**

### Toxicity

↓

**Sub-lytic AMP-like membrane perturbation**. The freed tail (charge +7 to +9) inserts into the outer leaflet of the lipid bilayer via the carpet/detergent model (no helix induction, charge-dependent, not sequence-specific). Insertion is preferential at lipid raft boundaries (where cholesterol's protective effect fails) and at regions of high membrane curvature (dendritic spine necks, radius 50–200 nm). **[E]** for AMP carpet mechanism in general; **[H]** for PrP N-terminal acting this way in vivo.

↓

**NMDA receptor modulation by mechanosensitivity**, not by glutamate excess:
- Membrane thinning/tension → NMDA is mechanosensitive → alleviates Mg²⁺ block → Ca²⁺ influx (Bhatt et al., PNAS 2007). **[E]**
- PLA2 activation by the membrane-inserted peptide → arachidonic acid → additional NMDA potentiation. **[E]** for the pathway; **[H]** for PrP triggering it.
- At sub-lytic doses, cationic peptides can also cause NMDA receptor endocytosis → reduced Ca²⁺ → neuroprotection (Bhatt 2015). This may explain dose-dependent signal inversion. **[E]** for CARPs; **[H]** for PrP.

↓

**Downstream cascades** (all **[E]**):
- Ca²⁺ → p38 MAPK / MK2-3 → actin collapse → spine retraction
- PERK → eIF2α-P → global translation shutdown
- zDHHC9/21 → PIKfyve → lysosomal trafficking failure → vacuoles (spongiform change)

↓

**Synaptic selectivity**: dendritic spines are preferentially damaged because they have the highest membrane curvature (most lipid packing defects → most AMP insertion) and are enriched in lipid rafts (where PrP^C resides). This predicts spine neck collapse before head loss — consistent with observed dendritic beading. **[H]**

↓

**Neuronal death**. Astrocytes propagate prions but do not degenerate — possibly because astrocytic membrane composition (different sphingolipid content, different raft architecture) is less susceptible to AMP-type perturbation. **[H]**

---

## 4. Functional Architecture by Domain

| Domain | Residues | Role in v5 | Evidence level |
|--------|----------|------------|----------------|
| **N-terminal tail** | 23–93 | Dual role: (1) drives LLPS upstream (polybasic motifs, including KKRPKP at 23–28 and octarepeat region 51–91 that binds Cu²⁺); (2) freed as the toxic AMP effector downstream. Net charge +7 to +9. Highest charge density among all amyloid flanking regions. | **[E]** for both activities independently |
| **Hydrophobic gatekeeper** | 112–133 | Controls the liquid-to-solid transition in condensates. Contains the palindrome AGAAAAGA. G127V "locks the gate open" (keeps condensates liquid, prevents fibril nucleation). α-cleavage at ~111 cuts here — fragments cannot form fibrils. Codon 129 (M/V polymorphism) is within this domain and modulates strain phenotype (D178N + M129 → FFI; D178N + V129 → CJD). | **[E]** for G127V; **[E]** for α-cleavage site; **[H]** for "gatekeeper" framing |
| **Globular domain** | ~134–225 | Structural core of the fibril (3 α-helices + 2 β-strands in PrP^C; PIRIBS β-sheet in PrP^Sc). Provides the auto-inhibitory contact surface that sequesters the N-terminal tail in monomeric PrP^C. Most pathogenic point mutations (17/30 cataloged) are in this domain, producing the widest range of phenotypes. | **[E]** |
| **GPI anchor** | C-terminal | Dual-function regulator. **Brake**: suppresses LLPS through 2D topological confinement — PrP released from membranes spontaneously transitions to insoluble aggregates (Tatzelt PNAS 2025). **Amplifier**: concentrates conversion products at the neuronal membrane, keeping AMP above the toxicity threshold locally. Without GPI: prion replicates freely (no LLPS brake) but tails are diluted in extracellular space → below threshold → no clinical disease (Chesebro 2005). | **[E]** for each role; **[H]** for the combined framing |

---

## 5. Specification Sheet

Any valid model must simultaneously satisfy these constraints:

| # | Constraint | Evidence | v5 explanation |
|---|-----------|----------|----------------|
| R1 | Requires PrP^C; knockout is resistant | Büeler 1992; Brandner 1996 | No substrate → no cascade |
| R2 | PrP^C depletion reverses early spongiform change | Mallucci 2003 | Stops AMP generation |
| R3 | PrP without GPI replicates without causing disease | Chesebro 2005 | No GPI = no local membrane concentration = below AMP threshold |
| R4 | Conversion in the membrane plane is necessary; binding alone is insufficient | Gatdula 2026 | Conversion in endosomal compartment; toxicity at membrane |
| R5 | Reproducible without PrP^Sc (anti-globular antibodies; FTgpi) | Sonati 2013; Aguzzi 2015 | Antibodies destabilize the auto-inhibitory lock → free the tail; FTgpi = constitutively active tail |
| R6 | Correlates poorly with total PrP^Sc load | Sandberg 2011 | Toxicity comes from the freed flanking region, not from the fibril itself |
| R7 | Infectivity and toxicity are decoupled | Sandberg 2011; Benilova 2020 | Infectivity = fibril replication; toxicity = membrane perturbation by freed tail. Different parameters |
| R8 | Depends on polybasic region 23–31 / KKRPKP | Solomon/Harris; Sci Adv 2026 | KKRPKP = charge center of the AMP (charge-dependent, not sequence-specific — CARPs show any sufficiently cationic peptide does the same) |
| R9 | Dose threshold exists | Sci Adv 2026 | Sub-lytic AMP threshold, documented for LL-37, melittin, CARPs |
| R10 | Partially shared with Aβ, tau, α-synuclein | Corbett 2020 — contested | General amyloid fuzzy-coat biology; Aβ binds any cationic membrane peptide (LL-37 does the same) |
| R11 | Different strains → different regional patterns in the same host | Tuzi 2008; Barria 2018 | Fibril conformation × regional cofactors (lipid raft composition, sialylation, endocytic routes) |
| R12 | Glycosylation ratio varies by strain, correlates with phenotype | Tuzi 2008 | Readout, not driver — completely deglycosylated PrP^Sc maintains strain-specific tropism |
| R13 | Astrocytes propagate without degenerating | Bhérer 2024 | Propagation ≠ toxicity; astrocyte membrane composition may resist AMP insertion |

---

## 6. Cross-Disciplinary Evidence

### 6.1 AMPs modulate neuronal ion channels at sub-lytic concentrations

Melittin potentiates AMPA receptors at 50 ng/mL (sub-lytic, receptor-subtype selective — ruling out crude disruption; Bhatt 1992). NMDA receptors are mechanosensitive to bilayer tension (Bhatt PNAS 2007). LL-37 modulates P2X7 at sub-lytic doses (Elssner JBC 2004). Cationic arginine-rich peptides (CARPs) block NMDA electrostatically — **charge-dependent, sequence-independent** (Front Neurol 2020). **[E]**

### 6.2 LL-37: clean precedent for dose-dependent signal inversion

LL-37 (human cathelicidin): antimicrobial at 0.02–16 µM; neuroprotective at physiological levels; causes Alzheimer-like phenotypes in mouse and monkey when accumulated (elevated Aβ, neurofibrillary tangles, neuronal death — Mol Psychiatry 2022). Binds Aβ. Cationic, amphipathic, membrane-associated. The structural parallel with PrP N-terminal is exact. **[E]**

### 6.3 LLPS → amyloid is a general mechanism

Validated for α-synuclein, tau, Aβ, FUS, TDP-43. Same sequence regions drive both LLPS and fibrilization. Two distinct thresholds: Csat for LLPS onset, a higher one for liquid → solid. For PrP: N-terminal polybasic motifs drive LLPS (Yi Liang / Cong Liu groups). Dog PrP residues Arg177/Asp159 slow LLPS and inhibit human PrP amyloid formation — a structural explanation for canine resistance independent of the N-terminal tail. **[E]**

### 6.4 NINJ1 and the fuzzy coat principle

NINJ1 (Nature 2024): auto-inhibited homodimer → re-oligomerizes into filaments exposing the membrane-rupturing face → "cookie-cutter" membrane rupture. Release mechanism = change of oligomeric state (not phosphorylation, not cleavage). Directly cytotoxic. Closest structural analogy to PrP conversion. **[E]**

The amyloid fuzzy coat (JBC 2023 review): fibrils universally consist of an ordered core flanked by disordered regions that mediate membrane binding/disruption. The structural logic of v5 is a general principle of amyloid biology, not a PrP-specific claim. **[E]**

### 6.5 Charge gradient correlates with disease severity

| Protein | Flanking charge | Membrane mechanism | AMP activity | Disease timescale |
|---------|-----------------|-------------------|-------------|-------------------|
| **PrP** | **+7 to +9** | Carpet (no helix) | **Yes** (Mangé 2009; prionins 2026) | Months |
| **α-synuclein** | +4 | Amphipathic α-helix | Yes (Park 2016) | Years |
| **Tau** | Polyampholyte | Electrostatic binding | Not demonstrated | Decades |
| **Aβ42** | −3 | Oligomeric pore (different mechanism) | No (anionic) | Decades |

PrP has the highest net positive charge of any amyloid flanking region. The severity gradient (months → years → decades) is consistent with charge, though confounded by replication speed and cellular localization. **[H]**

### 6.6 Heparan sulfate binds KKRPKP — reinterpretation of pentosan polysulfate

Three HS binding sites on PrP: residues 23–52, 53–93, 110–128 (JBC 2002). Polyanions neutralize AMPs by electrostatic binding (endogenous mechanism — 2-O-sulfated HS neutralizes LL-37 at nanomolar affinity). Pentosan polysulfate (PPS) inhibited toxicity in cells but stimulated conversion cell-free — a paradox resolved under v5: PPS neutralizes the freed AMP (toxicity), not the conversion. PPS failed clinically due to delivery (intraventricular infusion), not mechanism. **[H]**

---

## 7. Resolved and Open Gaps

| # | Gap | Status | v5 answer / candidate |
|---|-----|--------|----------------------|
| L1 | PrP^L identity in synaptotoxicity assay | **Candidate** | LLPS intermediates on the liquid → solid transition pathway; β-oligomers; sarkosyl-sensitive assemblies. PrP^L has never been positively identified — 5 candidates compete, none confirmed. |
| L2 | Dose free N-terminal fragment during infection | **Open** | Testable, now motivated by prionins. No one has measured this. |
| L3 | Infectivity titers in N-terminal-truncated mice | **Explained** | Shorter tail = less potent AMP (lower charge). Confound persists: deletions also impair conversion. |
| L4 | Why α-cleavage protects if the fragment is toxic | **Resolved** | Substrate destruction: α-cleavage at ~111 cuts within the fibril core region; neither C1 nor N1 can serve as conversion substrate. N1 transgenic does NOT protect in vivo (Bhullar 2020 — N1 fails ER translocation, never reaches extracellular space). Signal inversion is unnecessary. |
| L5 | What the tail contacts in the bilayer | **Resolved** | AMP carpet mechanism (charge-dependent, not receptor-mediated). NMDA modulated by bilayer mechanosensitivity + PLA2 → arachidonic acid. |
| L6 | Is 23–33 both binding site and effector? | **Resolved** | Yes, same site. Charge-dependent (+4 in 6 residues), not receptor-specific. CARPs demonstrate any sufficiently cationic peptide does the same. |
| L7 | Strain → clinical topography | **Strong candidate** | Tripartite interaction: fibril conformation × regional lipid raft composition (sphingomyelin, cholesterol, sialylation) × cell-type-specific endocytic routes. Glycosylation is readout, not driver. |
| L8 | First misfolded molecule in sporadic disease | **Resolved** | Convergence of three age-related vulnerabilities: (1) endosomal transit through MVB (pH 4.5, geometric concentration, GPI brake released on ILVs); (2) Cu²⁺ + H₂O₂ triggers liquid → solid transition (Evans Sci Adv 2023); (3) lipid raft degradation with aging weakens the GPI brake. ~1/million/year incidence reflects probability all three align. Deep sequencing of 205 sCJD cases found no somatic PRNP mutations (2024). |
| L9 | Physiological function of PrP^C | **Strong candidate** | Ancestral AMP (confirmed by prionins, Nature Microbiology 2026) + copper buffering at synapses via liquid LLPS condensates (Evans 2023 — copper keeps condensates liquid; proposed function = synaptic copper homeostasis). Prion disease = friendly fire — the endogenous AMP turns on its own membranes. Parallels Moir-Tanzi model for Aβ. |

---

## 8. Confidence Levels

| Hypothesis | Confidence | Basis |
|-----------|-----------|-------|
| N-terminal tail is the toxic effector | 90% | Direct demonstration in vivo (Sci Adv 2026) + convergence of 4 independent systems |
| The effector mechanism is AMP-like membrane perturbation | 75% | Prionins (Nat Microbiol 2026); LL-37 parallel; CARP charge-dependence. Not directly tested for PrP in membrane context |
| Dose threshold exists | 80% | Sci Adv 2026 + sub-lytic AMP literature |
| Conversion is the unlocking event | 70% | Cryo-EM structure + G127V protector; not directly tested |
| GPI is a dual-function regulator (brake + amplifier) | 65% | Each role demonstrated independently (Tatzelt 2025; Chesebro 2005); the combined framing is ours |
| NMDA modulation is via mechanosensitivity, not receptor binding | 60% | Established for melittin/CARPs; never tested for PrP N-terminal specifically |
| Sporadic disease initiates by endosomal LLPS + Cu/ROS | 55% | Each component has evidence; the convergence model is ours. Never proposed in literature |
| PrP^L = LLPS intermediates | 35% | Five candidates compete; none confirmed |
| Synaptic selectivity is via membrane curvature | 50% | AMP-curvature preference established (Vanni JBC 2008); spine curvature data exists; never connected for prions |
| Charge gradient explains severity difference across amyloidoses | 40% | Correlative; confounded by replication speed and localization |

---

## 9. Model Revision History

| Version | Proposal | What killed it |
|---------|----------|----------------|
| v1 | PrP^Sc aggregates are the toxic agent | Weak correlation; disease without plaques; purified prions not toxic |
| v2 | Multivalent avidity clusters PrP^C → unlocks tail | Species-incompatible PrP^Sc binds but is not toxic (Gatdula 2026) |
| v3 | Conversion destroys the lock; effector requires GPI | Soluble N-terminal fragment is toxic on its own (Sci Adv 2026) |
| v4 | Released tail is a specific effector with signal inversion | N1 does not protect in vivo (Bhullar 2020); toxicity is charge-dependent, not sequence-specific (CARPs); LLPS as upstream pathway not considered |
| **v5** | **Ancestral AMP released by general amyloid conversion; specificity in autocatalytic replication; tropism by lipid composition; sporadic initiation by endosomal convergence** | *Open* |

Ideas that survived all revisions: the N-terminal tail is the effector; a dose threshold exists; conversion is the unlocking event. What changed from v4 to v5: specificity migrated from the *toxicity mechanism* to the *replication mechanism*.

---

## 10. Honest Limitations

**On method.** This investigation's confidence grew monotonically — each search round found support. This is structurally suspect. Searches were directed by hypotheses the model already favored, creating confirmation bias. The model was corrected four times substantively (v1→v5) — a good sign — but was never abandoned.

**On search bias.** The cross-disciplinary findings (AMPs, LLPS, NINJ1, LL-37) came from fields that do not cite the prion field. This virtually guarantees that other reframings remain untried. Searches were in English and Chinese (translated); other languages may contain relevant findings.

**On negative results.** "Nobody did it" and "somebody did it and got nothing" are indistinguishable by search.

**On elegance.** The model grew more elegant with each iteration. All gaps found candidates. The therapeutic combination has a coherent narrative. Biological systems are rarely this clean. The ease with which pieces fit may mean selection of pieces that fit.

**On untested predictions.** The central predictions of v5 (memantine protects in vivo; HS neutralizes toxicity; AMPs rise during neuroinflammation; charge determines flank potency) are all testable but none has been tested. A model with elegant predictions and zero direct experimental validation is a well-dressed hypothesis, not a fact.

**On "friendly fire."** The evolutionary narrative (PrP as ancestral AMP → disease = friendly fire) is seductive and may have captured the investigation. The same narrative was proposed for Aβ (Moir-Tanzi) and is contested. Narrative is not evidence.
