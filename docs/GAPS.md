# Gaps — Resolved, Candidate, and Open

Status of all identified knowledge gaps, from the original dossier (L1–L9) through those discovered during the v5 investigation (L10–L15).

---

## Original gaps (L1–L9)

### L1. PrPL in synaptotoxicity assay — `CANDIDATE`

**Original question.** The protease-sensitive species that accumulates in Sandberg's phase 2 has never been tested in a synaptotoxicity assay.

**v5 status.** PrPL has never been positively identified. It is operationally defined by subtraction (total PrP − PrPC − PK-resistant PrPSc). Five candidates compete:

| Candidate | Source | Mechanism |
|-----------|--------|-----------|
| Released N-terminal tails | v5 hypothesis | AMP freed by fibril formation |
| β-rich oligomers (PrP*) | Bhérer/Bhullar groups | Soluble low-MW oligomers, ~20 monomers |
| Sarkosyl-sensitive assemblies | Benilova/Collinge 2020 | Sarkosyl destroys toxicity but not infectivity |
| Syntaxin-6-stabilized intermediates | Sangar et al., eLife 2024 | On-pathway transient intermediates |
| PK-resistant multimers | Bhérer 2018 (contradicts Sandberg) | SEC fractions 5–10 |

**Key constraint.** Benilova et al. (PNAS 2020): the toxic species is sarkosyl-sensitive; infectious prions are sarkosyl-resistant. Whatever PrPL is, it dissolves in sarkosyl. No mass spectrometry has been applied.

**What would resolve it.** Experiment E2 (quantify free N-terminal fragment throughout infection).

---

### L2. Dose free N-terminal fragment throughout infection — `OPEN`

**Original question.** Nobody has measured free N-terminal fragment (23–90) levels during prion infection.

**v5 status.** Still open. Now motivated by the prionins discovery (Nature Microbiology 2026) confirming that PrP fragments are functional AMPs. Experiment E2 directly addresses this.

---

### L3. Infectivity titers in truncated mice at terminal stage — `EXPLAINED (with caveat)`

**Original question.** PrPΔ23-88 mice still die of prion disease (161 days vs 50 for WT overexpressors). If the tail is the effector, why don't they survive longer?

**v5 explanation.** Shorter tail = less potent AMP (lower charge). But the confound persists: deletions also impair conversion efficiency. The two effects (less toxin and less conversion) cannot be cleanly separated in this model. Minikel (CureFFI 2015) notes that Δ23-88 mice still die, arguing downstream interventions are marginal vs anti-propagation.

**Caveat.** The deletion ladder is not a clean test of v5 because conversion and toxicity are coupled through the same sequence.

---

### L4. Why α-cleavage protects if the fragment is toxic — `RESOLVED`

**Original question.** N1 (23–111) was claimed to be neuroprotective (sequestering oligomers), but the same fragment appears toxic above a threshold.

**v5 resolution.** The paradox dissolves. α-cleavage protects by **substrate destruction**, not by releasing a protective fragment:

1. α-cleavage cuts at ~residue 111, within the fibril core (94–225).
2. Neither N1 (23–111) nor C1 (112–231) retains the full core sequence needed for conversion.
3. C1 on the membrane cannot be converted to PrPSc and may inhibit conversion (Bhérer et al., Prion 2012).
4. N1 transgenic mice show NO survival benefit (Bhullar et al., 2020, PMC7253391) — N1 fails ER translocation, retains its signal peptide, and accumulates in the cytosol.
5. The "dose-dependent signal inversion" hypothesis from v4 is unnecessary.

**Supporting evidence.** ADAM10 overexpression correlates with longer survival (Altmeppen et al., eLife 2015) — consistent with more substrate destruction, not more N1 release.

---

### L5. What the tail touches in the bilayer — `RESOLVED`

**Original question.** Does the freed N-terminal form a pore, modulate a channel, or act through a specific receptor?

**v5 resolution.** None of the above. It acts as a **carpet-model AMP**:

1. Polycationic peptide (charge +4 from KKRPKP alone, +7 to +9 total) inserts into the outer leaflet without forming a transmembrane pore.
2. Asymmetric insertion causes local membrane thinning and tension.
3. NMDA receptors are mechanosensitive to bilayer tension (Bhatt et al., PNAS 2007) — Mg²⁺ block relieved → Ca²⁺ influx.
4. PLA2 activation by the inserted peptide → arachidonic acid → additional NMDA potentiation.
5. No specific receptor is required. CARPs demonstrate NMDA modulation is charge-dependent, sequence-independent (Front Neurol 2020).
6. Melittin at sub-lytic concentrations (50 ng/mL) selectively potentiates AMPA-mediated Ca²⁺ influx with receptor-subtype selectivity (Bhatt et al., 1992).

---

### L6. Is 23–33 both the binding site and the effector? — `RESOLVED`

**Original question.** If the same residues (KKRPKP, 23–28) bind the membrane and cause toxicity, there is no "safe" separated target.

**v5 resolution.** Yes, the same site. The binding and effector functions are both charge-dependent, not receptor-mediated. KKRPKP contributes +4 in 6 residues — the highest charge density segment. This is the same charge density as melittin's C-terminal basic cluster (KRKRQQ). There is no separate receptor to target; the interaction is electrostatic with the anionic membrane surface.

**Therapeutic implication.** The target is charge neutralization (polyanions binding KKRPKP), not receptor blockade.

---

### L7. How strain conformation becomes clinical topography — `CANDIDATE (strong)`

**Original question.** No model explains why different prion strains target different brain regions in the same host.

**v5 candidate.** A three-way interaction:

1. **Fibril conformation (strain)** selects among regionally distributed cofactors (Bhérer et al., PLoS Pathog 2020 — cofactor and glycosylation preferences are strain-determined).
2. **Regional lipid composition** — sphingomyelin depletion → 4× more PrPSc; cholesterol depletion → less PrPSc (PMC9871914). Brain lipidome is region- and cell-type-specific (Nat Commun 2024, 419 species profiled).
3. **Regional sialylation** — thalamus PrPSc is less sialylated, replicates faster (Makarava et al., IJMS 2020). Lower sialylation = faster replication = greater vulnerability.

Glycosylation ratio is a readout, not a driver (Tuzi et al., J Virol 2008 — deglycosylated PrPSc maintains strain-specific tropism). PMCA with region-specific homogenates reproduces in vivo tropism (Barria et al., JBC 2018).

---

### L8. Origin of the first misfolded molecule in sporadic disease — `RESOLVED (hypothesis)`

**Original question.** 85% of prion cases are sporadic. What triggers the first conversion?

**v5 resolution.** Three age-related vulnerabilities converge:

1. **Endosomal transit** — the multivesicular body (MVB/late endosome) is the major intracellular site of prion conversion (J Cell Sci 2015, PMC4379730). pH 4.5–5.5 destabilizes PrP. On intraluminal vesicles, the GPI anchor points inward — PrP's ectodomain is exposed to the lumen, free from 2D membrane constraint. The topological brake is released.
2. **Copper + oxidative stress** — Cu²⁺ drives PrP phase separation at physiological ratios; copper alone keeps condensates liquid (Evans et al., Sci Adv 2023). Cu + H₂O₂ triggers the liquid-to-solid transition → amyloid with dityrosine cross-links. Age-related oxidative stress pushes the condensate past the transition.
3. **Lipid raft degradation with aging** — progressive loss of sphingolipids and cholesterol weakens the GPI brake (PMC7193971).

The ~1/million/year incidence reflects the probability that all three factors align in the same cell, at the same moment, with sufficient local concentration to nucleate a self-propagating seed.

**Supporting evidence.** Deep sequencing of PRNP in 205 sCJD cases found no enrichment of somatic mutations vs controls (PMC11328154, 2024). The Y145Stop variant attains self-seeding amyloid state via LLPS without the classic misfolding step (PNAS 2021). Bank vole PrP spontaneously generates transmissible prions (PNAS 2012).

---

### L9. Normal physiological function of PrPC — `CANDIDATE (strong)`

**Original question.** Without knowing PrPC's function, it is hard to know what breaking it means.

**v5 candidate.** Dual function:

1. **Ancestral antimicrobial peptide.** The N-terminal has direct antimicrobial activity against Gram+, Gram−, and fungi (Mangé/Pasupuleti et al., PLoS One 2009). 1,179 AMP fragments from prion-like proteins — "prionins" — confirmed by deep-learning screen and synthesis (Nature Microbiology 2026). PrP expression increases during wounding. The Moir–Tanzi model for Aβ (antimicrobial peptide → pathological when chronically activated) applies directly.

2. **Copper buffering at synapses via LLPS.** Cu²⁺ at physiological ratios drives PrP into liquid condensates that sequester free copper (Evans et al., Sci Adv 2023). This prevents Fenton-reaction toxicity from free Cu²⁺. The condensate is maintained in liquid state by the GPI-anchor topological brake. Prion disease disrupts both functions: the AMP is released against the host's own membranes (friendly fire), and copper buffering fails as condensates solidify.

---

## New gaps (L10–L15)

### L10. PrPL identity — `OPEN`

Five competing candidates (see L1), zero direct evidence for any. No mass spectrometry, no specific antibody, no isolation. Operationally defined by subtraction. The most basic question in the temporal model is unanswered.

**What would resolve it.** SEC fractionation of phase-2 brain homogenate → mass spectrometry on the protease-sensitive fraction → identify all PrP-derived species present.

---

### L11. Does LLPS initiate sporadic disease in vivo? — `OPEN`

The LLPS → sporadic disease hypothesis (L8) is mechanistically supported but has **never been explicitly proposed or tested**. The pieces are published under separate labels (membrane biophysics, copper biology, aging proteostasis) but no one has connected them.

**What would resolve it.** Experiment E8 (measure PrP Csat as a function of age/PTMs). Also: detect LLPS condensates in aging brain tissue (would require live imaging or cryo-ET in aged animals).

---

### L12. What are the endosomal cofactors? — `OPEN`

The endosome/MVB is the conversion site (J Cell Sci 2015), and region-specific PMCA shows non-PrPC cofactors contribute to tropism (Barria et al., JBC 2018 — correcting for PrPC levels reduces but does not eliminate the region-specific pattern, r = 0.65 → 0.40, still p = 0.01). But the cofactors themselves — lipids? RNA? heparan sulfate proteoglycans? — are unidentified.

**What would resolve it.** Proteomics/lipidomics on endosomal fractions from brain regions with different prion susceptibility, compared between strains.

---

### L13. Is the neuroinflammation → AMP feedback loop real? — `OPEN`

v5 predicts a positive feedback: freed N-terminal tail → membrane perturbation → neuroinflammation → microglial activation → upregulation of endogenous AMPs (LL-37/CRAMP) → additional membrane perturbation → more inflammation. If this loop is real, anti-inflammatory interventions that suppress AMP production would have a secondary neuroprotective benefit.

**What would resolve it.** Experiment E6 (AMP profiling throughout infection). If CRAMP levels rise with microglial activation during prion disease — which has never been measured — the loop is plausible.

---

### L14. Why do astrocytes propagate prions without dying? — `OPEN`

Astrocytes propagate prions (Bhérer et al., Biology 2024) but do not degenerate. Under v5, possible explanations:

1. Astrocyte membranes have different lipid composition (lower sphingolipid content, simpler glycolipids) — possibly less susceptible to AMP insertion at raft boundaries.
2. Astrocytes have higher antioxidant capacity — less ROS → condensates stay liquid → less solidification → less tail release.
3. Astrocytes express different NMDA receptor subunit compositions — possibly less mechanosensitive.
4. Astrocytes are the brain's primary cholesterol producers — their raft composition may resist AMP perturbation differently.

**What would resolve it.** Compare AMP-peptide toxicity assay (E7) in primary neurons vs primary astrocytes. If astrocytes are resistant to the same charge-gradient peptides that kill neurons, the difference is in the target cell, not the replication machinery.

---

### L15. Can CWD cross the species barrier to humans? — `OPEN`

Chronic wasting disease is the largest ongoing prion epidemic (cervids across North America). Under v5, the species barrier is about conversion efficiency (globular domain compatibility + cofactor availability), not toxicity (the N-terminal charge is identical across mammals at +6.5). The question is whether cervid PrPSc can template human PrPC conversion.

**What would resolve it.** In vitro: human PrPC in PMCA seeded with CWD PrPSc (has been done — very inefficient but not zero). In vivo: humanized transgenic mice inoculated with CWD (has been done — long incubation, subclinical infection in some models). Epidemiology: monitor CJD incidence in high-exposure populations (hunters, taxidermists) over decades.

**v5 note.** The identical N-terminal charge across species means that IF conversion occurs, the toxicity mechanism would be the same. The barrier is entirely at the conversion step.
