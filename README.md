# Prion Disease Neurotoxicity — Exploratory Investigation

An exploratory investigation into the mechanism of neurotoxicity in prion disease,
conducted by reverse engineering: starting from observed damage and reconstructing
the causal chain to the molecular property that initiates it.

## The v5 Model

Prion disease is a specific instance of a general amyloid toxicity mechanism — the
disordered N-terminal "fuzzy coat" freed by fibril formation perturbs neuronal membranes
like an antimicrobial peptide (AMP). Prion-specific features are not in the toxicity
(which is generic cationic-peptide biophysics) but in three properties: truly
autocatalytic self-templating replication (~1000× faster than tau), the GPI anchor as a
dual brake/amplifier, and strain conformation selecting among regionally distributed
cofactors to produce tropism.

## Status

**Checkpoint: August 28, 2026.** Model version 5. Not peer-reviewed. Not clinical guidance.

Investigation conducted via 24 cross-disciplinary literature searches spanning prion
biology, antimicrobial peptides, LLPS/phase separation, membrane biophysics, cell death
(MLKL, gasdermin, NINJ1), receptor signaling (Eph/ephrin), copper homeostasis, and the
full therapeutic pipeline. Three computational analyses completed; MD simulations in
progress.

## Repository Structure

```
cjd/
├── README.md
├── dossie-v5.html              # Full dossier (publishable artifact)
│
├── docs/
│   ├── MODEL.md                # v5 model description
│   ├── CROSS_DISCIPLINARY.md   # Findings from adjacent fields
│   ├── THERAPEUTICS.md         # Treatment landscape + v5 predictions
│   ├── EXPERIMENTS.md          # Proposed experiments
│   ├── LIMITATIONS.md          # Honest limitations and caveats
│   ├── analysis/
│   │   └── COMPUTATIONAL.md    # Computational analysis results
│   └── references/
│       └── BIBLIOGRAPHY.md     # Complete bibliography (~150 entries)
│
├── analysis_charge.py          # Cross-species PrP N-terminal charge analysis
├── analysis_kinetic.py         # Kinetic model of therapeutic window (ODE)
├── analysis_mutations.py       # PRNP mutations mapped to v5 domains
├── charge_analysis.json        # Output: charge data
├── kinetic_model.json          # Output: kinetic model curves
├── mutations_analysis.json     # Output: mutation mapping
│
├── md_setup.py                 # MD: peptide-in-water equilibration (OpenMM)
├── md_membrane.py              # MD: extended 10 ns simulation
├── md_output/                  # MD trajectories and structures
├── md_membrane_output/         # Extended MD output
│
├── runpod/                     # Cloud GPU simulation setup
│   ├── CHARMM_GUI_WORKFLOW.md  # Step-by-step membrane simulation guide
│   ├── run_production.py       # Production MD for CHARMM-GUI systems
│   ├── analyze_membrane.py     # Post-simulation analysis
│   ├── Dockerfile
│   └── requirements.txt
│
├── pyproject.toml
└── uv.lock
```

## Running the Analyses

Requires Python 3.14+ and [uv](https://docs.astral.sh/uv/).

```bash
# Cross-species charge analysis
uv run python analysis_charge.py

# Kinetic model (therapeutic window)
uv run python analysis_kinetic.py

# Mutation mapping to v5 domains
uv run python analysis_mutations.py

# MD simulation — peptide in water (requires OpenMM)
uv run python md_setup.py

# Extended MD — 10 ns (runs ~40 min on M1 Pro)
uv run python md_membrane.py --steps 5000000 --peptide wt
```

## Key Findings

- PrP N-terminal fragments are literal AMPs (prionins, Nature Microbiology 2026)
- The GPI anchor suppresses LLPS (Tatzelt, PNAS 2025) — dual brake/amplifier
- G127V protects by keeping LLPS condensates liquid, not by altering toxicity
- NMDA receptors are mechanosensitive to bilayer tension (no specific receptor needed)
- Memantine has never been tested in prion-infected mice (33-year gap since in vitro data)
- Octarepeat insertions correlate r = −0.989 with age of onset (charge dose-response)
- Combination therapy (ASO + NMDA antagonist + AMP neutralizer) is predicted but untested

## Disclaimer

This repository documents an exploratory investigation, not a systematic review.
It mixes established results, contested findings, and untested hypotheses — each
labeled accordingly in the documentation. References were checked against abstracts
and text excerpts; volume/page details should be verified before formal use. No
conclusion here substitutes primary literature. Not peer-reviewed. Not clinical guidance.
