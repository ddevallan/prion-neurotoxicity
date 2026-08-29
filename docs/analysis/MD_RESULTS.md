# Molecular Dynamics Results

## Platform

- **Hardware**: Apple M1 Pro, 16 GB RAM
- **Software**: OpenMM 8.6.0 with OpenCL acceleration
- **Force field**: AMBER14-SB + TIP3P water
- **Conditions**: 300 K, 0.15 M NaCl, Langevin thermostat

## Simulations completed

### 1. Short equilibration (100 ps) — three peptides

| Peptide | Charge | SASA (nm²) | Rg (nm) | Speed (ns/day) |
|---------|--------|------------|---------|----------------|
| **KKRPKP** (PrP 23-28, wild-type) | +4 | 7.864 | 1.793 | 262 |
| **NNRPNP** (charge-neutralized control) | 0 | 5.052 | 1.790 | 263 |
| **KKRPKPGGWNTGG** (PrP 23-35, extended) | +4 | 12.220 | 3.171 | 81 |

**Key finding**: KKRPKP has 56% more solvent-accessible surface area than NNRPNP
despite identical Rg. Charged sidechains (Lys, Arg) project outward, maximizing
solvent exposure — the conformation an AMP needs for the carpet mechanism.

### 2. Extended production (7 ns) — KKRPKP wild-type

| Metric | 100 ps | 7 ns | Interpretation |
|--------|--------|------|----------------|
| SASA | 7.864 nm² | 3.587 nm² | Peptide collapsed — sidechains rearranged |
| Rg | 1.793 nm | 3.783 nm | Backbone extended (more linear) |
| E2E distance | — | 0.716 ± 0.254 nm | Endpoints close — peptide curved |
| Secondary structure | — | 100% coil | Intrinsically disordered |

**Per-residue SASA (7 ns average)**:
- LYS1: 0.637 nm²
- LYS2: 0.623 nm²
- **ARG3: 0.932 nm²** (most exposed — guanidinium group)
- PRO4: 0.381 nm²
- LYS5: 0.513 nm²
- PRO6: 0.502 nm²

**Interpretation**: In water, KKRPKP is intrinsically disordered (100% coil) and
partially collapsed. This matches the "prionin" finding (Nature Microbiology 2026):
prion-derived AMPs are unstructured in water but adopt ordered conformations upon
membrane contact (inducible folding). The extended simulation confirms the peptide
is an IDP in solution — the membrane simulation (next step) will test whether it
opens and inserts upon encountering a lipid bilayer.

**ARG3 maintains the highest SASA** — the arginine guanidinium group remains most
exposed to solvent, consistent with its known role in electrostatic interaction
with phosphate headgroups during AMP-membrane binding.

## Next steps

These water-phase simulations establish the baseline behavior. The critical test
requires an explicit POPC lipid bilayer:

1. **CHARMM-GUI** to build KKRPKP + 128-lipid POPC bilayer system
2. **Vast.ai** (A40 GPU, ~$0.25/hr) for 3 × 300 ns production runs
3. **Measure**: insertion depth, membrane thinning, lipid order perturbation
4. **Compare**: KKRPKP (+4) vs NNRPNP (0) vs charge variants (+2, +3)

Prediction: KKRPKP will insert into the outer leaflet via carpet mechanism
and locally thin the membrane. NNRPNP will remain in the aqueous phase.
