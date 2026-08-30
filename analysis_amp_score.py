"""
Formal AMP (antimicrobial peptide) prediction scores for PrP fragments.
Compares with known AMPs and amyloid flanking regions.
Uses published sequence-based classification features.
"""
import numpy as np
import json

# ============================================================
# SCALES AND CONSTANTS
# ============================================================

HYDROPHOBIC_AA = set('AVILMFWP')

KD_SCALE = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'E': -3.5, 'Q': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2,
}

BOMAN_SCALE = {
    'A': -0.01, 'R': 0.04, 'N': 0.06, 'D': 0.15, 'C': -0.24,
    'E': 0.33, 'Q': -0.30, 'G': -0.01, 'H': 0.13, 'I': -0.80,
    'L': -0.56, 'K': -0.23, 'M': -0.28, 'F': -0.52, 'P': -0.09,
    'S': 0.11, 'T': 0.01, 'W': 0.07, 'Y': -0.17, 'V': -0.47,
}

# Instability index dipeptide weights (Guruprasad et al. 1990)
# Simplified: use average instability contribution per amino acid
INSTABILITY_WEIGHTS = {
    'A': 0.06, 'R': 0.27, 'N': 0.20, 'D': 0.44, 'C': 0.11,
    'E': 0.40, 'Q': 0.10, 'G': 0.16, 'H': 0.22, 'I': 0.04,
    'L': 0.05, 'K': 0.26, 'M': 0.08, 'F': 0.03, 'P': 0.32,
    'S': 0.16, 'T': 0.15, 'W': 0.05, 'Y': 0.09, 'V': 0.04,
}

PK_VALUES = {
    'K': 10.5, 'R': 12.4, 'H': 6.0,
    'D': 3.65, 'E': 4.25,
    'N_term': 9.69, 'C_term': 2.34,
}

# ============================================================
# FEATURE CALCULATORS
# ============================================================

def net_charge(seq, ph=7.4):
    charge = 0.0
    for aa in seq.upper():
        if aa in ('K', 'R'):
            pk = PK_VALUES[aa]
            charge += 1.0 / (1.0 + 10**(ph - pk))
        elif aa == 'H':
            pk = PK_VALUES['H']
            charge += 1.0 / (1.0 + 10**(ph - pk))
        elif aa in ('D', 'E'):
            pk = PK_VALUES[aa]
            charge -= 1.0 / (1.0 + 10**(pk - ph))
    # Terminal charges
    charge += 1.0 / (1.0 + 10**(ph - PK_VALUES['N_term']))
    charge -= 1.0 / (1.0 + 10**(PK_VALUES['C_term'] - ph))
    return round(charge, 2)

def hydrophobic_ratio(seq):
    count = sum(1 for aa in seq.upper() if aa in HYDROPHOBIC_AA)
    return round(count / len(seq), 3)

def gravy(seq):
    vals = [KD_SCALE.get(aa, 0) for aa in seq.upper()]
    return round(np.mean(vals), 3) if vals else 0

def boman_index(seq):
    vals = [BOMAN_SCALE.get(aa, 0) for aa in seq.upper()]
    return round(np.mean(vals), 3) if vals else 0

def hydrophobic_moment(seq, angle=100, window=11):
    if len(seq) < window:
        window = len(seq)
    moments = []
    for i in range(len(seq) - window + 1):
        win = seq[i:i+window]
        hx = sum(KD_SCALE.get(aa, 0) * np.cos(np.radians(angle * j)) for j, aa in enumerate(win))
        hy = sum(KD_SCALE.get(aa, 0) * np.sin(np.radians(angle * j)) for j, aa in enumerate(win))
        moments.append(np.sqrt(hx**2 + hy**2) / window)
    return round(max(moments), 3) if moments else 0

def charge_density(seq):
    return round(net_charge(seq) / len(seq), 4)

def instability_index(seq):
    vals = [INSTABILITY_WEIGHTS.get(aa, 0.1) for aa in seq.upper()]
    ii = (10.0 / len(seq)) * sum(vals) * len(seq)
    return round(ii, 1)

def isoelectric_point(seq):
    for ph_test in np.arange(0, 14, 0.01):
        c = net_charge(seq, ph=ph_test)
        if c <= 0:
            return round(ph_test, 2)
    return 14.0

def amp_likelihood_score(charge, hyd_ratio, amphipathicity, chg_density, gravy_val):
    score = 0
    if charge >= 2:
        score += 20
    if 0.30 <= hyd_ratio <= 0.60:
        score += 20
    if amphipathicity > 0.3:
        score += 20
    if chg_density > 0.05:
        score += 20
    if gravy_val < 0:
        score += 20
    return score

# ============================================================
# PEPTIDES
# ============================================================

PEPTIDES = [
    ("KKRPKP (PrP 23-28)", "KKRPKP", "PrP"),
    ("PrP 23-35", "KKRPKPGGWNTGG", "PrP"),
    ("PrP 23-93 (full N-term)", "KKRPKPGGWNTGGSRYPGQGSPGGNRYPPQGGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQGGGTHSQW", "PrP"),
    ("NNRPNP (neutral ctrl)", "NNRPNP", "Control"),
    ("LL-37 (cathelicidin)", "LLGDFFRKSKEKIGKEFKRIVQRIKDFLRNLVPRTES", "Known AMP"),
    ("Melittin (bee venom)", "GIGAVLKVLTTGLPALISWIKRKRQQ", "Known AMP"),
    ("Magainin 2 (frog)", "GIGKFLHSAKKFGKAFVGEIMNS", "Known AMP"),
    ("a-Syn 1-30", "MDVFMKGLSKAKEGVVAAAEKTKQGVAEAA", "Amyloid flank"),
    ("Tau 1-30", "MAEPRQEFEVMEDHAGTYGLGDRKDQGGYT", "Amyloid flank"),
]

# ============================================================
# ANALYSIS
# ============================================================

print("=" * 100)
print("FORMAL AMP PREDICTION SCORES")
print("=" * 100)

results = []

header = f"{'Peptide':<28} {'Len':>4} {'Chg':>6} {'H%':>6} {'GRAVY':>7} {'uH':>6} {'Boman':>7} {'ChgD':>7} {'pI':>6} {'AMP%':>5} {'Type':<12}"
print(f"\n{header}")
print("-" * 100)

for name, seq, ptype in PEPTIDES:
    chg = net_charge(seq)
    hr = hydrophobic_ratio(seq)
    g = gravy(seq)
    mu = hydrophobic_moment(seq)
    bi = boman_index(seq)
    cd = charge_density(seq)
    pi = isoelectric_point(seq)
    amp = amp_likelihood_score(chg, hr, mu, cd, g)

    results.append({
        "name": name, "sequence": seq, "type": ptype, "length": len(seq),
        "charge": chg, "hydrophobic_ratio": hr, "gravy": g,
        "amphipathicity": mu, "boman": bi, "charge_density": cd,
        "pI": pi, "amp_score": amp,
    })

    print(f"{name:<28} {len(seq):>4} {chg:>+6.1f} {hr:>6.1%} {g:>7.2f} {mu:>6.3f} {bi:>7.3f} {cd:>7.4f} {pi:>6.2f} {amp:>5} {ptype:<12}")

# Sort by AMP score
print(f"\n{'=' * 60}")
print("RANKED BY AMP LIKELIHOOD SCORE")
print(f"{'=' * 60}")
sorted_results = sorted(results, key=lambda x: x['amp_score'], reverse=True)
for i, r in enumerate(sorted_results):
    bar = "#" * (r['amp_score'] // 5)
    print(f"  {i+1}. {r['name']:<28} AMP={r['amp_score']:>3}% {bar}")

# Interpretation
print(f"\n{'=' * 60}")
print("INTERPRETATION")
print(f"{'=' * 60}")

prp_scores = [r for r in results if r['type'] == 'PrP']
amp_scores = [r for r in results if r['type'] == 'Known AMP']
ctrl_scores = [r for r in results if r['type'] == 'Control']
flank_scores = [r for r in results if r['type'] == 'Amyloid flank']

prp_mean = np.mean([r['amp_score'] for r in prp_scores])
amp_mean = np.mean([r['amp_score'] for r in amp_scores])
ctrl_mean = np.mean([r['amp_score'] for r in ctrl_scores]) if ctrl_scores else 0
flank_mean = np.mean([r['amp_score'] for r in flank_scores])

print(f"\n  Mean AMP scores by category:")
print(f"    Known AMPs:       {amp_mean:.0f}%")
print(f"    PrP fragments:    {prp_mean:.0f}%")
print(f"    Amyloid flanks:   {flank_mean:.0f}%")
print(f"    Neutral control:  {ctrl_mean:.0f}%")

print(f"\n  v5 model predictions:")
prp93 = next(r for r in results if '23-93' in r['name'])
nn = next(r for r in results if 'NNRPNP' in r['name'])
print(f"    PrP 23-93 AMP score ({prp93['amp_score']}%) vs NNRPNP ({nn['amp_score']}%): ", end="")
print("CONFIRMED" if prp93['amp_score'] > nn['amp_score'] else "REFUTED")

print(f"    PrP fragments in AMP range? ", end="")
print("YES" if prp_mean >= amp_mean * 0.6 else "NO")

print(f"    PrP > other amyloid flanks? ", end="")
print("CONFIRMED" if prp_mean > flank_mean else "REFUTED")

with open("/Users/allan/Projects/cjd/amp_scores.json", "w") as f:
    json.dump(results, f, indent=2)
print(f"\nSaved to amp_scores.json")
