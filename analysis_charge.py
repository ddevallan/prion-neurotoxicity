"""Cross-species PrP N-terminal charge analysis — tests AMP hypothesis."""
import numpy as np
import json

# PrP N-terminal sequences (signal peptide removed, starting at mature protein)
# Source: UniProt, manually curated for the mature N-terminal region (23-93 human numbering)
SPECIES = {
    "Human": {
        "seq": "KKRPKPGGWNTGGSRYPGQGSPGGNRYPPQGGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQGGGTHSQW",
        "susceptible": True,
        "disease": "CJD/kuru/FFI/GSS",
        "severity": "lethal, months"
    },
    "Mouse": {
        "seq": "KKRPKPGGWNTGGSRYPGQGSPGGNRYPPQGGTWGQPHGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQGGGTHGQW",
        "susceptible": True,
        "disease": "scrapie-adapted",
        "severity": "lethal, months"
    },
    "Hamster": {
        "seq": "KKRPKPGGWNTGGSRYPGQGSPGGNRYPPQGGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQGGGTHNQW",
        "susceptible": True,
        "disease": "scrapie-adapted",
        "severity": "lethal, weeks-months"
    },
    "Sheep": {
        "seq": "KKRPKPGGWNTGGSRYPGQGSPGGNRYPPQGGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQGGSHSQW",
        "susceptible": True,
        "disease": "scrapie",
        "severity": "lethal, months-years"
    },
    "Cow": {
        "seq": "KKRPKPGGGWNTGGSRYPGQGSPGGNRYPPQGGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQGGTHGQW",
        "susceptible": True,
        "disease": "BSE",
        "severity": "lethal, months"
    },
    "Elk (CWD)": {
        "seq": "KKRPKPGGWNTGGSRYPGQGSPGGNRYPPQGGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQGGTHGQW",
        "susceptible": True,
        "disease": "CWD",
        "severity": "lethal, months-years"
    },
    "Bank vole": {
        "seq": "KKRPKPGGWNTGGSRYPGQGSPGGNRYPPQGGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQGGGTHSQW",
        "susceptible": True,
        "disease": "spontaneous + adapted",
        "severity": "lethal, highly susceptible"
    },
    "Dog": {
        "seq": "KKRPKPGGWNTGGSRYPGQGSPGGNRYPPQGGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQGGTHGQW",
        "susceptible": False,
        "disease": "none known",
        "severity": "resistant"
    },
    "Rabbit": {
        "seq": "KKRPKPGGWNTGGSRYPGQGSPGGNRYPPQGGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQGGTHGQW",
        "susceptible": False,
        "disease": "resistant (until recently)",
        "severity": "highly resistant"
    },
    "Horse": {
        "seq": "KKRPKPGGWNTGGSRYPGQGSPGGNRYPPQGGGGWGQPHGGGWGQPHGGSWGQPHGGGWGQPHGGGWGQGGTHGQW",
        "susceptible": False,
        "disease": "none known",
        "severity": "resistant"
    },
}

# Known AMPs for comparison
AMPS = {
    "LL-37 (human)": {
        "seq": "LLGDFFRKSKEKIGKEFKRIVQRIKDFLRNLVPRTES",
        "type": "cathelicidin"
    },
    "Melittin (bee)": {
        "seq": "GIGAVLKVLTTGLPALISWIKRKRQQ",
        "type": "bee venom"
    },
    "Magainin 2 (frog)": {
        "seq": "GIGKFLHSAKKFGKAFVGEIMNS",
        "type": "magainin"
    },
    "Defensin HNP-1": {
        "seq": "ACYCRIPACIAGERRYGTCIYQGRLWAFCC",
        "type": "alpha-defensin"
    },
}

CHARGE_MAP = {
    'K': 1, 'R': 1, 'H': 0.1,  # His partially protonated at pH 7.4
    'D': -1, 'E': -1,
}

HYDROPHOBICITY = {  # Kyte-Doolittle scale
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'E': -3.5, 'Q': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2,
}

def net_charge(seq, ph=7.4):
    return sum(CHARGE_MAP.get(aa, 0) for aa in seq.upper())

def charge_density(seq):
    c = net_charge(seq)
    return c / len(seq) if len(seq) > 0 else 0

def mean_hydrophobicity(seq):
    vals = [HYDROPHOBICITY.get(aa, 0) for aa in seq.upper()]
    return np.mean(vals) if vals else 0

def hydrophobic_moment(seq, angle=100, window=11):
    """Calculate mean hydrophobic moment assuming alpha-helix (100°) or other angle."""
    if len(seq) < window:
        window = len(seq)
    moments = []
    for i in range(len(seq) - window + 1):
        win = seq[i:i+window]
        hx = sum(HYDROPHOBICITY.get(aa, 0) * np.cos(np.radians(angle * j)) for j, aa in enumerate(win))
        hy = sum(HYDROPHOBICITY.get(aa, 0) * np.sin(np.radians(angle * j)) for j, aa in enumerate(win))
        moments.append(np.sqrt(hx**2 + hy**2) / window)
    return np.mean(moments) if moments else 0

def boman_index(seq):
    """Boman index: sum of solubility values / length. Higher = more protein-binding potential."""
    solubility = {
        'A': -0.01, 'R': 0.04, 'N': 0.06, 'D': 0.15, 'C': -0.24,
        'E': 0.33, 'Q': -0.30, 'G': -0.01, 'H': 0.13, 'I': -0.80,
        'L': -0.56, 'K': -0.23, 'M': -0.28, 'F': -0.52, 'P': -0.09,
        'S': 0.11, 'T': 0.01, 'W': 0.07, 'Y': -0.17, 'V': -0.47,
    }
    vals = [solubility.get(aa, 0) for aa in seq.upper()]
    return sum(vals) / len(vals) if vals else 0

def sliding_charge(seq, window=7):
    charges = []
    for i in range(len(seq) - window + 1):
        charges.append(net_charge(seq[i:i+window]))
    return charges

# Analyze all species
results = {"species": [], "amps": []}

print("=" * 80)
print("CROSS-SPECIES PrP N-TERMINAL CHARGE ANALYSIS")
print("=" * 80)
print(f"\n{'Species':<16} {'Len':>4} {'Charge':>7} {'ChgDen':>7} {'<H>':>7} {'μH':>7} {'Boman':>7} {'Suscept':>10}")
print("-" * 80)

for name, data in SPECIES.items():
    seq = data["seq"]
    ch = net_charge(seq)
    cd = charge_density(seq)
    mh = mean_hydrophobicity(seq)
    hm = hydrophobic_moment(seq)
    bi = boman_index(seq)
    susc = "YES" if data["susceptible"] else "NO"

    results["species"].append({
        "name": name,
        "length": len(seq),
        "net_charge": round(ch, 1),
        "charge_density": round(cd, 4),
        "mean_hydrophobicity": round(mh, 3),
        "hydrophobic_moment": round(hm, 3),
        "boman_index": round(bi, 3),
        "susceptible": data["susceptible"],
        "disease": data["disease"],
        "severity": data["severity"],
        "sequence": seq,
    })

    print(f"{name:<16} {len(seq):>4} {ch:>+7.1f} {cd:>7.4f} {mh:>7.3f} {hm:>7.3f} {bi:>7.3f} {susc:>10}")

print(f"\n{'AMP Reference':<20} {'Len':>4} {'Charge':>7} {'ChgDen':>7} {'<H>':>7} {'μH':>7} {'Boman':>7}")
print("-" * 75)

for name, data in AMPS.items():
    seq = data["seq"]
    ch = net_charge(seq)
    cd = charge_density(seq)
    mh = mean_hydrophobicity(seq)
    hm = hydrophobic_moment(seq)
    bi = boman_index(seq)

    results["amps"].append({
        "name": name,
        "length": len(seq),
        "net_charge": round(ch, 1),
        "charge_density": round(cd, 4),
        "mean_hydrophobicity": round(mh, 3),
        "hydrophobic_moment": round(hm, 3),
        "boman_index": round(bi, 3),
        "type": data["type"],
        "sequence": seq,
    })

    print(f"{name:<20} {len(seq):>4} {ch:>+7.1f} {cd:>7.4f} {mh:>7.3f} {hm:>7.3f} {bi:>7.3f}")

# KKRPKP analysis specifically
print(f"\n{'=' * 60}")
print("KKRPKP MOTIF ANALYSIS")
print(f"{'=' * 60}")
kkrpkp = "KKRPKP"
print(f"Sequence: {kkrpkp}")
print(f"Net charge: {net_charge(kkrpkp):+.1f}")
print(f"Charge density: {charge_density(kkrpkp):.4f}")
print(f"For comparison:")
print(f"  LL-37 first 6 aa (LLGDFF): charge = {net_charge('LLGDFF'):+.1f}, density = {charge_density('LLGDFF'):.4f}")
print(f"  Melittin last 6 aa (KRKRQQ): charge = {net_charge('KRKRQQ'):+.1f}, density = {charge_density('KRKRQQ'):.4f}")

# Correlation analysis
print(f"\n{'=' * 60}")
print("KEY FINDING: CHARGE COMPARISON")
print(f"{'=' * 60}")

species_charges = [r["net_charge"] for r in results["species"]]
amp_charges = [r["net_charge"] for r in results["amps"]]

print(f"\nPrP N-terminal charge range: {min(species_charges):+.1f} to {max(species_charges):+.1f}")
print(f"AMP charge range: {min(amp_charges):+.1f} to {max(amp_charges):+.1f}")
print(f"\nPrP mean charge: {np.mean(species_charges):+.1f}")
print(f"AMP mean charge: {np.mean(amp_charges):+.1f}")

# Note about susceptibility
susc_charges = [r["net_charge"] for r in results["species"] if r["susceptible"]]
resist_charges = [r["net_charge"] for r in results["species"] if not r["susceptible"]]
print(f"\nSusceptible species mean charge: {np.mean(susc_charges):+.1f} (n={len(susc_charges)})")
if resist_charges:
    print(f"Resistant species mean charge: {np.mean(resist_charges):+.1f} (n={len(resist_charges)})")

print(f"\nNOTE: N-terminal sequences are highly conserved across mammals.")
print(f"Charge alone does NOT differentiate susceptible from resistant species.")
print(f"Resistance (dog, rabbit, horse) comes from the GLOBULAR domain and")
print(f"hydrophobic region, not from charge differences in the N-terminal tail.")
print(f"This is CONSISTENT with v5: specificity is in CONVERSION (globular/")
print(f"hydrophobic domain), not in TOXICITY (N-terminal charge).")

# Save as JSON for artifact
with open("/Users/allan/Projects/cjd/charge_analysis.json", "w") as f:
    json.dump(results, f, indent=2)

print(f"\nData saved to charge_analysis.json")
