"""Test: does PrP become more aggregation-prone at endosomal pH (4.5)?
Uses TANGO-like aggregation propensity scoring from sequence.
No external tools needed — pure analytical calculation."""
import numpy as np
import json

# PrP sequence (human, mature protein 23-230)
PRP_SEQ = "KKRPKPGGWNTGGSRYPGQGSPGGNRYPPQGGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQGGGTHSQWNKPSKPKTNMKHMAGAAAAGAVVGGLGGYMLGSAMSRPIIHFGSDYEDRYYRENMHRYPNQVYYRPMDEYSNQNNFVHDCVNITIKQHTVTTTTKGENFTETDVKMMERVVEQMCITQYERESQAYYQRGS"

# Amino acid properties relevant to aggregation
# Hydrophobicity (Kyte-Doolittle)
HYDRO = {'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'E':-3.5,'Q':-3.5,'G':-0.4,
         'H':-3.2,'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,
         'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}

# Beta-sheet propensity (Chou-Fasman)
BETA = {'A':0.83,'R':0.93,'N':0.89,'D':0.54,'C':1.19,'E':0.37,'Q':1.10,'G':0.75,
        'H':0.87,'I':1.60,'L':1.30,'K':0.74,'M':1.05,'F':1.38,'P':0.55,'S':0.75,
        'T':1.19,'W':1.37,'Y':1.47,'V':1.70}

# Charge at different pH
def charge_at_ph(seq, ph):
    """Calculate net charge per residue at given pH."""
    charge = 0
    for aa in seq:
        if aa == 'K':  # pKa 10.5
            charge += 1.0 / (1.0 + 10**(ph - 10.5))
        elif aa == 'R':  # pKa 12.5
            charge += 1.0 / (1.0 + 10**(ph - 12.5))
        elif aa == 'H':  # pKa 6.0
            charge += 1.0 / (1.0 + 10**(ph - 6.0))
        elif aa == 'D':  # pKa 3.9
            charge -= 1.0 / (1.0 + 10**(3.9 - ph))
        elif aa == 'E':  # pKa 4.1
            charge -= 1.0 / (1.0 + 10**(4.1 - ph))
    # Terminal charges
    charge += 1.0 / (1.0 + 10**(ph - 8.0))  # N-term
    charge -= 1.0 / (1.0 + 10**(3.1 - ph))  # C-term
    return charge

def aggregation_propensity(seq, ph, window=7):
    """Simple TANGO-like aggregation propensity score.
    High score = more likely to form beta-aggregates.
    Factors: hydrophobicity, beta-propensity, low charge."""
    scores = []
    for i in range(len(seq) - window + 1):
        win = seq[i:i+window]
        # Hydrophobic contribution (aggregation favored by hydrophobic regions)
        h = np.mean([HYDRO.get(aa, 0) for aa in win])
        # Beta-sheet propensity
        b = np.mean([BETA.get(aa, 1.0) for aa in win])
        # Charge penalty (high charge disfavors aggregation)
        q = abs(charge_at_ph(win, ph)) / len(win)
        # Score: high hydrophobicity + high beta + low charge = aggregation prone
        score = (h + 2) * b * max(0.1, 1.0 - q * 5)
        scores.append(score)
    return scores

# Analyze at both pH values
print("=" * 65)
print("PrP AGGREGATION PROPENSITY: pH 7.4 vs pH 4.5")
print("=" * 65)

results = {}
for ph in [7.4, 4.5]:
    scores = aggregation_propensity(PRP_SEQ, ph)
    net_charge = charge_at_ph(PRP_SEQ, ph)

    # Find aggregation-prone regions (top 10%)
    threshold = np.percentile(scores, 90)
    hotspots = []
    for i, s in enumerate(scores):
        if s >= threshold:
            hotspots.append((i+23, PRP_SEQ[i:i+7], round(s, 2)))

    # His protonation changes
    n_his = PRP_SEQ.count('H')
    his_charge_74 = n_his * (1.0 / (1.0 + 10**(7.4 - 6.0)))
    his_charge_45 = n_his * (1.0 / (1.0 + 10**(4.5 - 6.0)))

    results[f'pH_{ph}'] = {
        'net_charge': round(net_charge, 1),
        'mean_score': round(np.mean(scores), 3),
        'max_score': round(np.max(scores), 3),
        'n_hotspots': len(hotspots),
        'hotspots': hotspots[:10],
        'his_charge': round(his_charge_45 if ph == 4.5 else his_charge_74, 1),
    }

    print(f"\npH {ph}:")
    print(f"  Net charge: {net_charge:+.1f}")
    print(f"  His charge: {his_charge_45 if ph==4.5 else his_charge_74:.1f} / {n_his} His residues")
    print(f"  Mean aggregation score: {np.mean(scores):.3f}")
    print(f"  Max aggregation score: {np.max(scores):.3f}")
    print(f"  Aggregation hotspots (top 10%):")
    for pos, motif, score in hotspots[:10]:
        region = "N-terminal" if pos < 94 else "hydrophobic" if pos < 134 else "globular"
        print(f"    Res {pos}: {motif} score={score} [{region}]")

# Comparison
print(f"\n{'='*65}")
print("COMPARISON: pH 7.4 vs pH 4.5")
print(f"{'='*65}")
r74 = results['pH_7.4']
r45 = results['pH_4.5']
print(f"  Net charge:     {r74['net_charge']:+.1f} → {r45['net_charge']:+.1f} (more positive at low pH)")
print(f"  Mean agg score: {r74['mean_score']:.3f} → {r45['mean_score']:.3f} ({(r45['mean_score']/r74['mean_score']-1)*100:+.1f}%)")
print(f"  Max agg score:  {r74['max_score']:.3f} → {r45['max_score']:.3f}")
print(f"  Hotspots:       {r74['n_hotspots']} → {r45['n_hotspots']}")

change = "MORE" if r45['mean_score'] > r74['mean_score'] else "LESS"
print(f"\n  PrP is {change} aggregation-prone at endosomal pH 4.5")
print(f"  His residues become protonated (+charge): {r74['his_charge']:.1f} → {r45['his_charge']:.1f}")
print(f"  This INCREASES charge on the octarepeats (His-containing)")
print(f"  Higher charge REDUCES aggregation propensity locally")
print(f"  But INCREASES LLPS (more multivalent cation-pi interactions)")

with open('/Users/allan/Projects/cjd/ph_aggregation_results.json', 'w') as f:
    json.dump(results, f, indent=2, default=str)
print(f"\nSaved to ph_aggregation_results.json")
