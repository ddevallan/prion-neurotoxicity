"""
Sensitivity analysis of the v2 kinetic model.
Varies each parameter ±50% and measures impact on survival.
Identifies which intervention targets matter most.
"""
import numpy as np
from scipy.integrate import solve_ivp
import json

# Import model parameters from v2
K_CONV = 0.25
K_CLEAR_SC = 0.10
PMAX = 10.0
K_PROD = 1.0
K_CLEAR_C = 0.1
K_RELEASE = 0.12
K_DEGRADE_AMP = 0.15
CA_MECHANISM_A_WEIGHT = 0.35
CA_MECHANISM_B_WEIGHT = 0.40
AMP_THRESHOLD = 0.3
CA_MECHANISM_C_WEIGHT = 0.25
K_DAMAGE = 1.0
DAMAGE_LETHAL = 130.0
K_ROS_BASELINE = 0.02
K_ROS_AMP = 0.05

def model(t, y, params):
    S, P, A, ROS, D = y
    p = params

    S_ss = p['K_PROD'] / p['K_CLEAR_C']
    prpc_consumed = max(0, 1 - S / S_ss)
    conv = p['K_CONV'] * P * S / (1 + P / p['PMAX'])

    dS = p['K_PROD'] - conv - p['K_CLEAR_C'] * S
    dP = conv - p['K_CLEAR_SC'] * P
    dA = p['K_RELEASE'] * conv - p['K_DEGRADE_AMP'] * A
    dROS = p['K_ROS_BASELINE'] * prpc_consumed + p['K_ROS_AMP'] * max(0, A - p['AMP_THRESHOLD']) - 0.3 * ROS

    mech_a = p['CA_A'] * prpc_consumed
    mech_b = p['CA_B'] * max(0, A - p['AMP_THRESHOLD'])
    mech_c = p['CA_C'] * (P / p['PMAX'])

    dD = p['K_DAMAGE'] * (mech_a + mech_b + mech_c) + 0.15 * (P / p['PMAX']) + 0.1 * ROS

    return [dS, dP, dA, dROS, dD]

def get_death_day(params, t_max=400):
    S0 = params['K_PROD'] / params['K_CLEAR_C']
    y0 = [S0, 0.001, 0.0, 0.0, 0.0]
    sol = solve_ivp(model, (0, t_max), y0, t_eval=np.linspace(0, t_max, 2000),
                    args=(params,), method='RK45', max_step=0.5)
    death_idx = np.where(sol.y[4] > params['DAMAGE_LETHAL'])[0]
    return sol.t[death_idx[0]] if len(death_idx) > 0 else t_max

BASE_PARAMS = {
    'K_CONV': K_CONV, 'K_CLEAR_SC': K_CLEAR_SC, 'PMAX': PMAX,
    'K_PROD': K_PROD, 'K_CLEAR_C': K_CLEAR_C,
    'K_RELEASE': K_RELEASE, 'K_DEGRADE_AMP': K_DEGRADE_AMP,
    'CA_A': CA_MECHANISM_A_WEIGHT, 'CA_B': CA_MECHANISM_B_WEIGHT,
    'AMP_THRESHOLD': AMP_THRESHOLD, 'CA_C': CA_MECHANISM_C_WEIGHT,
    'K_DAMAGE': K_DAMAGE, 'DAMAGE_LETHAL': DAMAGE_LETHAL,
    'K_ROS_BASELINE': K_ROS_BASELINE, 'K_ROS_AMP': K_ROS_AMP,
}

PARAM_LABELS = {
    'K_CONV': 'Taxa de conversão',
    'K_CLEAR_SC': 'Clearance de PrPSc',
    'K_PROD': 'Produção de PrPC',
    'K_RELEASE': 'Liberação de AMP',
    'K_DEGRADE_AMP': 'Degradação de AMP',
    'CA_A': 'Peso mech A (perda PrPC→NMDA)',
    'CA_B': 'Peso mech B (AMP→membrana→NMDA)',
    'AMP_THRESHOLD': 'Limiar de AMP',
    'CA_C': 'Peso mech C (AMPA remodeling)',
    'K_DAMAGE': 'Taxa de dano',
    'K_ROS_BASELINE': 'ROS basal (cobre)',
    'K_ROS_AMP': 'ROS amplificado (membrana)',
}

# Parameters that represent THERAPEUTIC TARGETS
THERAPEUTIC_MAP = {
    'K_CONV': 'G127V / anti-conversão',
    'K_PROD': 'ASO / siRNA',
    'K_RELEASE': 'α-clivagem / ADAM10',
    'K_DEGRADE_AMP': 'HS analog / neutralização AMP',
    'CA_A': 'Memantina (parcial)',
    'CA_B': 'Memantina (parcial)',
    'CA_C': 'Perampanel',
    'K_ROS_BASELINE': 'NAC / scavenger ROS',
    'K_ROS_AMP': 'NAC / scavenger ROS',
    'AMP_THRESHOLD': '(propriedade intrínseca)',
}

print("=" * 80)
print("SENSITIVITY ANALYSIS — KINETIC MODEL v2")
print("Which parameters matter most for survival?")
print("=" * 80)

baseline = get_death_day(BASE_PARAMS)
print(f"\nBaseline death: day {baseline:.1f}")

results = []

# Vary each parameter ±50%
for param_name in PARAM_LABELS:
    base_val = BASE_PARAMS[param_name]

    # -50%
    p_low = dict(BASE_PARAMS)
    p_low[param_name] = base_val * 0.5
    death_low = get_death_day(p_low)

    # +50%
    p_high = dict(BASE_PARAMS)
    p_high[param_name] = base_val * 1.5
    death_high = get_death_day(p_high)

    # Sensitivity = (Δ death / Δ param) normalized
    sensitivity = (death_high - death_low) / baseline
    direction = "↑ param = lives longer" if death_high > death_low else "↑ param = dies sooner"

    results.append({
        'param': param_name,
        'label': PARAM_LABELS[param_name],
        'therapeutic': THERAPEUTIC_MAP.get(param_name, '—'),
        'base': base_val,
        'death_low': round(death_low, 1),
        'death_high': round(death_high, 1),
        'sensitivity': round(abs(sensitivity), 3),
        'direction': direction,
    })

# Sort by sensitivity
results.sort(key=lambda x: x['sensitivity'], reverse=True)

print(f"\n{'Rank':>4} {'Parameter':<35} {'−50%':>8} {'Base':>8} {'+50%':>8} {'|ΔS|':>7} {'Therapeutic target':<25}")
print("-" * 100)

for i, r in enumerate(results):
    d_low = f"{r['death_low']:.0f}" if r['death_low'] < 400 else ">400"
    d_high = f"{r['death_high']:.0f}" if r['death_high'] < 400 else ">400"
    print(f"{i+1:>4} {r['label']:<35} {d_low:>8} {baseline:>8.0f} {d_high:>8} {r['sensitivity']:>7.3f} {r['therapeutic']:<25}")

print(f"""
INTERPRETATION:

Parameters ranked by impact on survival (|ΔS| = normalized sensitivity):

Top tier (most impactful):
  → These are the targets that matter most therapeutically.
  → A 50% change in these parameters has the largest effect on survival.

Middle tier:
  → Relevant but secondary targets.

Bottom tier:
  → Changing these parameters by 50% barely affects survival.
  → Not worth targeting therapeutically.

KEY QUESTION: Is reducing CONVERSION (ASO/G127V) more impactful than
blocking TOXICITY (memantine/perampanel)?

If K_CONV and K_PROD dominate → ASO is the priority, neuroprotection is adjunct.
If CA_A/B/C and K_DAMAGE dominate → neuroprotection might be sufficient alone.
""")

# Tornado chart data
print(f"\n{'=' * 60}")
print("TORNADO CHART (visual)")
print(f"{'=' * 60}")
print(f"{'Parameter':<30} {'Effect':>50}")
print("-" * 80)

max_sens = max(r['sensitivity'] for r in results)
for r in results:
    bar_len = int(r['sensitivity'] / max_sens * 40)
    bar = "█" * bar_len
    print(f"{r['label']:<30} {bar} {r['sensitivity']:.3f}")

# Save
with open("/Users/allan/Projects/cjd/sensitivity_results.json", "w") as f:
    json.dump(results, f, indent=2)
print(f"\nResults saved to sensitivity_results.json")
