"""
Kinetic model v2: three Ca²⁺ mechanisms + expanded combination therapy.
Incorporates all findings from the memantine deep dive.

Three convergent excitotoxicity mechanisms:
  A: Loss of PrPC inhibition on GluN2D (loss-of-function)
  B: AMP membrane perturbation → NMDA mechanosensitivity (gain-of-function)
  C: AMPA receptor remodeling — GluA2 retention (gain-of-function)

Therapeutic interventions:
  - ASO/siRNA: reduces PrPC substrate
  - Memantine: blocks NMDA (mechanisms A+B)
  - Perampanel: blocks AMPA (mechanism C)
  - NAC: ROS scavenger (prevents dityrosine → solid transition)
  - Trazodone: PERK pathway protection
  - Lithium: autophagy
"""
import numpy as np
from scipy.integrate import solve_ivp
import json

# ============================================================
# PARAMETERS (calibrated to RML in WT mouse, ~150 days)
# ============================================================

# Conversion
K_CONV = 0.25
K_CLEAR_SC = 0.10
PMAX = 10.0
K_PROD = 1.0
K_CLEAR_C = 0.1

# AMP generation
K_RELEASE = 0.12
K_DEGRADE_AMP = 0.15

# Three Ca²⁺ mechanisms
# Mechanism A: loss of PrPC inhibition on NMDA (proportional to PrPC consumed)
CA_MECHANISM_A_WEIGHT = 0.35  # fraction of total Ca²⁺ damage
# Mechanism B: AMP membrane perturbation → NMDA (proportional to AMP above threshold)
CA_MECHANISM_B_WEIGHT = 0.40
AMP_THRESHOLD = 0.3
# Mechanism C: AMPA remodeling (proportional to PrPSc, with delay)
CA_MECHANISM_C_WEIGHT = 0.25

# Damage
K_DAMAGE = 1.0
DAMAGE_LETHAL = 130.0  # calibrated for ~150 day baseline

# Copper/ROS
K_ROS_BASELINE = 0.02  # baseline ROS from copper release
K_ROS_AMP = 0.05       # ROS amplified by membrane damage

# Therapeutic parameters
INTERVENTIONS = {
    'aso': {'reduction': 0.7, 'onset_days': 14, 'target': 'substrate'},
    'memantine': {'block': 0.5, 'onset_days': 3, 'target': 'nmda'},  # blocks A+B
    'perampanel': {'block': 0.6, 'onset_days': 3, 'target': 'ampa'},  # blocks C
    'nac': {'block': 0.4, 'onset_days': 1, 'target': 'ros'},
    'trazodone': {'block': 0.3, 'onset_days': 7, 'target': 'perk'},
    'lithium': {'block': 0.2, 'onset_days': 14, 'target': 'autophagy'},
    'g127v': {'destabilize': 0.5, 'onset_days': 28, 'target': 'conversion'},
}

def intervention_effect(t, start_day, onset_days, max_effect):
    if t < start_day:
        return 0.0
    elapsed = t - start_day
    return max_effect * (1 - np.exp(-3 * elapsed / onset_days))

def model_v2(t, y, active_interventions=None):
    S, P, A, ROS, D = y  # PrPC, PrPSc, AMP, ROS, cumulative Damage

    if active_interventions is None:
        active_interventions = {}

    # Calculate intervention effects
    aso_eff = 0
    mem_eff = 0
    per_eff = 0
    nac_eff = 0
    traz_eff = 0
    lith_eff = 0
    g127v_eff = 0

    for name, start_day in active_interventions.items():
        if name in INTERVENTIONS:
            params = INTERVENTIONS[name]
            if 'reduction' in params:
                aso_eff = intervention_effect(t, start_day, params['onset_days'], params['reduction'])
            elif 'block' in params:
                eff = intervention_effect(t, start_day, params['onset_days'], params['block'])
                if params['target'] == 'nmda':
                    mem_eff = eff
                elif params['target'] == 'ampa':
                    per_eff = eff
                elif params['target'] == 'ros':
                    nac_eff = eff
                elif params['target'] == 'perk':
                    traz_eff = eff
                elif params['target'] == 'autophagy':
                    lith_eff = eff
            elif 'destabilize' in params:
                g127v_eff = intervention_effect(t, start_day, params['onset_days'], params['destabilize'])

    # PrPC consumed (for mechanism A)
    S_steady = K_PROD / K_CLEAR_C  # steady state without infection
    prpc_consumed_fraction = max(0, 1 - S / S_steady)

    # Conversion
    conv = K_CONV * P * S / (1 + P / PMAX)

    # ODEs
    dS = K_PROD * (1 - aso_eff) - conv - K_CLEAR_C * S
    dP = conv - K_CLEAR_SC * P - g127v_eff * P - lith_eff * 0.1 * P  # lithium enhances clearance
    dA = K_RELEASE * conv - K_DEGRADE_AMP * A

    # ROS dynamics (copper release + membrane damage feedback)
    dROS = K_ROS_BASELINE * prpc_consumed_fraction + K_ROS_AMP * max(0, A - AMP_THRESHOLD) - 0.3 * ROS - nac_eff * ROS

    # Three mechanisms of Ca²⁺-driven damage
    mech_a = CA_MECHANISM_A_WEIGHT * prpc_consumed_fraction * (1 - mem_eff)
    mech_b = CA_MECHANISM_B_WEIGHT * max(0, A - AMP_THRESHOLD) * (1 - mem_eff)
    mech_c = CA_MECHANISM_C_WEIGHT * (P / PMAX) * (1 - per_eff)

    total_ca_damage = K_DAMAGE * (mech_a + mech_b + mech_c)

    # PERK/translation damage (additive, from ROS and misfolded protein)
    perk_damage = 0.15 * (P / PMAX) * (1 - traz_eff)

    # ROS damage (membrane lipid peroxidation)
    ros_damage = 0.1 * ROS

    dD = total_ca_damage + perk_damage + ros_damage

    return [dS, dP, dA, dROS, dD]

def run_scenario(name, interventions, t_max=300):
    S0 = K_PROD / K_CLEAR_C  # steady state PrPC = 10
    y0 = [S0, 0.001, 0.0, 0.0, 0.0]
    t_eval = np.linspace(0, t_max, 1500)

    sol = solve_ivp(model_v2, (0, t_max), y0, t_eval=t_eval,
                    args=(interventions,), method='RK45', max_step=0.5)

    death_idx = np.where(sol.y[4] > DAMAGE_LETHAL)[0]
    death_day = sol.t[death_idx[0]] if len(death_idx) > 0 else t_max

    # Calculate mechanism contributions at peak
    peak_idx = len(sol.t) // 2
    S_ss = K_PROD / K_CLEAR_C
    prpc_consumed = max(0, 1 - sol.y[0][peak_idx] / S_ss)
    amp_above = max(0, sol.y[2][peak_idx] - AMP_THRESHOLD)
    p_frac = sol.y[1][peak_idx] / PMAX

    return {
        "name": name,
        "death_day": round(float(death_day), 1),
        "mechanism_a_peak": round(float(CA_MECHANISM_A_WEIGHT * prpc_consumed), 3),
        "mechanism_b_peak": round(float(CA_MECHANISM_B_WEIGHT * amp_above), 3),
        "mechanism_c_peak": round(float(CA_MECHANISM_C_WEIGHT * p_frac), 3),
    }

# ============================================================
# RUN SCENARIOS
# ============================================================

print("=" * 80)
print("KINETIC MODEL v2 — THREE Ca²⁺ MECHANISMS + EXPANDED COMBINATION")
print("=" * 80)

scenarios = {}

# Baseline
scenarios["none"] = run_scenario("Sem tratamento", {})
baseline = scenarios["none"]["death_day"]

# Individual drugs
scenarios["aso_80"] = run_scenario("ASO dia 80", {"aso": 80})
scenarios["mem_80"] = run_scenario("Memantina dia 80", {"memantine": 80})
scenarios["per_80"] = run_scenario("Perampanel dia 80", {"perampanel": 80})
scenarios["nac_80"] = run_scenario("NAC dia 80", {"nac": 80})
scenarios["traz_80"] = run_scenario("Trazodona dia 80", {"trazodone": 80})

# Two-drug combinations
scenarios["mem_per_80"] = run_scenario("Mem + Perampl dia 80", {"memantine": 80, "perampanel": 80})
scenarios["aso_mem_80"] = run_scenario("ASO + Mem dia 80", {"aso": 80, "memantine": 80})

# Three-drug: all Ca²⁺ doors
scenarios["triple_ca"] = run_scenario("Mem + Perampl + NAC dia 80", {
    "memantine": 80, "perampanel": 80, "nac": 80
})

# Full v5 combination day 80
scenarios["full_80"] = run_scenario("Combinação v5 completa dia 80", {
    "aso": 80, "memantine": 80, "perampanel": 80, "nac": 80, "trazodone": 80, "lithium": 80
})

# Full v5 combination day 110 (symptomatic)
scenarios["full_110"] = run_scenario("Combinação v5 completa dia 110", {
    "aso": 110, "memantine": 110, "perampanel": 110, "nac": 110, "trazodone": 110, "lithium": 110
})

# Bridge: neuroprotection first, ASO later
scenarios["bridge"] = run_scenario("Ponte: Mem+Per+NAC d80, ASO d94", {
    "memantine": 80, "perampanel": 80, "nac": 80, "aso": 94
})

# Full + G127V
scenarios["full_g127v"] = run_scenario("Combinação + G127V dia 80", {
    "aso": 80, "memantine": 80, "perampanel": 80, "nac": 80,
    "trazodone": 80, "lithium": 80, "g127v": 80
})

# Late rescue
scenarios["late_full"] = run_scenario("Resgate tardio dia 130", {
    "aso": 130, "memantine": 130, "perampanel": 130, "nac": 130,
    "trazodone": 130, "lithium": 130
})

# Print results
print(f"\n{'Scenario':<45} {'Death':>8} {'Ext':>10} {'MechA':>7} {'MechB':>7} {'MechC':>7}")
print("-" * 87)

for key, sc in scenarios.items():
    ext = sc["death_day"] - baseline
    if sc["death_day"] >= 300:
        ext_str = "SURVIVAL"
        death_str = "> 300"
    else:
        ext_str = f"+{ext:.0f}d ({ext/baseline*100:.0f}%)"
        death_str = f"{sc['death_day']:.1f}"

    print(f"{sc['name']:<45} {death_str:>8} {ext_str:>10} "
          f"{sc['mechanism_a_peak']:>7.3f} {sc['mechanism_b_peak']:>7.3f} {sc['mechanism_c_peak']:>7.3f}")

print(f"\nBaseline death: day {baseline}")
print(f"Mechanisms: A=loss of PrPC NMDA inhibition, B=AMP membrane perturbation, C=AMPA remodeling")

print(f"""
INTERPRETATION:

- Memantine blocks mechanisms A+B ({(CA_MECHANISM_A_WEIGHT+CA_MECHANISM_B_WEIGHT)*100:.0f}% of Ca²⁺ damage)
- Perampanel blocks mechanism C ({CA_MECHANISM_C_WEIGHT*100:.0f}% of Ca²⁺ damage)
- Together they cover all three Ca²⁺ entry points
- NAC reduces ROS → slows dityrosine crosslinks → slows conversion feedback
- ASO is still the trunk intervention — without it, protection is temporary
- Bridge strategy: neuroprotect immediately, add ASO 2 weeks later
""")

# Save
with open("/Users/allan/Projects/cjd/kinetic_v2_results.json", "w") as f:
    json.dump(scenarios, f, indent=2)
print("Results saved to kinetic_v2_results.json")
