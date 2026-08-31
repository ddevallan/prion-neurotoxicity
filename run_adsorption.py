"""Driver: unbiased adsorption replicas across GPUs, pooled and compared.

Each system gets independent replicas; the pooled report gives bound fraction,
insertion depth, and POPS-vs-POPC contact enrichment with error bars.
"""
import subprocess, os, sys, json, glob, time
import numpy as np
from concurrent.futures import ThreadPoolExecutor

BASE = '/workspace/prion-neurotoxicity' if os.path.isdir('/workspace/prion-neurotoxicity') \
       else '/Users/allan/Projects/cjd'
OUT_DIR = '/workspace/adsorption' if os.path.isdir('/workspace') else \
          '/Users/allan/Projects/cjd/results_adsorption'
os.makedirs(OUT_DIR, exist_ok=True)

# system -> (n_replicas, production ns, concurrent jobs per GPU)
# The PrP system is 79k atoms versus 35k for the hexapeptides, so it gets a
# shorter run and fewer concurrent jobs per GPU.
PLAN = json.loads(os.environ.get('PLAN', json.dumps({
    'KKRPKP_water': [4, 40, 3],
    'PrP_23_93':    [4, 40, 2],
})))

try:
    n_gpus = len(subprocess.check_output(['nvidia-smi', '--list-gpus'],
                                         text=True).strip().split('\n'))
except Exception:
    n_gpus = 1

jobs = []
for sysname, (nrep, prod_ns, per_gpu) in PLAN.items():
    for r in range(nrep):
        jobs.append((sysname, r, prod_ns))
slots = n_gpus * max(p[2] for p in PLAN.values())
print(f"{n_gpus} GPU(s), {len(jobs)} jobs, {slots} slots", flush=True)
for s, (n, ns, pg) in PLAN.items():
    print(f"  {s}: {n} replicas x {ns} ns", flush=True)
print(flush=True)


def run(i):
    sysname, rep, prod_ns = jobs[i]
    gpu = str(i % n_gpus)
    tag = f"{sysname}_rep{rep}"
    if os.path.exists(f"{OUT_DIR}/{tag}.json"):
        print(f"[skip] {tag}", flush=True)
        return
    t0 = time.time()
    env = {**os.environ, 'PYTHONUNBUFFERED': '1',
           'PROD_NS': str(prod_ns), 'EQUIL_NS': '2'}
    with open(f"{OUT_DIR}/{tag}.log", 'w') as log:
        rc = subprocess.call([sys.executable, '-u', f'{BASE}/adsorption_md.py',
                              sysname, str(rep), gpu],
                             stdout=log, stderr=subprocess.STDOUT, env=env)
    print(f"[{'ok' if rc == 0 else 'FAIL'}] {tag} gpu{gpu} "
          f"({(time.time()-t0)/60:.0f} min)", flush=True)


with ThreadPoolExecutor(max_workers=slots) as ex:
    list(ex.map(run, range(len(jobs))))

# ---- pooled report --------------------------------------------------------
print("\n" + "=" * 70)
print("ADSORPTION RESULTS")
print("=" * 70)
summary = {}
for sysname in PLAN:
    files = sorted(glob.glob(f"{OUT_DIR}/{sysname}_rep*.json"))
    if not files:
        print(f"\n{sysname}: no results")
        continue
    d = [json.load(open(f)) for f in files]
    bf = np.array([x['bound_fraction'] for x in d])
    zr = np.array([x['z_rel_mean_last10ns'] for x in d])
    en = np.array([x['pops_enrichment'] for x in d if x['pops_enrichment']])
    fc = [x['first_contact_ns'] for x in d]
    fcv = np.array([f for f in fc if f is not None])

    def pm(a):
        return (f"{a.mean():+.3f} +/- {a.std(ddof=1)/np.sqrt(len(a)):.3f}"
                if len(a) > 1 else f"{a.mean():+.3f}")

    print(f"\n{sysname}  (n={len(d)})")
    print(f"  bound fraction    : {pm(bf)}")
    print(f"  bound in {len(fcv)}/{len(d)} replicas"
          + (f", first contact {fcv.mean():.1f} ns" if len(fcv) else ""))
    print(f"  z_rel last 10 ns  : {pm(zr)} nm  (negative = below phosphates)")
    if len(en):
        print(f"  POPS enrichment   : {pm(en)}x  "
              f"(1.0 = no preference over the 20% baseline)")
    summary[sysname] = {
        'n': len(d),
        'bound_fraction_mean': float(bf.mean()),
        'bound_fraction_sem': float(bf.std(ddof=1)/np.sqrt(len(bf))) if len(bf) > 1 else None,
        'n_bound_replicas': int(len(fcv)),
        'first_contact_ns_mean': float(fcv.mean()) if len(fcv) else None,
        'z_rel_mean': float(zr.mean()),
        'z_rel_sem': float(zr.std(ddof=1)/np.sqrt(len(zr))) if len(zr) > 1 else None,
        'pops_enrichment_mean': float(en.mean()) if len(en) else None,
        'pops_enrichment_sem': float(en.std(ddof=1)/np.sqrt(len(en))) if len(en) > 1 else None,
        'per_replica': {'bound_fraction': bf.tolist(), 'z_rel': zr.tolist(),
                        'enrichment': en.tolist(), 'first_contact_ns': fc},
    }

# Pairwise test between any two systems present (e.g. KKRPKP vs NNRPNP).
keys = [k for k in PLAN if k in summary]
if len(keys) >= 2:
    from itertools import combinations
    from scipy import stats
    summary['comparisons'] = {}
    for a, b in combinations(keys, 2):
        za = np.array(summary[a]['per_replica']['z_rel'])
        zb = np.array(summary[b]['per_replica']['z_rel'])
        if len(za) > 1 and len(zb) > 1:
            t, p = stats.ttest_ind(za, zb, equal_var=False)
            print(f"\n{a} vs {b}: z_rel difference "
                  f"{za.mean()-zb.mean():+.3f} nm, Welch p={p:.4f}")
            summary['comparisons'][f'{a}_vs_{b}'] = {
                'z_rel_difference': float(za.mean()-zb.mean()),
                'welch_p': float(p)}

with open(f"{OUT_DIR}/SUMMARY.json", 'w') as f:
    json.dump(summary, f, indent=2)
print(f"\nWrote {OUT_DIR}/SUMMARY.json")
