"""Driver: membrane-perturbation replicas across GPUs, with bare-membrane controls.

One replica per GPU (validated by bench_concurrency.py: 841 ns/day solo vs
825 at 2 concurrent — stacking does not help on RTX 4090 with these systems).
"""
import subprocess, os, sys, json, glob, time
from concurrent.futures import ThreadPoolExecutor

BASE = '/workspace/prion-neurotoxicity' if os.path.isdir('/workspace/prion-neurotoxicity') \
       else '/Users/allan/Projects/cjd'
OUT_DIR = '/workspace/perturbation' if os.path.isdir('/workspace') else \
          '/Users/allan/Projects/cjd/results_perturbation'
os.makedirs(OUT_DIR, exist_ok=True)

# system -> (n_replicas, production ns)
PLAN = json.loads(os.environ.get('PLAN', json.dumps({
    'KKRPKP_water': [3, 60],
    'NNRPNP_water': [3, 60],
    'PrP_23_93':    [3, 60],
    'PrP_dCC1':     [3, 60],
    'bare_66':      [2, 60],
    'bare_85':      [2, 60],
})))

try:
    n_gpus = len(subprocess.check_output(['nvidia-smi', '--list-gpus'],
                                         text=True).strip().split('\n'))
except Exception:
    n_gpus = 1

jobs = [(s, r) for s, (nrep, _) in PLAN.items() for r in range(nrep)]
print(f"{n_gpus} GPU(s), {len(jobs)} jobs (1 per GPU, sequential waves)", flush=True)
for s, (n, ns) in PLAN.items():
    print(f"  {s}: {n} replicas x {ns} ns", flush=True)
print(flush=True)


def run(job_idx):
    system, rep = jobs[job_idx]
    gpu = str(job_idx % n_gpus)
    tag = f"{system}_rep{rep}"
    if os.path.exists(f"{OUT_DIR}/{tag}.json"):
        print(f"[skip] {tag}", flush=True)
        return
    prod_ns = PLAN[system][1]
    t0 = time.time()
    env = {**os.environ, 'PYTHONUNBUFFERED': '1',
           'PROD_NS': str(prod_ns), 'PLATFORM': 'CUDA'}
    with open(f"{OUT_DIR}/{tag}.log", 'w') as log:
        rc = subprocess.call([sys.executable, '-u',
                              f'{BASE}/membrane_perturbation_md.py',
                              system, str(rep), gpu],
                             stdout=log, stderr=subprocess.STDOUT, env=env)
    print(f"[{'ok' if rc == 0 else 'FAIL'}] {tag} gpu{gpu} "
          f"({(time.time()-t0)/60:.0f} min)", flush=True)


with ThreadPoolExecutor(max_workers=n_gpus) as ex:
    list(ex.map(run, range(len(jobs))))

# ---- pooled report --------------------------------------------------------
print("\n" + "=" * 70)
print("PERTURBATION RESULTS")
print("=" * 70)

import numpy as np

summary = {}
for sysname in PLAN:
    files = sorted(glob.glob(f"{OUT_DIR}/{sysname}_rep*.json"))
    if not files:
        print(f"\n{sysname}: no results")
        continue
    d = [json.load(open(f)) for f in files]

    def stat(key):
        v = np.array([x[key] for x in d if x.get(key) is not None], float)
        v = v[np.isfinite(v)]
        if len(v) == 0:
            return None, None
        return float(v.mean()), float(v.std(ddof=1)/np.sqrt(len(v))) if len(v) > 1 else None

    bf_m, bf_e = stat('bound_fraction')
    th_m, th_e = stat('thinning_nm')
    ds_m, ds_e = stat('disordering')
    en = [x['pops_enrichment'] for x in d if x.get('pops_enrichment')]
    en_m = float(np.mean(en)) if en else None
    en_e = float(np.std(en, ddof=1)/np.sqrt(len(en))) if len(en) > 1 else None
    fc = [x['first_contact_ns'] for x in d]

    print(f"\n{sysname}  (n={len(d)})")
    if bf_m is not None:
        print(f"  bound fraction     : {bf_m:.3f}" + (f" +/- {bf_e:.3f}" if bf_e else ""))
    if th_m is not None:
        print(f"  thinning (loc-dis) : {th_m*10:+.2f}" + (f" +/- {th_e*10:.2f} A" if th_e else " A"))
    if ds_m is not None:
        print(f"  disordering dS_CD  : {ds_m:+.4f}" + (f" +/- {ds_e:.4f}" if ds_e else ""))
    if en_m is not None:
        print(f"  POPS enrichment    : {en_m:.2f}" + (f" +/- {en_e:.2f}x" if en_e else "x"))
    print(f"  first contact      : {fc}")

    summary[sysname] = {
        'n': len(d), 'bound_fraction': bf_m, 'bound_fraction_sem': bf_e,
        'thinning_nm': th_m, 'thinning_sem': th_e,
        'disordering': ds_m, 'disordering_sem': ds_e,
        'pops_enrichment': en_m, 'pops_enrichment_sem': en_e,
        'first_contact_ns': fc,
    }

with open(f"{OUT_DIR}/SUMMARY.json", 'w') as f:
    json.dump(summary, f, indent=2)
print(f"\nWrote {OUT_DIR}/SUMMARY.json")
