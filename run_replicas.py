"""Driver: run steered-MD replicas across GPUs and pool the statistics.

A 35k-atom PME system does not saturate a modern GPU — it uses well under 1 GB
of VRAM and leaves SMs idle. Running several replicas concurrently per GPU
therefore buys most of a linear speedup, so the slot count is (gpus x
CONCURRENT_PER_GPU) rather than one job per GPU.
"""
import subprocess, os, sys, json, glob, time
import numpy as np
from concurrent.futures import ThreadPoolExecutor

N_REPLICAS = int(os.environ.get('N_REPLICAS', 8))
CONCURRENT_PER_GPU = int(os.environ.get('CONCURRENT_PER_GPU', 3))
SYSTEMS = ['KKRPKP', 'NNRPNP']

BASE = '/workspace/prion-neurotoxicity' if os.path.isdir('/workspace/prion-neurotoxicity') \
       else '/Users/allan/Projects/cjd'
OUT_DIR = '/workspace/replicas' if os.path.isdir('/workspace') else \
          '/Users/allan/Projects/cjd/results_replicas'
os.makedirs(OUT_DIR, exist_ok=True)

try:
    n_gpus = len(subprocess.check_output(
        ['nvidia-smi', '--list-gpus'], text=True).strip().split('\n'))
except Exception:
    n_gpus = 1
n_slots = n_gpus * CONCURRENT_PER_GPU
jobs = [(s, r) for s in SYSTEMS for r in range(N_REPLICAS)]
print(f"{n_gpus} GPU(s), {CONCURRENT_PER_GPU} concurrent each = {n_slots} slots")
print(f"{len(jobs)} jobs: {N_REPLICAS} replicas x {len(SYSTEMS)} systems\n", flush=True)


def run(job_idx):
    system, rep = jobs[job_idx]
    gpu = str(job_idx % n_gpus)
    tag = f"{system}_rep{rep}"
    if os.path.exists(f"{OUT_DIR}/{tag}.json"):
        print(f"[skip] {tag} already done", flush=True)
        return
    t0 = time.time()
    with open(f"{OUT_DIR}/{tag}.log", 'w') as log:
        rc = subprocess.call(
            [sys.executable, '-u', f'{BASE}/steered_replicas.py', system, str(rep), gpu],
            stdout=log, stderr=subprocess.STDOUT,
            env={**os.environ, 'PYTHONUNBUFFERED': '1'})
    print(f"[{'ok' if rc == 0 else 'FAIL'}] {tag} gpu{gpu} "
          f"({(time.time()-t0)/60:.0f} min)", flush=True)


with ThreadPoolExecutor(max_workers=n_slots) as ex:
    list(ex.map(run, range(len(jobs))))

# ---- Pool and test --------------------------------------------------------
print("\n" + "=" * 68)
print("REPLICA STATISTICS")
print("=" * 68)

data = {}
for s in SYSTEMS:
    works, dzs = [], []
    for f in sorted(glob.glob(f"{OUT_DIR}/{s}_rep*.json")):
        d = json.load(open(f))
        works.append(d['work_kj_mol'])
        dzs.append(d['dz_mean_nm'])
    data[s] = {'work': np.array(works), 'dz': np.array(dzs)}

kT = 8.314e-3 * 303.15
summary = {'kT_kj_mol': kT, 'systems': {}}
for s in SYSTEMS:
    w = data[s]['work']
    if len(w) == 0:
        print(f"{s}: no results")
        continue
    sem = w.std(ddof=1) / np.sqrt(len(w)) if len(w) > 1 else float('nan')
    print(f"\n{s}  (n={len(w)})")
    print(f"  Work  = {w.mean():+8.2f} +/- {sem:.2f} kJ/mol (SEM)")
    print(f"          {w.mean()/kT:+8.1f} kT")
    print(f"  Range = [{w.min():+.2f}, {w.max():+.2f}]")
    print(f"  dz    = {data[s]['dz'].mean():+.4f} nm")
    summary['systems'][s] = {
        'n': len(w), 'work_mean': float(w.mean()), 'work_sem': float(sem),
        'work_std': float(w.std(ddof=1)) if len(w) > 1 else None,
        'work_kT': float(w.mean()/kT),
        'dz_mean': float(data[s]['dz'].mean()),
        'works': [float(x) for x in w],
    }

a, b = data['KKRPKP']['work'], data['NNRPNP']['work']
if len(a) > 1 and len(b) > 1:
    from scipy import stats
    t, p = stats.ttest_ind(a, b, equal_var=False)
    u, pu = stats.mannwhitneyu(a, b, alternative='two-sided')
    pooled = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
    d = (a.mean() - b.mean()) / pooled if pooled > 0 else float('nan')
    print(f"\n{'-'*68}")
    print(f"KKRPKP (+4) vs NNRPNP (0)")
    print(f"  Difference : {a.mean()-b.mean():+.2f} kJ/mol "
          f"({(a.mean()-b.mean())/kT:+.1f} kT)")
    print(f"  Welch t    : t={t:.3f}, p={p:.5f}")
    print(f"  Mann-Whitney: U={u:.1f}, p={pu:.5f}")
    print(f"  Cohen's d  : {d:.2f}")
    verdict = ("SIGNIFICANT: charge changes membrane affinity"
               if p < 0.05 else
               "NOT SIGNIFICANT: the single-trajectory difference does not survive replication")
    print(f"\n  => {verdict}")
    summary['comparison'] = {
        'difference_kj_mol': float(a.mean()-b.mean()),
        'difference_kT': float((a.mean()-b.mean())/kT),
        'welch_t': float(t), 'welch_p': float(p),
        'mannwhitney_u': float(u), 'mannwhitney_p': float(pu),
        'cohens_d': float(d), 'verdict': verdict,
    }

with open(f"{OUT_DIR}/SUMMARY.json", 'w') as f:
    json.dump(summary, f, indent=2)
print(f"\nWrote {OUT_DIR}/SUMMARY.json")
