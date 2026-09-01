"""Long equilibration: same script, longer production.

The |S_CD| was stuck at 0.36 in all 16 runs (60 ns). This run extends to
200 ns to determine (a) whether S_CD converges and (b) what the converged
baseline properties are.

Only two systems: bare_66 (the baseline) and KKRPKP_water rep1 (the replica
that was bound 100% and showed -1.0 A thinning). If the thinning persists
with a converged bilayer, the result stands. If it vanishes, it was an
equilibration artefact.
"""
import subprocess, os, sys, json, time

BASE = '/workspace/prion-neurotoxicity' if os.path.isdir('/workspace/prion-neurotoxicity') \
       else '/Users/allan/Projects/cjd'
OUT_DIR = '/workspace/long_equil' if os.path.isdir('/workspace') else \
          '/Users/allan/Projects/cjd/results_long_equil'
os.makedirs(OUT_DIR, exist_ok=True)

PROD_NS = float(os.environ.get('PROD_NS', 200))

JOBS = [
    ('bare_66', 0),
    ('KKRPKP_water', 1),
]

try:
    n_gpus = len(subprocess.check_output(['nvidia-smi', '--list-gpus'],
                                         text=True).strip().split('\n'))
except Exception:
    n_gpus = 1

print(f"{n_gpus} GPU(s), {len(JOBS)} jobs, {PROD_NS} ns each", flush=True)
from concurrent.futures import ThreadPoolExecutor


def run(i):
    system, rep = JOBS[i]
    gpu = str(i % n_gpus)
    tag = f"{system}_rep{rep}"
    out = f"{OUT_DIR}/{tag}.json"
    if os.path.exists(out):
        print(f"[skip] {tag}", flush=True)
        return
    t0 = time.time()
    env = {**os.environ, 'PYTHONUNBUFFERED': '1', 'PROD_NS': str(PROD_NS)}
    if 'PLATFORM' not in os.environ:
        env['PLATFORM'] = 'CUDA'
    with open(f"{OUT_DIR}/{tag}.log", 'w') as log:
        rc = subprocess.call([sys.executable, '-u',
                              f'{BASE}/membrane_perturbation_md.py',
                              system, str(rep), gpu],
                             stdout=log, stderr=subprocess.STDOUT, env=env)
    # Move the result to our output dir
    src = f"{BASE}/results_perturbation/{tag}.json"
    if os.path.exists(src):
        import shutil
        shutil.move(src, out)
    print(f"[{'ok' if rc == 0 else 'FAIL'}] {tag} gpu{gpu} "
          f"({(time.time()-t0)/60:.0f} min)", flush=True)


with ThreadPoolExecutor(max_workers=n_gpus) as ex:
    list(ex.map(run, range(len(JOBS))))

# Report
import numpy as np
print("\n" + "=" * 70)
print("LONG EQUILIBRATION RESULTS")
print("=" * 70)
for system, rep in JOBS:
    tag = f"{system}_rep{rep}"
    f = f"{OUT_DIR}/{tag}.json"
    if not os.path.exists(f):
        print(f"\n{tag}: MISSING"); continue
    d = json.load(open(f))
    tr = d['trajectory']
    t = np.array([x['t_ns'] for x in tr])
    scd = np.array([x['scd_global'] for x in tr], dtype=float)
    apl = np.array([x['apl'] for x in tr], dtype=float)
    th = np.array([x['thickness_global'] for x in tr], dtype=float)
    ok = np.isfinite(scd) & np.isfinite(apl) & np.isfinite(th)
    # Report in 50 ns windows
    print(f"\n{tag}:")
    print(f"  {'window':>12}  {'|S_CD|':>8}  {'APL (A²)':>10}  {'thick (A)':>10}")
    for lo in range(0, int(PROD_NS), max(1, int(PROD_NS/6))):
        hi = lo + 50
        mask = ok & (t >= lo) & (t < hi)
        if mask.sum() < 10: continue
        print(f"  {lo:>3}-{hi:>3} ns    {scd[mask].mean():.4f}    {apl[mask].mean():>8.1f}    {th[mask].mean()*10:>8.1f}")

print(f"\nWrote results to {OUT_DIR}/")
