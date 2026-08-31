"""How many simulations should share one GPU?

Every run so far used 2-3 jobs per GPU by guess. The question has direct cost
consequences: if 4 concurrent jobs give more aggregate throughput than 2, the
same rental buys twice the sampling. If they give less, we have been paying
for contention.

Measures wall-clock ns/day per job and the aggregate, for 1, 2, 3 and 4
concurrent jobs pinned to a single GPU, on the real system.
"""
import json, os, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor

BASE = os.path.dirname(os.path.abspath(__file__))
SYSTEM = os.environ.get('BENCH_SYSTEM', 'KKRPKP_water')
STEPS = int(os.environ.get('BENCH_STEPS', 4000))
LEVELS = [int(x) for x in os.environ.get('BENCH_LEVELS', '1,2,3,4').split(',')]

WORKER = r'''
import openmm as mm, openmm.app as app, openmm.unit as u
import json, glob, os, sys, time
sysname, gpu = sys.argv[1], sys.argv[2]
steps = int(sys.argv[3])
PATHS = {'KKRPKP_water':'charmm_gui_kkrpkp_water/openmm',
         'NNRPNP_water':'charmm_gui_nnrpnp_water/openmm',
         'PrP_23_93':'charmm_gui_prp2393/openmm'}
for b in ['/workspace/prion-neurotoxicity','/Users/allan/Projects/cjd']:
    p = os.path.join(b, PATHS[sysname])
    if os.path.isdir(p):
        os.chdir(p); break
psf = app.CharmmPsfFile('step5_input.psf'); pdb = app.PDBFile('step5_input.pdb')
dims = json.load(open('sysinfo.dat'))['dimensions']
psf.setBox(*[d*u.angstroms for d in dims[:3]])
tp = '../toppar'
pf = (sorted(glob.glob(tp+'/*.rtf')) + sorted(glob.glob(tp+'/*.prm')) +
      sorted(glob.glob(tp+'/*.str')) + ['toppar.str'])
system = psf.createSystem(app.CharmmParameterSet(*pf), nonbondedMethod=app.PME,
                          nonbondedCutoff=1.2*u.nanometers,
                          switchDistance=1.0*u.nanometers,
                          constraints=app.HBonds, hydrogenMass=4.0*u.amu)
system.addForce(mm.MonteCarloMembraneBarostat(
    1.0*u.bar, 0.0*u.bar*u.nanometers, 303.15*u.kelvin,
    mm.MonteCarloMembraneBarostat.XYIsotropic,
    mm.MonteCarloMembraneBarostat.ZFree))
sim = app.Simulation(psf.topology, system,
    mm.LangevinMiddleIntegrator(303.15*u.kelvin, 1/u.picosecond, 0.004*u.picoseconds),
    mm.Platform.getPlatformByName('CUDA'),
    {'Precision':'mixed','DeviceIndex':gpu})
sim.context.setPositions(pdb.positions)
sim.minimizeEnergy(maxIterations=100)
sim.step(500)                      # warm up: kernel compilation, allocation
t = time.time(); sim.step(steps); el = time.time() - t
print(json.dumps({'ns_per_day': steps*0.004/1000/(el/86400)}))
'''

worker_path = os.path.join(BASE, '_bench_worker.py')
with open(worker_path, 'w') as f:
    f.write(WORKER)


def one(_):
    out = subprocess.run([sys.executable, worker_path, SYSTEM, '0', str(STEPS)],
                         capture_output=True, text=True)
    for line in out.stdout.strip().splitlines()[::-1]:
        try:
            return json.loads(line)['ns_per_day']
        except Exception:
            continue
    print(out.stderr[-400:], file=sys.stderr)
    return None


print(f"system {SYSTEM}, {STEPS} steps at 4 fs, single GPU\n")
print(f"{'jobs':<6}{'ns/day each':>14}{'aggregate':>12}{'efficiency':>12}")
results = {}
solo = None
for n in LEVELS:
    with ThreadPoolExecutor(max_workers=n) as ex:
        vals = [v for v in ex.map(one, range(n)) if v]
    if not vals:
        print(f"{n:<6}{'FAILED':>14}")
        continue
    each, agg = sum(vals)/len(vals), sum(vals)
    solo = solo or agg
    results[n] = {'each': each, 'aggregate': agg, 'efficiency': agg/solo}
    print(f"{n:<6}{each:>14.0f}{agg:>12.0f}{agg/solo:>11.2f}x")

best = max(results, key=lambda k: results[k]['aggregate']) if results else None
if best:
    print(f"\nBest aggregate throughput: {best} concurrent job(s) "
          f"({results[best]['aggregate']:.0f} ns/day, "
          f"{results[best]['efficiency']:.2f}x over one)")
with open(os.path.join(BASE, 'results_perturbation', 'BENCH_CONCURRENCY.json'), 'w') as f:
    json.dump({'system': SYSTEM, 'steps': STEPS, 'results': results, 'best': best}, f, indent=2)
os.remove(worker_path)
