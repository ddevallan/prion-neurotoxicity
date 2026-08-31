"""Does the peptide perturb the bilayer, or only stick to it?

This supersedes adsorption_md.py, which answered a narrower question and had
two methodological faults:

1. No restrained equilibration. It ran minimise + 2 ns of free MD, so the
   peptide was loose while the membrane was still relaxing and bound before
   production began. Reported first-contact times of 0.02 ns were an artefact
   of that: the peptide was already attached when the clock started. Here the
   peptide is held in bulk water through CHARMM-GUI's six-stage schedule
   (4000 -> 0 kJ/mol/nm^2) and released only at production, so time-to-bind
   means what it says.

2. Only scalars were saved, so nothing about the membrane could be recovered
   afterwards. The observable the model actually needs is whether binding
   thins or disorders the bilayer -- Kloda 2007 has NMDA responding to
   bilayer tension, so adsorption alone does not close the mechanism. Those
   metrics are now computed inline.

Perturbation is measured two ways, because neither control is clean alone:

  SPATIAL  local (<1.5 nm of the peptide) vs distal (>2.5 nm), same frame.
  TEMPORAL global bilayer properties before vs after the peptide binds.

Three limits are worth stating before any number is read:

1. The distal reference is not unperturbed. In the 6.6 nm hexapeptide box the
   farthest a lipid can sit from the peptide is 3.3 nm, so the distal annulus
   spans only 2.5-3.3 nm; peptide-induced perturbation typically reaches
   1-3 nm. The PrP box (8.5 nm, 4.2 nm max) is better but not clean. So this
   measures a LOWER BOUND: a positive thinning signal is real and probably
   understated, while a null result is ambiguous -- the reference may be just
   as perturbed as the test region. The temporal control exists to cover that
   case.

2. One peptide per box is P:L = 1:130 (hexapeptides) and 1:216 (PrP).
   Carpet-model thinning in experiments is usually seen at 1:100 to 1:10.
   We may be below the threshold at which a single peptide measurably thins a
   bilayer, so again a null result is weak evidence against the mechanism.

3. 1.9 ns of restrained equilibration relaxes clashes; it does not converge
   the area per lipid, which typically needs tens of ns. Area and order
   parameter are therefore tracked through production and only the last 20 ns
   is averaged, with the trace kept so convergence can be checked rather than
   assumed. A CHARMM-GUI starting structure has |S_CD| flat near 0.36 because
   the library tails are built extended; a properly equilibrated POPC bilayer
   relaxes to a plateau near 0.20 that falls toward the chain end. That
   relaxation is the check that equilibration worked.

Hydrogen mass repartitioning allows a 4 fs step (Balusek et al., JCTC 2019,
validated for CHARMM36 lipids), roughly doubling sampling per GPU-hour.

Usage: python membrane_perturbation_md.py <system> <replica> <gpu>
"""
import openmm as mm
import openmm.app as app
import openmm.unit as unit
import numpy as np
import json, sys, os, time, glob

system_name = sys.argv[1] if len(sys.argv) > 1 else 'KKRPKP_water'
replica = int(sys.argv[2]) if len(sys.argv) > 2 else 0
gpu_index = sys.argv[3] if len(sys.argv) > 3 else '0'

PATHS = {
    'KKRPKP_water': 'charmm_gui_kkrpkp_water/openmm',
    'NNRPNP_water': 'charmm_gui_nnrpnp_water/openmm',
    'PrP_23_93': 'charmm_gui_prp2393/openmm',
    'PrP_dCC1': 'charmm_gui_prp_dcc1/openmm',
}
for base in ['/workspace/prion-neurotoxicity', '/Users/allan/Projects/cjd']:
    p = os.path.join(base, PATHS[system_name])
    if os.path.isdir(p):
        os.chdir(p)
        break
else:
    sys.exit(f"ERROR: no system directory for {system_name}")

OUT_DIR = '/workspace/perturbation' if os.path.isdir('/workspace') else \
          '/Users/allan/Projects/cjd/results_perturbation'
os.makedirs(OUT_DIR, exist_ok=True)
tag = f"{system_name}_rep{replica}"
SEED = 30000 + replica * 31 + (abs(hash(system_name)) % 977)

PROD_NS = float(os.environ.get('PROD_NS', 60))
DT = 0.004                      # 4 fs, enabled by hydrogen mass repartitioning
SAMPLE_PS = float(os.environ.get('SAMPLE_PS', 20))
# TEST=1 shrinks everything so the whole path -- equilibration schedule,
# restraint release, every metric, the JSON -- can be exercised locally in
# minutes before any GPU time is spent.
TEST = os.environ.get('TEST') == '1'
if TEST:
    PROD_NS, SAMPLE_PS = float(os.environ.get('PROD_NS', 0.1)), 5.0
print(f"=== {tag} | GPU {gpu_index} | seed {SEED} | {PROD_NS} ns ===", flush=True)

psf = app.CharmmPsfFile('step5_input.psf')
pdb = app.PDBFile('step5_input.pdb')
with open('sysinfo.dat') as f:
    dims = json.load(f)['dimensions']
psf.setBox(dims[0]*unit.angstroms, dims[1]*unit.angstroms, dims[2]*unit.angstroms)

toppar = os.path.join(os.path.dirname(os.getcwd()), 'toppar')
pfiles = (sorted(glob.glob(f'{toppar}/*.rtf')) + sorted(glob.glob(f'{toppar}/*.prm')) +
          sorted(glob.glob(f'{toppar}/*.str')) + ['toppar.str'])
params = app.CharmmParameterSet(*pfiles)

system = psf.createSystem(params, nonbondedMethod=app.PME,
                          nonbondedCutoff=1.2*unit.nanometers,
                          switchDistance=1.0*unit.nanometers,
                          constraints=app.HBonds,
                          hydrogenMass=4.0*unit.amu)
baro = mm.MonteCarloMembraneBarostat(
    1.0*unit.bar, 0.0*unit.bar*unit.nanometers, 303.15*unit.kelvin,
    mm.MonteCarloMembraneBarostat.XYIsotropic,
    mm.MonteCarloMembraneBarostat.ZFree)
baro_idx = system.addForce(baro)

# ---- selections ----------------------------------------------------------
PROTEIN_RES = {'ALA','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS','HSD','HSE','HSP',
               'ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'}
WATER_RES, WATER_O = {'HOH','TIP3','WAT','SOL'}, {'O','OH2','OW'}
# Order parameter is taken on the SATURATED sn-1 palmitoyl chain, named C3x in
# CHARMM36 (C31 is the ester carbonyl, so the chain proper runs C32-C316). The
# C2x series is the sn-2 oleoyl chain, whose double bond puts a kink in the
# profile and makes S_CD harder to read; an earlier version selected it by
# mistake while the comment claimed otherwise.
SN1 = {f'C3{i}' for i in range(2, 17)}

pep, pep_bb, pep_sc, popc, pops, phos, wat_o, tail = [], [], [], [], [], [], [], []
lipid_res = {}
for a in psf.topology.atoms():
    rn, nm, i = a.residue.name, a.name, a.index
    if rn in WATER_RES:
        if nm in WATER_O:
            wat_o.append(i)
        continue
    if rn in PROTEIN_RES:
        if nm.startswith('H'):
            continue
        pep.append(i)
        (pep_bb if nm in ('N', 'CA', 'C', 'O') else pep_sc).append(i)
    elif rn in ('POPC', 'POPS'):
        if nm.startswith('H'):
            continue
        (popc if rn == 'POPC' else pops).append(i)
        lipid_res.setdefault(a.residue.index, []).append(i)
        if nm == 'P':
            phos.append(i)
        if nm in SN1:
            tail.append((a.residue.index, nm, i))

pep = np.array(pep, int); pep_bb = np.array(pep_bb, int); pep_sc = np.array(pep_sc, int)
popc = np.array(popc, int); pops = np.array(pops, int)
phos = np.array(phos, int); wat_o = np.array(wat_o, int)
lipid = np.concatenate([popc, pops])
phos_set = set(phos.tolist())
is_pops = np.concatenate([np.zeros(len(popc), bool), np.ones(len(pops), bool)])
print(f"peptide {len(pep)} | POPC {len(popc)} | POPS {len(pops)} | "
      f"P {len(phos)} | water {len(wat_o)} | sn1 carbons {len(tail)}", flush=True)
for nm_, sel in [('peptide', pep), ('P', phos), ('water', wat_o), ('lipid', lipid)]:
    if len(sel) == 0:
        sys.exit(f"ERROR: {nm_} selection empty")

# sn-1 carbons indexed per lipid residue, so the order parameter can be split
# by each lipid's own position relative to the peptide.
# S_CD is defined on the C-H vector, not on consecutive carbons. This is an
# all-atom model, so the hydrogens are present and the real quantity is
# available; using C-C vectors would give S_CC, a different number that only
# maps onto S_CD through a transformation. Build carbon -> bonded-hydrogen
# from the topology.
_bonded = {}
for _b in psf.topology.bonds():
    _bonded.setdefault(_b[0].index, []).append(_b[1])
    _bonded.setdefault(_b[1].index, []).append(_b[0])

res_sn1 = {}          # residue -> {carbon name: (C index, [H indices])}
for ridx, nm, i in tail:
    hs = [x.index for x in _bonded.get(i, []) if x.name.startswith('H')]
    if hs:
        res_sn1.setdefault(ridx, {})[nm] = (i, hs)
SN1_CARBONS = sorted({nm for d in res_sn1.values() for nm in d},
                     key=lambda x: int(x[2:]))
_nh = sum(len(h) for d in res_sn1.values() for _, h in d.values())
print(f"sn-1: {len(SN1_CARBONS)} carbons with H, {_nh} C-H vectors "
      f"over {len(res_sn1)} lipids", flush=True)

# ---- perturbation metrics ------------------------------------------------
def leaflet_split(zp):
    m = zp.mean()
    return zp > m, zp <= m


def frame(pos, box):
    P = pos[pep]
    com = P.mean(axis=0)
    zp = pos[phos, 2]
    up_mask, lo_mask = leaflet_split(zp)
    z_up, z_lo = zp[up_mask].mean(), zp[lo_mask].mean()

    # contacts, restricted to a neighbourhood so the cost stays flat
    lo_b, hi_b = P.min(axis=0) - 1.2, P.max(axis=0) + 1.2
    lc = pos[lipid]
    near = ((lc[:, 0] > lo_b[0]) & (lc[:, 0] < hi_b[0]) &
            (lc[:, 1] > lo_b[1]) & (lc[:, 1] < hi_b[1]) &
            (lc[:, 2] > lo_b[2]) & (lc[:, 2] < hi_b[2]))
    if near.any():
        d = np.linalg.norm(lc[near][None, :, :] - P[:, None, :], axis=2)
        dmin = float(d.min())
        close = (d < 0.4).any(axis=0)
        ps = is_pops[near]
        n_ps, n_pc = int(close[ps].sum()), int(close[~ps].sum())
    else:
        dmin, n_ps, n_pc = float('inf'), 0, 0

    # LOCAL vs DISTAL: upper-leaflet phosphates within 1.5 nm of the peptide
    # COM in xy, against those beyond 2.5 nm. Same frame, same bilayer.
    upper_lipids = {ri: idx for ri, idx in lipid_res.items()
                    if pos[idx][:, 2].mean() > zp.mean()}
    pxy = pos[phos][:, :2]
    r = np.linalg.norm(pxy - com[:2], axis=1)
    loc = up_mask & (r < 1.5)
    dis = up_mask & (r > 2.5)
    th_loc = float(z_up - z_lo) if loc.sum() < 3 else \
        float(pos[phos][loc, 2].mean() - z_lo)
    th_dis = float(z_up - z_lo) if dis.sum() < 3 else \
        float(pos[phos][dis, 2].mean() - z_lo)

    # Order parameter S_CD = <3cos^2(theta)-1>/2 over consecutive sn-1 carbons,
    # split local vs distal by each lipid's own xy position so it answers the
    # same question the thickness does: is the bilayer under the peptide
    # disordered relative to the rest of the same bilayer?
    lip_r = {ri: np.linalg.norm(pos[idx][:, :2].mean(axis=0) - com[:2])
             for ri, idx in upper_lipids.items()}

    def scd(sel):
        """|S_CD| averaged over the sn-1 chain. theta is the angle between each
        C-H bond and the bilayer normal (z). Reported as the absolute value,
        the convention in the lipid literature, where the plateau sits near
        0.20 for POPC."""
        vals = []
        for cname in SN1_CARBONS:
            ci, hi = [], []
            for ri in sel:
                e = res_sn1.get(ri, {}).get(cname)
                if e:
                    for h in e[1]:
                        ci.append(e[0])
                        hi.append(h)
            if not ci:
                continue
            v = pos[hi] - pos[ci]
            cos2 = (v[:, 2]**2) / (v**2).sum(axis=1)
            vals.append(float(np.abs((3*cos2 - 1) / 2).mean()))
        return float(np.mean(vals)) if vals else float('nan')

    # Radial shells rather than one local/distal split. Franco et al. 2022
    # (BP100 at P:L 1:128, PMID 35139324) found shell analysis essential:
    # global averages masked a local effect that was plainly there. Huang's
    # two-state model puts P/L* near 1/50, so at our 1:130 the GLOBAL thinning
    # expected is under 0.5 A, inside thermal noise -- only the first shells
    # can show anything.
    shells = [(0.0, 1.0), (1.0, 2.0), (2.0, 3.0), (3.0, 99.0)]
    prof = []
    for lo_r, hi_r in shells:
        sel = [ri for ri, d in lip_r.items() if lo_r <= d < hi_r]
        if len(sel) >= 3:
            zs = np.concatenate([[pos[i][2] for i in lipid_res[ri]
                                  if i in set(phos.tolist())] for ri in sel]) \
                 if False else np.array([pos[[i for i in lipid_res[ri]
                                              if i in phos_set]][:, 2].mean()
                                         for ri in sel])
            prof.append({'r': (lo_r + min(hi_r, 4.0)) / 2, 'n': len(sel),
                         'thickness': float(zs.mean() - z_lo),
                         'scd': scd(sel)})
        else:
            prof.append({'r': (lo_r + min(hi_r, 4.0)) / 2, 'n': len(sel),
                         'thickness': float('nan'), 'scd': float('nan')})

    loc_res = [ri for ri, d in lip_r.items() if d < 1.5]
    dis_res = [ri for ri, d in lip_r.items() if d > 2.5]
    scd_loc = scd(loc_res) if len(loc_res) >= 3 else float('nan')
    scd_dis = scd(dis_res) if len(dis_res) >= 3 else float('nan')
    # Whole-leaflet values, for the temporal control and for the equilibration
    # check described at the top of this file.
    scd_all = scd(list(upper_lipids))
    th_all = float(z_up - z_lo)

    return {'z_com': float(com[2]), 'z_rel': float(com[2] - z_up),
            'min_dist': dmin, 'contacts_popc': n_pc, 'contacts_pops': n_ps,
            'thickness_local': th_loc, 'thickness_distal': th_dis,
            'thinning': th_dis - th_loc,
            # Live box vectors: under NPT the box breathes, so the initial
            # dimensions would give a constant, wrong area per lipid.
            'apl': float(box[0][0]*box[1][1]*100.0 / max(up_mask.sum(), 1)),
            'scd_local': scd_loc, 'scd_distal': scd_dis,
            'scd_global': scd_all, 'thickness_global': th_all,
            'disordering': (scd_dis - scd_loc) if np.isfinite(scd_loc) and
                           np.isfinite(scd_dis) else float('nan'),
            'n_local_lipids': len(loc_res), 'n_local_P': int(loc.sum()),
            'radial': prof}


# ---- restrained equilibration, CHARMM-GUI's schedule ---------------------
pos0 = pdb.getPositions(asNumpy=True).value_in_unit(unit.nanometers)

# The peptide is held in 3D: it must still be unbound when production starts,
# or time-to-bind is measured from an already-attached state.
rest = mm.CustomExternalForce('k_rest*periodicdistance(x,y,z,x0,y0,z0)^2')
rest.addGlobalParameter('k_rest', 0.0)
for p_ in ('x0', 'y0', 'z0'):
    rest.addPerParticleParameter(p_)
for i in pep:
    rest.addParticle(int(i), pos0[i])
system.addForce(rest)

# Phosphates are held in z only, as CHARMM-GUI's membrane_restraint.str does.
# Pinning them in x and y as well would stop lateral diffusion and prevent the
# bilayer from relaxing its area under the barostat -- the one thing this
# equilibration exists to allow.
zrest = mm.CustomExternalForce('k_zrest*(z-z0)^2')
zrest.addGlobalParameter('k_zrest', 0.0)
zrest.addPerParticleParameter('z0')
for i in phos:
    zrest.addParticle(int(i), [pos0[i][2]])
system.addForce(zrest)

def pick_platform():
    for name in (os.environ.get('PLATFORM') or 'CUDA,OpenCL,CPU').split(','):
        try:
            pl = mm.Platform.getPlatformByName(name.strip())
            if name.strip() == 'CUDA':
                return pl, {'Precision': 'mixed', 'DeviceIndex': gpu_index}
            if name.strip() == 'OpenCL':
                # OpenCLPlatformIndex must be given explicitly: without it the
                # automatic device search fails on Apple Silicon even though
                # the GPU is usable ("No compatible OpenCL platform").
                # Apple Silicon GPUs have no fp64, so 'mixed' and 'double'
                # both fail here; only 'single' initialises. That is fine for
                # a smoke test, which checks code paths rather than numerics --
                # production on CUDA keeps mixed precision.
                return pl, {'Precision': 'single',
                            'OpenCLPlatformIndex':
                                os.environ.get('OPENCL_PLATFORM', '0')}
            pl.setPropertyDefaultValue('Threads',
                                       os.environ.get('THREADS', str(os.cpu_count() or 4)))
            return pl, {}
        except Exception:
            continue
    raise RuntimeError('no usable platform')


platform, props = pick_platform()
print(f"Platform: {platform.getName()}", flush=True)
integ = mm.LangevinMiddleIntegrator(303.15*unit.kelvin, 1.0/unit.picosecond,
                                    0.001*unit.picoseconds)
integ.setRandomNumberSeed(SEED)
sim = app.Simulation(psf.topology, system, integ, platform, props)
sim.context.setPositions(pdb.positions)

print("Minimising...", flush=True)
sim.minimizeEnergy(maxIterations=10000)
sim.context.setVelocitiesToTemperature(303.15*unit.kelvin, SEED)

# (force constant kJ/mol/nm^2, ns, timestep ps, barostat on)
SCHEDULE = [(4000, 0.125, 0.001, False), (2000, 0.125, 0.001, False),
            (1000, 0.125, 0.001, True),  (500, 0.5, 0.002, True),
            (200, 0.5, 0.002, True),     (50, 0.5, 0.004, True)]
if TEST:      # same six stages, same force constants, much shorter
    _div = float(os.environ.get('TEST_DIV', 40))
    SCHEDULE = [(k, ns/_div, dt, b) for k, ns, dt, b in SCHEDULE]
for n, (k, ns, dt, baro_on) in enumerate(SCHEDULE, 1):
    sim.context.setParameter('k_rest', k)
    sim.context.setParameter('k_zrest', k / 4.0)   # lipids relax faster than the peptide
    sim.integrator.setStepSize(dt*unit.picoseconds)
    sim.context.setParameter(baro.Pressure(), (1.0 if baro_on else 0.0)*unit.bar)
    sim.step(int(ns*1000/dt))
    print(f"  equil {n}/6: k={k:>5} {ns} ns @ {dt*1000:.0f} fs "
          f"{'NPT' if baro_on else 'NVT'}", flush=True)
sim.context.setParameter('k_rest', 0.0)          # release everything
sim.context.setParameter('k_zrest', 0.0)
sim.context.setParameter(baro.Pressure(), 1.0*unit.bar)
sim.integrator.setStepSize(DT*unit.picoseconds)
_st = sim.context.getState(getPositions=True)
_p = _st.getPositions(asNumpy=True).value_in_unit(unit.nanometers)
_f = frame(_p, _st.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometers))
print(f"After equilibration: |S_CD| global {_f['scd_global']:.3f} "
      f"(start ~0.36; a relaxed POPC bilayer sits near 0.20), "
      f"thickness {_f['thickness_global']*10:.1f} A, "
      f"APL {_f['apl']:.1f} A^2, peptide dmin {_f['min_dist']:.2f} nm", flush=True)
if _f['min_dist'] < 0.4:
    print("  WARNING: peptide is already in contact before production; "
          "time-to-bind will not be meaningful.", flush=True)
print("Restraints released; production begins.", flush=True)


n_steps = int(PROD_NS * 1000 / DT)
every = int(SAMPLE_PS / DT)
n_frames = n_steps // every
traj = []
t0 = time.time()
for i in range(n_frames):
    sim.step(every)
    st = sim.context.getState(getPositions=True)
    pos = st.getPositions(asNumpy=True).value_in_unit(unit.nanometers)
    box = st.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometers)
    f = frame(pos, box)
    f['t_ns'] = (i + 1) * SAMPLE_PS / 1000
    traj.append(f)
    if i % 100 == 0:
        el = time.time() - t0
        print(f"  {f['t_ns']:6.1f} ns  dmin={f['min_dist']:.2f} "
              f"thin={f['thinning']*10:+.2f} A  dS={f['disordering']:+.3f}  "
              f"{f['t_ns']/(el/86400) if el else 0:.0f} ns/day", flush=True)

# ---- summary with block-averaged error -----------------------------------
t = np.array([f['t_ns'] for f in traj])
late = t > t.max() - 20


def blocked(key, n_blocks=5):
    v = np.array([f[key] for f in traj])[late]
    if len(v) < n_blocks:
        return float(v.mean()), float('nan')
    b = np.array([x.mean() for x in np.array_split(v, n_blocks)])
    return float(b.mean()), float(b.std(ddof=1)/np.sqrt(n_blocks))


dmin = np.array([f['min_dist'] for f in traj])
bound = dmin < 0.4
thin_m, thin_e = blocked('thinning')
scd_m, scd_e = blocked('disordering')
cpc = np.array([f['contacts_popc'] for f in traj])[late]
cps = np.array([f['contacts_pops'] for f in traj])[late]
tot = cpc + cps
frac = len(pops) / (len(pops) + len(popc))
enrich = float(np.mean((cps[tot > 0]/tot[tot > 0])/frac)) if (tot > 0).any() else None

result = {
    'system': system_name, 'replica': replica, 'seed': SEED,
    'prod_ns': PROD_NS, 'dt_fs': DT*1000, 'hmr': True,
    'equilibration': 'CHARMM-GUI 6-stage restrained, peptide held unbound',
    'bound_fraction': float(bound.mean()),
    'first_contact_ns': float(t[bound][0]) if bound.any() else None,
    'thinning_nm': thin_m, 'thinning_sem': thin_e,
    'disordering': scd_m, 'disordering_sem': scd_e,
    'pops_enrichment': enrich,
    'contacts_popc_mean': float(cpc.mean()), 'contacts_pops_mean': float(cps.mean()),
    'trajectory': traj,
}
with open(f"{OUT_DIR}/{tag}.json", 'w') as f:
    json.dump(result, f, indent=2)
print(f"\n{tag}: bound {bound.mean()*100:.0f}%, first contact "
      f"{result['first_contact_ns']} ns")
print(f"  thinning {thin_m*10:+.2f} +/- {thin_e*10:.2f} A  (local vs distal)")
print(f"  disordering (S_CD distal - local) {scd_m:+.4f} +/- {scd_e:.4f}")
print(f"  POPS enrichment {enrich if enrich else float('nan'):.2f}x")
print(f"  {(time.time()-t0)/60:.0f} min")
