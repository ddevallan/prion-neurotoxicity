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

Perturbation is measured LOCAL vs DISTAL within the same frame: the patch
under the peptide against the rest of the same bilayer. That is a built-in
control -- same lipids, same frame, same thermostat -- and avoids comparing
against a separate bare-membrane run.

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
SAMPLE_PS = 20.0
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
is_pops = np.concatenate([np.zeros(len(popc), bool), np.ones(len(pops), bool)])
print(f"peptide {len(pep)} | POPC {len(popc)} | POPS {len(pops)} | "
      f"P {len(phos)} | water {len(wat_o)} | sn1 carbons {len(tail)}", flush=True)
for nm_, sel in [('peptide', pep), ('P', phos), ('water', wat_o), ('lipid', lipid)]:
    if len(sel) == 0:
        sys.exit(f"ERROR: {nm_} selection empty")

# sn-1 carbons indexed per lipid residue, so the order parameter can be split
# by each lipid's own position relative to the peptide.
res_sn1 = {}
for ridx, nm, i in tail:
    res_sn1.setdefault(ridx, {})[nm] = i
_names = sorted({nm for _, nm, _ in tail}, key=lambda x: int(x[2:]))
SN1_PAIRS = list(zip(_names[:-1], _names[1:]))
print(f"sn-1 chain: {len(_names)} carbons, {len(SN1_PAIRS)} C-C vectors", flush=True)

# ---- restrained equilibration, CHARMM-GUI's schedule ---------------------
rest = mm.CustomExternalForce('k_rest*periodicdistance(x,y,z,x0,y0,z0)^2')
rest.addGlobalParameter('k_rest', 0.0)
for p_ in ('x0', 'y0', 'z0'):
    rest.addPerParticleParameter(p_)
pos0 = pdb.getPositions(asNumpy=True).value_in_unit(unit.nanometers)
# Restrain the peptide and the lipid phosphates. Holding the peptide is the
# point: it must still be unbound when production starts, or time-to-bind is
# measured from an already-attached state.
for i in np.concatenate([pep, phos]):
    rest.addParticle(int(i), pos0[i])
system.addForce(rest)

platform = mm.Platform.getPlatformByName('CUDA')
props = {'Precision': 'mixed', 'DeviceIndex': gpu_index}
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
for n, (k, ns, dt, baro_on) in enumerate(SCHEDULE, 1):
    sim.context.setParameter('k_rest', k)
    sim.integrator.setStepSize(dt*unit.picoseconds)
    sim.context.setParameter(baro.Pressure(), (1.0 if baro_on else 0.0)*unit.bar)
    sim.step(int(ns*1000/dt))
    print(f"  equil {n}/6: k={k:>5} {ns} ns @ {dt*1000:.0f} fs "
          f"{'NPT' if baro_on else 'NVT'}", flush=True)
sim.context.setParameter('k_rest', 0.0)          # release everything
sim.context.setParameter(baro.Pressure(), 1.0*unit.bar)
sim.integrator.setStepSize(DT*unit.picoseconds)
print("Restraints released; production begins with the peptide unbound.", flush=True)


# ---- perturbation metrics ------------------------------------------------
def leaflet_split(zp):
    m = zp.mean()
    return zp > m, zp <= m


def frame(pos):
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
        vals = []
        for a_, b_ in SN1_PAIRS:
            va, vb = [], []
            for ri in sel:
                ia, ib = res_sn1.get(ri, {}).get(a_), res_sn1.get(ri, {}).get(b_)
                if ia is not None and ib is not None:
                    va.append(ia)
                    vb.append(ib)
            if not va:
                continue
            v = pos[vb] - pos[va]
            cos2 = (v[:, 2]**2) / (v**2).sum(axis=1)
            vals.append(float((3*cos2.mean() - 1) / 2))
        return float(np.mean(vals)) if vals else float('nan')

    loc_res = [ri for ri, d in lip_r.items() if d < 1.5]
    dis_res = [ri for ri, d in lip_r.items() if d > 2.5]
    scd_loc = scd(loc_res) if len(loc_res) >= 3 else float('nan')
    scd_dis = scd(dis_res) if len(dis_res) >= 3 else float('nan')

    return {'z_com': float(com[2]), 'z_rel': float(com[2] - z_up),
            'min_dist': dmin, 'contacts_popc': n_pc, 'contacts_pops': n_ps,
            'thickness_local': th_loc, 'thickness_distal': th_dis,
            'thinning': th_dis - th_loc,
            'apl': float(dims[0]*dims[1]/100.0 / max(up_mask.sum(), 1)),
            'scd_local': scd_loc, 'scd_distal': scd_dis,
            'disordering': (scd_dis - scd_loc) if np.isfinite(scd_loc) and
                           np.isfinite(scd_dis) else float('nan'),
            'n_local_lipids': len(loc_res), 'n_local_P': int(loc.sum())}


n_steps = int(PROD_NS * 1000 / DT)
every = int(SAMPLE_PS / DT)
n_frames = n_steps // every
traj = []
t0 = time.time()
for i in range(n_frames):
    sim.step(every)
    pos = sim.context.getState(getPositions=True).getPositions(
        asNumpy=True).value_in_unit(unit.nanometers)
    f = frame(pos)
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
