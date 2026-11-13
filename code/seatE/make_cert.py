"""make_cert.py -- assemble seatB2's cert.json (NOTES-B2 sec. 7 schema) from the per-sub-square certificates of production A
(prod_A.py: cmtoric scalar layer + certify_b3) and, when available, production B (tensor sups: kappa, K_1, K_2, Rc_1..Rc_3, R0_1, R0_2).
python make_cert.py PRODA_DIR OUT_cert.json [PRODB_DIR]
Rules: every number written is a certified bound in the stated direction, valid for all a in the a* ball, rounded OUTWARD
(lower bounds x (1 - 1e-12), upper bounds x (1 + 1e-12) + 1e-300).  Fields not yet produced are written as "MISSING" and the file
is then marked certified = false.  lambda1 is NOT a seatE quantity: it is copied from seatB2/cert_sampled.json with its source noted.
Global checks: all 5 x 16 sub-squares present with status ok; sum of the volume pieces encloses (2 pi)^2 Area(P)."""
import sys, os, json, glob, math
from flint import arb, arb_mat, ctx
import common, cmtoric

ctx.prec = 256
d, outf = sys.argv[1], sys.argv[2]
pb = None
if len(sys.argv) > 3:
    RB = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(sys.argv[3], 'clwB_e*_*_*.json')))]
    okB = [r for r in RB if r.get('status') == 'ok']
    if len(okB) == 80:
        pb = {key: max(r['cert_a'][src] for r in okB) for key, src in
              (('kappa', 'kappa'), ('K1', 'K1'), ('K2', 'K2'), ('Rc1', 'Rc1'), ('Rc2', 'Rc2'), ('Rc3', 'Rc3'), ('R0_1', 'R0c1'), ('R0_2', 'R0c2'))}
    print('production B: %d/%d sub-squares ok' % (len(okB), len(RB)))
files = sorted(glob.glob(os.path.join(d, 'clw_e*_*_*.json')))
R = [json.load(open(f)) for f in files]
have = {(r['k'], r['i'], r['j']) for r in R}
need = {(k, i, j) for k in range(5) for i in range(4) for j in range(4)}
missing_subs = sorted(need - have)
up = lambda x: repr(float(x) * (1 + 1e-12) + 1e-300) if x >= 0 else repr(float(x) * (1 - 1e-12))
lo = lambda x: repr(float(x) * (1 - 1e-12) - 1e-300) if x >= 0 else repr(float(x) * (1 + 1e-12))
A = [r['A'] for r in R]
B3 = [r['B3'] for r in R]
bad = [(r['k'], r['i'], r['j']) for r in R if r['B3'].get('status') != 'ok' or r['A']['cert']['ric_lower'] is None or r['A']['cert']['ric_upper'] is None]
prob = cmtoric.AProb('clw', None, prec=256)
P = prob.P
V = (2 * arb.pi()) ** 2 * prob.area
# Gram of e = (1, x1, x2) in dx (exact Gauss on the bilinear elements)
GP = cmtoric.gram_P(P)
GPi = GP.inv()
fro = arb(0)
for a in range(3):
    for b in range(3):
        fro += GPi[a, b] ** 2
GPinv_norm = fro.sqrt()


def ball(entry):
    return arb(float(entry['mid']), float(entry['rad']))


def ssum(name):
    tot = arb(0)
    for a in A:
        tot += ball(a['integrals'][name])
    return tot
G = {}
for j, key in ((1, 'Ux'), (2, 'AUx')):
    M = arb_mat(3, 3)
    for (a_, b_) in ((0, 0), (0, 1), (1, 1)):
        v = ssum('%s_%d%d' % (key, a_, b_))
        M[a_ + 1, b_ + 1] = v; M[b_ + 1, a_ + 1] = v
    G[j] = M
kj = {}
certified_ok_kj = True
for j in (1, 2):
    tr = arb(0)
    for a_ in range(3):
        for b_ in range(3):
            tr += GPi[b_, a_] * G[j][a_, b_]
    if not (tr > 0):                                  # BAR-MET5-amend1: fail-safe instead of 0
        certified_ok_kj = False
        kj[j] = arb(float('inf'))
    else:
        kj[j] = tr.sqrt()                             # lambda_max(G_P^-1 G_j) <= tr (G_j PSD, G_P PD)
kappa_A = 1 + kj[1] + kj[2]
vol_sum = ssum('vol')
cert = {}
cert['label'] = 'seatE production A (arb Chebyshev models, a* ball) + certify_b3' + (' + production B' if pb else ' ; production-B fields MISSING')
mu = min(a['cert']['ric_lower'] for a in A)
Rp = max(a['cert']['ric_upper'] for a in A)
cert['mu'] = lo(mu)
cert['Rm'] = '0' if mu > 0 else 'MISSING'
cert['Rp'] = up(Rp); cert['Rop'] = up(max(Rp, abs(mu)))
K0 = max(a['cert']['K0'] for a in A)
R00 = max(a['cert']['R0_0'] for a in A)
S_max = max(a['cert']['S_max'] for a in A); S_min = min(a['cert']['S_min'] for a in A)
rho_sup = max(a['sup']['rho'] for a in A)
drho = max(a['sup']['drho'] for a in A)
dS = max(a['cert']['dS'] for a in A); ndS = max(a['cert']['ndS'] for a in A); n2dS = max(a['cert']['n2dS'] for a in A)
rj = [max(b['b3']['rho%d' % j]['sup'] for b in B3 if b.get('status') == 'ok') for j in range(4)] if not bad else [float('nan')] * 4
cert['kappa'] = up(pb['kappa']) if pb else 'MISSING'
cert['K'] = [up(K0)] + ([up(pb['K1']), up(pb['K2'])] if pb else ['MISSING', 'MISSING'])
cert['Rc'] = [up(pb['Rc%d' % q]) for q in (1, 2, 3)] if pb else ['MISSING'] * 3
cert['R0'] = [up(R00)] + ([up(pb['R0_1']), up(pb['R0_2'])] if pb else ['MISSING', 'MISSING'])
# |nabla^q ds| <= |nabla^q dS| + |nabla^{q+1} rho|  (s = S + rho)
cert['S_d'] = [up(dS + max(drho, rj[1])), up(ndS + rj[2]), up(n2dS + rj[3])]
cert['S_max'] = up(S_max); cert['S_min'] = lo(S_min)
cert['rho_sup'] = up(rho_sup)
gS = (prob.Sabc[1] ** 2 + prob.Sabc[2] ** 2).sqrt()
cert['gradS_x'] = up(float(gS.upper()))
sampled = json.load(open(os.path.join(common.HERE, '..', 'seatB2', 'cert_sampled.json')))
_g4 = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'seatG45', 'G4_eps_m48.json')))
_l1 = float(_g4['RESULT']['lambda1_M84_lower'])
import math as _m
_l1_down = _m.floor(_l1 * 1e11) / 1e11          # rounded DOWN to 12 significant digits
assert _l1_down <= _l1
cert['lambda1'] = repr(_l1_down)
cert['lambda1_note'] = 'BAR-MET5-amend1: certified lower bound lambda1_M84_lower from seatG45/G4_eps_m48.json, rounded down'
cert['V'] = [lo(float(V.lower())), up(float(V.upper()))]
cert['kappa_A'] = up(float(kappa_A.upper()))
cert['GPinv_norm'] = up(float(GPinv_norm.upper()))
cert['GP_diag'] = [up(float(GP[0, 0].upper())), up(float(GP[1, 1].upper())), up(float(GP[2, 2].upper()))]
cert['U_sup'] = up(max(a['cert']['U_sup'] for a in A))
cert['rho_L2'] = up(math.sqrt(sum(a['integrals']['rho2']['upper'] for a in A)))
cert['lap_rho_L2'] = up(math.sqrt(sum(a['integrals']['lap2']['upper'] for a in A)))
cert['rho_H3'] = up(float(V.upper().sqrt()) * sum(rj)) if not bad else 'MISSING'
cert['rho_H3_note'] = '||rho||_H3 = sum_j ||nabla^j rho||_L2(M) <= sqrt(V) sum_j sup|nabla^j rho| (certify_b3 interpolation certificates)'
cert['w_prime'] = sampled['w_prime']; cert['d_target'] = sampled['d_target']
cert['K1_admissible'] = dict(detm_inf=lo(min(a['cert']['detm_inf'] for a in A)), m11_inf=lo(min(a['cert']['m11_inf'] for a in A)),
                             m22_inf=lo(min(a['cert']['m22_inf'] for a in A)), Gjac_inf=lo(min(a['Gjac_inf'] for a in A)))
cert['checks'] = dict(subsquares=len(R), missing_subsquares=missing_subs, bad_subsquares=bad,
                      vol_sum=vol_sum.str(15, radius=True), V=V.str(15, radius=True), vol_consistent=bool(vol_sum.overlaps(V)),
                      k1=kj[1].str(10, radius=True), k2=kj[2].str(10, radius=True))
complete = certified_ok_kj and not missing_subs and not bad and pb is not None and bool(vol_sum.overlaps(V)) and mu > 0 and \
    all(float(v) > 0 for v in [cert['K1_admissible'][x] for x in ('detm_inf', 'm11_inf', 'm22_inf', 'Gjac_inf')])
cert['certified'] = bool(complete)
json.dump(cert, open(outf, 'w'), indent=1)
print(json.dumps(cert, indent=1))
