"""prod_A.py -- production A (scalar layer + rho-derivative certificates), one vertex element per process.
python prod_A.py ELEMENT OUTDIR [case] [subs]    case: clw (default) | square | page | null30;  subs: 'all' (16 sub-squares of side 1/4)
or a comma list like 0:0,2:2 (sub-square indices i:j -> [i/4,(i+1)/4] x [j/4,(j+1)/4]).
Per sub-square: (A) cmtoric.certify_sub (arb Chebyshev models, a* ball); (B3) certify.certify_b3 (|nabla^j rho|, j <= 3; float-disc
torus crude bounds + first-kind arb nodes).  One JSON per sub-square in OUTDIR (checkpoint: existing files are skipped)."""
import sys, os, json, time, math
import numpy as np
from flint import ctx, arb
import common, cm, cmtoric, fbox, jets, certify

k = int(sys.argv[1]); outdir = sys.argv[2]
case = sys.argv[3] if len(sys.argv) > 3 else 'clw'
subs = sys.argv[4] if len(sys.argv) > 4 else 'all'
DEG, PREC, H, RHO = 28, 192, 0.25, 2.0
os.makedirs(outdir, exist_ok=True)
cm.set_degree(DEG)
SR = common.SEATR
src = {'clw': ('clw', os.path.join(SR, 'sol_clw_B448_M84.json')), 'page': ('page', os.path.join(SR, 'sol_page_B512_M56.json')),
       'square': ('square', None), 'triangle': ('triangle', None), 'null30': ('clw', None)}
kind, sol = src[case]
if case == 'null30':
    sys.path.insert(0, common.SEATA4)
    import clw_metric as Wm
    d30 = np.load(os.path.join(common.HERE, '..', 'seatW', 'bg_bach_M30.npz'))
    prob = cmtoric.AProb(kind, None, prec=PREC, Cfloat=Wm.coef_matrix(d30['a'], int(d30['M'])))
else:
    prob = cmtoric.AProb(kind, sol, prec=PREC)
bank = fbox.PatchBank(prob, h=0.1, rp=0.4, m=120, R=1.2, kd=12) if prob.Cint is not None else None   # amend1: rp 0.2 -> 0.4 (analyticity cover reaches 0.342)
model = jets.model_from_prob(prob, bank, cpow=0)
log = lambda s: print(s, flush=True)
ij = [(i, j) for i in range(4) for j in range(4)] if subs == 'all' else [tuple(int(x) for x in p.split(':')) for p in subs.split(',')]
log('prod_A element %d case %s: %d sub-squares, degree %d, prec %d, a* ball %s' % (k, case, len(ij), DEG, PREC,
    prob.astar.str(30, radius=True) if hasattr(prob, 'astar') else '-'))
t00 = time.time()
for (i, j) in ij:
    s0, t0 = i * H, j * H
    fn = os.path.join(outdir, '%s_e%d_%d_%d.json' % (case, k, i, j))
    if os.path.exists(fn):
        log('skip %s (exists)' % fn); continue
    ctx.prec = PREC
    rA = cmtoric.certify_sub(prob, k, s0, t0, H, log=log)
    chart = 'cart' if (i == 0 and j == 0) else ('rodA' if i == 0 else ('rodB' if j == 0 else 'polar'))
    sub = certify.Sub(prob, model, k, s0, t0, H, chart)
    rB = certify.certify_b3(sub, RHO, 100 if chart != 'polar' else 200, log=log)
    if rB.get('status', '').startswith('FAIL torus') and chart == 'polar':
        # amend1: near the element centroid the polar torus boxes inflate to inf; the symplectic chart (all l_r large there,
        # x(s,t) polynomial) gives finite crude bounds.  Same certificate, different chart for the crude bound only.
        log('  b3: polar torus infinite -> retry in xchart')
        sub = certify.Sub(prob, model, k, s0, t0, H, 'xchart')
        rB = certify.certify_b3(sub, RHO, 200, log=log)
        if rB.get('status', '').startswith('FAIL torus'):
            # amend3 (BAR-MET3-amend3): xchart also infinite -> 4 quarter sub-squares (side H/2), each with the same rule
            # (polar, xchart retry) and settings; merged conservatively: sup|nabla^j rho| = max over the quarters.
            log('  b3: xchart torus also infinite -> split into 4 quarter sub-squares')
            pieces = []
            for (a, b) in ((0, 0), (0, 1), (1, 0), (1, 1)):
                qs0, qt0 = s0 + a * (H / 2), t0 + b * (H / 2)
                log('  b3 -- quarter [%d,%d] s0=%.4f t0=%.4f h=%.4f' % (a, b, qs0, qt0, H / 2))
                q = certify.certify_b3(certify.Sub(prob, model, k, qs0, qt0, H / 2, 'polar'), RHO, 200, log=log)
                if q.get('status', '').startswith('FAIL torus'):
                    log('  b3: polar torus infinite -> retry in xchart')
                    q = certify.certify_b3(certify.Sub(prob, model, k, qs0, qt0, H / 2, 'xchart'), RHO, 200, log=log)
                pieces.append(q)
            if all(q.get('status') == 'ok' and all(math.isfinite(q['b3']['rho%d' % jj]['sup']) for jj in range(4)) for q in pieces):
                b3 = {'rho%d' % jj: dict(sup=max(q['b3']['rho%d' % jj]['sup'] for q in pieces),
                                         interp_err=max(q['b3']['rho%d' % jj]['interp_err'] for q in pieces),
                                         node_max=max(q['b3']['rho%d' % jj]['node_max'] for q in pieces)) for jj in range(4)}
                rB = dict(k=k, s0=s0, t0=t0, h=H, chart='quarters', status='ok', b3=b3, n_rho=max(q['n_rho'] for q in pieces), quarters=pieces)
                log('  b3 quarters merged: sup|nabla^j rho| <= %s' % ' '.join('%.2e' % b3['rho%d' % jj]['sup'] for jj in range(4)))
            else:
                rB = dict(k=k, s0=s0, t0=t0, h=H, chart='quarters', status='FAIL quarter', quarters=pieces)
    json.dump(dict(case=case, k=k, i=i, j=j, A=rA, B3=rB), open(fn, 'w'), indent=1, default=float)
log('element %d done in %.0fs' % (k, time.time() - t00))
