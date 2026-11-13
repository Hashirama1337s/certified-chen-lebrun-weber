"""diag_e1_13.py -- amend3 diagnosis (read-only, frozen amend2 code): the e1 [1,3] crash of production B.
(1) polar torus: which crude bounds are infinite;  (2) xchart torus, chunk by chunk (same chunks as certify.torus_max):
per box, lower bounds of |l_r| (all r, as metric_x uses them) and of |det H| (order 0), radius of x, non-finite entries of the
order-0 metric, and whether np.linalg.inv raises on the chunk."""
import sys, os, math, time
import numpy as np
import common, cmtoric, fbox, jets, jets2, certify
from common import Disc

k, i, j = 1, 1, 3
H, RHO = 0.25, 2.0
prob = cmtoric.AProb('clw', os.path.join(common.SEATR, 'sol_clw_B448_M84.json'), prec=192)
bank = fbox.PatchBank(prob, h=0.1, rp=0.4, m=120, R=1.2, kd=12)
model = jets.model_from_prob(prob, bank, cpow=0)
t0 = time.time()
do_polar = len(sys.argv) > 1 and sys.argv[1] == 'polar'
if do_polar:
    sub = certify.Sub(prob, model, k, i * H, j * H, H, 'polar')
    M, nb, rb = certify.torus_max(sub, RHO, 200)
    print('polar torus (%.0fs): infinite keys %s' % (time.time() - t0, [kk for kk, v in M.items() if not math.isfinite(v)]), flush=True)
    print('  finite maxima: ' + ' '.join('%s<=%.1e' % (kk, v) for kk, v in M.items()), flush=True)
sub = certify.Sub(prob, model, k, i * H, j * H, H, 'xchart')
U = certify.ellipse_arcs(RHO, 100, True); V = certify.ellipse_arcs(RHO, 200, False)
Uu, Vv = certify.prod_boxes(U, V)
S, T = sub.st(Uu, Vv)
print('xchart torus: %d boxes' % len(S.m), flush=True)
nbad_total = 0; nraise = 0
for c in range(0, len(S.m), 512):
    Sc, Tc = S[c:c + 512], T[c:c + 512]
    X = sub.chartV(Sc, Tc)
    lrmin = np.full(len(Sc.m), np.inf)
    for r in range(model.d):
        nu = model.nu[r]
        lr = X[0] * float(nu[0]) + X[1] * float(nu[1]) + Disc(np.full(len(Sc.m), model.lam[r][0]), np.full(len(Sc.m), model.lam[r][1]))
        lrmin = np.minimum(lrmin, np.abs(lr.m) - lr.r)
    with np.errstate(all='ignore'):
        g, _, _ = jets2.metric_x(model, X, 5)
    g0m = g.m[:, 0]; g0r = g.r[:, 0]
    Hm = g0m[:, :2, :2]; Hr = g0r[:, :2, :2]
    dH = Hm[:, 0, 0] * Hm[:, 1, 1] - Hm[:, 0, 1] * Hm[:, 1, 0]
    nonfin = ~np.isfinite(g0m).all(axis=(1, 2)) | ~np.isfinite(g0r).all(axis=(1, 2))
    G = np.moveaxis(Disc(np.moveaxis(g0m, 0, -1)).m, -1, 0)
    sing = np.array([np.linalg.matrix_rank(G[b]) < 4 if np.isfinite(G[b]).all() else True for b in range(len(G))])
    try:
        np.linalg.inv(g0m); raised = False
    except np.linalg.LinAlgError:
        raised = True
    bad = nonfin | sing | (lrmin <= 0)
    if raised:
        nraise += 1
    if bad.any() or raised:
        nbad_total += int(bad.sum())
        for b in np.nonzero(bad)[0][:6]:
            print('  chunk %d box %d: s=%s t=%s |x| rad %.2e/%.2e  min_r(|l_r| lower)=%.3e  H0 mid det=%s  H0 rad max=%.2e  nonfinite=%s singular=%s' % (
                c // 512, b, np.round(Sc.m[b], 4), np.round(Tc.m[b], 4), X[0].r[b], X[1].r[b], lrmin[b], dH[b], np.max(Hr[b]), nonfin[b], sing[b]), flush=True)
        print('  chunk %d: %d bad boxes, inv raises=%s, min |l_r| lower over chunk %.3e' % (c // 512, int(bad.sum()), raised, lrmin.min()), flush=True)
print('done (%.0fs): %d bad boxes, %d chunks where inv raises' % (time.time() - t0, nbad_total, nraise), flush=True)
