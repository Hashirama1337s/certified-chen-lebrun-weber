"""prod_A_refine.py -- refinement of the production-A scalar layer (BAR-MET5): re-certify chosen sub-squares as n x n sub-sub-squares
with the SAME frozen code (cmtoric.certify_sub, degree 28, prec 192, K = 4, same AProb / a* ball) and merge conservatively.
python prod_A_refine.py SRC_DIR DST_DIR LIST [n]      LIST = k:i:j,k:i:j,...   n = 2 (default; 4 sub-sub-squares) or 4 (16)
For each listed sub-square: the scalar record 'A' of SRC_DIR/clw_e{k}_{i}_{j}.json is replaced by the merge of the n^2 children;
'B3' (rho-derivative certificates) and every other field are copied unchanged.  Written to DST_DIR (checkpoint: a DST file that already
carries A.refined == n is skipped).  Nothing in SRC_DIR is modified.
Merge (every child certifies its own closed sub-sub-square; the children cover the sub-square):
  integrals (vol, rho2, lap2, drho2, Ux_*, AUx_*): SUM of the child balls (arb), upper = sum of child uppers;
  upper bounds (sups, maxima of bounds, remainders, ric_upper, K0, Rc0, R0_0, s_sup, dS, ndS, n2dS, U_sup, S_max): MAX;
  lower bounds (ric_lower, s_inf, detm/m11/m22_inf, S_min, Gjac_inf, rho_sup_lower, bounds[..][0]): MIN; any None -> None."""
import sys, os, json, time, math
from flint import ctx, arb
import common, cm, cmtoric

DEG, PREC, H = 28, 192, 0.25                     # identical to prod_A.py
src, dst, lst = sys.argv[1], sys.argv[2], sys.argv[3]
n = int(sys.argv[4]) if len(sys.argv) > 4 else 2
os.makedirs(dst, exist_ok=True)
cm.set_degree(DEG)
prob = cmtoric.AProb('clw', os.path.join(common.SEATR, 'sol_clw_B448_M84.json'), prec=PREC)
log = lambda s: print(s, flush=True)
MAXC = ('K0', 'Rc0', 'R0_0', 's_sup', 'dS', 'ndS', 'n2dS', 'U_sup', 'S_max', 'ric_upper')
MINC = ('s_inf', 'detm_inf', 'm11_inf', 'm22_inf', 'S_min', 'ric_lower')


def mn(vals):
    return None if any(v is None for v in vals) else min(vals)


def mx(vals):
    return None if any(v is None for v in vals) else max(vals)


def merge(ch, k, s0, t0, h):
    out = dict(k=k, s0=s0, t0=t0, h=h, K=ch[0]['K'], degree=ch[0]['degree'], prec=ch[0]['prec'], refined=n)
    out['bounds'] = {kk: [min(c['bounds'][kk][0] for c in ch), max(c['bounds'][kk][1] for c in ch)] for kk in ch[0]['bounds']}
    out['rem'] = {kk: max(c['rem'][kk] for c in ch) for kk in ch[0]['rem']}
    out['ric_lower'] = mn([c['ric_lower'] for c in ch]); out['ric_upper'] = mx([c['ric_upper'] for c in ch])
    out['Gjac_inf'] = min(c['Gjac_inf'] for c in ch)
    old = ctx.prec
    ctx.prec = 256
    ints = {}
    for name in ch[0]['integrals']:
        S = arb(0)
        for c in ch:
            e = c['integrals'][name]
            S += arb(float(e['mid']), float(e['rad']))
        mf = float(S.mid())
        ints[name] = dict(mid=repr(mf), rad=(float(S.rad()) + abs(float(S.mid() - arb(mf))) + abs(mf) * 2.0 ** -52) * (1 + 1e-9) + 1e-300,
                          upper=sum(c['integrals'][name]['upper'] for c in ch) * (1 + 1e-12) + 1e-300)
    ctx.prec = old
    out['integrals'] = ints
    out['rho_sup_lower'] = min(c['rho_sup_lower'] for c in ch)
    out['sup'] = {kk: max(c['sup'][kk] for c in ch) for kk in ch[0]['sup']}
    cert = {kk: mx([c['cert'][kk] for c in ch]) for kk in MAXC}
    cert.update({kk: mn([c['cert'][kk] for c in ch]) for kk in MINC})
    out['cert'] = cert
    out['secs_build'] = sum(c['secs_build'] for c in ch); out['secs'] = sum(c['secs'] for c in ch)
    out['children'] = ch
    return out


for tok in lst.split(','):
    k, i, j = (int(x) for x in tok.split(':'))
    name = 'clw_e%d_%d_%d.json' % (k, i, j)
    fo = os.path.join(dst, name)
    if os.path.exists(fo) and json.load(open(fo))['A'].get('refined') == n:
        log('skip %s (refined %d exists)' % (fo, n)); continue
    r = json.load(open(os.path.join(src, name)))
    s0, t0 = i * H, j * H
    hc = H / n
    log('== refine e%d [%d,%d] into %d x %d (h %.5f)' % (k, i, j, n, n, hc))
    t00 = time.time()
    ch = []
    for a in range(n):
        for b in range(n):
            ctx.prec = PREC
            ch.append(cmtoric.certify_sub(prob, k, s0 + a * hc, t0 + b * hc, hc, log=log))
    A = merge(ch, k, s0, t0, H)
    o = A['integrals']; p = r['A']['integrals']
    log('  merged (%.0fs): lap2 upper %.3e (was %.3e)  rho2 upper %.3e (was %.3e)  drho2 upper %.3e (was %.3e)  vol %.6f (was %.6f)' % (
        time.time() - t00, o['lap2']['upper'], p['lap2']['upper'], o['rho2']['upper'], p['rho2']['upper'],
        o['drho2']['upper'], p['drho2']['upper'], float(o['vol']['mid']), float(p['vol']['mid'])))
    r['A'] = A
    json.dump(r, open(fo, 'w'), indent=1, default=float)
log('refine done')
