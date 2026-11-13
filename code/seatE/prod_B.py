"""prod_B.py -- production B (tensor sups by the float-disc engine + Chebyshev interpolation), one vertex element per process.
python prod_B.py ELEMENT OUTDIR [subs]   (subs: 'all' or 'i:j,...').  Same chart rule and checkpointing as prod_A.py.
amend3 (BAR-MET3-amend3): if the xchart retry is also infinite, the sub-square is split into 4 quarter sub-squares (side H/2), each
certified with the same rule (polar, xchart retry), rho, n_arc, tol_rel, n_max; merged conservatively (max of sups, min of infs).
amend4 (BAR-MET3-amend4): both triggers fire only on TARGETED crude bounds; status ok only with finite targeted values."""
import sys, os, json, time, math
from flint import ctx
import common, cmtoric, fbox, jets, certify

k = int(sys.argv[1]); outdir = sys.argv[2]; subs = sys.argv[3] if len(sys.argv) > 3 else 'all'
H, RHO = 0.25, 2.0
os.makedirs(outdir, exist_ok=True)
prob = cmtoric.AProb('clw', os.path.join(common.SEATR, 'sol_clw_B448_M84.json'), prec=192)
bank = fbox.PatchBank(prob, h=0.1, rp=0.4, m=120, R=1.2, kd=12)   # amend1 (as BAR-MET2-amend1): analyticity cover reaches 0.342
model = jets.model_from_prob(prob, bank, cpow=0)
log = lambda s: print(s, flush=True)
ij = [(i, j) for i in range(4) for j in range(4)] if subs == 'all' else [tuple(int(x) for x in p.split(':')) for p in subs.split(',')]
INF_KEYS = ('s_inf', 'ric_lower')                    # lower bounds: merged by min; every other cert_a entry is an upper bound (max)


TARGET = ('K1', 'K2', 'Rc1', 'Rc2', 'Rc3', 'R0c1', 'R0c2', 'Q0', 's')   # crude bounds behind K1, K2, Rc1-3, R0_1-2, the Ricci bounds
NEED = ['K1', 'K2', 'kappa', 'Rc1', 'Rc2', 'Rc3', 'R0c1', 'R0c2', 'ric_lower']   # amend4: targeted certificate values
TRACES = ('T8', 'T16', 'T32', 'T64')                                       # kappa = min over these: needs ANY one finite


def finite(M):
    """amend4 (BAR-MET3-amend4): the chart switch and the split fire only when a TARGETED crude bound is infinite.
    Untargeted overflow (e.g. T64 = inf with T8 finite; K0, K3, S0-2, ...) no longer triggers; those entries are certified as before
    (inf stays inf, never reported as finite)."""
    return all(math.isfinite(M[kk]) for kk in TARGET) and any(math.isfinite(M[kk]) for kk in TRACES)


def plan_chart(s0, t0, h, chart, narc):
    """the amend1/amend2 rule, unchanged: plan in the chart; polar with an infinite crude bound -> re-plan in xchart."""
    sub = certify.Sub(prob, model, k, s0, t0, h, chart)
    pl = certify.plan(sub, RHO, narc, tol_rel=0.02, n_max=200, log=log)
    if pl['ok'] and chart == 'polar' and not finite(pl['M']):
        # amend1: polar torus crude bounds infinite near the element centroid -> symplectic chart (as BAR-MET2-amend1)
        log('  polar torus infinite -> retry in xchart')
        chart = 'xchart'
        sub = certify.Sub(prob, model, k, s0, t0, h, chart)
        pl = certify.plan(sub, RHO, narc, tol_rel=0.02, n_max=200, log=log)
    return sub, pl, chart


def cert_piece(sub, pl, chart, narc):
    if not pl['ok']:
        return dict(status='FAIL analyticity', den=pl['den'])
    r = certify.certify(sub, RHO, narc, pl['n'], do_b=False, log=log, den=pl['den'], M=pl['M'], nbox=pl['nbox'], rbox=pl['rbox'])
    return dict(status=r['status'], chart=chart, cert_a=r['cert_a'], n=r['n'], M=r['M'])


def merge(pieces):
    """conservative merge of the quarter certificates: upper bounds by max, lower bounds by min; None (no certificate) or NaN
    in any piece -> None / FAIL."""
    if any(p['status'] != 'ok' for p in pieces):
        return dict(status='FAIL quarter (%s)' % ', '.join(p['status'] for p in pieces), chart='quarters', quarters=pieces)
    keys = pieces[0]['cert_a'].keys()
    ca = {}
    for kk in keys:
        vals = [p['cert_a'].get(kk) for p in pieces]
        if any(v is None for v in vals):
            ca[kk] = None
        else:
            ca[kk] = min(vals) if kk in INF_KEYS else max(vals)
    nan = any(isinstance(v, float) and math.isnan(v) for p in pieces for v in p['cert_a'].values())
    M = {kk: max(p['M'][kk] for p in pieces) for kk in pieces[0]['M']}
    # a split exists to produce usable bounds: the merged sub-square is 'ok' only if every NEEDED value is finite (kappa_T16..T64 may be
    # inf, as before, since kappa is their min); otherwise it is reported as a failure, never as an 'ok' with inf bounds
    need = NEED
    unusable = [kk for kk in need if ca.get(kk) is None or not math.isfinite(ca[kk])]
    st = 'FAIL nan' if nan else ('FAIL quarter infinite (%s)' % ','.join(unusable) if unusable else 'ok')
    return dict(status=st, chart='quarters', cert_a=ca, n=max(p['n'] for p in pieces), M=M, quarters=pieces)


t00 = time.time()
for (i, j) in ij:
    fn = os.path.join(outdir, 'clwB_e%d_%d_%d.json' % (k, i, j))
    if os.path.exists(fn):
        log('skip %s' % fn); continue
    chart = 'cart' if (i == 0 and j == 0) else ('rodA' if i == 0 else ('rodB' if j == 0 else 'polar'))
    narc = 100 if chart != 'polar' else 200
    log('== e%d [%d,%d] chart %s' % (k, i, j, chart))
    sub, pl, ch = plan_chart(i * H, j * H, H, chart, narc)
    if pl['ok'] and ch == 'xchart' and not finite(pl['M']):
        # amend3: the xchart retry is also infinite -> 4 quarter sub-squares, same rule and settings, conservative merge
        log('  xchart torus also infinite -> split into 4 quarter sub-squares')
        pieces = []
        for (a, b) in ((0, 0), (0, 1), (1, 0), (1, 1)):
            s0, t0 = i * H + a * (H / 2), j * H + b * (H / 2)
            log('  -- quarter [%d,%d] s0=%.4f t0=%.4f h=%.4f chart %s' % (a, b, s0, t0, H / 2, chart))
            qs, qpl, qch = plan_chart(s0, t0, H / 2, chart, narc)
            pc = cert_piece(qs, qpl, qch, narc)
            pc.update(s0=s0, t0=t0, h=H / 2)
            pieces.append(pc)
        rec = merge(pieces)
    else:
        rec = cert_piece(sub, pl, ch, narc)
    if rec['status'] == 'FAIL analyticity':
        json.dump(dict(k=k, i=i, j=j, status='FAIL analyticity', den=rec['den']), open(fn, 'w'), default=float); continue
    if rec['status'] == 'ok' and any(isinstance(v, float) and math.isnan(v) for v in rec['cert_a'].values()):
        rec['status'] = 'FAIL nan'                                            # amend3: never report a NaN as a bound
    if rec['status'] == 'ok':                                                 # amend4: never 'ok' with an unusable targeted value
        bad_t = [kk for kk in NEED if rec['cert_a'].get(kk) is None or not math.isfinite(rec['cert_a'][kk])]
        if bad_t:
            rec['status'] = 'FAIL targeted infinite (%s)' % ','.join(bad_t)
    out = dict(k=k, i=i, j=j, status=rec['status'], chart=rec['chart'], cert_a=rec.get('cert_a'), n=rec.get('n'), M=rec.get('M'))
    if 'quarters' in rec:
        out['quarters'] = rec['quarters']
    json.dump(out, open(fn, 'w'), indent=1, default=float)
log('element %d done in %.0fs' % (k, time.time() - t00))
