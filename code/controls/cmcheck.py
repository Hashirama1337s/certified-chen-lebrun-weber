"""BAR-CTRL1 containment check (check logic as in the referee's REFEREE-CM run).
CM model value P(pt) must satisfy |P - truth|_lower <= rem at every grid point; truth = exact square values, or arbrho for CLW.
python cmcheck.py CASE DEG NG OUT.json
  env CM_MUT = none     -> the SHIPPING pinned seatE/cm.py is imported (hash-checked)
  env CM_MUT = notail | noinvrem -> the mutated COPY in controls/mut/cm.py is imported (pinned file untouched)"""
import sys, os, time, json, hashlib
HERE = os.path.dirname(os.path.abspath(__file__))
SE = os.path.join(HERE, '..', 'seatE')
PINS = {'cm.py': '8697432728c08f732418af0f17b4a3cc123096cbd796c2f9dad989c040542025',
        'cmtoric.py': 'b690368898de618d2e54a93a56005c6b7ed1c263025e6b4021ced0675fafb13f',
        'arbrho.py': 'add126b7bcc071db59f0aa35fc40c4f5e8b0c8e1a1dfcc07d1c34d1784909b2a',
        'common.py': '2c558a5f8e43e739c694f4f24b39f01dd5b511e67a941e48de3006d4bbe113b2'}


def pin_check():
    for f, h in PINS.items():
        assert hashlib.sha256(open(os.path.join(SE, f), 'rb').read()).hexdigest() == h, 'pinned file changed: ' + f
    return True


pin_check()
MUTENV = os.environ.get('CM_MUT', 'none')
sys.path.insert(0, SE)
if MUTENV != 'none':
    sys.path.insert(0, os.path.join(HERE, 'mut'))      # mutated copy wins
import cm
if MUTENV == 'none':
    assert os.path.normcase(os.path.dirname(os.path.abspath(cm.__file__))) == os.path.normcase(os.path.abspath(SE)), cm.__file__
    assert not hasattr(cm, 'MUT')
else:
    assert os.path.normcase(os.path.dirname(os.path.abspath(cm.__file__))) == os.path.normcase(os.path.join(HERE, 'mut')), cm.__file__
    assert cm.MUT == MUTENV
from flint import arb, ctx
import common, cmtoric

case, deg, NG, OUT = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
cm.set_degree(deg)
PREC = 192
print('MUT', MUTENV, 'case', case, 'deg', deg, 'NG', NG, 'cm from', cm.__file__, flush=True)
T0 = time.time()


def peval(F, u, v):
    C = F.cheb()
    D = cm.D
    Tu = [arb(1), u]; Tv = [arb(1), v]
    for k in range(2, D + 1):
        Tu.append(2 * u * Tu[-1] - Tu[-2]); Tv.append(2 * v * Tv[-1] - Tv[-2])
    tot = arb(0)
    for i in range(D + 1):
        for j in range(D + 1):
            tot += C[i][j] * Tu[i] * Tv[j]
    return tot


rows = []


def check(sub, name, F, pts, truth):
    worst = -1e300; nfail = 0
    r = F.rem
    for (u, v), tr in zip(pts, truth):
        p = peval(F, u, v)
        err = (p - tr).abs_lower()                 # certain lower bound of |P - truth|
        ratio = float(err / r) if r > 0 else (float('inf') if err > 0 else 0.0)
        worst = max(worst, ratio)
        if err > r:
            nfail += 1
    print('  %-6s rem %.3e  max(|P-truth|_lower / rem) = %.3e  violations %d/%d' % (name, float(r.abs_upper()), worst, nfail, len(pts)), flush=True)
    rows.append(dict(sub=sub, quantity=name, rem=float(r.abs_upper()), max_ratio=worst, violations=nfail, points=len(pts)))
    return nfail


if case == 'square':
    prob = cmtoric.AProb('square', None, prec=PREC)
    subs = [(0, 0.0, 0.0, 0.25), (0, 0.5, 0.5, 0.25)]
else:
    prob = cmtoric.AProb('clw', os.path.join(common.SEATR, 'sol_clw_B448_M84.json'), prec=PREC)
    subs = [(0, 0.5, 0.5, 0.25)]
tot_fail = 0
bounds_s = []
for (k, s0, t0, h) in subs:
    ctx.prec = PREC
    t1 = time.time()
    B_ = cmtoric.build(prob, k, s0, t0, h, K=4)
    Q = cmtoric.quantities(B_)
    sub = 'e%d [%g,%g]x[%g,%g]' % (k, s0, s0 + h, t0, t0 + h)
    print('sub %s  build+quant %.0fs' % (sub, time.time() - t1), flush=True)
    g = [-1 + 2 * (i + 0.5) / NG for i in range(NG)] + [-1.0, 1.0]
    pts = [(arb(a), arb(b)) for a in g for b in g]
    if case == 'square':
        n = len(pts)
        tot_fail += check(sub, 's', Q['s'], pts, [arb(4)] * n)
        tot_fail += check(sub, 'rho', Q['rho'], pts, [arb(0)] * n)
        tot_fail += check(sub, 'Riem2', Q['Riem2'], pts, [arb(8)] * n)
        tot_fail += check(sub, 'Ric2', Q['Ric2'], pts, [arb(4)] * n)
        tot_fail += check(sub, 'detA', Q['detA'], pts, [arb(1)] * n)
        b = Q['s'].bounds(); print('  bounds(s) = [%.15f, %.15f]  (square formula control sees only this)' % b)
        bounds_s.append([sub, b[0], b[1]])
    else:
        import arbrho
        from met1_core import obj
        ss = [arb(s0) + arb(h) / 2 * (1 + u) for u, v in pts]
        tt = [arb(t0) + arb(h) / 2 * (1 + v) for u, v in pts]
        ctx.prec = PREC
        t2 = time.time()
        gj, psi, Sc = arbrho.build(prob, k, obj(ss), obj(tt), K=4)
        print('  arbrho %.0fs' % (time.time() - t2), flush=True)
        sc = [Sc.c[0][i] for i in range(len(pts))]
        ps = [psi.c[0][i] for i in range(len(pts))]
        print('  sample: CM s(P) %s  arbrho Scal %s' % (peval(Q['s'], *pts[0]).str(25), sc[0].str(25)))
        tot_fail += check(sub, 's', Q['s'], pts, sc)
        tot_fail += check(sub, 'rho', Q['rho'], pts, ps)
pin_check()
wall = time.time() - T0
print('TOTAL violations', tot_fail, 'wall %.0fs' % wall)
json.dump(dict(case=case, deg=deg, mut=MUTENV, NG=NG, prec=PREC, points_per_sub=(NG + 2) ** 2, total_violations=tot_fail,
               max_ratio=max(r['max_ratio'] for r in rows), wall_s=round(wall, 1), rows=rows, bounds_s=bounds_s,
               cm_file=os.path.relpath(cm.__file__, os.path.join(HERE, '..')), pinned_files_unchanged=True),
          open(OUT, 'w'), indent=1)
