"""audit_torus.py -- amend3 audit (BAR-MET3-amend3): were any crude bounds of the already-written sub-squares affected by the
amend2 NaN drop (Python max() keeps the old value when the new one is NaN)?
python audit_torus.py A|B ELEMENT DIR [subs]
Recomputes ONLY the torus crude bounds M of each written sub-square (status ok) with the amend3 code (identical to amend2 except
NaN -> inf) in the chart recorded in its file, and compares with the M stored in the file, key by key, EXACTLY.
  all equal  -> no NaN occurred there: the stored certificate is what amend3 computes (bit-identical), it stands.
  any key inf now but finite in the file -> a NaN was dropped: the sub-square must be deleted and recomputed under amend3.
  any other difference -> reported (nondeterminism); the lead decides.
Files with chart 'quarters' (written by amend3 itself) are skipped."""
import sys, os, json, math, glob
import numpy as np
import common, cmtoric, fbox, jets, certify

part, k, d = sys.argv[1], int(sys.argv[2]), sys.argv[3]
subs = sys.argv[4] if len(sys.argv) > 4 else 'all'
H, RHO = 0.25, 2.0
prob = cmtoric.AProb('clw', os.path.join(common.SEATR, 'sol_clw_B448_M84.json'), prec=192)
bank = fbox.PatchBank(prob, h=0.1, rp=0.4, m=120, R=1.2, kd=12)
model = jets.model_from_prob(prob, bank, cpow=0)
pat = ('clw_e%d_*_*.json' if part == 'A' else 'clwB_e%d_*_*.json') % k
files = sorted(glob.glob(os.path.join(d, pat)))
bad = 0
for fn in files:
    r = json.load(open(fn))
    i, j = r['i'], r['j']
    if subs != 'all' and '%d:%d' % (i, j) not in subs.split(','):
        continue
    rec = r['B3'] if part == 'A' else r
    if rec.get('status') != 'ok' or rec.get('chart') == 'quarters':
        print('%s: skipped (status %s, chart %s)' % (os.path.basename(fn), rec.get('status'), rec.get('chart')), flush=True); continue
    chart = rec['chart']
    narc = 100 if chart not in ('polar', 'xchart') else 200
    sub = certify.Sub(prob, model, k, i * H, j * H, H, chart)
    if part == 'B':
        M, _, _ = certify.torus_max(sub, RHO, narc)
    else:
        U = certify.ellipse_arcs(RHO, narc // 2, True); V = certify.ellipse_arcs(RHO, narc, False)
        S, T = sub.st(*certify.prod_boxes(U, V))
        M = {}
        for o in sub.run(S, T, want=('abreu',)):
            for jj in range(4):
                M['rho%d' % jj] = max(M.get('rho%d' % jj, 0.0), certify._ub(o['abreu'][jj]))
    old = rec['M']
    diff = {kk: (old.get(kk), v) for kk, v in M.items() if not (old.get(kk) == v)}
    dropped = [kk for kk, (a, b) in diff.items() if b == math.inf and (a is not None and math.isfinite(a))]
    verdict = 'EQUAL (stands)' if not diff else ('NaN DROPPED in %s -> RECOMPUTE' % dropped if dropped else 'DIFFERENT %s' % diff)
    bad += bool(diff)
    print('%s chart %s: %s' % (os.path.basename(fn), chart, verdict), flush=True)
print('audit %s e%d: %d file(s) differ' % (part, k, bad), flush=True)
